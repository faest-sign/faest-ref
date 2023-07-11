/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "aes.h"

#include "fields.h"
#include "compat.h"
#include "utils.h"

#if defined(HAVE_OPENSSL)
#include <openssl/evp.h>
#endif
#include <string.h>

#define ROUNDS_128 10
#define ROUNDS_192 12
#define ROUNDS_256 14

#define KEY_WORDS_128 4
#define KEY_WORDS_192 6
#define KEY_WORDS_256 8

#define AES_BLOCK_WORDS 4
#define RIJNDAEL_BLOCK_WORDS_192 6
#define RIJNDAEL_BLOCK_WORDS_256 8

static const bf8_t round_constants[30] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a,
    0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91,
};

static int contains_zero(const bf8_t* block) {
  return !block[0] | !block[1] | !block[2] | !block[3];
}

static bf8_t compute_sbox(bf8_t in) {
  bf8_t t  = bf8_inv(in);
  bf8_t t0 = set_bit(
      get_bit(t, 0) ^ get_bit(t, 4) ^ get_bit(t, 5) ^ get_bit(t, 6) ^ get_bit(t, 7) ^ 0x01, 0);
  bf8_t t1 = set_bit(
      get_bit(t, 0) ^ get_bit(t, 1) ^ get_bit(t, 5) ^ get_bit(t, 6) ^ get_bit(t, 7) ^ 0x01, 1);
  bf8_t t2 =
      set_bit(get_bit(t, 0) ^ get_bit(t, 1) ^ get_bit(t, 2) ^ get_bit(t, 6) ^ get_bit(t, 7), 2);
  bf8_t t3 =
      set_bit(get_bit(t, 0) ^ get_bit(t, 1) ^ get_bit(t, 2) ^ get_bit(t, 3) ^ get_bit(t, 7), 3);
  bf8_t t4 =
      set_bit(get_bit(t, 0) ^ get_bit(t, 1) ^ get_bit(t, 2) ^ get_bit(t, 3) ^ get_bit(t, 4), 4);
  bf8_t t5 = set_bit(
      get_bit(t, 1) ^ get_bit(t, 2) ^ get_bit(t, 3) ^ get_bit(t, 4) ^ get_bit(t, 5) ^ 0x01, 5);
  bf8_t t6 = set_bit(
      get_bit(t, 2) ^ get_bit(t, 3) ^ get_bit(t, 4) ^ get_bit(t, 5) ^ get_bit(t, 6) ^ 0x01, 6);
  bf8_t t7 =
      set_bit(get_bit(t, 3) ^ get_bit(t, 4) ^ get_bit(t, 5) ^ get_bit(t, 6) ^ get_bit(t, 7), 7);
  return t0 ^ t1 ^ t2 ^ t3 ^ t4 ^ t5 ^ t6 ^ t7;
}

void aes_increment_iv(uint8_t* iv) {
  for (unsigned int i = 16; i > 0; i--) {
    if (iv[i - 1] == 0xff) {
      iv[i - 1] = 0x00;
      continue;
    }
    iv[i - 1] += 0x01;
    break;
  }
}

// ## AES ##
// Round Functions
static void add_round_key(unsigned int round, aes_block_t state, const aes_round_keys_t* round_key,
                          unsigned int block_words) {
  for (unsigned int c = 0; c < block_words; c++) {
    xor_u8_array(&state[c][0], &round_key->round_keys[round][c][0], &state[c][0], AES_NR);
  }
}

static int sub_bytes(aes_block_t state, unsigned int block_words) {
  int ret = 0;

  for (unsigned int c = 0; c < block_words; c++) {
    ret |= contains_zero(&state[c][0]);
    for (unsigned int r = 0; r < AES_NR; r++) {
      state[c][r] = compute_sbox(state[c][r]);
    }
  }

  return ret;
}

static void shift_row(aes_block_t state, unsigned int block_words) {
  aes_block_t new_state;
  switch (block_words) {
  case 4:
  case 6:
    for (unsigned int i = 0; i < block_words; ++i) {
      new_state[i][0] = state[i][0];
      new_state[i][1] = state[(i + 1) % block_words][1];
      new_state[i][2] = state[(i + 2) % block_words][2];
      new_state[i][3] = state[(i + 3) % block_words][3];
    }
    break;
  case 8:
    for (unsigned int i = 0; i < block_words; i++) {
      new_state[i][0] = state[i][0];
      new_state[i][1] = state[(i + 1) % 8][1];
      new_state[i][2] = state[(i + 3) % 8][2];
      new_state[i][3] = state[(i + 4) % 8][3];
    }
    break;
  }

  for (unsigned int i = 0; i < block_words; ++i) {
    memcpy(&state[i][0], &new_state[i][0], AES_NR);
  }
}

static void mix_column(aes_block_t state, unsigned int block_words) {
  for (unsigned int c = 0; c < block_words; c++) {
    bf8_t tmp = bf8_mul(state[c][0], 0x02) ^ bf8_mul(state[c][1], 0x03) ^ state[c][2] ^ state[c][3];
    bf8_t tmp_1 =
        state[c][0] ^ bf8_mul(state[c][1], 0x02) ^ bf8_mul(state[c][2], 0x03) ^ state[c][3];
    bf8_t tmp_2 =
        state[c][0] ^ state[c][1] ^ bf8_mul(state[c][2], 0x02) ^ bf8_mul(state[c][3], 0x03);
    bf8_t tmp_3 =
        bf8_mul(state[c][0], 0x03) ^ state[c][1] ^ state[c][2] ^ bf8_mul(state[c][3], 0x02);

    state[c][0] = tmp;
    state[c][1] = tmp_1;
    state[c][2] = tmp_2;
    state[c][3] = tmp_3;
  }
}

// Key Expansion functions
static void sub_words(bf8_t* words) {
  words[0] = compute_sbox(words[0]);
  words[1] = compute_sbox(words[1]);
  words[2] = compute_sbox(words[2]);
  words[3] = compute_sbox(words[3]);
}

static void rot_word(bf8_t* words) {
  bf8_t tmp = words[0];
  words[0]  = words[1];
  words[1]  = words[2];
  words[2]  = words[3];
  words[3]  = tmp;
}

int expand_key(aes_round_keys_t* round_keys, const uint8_t* key, unsigned int key_words,
               unsigned int block_words, unsigned int num_rounds) {
  int ret = 0;

  for (unsigned int k = 0; k < key_words; k++) {
    round_keys->round_keys[k / block_words][k % block_words][0] = bf8_load(&key[4 * k]);
    round_keys->round_keys[k / block_words][k % block_words][1] = bf8_load(&key[(4 * k) + 1]);
    round_keys->round_keys[k / block_words][k % block_words][2] = bf8_load(&key[(4 * k) + 2]);
    round_keys->round_keys[k / block_words][k % block_words][3] = bf8_load(&key[(4 * k) + 3]);
  }

  for (unsigned int k = key_words; k < block_words * (num_rounds + 1); ++k) {
    bf8_t tmp[AES_NR];
    memcpy(tmp, round_keys->round_keys[(k - 1) / block_words][(k - 1) % block_words], sizeof(tmp));

    if (k % key_words == 0) {
      rot_word(tmp);
      ret |= contains_zero(tmp);
      sub_words(tmp);
      tmp[0] ^= round_constants[(k / key_words) - 1];
    }

    if (key_words > 6 && (k % key_words) == 4) {
      ret |= contains_zero(tmp);
      sub_words(tmp);
    }

    unsigned int m = k - key_words;
    round_keys->round_keys[k / block_words][k % block_words][0] =
        round_keys->round_keys[m / block_words][m % block_words][0] ^ tmp[0];
    round_keys->round_keys[k / block_words][k % block_words][1] =
        round_keys->round_keys[m / block_words][m % block_words][1] ^ tmp[1];
    round_keys->round_keys[k / block_words][k % block_words][2] =
        round_keys->round_keys[m / block_words][m % block_words][2] ^ tmp[2];
    round_keys->round_keys[k / block_words][k % block_words][3] =
        round_keys->round_keys[m / block_words][m % block_words][3] ^ tmp[3];
  }

  return ret;
}

// Calling Functions

int aes128_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, KEY_WORDS_128, AES_BLOCK_WORDS, ROUNDS_128);
}

int aes192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, KEY_WORDS_192, AES_BLOCK_WORDS, ROUNDS_192);
}

int aes256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, KEY_WORDS_256, AES_BLOCK_WORDS, ROUNDS_256);
}

int rijndael192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, KEY_WORDS_192, RIJNDAEL_BLOCK_WORDS_192, ROUNDS_192);
}

int rijndael256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, KEY_WORDS_256, RIJNDAEL_BLOCK_WORDS_256, ROUNDS_256);
}

static void load_state(aes_block_t state, const uint8_t* src, unsigned int block_words) {
  for (unsigned int i = 0; i != block_words * 4; ++i) {
    state[i / 4][i % 4] = bf8_load(&src[i]);
  }
}

static void store_state(uint8_t* dst, aes_block_t state, unsigned int block_words) {
  for (unsigned int i = 0; i != block_words * 4; ++i) {
    bf8_store(&dst[i], state[i / 4][i % 4]);
  }
}

static int aes_encrypt(const aes_round_keys_t* keys, aes_block_t state, unsigned int block_words,
                       unsigned int num_rounds) {
  int ret = 0;

  // first round
  add_round_key(0, state, keys, block_words);

  for (unsigned int round = 1; round < num_rounds; ++round) {
    ret |= sub_bytes(state, block_words);
    shift_row(state, block_words);
    mix_column(state, block_words);
    add_round_key(round, state, keys, block_words);
  }

  // last round
  ret |= sub_bytes(state, block_words);
  shift_row(state, block_words);
  add_round_key(num_rounds, state, keys, block_words);

  return ret;
}

int aes128_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, AES_BLOCK_WORDS);
  const int ret = aes_encrypt(key, state, AES_BLOCK_WORDS, ROUNDS_128);
  store_state(ciphertext, state, AES_BLOCK_WORDS);
  return ret;
}

int aes192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, AES_BLOCK_WORDS);
  const int ret = aes_encrypt(key, state, AES_BLOCK_WORDS, ROUNDS_192);
  store_state(ciphertext, state, AES_BLOCK_WORDS);
  return ret;
}

int aes256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, AES_BLOCK_WORDS);
  const int ret = aes_encrypt(key, state, AES_BLOCK_WORDS, ROUNDS_256);
  store_state(ciphertext, state, AES_BLOCK_WORDS);
  return ret;
}

int rijndael192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                              uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, RIJNDAEL_BLOCK_WORDS_192);
  const int ret = aes_encrypt(key, state, RIJNDAEL_BLOCK_WORDS_192, ROUNDS_192);
  store_state(ciphertext, state, RIJNDAEL_BLOCK_WORDS_192);
  return ret;
}

int rijndael256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                              uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, RIJNDAEL_BLOCK_WORDS_256);
  const int ret = aes_encrypt(key, state, RIJNDAEL_BLOCK_WORDS_256, ROUNDS_256);
  store_state(ciphertext, state, RIJNDAEL_BLOCK_WORDS_256);
  return ret;
}

void prg(const uint8_t* key, const uint8_t* iv, uint8_t* out, unsigned int seclvl, size_t outlen) {
#if !defined(HAVE_OPENSSL)
  uint8_t internal_iv[16];
  memcpy(internal_iv, iv, sizeof(internal_iv));

  aes_round_keys_t round_key;

  switch (seclvl) {
  case 256:
    aes256_init_round_keys(&round_key, key);
    for (; outlen >= 16; outlen -= 16, out += 16) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_256);
      store_state(out, state, AES_BLOCK_WORDS);
      aes_increment_iv(internal_iv);
    }
    if (outlen) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_256);
      uint8_t tmp[16];
      store_state(tmp, state, AES_BLOCK_WORDS);
      memcpy(out, tmp, outlen);
    }
    return;
  case 192:
    aes192_init_round_keys(&round_key, key);
    for (; outlen >= 16; outlen -= 16, out += 16) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_192);
      store_state(out, state, AES_BLOCK_WORDS);
      aes_increment_iv(internal_iv);
    }
    if (outlen) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_192);
      uint8_t tmp[16];
      store_state(tmp, state, AES_BLOCK_WORDS);
      memcpy(out, tmp, outlen);
    }
    return;
  default:
    aes128_init_round_keys(&round_key, key);
    for (; outlen >= 16; outlen -= 16, out += 16) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_128);
      store_state(out, state, AES_BLOCK_WORDS);
      aes_increment_iv(internal_iv);
    }
    if (outlen) {
      aes_block_t state;
      load_state(state, internal_iv, AES_BLOCK_WORDS);
      aes_encrypt(&round_key, state, AES_BLOCK_WORDS, ROUNDS_128);
      uint8_t tmp[16];
      store_state(tmp, state, AES_BLOCK_WORDS);
      memcpy(out, tmp, outlen);
    }
    return;
  }
#else
  const EVP_CIPHER* cipher;
  switch (seclvl) {
  case 256:
    cipher = EVP_aes_256_ctr();
    break;
  case 192:
    cipher = EVP_aes_192_ctr();
    break;
  default:
    cipher = EVP_aes_128_ctr();
    break;
  }

  EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
  assert(ctx);

  EVP_EncryptInit_ex(ctx, cipher, NULL, key, iv);

  static const uint8_t plaintext[16] = {0};

  int len = 0;
  for (size_t idx = 0; idx < outlen / 16; idx += 1, out += 16) {
    EVP_EncryptUpdate(ctx, out, &len, plaintext, sizeof(plaintext));
  }
  if (outlen % 16) {
    EVP_EncryptUpdate(ctx, out, &len, plaintext, outlen % 16);
  }
  EVP_EncryptFinal_ex(ctx, out, &len);
  EVP_CIPHER_CTX_free(ctx);
#endif
}

uint8_t* aes_extend_witness(const uint8_t* key, const uint8_t* in, const faest_paramset_t* params) {
  const unsigned int lambda     = params->faest_param.lambda;
  const unsigned int l          = params->faest_param.l;
  const unsigned int L_ke       = params->faest_param.Lke;
  const unsigned int S_ke       = params->faest_param.Ske;
  const unsigned int num_rounds = params->faest_param.R;

  uint8_t* w           = malloc((l + 7) / 8);
  uint8_t* const w_out = w;

  unsigned int block_words = AES_BLOCK_WORDS;
  unsigned int beta        = 1;
  switch (params->faest_paramid) {
  case FAEST_192F:
  case FAEST_192S:
  case FAEST_256F:
  case FAEST_256S:
    beta = 2;
    break;
  case FAEST_EM_192F:
  case FAEST_EM_192S:
    block_words = RIJNDAEL_BLOCK_WORDS_192;
    break;
  case FAEST_EM_256F:
  case FAEST_EM_256S:
    block_words = RIJNDAEL_BLOCK_WORDS_256;
    break;
  default:
    break;
  }

  if (!L_ke) {
    // switch input and key for EM
    const uint8_t* tmp = key;
    key                = in;
    in                 = tmp;
  }

  // Step 3
  aes_round_keys_t round_keys;
  switch (lambda) {
  case 256:
    if (block_words == RIJNDAEL_BLOCK_WORDS_256) {
      rijndael256_init_round_keys(&round_keys, key);
    } else {
      aes256_init_round_keys(&round_keys, key);
    }
    break;
  case 192:
    if (block_words == RIJNDAEL_BLOCK_WORDS_192) {
      rijndael192_init_round_keys(&round_keys, key);
    } else {
      aes192_init_round_keys(&round_keys, key);
    }
    break;
  default:
    aes128_init_round_keys(&round_keys, key);
    break;
  }

  // Step 4
  if (L_ke > 0) {
    // Key schedule constraints only needed for normal AES, not EM variant.
    for (unsigned int i = 0; i != params->faest_param.Nwd; ++i) {
      memcpy(w, round_keys.round_keys[i / 4][i % 4], sizeof(aes_word_t));
      w += sizeof(aes_word_t);
    }

    for (unsigned int j = 0, ik = params->faest_param.Nwd; j < S_ke / 4; ++j) {
      memcpy(w, round_keys.round_keys[ik / 4][ik % 4], sizeof(aes_word_t));
      w += sizeof(aes_word_t);
      ik += lambda == 192 ? 6 : 4;
    }
  } else {
    // saving the OWF key to the extended witness
    memcpy(w, in, lambda / 8);
    w += lambda / 8;
  }

  // Step 10
  for (unsigned b = 0; b < beta; ++b, in += sizeof(aes_word_t) * block_words) {
    // Step 12
    aes_block_t state;
    load_state(state, in, block_words);

    // Step 13
    add_round_key(0, state, &round_keys, block_words);

    for (unsigned int round = 1; round < num_rounds; ++round) {
      // Step 15
      sub_bytes(state, block_words);
      // Step 16
      shift_row(state, block_words);
      // Step 17
      store_state(w, state, block_words);
      w += sizeof(aes_word_t) * block_words;
      // Step 18
      mix_column(state, block_words);
      // Step 19
      add_round_key(round, state, &round_keys, block_words);
    }
    // last round is not commited to, so not computed
  }

  return w_out;
}
