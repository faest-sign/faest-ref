/*
 *  SPDX-License-Identifier: MIT
 */

// AES CTR mode 128/192/256
// Reference - https://nvlpubs.nist.gov/nistpubs/fips/nist.fips.197.pdf
// Tested against Appendix B - Cipher Example

#include "aes.h"

#include "fields.h"
#include "compat.h"

#include <string.h>

// # of rows
#define AES_NR 4

// block with 4 (AES) up to 8 (Rijndael-256) units
typedef aes_word_t aes_block_t[8];

static const bf8_t round_constants[30] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a,
    0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91};

ATTR_CONST static inline bf8_t get_bit(bf8_t in, uint8_t index) {
  return (in >> index) & 0x01;
}

ATTR_CONST static inline bf8_t set_bit(bf8_t in, uint8_t index) {
  return (in << index);
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
    for (unsigned int r = 0; r < AES_NR; r++) {
      state[c][r] ^= round_key->round_keys[round][c][r];
    }
  }
}

static void sub_bytes(aes_block_t state, unsigned int block_words) {
  for (unsigned int c = 0; c < block_words; c++) {
    for (unsigned int r = 0; r < AES_NR; r++) {
      state[c][r] = compute_sbox(state[c][r]);
    }
  }
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
    state[i][0] = new_state[i][0];
    state[i][1] = new_state[i][1];
    state[i][2] = new_state[i][2];
    state[i][3] = new_state[i][3];
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

static bool expand_key(aes_round_keys_t* round_keys, const uint8_t* key, unsigned int key_words,
                       unsigned int block_words, unsigned int num_rounds) {
  int ret             = 0;
  const bf8_t zero[4] = {0};

  for (unsigned int k = 0; k < key_words; k++) {
    round_keys->round_keys[k / block_words][k % block_words][0] = bf8_load(&key[4 * k]);
    round_keys->round_keys[k / block_words][k % block_words][1] = bf8_load(&key[(4 * k) + 1]);
    round_keys->round_keys[k / block_words][k % block_words][2] = bf8_load(&key[(4 * k) + 2]);
    round_keys->round_keys[k / block_words][k % block_words][3] = bf8_load(&key[(4 * k) + 3]);
  }

  for (unsigned int k = key_words; k < block_words * (num_rounds + 1); ++k) {
    bf8_t tmp[4];
    memcpy(tmp, round_keys->round_keys[(k - 1) / block_words][(k - 1) % block_words], sizeof(tmp));

    if (k % key_words == 0) {
      rot_word(tmp);
      ret |= faest_timingsafe_bcmp(tmp, zero, sizeof(tmp));
      sub_words(tmp);
      tmp[0] ^= round_constants[(k / key_words) - 1];
    }

    if (key_words == 8 && (k % key_words) == 4) {
      ret |= faest_timingsafe_bcmp(tmp, zero, sizeof(tmp));
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

  return ret == 0;
}

// Calling Functions

bool aes128_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, 4, 4, 10);
}

bool aes192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, 6, 4, 12);
}

bool aes256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, 8, 4, 14);
}

bool rijndael192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, 6, 6, 12);
}

bool rijndael256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  return expand_key(round_key, key, 8, 8, 14);
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

static void aes_encrypt(const aes_round_keys_t* keys, aes_block_t state, unsigned int block_words,
                        unsigned int num_rounds) {
  // first round
  add_round_key(0, state, keys, block_words);

  for (unsigned int round = 1; round < num_rounds; ++round) {
    sub_bytes(state, block_words);
    shift_row(state, block_words);
    mix_column(state, block_words);
    add_round_key(round, state, keys, block_words);
  }

  // last round
  sub_bytes(state, block_words);
  shift_row(state, block_words);
  add_round_key(num_rounds, state, keys, block_words);
}

void aes128_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                          uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, 4);
  aes_encrypt(key, state, 4, 10);
  store_state(ciphertext, state, 4);
}

void aes192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                          uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, 4);
  aes_encrypt(key, state, 4, 12);
  store_state(ciphertext, state, 4);
}

void aes256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                          uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, 4);
  aes_encrypt(key, state, 4, 14);
  store_state(ciphertext, state, 4);
}

void rijndael192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                               uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, 6);
  aes_encrypt(key, state, 6, 12);
  store_state(ciphertext, state, 6);
}

void rijndael256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                               uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, 8);
  aes_encrypt(key, state, 8, 14);
  store_state(ciphertext, state, 8);
}

void aes128_ctr_encrypt(const aes_round_keys_t* key, const uint8_t* iv, const uint8_t* plaintext,
                        uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, iv, 4);
  aes_encrypt(key, state, 4, 10);

  for (unsigned int i = 0; i < 16; ++i) {
    ciphertext[i] = plaintext[i] ^ state[i / 4][i % 4];
  }
}

void aes192_ctr_encrypt(const aes_round_keys_t* key, const uint8_t* iv, const uint8_t* plaintext,
                        uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, iv, 4);
  aes_encrypt(key, state, 4, 12);

  for (unsigned int i = 0; i < 16; ++i) {
    ciphertext[i] = plaintext[i] ^ state[i / 4][i % 4];
  }
}

void aes256_ctr_encrypt(const aes_round_keys_t* key, const uint8_t* iv, const uint8_t* plaintext,
                        uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, iv, 4);
  aes_encrypt(key, state, 4, 14);

  for (unsigned int i = 0; i < 16; ++i) {
    ciphertext[i] = plaintext[i] ^ state[i / 4][i % 4];
  }
}

void aes_prg(const uint8_t* key, uint8_t* iv, uint8_t* out, uint16_t seclvl) {
  aes_round_keys_t round_key;
  aes_block_t state;
  load_state(state, iv, 4);

  switch (seclvl) {
  case 256:
    aes256_init_round_keys(&round_key, key);
    for (uint64_t i = 0; i < 4; i++) {
      aes_encrypt(&round_key, state, 4, 14);
      store_state(out + i * 16, state, 4);
      aes_increment_iv(iv);
    }
    return;
  case 192:
    aes192_init_round_keys(&round_key, key);
    for (uint64_t i = 0; i < 3; i++) {
      aes_encrypt(&round_key, state, 4, 12);
      store_state(out + i * 16, state, 4);
      aes_increment_iv(iv);
    }
    return;
  default:
    aes128_init_round_keys(&round_key, key);
    for (uint64_t i = 0; i < 2; i++) {
      aes_encrypt(&round_key, state, 4, 10);
      store_state(out + i * 16, state, 4);
      aes_increment_iv(iv);
    }
    return;
  }
}
