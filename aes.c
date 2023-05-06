/*
 *  SPDX-License-Identifier: MIT
 */

// AES CTR mode 128/192/256
// Reference - https://nvlpubs.nist.gov/nistpubs/fips/nist.fips.197.pdf
// Tested against Appendix B - Cipher Example

#include "aes.h"

#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>

#include "fields.h"

// # of columns
#define AES_NC 4
// # of rows
#define AES_NR 4
// block size
#define AES_BLOCK_SIZE (AES_NC * AES_NR)

typedef bf8_t state_t[AES_NC][AES_NR];

static const bf8_t round_constants[10] = {0x01, 0x02, 0x04, 0x08, 0x10,
                                          0x20, 0x40, 0x80, 0x1b, 0x36};

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

void increment_iv(uint8_t* iv) {
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
static void add_round_key(unsigned int round, state_t state, const aes_round_key_t* round_key) {
  for (unsigned int c = 0; c < AES_NC; c++) {
    for (unsigned int r = 0; r < AES_NR; r++) {
      state[c][r] ^= round_key[round][c][r];
    }
  }
}
static void sub_bytes(state_t state) {
  for (unsigned int c = 0; c < AES_NC; c++) {
    for (unsigned int r = 0; r < AES_NR; r++) {
      state[c][r] = compute_sbox(state[c][r]);
    }
  }
}
static void shift_row(state_t state) {
  bf8_t tmp   = state[0][1];
  state[0][1] = state[1][1];
  state[1][1] = state[2][1];
  state[2][1] = state[3][1];
  state[3][1] = tmp;

  bf8_t tmp_20 = state[0][2];
  bf8_t tmp_21 = state[1][2];
  state[0][2]  = state[2][2];
  state[1][2]  = state[3][2];
  state[2][2]  = tmp_20;
  state[3][2]  = tmp_21;

  tmp         = state[0][3];
  state[0][3] = state[3][3];
  state[3][3] = state[2][3];
  state[2][3] = state[1][3];
  state[1][3] = tmp;
}

static void mix_column(state_t state) {
  for (unsigned int c = 0; c < AES_NC; c++) {
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

static void expand_key(aes_round_key_t* round_keys, const uint8_t* key, unsigned int kwords,
                       unsigned int nround, bool do_sub_words) {
  for (unsigned int k = 0; k < kwords; k++) {
    round_keys[0][k][0] = bf8_load(&key[4 * k]);
    round_keys[0][k][1] = bf8_load(&key[(4 * k) + 1]);
    round_keys[0][k][2] = bf8_load(&key[(4 * k) + 2]);
    round_keys[0][k][3] = bf8_load(&key[(4 * k) + 3]);
  }

  for (unsigned int k = kwords; k < AES_NC * (nround + 1); ++k) {
    bf8_t tmp[4];
    memcpy(tmp, round_keys[(k - 1) / 4][(k - 1) % 4], sizeof(tmp));

    if (k % kwords == 0) {
      rot_word(tmp);
      sub_words(tmp);
      tmp[0] ^= round_constants[(k / kwords) - 1];
    }

    if (do_sub_words && (k % kwords) == 4) {
      sub_words(tmp);
    }

    unsigned int m              = k - kwords;
    round_keys[k / 4][k % 4][0] = round_keys[m / 4][m % 4][0] ^ tmp[0];
    round_keys[k / 4][k % 4][1] = round_keys[m / 4][m % 4][1] ^ tmp[1];
    round_keys[k / 4][k % 4][2] = round_keys[m / 4][m % 4][2] ^ tmp[2];
    round_keys[k / 4][k % 4][3] = round_keys[m / 4][m % 4][3] ^ tmp[3];
  }
}

// Calling Functions

void aes128_init_round_keys(aes128_round_keys_t* round_key, const uint8_t* key) {
  expand_key(round_key->keys, key, 4, 10, false);
}

void aes192_init_round_keys(aes192_round_keys_t* round_key, const uint8_t* key) {
  expand_key(round_key->keys, key, 6, 12, false);
}

void aes256_init_round_keys(aes256_round_keys_t* round_key, const uint8_t* key) {
  expand_key(round_key->keys, key, 8, 14, true);
}

static void load_state(state_t iv, const uint8_t* src) {
  for (unsigned int i = 0; i != 16; ++i) {
    iv[i / 4][i % 4] = bf8_load(&src[i]);
  }
}

static void aes_encrypt(const aes_round_key_t* keys, state_t state, unsigned int num_rounds) {
  // first round
  add_round_key(0, state, keys);

  for (unsigned int round = 1; round < num_rounds; ++round) {
    sub_bytes(state);
    shift_row(state);
    mix_column(state);
    add_round_key(round, state, keys);
  }

  // last round
  sub_bytes(state);
  shift_row(state);
  add_round_key(num_rounds, state, keys);
}

void aes128_ctr_encrypt(const uint8_t* key, const uint8_t* iv,
                        const uint8_t* plaintext, uint8_t* ciphertext) {
  state_t state;
  load_state(state, iv);
  aes128_round_keys_t round_key;
  aes128_init_round_keys(&round_key, key);
  aes_encrypt(round_key.keys, state, 10);

  for (unsigned int i = 0; i < 16; ++i) {
    ciphertext[i] = plaintext[i] ^ state[i / 4][i % 4];
  }
}
void aes192_ctr_encrypt(const uint8_t* key, const uint8_t* iv,
                        const uint8_t* plaintext, uint8_t* ciphertext) {
  state_t state;
  load_state(state, iv);
  aes192_round_keys_t round_key;
  aes192_init_round_keys(&round_key, key);
  aes_encrypt(round_key.keys, state, 12);

  for (unsigned int i = 0; i < 16; ++i) {
    ciphertext[i] = plaintext[i] ^ state[i / 4][i % 4];
  }
}
void aes256_ctr_encrypt(const uint8_t* key, const uint8_t* iv,
                        const uint8_t* plaintext, uint8_t* ciphertext) {
  state_t state;
  load_state(state, iv);
  aes256_round_keys_t round_key;
  aes256_init_round_keys(&round_key, key);
  aes_encrypt(round_key.keys, state, 14);

  for (unsigned int i = 0; i < 16; ++i) {
    ciphertext[i] = plaintext[i] ^ state[i / 4][i % 4];
  }
}


void aes_prg(const uint8_t* key, uint8_t* iv, uint8_t* out, uint16_t seclvl) {
  uint32_t outlenbits = seclvl*2;
  uint64_t k = 0;
  switch (seclvl) {
    case 256:
      for(uint64_t i = 0; i < ceil(outlenbits/256); i++) { 
        state_t state;
        load_state(state, iv);
        aes256_round_keys_t round_key;
        aes256_init_round_keys(&round_key, key);
        aes_encrypt(round_key.keys, state, 14);
        for (unsigned int j = 0; j < 16; ++j) {
          out[k] = state[j / 4][j % 4];
          k++;
        }
        increment_iv(iv);
      }
      break;
    case 192:
      for(uint64_t i = 0; i < ceil(outlenbits/192); i++) { 
        state_t state;
        load_state(state, iv);
        aes192_round_keys_t round_key;
        aes192_init_round_keys(&round_key, key);
        aes_encrypt(round_key.keys, state, 12);
        for (unsigned int j = 0; j < 16; ++j) {
          out[k] = state[j / 4][j % 4];
          k++;
        }
        increment_iv(iv);
      }
      break;
    default:
      for(uint64_t i = 0; i < ceil(outlenbits/128); i++) { 
        state_t state;
        load_state(state, iv);
        aes128_round_keys_t round_key;
        aes128_init_round_keys(&round_key, key);
        aes_encrypt(round_key.keys, state, 10);
        for (unsigned int j = 0; j < 16; ++j) {
          out[k] = state[j / 4][j % 4];
          k++;
        }
        increment_iv(iv);
      }
      break;
  }
}