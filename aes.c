/*
 *  SPDX-License-Identifier: MIT
 */

// AES CTR mode 128/192/256
// Reference - https://nvlpubs.nist.gov/nistpubs/fips/nist.fips.197.pdf
// Tested against Appendix B - Cipher Example

#include "aes.h"

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>

#include "fields.h"

typedef bf8_t state_t[4][4];

static const uint8_t nc    = 4;       // # of columns
static const uint8_t nr    = 4;       // # of rows
static const uint8_t bsize = nc * nr; // block size

uint8_t kwords; // key words
uint8_t nround; // no round
uint16_t seclv; // 128,192,256

ATTR_CONST static bf8_t get_bit(bf8_t in, uint8_t index) {
  return (in >> index) & 0x01;
}

ATTR_CONST static bf8_t set_bit(bf8_t in, uint8_t index) {
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

void aes_increment_iv(bf8_t* iv) {
  for (uint8_t i = 16; i > 0; i--) {
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
static void add_round_key(bf8_t round, state_t* state, bf8_t* roundKey) {
  for (uint8_t c = 0; c < nc; c++) {
    for (uint8_t r = 0; r < nr; r++) {
      (*state)[c][r] ^= roundKey[(nc * round * 4) + (c * nc) + r];
    }
  }
}
static void sub_bytes(state_t* state) {
  for (uint8_t c = 0; c < nc; c++) {
    for (uint8_t r = 0; r < nr; r++) {
      (*state)[c][r] = compute_sbox((*state)[c][r]);
    }
  }
}
static void shift_row(state_t* state) {
  bf8_t tmp      = (*state)[0][1];
  (*state)[0][1] = (*state)[1][1];
  (*state)[1][1] = (*state)[2][1];
  (*state)[2][1] = (*state)[3][1];
  (*state)[3][1] = tmp;

  bf8_t tmp_20   = (*state)[0][2];
  bf8_t tmp_21   = (*state)[1][2];
  (*state)[0][2] = (*state)[2][2];
  (*state)[1][2] = (*state)[3][2];
  (*state)[2][2] = tmp_20;
  (*state)[3][2] = tmp_21;

  tmp            = (*state)[0][3];
  (*state)[0][3] = (*state)[3][3];
  (*state)[3][3] = (*state)[2][3];
  (*state)[2][3] = (*state)[1][3];
  (*state)[1][3] = tmp;
}

static void mix_column(state_t* state) {
  for (uint8_t c = 0; c < nc; c++) {
    bf8_t tmp = bf8_mul((*state)[c][0], 0x02) ^ bf8_mul((*state)[c][1], 0x03) ^ (*state)[c][2] ^
                (*state)[c][3];
    bf8_t tmp_1 = (*state)[c][0] ^ bf8_mul((*state)[c][1], 0x02) ^ bf8_mul((*state)[c][2], 0x03) ^
                  (*state)[c][3];
    bf8_t tmp_2 = (*state)[c][0] ^ (*state)[c][1] ^ bf8_mul((*state)[c][2], 0x02) ^
                  bf8_mul((*state)[c][3], 0x03);
    bf8_t tmp_3 = bf8_mul((*state)[c][0], 0x03) ^ (*state)[c][1] ^ (*state)[c][2] ^
                  bf8_mul((*state)[c][3], 0x02);

    (*state)[c][0] = tmp;
    (*state)[c][1] = tmp_1;
    (*state)[c][2] = tmp_2;
    (*state)[c][3] = tmp_3;
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

static bf8_t round_const(uint8_t in) {
  bf8_t tmp = 1;
  for (uint8_t i = 1; i < in; i++) {
    tmp = bf8_mul(tmp, 0x02);
  }
  return tmp;
}

// Main Functions
static void key_expansion(bf8_t* key, bf8_t* roundKey) {
  for (uint8_t k = 0; k < kwords; k++) {
    roundKey[4 * k]       = key[4 * k];
    roundKey[(4 * k) + 1] = key[(4 * k) + 1];
    roundKey[(4 * k) + 2] = key[(4 * k) + 2];
    roundKey[(4 * k) + 3] = key[(4 * k) + 3];
  }
  bf8_t tmp[4];
  for (uint8_t k = kwords; k < nc * (nround + 1); k++) {
    uint8_t j = (k - 1) * 4;
    tmp[0]    = roundKey[j];
    tmp[1]    = roundKey[j + 1];
    tmp[2]    = roundKey[j + 2];
    tmp[3]    = roundKey[j + 3];

    if (k % kwords == 0) {
      rot_word(tmp);
      sub_words(tmp);
      tmp[0] ^= round_const(k / kwords);
    }

    if (seclv == 256 && (k % kwords) == 4) {
      sub_words(tmp);
    }

    j               = k * 4;
    uint8_t m       = (k - kwords) * 4;
    roundKey[j]     = roundKey[m] ^ tmp[0];
    roundKey[j + 1] = roundKey[m + 1] ^ tmp[1];
    roundKey[j + 2] = roundKey[m + 2] ^ tmp[2];
    roundKey[j + 3] = roundKey[m + 3] ^ tmp[3];
  }
}

static void cipher(state_t* state, bf8_t* roundKey) {
  uint8_t round = 0;

  // first round
  add_round_key(round, state, roundKey);

  for (round = 1; round < nround; round++) {
    sub_bytes(state);
    shift_row(state);
    mix_column(state);
    add_round_key(round, state, roundKey);
  }

  // last round
  sub_bytes(state);
  shift_row(state);
  add_round_key(round, state, roundKey);
}

// Calling Functions
void aes_ctr_encrypt(bf8_t* key, bf8_t* iv, bf8_t* plaintext, bf8_t* output, uint16_t seclv_) {
  seclv = seclv_;
  if (seclv_ == 256) {
    kwords = 8;
    nround = 14;
    bf8_t round_key[16 * (14 + 1)];
    key_expansion(key, round_key);
    cipher((state_t*)iv, round_key);
    for (uint8_t i = 0; i < 16; i++) {
      output[i] = plaintext[i] ^ iv[i];
    }
  } else if (seclv_ == 192) {
    kwords = 6;
    nround = 12;
    bf8_t round_key[16 * (12 + 1)];
    key_expansion(key, round_key);
    cipher((state_t*)iv, round_key);
    for (uint8_t i = 0; i < 16; i++) {
      output[i] = plaintext[i] ^ iv[i];
    }
  } else {
    kwords = 4;
    nround = 10;
    bf8_t round_key[16 * (10 + 1)];
    key_expansion(key, round_key);
    cipher((state_t*)iv, round_key);
    for (uint8_t i = 0; i < 16; i++) {
      output[i] = plaintext[i] ^ iv[i];
    }
  }
}
