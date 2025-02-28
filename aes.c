/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#else
#include "macros.h"

#if defined(__x86_64__) || defined(__i386__) || defined(_M_IX86) || defined(_M_AMD64)
#if __has_include(<wmmintrin.h>)
#define HAVE_AESNI
#endif
#if __has_include(<immintrin.h>)
#define HAVE_AVX2
#endif

#if GNUC_CHECK(9, 0) || CLANG_CHECK(8)
#define HAVE_MM_LOADU_SI64
#endif
#endif
#endif

#include "aes.h"

#include "cpu.h"
#include "fields.h"
#include "compat.h"
#include "utils.h"

#if defined(HAVE_AESNI)
#include <wmmintrin.h>
#if defined(HAVE_AVX2)
#include <immintrin.h>
#endif
#endif
#if defined(HAVE_OPENSSL)
#include <openssl/evp.h>
#elif defined(_WIN32)
#include <windows.h>
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

static bf8_t compute_sbox(bf8_t in) {
  bf8_t t  = bf8_inv(in);
  bf8_t t0 = set_bit(parity8(t & (1 | (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7))), 0);
  t0 ^= set_bit(parity8(t & (1 | (1 << 1) | (1 << 5) | (1 << 6) | (1 << 7))), 1);
  t0 ^= set_bit(parity8(t & (1 | (1 << 1) | (1 << 2) | (1 << 6) | (1 << 7))), 2);
  t0 ^= set_bit(parity8(t & (1 | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 7))), 3);
  t0 ^= set_bit(parity8(t & (1 | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))), 4);
  t0 ^= set_bit(parity8(t & ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5))), 5);
  t0 ^= set_bit(parity8(t & ((1 << 2) | (1 << 3) | (1 << 4) | (1 << 5) | (1 << 6))), 6);
  t0 ^= set_bit(parity8(t & ((1 << 3) | (1 << 4) | (1 << 5) | (1 << 6) | (1 << 7))), 7);
  return t0 ^ (1 | (1 << 1) | (1 << 5) | (1 << 6));
}

#if defined(FAEST_TESTS)
void aes_increment_iv(uint8_t* iv)
#else
static inline void aes_increment_iv(uint8_t* iv)
#endif
{
  uint32_t iv0;
  memcpy(&iv0, iv, sizeof(uint32_t));
  iv0 = htole32(le32toh(iv0) + 1);
  memcpy(iv, &iv0, sizeof(uint32_t));
}

// ## AES ##
// Round Functions
static void add_round_key(unsigned int round, aes_block_t state, const aes_round_keys_t* round_key,
                          unsigned int block_words) {
  for (unsigned int c = 0; c < block_words; c++) {
    xor_u8_array(&state[c][0], &round_key->round_keys[round][c][0], &state[c][0], AES_NR);
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
    memcpy(&state[i][0], &new_state[i][0], AES_NR);
  }
}

static void mix_column(aes_block_t state, unsigned int block_words) {
  for (unsigned int c = 0; c < block_words; c++) {
    bf8_t tmp[4];
    tmp[0] = bf8_mul(state[c][0], 0x02) ^ bf8_mul(state[c][1], 0x03) ^ state[c][2] ^ state[c][3];
    tmp[1] = state[c][0] ^ bf8_mul(state[c][1], 0x02) ^ bf8_mul(state[c][2], 0x03) ^ state[c][3];
    tmp[2] = state[c][0] ^ state[c][1] ^ bf8_mul(state[c][2], 0x02) ^ bf8_mul(state[c][3], 0x03);
    tmp[3] = bf8_mul(state[c][0], 0x03) ^ state[c][1] ^ state[c][2] ^ bf8_mul(state[c][3], 0x02);

    memcpy(state[c], tmp, sizeof(tmp));
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
#if 0
  bf8_t tmp = words[0];
  words[0]  = words[1];
  words[1]  = words[2];
  words[2]  = words[3];
  words[3]  = tmp;
#else
  // in the most ideal case, this generates a simple rord instruction
  uint32_t w;
  memcpy(&w, words, sizeof(w));
  w = htole32(rotr32(le32toh(w), 8));
  memcpy(words, &w, sizeof(w));
#endif
}

void expand_key(aes_round_keys_t* round_keys, const uint8_t* key, unsigned int key_words,
                unsigned int block_words, unsigned int num_rounds) {
  for (unsigned int k = 0; k < key_words; k++) {
    memcpy(round_keys->round_keys[k / block_words][k % block_words], &key[4 * k], 4);
  }

  for (unsigned int k = key_words; k < block_words * (num_rounds + 1); ++k) {
    bf8_t tmp[AES_NR];
    memcpy(tmp, round_keys->round_keys[(k - 1) / block_words][(k - 1) % block_words], sizeof(tmp));

    if (k % key_words == 0) {
      rot_word(tmp);
      sub_words(tmp);
      tmp[0] ^= round_constants[(k / key_words) - 1];
    }

    if (key_words > 6 && (k % key_words) == 4) {
      sub_words(tmp);
    }

    const unsigned int m = k - key_words;
    xor_u8_array(round_keys->round_keys[m / block_words][m % block_words], tmp,
                 round_keys->round_keys[k / block_words][k % block_words], 4);
  }
}

// Calling Functions

void aes128_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  expand_key(round_key, key, KEY_WORDS_128, AES_BLOCK_WORDS, ROUNDS_128);
}

void aes192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  expand_key(round_key, key, KEY_WORDS_192, AES_BLOCK_WORDS, ROUNDS_192);
}

void aes256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  expand_key(round_key, key, KEY_WORDS_256, AES_BLOCK_WORDS, ROUNDS_256);
}

void rijndael192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  expand_key(round_key, key, KEY_WORDS_192, RIJNDAEL_BLOCK_WORDS_192, ROUNDS_192);
}

void rijndael256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key) {
  expand_key(round_key, key, KEY_WORDS_256, RIJNDAEL_BLOCK_WORDS_256, ROUNDS_256);
}

static void load_state(aes_block_t state, const uint8_t* src, unsigned int block_words) {
  for (unsigned int i = 0; i != block_words * 4; ++i) {
    state[i / 4][i % 4] = bf8_load(&src[i]);
  }
}

ATTR_CONST static bf8_t bf_exp_238(bf8_t x) {
  // 238 == 0b11101110
  bf8_t y = bf8_square(x); // x^2
  x       = bf8_square(y); // x^4
  y       = bf8_mul(x, y);
  x       = bf8_square(x); // x^8
  y       = bf8_mul(x, y);
  x       = bf8_square(x); // x^16
  x       = bf8_square(x); // x^32
  y       = bf8_mul(x, y);
  x       = bf8_square(x); // x^64
  y       = bf8_mul(x, y);
  x       = bf8_square(x); // x^128
  return bf8_mul(x, y);
}

#if !defined(FAEST_TESTS)
static
#endif
    uint8_t
    invnorm(uint8_t in) {
  // instead of computing in^(-17), we calculate in^238
  in = bf_exp_238(in);
  return set_bit(get_bit(in, 0), 0) ^ set_bit(get_bit(in, 6), 1) ^ set_bit(get_bit(in, 7), 2) ^
         set_bit(get_bit(in, 2), 3);
}

static uint8_t* store_invnorm_state(uint8_t* dst, aes_block_t state, unsigned int block_words) {
  for (unsigned int i = 0; i != block_words * 4; i += 2, ++dst) {
    uint8_t normstate_lo = invnorm(state[i / 4][i % 4]);
    uint8_t normstate_hi = invnorm(state[i / 4][(i + 1) % 4]);
    bf8_store(dst, (normstate_hi << 4) | normstate_lo);
  }
  return dst;
}

static uint8_t* store_state(uint8_t* dst, aes_block_t state, unsigned int block_words) {
  for (unsigned int i = 0; i != block_words * 4; ++i, ++dst) {
    bf8_store(dst, state[i / 4][i % 4]);
  }
  return dst;
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
  load_state(state, plaintext, AES_BLOCK_WORDS);
  aes_encrypt(key, state, AES_BLOCK_WORDS, ROUNDS_128);
  store_state(ciphertext, state, AES_BLOCK_WORDS);
}

void aes192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                          uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, AES_BLOCK_WORDS);
  aes_encrypt(key, state, AES_BLOCK_WORDS, ROUNDS_192);
  store_state(ciphertext, state, AES_BLOCK_WORDS);
}

void aes256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                          uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, AES_BLOCK_WORDS);
  aes_encrypt(key, state, AES_BLOCK_WORDS, ROUNDS_256);
  store_state(ciphertext, state, AES_BLOCK_WORDS);
}

void rijndael192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                               uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, RIJNDAEL_BLOCK_WORDS_192);
  aes_encrypt(key, state, RIJNDAEL_BLOCK_WORDS_192, ROUNDS_192);
  store_state(ciphertext, state, RIJNDAEL_BLOCK_WORDS_192);
}

void rijndael256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                               uint8_t* ciphertext) {
  aes_block_t state;
  load_state(state, plaintext, RIJNDAEL_BLOCK_WORDS_256);
  aes_encrypt(key, state, RIJNDAEL_BLOCK_WORDS_256, ROUNDS_256);
  store_state(ciphertext, state, RIJNDAEL_BLOCK_WORDS_256);
}

static void add_to_upper_word(uint8_t* iv, uint32_t tweak) {
  uint32_t iv3;
  memcpy(&iv3, iv + IV_SIZE - sizeof(uint32_t), sizeof(uint32_t));
  iv3 = htole32(le32toh(iv3) + tweak);
  memcpy(iv + IV_SIZE - sizeof(uint32_t), &iv3, sizeof(uint32_t));
}

#if defined(HAVE_AESNI)
#if defined(_MSC_VER)
// workarounds for MSVC
#include <immintrin.h>

#define __m128i_u __m128i
#elif !GNUC_CHECK(7, 0) && !CLANG_CHECK(9)
// workaround for gcc and clang
#define __m128i_u __m128i
#endif

#if !defined(HAVE_MM_LOADU_SI64)
ATTR_ALWAYS_INLINE static inline __m128i _mm_loadu_si64(const void* src) {
  uint64_t u0;
  memcpy(&u0, src, sizeof(u0));
#if !defined(_MSC_VER) || defined(__x86_64__)
  return _mm_set_epi64x(0, u0);
#else
  // MS VC for x86 (untested)
  union {
    uint64_t q;
    uint32_t r[2];
  } u;
  u.q = u0;
  return _mm_setr_epi32(u.r[0], u.r[1], 0, 0);
#endif
}
#endif

ATTR_TARGET_AESNI ATTR_ALWAYS_INLINE static inline __m128i sse2_increment_iv(__m128i iv) {
  return _mm_add_epi32(iv, _mm_set_epi32(0, 0, 0, 1));
}

ATTR_TARGET_AESNI static inline __m128i aes128_keyexpand(__m128i key) {
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  return _mm_xor_si128(key, _mm_slli_si128(key, 4));
}

ATTR_TARGET_AESNI static inline __m128i aes192_keyexpand_2(__m128i key, __m128i key2) {
  key  = _mm_shuffle_epi32(key, 0xff);
  key2 = _mm_xor_si128(key2, _mm_slli_si128(key2, 4));
  return _mm_xor_si128(key, key2);
}

#define KEYEXP128_H(K1, K2, I, S)                                                                  \
  _mm_xor_si128(aes128_keyexpand(K1), _mm_shuffle_epi32(_mm_aeskeygenassist_si128(K2, I), S))

#define KEYEXP128(K, I) KEYEXP128_H(K, K, I, 0xff)
#define KEYEXP192(K1, K2, I) KEYEXP128_H(K1, K2, I, 0x55)
#define KEYEXP192_2(K1, K2) aes192_keyexpand_2(K1, K2)
#define KEYEXP256(K1, K2, I) KEYEXP128_H(K1, K2, I, 0xff)
#define KEYEXP256_2(K1, K2) KEYEXP128_H(K1, K2, 0x00, 0xaa)

ATTR_TARGET_AESNI static void prg_aesni_128(const uint8_t* key, const uint8_t* iv, uint8_t* out,
                                            size_t outlen) {
  __m128i rk[ROUNDS_128 + 1];
  /* 128 bit key setup */
  rk[0]  = _mm_loadu_si128((const __m128i_u*)key);
  rk[1]  = KEYEXP128(rk[0], 0x01);
  rk[2]  = KEYEXP128(rk[1], 0x02);
  rk[3]  = KEYEXP128(rk[2], 0x04);
  rk[4]  = KEYEXP128(rk[3], 0x08);
  rk[5]  = KEYEXP128(rk[4], 0x10);
  rk[6]  = KEYEXP128(rk[5], 0x20);
  rk[7]  = KEYEXP128(rk[6], 0x40);
  rk[8]  = KEYEXP128(rk[7], 0x80);
  rk[9]  = KEYEXP128(rk[8], 0x1B);
  rk[10] = KEYEXP128(rk[9], 0x36);

  __m128i miv = _mm_loadu_si128((const __m128i_u*)iv);
  for (size_t idx = 0; idx < outlen / IV_SIZE; idx += 1, out += IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_128; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_128]);
    _mm_storeu_si128((__m128i_u*)out, m);
    miv = sse2_increment_iv(miv);
  }

  if (outlen % IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_128; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_128]);

    uint8_t last_block[IV_SIZE];
    _mm_storeu_si128((__m128i_u*)last_block, m);

    memcpy(out, last_block, outlen % IV_SIZE);
  }
}

ATTR_TARGET_AESNI static void prg_aesni_192(const uint8_t* key, uint8_t* iv, uint8_t* out,
                                            size_t outlen) {
  __m128i rk[ROUNDS_192 + 1];
  /* 192 bit key setup */
  __m128i temp[2];
  rk[0]   = _mm_loadu_si128((const __m128i_u*)key);
  rk[1]   = _mm_loadu_si64((const __m128i_u*)(key + 16));
  temp[0] = KEYEXP192(rk[0], rk[1], 0x01);
  temp[1] = KEYEXP192_2(temp[0], rk[1]);
  rk[1]   = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(rk[1]), _mm_castsi128_pd(temp[0]), 0));
  rk[2] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp[0]), _mm_castsi128_pd(temp[1]), 1));
  rk[3] = KEYEXP192(temp[0], temp[1], 0x02);
  rk[4] = KEYEXP192_2(rk[3], temp[1]);
  temp[0] = KEYEXP192(rk[3], rk[4], 0x04);
  temp[1] = KEYEXP192_2(temp[0], rk[4]);
  rk[4]   = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(rk[4]), _mm_castsi128_pd(temp[0]), 0));
  rk[5] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp[0]), _mm_castsi128_pd(temp[1]), 1));
  rk[6] = KEYEXP192(temp[0], temp[1], 0x08);
  rk[7] = KEYEXP192_2(rk[6], temp[1]);
  temp[0] = KEYEXP192(rk[6], rk[7], 0x10);
  temp[1] = KEYEXP192_2(temp[0], rk[7]);
  rk[7]   = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(rk[7]), _mm_castsi128_pd(temp[0]), 0));
  rk[8] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp[0]), _mm_castsi128_pd(temp[1]), 1));
  rk[9] = KEYEXP192(temp[0], temp[1], 0x20);
  rk[10]  = KEYEXP192_2(rk[9], temp[1]);
  temp[0] = KEYEXP192(rk[9], rk[10], 0x40);
  temp[1] = KEYEXP192_2(temp[0], rk[10]);
  rk[10] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(rk[10]), _mm_castsi128_pd(temp[0]), 0));
  rk[11] =
      _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp[0]), _mm_castsi128_pd(temp[1]), 1));
  rk[12] = KEYEXP192(temp[0], temp[1], 0x80);

  __m128i miv = _mm_loadu_si128((const __m128i_u*)iv);
  for (size_t idx = 0; idx < outlen / IV_SIZE; idx += 1, out += IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_192; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_192]);
    _mm_storeu_si128((__m128i_u*)out, m);
    miv = sse2_increment_iv(miv);
  }

  if (outlen % IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_192; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_192]);

    uint8_t last_block[IV_SIZE];
    _mm_storeu_si128((__m128i_u*)last_block, m);

    memcpy(out, last_block, outlen % IV_SIZE);
  }
}

ATTR_TARGET_AESNI static void prg_aesni_256(const uint8_t* key, uint8_t* iv, uint8_t* out,
                                            size_t outlen) {
  __m128i rk[ROUNDS_256 + 1];
  /* 256 bit key setup */
  rk[0]  = _mm_loadu_si128((const __m128i_u*)key);
  rk[1]  = _mm_loadu_si128((const __m128i_u*)(key + 16));
  rk[2]  = KEYEXP256(rk[0], rk[1], 0x01);
  rk[3]  = KEYEXP256_2(rk[1], rk[2]);
  rk[4]  = KEYEXP256(rk[2], rk[3], 0x02);
  rk[5]  = KEYEXP256_2(rk[3], rk[4]);
  rk[6]  = KEYEXP256(rk[4], rk[5], 0x04);
  rk[7]  = KEYEXP256_2(rk[5], rk[6]);
  rk[8]  = KEYEXP256(rk[6], rk[7], 0x08);
  rk[9]  = KEYEXP256_2(rk[7], rk[8]);
  rk[10] = KEYEXP256(rk[8], rk[9], 0x10);
  rk[11] = KEYEXP256_2(rk[9], rk[10]);
  rk[12] = KEYEXP256(rk[10], rk[11], 0x20);
  rk[13] = KEYEXP256_2(rk[11], rk[12]);
  rk[14] = KEYEXP256(rk[12], rk[13], 0x40);

  __m128i miv = _mm_loadu_si128((const __m128i_u*)iv);
  for (size_t idx = 0; idx < outlen / IV_SIZE; idx += 1, out += IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_256; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_256]);
    _mm_storeu_si128((__m128i_u*)out, m);
    miv = sse2_increment_iv(miv);
  }

  if (outlen % IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_256; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_256]);

    uint8_t last_block[IV_SIZE];
    _mm_storeu_si128((__m128i_u*)last_block, m);

    memcpy(out, last_block, outlen % IV_SIZE);
  }
}

#if defined(HAVE_AVX2)
ATTR_TARGET_AESNI_AVX static void prg_aesni_avx_128(const uint8_t* key, const uint8_t* iv,
                                                    uint8_t* out, size_t outlen) {
  __m128i rk[ROUNDS_128 + 1];
  /* 128 bit key setup */
  rk[0]  = _mm_loadu_si128((const __m128i_u*)key);
  rk[1]  = KEYEXP128(rk[0], 0x01);
  rk[2]  = KEYEXP128(rk[1], 0x02);
  rk[3]  = KEYEXP128(rk[2], 0x04);
  rk[4]  = KEYEXP128(rk[3], 0x08);
  rk[5]  = KEYEXP128(rk[4], 0x10);
  rk[6]  = KEYEXP128(rk[5], 0x20);
  rk[7]  = KEYEXP128(rk[6], 0x40);
  rk[8]  = KEYEXP128(rk[7], 0x80);
  rk[9]  = KEYEXP128(rk[8], 0x1B);
  rk[10] = KEYEXP128(rk[9], 0x36);

  __m128i miv = _mm_loadu_si128((const __m128i_u*)iv);
  for (size_t idx = 0; idx < outlen / IV_SIZE; idx += 1, out += IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_128; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_128]);
    _mm_storeu_si128((__m128i_u*)out, m);
    miv = sse2_increment_iv(miv);
  }

  if (outlen % IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_128; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_128]);

    uint8_t last_block[IV_SIZE];
    _mm_storeu_si128((__m128i_u*)last_block, m);

    memcpy(out, last_block, outlen % IV_SIZE);
  }
}

ATTR_TARGET_AESNI_AVX static void prg_aesni_avx_192(const uint8_t* key, uint8_t* iv, uint8_t* out,
                                                    size_t outlen) {
  __m128i rk[ROUNDS_192 + 1];
  /* 192 bit key setup */
  __m128i temp[2];
  rk[0]   = _mm_loadu_si128((const __m128i_u*)key);
  rk[1]   = _mm_loadu_si64((const __m128i_u*)(key + 16));
  temp[0] = KEYEXP192(rk[0], rk[1], 0x01);
  temp[1] = KEYEXP192_2(temp[0], rk[1]);
  rk[1]   = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(rk[1]), _mm_castsi128_pd(temp[0]), 0));
  rk[2] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp[0]), _mm_castsi128_pd(temp[1]), 1));
  rk[3] = KEYEXP192(temp[0], temp[1], 0x02);
  rk[4] = KEYEXP192_2(rk[3], temp[1]);
  temp[0] = KEYEXP192(rk[3], rk[4], 0x04);
  temp[1] = KEYEXP192_2(temp[0], rk[4]);
  rk[4]   = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(rk[4]), _mm_castsi128_pd(temp[0]), 0));
  rk[5] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp[0]), _mm_castsi128_pd(temp[1]), 1));
  rk[6] = KEYEXP192(temp[0], temp[1], 0x08);
  rk[7] = KEYEXP192_2(rk[6], temp[1]);
  temp[0] = KEYEXP192(rk[6], rk[7], 0x10);
  temp[1] = KEYEXP192_2(temp[0], rk[7]);
  rk[7]   = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(rk[7]), _mm_castsi128_pd(temp[0]), 0));
  rk[8] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp[0]), _mm_castsi128_pd(temp[1]), 1));
  rk[9] = KEYEXP192(temp[0], temp[1], 0x20);
  rk[10]  = KEYEXP192_2(rk[9], temp[1]);
  temp[0] = KEYEXP192(rk[9], rk[10], 0x40);
  temp[1] = KEYEXP192_2(temp[0], rk[10]);
  rk[10] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(rk[10]), _mm_castsi128_pd(temp[0]), 0));
  rk[11] =
      _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp[0]), _mm_castsi128_pd(temp[1]), 1));
  rk[12] = KEYEXP192(temp[0], temp[1], 0x80);

  __m128i miv = _mm_loadu_si128((const __m128i_u*)iv);
  for (size_t idx = 0; idx < outlen / IV_SIZE; idx += 1, out += IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_192; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_192]);
    _mm_storeu_si128((__m128i_u*)out, m);
    miv = sse2_increment_iv(miv);
  }

  if (outlen % IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_192; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_192]);

    uint8_t last_block[IV_SIZE];
    _mm_storeu_si128((__m128i_u*)last_block, m);

    memcpy(out, last_block, outlen % IV_SIZE);
  }
}

ATTR_TARGET_AESNI_AVX static void prg_aesni_avx_256(const uint8_t* key, uint8_t* iv, uint8_t* out,
                                                    size_t outlen) {
  __m128i rk[ROUNDS_256 + 1];
  /* 256 bit key setup */
  rk[0]  = _mm_loadu_si128((const __m128i_u*)key);
  rk[1]  = _mm_loadu_si128((const __m128i_u*)(key + 16));
  rk[2]  = KEYEXP256(rk[0], rk[1], 0x01);
  rk[3]  = KEYEXP256_2(rk[1], rk[2]);
  rk[4]  = KEYEXP256(rk[2], rk[3], 0x02);
  rk[5]  = KEYEXP256_2(rk[3], rk[4]);
  rk[6]  = KEYEXP256(rk[4], rk[5], 0x04);
  rk[7]  = KEYEXP256_2(rk[5], rk[6]);
  rk[8]  = KEYEXP256(rk[6], rk[7], 0x08);
  rk[9]  = KEYEXP256_2(rk[7], rk[8]);
  rk[10] = KEYEXP256(rk[8], rk[9], 0x10);
  rk[11] = KEYEXP256_2(rk[9], rk[10]);
  rk[12] = KEYEXP256(rk[10], rk[11], 0x20);
  rk[13] = KEYEXP256_2(rk[11], rk[12]);
  rk[14] = KEYEXP256(rk[12], rk[13], 0x40);

  __m128i miv = _mm_loadu_si128((const __m128i_u*)iv);
  for (size_t idx = 0; idx < outlen / IV_SIZE; idx += 1, out += IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_256; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_256]);
    _mm_storeu_si128((__m128i_u*)out, m);
    miv = sse2_increment_iv(miv);
  }

  if (outlen % IV_SIZE) {
    __m128i m = _mm_xor_si128(miv, rk[0]);
    for (unsigned int round = 1; round != ROUNDS_256; ++round) {
      m = _mm_aesenc_si128(m, rk[round]);
    }
    m = _mm_aesenclast_si128(m, rk[ROUNDS_256]);

    uint8_t last_block[IV_SIZE];
    _mm_storeu_si128((__m128i_u*)last_block, m);

    memcpy(out, last_block, outlen % IV_SIZE);
  }
}
#endif
#endif

void prg(const uint8_t* key, const uint8_t* iv, uint32_t tweak, uint8_t* out, unsigned int seclvl,
         size_t outlen) {
  uint8_t internal_iv[IV_SIZE];
  memcpy(internal_iv, iv, IV_SIZE);
  add_to_upper_word(internal_iv, tweak);

#if defined(HAVE_AESNI)
// use AES-NI if possible
#if defined(HAVE_AVX2)
  if (CPU_SUPPORTS_AESNI_AVX) {
    switch (seclvl) {
    case 256:
      prg_aesni_avx_256(key, internal_iv, out, outlen);
      break;
    case 192:
      prg_aesni_avx_192(key, internal_iv, out, outlen);
      break;
    default:
      prg_aesni_avx_128(key, internal_iv, out, outlen);
      break;
    }
    return;
  }
#endif

  if (CPU_SUPPORTS_AESNI) {
    switch (seclvl) {
    case 256:
      prg_aesni_256(key, internal_iv, out, outlen);
      break;
    case 192:
      prg_aesni_192(key, internal_iv, out, outlen);
      break;
    default:
      prg_aesni_128(key, internal_iv, out, outlen);
      break;
    }
    return;
  }
#endif

#if defined(HAVE_OPENSSL)
  const EVP_CIPHER* cipher;
  switch (seclvl) {
  case 256:
    cipher = EVP_aes_256_ecb();
    break;
  case 192:
    cipher = EVP_aes_192_ecb();
    break;
  default:
    cipher = EVP_aes_128_ecb();
    break;
  }
  assert(cipher);

  EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
  assert(ctx);

  EVP_EncryptInit_ex(ctx, cipher, NULL, key, NULL);

  int len = 0;
  for (size_t idx = 0; idx < outlen / IV_SIZE; idx += 1, out += IV_SIZE) {
    int ret = EVP_EncryptUpdate(ctx, out, &len, internal_iv, IV_SIZE);
    assert(ret == 1);
    assert(len == (int)IV_SIZE);
    (void)ret;
    (void)len;
    aes_increment_iv(internal_iv);
  }
  if (outlen % IV_SIZE) {
    uint8_t last_block[IV_SIZE];
    int ret = EVP_EncryptUpdate(ctx, last_block, &len, internal_iv, IV_SIZE);
    assert(ret == 1);
    assert(len == (int)IV_SIZE);
    (void)ret;
    (void)len;
    memcpy(out, last_block, outlen % IV_SIZE);
  }
  EVP_CIPHER_CTX_free(ctx);
#elif defined(_WIN32)
  BCRYPT_ALG_HANDLE aes_handle = NULL;
  BCRYPT_KEY_HANDLE key_handle = NULL;

  NTSTATUS ret = BCryptOpenAlgorithmProvider(&aes_handle, BCRYPT_AES_ALGORITHM, NULL, 0);
  assert(BCRYPT_SUCCESS(ret));
  ret = BCryptGenerateSymmetricKey(aes_handle, &key_handle, NULL, 0, (PUCHAR)key, seclvl / 8, 0);
  assert(BCRYPT_SUCCESS(ret));

  ret = BCryptSetProperty(aes_handle, BCRYPT_CHAINING_MODE, (PUCHAR)BCRYPT_CHAIN_MODE_ECB,
                          (ULONG)((wcslen(BCRYPT_CHAIN_MODE_ECB) + 1) * sizeof(wchar_t)), 0);
  assert(BCRYPT_SUCCESS(ret));

  for (size_t idx = 0; idx < outlen / IV_SIZE; idx += 1, out += IV_SIZE) {
    ULONG len = 0;
    ret = BCryptEncrypt(key_handle, internal_iv, IV_SIZE, NULL, NULL, 0, out, IV_SIZE, &len, 0);
    assert(BCRYPT_SUCCESS(ret));
    assert(len == (ULONG)IV_SIZE);
    (void)ret;
    (void)len;
    aes_increment_iv(internal_iv);
  }
  if (outlen % IV_SIZE) {
    uint8_t last_block[IV_SIZE];
    ULONG len = 0;
    ret = BCryptEncrypt(key_handle, internal_iv, IV_SIZE, NULL, NULL, 0, last_block, IV_SIZE, &len,
                        0);
    assert(BCRYPT_SUCCESS(ret));
    assert(len == (ULONG)IV_SIZE);
    (void)ret;
    (void)len;
    memcpy(out, last_block, outlen % IV_SIZE);
  }

  BCryptDestroyKey(key_handle);
  BCryptCloseAlgorithmProvider(aes_handle, 0);
#else
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
#endif
}

void aes_extend_witness(uint8_t* w, const uint8_t* key, const uint8_t* in,
                        const faest_paramset_t* params) {
  const unsigned int lambda      = params->lambda;
  const unsigned int num_rounds  = params->R;
  const unsigned int blocksize   = 32 * params->Nst;
  const unsigned int beta        = (lambda + blocksize - 1) / blocksize;
  const unsigned int block_words = blocksize / 32;
  const bool is_em               = faest_is_em(params);

#if !defined(NDEBUG)
  uint8_t* const w_out = w;
#endif

  if (is_em) {
    // switch input and key for EM
    const uint8_t* tmp = key;
    key                = in;
    in                 = tmp;
  }

  // Step 3
  aes_round_keys_t round_keys;
  switch (lambda) {
  case 256:
    if (is_em) {
      // for scan-build
      assert(block_words == RIJNDAEL_BLOCK_WORDS_256);
      rijndael256_init_round_keys(&round_keys, key);
    } else {
      aes256_init_round_keys(&round_keys, key);
    }
    break;
  case 192:
    if (is_em) {
      // for scan-build
      assert(block_words == RIJNDAEL_BLOCK_WORDS_192);
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
  if (!is_em) {
    const unsigned int nk = lambda / 32;

    // Key schedule constraints only needed for normal AES, not EM variant.
    for (unsigned int i = 0; i != nk; ++i) {
      memcpy(w, round_keys.round_keys[i / 4][i % 4], sizeof(aes_word_t));
      w += sizeof(aes_word_t);
    }

    const unsigned int S_ke = params->Ske;
    for (unsigned int j = 0, ik = nk; j < S_ke / 4; ++j) {
      memcpy(w, round_keys.round_keys[ik / 4][ik % 4], sizeof(aes_word_t));
      w += sizeof(aes_word_t);
      ik += lambda == 192 ? 6 : 4;
    }
  } else {
    const unsigned int lambda_bytes = lambda / 8;

    // saving the OWF key to the extended witness
    memcpy(w, in, lambda_bytes);
    w += lambda_bytes;
  }

  assert(w - w_out == params->Lke / 8);

  // Step 10
  // common part for AES-128, EM-128, EM-192, EM-256, first part for AES-192 and AES-256
  {
    // Step 12
    aes_block_t state;
    load_state(state, in, block_words);

    // Step 13
    add_round_key(0, state, &round_keys, block_words);

    for (unsigned int round = 1; round < num_rounds; ++round) {
      if (round % 2 == 1) {
        // save inverse norm of the S-box inputs, in coloumn major order
        w = store_invnorm_state(w, state, block_words);
      }
      // Step 15
      sub_bytes(state, block_words);
      // Step 16
      shift_row(state, block_words);
      if (round % 2 == 0) {
        // Step 17
        w = store_state(w, state, block_words);
      }
      // Step 18
      mix_column(state, block_words);
      // Step 19
      add_round_key(round, state, &round_keys, block_words);
    }
    // last round is not commited to, so not computed
  }

  if (beta == 2) {
    // AES-192 and AES-256
    uint8_t buf[16];
    memcpy(buf, in, sizeof(buf));
    buf[0] ^= 0x1;

    // Step 12
    aes_block_t state;
    load_state(state, buf, block_words);

    // Step 13
    add_round_key(0, state, &round_keys, block_words);

    for (unsigned int round = 1; round < num_rounds; ++round) {
      if (round % 2 == 1) {
        // save inverse norm of the S-box inputs, in coloumn major order
        w = store_invnorm_state(w, state, block_words);
      }
      // Step 15
      sub_bytes(state, block_words);
      // Step 16
      shift_row(state, block_words);
      if (round % 2 == 0) {
        // Step 17
        w = store_state(w, state, block_words);
      }
      // Step 18
      mix_column(state, block_words);
      // Step 19
      add_round_key(round, state, &round_keys, block_words);
    }
    // last round is not commited to, so not computed
  }

  assert(w - w_out == params->l / 8);
}
