/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_AESNI_H
#define FAEST_AESNI_H

#include "macros.h"
#include "compat.h"

#include <wmmintrin.h>
#if defined(HAVE_AVX2)
#include <immintrin.h>
#endif

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

#endif
