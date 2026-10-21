/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_FIELDS_SSE2_H
#define FAEST_FIELDS_SSE2_H

// Included by fields.h after the common field types and sizes.
#include <emmintrin.h>

// Extract a polynomial limb without aliasing the register storage as uint64_t.
ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline uint64_t bf_word(__m128i value,
                                                                             unsigned int lane) {
  if (lane) {
    value = _mm_srli_si128(value, 8);
  }
  return (uint64_t)_mm_cvtsi128_si64(value);
}

#define BF_VALUE(v, i) bf_word((v).values[(i) / 2], (i) % 2)

typedef struct {
  __m128i values[1];
} bf128_t;

typedef struct {
  __m128i values[2];
} bf192_t;

typedef struct {
  __m128i values[2];
} bf256_t;

typedef struct {
  __m128i values[3];
} bf384_t;

typedef struct {
  __m128i values[5];
} bf576_t;

typedef struct {
  __m128i values[6];
} bf768_t;

#define BF128_ALIGN 16
#define BF192_ALIGN 16
#define BF256_ALIGN 16

// GF(2^128) implementation

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf128_load(bf128_t* dst, const uint8_t* src) {
  dst->values[0] = _mm_loadu_si128((const __m128i*)(src + 0));
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf128_store(uint8_t* dst,
                                                                  const bf128_t* src) {
  _mm_storeu_si128((__m128i*)(dst + 0), src->values[0]);
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf128_t bf128_zero(void) {
  bf128_t ret;
  ret.values[0] = _mm_setzero_si128();
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf128_t bf128_from_bit(uint8_t bit) {
  bf128_t ret;
  ret.values[0] = _mm_set_epi64x(0, bit & 1u);
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf128_t bf128_one(void) {
  bf128_t ret;
  ret.values[0] = _mm_set_epi64x(0, 1);
  return ret;
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf128_add(bf128_t* dst, const bf128_t* lhs,
                                                                const bf128_t* rhs) {
  dst->values[0] = _mm_xor_si128(lhs->values[0], rhs->values[0]);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf128_add_inplace(bf128_t* lhs,
                                                                        const bf128_t* rhs) {
  lhs->values[0] = _mm_xor_si128(lhs->values[0], rhs->values[0]);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void
bf128_mul_bit(bf128_t* dst, const bf128_t* lhs, uint8_t rhs) {
  const __m128i mask = _mm_set1_epi32(-(int)(rhs & 1u));
  dst->values[0]     = _mm_and_si128(lhs->values[0], mask);
}

// Multiply lhs by each bit of rhs, from least to most significant.
void bf128_mul_bit_8(bf128_t* dst, const bf128_t* lhs, uint8_t rhs);
void bf128_mul(bf128_t* dst, const bf128_t* lhs, const bf128_t* rhs);
void bf128_square(bf128_t* dst, const bf128_t* lhs);
void bf128_mul_inplace(bf128_t* lhs, const bf128_t* rhs);
void bf128_mul_64_inplace(bf128_t* lhs, bf64_t rhs);

// GF(2^192) implementation

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf192_load(bf192_t* dst, const uint8_t* src) {
  dst->values[0] = _mm_loadu_si128((const __m128i*)(src + 0));
  dst->values[1] = _mm_loadl_epi64((const __m128i*)(src + 16));
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf192_store(uint8_t* dst,
                                                                  const bf192_t* src) {
  _mm_storeu_si128((__m128i*)(dst + 0), src->values[0]);
  _mm_storel_epi64((__m128i*)(dst + 16), src->values[1]);
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf192_t bf192_zero(void) {
  bf192_t ret;
  for (unsigned int i = 0; i != 2; ++i) {
    ret.values[i] = _mm_setzero_si128();
  }
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf192_t bf192_from_bit(uint8_t bit) {
  bf192_t ret   = bf192_zero();
  ret.values[0] = _mm_set_epi64x(0, bit & 1u);
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf192_t bf192_one(void) {
  return bf192_from_bit(1);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf192_add(bf192_t* dst, const bf192_t* lhs,
                                                                const bf192_t* rhs) {
  for (unsigned int i = 0; i != 2; ++i) {
    dst->values[i] = _mm_xor_si128(lhs->values[i], rhs->values[i]);
  }
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf192_add_inplace(bf192_t* lhs,
                                                                        const bf192_t* rhs) {
  bf192_add(lhs, lhs, rhs);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void
bf192_mul_bit(bf192_t* dst, const bf192_t* lhs, uint8_t rhs) {
  const __m128i mask = _mm_set1_epi32(-(int)(rhs & 1u));
  for (unsigned int i = 0; i != 2; ++i) {
    dst->values[i] = _mm_and_si128(lhs->values[i], mask);
  }
}

// Multiply lhs by each bit of rhs, from least to most significant.
void bf192_mul_bit_8(bf192_t* dst, const bf192_t* lhs, uint8_t rhs);
void bf192_mul(bf192_t* dst, const bf192_t* lhs, const bf192_t* rhs);
void bf192_square(bf192_t* dst, const bf192_t* lhs);
void bf192_mul_inplace(bf192_t* lhs, const bf192_t* rhs);
void bf192_mul_64_inplace(bf192_t* lhs, bf64_t rhs);

// GF(2^256) implementation

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf256_load(bf256_t* dst, const uint8_t* src) {
  dst->values[0] = _mm_loadu_si128((const __m128i*)(src + 0));
  dst->values[1] = _mm_loadu_si128((const __m128i*)(src + 16));
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf256_store(uint8_t* dst,
                                                                  const bf256_t* src) {
  _mm_storeu_si128((__m128i*)(dst + 0), src->values[0]);
  _mm_storeu_si128((__m128i*)(dst + 16), src->values[1]);
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf256_t bf256_zero(void) {
  bf256_t ret;
  for (unsigned int i = 0; i != 2; ++i) {
    ret.values[i] = _mm_setzero_si128();
  }
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf256_t bf256_from_bit(uint8_t bit) {
  bf256_t ret;
  ret.values[0] = _mm_set_epi64x(0, bit & 1u);
  ret.values[1] = _mm_setzero_si128();
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf256_t bf256_one(void) {
  return bf256_from_bit(1);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf256_add(bf256_t* dst, const bf256_t* lhs,
                                                                const bf256_t* rhs) {
  for (unsigned int i = 0; i != 2; ++i) {
    dst->values[i] = _mm_xor_si128(lhs->values[i], rhs->values[i]);
  }
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf256_add_inplace(bf256_t* lhs,
                                                                        const bf256_t* rhs) {
  bf256_add(lhs, lhs, rhs);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void
bf256_mul_bit(bf256_t* dst, const bf256_t* lhs, uint8_t rhs) {
  const __m128i mask = _mm_set1_epi32(-(int)(rhs & 1u));
  for (unsigned int i = 0; i != 2; ++i) {
    dst->values[i] = _mm_and_si128(lhs->values[i], mask);
  }
}

// Multiply lhs by each bit of rhs, from least to most significant.
void bf256_mul_bit_8(bf256_t* dst, const bf256_t* lhs, uint8_t rhs);
void bf256_mul(bf256_t* dst, const bf256_t* lhs, const bf256_t* rhs);
void bf256_square(bf256_t* dst, const bf256_t* lhs);
void bf256_mul_inplace(bf256_t* lhs, const bf256_t* rhs);
void bf256_mul_64_inplace(bf256_t* lhs, bf64_t rhs);

// GF(2^384) implementation

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf384_load(bf384_t* dst, const uint8_t* src) {
  dst->values[0] = _mm_loadu_si128((const __m128i*)(src + 0));
  dst->values[1] = _mm_loadu_si128((const __m128i*)(src + 16));
  dst->values[2] = _mm_loadu_si128((const __m128i*)(src + 32));
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf384_store(uint8_t* dst,
                                                                  const bf384_t* src) {
  _mm_storeu_si128((__m128i*)(dst + 0), src->values[0]);
  _mm_storeu_si128((__m128i*)(dst + 16), src->values[1]);
  _mm_storeu_si128((__m128i*)(dst + 32), src->values[2]);
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf384_t bf384_zero(void) {
  bf384_t ret;
  for (unsigned int i = 0; i != 3; ++i) {
    ret.values[i] = _mm_setzero_si128();
  }
  return ret;
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf384_add(bf384_t* dst, const bf384_t* lhs,
                                                                const bf384_t* rhs) {
  for (unsigned int i = 0; i != 3; ++i) {
    dst->values[i] = _mm_xor_si128(lhs->values[i], rhs->values[i]);
  }
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf384_add_inplace(bf384_t* lhs,
                                                                        const bf384_t* rhs) {
  bf384_add(lhs, lhs, rhs);
}

// GF(2^576) implementation

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf576_load(bf576_t* dst, const uint8_t* src) {
  dst->values[0] = _mm_loadu_si128((const __m128i*)(src + 0));
  dst->values[1] = _mm_loadu_si128((const __m128i*)(src + 16));
  dst->values[2] = _mm_loadu_si128((const __m128i*)(src + 32));
  dst->values[3] = _mm_loadu_si128((const __m128i*)(src + 48));
  dst->values[4] = _mm_loadl_epi64((const __m128i*)(src + 64));
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf576_store(uint8_t* dst,
                                                                  const bf576_t* src) {
  _mm_storeu_si128((__m128i*)(dst + 0), src->values[0]);
  _mm_storeu_si128((__m128i*)(dst + 16), src->values[1]);
  _mm_storeu_si128((__m128i*)(dst + 32), src->values[2]);
  _mm_storeu_si128((__m128i*)(dst + 48), src->values[3]);
  _mm_storel_epi64((__m128i*)(dst + 64), src->values[4]);
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf576_t bf576_zero(void) {
  bf576_t ret;
  for (unsigned int i = 0; i != 5; ++i) {
    ret.values[i] = _mm_setzero_si128();
  }
  return ret;
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf576_add(bf576_t* dst, const bf576_t* lhs,
                                                                const bf576_t* rhs) {
  for (unsigned int i = 0; i != 5; ++i) {
    dst->values[i] = _mm_xor_si128(lhs->values[i], rhs->values[i]);
  }
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf576_add_inplace(bf576_t* lhs,
                                                                        const bf576_t* rhs) {
  bf576_add(lhs, lhs, rhs);
}

// GF(2^768) implementation

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf768_load(bf768_t* dst, const uint8_t* src) {
  dst->values[0] = _mm_loadu_si128((const __m128i*)(src + 0));
  dst->values[1] = _mm_loadu_si128((const __m128i*)(src + 16));
  dst->values[2] = _mm_loadu_si128((const __m128i*)(src + 32));
  dst->values[3] = _mm_loadu_si128((const __m128i*)(src + 48));
  dst->values[4] = _mm_loadu_si128((const __m128i*)(src + 64));
  dst->values[5] = _mm_loadu_si128((const __m128i*)(src + 80));
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf768_store(uint8_t* dst,
                                                                  const bf768_t* src) {
  _mm_storeu_si128((__m128i*)(dst + 0), src->values[0]);
  _mm_storeu_si128((__m128i*)(dst + 16), src->values[1]);
  _mm_storeu_si128((__m128i*)(dst + 32), src->values[2]);
  _mm_storeu_si128((__m128i*)(dst + 48), src->values[3]);
  _mm_storeu_si128((__m128i*)(dst + 64), src->values[4]);
  _mm_storeu_si128((__m128i*)(dst + 80), src->values[5]);
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf768_t bf768_zero(void) {
  bf768_t ret;
  for (unsigned int i = 0; i != 6; ++i) {
    ret.values[i] = _mm_setzero_si128();
  }
  return ret;
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf768_add(bf768_t* dst, const bf768_t* lhs,
                                                                const bf768_t* rhs) {
  for (unsigned int i = 0; i != 6; ++i) {
    dst->values[i] = _mm_xor_si128(lhs->values[i], rhs->values[i]);
  }
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf768_add_inplace(bf768_t* lhs,
                                                                        const bf768_t* rhs) {
  for (unsigned int i = 0; i != 6; ++i) {
    lhs->values[i] = _mm_xor_si128(lhs->values[i], rhs->values[i]);
  }
}

#endif
