/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_FIELDS_H
#define FAEST_FIELDS_H

#include "macros.h"
#include "endian_compat.h"

#include <stdint.h>
#include <string.h>

FAEST_BEGIN_C_DECL

typedef uint8_t bf8_t;
typedef uint64_t bf64_t;

#if defined(HAVE_ATTR_VECTOR_SIZE)
#define BF_VALUE(v, i) ((v)[i])

typedef uint64_t bf128_t ATTR_VECTOR_SIZE(16);
typedef uint64_t bf192_t ATTR_VECTOR_SIZE(32);
typedef uint64_t bf256_t ATTR_VECTOR_SIZE(32);

#define BF128_ALIGN 16
#define BF192_ALIGN 32
#define BF256_ALIGN 32

#define BF128C(x0, x1)                                                                             \
  { x0, x1 }
#define BF192C(x0, x1, x2)                                                                         \
  { x0, x1, x2, UINT64_C(0) }
#define BF256C(x0, x1, x2, x3)                                                                     \
  { x0, x1, x2, x3 }
#else

#define BF_VALUE(v, i) ((v).values[i])

typedef struct {
  uint64_t values[2];
} bf128_t;

typedef struct {
  uint64_t values[3];
} bf192_t;

typedef struct {
  uint64_t values[4];
} bf256_t;

#define BF128C(x0, x1)                                                                             \
  {                                                                                                \
    { x0, x1 }                                                                                     \
  }
#define BF192C(x0, x1, x2)                                                                         \
  {                                                                                                \
    { x0, x1, x2 }                                                                                 \
  }
#define BF256C(x0, x1, x2, x3)                                                                     \
  {                                                                                                \
    { x0, x1, x2, x3 }                                                                             \
  }

#define BF128_ALIGN 16
#define BF192_ALIGN 16
#define BF256_ALIGN 32
#endif

// Placed after bitfield types
#include "vbb.h"


#define BF128_NUM_BYTES (128 / 8)
#define BF192_NUM_BYTES (192 / 8)
#define BF256_NUM_BYTES (256 / 8)

// GF(2^8) implementation

ATTR_PURE ATTR_ALWAYS_INLINE static inline bf8_t bf8_load(const uint8_t* src) {
  return *src;
}

ATTR_ALWAYS_INLINE static inline void bf8_store(uint8_t* dst, bf8_t src) {
  *dst = src;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf8_t bf8_zero(void) {
  return 0;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf8_t bf8_one(void) {
  return 1;
}

bf8_t bf8_rand(void);

ATTR_CONST ATTR_ALWAYS_INLINE inline bf8_t bf8_add(bf8_t lhs, bf8_t rhs) {
  return lhs ^ rhs;
}

ATTR_CONST bf8_t bf8_mul(bf8_t lhs, bf8_t rhs);
ATTR_CONST bf8_t bf8_inv(bf8_t lhs);

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf8_t bf8_from_bit(uint8_t bit) {
  return bit & 1;
}

// GF(2^64) implementation

ATTR_PURE ATTR_ALWAYS_INLINE inline bf64_t bf64_load(const uint8_t* src) {
  bf64_t ret;
  memcpy(&ret, src, sizeof(ret));
#if defined(FAEST_IS_BIG_ENDIAN)
  ret = le64toh(ret);
#endif
  return ret;
}

ATTR_ALWAYS_INLINE inline void bf64_store(uint8_t* dst, bf64_t src) {
#if defined(FAEST_IS_BIG_ENDIAN)
  src = htole64(src);
#endif
  memcpy(dst, &src, sizeof(src));
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_zero(void) {
  return 0;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_one(void) {
  return 1;
}

bf64_t bf64_rand(void);

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_add(bf64_t lhs, bf64_t rhs) {
  return lhs ^ rhs;
}

ATTR_CONST bf64_t bf64_mul(bf64_t lhs, bf64_t rhs);

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_from_bit(uint8_t bit) {
  return bit & 1;
}

// GF(2^128) implementation

ATTR_PURE ATTR_ALWAYS_INLINE static inline bf128_t bf128_load(const uint8_t* src) {
  bf128_t ret;
#if defined(FAEST_IS_BIG_ENDIAN)
  for (unsigned int i = 0; i != BF128_NUM_BYTES / sizeof(uint64_t); ++i, src += sizeof(uint64_t)) {
    memcpy(&BF_VALUE(ret, i), src, sizeof(uint64_t));
    BF_VALUE(ret, i) = le64toh(BF_VALUE(ret, i));
  }
#else
  memcpy(&ret, src, BF128_NUM_BYTES);
#endif
  return ret;
}

ATTR_ALWAYS_INLINE static inline void bf128_store(uint8_t* dst, bf128_t src) {
#if defined(FAEST_IS_BIG_ENDIAN)
  for (unsigned int i = 0; i != BF128_NUM_BYTES / sizeof(uint64_t); ++i, dst += sizeof(uint64_t)) {
    uint64_t tmp = htole64(BF_VALUE(src, i));
    memcpy(dst, &tmp, sizeof(tmp));
  }
#else
  memcpy(dst, &src, BF128_NUM_BYTES);
#endif
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_from_bf64(bf64_t src) {
  bf128_t ret      = BF128C(0, 0);
  BF_VALUE(ret, 0) = src;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_from_bf8(bf8_t src) {
  bf128_t ret      = BF128C(0, 0);
  BF_VALUE(ret, 0) = src;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_from_bit(uint8_t bit) {
  return bf128_from_bf8(bit & 1);
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_zero(void) {
  const bf128_t ret = BF128C(0, 0);
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_one(void) {
  const bf128_t ret = BF128C(1, 0);
  return ret;
}

ATTR_PURE bf128_t bf128_byte_combine(const bf128_t* x);
ATTR_PURE bf128_t bf128_byte_combine_bits(uint8_t x);
bf128_t bf128_rand(void);

#if defined(HAVE_ATTR_VECTOR_SIZE)
#define bf128_add(lhs, rhs) ((lhs) ^ (rhs))
#else
ATTR_CONST static inline bf128_t bf128_add(bf128_t lhs, bf128_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}
#endif

ATTR_CONST bf128_t bf128_mul(bf128_t lhs, bf128_t rhs);
ATTR_CONST bf128_t bf128_mul_64(bf128_t lhs, bf64_t rhs);
#if defined(HAVE_ATTR_VECTOR_SIZE)
#define bf128_mul_bit(lhs, rhs) ((lhs) & -((uint64_t)(rhs)&1))
#else
ATTR_CONST bf128_t bf128_mul_bit(bf128_t lhs, uint8_t rhs);
#endif
ATTR_PURE bf128_t bf128_sum_poly(const bf128_t* xs);

// GF(2^192) implemenation

ATTR_PURE ATTR_ALWAYS_INLINE static inline bf192_t bf192_load(const uint8_t* src) {
  bf192_t ret;
#if defined(FAEST_IS_BIG_ENDIAN)
  for (unsigned int i = 0; i != BF192_NUM_BYTES / sizeof(uint64_t); ++i, src += sizeof(uint64_t)) {
    memcpy(&BF_VALUE(ret, i), src, sizeof(uint64_t));
    BF_VALUE(ret, i) = le64toh(BF_VALUE(ret, i));
  }
#else
  memcpy(&ret, src, BF192_NUM_BYTES);
#endif
#if defined(HAVE_ATTR_VECTOR_SIZE)
  BF_VALUE(ret, 3) = 0;
#endif
  return ret;
}

ATTR_ALWAYS_INLINE static inline void bf192_store(uint8_t* dst, bf192_t src) {
#if defined(FAEST_IS_BIG_ENDIAN)
  for (unsigned int i = 0; i != BF192_NUM_BYTES / sizeof(uint64_t); ++i, dst += sizeof(uint64_t)) {
    uint64_t tmp = htole64(BF_VALUE(src, i));
    memcpy(dst, &tmp, sizeof(tmp));
  }
#else
  memcpy(dst, &src, BF192_NUM_BYTES);
#endif
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_from_bf64(bf64_t src) {
  bf192_t ret      = BF192C(0, 0, 0);
  BF_VALUE(ret, 0) = src;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_from_bf8(bf8_t src) {
  bf192_t ret      = BF192C(0, 0, 0);
  BF_VALUE(ret, 0) = src;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_from_bit(uint8_t bit) {
  return bf192_from_bf8(bit & 1);
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_zero(void) {
  const bf192_t ret = BF192C(0, 0, 0);
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_one(void) {
  const bf192_t ret = BF192C(1, 0, 0);
  return ret;
}

ATTR_PURE bf192_t bf192_byte_combine(const bf192_t* x);
ATTR_PURE bf192_t bf192_byte_combine_bits(uint8_t x);
bf192_t bf192_rand(void);

#if defined(HAVE_ATTR_VECTOR_SIZE)
#define bf192_add(lhs, rhs) ((lhs) ^ (rhs))
#else
ATTR_CONST static inline bf192_t bf192_add(bf192_t lhs, bf192_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}
#endif

ATTR_CONST bf192_t bf192_mul(bf192_t lhs, bf192_t rhs);
ATTR_CONST bf192_t bf192_mul_64(bf192_t lhs, bf64_t rhs);
#if defined(HAVE_ATTR_VECTOR_SIZE)
#define bf192_mul_bit(lhs, rhs) ((lhs) & -((uint64_t)(rhs)&1))
#else
ATTR_CONST bf192_t bf192_mul_bit(bf192_t lhs, uint8_t rhs);
#endif
ATTR_PURE bf192_t bf192_sum_poly(const bf192_t* xs);

// GF(2^256) implementation

ATTR_PURE ATTR_ALWAYS_INLINE static inline bf256_t bf256_load(const uint8_t* src) {
  bf256_t ret;
#if defined(FAEST_IS_BIG_ENDIAN)
  for (unsigned int i = 0; i != BF256_NUM_BYTES / sizeof(uint64_t); ++i, src += sizeof(uint64_t)) {
    memcpy(&BF_VALUE(ret, i), src, sizeof(uint64_t));
    BF_VALUE(ret, i) = le64toh(BF_VALUE(ret, i));
  }
#else
  memcpy(&ret, src, BF256_NUM_BYTES);
#endif
  return ret;
}

ATTR_ALWAYS_INLINE static inline void bf256_store(uint8_t* dst, bf256_t src) {
#if defined(FAEST_IS_BIG_ENDIAN)
  for (unsigned int i = 0; i != BF256_NUM_BYTES / sizeof(uint64_t); ++i, dst += sizeof(uint64_t)) {
    uint64_t tmp = htole64(BF_VALUE(src, i));
    memcpy(dst, &tmp, sizeof(tmp));
  }
#else
  memcpy(dst, &src, BF256_NUM_BYTES);
#endif
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_from_bf64(bf64_t src) {
  bf256_t ret      = BF256C(0, 0, 0, 0);
  BF_VALUE(ret, 0) = src;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_from_bf8(bf8_t src) {
  bf256_t ret      = BF256C(0, 0, 0, 0);
  BF_VALUE(ret, 0) = src;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_from_bit(uint8_t bit) {
  return bf256_from_bf8(bit & 1);
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_zero(void) {
  const bf256_t ret = BF256C(0, 0, 0, 0);
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_one(void) {
  const bf256_t ret = BF256C(1, 0, 0, 0);
  return ret;
}

ATTR_PURE bf256_t bf256_byte_combine(const bf256_t* x);
ATTR_PURE bf256_t bf256_byte_combine_bits(uint8_t x);
bf256_t bf256_rand(void);

#if defined(HAVE_ATTR_VECTOR_SIZE)
#define bf256_add(lhs, rhs) ((lhs) ^ (rhs))
#else
ATTR_CONST static inline bf256_t bf256_add(bf256_t lhs, bf256_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}
#endif

ATTR_CONST bf256_t bf256_mul(bf256_t lhs, bf256_t rhs);
ATTR_CONST bf256_t bf256_mul_64(bf256_t lhs, bf64_t rhs);
#if defined(HAVE_ATTR_VECTOR_SIZE)
#define bf256_mul_bit(lhs, rhs) ((lhs) & -((uint64_t)(rhs)&1))
#else
ATTR_CONST bf256_t bf256_mul_bit(bf256_t lhs, uint8_t rhs);
#endif
ATTR_PURE bf256_t bf256_sum_poly(const bf256_t* xs);
bf256_t bf256_sum_poly_vbb(vbb_t* vbb, unsigned int offset);
bf256_t bf256_byte_combine_vbb(vbb_t* vbb, unsigned int offset);
bf192_t bf192_sum_poly_vbb(vbb_t* vbb, unsigned int offset);
bf192_t bf192_byte_combine_vbb(vbb_t* vbb, unsigned int offset);
bf128_t bf128_sum_poly_vbb(vbb_t* vbb, unsigned int offset);
bf128_t bf128_byte_combine_vbb(vbb_t* vbb, unsigned int offset);

bf128_t bf128_byte_combine_vk(vbb_t* vbb, unsigned int offset);
bf192_t bf192_byte_combine_vk(vbb_t* vbb, unsigned int offset);
bf256_t bf256_byte_combine_vk(vbb_t* vbb, unsigned int offset);
FAEST_END_C_DECL

#endif
