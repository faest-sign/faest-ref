/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_FIELDS_NEON_H
#define FAEST_FIELDS_NEON_H

// Included by fields.h after the common field types and sizes.
// The NEON backend assumes little-endian targets.
#if defined(_MSC_VER) && !defined(__clang__) && (defined(_M_ARM64) || defined(_M_ARM64EC))
#include <arm64_neon.h>
#else
#include <arm_neon.h>
#endif

// Form a register with a zero high limb, including padding in 192- and 576-bit values.
ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline uint64x2_t bf_neon_low64(uint64_t low) {
  return vcombine_u64(vcreate_u64(low), vdup_n_u64(0));
}

// Extract a polynomial limb without aliasing the register storage as uint64_t.
ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline uint64_t bf_word(uint64x2_t value,
                                                                             unsigned int lane) {
  return lane ? vgetq_lane_u64(value, 1) : vgetq_lane_u64(value, 0);
}

// Byte-vector accesses permit unaligned serialized input and output.
ATTR_PURE ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline uint64x2_t
bf_neon_load(const uint8_t* src) {
  return vreinterpretq_u64_u8(vld1q_u8(src));
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf_neon_store(uint8_t* dst,
                                                                    uint64x2_t value) {
  vst1q_u8(dst, vreinterpretq_u8_u64(value));
}

#define BF_VALUE(v, i) bf_word((v).values[(i) / 2], (i) % 2)

typedef struct {
  uint64x2_t values[1];
} bf128_t;

typedef struct {
  uint64x2_t values[2];
} bf192_t;

typedef struct {
  uint64x2_t values[2];
} bf256_t;

typedef struct {
  uint64x2_t values[3];
} bf384_t;

typedef struct {
  uint64x2_t values[5];
} bf576_t;

typedef struct {
  uint64x2_t values[6];
} bf768_t;

#define BF128_ALIGN 16
#define BF192_ALIGN 16
#define BF256_ALIGN 16

// GF(2^128) implementation

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf128_load(bf128_t* dst, const uint8_t* src) {
  dst->values[0] = bf_neon_load(src + 0);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf128_store(uint8_t* dst,
                                                                  const bf128_t* src) {
  bf_neon_store(dst + 0, src->values[0]);
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf128_t bf128_zero(void) {
  bf128_t ret;
  ret.values[0] = vdupq_n_u64(0);
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf128_t bf128_from_bit(uint8_t bit) {
  bf128_t ret;
  ret.values[0] = bf_neon_low64(bit & 1u);
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf128_t bf128_one(void) {
  bf128_t ret;
  ret.values[0] = bf_neon_low64(1);
  return ret;
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf128_add(bf128_t* dst, const bf128_t* lhs,
                                                                const bf128_t* rhs) {
  dst->values[0] = veorq_u64(lhs->values[0], rhs->values[0]);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf128_add_inplace(bf128_t* lhs,
                                                                        const bf128_t* rhs) {
  lhs->values[0] = veorq_u64(lhs->values[0], rhs->values[0]);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void
bf128_mul_bit(bf128_t* dst, const bf128_t* lhs, uint8_t rhs) {
  const uint64x2_t mask = vdupq_n_u64(UINT64_C(0) - (rhs & 1u));
  dst->values[0]        = vandq_u64(lhs->values[0], mask);
}

// Multiply lhs by each bit of rhs, from least to most significant.
void bf128_mul_bit_8(bf128_t* dst, const bf128_t* lhs, uint8_t rhs);
void bf128_mul(bf128_t* dst, const bf128_t* lhs, const bf128_t* rhs);
void bf128_square(bf128_t* dst, const bf128_t* lhs);
void bf128_mul_inplace(bf128_t* lhs, const bf128_t* rhs);
void bf128_mul_64_inplace(bf128_t* lhs, bf64_t rhs);

// GF(2^192) implementation

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf192_load(bf192_t* dst, const uint8_t* src) {
  dst->values[0] = bf_neon_load(src + 0);
  dst->values[1] = vcombine_u64(vreinterpret_u64_u8(vld1_u8(src + 16)), vdup_n_u64(0));
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf192_store(uint8_t* dst,
                                                                  const bf192_t* src) {
  bf_neon_store(dst + 0, src->values[0]);
  vst1_u8(dst + 16, vreinterpret_u8_u64(vget_low_u64(src->values[1])));
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf192_t bf192_zero(void) {
  bf192_t ret;
  for (unsigned int i = 0; i != 2; ++i) {
    ret.values[i] = vdupq_n_u64(0);
  }
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf192_t bf192_from_bit(uint8_t bit) {
  bf192_t ret   = bf192_zero();
  ret.values[0] = bf_neon_low64(bit & 1u);
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf192_t bf192_one(void) {
  return bf192_from_bit(1);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf192_add(bf192_t* dst, const bf192_t* lhs,
                                                                const bf192_t* rhs) {
  for (unsigned int i = 0; i != 2; ++i) {
    dst->values[i] = veorq_u64(lhs->values[i], rhs->values[i]);
  }
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf192_add_inplace(bf192_t* lhs,
                                                                        const bf192_t* rhs) {
  bf192_add(lhs, lhs, rhs);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void
bf192_mul_bit(bf192_t* dst, const bf192_t* lhs, uint8_t rhs) {
  const uint64x2_t mask = vdupq_n_u64(UINT64_C(0) - (rhs & 1u));
  for (unsigned int i = 0; i != 2; ++i) {
    dst->values[i] = vandq_u64(lhs->values[i], mask);
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
  dst->values[0] = bf_neon_load(src + 0);
  dst->values[1] = bf_neon_load(src + 16);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf256_store(uint8_t* dst,
                                                                  const bf256_t* src) {
  bf_neon_store(dst + 0, src->values[0]);
  bf_neon_store(dst + 16, src->values[1]);
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf256_t bf256_zero(void) {
  bf256_t ret;
  for (unsigned int i = 0; i != 2; ++i) {
    ret.values[i] = vdupq_n_u64(0);
  }
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf256_t bf256_from_bit(uint8_t bit) {
  bf256_t ret;
  ret.values[0] = bf_neon_low64(bit & 1u);
  ret.values[1] = vdupq_n_u64(0);
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf256_t bf256_one(void) {
  return bf256_from_bit(1);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf256_add(bf256_t* dst, const bf256_t* lhs,
                                                                const bf256_t* rhs) {
  for (unsigned int i = 0; i != 2; ++i) {
    dst->values[i] = veorq_u64(lhs->values[i], rhs->values[i]);
  }
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf256_add_inplace(bf256_t* lhs,
                                                                        const bf256_t* rhs) {
  bf256_add(lhs, lhs, rhs);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void
bf256_mul_bit(bf256_t* dst, const bf256_t* lhs, uint8_t rhs) {
  const uint64x2_t mask = vdupq_n_u64(UINT64_C(0) - (rhs & 1u));
  for (unsigned int i = 0; i != 2; ++i) {
    dst->values[i] = vandq_u64(lhs->values[i], mask);
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
  dst->values[0] = bf_neon_load(src + 0);
  dst->values[1] = bf_neon_load(src + 16);
  dst->values[2] = bf_neon_load(src + 32);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf384_store(uint8_t* dst,
                                                                  const bf384_t* src) {
  bf_neon_store(dst + 0, src->values[0]);
  bf_neon_store(dst + 16, src->values[1]);
  bf_neon_store(dst + 32, src->values[2]);
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf384_t bf384_zero(void) {
  bf384_t ret;
  for (unsigned int i = 0; i != 3; ++i) {
    ret.values[i] = vdupq_n_u64(0);
  }
  return ret;
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf384_add(bf384_t* dst, const bf384_t* lhs,
                                                                const bf384_t* rhs) {
  for (unsigned int i = 0; i != 3; ++i) {
    dst->values[i] = veorq_u64(lhs->values[i], rhs->values[i]);
  }
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf384_add_inplace(bf384_t* lhs,
                                                                        const bf384_t* rhs) {
  bf384_add(lhs, lhs, rhs);
}

// GF(2^576) implementation

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf576_load(bf576_t* dst, const uint8_t* src) {
  dst->values[0] = bf_neon_load(src + 0);
  dst->values[1] = bf_neon_load(src + 16);
  dst->values[2] = bf_neon_load(src + 32);
  dst->values[3] = bf_neon_load(src + 48);
  dst->values[4] = vcombine_u64(vreinterpret_u64_u8(vld1_u8(src + 64)), vdup_n_u64(0));
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf576_store(uint8_t* dst,
                                                                  const bf576_t* src) {
  bf_neon_store(dst + 0, src->values[0]);
  bf_neon_store(dst + 16, src->values[1]);
  bf_neon_store(dst + 32, src->values[2]);
  bf_neon_store(dst + 48, src->values[3]);
  vst1_u8(dst + 64, vreinterpret_u8_u64(vget_low_u64(src->values[4])));
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf576_t bf576_zero(void) {
  bf576_t ret;
  for (unsigned int i = 0; i != 5; ++i) {
    ret.values[i] = vdupq_n_u64(0);
  }
  return ret;
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf576_add(bf576_t* dst, const bf576_t* lhs,
                                                                const bf576_t* rhs) {
  for (unsigned int i = 0; i != 5; ++i) {
    dst->values[i] = veorq_u64(lhs->values[i], rhs->values[i]);
  }
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf576_add_inplace(bf576_t* lhs,
                                                                        const bf576_t* rhs) {
  bf576_add(lhs, lhs, rhs);
}

// GF(2^768) implementation

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf768_load(bf768_t* dst, const uint8_t* src) {
  dst->values[0] = bf_neon_load(src + 0);
  dst->values[1] = bf_neon_load(src + 16);
  dst->values[2] = bf_neon_load(src + 32);
  dst->values[3] = bf_neon_load(src + 48);
  dst->values[4] = bf_neon_load(src + 64);
  dst->values[5] = bf_neon_load(src + 80);
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf768_store(uint8_t* dst,
                                                                  const bf768_t* src) {
  bf_neon_store(dst + 0, src->values[0]);
  bf_neon_store(dst + 16, src->values[1]);
  bf_neon_store(dst + 32, src->values[2]);
  bf_neon_store(dst + 48, src->values[3]);
  bf_neon_store(dst + 64, src->values[4]);
  bf_neon_store(dst + 80, src->values[5]);
}

ATTR_CONST ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline bf768_t bf768_zero(void) {
  bf768_t ret;
  for (unsigned int i = 0; i != 6; ++i) {
    ret.values[i] = vdupq_n_u64(0);
  }
  return ret;
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf768_add(bf768_t* dst, const bf768_t* lhs,
                                                                const bf768_t* rhs) {
  for (unsigned int i = 0; i != 6; ++i) {
    dst->values[i] = veorq_u64(lhs->values[i], rhs->values[i]);
  }
}

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void bf768_add_inplace(bf768_t* lhs,
                                                                        const bf768_t* rhs) {
  for (unsigned int i = 0; i != 6; ++i) {
    lhs->values[i] = veorq_u64(lhs->values[i], rhs->values[i]);
  }
}

#endif
