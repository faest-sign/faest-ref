/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "fields.h"

#if FAEST_FIELDS == FAEST_FIELDS_USE_NEON
// Insert a zero bit between the low 32 bits in each 64-bit lane.
// Numeric shifts keep polynomial bit ordering independent of target endianness.
static inline uint64x2_t bf_square_spread(uint64x2_t value) {
  value = vandq_u64(value, vdupq_n_u64(UINT64_C(0x00000000ffffffff)));
  value = vandq_u64(vorrq_u64(value, vshlq_n_u64(value, 16)),
                    vdupq_n_u64(UINT64_C(0x0000ffff0000ffff)));
  value =
      vandq_u64(vorrq_u64(value, vshlq_n_u64(value, 8)), vdupq_n_u64(UINT64_C(0x00ff00ff00ff00ff)));
  value =
      vandq_u64(vorrq_u64(value, vshlq_n_u64(value, 4)), vdupq_n_u64(UINT64_C(0x0f0f0f0f0f0f0f0f)));
  value =
      vandq_u64(vorrq_u64(value, vshlq_n_u64(value, 2)), vdupq_n_u64(UINT64_C(0x3333333333333333)));
  return vandq_u64(vorrq_u64(value, vshlq_n_u64(value, 1)),
                   vdupq_n_u64(UINT64_C(0x5555555555555555)));
}

// Expand each 64-bit polynomial limb into one 128-bit register.
static inline bf256_t bf_square_expand(uint64x2_t value) {
  const uint64x2_t upper = vshrq_n_u64(value, 32);
  bf256_t square;
  square.values[0] = bf_square_spread(vcombine_u64(vget_low_u64(value), vget_low_u64(upper)));
  square.values[1] = bf_square_spread(vcombine_u64(vget_high_u64(value), vget_high_u64(upper)));
  return square;
}

static inline uint64x2_t bf_square_fold_7_2_1(uint64x2_t high) {
  return veorq_u64(veorq_u64(high, vshlq_n_u64(high, 1)),
                   veorq_u64(vshlq_n_u64(high, 2), vshlq_n_u64(high, 7)));
}

static inline uint64x2_t bf_square_carry_7_2_1(uint64x2_t high) {
  return veorq_u64(veorq_u64(vshrq_n_u64(high, 63), vshrq_n_u64(high, 62)), vshrq_n_u64(high, 57));
}

static inline uint64x2_t bf_square_fold_10_5_2(uint64x2_t high) {
  return veorq_u64(veorq_u64(high, vshlq_n_u64(high, 2)),
                   veorq_u64(vshlq_n_u64(high, 5), vshlq_n_u64(high, 10)));
}

static inline uint64x2_t bf_square_carry_10_5_2(uint64x2_t high) {
  return veorq_u64(veorq_u64(vshrq_n_u64(high, 62), vshrq_n_u64(high, 59)), vshrq_n_u64(high, 54));
}

// Fold pairs of upper limbs and their carries into the lower half. The carry from the
// highest limb crosses the field boundary, so fold it once more into the lowest limb.
static inline void bf128_reduce(bf128_t* dst, uint64x2_t low, uint64x2_t high) {
  const uint64x2_t carry = bf_square_carry_7_2_1(high);

  low            = veorq_u64(low, bf_square_fold_7_2_1(high));
  low            = veorq_u64(low, vextq_u64(vdupq_n_u64(0), carry, 1));
  dst->values[0] = veorq_u64(low, bf_square_fold_7_2_1(vextq_u64(carry, vdupq_n_u64(0), 1)));
}

static inline void bf192_reduce(bf192_t* dst, uint64x2_t square0, uint64x2_t square1,
                                uint64x2_t square2) {
  const uint64x2_t high0  = vextq_u64(square1, square2, 1);
  const uint64x2_t high1  = vextq_u64(square2, vdupq_n_u64(0), 1);
  const uint64x2_t carry0 = bf_square_carry_7_2_1(high0);
  const uint64x2_t carry1 = bf_square_carry_7_2_1(high1);

  uint64x2_t low0 = veorq_u64(square0, bf_square_fold_7_2_1(high0));
  uint64x2_t low1 =
      veorq_u64(vcombine_u64(vget_low_u64(square1), vdup_n_u64(0)), bf_square_fold_7_2_1(high1));
  low0 = veorq_u64(low0, vextq_u64(vdupq_n_u64(0), carry0, 1));
  low0 = veorq_u64(low0, bf_square_fold_7_2_1(carry1));
  low1 = veorq_u64(low1, vextq_u64(carry0, vdupq_n_u64(0), 1));

  dst->values[0] = low0;
  dst->values[1] = vcombine_u64(vget_low_u64(low1), vdup_n_u64(0));
}

static inline void bf256_reduce(bf256_t* dst, uint64x2_t low0, uint64x2_t low1, uint64x2_t high0,
                                uint64x2_t high1) {
  const uint64x2_t carry0 = bf_square_carry_10_5_2(high0);
  const uint64x2_t carry1 = bf_square_carry_10_5_2(high1);

  low0 = veorq_u64(low0, bf_square_fold_10_5_2(high0));
  low1 = veorq_u64(low1, bf_square_fold_10_5_2(high1));
  low0 = veorq_u64(low0, vextq_u64(vdupq_n_u64(0), carry0, 1));
  low0 = veorq_u64(low0, bf_square_fold_10_5_2(vextq_u64(carry1, vdupq_n_u64(0), 1)));
  low1 = veorq_u64(low1, vextq_u64(vdupq_n_u64(0), carry1, 1));
  low1 = veorq_u64(low1, vextq_u64(carry0, vdupq_n_u64(0), 1));

  dst->values[0] = low0;
  dst->values[1] = low1;
}

static inline uint64x2_t bf_shift_left_1(uint64x2_t value) {
  return vorrq_u64(vshlq_n_u64(value, 1), vextq_u64(vdupq_n_u64(0), vshrq_n_u64(value, 63), 1));
}

static inline uint64x2_t bf_shift_left_1_high(uint64x2_t high, uint64x2_t low) {
  // Combine the preceding limb for each lane before extracting its carry bit.
  const uint64x2_t preceding = vextq_u64(low, high, 1);
  return vorrq_u64(vshlq_n_u64(high, 1), vshrq_n_u64(preceding, 63));
}

// Double a polynomial and fold the overflow into its low register.
static void bf_dbl(uint64x2_t* value, unsigned int words, uint64_t modulus) {
  const unsigned int count = (words + 1) / 2;
  const uint64_t overflow  = bf_word(value[count - 1], (words - 1) % 2) >> 63;
  for (unsigned int i = count - 1; i > 0; --i) {
    value[i] = bf_shift_left_1_high(value[i], value[i - 1]);
  }
  value[0] = veorq_u64(bf_shift_left_1(value[0]), bf_neon_low64(modulus & -overflow));
  if (words & 1) {
    value[count - 1] = vcombine_u64(vget_low_u64(value[count - 1]), vdup_n_u64(0));
  }
}

// GF(2^128) implementation

// Process four multiplier bits at a time using four fixed multiples of the input.
// Selection uses masks rather than a table indexed by multiplier bits.
// Both sizes are constant at each call site, allowing the compiler to specialize the loops.
void bf128_mul(bf128_t* dst, const bf128_t* lhs, const bf128_t* rhs) {
  uint64x2_t input_low[4];
  input_low[0] = lhs->values[0];

  const uint64x2_t mod = bf_neon_low64(BF128_MODULUS);
  for (unsigned int i = 1; i < 4; ++i) {
    const uint64x2_t top      = input_low[i - 1];
    const uint64x2_t overflow = vdupq_n_u64(UINT64_C(0) - (bf_word(top, 1) >> 63));
    input_low[i] = veorq_u64(bf_shift_left_1(input_low[i - 1]), vandq_u64(overflow, mod));
  }

  uint64x2_t result_low = vdupq_n_u64(0);
  for (unsigned int word = 2; word-- > 0;) {
    uint64_t multiplier = bf_word(rhs->values[0], word % 2);
    for (unsigned int nibble = 0; nibble < 16; ++nibble, multiplier <<= 4) {
      uint64x2_t overflow = vextq_u64(vshrq_n_u64(result_low, 60), vdupq_n_u64(0), 1);

      result_low = vorrq_u64(vshlq_n_u64(result_low, 4),
                             vextq_u64(vdupq_n_u64(0), vshrq_n_u64(result_low, 60), 1));
      // Only four bits overflow, so folding them cannot cross another limb boundary.
      result_low = veorq_u64(result_low, bf_square_fold_7_2_1(overflow));
      for (unsigned int bit = 0; bit < 4; ++bit) {
        const uint64x2_t mask = vdupq_n_u64(UINT64_C(0) - ((multiplier >> (60 + bit)) & 1));
        result_low            = veorq_u64(result_low, vandq_u64(input_low[bit], mask));
      }
    }
  }
  dst->values[0] = result_low;
}

void bf128_mul_inplace(bf128_t* lhs, const bf128_t* rhs) {
  bf128_mul(lhs, lhs, rhs);
}

void bf128_mul_64_inplace(bf128_t* lhs, bf64_t rhs) {
  uint64x2_t input_low[4];
  input_low[0] = lhs->values[0];

  const uint64x2_t mod = bf_neon_low64(BF128_MODULUS);
  for (unsigned int i = 1; i < 4; ++i) {
    const uint64x2_t top      = input_low[i - 1];
    const uint64x2_t overflow = vdupq_n_u64(UINT64_C(0) - (bf_word(top, 1) >> 63));
    input_low[i] = veorq_u64(bf_shift_left_1(input_low[i - 1]), vandq_u64(overflow, mod));
  }

  uint64x2_t result_low = vdupq_n_u64(0);
  for (unsigned int nibble = 0; nibble < 16; ++nibble, rhs <<= 4) {
    uint64x2_t overflow = vextq_u64(vshrq_n_u64(result_low, 60), vdupq_n_u64(0), 1);

    result_low = vorrq_u64(vshlq_n_u64(result_low, 4),
                           vextq_u64(vdupq_n_u64(0), vshrq_n_u64(result_low, 60), 1));
    // Only four bits overflow, so folding them cannot cross another limb boundary.
    result_low = veorq_u64(result_low, bf_square_fold_7_2_1(overflow));
    for (unsigned int bit = 0; bit < 4; ++bit) {
      const uint64x2_t mask = vdupq_n_u64(UINT64_C(0) - ((rhs >> (60 + bit)) & 1));
      result_low            = veorq_u64(result_low, vandq_u64(input_low[bit], mask));
    }
  }
  lhs->values[0] = result_low;
}

void bf128_square(bf128_t* dst, const bf128_t* lhs) {
  const bf256_t square = bf_square_expand(lhs->values[0]);
  bf128_reduce(dst, square.values[0], square.values[1]);
}

#if defined(FAEST_TESTS)
void bf128_sum_poly(bf128_t* dst, const bf128_t* xs) {
  *dst = xs[128 - 1];
  for (unsigned int i = 1; i != 128; ++i) {
    bf_dbl(dst->values, 2, BF128_MODULUS);
    bf128_add_inplace(dst, &xs[128 - 1 - i]);
  }
}
#endif

void bf128_sq_bit(bf128_t* out_tag, const bf128_t* in_tag) {
  bf128_add(&out_tag[7], &in_tag[6], &in_tag[7]);
  bf128_add(&out_tag[0], &in_tag[0], &in_tag[4]);
  bf128_add_inplace(&out_tag[0], &in_tag[6]);
  bf128_add(&out_tag[1], &in_tag[4], &out_tag[7]);
  bf128_add(&out_tag[2], &in_tag[1], &in_tag[5]);
  bf128_add(&out_tag[3], &out_tag[1], &in_tag[5]);
  bf128_add(&out_tag[4], &in_tag[2], &in_tag[4]);
  bf128_add_inplace(&out_tag[4], &in_tag[7]);
  bf128_add(&out_tag[5], &in_tag[5], &in_tag[6]);
  bf128_add(&out_tag[6], &in_tag[3], &in_tag[5]);
}

// GF(2^192) implementation

void bf192_mul(bf192_t* dst, const bf192_t* lhs, const bf192_t* rhs) {
  uint64x2_t input_low[4];
  uint64x2_t input_high[4];
  input_low[0]  = lhs->values[0];
  input_high[0] = lhs->values[1];

  const uint64x2_t mod = bf_neon_low64(BF192_MODULUS);
  for (unsigned int i = 1; i < 4; ++i) {
    const uint64x2_t top      = input_high[i - 1];
    const uint64x2_t overflow = vdupq_n_u64(UINT64_C(0) - (bf_word(top, 0) >> 63));
    input_low[i]  = veorq_u64(bf_shift_left_1(input_low[i - 1]), vandq_u64(overflow, mod));
    input_high[i] = vorrq_u64(vshlq_n_u64(input_high[i - 1], 1),
                              vextq_u64(vshrq_n_u64(input_low[i - 1], 63), vdupq_n_u64(0), 1));
  }

  uint64x2_t result_low  = vdupq_n_u64(0);
  uint64x2_t result_high = vdupq_n_u64(0);
  for (unsigned int word = 3; word-- > 0;) {
    uint64_t multiplier = bf_word(rhs->values[word / 2], word % 2);
    for (unsigned int nibble = 0; nibble < 16; ++nibble, multiplier <<= 4) {
      uint64x2_t overflow = vshrq_n_u64(result_high, 60);
      result_high         = vorrq_u64(vshlq_n_u64(result_high, 4),
                                      vextq_u64(vshrq_n_u64(result_low, 60), vdupq_n_u64(0), 1));

      result_low = vorrq_u64(vshlq_n_u64(result_low, 4),
                             vextq_u64(vdupq_n_u64(0), vshrq_n_u64(result_low, 60), 1));
      // Only four bits overflow, so folding them cannot cross another limb boundary.
      result_low = veorq_u64(result_low, bf_square_fold_7_2_1(overflow));
      for (unsigned int bit = 0; bit < 4; ++bit) {
        const uint64x2_t mask = vdupq_n_u64(UINT64_C(0) - ((multiplier >> (60 + bit)) & 1));
        result_low            = veorq_u64(result_low, vandq_u64(input_low[bit], mask));
        result_high           = veorq_u64(result_high, vandq_u64(input_high[bit], mask));
      }
    }
  }
  dst->values[0] = result_low;
  dst->values[1] = result_high;
}

void bf192_mul_inplace(bf192_t* lhs, const bf192_t* rhs) {
  bf192_mul(lhs, lhs, rhs);
}

void bf192_mul_64_inplace(bf192_t* lhs, bf64_t rhs) {
  uint64x2_t input_low[4];
  uint64x2_t input_high[4];
  input_low[0]  = lhs->values[0];
  input_high[0] = lhs->values[1];

  const uint64x2_t mod = bf_neon_low64(BF192_MODULUS);
  for (unsigned int i = 1; i < 4; ++i) {
    const uint64x2_t top      = input_high[i - 1];
    const uint64x2_t overflow = vdupq_n_u64(UINT64_C(0) - (bf_word(top, 0) >> 63));
    input_low[i]  = veorq_u64(bf_shift_left_1(input_low[i - 1]), vandq_u64(overflow, mod));
    input_high[i] = vorrq_u64(vshlq_n_u64(input_high[i - 1], 1),
                              vextq_u64(vshrq_n_u64(input_low[i - 1], 63), vdupq_n_u64(0), 1));
  }

  uint64x2_t result_low  = vdupq_n_u64(0);
  uint64x2_t result_high = vdupq_n_u64(0);
  for (unsigned int nibble = 0; nibble < 16; ++nibble, rhs <<= 4) {
    uint64x2_t overflow = vshrq_n_u64(result_high, 60);
    result_high         = vorrq_u64(vshlq_n_u64(result_high, 4),
                                    vextq_u64(vshrq_n_u64(result_low, 60), vdupq_n_u64(0), 1));

    result_low = vorrq_u64(vshlq_n_u64(result_low, 4),
                           vextq_u64(vdupq_n_u64(0), vshrq_n_u64(result_low, 60), 1));
    // Only four bits overflow, so folding them cannot cross another limb boundary.
    result_low = veorq_u64(result_low, bf_square_fold_7_2_1(overflow));
    for (unsigned int bit = 0; bit < 4; ++bit) {
      const uint64x2_t mask = vdupq_n_u64(UINT64_C(0) - ((rhs >> (60 + bit)) & 1));
      result_low            = veorq_u64(result_low, vandq_u64(input_low[bit], mask));
      result_high           = veorq_u64(result_high, vandq_u64(input_high[bit], mask));
    }
  }
  lhs->values[0] = result_low;
  lhs->values[1] = result_high;
}

void bf192_square(bf192_t* dst, const bf192_t* lhs) {
  const bf256_t square0 = bf_square_expand(lhs->values[0]);
  const bf256_t square1 = bf_square_expand(lhs->values[1]);
  bf192_reduce(dst, square0.values[0], square0.values[1], square1.values[0]);
}

#if defined(FAEST_TESTS)
void bf192_sum_poly(bf192_t* dst, const bf192_t* xs) {
  *dst = xs[192 - 1];
  for (unsigned int i = 1; i != 192; ++i) {
    bf_dbl(dst->values, 3, BF192_MODULUS);
    bf192_add_inplace(dst, &xs[192 - 1 - i]);
  }
}
#endif

void bf192_sq_bit(bf192_t* out_tag, const bf192_t* in_tag) {
  bf192_add(&out_tag[7], &in_tag[6], &in_tag[7]);
  bf192_add(&out_tag[0], &in_tag[0], &in_tag[4]);
  bf192_add_inplace(&out_tag[0], &in_tag[6]);
  bf192_add(&out_tag[1], &in_tag[4], &out_tag[7]);
  bf192_add(&out_tag[2], &in_tag[1], &in_tag[5]);
  bf192_add(&out_tag[3], &out_tag[1], &in_tag[5]);
  bf192_add(&out_tag[4], &in_tag[2], &in_tag[4]);
  bf192_add_inplace(&out_tag[4], &in_tag[7]);
  bf192_add(&out_tag[5], &in_tag[5], &in_tag[6]);
  bf192_add(&out_tag[6], &in_tag[3], &in_tag[5]);
}

// GF(2^256) implementation

void bf256_mul(bf256_t* dst, const bf256_t* lhs, const bf256_t* rhs) {
  uint64x2_t input_low[4];
  uint64x2_t input_high[4];
  input_low[0]  = lhs->values[0];
  input_high[0] = lhs->values[1];

  const uint64x2_t mod = bf_neon_low64(BF256_MODULUS);
  for (unsigned int i = 1; i < 4; ++i) {
    const uint64x2_t top      = input_high[i - 1];
    const uint64x2_t overflow = vdupq_n_u64(UINT64_C(0) - (bf_word(top, 1) >> 63));
    input_low[i]  = veorq_u64(bf_shift_left_1(input_low[i - 1]), vandq_u64(overflow, mod));
    input_high[i] = bf_shift_left_1_high(input_high[i - 1], input_low[i - 1]);
  }

  uint64x2_t result_low  = vdupq_n_u64(0);
  uint64x2_t result_high = vdupq_n_u64(0);
  for (unsigned int word = 4; word-- > 0;) {
    uint64_t multiplier = bf_word(rhs->values[word / 2], word % 2);
    for (unsigned int nibble = 0; nibble < 16; ++nibble, multiplier <<= 4) {
      uint64x2_t overflow        = vextq_u64(vshrq_n_u64(result_high, 60), vdupq_n_u64(0), 1);
      const uint64x2_t preceding = vextq_u64(result_low, result_high, 1);
      result_high = vorrq_u64(vshlq_n_u64(result_high, 4), vshrq_n_u64(preceding, 60));
      result_low  = vorrq_u64(vshlq_n_u64(result_low, 4),
                              vextq_u64(vdupq_n_u64(0), vshrq_n_u64(result_low, 60), 1));
      // Only four bits overflow, so folding them cannot cross another limb boundary.
      result_low = veorq_u64(result_low, bf_square_fold_10_5_2(overflow));
      for (unsigned int bit = 0; bit < 4; ++bit) {
        const uint64x2_t mask = vdupq_n_u64(UINT64_C(0) - ((multiplier >> (60 + bit)) & 1));
        result_low            = veorq_u64(result_low, vandq_u64(input_low[bit], mask));
        result_high           = veorq_u64(result_high, vandq_u64(input_high[bit], mask));
      }
    }
  }
  dst->values[0] = result_low;
  dst->values[1] = result_high;
}

void bf256_mul_inplace(bf256_t* lhs, const bf256_t* rhs) {
  bf256_mul(lhs, lhs, rhs);
}

void bf256_mul_64_inplace(bf256_t* lhs, bf64_t rhs) {
  uint64x2_t input_low[4];
  uint64x2_t input_high[4];
  input_low[0]  = lhs->values[0];
  input_high[0] = lhs->values[1];

  const uint64x2_t mod = bf_neon_low64(BF256_MODULUS);
  for (unsigned int i = 1; i < 4; ++i) {
    const uint64x2_t top      = input_high[i - 1];
    const uint64x2_t overflow = vdupq_n_u64(UINT64_C(0) - (bf_word(top, 1) >> 63));
    input_low[i]  = veorq_u64(bf_shift_left_1(input_low[i - 1]), vandq_u64(overflow, mod));
    input_high[i] = bf_shift_left_1_high(input_high[i - 1], input_low[i - 1]);
  }

  uint64x2_t result_low  = vdupq_n_u64(0);
  uint64x2_t result_high = vdupq_n_u64(0);
  for (unsigned int nibble = 0; nibble < 16; ++nibble, rhs <<= 4) {
    uint64x2_t overflow        = vextq_u64(vshrq_n_u64(result_high, 60), vdupq_n_u64(0), 1);
    const uint64x2_t preceding = vextq_u64(result_low, result_high, 1);
    result_high                = vorrq_u64(vshlq_n_u64(result_high, 4), vshrq_n_u64(preceding, 60));
    result_low                 = vorrq_u64(vshlq_n_u64(result_low, 4),
                                           vextq_u64(vdupq_n_u64(0), vshrq_n_u64(result_low, 60), 1));
    // Only four bits overflow, so folding them cannot cross another limb boundary.
    result_low = veorq_u64(result_low, bf_square_fold_10_5_2(overflow));
    for (unsigned int bit = 0; bit < 4; ++bit) {
      const uint64x2_t mask = vdupq_n_u64(UINT64_C(0) - ((rhs >> (60 + bit)) & 1));
      result_low            = veorq_u64(result_low, vandq_u64(input_low[bit], mask));
      result_high           = veorq_u64(result_high, vandq_u64(input_high[bit], mask));
    }
  }
  lhs->values[0] = result_low;
  lhs->values[1] = result_high;
}

void bf256_square(bf256_t* dst, const bf256_t* lhs) {
  const bf256_t square0 = bf_square_expand(lhs->values[0]);
  const bf256_t square1 = bf_square_expand(lhs->values[1]);
  bf256_reduce(dst, square0.values[0], square0.values[1], square1.values[0], square1.values[1]);
}

#if defined(FAEST_TESTS)
void bf256_sum_poly(bf256_t* dst, const bf256_t* xs) {
  *dst = xs[256 - 1];
  for (unsigned int i = 1; i != 256; ++i) {
    bf_dbl(dst->values, 4, BF256_MODULUS);
    bf256_add_inplace(dst, &xs[256 - 1 - i]);
  }
}
#endif

void bf256_sq_bit(bf256_t* out_tag, const bf256_t* in_tag) {
  bf256_add(&out_tag[7], &in_tag[6], &in_tag[7]);
  bf256_add(&out_tag[0], &in_tag[0], &in_tag[4]);
  bf256_add_inplace(&out_tag[0], &in_tag[6]);
  bf256_add(&out_tag[1], &in_tag[4], &out_tag[7]);
  bf256_add(&out_tag[2], &in_tag[1], &in_tag[5]);
  bf256_add(&out_tag[3], &out_tag[1], &in_tag[5]);
  bf256_add(&out_tag[4], &in_tag[2], &in_tag[4]);
  bf256_add_inplace(&out_tag[4], &in_tag[7]);
  bf256_add(&out_tag[5], &in_tag[5], &in_tag[6]);
  bf256_add(&out_tag[6], &in_tag[3], &in_tag[5]);
}

// Multiply a large field element by a smaller field element using register masks.
static void bf_mul_large(uint64x2_t* lhs, const uint64x2_t* rhs, unsigned int words,
                         unsigned int rhs_bits, uint64_t modulus) {
  const unsigned int count = (words + 1) / 2;
  uint64x2_t input[6];
  uint64x2_t result[6];
  for (unsigned int i = 0; i != count; ++i) {
    input[i]  = lhs[i];
    result[i] = vdupq_n_u64(0);
  }
  for (unsigned int bit = 0; bit != rhs_bits; ++bit) {
    const uint64_t word   = bf_word(rhs[bit / 128], (bit / 64) % 2);
    const uint64x2_t mask = vdupq_n_u64(UINT64_C(0) - ((word >> (bit % 64)) & 1));
    for (unsigned int i = 0; i != count; ++i) {
      result[i] = veorq_u64(result[i], vandq_u64(input[i], mask));
    }
    bf_dbl(input, words, modulus);
  }
  for (unsigned int i = 0; i != count; ++i) {
    lhs[i] = result[i];
  }
}

void bf384_mul_128_inplace(bf384_t* lhs, const bf128_t* rhs) {
  bf_mul_large(lhs->values, rhs->values, 6, 128, BF384_MODULUS);
}

void bf576_mul_192_inplace(bf576_t* lhs, const bf192_t* rhs) {
  bf_mul_large(lhs->values, rhs->values, 9, 192, BF576_MODULUS);
}

void bf768_mul_256_inplace(bf768_t* lhs, const bf256_t* rhs) {
  bf_mul_large(lhs->values, rhs->values, 12, 256, BF768_MODULUS);
}

#endif
