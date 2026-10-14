/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "fields.h"

#if FAEST_FIELDS == FAEST_FIELDS_USE_SSE2
#include <emmintrin.h>

ATTR_TARGET_SSE2 static inline void bf_square_expand(uint64_t* dst, uint64_t value) {
  __m128i expanded = _mm_set_epi32(0, 0, (int32_t)(value >> 32), (int32_t)value);
  // Insert a zero byte after each input byte, then spread the bits within each 16-bit lane.
  expanded = _mm_unpacklo_epi8(expanded, _mm_setzero_si128());
  expanded =
      _mm_and_si128(_mm_or_si128(expanded, _mm_slli_epi16(expanded, 4)), _mm_set1_epi16(0x0f0f));
  expanded =
      _mm_and_si128(_mm_or_si128(expanded, _mm_slli_epi16(expanded, 2)), _mm_set1_epi16(0x3333));
  expanded =
      _mm_and_si128(_mm_or_si128(expanded, _mm_slli_epi16(expanded, 1)), _mm_set1_epi16(0x5555));
  _mm_store_si128((__m128i*)dst, expanded);
}

ATTR_TARGET_SSE2 static inline __m128i bf_square_fold_7_2_1(__m128i high) {
  return _mm_xor_si128(_mm_xor_si128(high, _mm_slli_epi64(high, 1)),
                       _mm_xor_si128(_mm_slli_epi64(high, 2), _mm_slli_epi64(high, 7)));
}

ATTR_TARGET_SSE2 static inline __m128i bf_square_carry_7_2_1(__m128i high) {
  return _mm_xor_si128(_mm_xor_si128(_mm_srli_epi64(high, 63), _mm_srli_epi64(high, 62)),
                       _mm_srli_epi64(high, 57));
}

ATTR_TARGET_SSE2 static inline __m128i bf_square_fold_10_5_2(__m128i high) {
  return _mm_xor_si128(_mm_xor_si128(high, _mm_slli_epi64(high, 2)),
                       _mm_xor_si128(_mm_slli_epi64(high, 5), _mm_slli_epi64(high, 10)));
}

ATTR_TARGET_SSE2 static inline __m128i bf_square_carry_10_5_2(__m128i high) {
  return _mm_xor_si128(_mm_xor_si128(_mm_srli_epi64(high, 62), _mm_srli_epi64(high, 59)),
                       _mm_srli_epi64(high, 54));
}

// Fold pairs of upper limbs and their carries into the lower half. The carry from the
// highest limb crosses the field boundary, so fold it once more into the lowest limb.
ATTR_TARGET_SSE2 static void bf_square_reduce_7_2_1_2(uint64_t* square) {
  const __m128i high  = _mm_load_si128((const __m128i*)(square + 2));
  const __m128i carry = bf_square_carry_7_2_1(high);

  __m128i low = _mm_xor_si128(_mm_load_si128((const __m128i*)square), bf_square_fold_7_2_1(high));
  low         = _mm_xor_si128(low, _mm_slli_si128(carry, 8));
  low         = _mm_xor_si128(low, bf_square_fold_7_2_1(_mm_srli_si128(carry, 8)));
  _mm_store_si128((__m128i*)square, low);
}

ATTR_TARGET_SSE2 static void bf_square_reduce_7_2_1_3(uint64_t* square) {
  const __m128i high0  = _mm_loadu_si128((const __m128i*)(square + 3));
  const __m128i high1  = _mm_loadl_epi64((const __m128i*)(square + 5));
  const __m128i carry0 = bf_square_carry_7_2_1(high0);
  const __m128i carry1 = bf_square_carry_7_2_1(high1);

  __m128i low0 = _mm_xor_si128(_mm_load_si128((const __m128i*)square), bf_square_fold_7_2_1(high0));
  __m128i low1 =
      _mm_xor_si128(_mm_loadl_epi64((const __m128i*)(square + 2)), bf_square_fold_7_2_1(high1));
  low0 = _mm_xor_si128(low0, _mm_slli_si128(carry0, 8));
  low0 = _mm_xor_si128(low0, bf_square_fold_7_2_1(carry1));
  low1 = _mm_xor_si128(low1, _mm_srli_si128(carry0, 8));
  _mm_store_si128((__m128i*)square, low0);
  _mm_storel_epi64((__m128i*)(square + 2), low1);
}

ATTR_TARGET_SSE2 static void bf_square_reduce_10_5_2(uint64_t* square) {
  const __m128i high0  = _mm_load_si128((const __m128i*)(square + 4));
  const __m128i high1  = _mm_load_si128((const __m128i*)(square + 6));
  const __m128i carry0 = bf_square_carry_10_5_2(high0);
  const __m128i carry1 = bf_square_carry_10_5_2(high1);

  __m128i low0 =
      _mm_xor_si128(_mm_load_si128((const __m128i*)square), bf_square_fold_10_5_2(high0));
  __m128i low1 =
      _mm_xor_si128(_mm_load_si128((const __m128i*)(square + 2)), bf_square_fold_10_5_2(high1));
  low0 = _mm_xor_si128(low0, _mm_slli_si128(carry0, 8));
  low0 = _mm_xor_si128(low0, bf_square_fold_10_5_2(_mm_srli_si128(carry1, 8)));
  low1 = _mm_xor_si128(low1, _mm_slli_si128(carry1, 8));
  low1 = _mm_xor_si128(low1, _mm_srli_si128(carry0, 8));
  _mm_store_si128((__m128i*)square, low0);
  _mm_store_si128((__m128i*)(square + 2), low1);
}

ATTR_TARGET_SSE2 static inline __m128i bf_shift_left_1_sse2(__m128i value) {
  return _mm_or_si128(_mm_slli_epi64(value, 1), _mm_slli_si128(_mm_srli_epi64(value, 63), 8));
}

ATTR_TARGET_SSE2 static inline __m128i bf_shift_left_1_high_sse2(__m128i high, __m128i low) {
  // Combine the preceding limb for each lane before extracting its carry bit.
  const __m128i preceding =
      _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(low), _mm_castsi128_pd(high), 1));
  return _mm_or_si128(_mm_slli_epi64(high, 1), _mm_srli_epi64(preceding, 63));
}

// Process four multiplier bits at a time using four fixed multiples of the input.
// Selection uses masks rather than a table indexed by multiplier bits.
// Both sizes are constant at each call site, allowing the compiler to specialize the loops.
ATTR_TARGET_SSE2 ATTR_ALWAYS_INLINE static inline void
bf128_mul_sse2(bf128_t* dst, const bf128_t* lhs, const __m128i* rhs, unsigned int words) {
  __m128i input_low[4];
  input_low[0] = lhs->values[0];

  const __m128i mod = _mm_set_epi32(0, 0, 0, (int32_t)BF128_MODULUS);
  for (unsigned int i = 1; i < 4; ++i) {
    const __m128i top      = input_low[i - 1];
    const __m128i overflow = _mm_srai_epi32(_mm_shuffle_epi32(top, _MM_SHUFFLE(3, 3, 3, 3)), 31);
    input_low[i] =
        _mm_xor_si128(bf_shift_left_1_sse2(input_low[i - 1]), _mm_and_si128(overflow, mod));
  }

  __m128i result_low = _mm_setzero_si128();
  for (unsigned int word = words; word-- > 0;) {
    uint64_t multiplier = bf_sse2_word(rhs[word / 2], word % 2);
    for (unsigned int nibble = 0; nibble < 16; ++nibble, multiplier <<= 4) {
      __m128i overflow = _mm_srli_si128(_mm_srli_epi64(result_low, 60), 8);

      result_low = _mm_or_si128(_mm_slli_epi64(result_low, 4),
                                _mm_slli_si128(_mm_srli_epi64(result_low, 60), 8));
      // Only four bits overflow, so folding them cannot cross another limb boundary.
      result_low = _mm_xor_si128(result_low, bf_square_fold_7_2_1(overflow));
      for (unsigned int bit = 0; bit < 4; ++bit) {
        const __m128i mask = _mm_set1_epi32(-(int)((multiplier >> (60 + bit)) & 1));
        result_low         = _mm_xor_si128(result_low, _mm_and_si128(input_low[bit], mask));
      }
    }
  }
  dst->values[0] = result_low;
}

ATTR_TARGET_SSE2 ATTR_ALWAYS_INLINE static inline void
bf192_mul_sse2(bf192_t* dst, const bf192_t* lhs, const __m128i* rhs, unsigned int words) {
  __m128i input_low[4];
  __m128i input_high[4];
  input_low[0]  = lhs->values[0];
  input_high[0] = lhs->values[1];

  const __m128i mod = _mm_set_epi32(0, 0, 0, (int32_t)BF192_MODULUS);
  for (unsigned int i = 1; i < 4; ++i) {
    const __m128i top      = input_high[i - 1];
    const __m128i overflow = _mm_srai_epi32(_mm_shuffle_epi32(top, _MM_SHUFFLE(1, 1, 1, 1)), 31);
    input_low[i] =
        _mm_xor_si128(bf_shift_left_1_sse2(input_low[i - 1]), _mm_and_si128(overflow, mod));
    input_high[i] = _mm_or_si128(_mm_slli_epi64(input_high[i - 1], 1),
                                 _mm_srli_si128(_mm_srli_epi64(input_low[i - 1], 63), 8));
  }

  __m128i result_low  = _mm_setzero_si128();
  __m128i result_high = _mm_setzero_si128();
  for (unsigned int word = words; word-- > 0;) {
    uint64_t multiplier = bf_sse2_word(rhs[word / 2], word % 2);
    for (unsigned int nibble = 0; nibble < 16; ++nibble, multiplier <<= 4) {
      __m128i overflow = _mm_srli_epi64(result_high, 60);
      result_high      = _mm_or_si128(_mm_slli_epi64(result_high, 4),
                                      _mm_srli_si128(_mm_srli_epi64(result_low, 60), 8));

      result_low = _mm_or_si128(_mm_slli_epi64(result_low, 4),
                                _mm_slli_si128(_mm_srli_epi64(result_low, 60), 8));
      // Only four bits overflow, so folding them cannot cross another limb boundary.
      result_low = _mm_xor_si128(result_low, bf_square_fold_7_2_1(overflow));
      for (unsigned int bit = 0; bit < 4; ++bit) {
        const __m128i mask = _mm_set1_epi32(-(int)((multiplier >> (60 + bit)) & 1));
        result_low         = _mm_xor_si128(result_low, _mm_and_si128(input_low[bit], mask));
        result_high        = _mm_xor_si128(result_high, _mm_and_si128(input_high[bit], mask));
      }
    }
  }
  dst->values[0] = result_low;
  dst->values[1] = result_high;
}

ATTR_TARGET_SSE2 ATTR_ALWAYS_INLINE static inline void
bf256_mul_sse2(bf256_t* dst, const bf256_t* lhs, const __m128i* rhs, unsigned int words) {
  __m128i input_low[4];
  __m128i input_high[4];
  input_low[0]  = lhs->values[0];
  input_high[0] = lhs->values[1];

  const __m128i mod = _mm_set_epi32(0, 0, 0, BF256_MODULUS);
  for (unsigned int i = 1; i < 4; ++i) {
    const __m128i top      = input_high[i - 1];
    const __m128i overflow = _mm_srai_epi32(_mm_shuffle_epi32(top, _MM_SHUFFLE(3, 3, 3, 3)), 31);
    input_low[i] =
        _mm_xor_si128(bf_shift_left_1_sse2(input_low[i - 1]), _mm_and_si128(overflow, mod));
    input_high[i] = bf_shift_left_1_high_sse2(input_high[i - 1], input_low[i - 1]);
  }

  __m128i result_low  = _mm_setzero_si128();
  __m128i result_high = _mm_setzero_si128();
  for (unsigned int word = words; word-- > 0;) {
    uint64_t multiplier = bf_sse2_word(rhs[word / 2], word % 2);
    for (unsigned int nibble = 0; nibble < 16; ++nibble, multiplier <<= 4) {
      __m128i overflow        = _mm_srli_si128(_mm_srli_epi64(result_high, 60), 8);
      const __m128i preceding = _mm_castpd_si128(
          _mm_shuffle_pd(_mm_castsi128_pd(result_low), _mm_castsi128_pd(result_high), 1));
      result_high = _mm_or_si128(_mm_slli_epi64(result_high, 4), _mm_srli_epi64(preceding, 60));
      result_low  = _mm_or_si128(_mm_slli_epi64(result_low, 4),
                                 _mm_slli_si128(_mm_srli_epi64(result_low, 60), 8));
      // Only four bits overflow, so folding them cannot cross another limb boundary.
      result_low = _mm_xor_si128(result_low, bf_square_fold_10_5_2(overflow));
      for (unsigned int bit = 0; bit < 4; ++bit) {
        const __m128i mask = _mm_set1_epi32(-(int)((multiplier >> (60 + bit)) & 1));
        result_low         = _mm_xor_si128(result_low, _mm_and_si128(input_low[bit], mask));
        result_high        = _mm_xor_si128(result_high, _mm_and_si128(input_high[bit], mask));
      }
    }
  }
  dst->values[0] = result_low;
  dst->values[1] = result_high;
}

void bf128_mul(bf128_t* dst, const bf128_t* lhs, const bf128_t* rhs) {
  bf128_mul_sse2(dst, lhs, rhs->values, 2);
}

void bf128_square(bf128_t* dst, const bf128_t* lhs) {
  ATTR_ALIGNED(16) uint64_t square[4];
  for (unsigned int i = 0; i != 2; ++i) {
    bf_square_expand(square + 2 * i, BF_VALUE(*lhs, i));
  }
  bf_square_reduce_7_2_1_2(square);

  dst->values[0] = _mm_load_si128((const __m128i*)square);
}

void bf128_mul_inplace(bf128_t* lhs, const bf128_t* rhs) {
  bf128_mul(lhs, lhs, rhs);
}

void bf128_mul_64_inplace(bf128_t* lhs, bf64_t rhs) {
  const __m128i multiplier = _mm_set_epi64x(0, (long long)rhs);
  bf128_mul_sse2(lhs, lhs, &multiplier, 1);
}

void bf192_mul(bf192_t* dst, const bf192_t* lhs, const bf192_t* rhs) {
  bf192_mul_sse2(dst, lhs, rhs->values, 3);
}

void bf192_square(bf192_t* dst, const bf192_t* lhs) {
  ATTR_ALIGNED(16) uint64_t square[6];
  for (unsigned int i = 0; i != 3; ++i) {
    bf_square_expand(square + 2 * i, BF_VALUE(*lhs, i));
  }
  bf_square_reduce_7_2_1_3(square);

  dst->values[0] = _mm_load_si128((const __m128i*)square);
  dst->values[1] = _mm_loadl_epi64((const __m128i*)(square + 2));
}

void bf192_mul_inplace(bf192_t* lhs, const bf192_t* rhs) {
  bf192_mul(lhs, lhs, rhs);
}

void bf192_mul_64_inplace(bf192_t* lhs, bf64_t rhs) {
  const __m128i multiplier = _mm_set_epi64x(0, (long long)rhs);
  bf192_mul_sse2(lhs, lhs, &multiplier, 1);
}

void bf256_mul(bf256_t* dst, const bf256_t* lhs, const bf256_t* rhs) {
  bf256_mul_sse2(dst, lhs, rhs->values, 4);
}

void bf256_square(bf256_t* dst, const bf256_t* lhs) {
  ATTR_ALIGNED(16) uint64_t square[8];
  for (unsigned int i = 0; i != 4; ++i) {
    bf_square_expand(square + 2 * i, BF_VALUE(*lhs, i));
  }
  bf_square_reduce_10_5_2(square);

  dst->values[0] = _mm_load_si128((const __m128i*)square);
  dst->values[1] = _mm_load_si128((const __m128i*)(square + 2));
}

void bf256_mul_inplace(bf256_t* lhs, const bf256_t* rhs) {
  bf256_mul(lhs, lhs, rhs);
}

void bf256_mul_64_inplace(bf256_t* lhs, bf64_t rhs) {
  const __m128i multiplier = _mm_set_epi64x(0, (long long)rhs);
  bf256_mul_sse2(lhs, lhs, &multiplier, 1);
}

// Double a polynomial and fold the overflow into its low register.
ATTR_TARGET_SSE2 static void bf_sse2_dbl(__m128i* value, unsigned int words, uint64_t modulus) {
  const unsigned int count = (words + 1) / 2;
  const uint64_t overflow  = bf_sse2_word(value[count - 1], (words - 1) % 2) >> 63;
  for (unsigned int i = count - 1; i > 0; --i) {
    value[i] = bf_shift_left_1_high_sse2(value[i], value[i - 1]);
  }
  value[0] = _mm_xor_si128(bf_shift_left_1_sse2(value[0]),
                           _mm_set_epi64x(0, (long long)(modulus & -overflow)));
  if (words & 1) {
    value[count - 1] = _mm_unpacklo_epi64(value[count - 1], _mm_setzero_si128());
  }
}

// Multiply a large field element by a smaller field element using register masks.
ATTR_TARGET_SSE2 static void bf_sse2_mul_large(__m128i* lhs, const __m128i* rhs, unsigned int words,
                                               unsigned int rhs_bits, uint64_t modulus) {
  const unsigned int count = (words + 1) / 2;
  __m128i input[6];
  __m128i result[6];
  for (unsigned int i = 0; i != count; ++i) {
    input[i]  = lhs[i];
    result[i] = _mm_setzero_si128();
  }
  for (unsigned int bit = 0; bit != rhs_bits; ++bit) {
    const uint64_t word = bf_sse2_word(rhs[bit / 128], (bit / 64) % 2);
    const __m128i mask  = _mm_set1_epi32(-(int)((word >> (bit % 64)) & 1));
    for (unsigned int i = 0; i != count; ++i) {
      result[i] = _mm_xor_si128(result[i], _mm_and_si128(input[i], mask));
    }
    bf_sse2_dbl(input, words, modulus);
  }
  for (unsigned int i = 0; i != count; ++i) {
    lhs[i] = result[i];
  }
}

void bf384_mul_128_inplace(bf384_t* lhs, const bf128_t* rhs) {
  bf_sse2_mul_large(lhs->values, rhs->values, 6, 128, BF384_MODULUS);
}

void bf576_mul_192_inplace(bf576_t* lhs, const bf192_t* rhs) {
  bf_sse2_mul_large(lhs->values, rhs->values, 9, 192, BF576_MODULUS);
}

void bf768_mul_256_inplace(bf768_t* lhs, const bf256_t* rhs) {
  bf_sse2_mul_large(lhs->values, rhs->values, 12, 256, BF768_MODULUS);
}

#if defined(FAEST_TESTS)
void bf128_sum_poly(bf128_t* dst, const bf128_t* xs) {
  *dst = xs[128 - 1];
  for (unsigned int i = 1; i != 128; ++i) {
    bf_sse2_dbl(dst->values, 2, BF128_MODULUS);
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

#if defined(FAEST_TESTS)
void bf192_sum_poly(bf192_t* dst, const bf192_t* xs) {
  *dst = xs[192 - 1];
  for (unsigned int i = 1; i != 192; ++i) {
    bf_sse2_dbl(dst->values, 3, BF192_MODULUS);
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

#if defined(FAEST_TESTS)
void bf256_sum_poly(bf256_t* dst, const bf256_t* xs) {
  *dst = xs[256 - 1];
  for (unsigned int i = 1; i != 256; ++i) {
    bf_sse2_dbl(dst->values, 4, BF256_MODULUS);
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

#endif
