/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "fields.h"

#if FAEST_FIELDS == FAEST_FIELDS_USE_AVX2
#define TARGET ATTR_TARGET_AVX2 ATTR_TARGET_PCLMUL

// Move each 64-bit limb to its successor across the 128-bit lane boundary.
TARGET static inline __m256i bf_shift_word(__m256i value) {
  return _mm256_and_si256(_mm256_permute4x64_epi64(value, _MM_SHUFFLE(2, 1, 0, 0)),
                          _mm256_set_epi64x(-1, -1, -1, 0));
}

TARGET static inline __m256i bf_keep_words(__m256i value, unsigned int words) {
  const unsigned int remaining = words % 4;
  if (remaining) {
    const __m256i mask =
        _mm256_set_epi64x(0, -(long long)(remaining > 2), -(long long)(remaining > 1), -1);
    value = _mm256_and_si256(value, mask);
  }
  return value;
}

// At most four bits overflow. Multiply them by the low part of the modulus.
TARGET static inline uint64_t bf_fold_overflow(uint64_t overflow, uint64_t modulus) {
  return (modulus & -(overflow & 1)) ^ ((modulus << 1) & -((overflow >> 1) & 1)) ^
         ((modulus << 2) & -((overflow >> 2) & 1)) ^ ((modulus << 3) & -((overflow >> 3) & 1));
}

TARGET static inline void bf_shift_reduce(__m256i* value, unsigned int words, uint64_t modulus,
                                          int bits) {
  const unsigned int count = (words + 3) / 4;
  const uint64_t overflow  = bf_word(value[count - 1], (words - 1) % 4) >> (64 - bits);
  for (unsigned int i = count; i-- > 0;) {
    __m256i carry = bf_shift_word(_mm256_srli_epi64(value[i], 64 - bits));
    if (i) {
      const __m256i previous = _mm256_permute4x64_epi64(value[i - 1], _MM_SHUFFLE(3, 3, 3, 3));
      carry = _mm256_blend_epi32(carry, _mm256_srli_epi64(previous, 64 - bits), 0x03);
    }
    value[i] = _mm256_or_si256(_mm256_slli_epi64(value[i], bits), carry);
  }
  value[count - 1] = bf_keep_words(value[count - 1], words);
  value[0]         = _mm256_xor_si256(
      value[0], _mm256_set_epi64x(0, 0, 0, (long long)bf_fold_overflow(overflow, modulus)));
}

// Select a 128-bit polynomial block without aliasing the register storage.
TARGET static inline __m128i bf_half(const __m256i* value, unsigned int block) {
  return block & 1 ? _mm256_extracti128_si256(value[block / 2], 1)
                   : _mm256_castsi256_si128(value[block / 2]);
}

// Karatsuba multiplication of two full 128-bit polynomials.
TARGET static inline void bf_clmul128(__m128i* low, __m128i* high, __m128i lhs, __m128i rhs) {
  const __m128i low_product  = _mm_clmulepi64_si128(lhs, rhs, 0x00);
  const __m128i high_product = _mm_clmulepi64_si128(lhs, rhs, 0x11);
  const __m128i lhs_sum      = _mm_xor_si128(lhs, _mm_shuffle_epi32(lhs, _MM_SHUFFLE(1, 0, 3, 2)));
  const __m128i rhs_sum      = _mm_xor_si128(rhs, _mm_shuffle_epi32(rhs, _MM_SHUFFLE(1, 0, 3, 2)));
  const __m128i middle       = _mm_xor_si128(_mm_clmulepi64_si128(lhs_sum, rhs_sum, 0x00),
                                             _mm_xor_si128(low_product, high_product));
  *low                       = _mm_xor_si128(low_product, _mm_slli_si128(middle, 8));
  *high                      = _mm_xor_si128(high_product, _mm_srli_si128(middle, 8));
}

TARGET static inline void bf_clmul128_64(__m128i* low, __m128i* high, __m128i lhs, __m128i rhs) {
  const __m128i low_product = _mm_clmulepi64_si128(lhs, rhs, 0x00);
  const __m128i upper       = _mm_clmulepi64_si128(lhs, rhs, 0x01);
  *low                      = _mm_xor_si128(low_product, _mm_slli_si128(upper, 8));
  *high                     = _mm_srli_si128(upper, 8);
}

// Partial blocks are needed only by the generic large-field multiplier.
TARGET static inline void bf_clmul_partial(__m128i* low, __m128i* high, __m128i lhs, __m128i rhs,
                                           bool lhs_short, bool rhs_short) {
  if (lhs_short && rhs_short) {
    *low  = _mm_clmulepi64_si128(lhs, rhs, 0x00);
    *high = _mm_setzero_si128();
  } else if (rhs_short) {
    bf_clmul128_64(low, high, lhs, rhs);
  } else if (lhs_short) {
    bf_clmul128_64(low, high, rhs, lhs);
  } else {
    bf_clmul128(low, high, lhs, rhs);
  }
}

TARGET static void bf_mul_large(__m256i* dst, const __m256i* lhs, const __m256i* rhs,
                                unsigned int words, unsigned int rhs_words, uint64_t modulus) {
  const unsigned int blocks     = (words + 1) / 2;
  const unsigned int rhs_blocks = (rhs_words + 1) / 2;
  // The largest product is 768 by 256 bits, or eight 128-bit blocks.
  __m128i product[8];
  for (unsigned int i = 0; i != blocks + rhs_blocks; ++i) {
    product[i] = _mm_setzero_si128();
  }
  for (unsigned int i = 0; i != blocks; ++i) {
    __m128i left = bf_half(lhs, i);
    if ((words & 1) && i == blocks - 1) {
      left = _mm_unpacklo_epi64(left, _mm_setzero_si128());
    }
    for (unsigned int j = 0; j != rhs_blocks; ++j) {
      __m128i right = bf_half(rhs, j);
      if ((rhs_words & 1) && j == rhs_blocks - 1) {
        right = _mm_unpacklo_epi64(right, _mm_setzero_si128());
      }
      __m128i low, high;
      bf_clmul_partial(&low, &high, left, right, (words & 1) && i == blocks - 1,
                       (rhs_words & 1) && j == rhs_blocks - 1);
      product[i + j]     = _mm_xor_si128(product[i + j], low);
      product[i + j + 1] = _mm_xor_si128(product[i + j + 1], high);
    }
  }

  // Fold each upper 64-bit limb with the low polynomial of the field modulus.
  // Descending order also folds any carry that crosses the field boundary again.
  const __m128i mod = _mm_set_epi64x(0, (long long)modulus);
  for (unsigned int i = words + rhs_words; i-- > words;) {
    const __m128i high     = i & 1 ? _mm_srli_si128(product[i / 2], 8) : product[i / 2];
    const __m128i folded   = _mm_clmulepi64_si128(high, mod, 0x00);
    const unsigned int low = i - words;
    if (low & 1) {
      product[low / 2]     = _mm_xor_si128(product[low / 2], _mm_slli_si128(folded, 8));
      product[low / 2 + 1] = _mm_xor_si128(product[low / 2 + 1], _mm_srli_si128(folded, 8));
    } else {
      product[low / 2] = _mm_xor_si128(product[low / 2], folded);
    }
  }
  const unsigned int count = (words + 3) / 4;
  for (unsigned int i = 0; i != count; ++i) {
    // For a partial last block, do not read beyond the initialized product.
    const __m128i high = 4 * i + 2 < words ? product[2 * i + 1] : _mm_setzero_si128();
    dst[i]             = _mm256_inserti128_si256(_mm256_castsi128_si256(product[2 * i]), high, 1);
  }
  dst[count - 1] = bf_keep_words(dst[count - 1], words);
}

// In characteristic two the cross terms cancel: square each 64-bit half independently.
TARGET static inline __m256i bf_square_expand(__m128i value) {
  const __m128i low  = _mm_clmulepi64_si128(value, value, 0x00);
  const __m128i high = _mm_clmulepi64_si128(value, value, 0x11);
  return _mm256_inserti128_si256(_mm256_castsi128_si256(low), high, 1);
}

TARGET static inline __m256i bf_fold_7_2_1(__m256i high) {
  return _mm256_xor_si256(_mm256_xor_si256(high, _mm256_slli_epi64(high, 7)),
                          _mm256_xor_si256(_mm256_slli_epi64(high, 2), _mm256_slli_epi64(high, 1)));
}

TARGET static inline __m256i bf_carry_7_2_1(__m256i high) {
  return _mm256_xor_si256(
      _mm256_xor_si256(_mm256_srli_epi64(high, 57), _mm256_srli_epi64(high, 62)),
      _mm256_srli_epi64(high, 63));
}

TARGET static inline __m256i bf_fold_10_5_2(__m256i high) {
  return _mm256_xor_si256(_mm256_xor_si256(high, _mm256_slli_epi64(high, 10)),
                          _mm256_xor_si256(_mm256_slli_epi64(high, 5), _mm256_slli_epi64(high, 2)));
}

TARGET static inline __m256i bf_carry_10_5_2(__m256i high) {
  return _mm256_xor_si256(
      _mm256_xor_si256(_mm256_srli_epi64(high, 54), _mm256_srli_epi64(high, 59)),
      _mm256_srli_epi64(high, 62));
}

// Each reduction has a fixed field boundary and clears any representation padding.
TARGET static inline void bf128_reduce(bf128_t* dst, __m128i low, __m128i high) {
  const __m256i upper    = _mm256_zextsi128_si256(high);
  const __m256i carry    = bf_carry_7_2_1(upper);
  __m256i result         = _mm256_xor_si256(_mm256_zextsi128_si256(low), bf_fold_7_2_1(upper));
  result                 = _mm256_xor_si256(result, bf_shift_word(carry));
  const __m256i overflow = _mm256_zextsi128_si256(_mm_srli_si128(_mm256_castsi256_si128(carry), 8));
  result                 = _mm256_xor_si256(result, bf_fold_7_2_1(overflow));
  dst->values[0]         = _mm256_zextsi128_si256(_mm256_castsi256_si128(result));
}

TARGET static inline void bf192_reduce(bf192_t* dst, __m128i product0, __m128i product1,
                                       __m128i product2) {
  // Split the 384-bit product at bit 192, halfway through product1.
  const __m128i upper0 =
      _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(product1), _mm_castsi128_pd(product2), 1));
  const __m256i high =
      _mm256_inserti128_si256(_mm256_castsi128_si256(upper0), _mm_srli_si128(product2, 8), 1);
  __m256i low            = _mm256_inserti128_si256(_mm256_castsi128_si256(product0),
                                                   _mm_unpacklo_epi64(product1, _mm_setzero_si128()), 1);
  const __m256i carry    = bf_carry_7_2_1(high);
  low                    = _mm256_xor_si256(low, bf_fold_7_2_1(high));
  low                    = _mm256_xor_si256(low, bf_shift_word(carry));
  const __m256i overflow = _mm256_zextsi128_si256(
      _mm_unpacklo_epi64(_mm256_extracti128_si256(carry, 1), _mm_setzero_si128()));
  low            = _mm256_xor_si256(low, bf_fold_7_2_1(overflow));
  dst->values[0] = _mm256_and_si256(low, _mm256_set_epi64x(0, -1, -1, -1));
}

TARGET static inline void bf256_reduce(bf256_t* dst, __m256i low, __m256i high) {
  const __m256i carry = bf_carry_10_5_2(high);
  low                 = _mm256_xor_si256(low, bf_fold_10_5_2(high));
  low                 = _mm256_xor_si256(low, bf_shift_word(carry));
  const __m256i overflow =
      _mm256_zextsi128_si256(_mm_srli_si128(_mm256_extracti128_si256(carry, 1), 8));
  dst->values[0] = _mm256_xor_si256(low, bf_fold_10_5_2(overflow));
}

// GF(2^128) implementation

TARGET void bf128_mul(bf128_t* dst, const bf128_t* lhs, const bf128_t* rhs) {
  const __m128i left  = _mm256_castsi256_si128(lhs->values[0]);
  const __m128i right = _mm256_castsi256_si128(rhs->values[0]);
  __m128i low, high;
  bf_clmul128(&low, &high, left, right);
  bf128_reduce(dst, low, high);
}

TARGET void bf128_mul_inplace(bf128_t* lhs, const bf128_t* rhs) {
  bf128_mul(lhs, lhs, rhs);
}

TARGET void bf128_mul_64_inplace(bf128_t* lhs, bf64_t rhs) {
  const __m128i left  = _mm256_castsi256_si128(lhs->values[0]);
  const __m128i right = _mm_set_epi64x(0, (long long)rhs);
  __m128i low, high;
  bf_clmul128_64(&low, &high, left, right);
  bf128_reduce(lhs, low, high);
}

TARGET void bf128_square(bf128_t* dst, const bf128_t* lhs) {
  const __m256i square = bf_square_expand(_mm256_castsi256_si128(lhs->values[0]));
  bf128_reduce(dst, _mm256_castsi256_si128(square), _mm256_extracti128_si256(square, 1));
}

#if defined(FAEST_TESTS)
TARGET void bf128_sum_poly(bf128_t* dst, const bf128_t* xs) {
  *dst = xs[128 - 1];
  for (unsigned int i = 1; i != 128; ++i) {
    bf_shift_reduce(dst->values, 2, BF128_MODULUS, 1);
    bf128_add_inplace(dst, &xs[128 - 1 - i]);
  }
}
#endif

TARGET void bf128_sq_bit(bf128_t* out_tag, const bf128_t* in_tag) {
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

TARGET void bf192_mul(bf192_t* dst, const bf192_t* lhs, const bf192_t* rhs) {
  const __m128i a0 = _mm256_castsi256_si128(lhs->values[0]);
  const __m128i a1 = _mm_srli_si128(a0, 8);
  const __m128i a2 = _mm256_extracti128_si256(lhs->values[0], 1);
  const __m128i b0 = _mm256_castsi256_si128(rhs->values[0]);
  const __m128i b1 = _mm_srli_si128(b0, 8);
  const __m128i b2 = _mm256_extracti128_si256(rhs->values[0], 1);

  // Three diagonal products and three pairwise Karatsuba cross products.
  const __m128i d0 = _mm_clmulepi64_si128(a0, b0, 0x00);
  const __m128i d1 = _mm_clmulepi64_si128(a0, b0, 0x11);
  const __m128i d2 = _mm_clmulepi64_si128(a2, b2, 0x00);
  const __m128i c01 =
      _mm_xor_si128(_mm_clmulepi64_si128(_mm_xor_si128(a0, a1), _mm_xor_si128(b0, b1), 0x00),
                    _mm_xor_si128(d0, d1));
  const __m128i c02 =
      _mm_xor_si128(_mm_clmulepi64_si128(_mm_xor_si128(a0, a2), _mm_xor_si128(b0, b2), 0x00),
                    _mm_xor_si128(d0, d2));
  const __m128i c12 =
      _mm_xor_si128(_mm_clmulepi64_si128(_mm_xor_si128(a1, a2), _mm_xor_si128(b1, b2), 0x00),
                    _mm_xor_si128(d1, d2));
  const __m128i product0 = _mm_xor_si128(d0, _mm_slli_si128(c01, 8));
  const __m128i product1 = _mm_xor_si128(
      _mm_xor_si128(d1, c02), _mm_xor_si128(_mm_srli_si128(c01, 8), _mm_slli_si128(c12, 8)));
  const __m128i product2 = _mm_xor_si128(d2, _mm_srli_si128(c12, 8));
  bf192_reduce(dst, product0, product1, product2);
}

TARGET void bf192_mul_inplace(bf192_t* lhs, const bf192_t* rhs) {
  bf192_mul(lhs, lhs, rhs);
}

TARGET void bf192_mul_64_inplace(bf192_t* lhs, bf64_t rhs) {
  const __m128i left0 = _mm256_castsi256_si128(lhs->values[0]);
  const __m128i left1 = _mm256_extracti128_si256(lhs->values[0], 1);
  const __m128i right = _mm_set_epi64x(0, (long long)rhs);
  __m128i low, high;
  bf_clmul128_64(&low, &high, left0, right);
  high = _mm_xor_si128(high, _mm_clmulepi64_si128(left1, right, 0x00));
  bf192_reduce(lhs, low, high, _mm_setzero_si128());
}

TARGET void bf192_square(bf192_t* dst, const bf192_t* lhs) {
  const __m256i square = bf_square_expand(_mm256_castsi256_si128(lhs->values[0]));
  const __m128i upper  = _mm256_extracti128_si256(lhs->values[0], 1);
  const __m128i high   = _mm_clmulepi64_si128(upper, upper, 0x00);
  bf192_reduce(dst, _mm256_castsi256_si128(square), _mm256_extracti128_si256(square, 1), high);
}

#if defined(FAEST_TESTS)
TARGET void bf192_sum_poly(bf192_t* dst, const bf192_t* xs) {
  *dst = xs[192 - 1];
  for (unsigned int i = 1; i != 192; ++i) {
    bf_shift_reduce(dst->values, 3, BF192_MODULUS, 1);
    bf192_add_inplace(dst, &xs[192 - 1 - i]);
  }
}
#endif

TARGET void bf192_sq_bit(bf192_t* out_tag, const bf192_t* in_tag) {
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

TARGET void bf256_mul(bf256_t* dst, const bf256_t* lhs, const bf256_t* rhs) {
  const __m128i left0  = _mm256_castsi256_si128(lhs->values[0]);
  const __m128i left1  = _mm256_extracti128_si256(lhs->values[0], 1);
  const __m128i right0 = _mm256_castsi256_si128(rhs->values[0]);
  const __m128i right1 = _mm256_extracti128_si256(rhs->values[0], 1);

  // Karatsuba at both the 128-bit and 256-bit levels: nine PCLMUL operations.
  __m128i low0, low1, high0, high1, middle0, middle1;
  bf_clmul128(&low0, &low1, left0, right0);
  bf_clmul128(&high0, &high1, left1, right1);
  bf_clmul128(&middle0, &middle1, _mm_xor_si128(left0, left1), _mm_xor_si128(right0, right1));
  middle0 = _mm_xor_si128(middle0, _mm_xor_si128(low0, high0));
  middle1 = _mm_xor_si128(middle1, _mm_xor_si128(low1, high1));
  const __m256i low =
      _mm256_inserti128_si256(_mm256_castsi128_si256(low0), _mm_xor_si128(low1, middle0), 1);
  const __m256i high =
      _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_xor_si128(high0, middle1)), high1, 1);
  bf256_reduce(dst, low, high);
}

TARGET void bf256_mul_inplace(bf256_t* lhs, const bf256_t* rhs) {
  bf256_mul(lhs, lhs, rhs);
}

TARGET void bf256_mul_64_inplace(bf256_t* lhs, bf64_t rhs) {
  const __m128i left0 = _mm256_castsi256_si128(lhs->values[0]);
  const __m128i left1 = _mm256_extracti128_si256(lhs->values[0], 1);
  const __m128i right = _mm_set_epi64x(0, (long long)rhs);
  __m128i low0, low1, high0, high1;
  bf_clmul128_64(&low0, &low1, left0, right);
  bf_clmul128_64(&high0, &high1, left1, right);
  const __m256i low =
      _mm256_inserti128_si256(_mm256_castsi128_si256(low0), _mm_xor_si128(low1, high0), 1);
  bf256_reduce(lhs, low, _mm256_zextsi128_si256(high1));
}

TARGET void bf256_square(bf256_t* dst, const bf256_t* lhs) {
  const __m256i low  = bf_square_expand(_mm256_castsi256_si128(lhs->values[0]));
  const __m256i high = bf_square_expand(_mm256_extracti128_si256(lhs->values[0], 1));
  bf256_reduce(dst, low, high);
}

#if defined(FAEST_TESTS)
TARGET void bf256_sum_poly(bf256_t* dst, const bf256_t* xs) {
  *dst = xs[256 - 1];
  for (unsigned int i = 1; i != 256; ++i) {
    bf_shift_reduce(dst->values, 4, BF256_MODULUS, 1);
    bf256_add_inplace(dst, &xs[256 - 1 - i]);
  }
}
#endif

TARGET void bf256_sq_bit(bf256_t* out_tag, const bf256_t* in_tag) {
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

TARGET void bf384_mul_128_inplace(bf384_t* lhs, const bf128_t* rhs) {
  bf_mul_large(lhs->values, lhs->values, rhs->values, 6, 2, BF384_MODULUS);
}

TARGET void bf576_mul_192_inplace(bf576_t* lhs, const bf192_t* rhs) {
  bf_mul_large(lhs->values, lhs->values, rhs->values, 9, 3, BF576_MODULUS);
}

TARGET void bf768_mul_256_inplace(bf768_t* lhs, const bf256_t* rhs) {
  bf_mul_large(lhs->values, lhs->values, rhs->values, 12, 4, BF768_MODULUS);
}

#endif
