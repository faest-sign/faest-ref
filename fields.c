/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "fields.h"
#include "randomness.h"
#include "utils.h"

// GF(2^8) with X^8 + X^4 + X^3 + X^1 + 1
#define bf8_modulus (UINT8_C((1 << 4) | (1 << 3) | (1 << 1) | 1))
// GF(2^64) with X^64 + X^4 + X^3 + X^1 + 1
#define bf64_modulus (UINT64_C((1 << 4) | (1 << 3) | (1 << 1) | 1))
// GF(2^128) with X^128 + X^7 + X^2 + X^1 + 1
#define bf128_modulus (UINT64_C((1 << 7) | (1 << 2) | (1 << 1) | 1))
// GF(2^192) with X^192 + X^7 + X^2 + X^1 + 1
#define bf192_modulus (UINT64_C((1 << 7) | (1 << 2) | (1 << 1) | 1))
// GF(2^256) with X^256 + X^10 + X^5 + X^2 + 1
#define bf256_modulus (UINT64_C((1 << 10) | (1 << 5) | (1 << 2) | 1))
// GF(2^384) with X^384 + X^12 + X^3 + X^2 + 1
#define bf384_modulus (UINT64_C((1 << 12) | (1 << 3) | (1 << 2) | 1))
// GF(2^576) with X^576 + X^13 + X^4 + X^3 + 1
#define bf576_modulus (UINT64_C((1 << 13) | (1 << 4) | (1 << 3) | 1))
// GF(2^768) with X^768 + X^19 + X^17 + X^4 + 1
#define bf768_modulus (UINT64_C((1 << 19) | (1 << 17) | (1 << 4) | 1))

#define U64C(x0, x1, x2, x3, x4, x5, x6, x7)                                                       \
  ((UINT64_C(x7) << 56) | (UINT64_C(x6) << 48) | (UINT64_C(x5) << 40) | (UINT64_C(x4) << 32) |     \
   (UINT64_C(x3) << 24) | (UINT64_C(x2) << 16) | (UINT64_C(x1) << 8) | UINT64_C(x0))

ATTR_CONST uint8_t bits_sq(uint8_t x) {
  uint8_t res = set_bit(get_bit(x, 0) ^ get_bit(x, 4) ^ get_bit(x, 6), 0);
  res |= set_bit(get_bit(x, 4) ^ get_bit(x, 6) ^ get_bit(x, 7), 1);
  res |= set_bit(get_bit(x, 1) ^ get_bit(x, 5), 2);
  res |= set_bit(get_bit(x, 4) ^ get_bit(x, 5) ^ get_bit(x, 6) ^ get_bit(x, 7), 3);
  res |= set_bit(get_bit(x, 2) ^ get_bit(x, 4) ^ get_bit(x, 7), 4);
  res |= set_bit(get_bit(x, 5) ^ get_bit(x, 6), 5);
  res |= set_bit(get_bit(x, 3) ^ get_bit(x, 5), 6);
  res |= set_bit(get_bit(x, 6) ^ get_bit(x, 7), 7);
  return res;
}

// GF(2^8) implementation

bf8_t bf8_rand(void) {
  bf8_t ret;
  rand_bytes(&ret, sizeof(ret));
  return ret;
}

bf8_t bf8_mul(bf8_t lhs, bf8_t rhs) {
  bf8_t result = -(rhs & 1) & lhs;
  for (unsigned int idx = 1; idx < 8; ++idx) {
    const uint8_t mask = -((lhs >> 7) & 1);
    lhs                = (lhs << 1) ^ (mask & bf8_modulus);
    result ^= -((rhs >> idx) & 1) & lhs;
  }
  return result;
}

bf8_t bf8_square(bf8_t lhs) {
  bf8_t result = -(lhs & 1) & lhs;
  bf8_t rhs    = lhs;
  for (unsigned int idx = 1; idx < 8; ++idx) {
    const uint8_t mask = -((lhs >> 7) & 1);
    lhs                = (lhs << 1) ^ (mask & bf8_modulus);
    result ^= -((rhs >> idx) & 1) & lhs;
  }
  return result;
}

bf8_t bf8_inv(bf8_t in) {
  const bf8_t t2   = bf8_square(in);
  const bf8_t t3   = bf8_mul(in, t2);
  const bf8_t t5   = bf8_mul(t3, t2);
  const bf8_t t7   = bf8_mul(t5, t2);
  const bf8_t t14  = bf8_square(t7);
  const bf8_t t28  = bf8_square(t14);
  const bf8_t t56  = bf8_square(t28);
  const bf8_t t63  = bf8_mul(t56, t7);
  const bf8_t t126 = bf8_square(t63);
  const bf8_t t252 = bf8_square(t126);
  return bf8_mul(t252, t2);
}

// GF(2^64) implementation

bf64_t bf64_rand(void) {
  bf64_t ret;
  rand_bytes((uint8_t*)&ret, sizeof(ret));
  return ret;
}

bf64_t bf64_mul(bf64_t lhs, bf64_t rhs) {
  bf64_t result = (-(rhs & 1)) & lhs;
  for (unsigned int idx = 1; idx != 64; ++idx) {
    const uint64_t mask = -((lhs >> 63) & 1);
    lhs                 = (lhs << 1) ^ (mask & bf64_modulus);
    result ^= (-((rhs >> idx) & 1)) & lhs;
  }
  return result;
}

#define bf64_bit_to_mask(value, bit) -((((uint64_t)(value)) >> (bit)) & 1)

// GF(2^128) implementation

static const bf128_t bf128_alpha[7] = {
    BF128C(U64C(0x0d, 0xce, 0x60, 0x55, 0xac, 0xe8, 0x3f, 0xa1),
           U64C(0x1c, 0x9a, 0x97, 0xa9, 0x55, 0x85, 0x3d, 0x05)),
    BF128C(U64C(0xe1, 0xae, 0x88, 0x34, 0xca, 0x59, 0x77, 0xec),
           U64C(0x84, 0xbb, 0xbf, 0x9c, 0x43, 0xb7, 0xf4, 0x4c)),
    BF128C(U64C(0xa8, 0x46, 0x39, 0x36, 0xae, 0x02, 0xcf, 0xbf),
           U64C(0xc6, 0xd2, 0x51, 0x7d, 0x4f, 0x60, 0xad, 0x35)),
    BF128C(U64C(0x49, 0x98, 0x2e, 0x3c, 0x48, 0x30, 0x83, 0x6b),
           U64C(0xfe, 0x22, 0xa2, 0x40, 0x46, 0x36, 0xcb, 0x0d)),
    BF128C(U64C(0xb4, 0x82, 0x1b, 0x7b, 0x27, 0x49, 0x2b, 0x25),
           U64C(0xa5, 0xde, 0x88, 0x1a, 0xe1, 0x10, 0x98, 0x54)),
    BF128C(U64C(0x22, 0xff, 0x21, 0x25, 0xef, 0xf2, 0x2b, 0xc7),
           U64C(0x75, 0x1f, 0x0c, 0x6c, 0x68, 0xa5, 0x81, 0xd6)),
    BF128C(U64C(0xbc, 0xf9, 0x36, 0xe1, 0x94, 0x8e, 0x7a, 0x7a),
           U64C(0xe0, 0x8f, 0xb7, 0x4f, 0x1a, 0x31, 0x50, 0x09)),
};

#if defined(FAEST_TESTS)
bf128_t bf128_get_alpha(unsigned int idx) {
  return bf128_alpha[idx];
}
#endif

bf128_t bf128_byte_combine(const bf128_t* x) {
  bf128_t bf_out = x[0];
  for (unsigned int i = 1; i < 8; ++i) {
    bf_out = bf128_add(bf_out, bf128_mul(x[i], bf128_alpha[i - 1]));
  }
  return bf_out;
}

void bf128_sq_bit(bf128_t* out_tag, const bf128_t* in_tag) {
  out_tag[0] = bf128_add(in_tag[0], bf128_add(in_tag[4], in_tag[6]));
  out_tag[1] = bf128_add(in_tag[4], bf128_add(in_tag[6], in_tag[7]));
  out_tag[2] = bf128_add(in_tag[1], in_tag[5]);
  out_tag[3] = bf128_add(bf128_add(in_tag[4], in_tag[5]), bf128_add(in_tag[6], in_tag[7]));
  out_tag[4] = bf128_add(in_tag[2], bf128_add(in_tag[4], in_tag[7]));
  out_tag[5] = bf128_add(in_tag[5], in_tag[6]);
  out_tag[6] = bf128_add(in_tag[3], in_tag[5]);
  out_tag[7] = bf128_add(in_tag[6], in_tag[7]);
}

void bf128_sq_bit_inplace(bf128_t* tag) {
  tag[0]           = bf128_add(tag[0], bf128_add(tag[4], tag[6]));
  const bf128_t i1 = tag[1];
  tag[1]           = bf128_add(tag[4], bf128_add(tag[6], tag[7]));
  const bf128_t i2 = tag[2];
  tag[2]           = bf128_add(i1, tag[5]);
  const bf128_t i3 = tag[3];
  tag[3]           = bf128_add(bf128_add(tag[4], tag[5]), bf128_add(tag[6], tag[7]));
  tag[4]           = bf128_add(i2, bf128_add(tag[4], tag[7]));
  const bf128_t i5 = tag[5];
  tag[5]           = bf128_add(tag[5], tag[6]);
  const bf128_t i6 = tag[6];
  tag[6]           = bf128_add(i3, i5);
  tag[7]           = bf128_add(i6, tag[7]);
}

bf128_t bf128_byte_combine_sq(const bf128_t* x) {
  bf128_t bf_tmp[8];
  bf128_sq_bit(bf_tmp, x);
  return bf128_byte_combine(bf_tmp);
}

bf128_t bf128_byte_combine_bits(uint8_t x) {
#if defined(HAVE_ATTR_VECTOR_SIZE)
  return bf128_from_bit(get_bit(x, 0)) ^ bf128_mul_bit(bf128_alpha[1 - 1], get_bit(x, 1)) ^
         bf128_mul_bit(bf128_alpha[2 - 1], get_bit(x, 2)) ^
         bf128_mul_bit(bf128_alpha[3 - 1], get_bit(x, 3)) ^
         bf128_mul_bit(bf128_alpha[4 - 1], get_bit(x, 4)) ^
         bf128_mul_bit(bf128_alpha[5 - 1], get_bit(x, 5)) ^
         bf128_mul_bit(bf128_alpha[6 - 1], get_bit(x, 6)) ^
         bf128_mul_bit(bf128_alpha[7 - 1], get_bit(x, 7));
#else
  bf128_t bf_out = bf128_from_bit(get_bit(x, 0));
  for (unsigned int i = 1; i < 8; ++i) {
    bf_out = bf128_add(bf_out, bf128_mul_bit(bf128_alpha[i - 1], get_bit(x, i)));
  }
  return bf_out;
#endif
}

bf128_t bf128_byte_combine_bits_sq(uint8_t x) {
  return bf128_byte_combine_bits(bits_sq(x));
}

bf128_t bf128_rand(void) {
  uint8_t buf[BF128_NUM_BYTES];
  rand_bytes(buf, sizeof(buf));
  return bf128_load(buf);
}

#if defined(HAVE_ATTR_VECTOR_SIZE)
#define bf128_and_64(lhs, rhs) ((lhs) & (rhs))
#else
ATTR_CONST
static inline bf128_t bf128_and_64(bf128_t lhs, bf64_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] &= rhs;
  }
  return lhs;
}
#endif

#if defined(HAVE_ATTR_VECTOR_SIZE)
#if __has_builtin(__builtin_shufflevector)
#define bf128_shift_right_64(v1) __builtin_shufflevector((v1), bf128_zero(), 2, 0)
#else
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_shift_right_64(bf128_t v1) {
  bf128_t ret;
  BF_VALUE(ret, 0) = 0;
  BF_VALUE(ret, 1) = BF_VALUE(v1, 0);
  return ret;
}
#endif

#define bf128_shift_left_1(value) ((value << 1) | bf128_shift_right_64(value >> 63))
#else
ATTR_CONST
static inline bf128_t bf128_shift_left_1(bf128_t value) {
  value.values[1] = (value.values[1] << 1) | (value.values[0] >> 63);
  value.values[0] = value.values[0] << 1;
  return value;
}
#endif

ATTR_CONST
static inline uint64_t bf128_bit_to_uint64_mask(bf128_t value, unsigned int bit) {
  const unsigned int byte_idx = bit / 64;
  const unsigned int bit_idx  = bit % 64;

  return -((BF_VALUE(value, byte_idx) >> bit_idx) & 1);
}

bf128_t bf128_mul(bf128_t lhs, bf128_t rhs) {
  bf128_t result = bf128_and_64(lhs, bf128_bit_to_uint64_mask(rhs, 0));
  for (unsigned int idx = 1; idx != 128; ++idx) {
    const uint64_t mask = bf128_bit_to_uint64_mask(lhs, 128 - 1);
    lhs                 = bf128_shift_left_1(lhs);
    BF_VALUE(lhs, 0) ^= (mask & bf128_modulus);

    result = bf128_add(result, bf128_and_64(lhs, bf128_bit_to_uint64_mask(rhs, idx)));
  }
  return result;
}

bf128_t bf128_mul_64(bf128_t lhs, bf64_t rhs) {
  bf128_t result = bf128_and_64(lhs, bf64_bit_to_mask(rhs, 0));
  for (unsigned int idx = 1; idx != 64; ++idx) {
    const uint64_t mask = bf128_bit_to_uint64_mask(lhs, 128 - 1);
    lhs                 = bf128_shift_left_1(lhs);
    BF_VALUE(lhs, 0) ^= (mask & bf128_modulus);

    result = bf128_add(result, bf128_and_64(lhs, bf64_bit_to_mask(rhs, idx)));
  }
  return result;
}

#if !defined(HAVE_ATTR_VECTOR_SIZE)
bf128_t bf128_mul_bit(bf128_t lhs, uint8_t rhs) {
  return bf128_and_64(lhs, -((uint64_t)rhs & 1));
}
#endif

ATTR_CONST static inline bf128_t bf128_dbl(bf128_t lhs) {
  uint64_t mask = bf128_bit_to_uint64_mask(lhs, 128 - 1);
  lhs           = bf128_shift_left_1(lhs);
  BF_VALUE(lhs, 0) ^= (mask & bf128_modulus);

  return lhs;
}

bf128_t bf128_sum_poly(const bf128_t* xs) {
  bf128_t ret = xs[128 - 1];
  for (size_t i = 1; i < 128; ++i) {
    ret = bf128_add(bf128_dbl(ret), xs[128 - 1 - i]);
  }
  return ret;
}

bf128_t bf128_sum_poly_bits(const uint8_t* xs) {
  bf128_t ret = bf128_from_bit(ptr_get_bit(xs, 128 - 1));
  for (size_t i = 1; i < 128; ++i) {
    ret = bf128_add(bf128_dbl(ret), bf128_from_bit(ptr_get_bit(xs, 128 - 1 - i)));
  }
  return ret;
}

// GF(2^192) implementation

static const bf192_t bf192_alpha[7] = {
    BF192C(U64C(0x63, 0x97, 0x38, 0x6f, 0xd5, 0xa3, 0xc8, 0xcc),
           U64C(0xea, 0xbd, 0x6e, 0x96, 0x6c, 0xd7, 0x65, 0xe6),
           U64C(0x62, 0x36, 0x6b, 0x0e, 0x14, 0xc8, 0x0b, 0x31)),
    BF192C(U64C(0xbb, 0x50, 0xf4, 0x7c, 0x9e, 0x61, 0x33, 0xb2),
           U64C(0x26, 0x3f, 0x63, 0xd5, 0x19, 0x1f, 0xf6, 0x7b),
           U64C(0x34, 0xdb, 0x91, 0xd4, 0x26, 0x37, 0x93, 0xda)),
    BF192C(U64C(0x0d, 0x8a, 0x39, 0xf5, 0x13, 0x2c, 0x6d, 0x9c),
           U64C(0x19, 0x8d, 0x32, 0x06, 0x77, 0xe3, 0x32, 0x82),
           U64C(0xf6, 0x4e, 0x75, 0x3c, 0x70, 0x0d, 0x3b, 0x0c)),
    BF192C(U64C(0x5d, 0xf7, 0x2b, 0xbd, 0x7c, 0x74, 0x20, 0xdd),
           U64C(0x2e, 0xd2, 0x58, 0x00, 0xab, 0x42, 0x55, 0x7a),
           U64C(0x51, 0x12, 0xbc, 0x94, 0x9c, 0x51, 0xec, 0x45)),
    BF192C(U64C(0xf8, 0x2b, 0xce, 0x8a, 0xe2, 0x0c, 0xd5, 0xd8),
           U64C(0x84, 0xbe, 0xde, 0x67, 0xb7, 0x8c, 0x16, 0x08),
           U64C(0x45, 0x70, 0xa6, 0x4b, 0x6a, 0x14, 0x7d, 0xd6)),
    BF192C(U64C(0xba, 0xe1, 0xd5, 0xee, 0x76, 0x9c, 0x0f, 0x97),
           U64C(0x48, 0x20, 0xd7, 0x5f, 0xae, 0xf7, 0xea, 0xf3),
           U64C(0x43, 0xea, 0x6c, 0x69, 0x5f, 0xbd, 0xa6, 0x29)),
    BF192C(U64C(0x71, 0x85, 0x06, 0x65, 0xc2, 0x5d, 0x94, 0xf5),
           U64C(0xd3, 0xe9, 0x06, 0x39, 0x62, 0xfd, 0x19, 0x60),
           U64C(0xb0, 0xc4, 0x87, 0x0f, 0x54, 0x56, 0x7c, 0xc7)),
};

#if defined(FAEST_TESTS)
bf192_t bf192_get_alpha(unsigned int idx) {
  return bf192_alpha[idx];
}
#endif

bf192_t bf192_byte_combine(const bf192_t* x) {
  bf192_t bf_out = x[0];
  for (unsigned int i = 1; i < 8; ++i) {
    bf_out = bf192_add(bf_out, bf192_mul(x[i], bf192_alpha[i - 1]));
  }
  return bf_out;
}

void bf192_sq_bit(bf192_t* out_tag, const bf192_t* in_tag) {
  out_tag[0] = bf192_add(in_tag[0], bf192_add(in_tag[4], in_tag[6]));
  out_tag[1] = bf192_add(in_tag[4], bf192_add(in_tag[6], in_tag[7]));
  out_tag[2] = bf192_add(in_tag[1], in_tag[5]);
  out_tag[3] = bf192_add(bf192_add(in_tag[4], in_tag[5]), bf192_add(in_tag[6], in_tag[7]));
  out_tag[4] = bf192_add(in_tag[2], bf192_add(in_tag[4], in_tag[7]));
  out_tag[5] = bf192_add(in_tag[5], in_tag[6]);
  out_tag[6] = bf192_add(in_tag[3], in_tag[5]);
  out_tag[7] = bf192_add(in_tag[6], in_tag[7]);
}

void bf192_sq_bit_inplace(bf192_t* tag) {
  tag[0]           = bf192_add(tag[0], bf192_add(tag[4], tag[6]));
  const bf192_t i1 = tag[1];
  tag[1]           = bf192_add(tag[4], bf192_add(tag[6], tag[7]));
  const bf192_t i2 = tag[2];
  tag[2]           = bf192_add(i1, tag[5]);
  const bf192_t i3 = tag[3];
  tag[3]           = bf192_add(bf192_add(tag[4], tag[5]), bf192_add(tag[6], tag[7]));
  tag[4]           = bf192_add(i2, bf192_add(tag[4], tag[7]));
  const bf192_t i5 = tag[5];
  tag[5]           = bf192_add(tag[5], tag[6]);
  const bf192_t i6 = tag[6];
  tag[6]           = bf192_add(i3, i5);
  tag[7]           = bf192_add(i6, tag[7]);
}

bf192_t bf192_byte_combine_sq(const bf192_t* x) {
  bf192_t bf_tmp[8];
  bf192_sq_bit(bf_tmp, x);
  return bf192_byte_combine(bf_tmp);
}

bf192_t bf192_byte_combine_bits(uint8_t x) {
#if defined(HAVE_ATTR_VECTOR_SIZE)
  return bf192_from_bit(get_bit(x, 0)) ^ bf192_mul_bit(bf192_alpha[1 - 1], get_bit(x, 1)) ^
         bf192_mul_bit(bf192_alpha[2 - 1], get_bit(x, 2)) ^
         bf192_mul_bit(bf192_alpha[3 - 1], get_bit(x, 3)) ^
         bf192_mul_bit(bf192_alpha[4 - 1], get_bit(x, 4)) ^
         bf192_mul_bit(bf192_alpha[5 - 1], get_bit(x, 5)) ^
         bf192_mul_bit(bf192_alpha[6 - 1], get_bit(x, 6)) ^
         bf192_mul_bit(bf192_alpha[7 - 1], get_bit(x, 7));
#else
  bf192_t bf_out = bf192_from_bit(get_bit(x, 0));
  for (unsigned int i = 1; i < 8; ++i) {
    bf_out = bf192_add(bf_out, bf192_mul_bit(bf192_alpha[i - 1], get_bit(x, i)));
  }
  return bf_out;
#endif
}

bf192_t bf192_byte_combine_bits_sq(uint8_t x) {
  return bf192_byte_combine_bits(bits_sq(x));
}

bf192_t bf192_rand(void) {
  uint8_t buf[BF192_NUM_BYTES];
  rand_bytes(buf, sizeof(buf));
  return bf192_load(buf);
}

#if defined(HAVE_ATTR_VECTOR_SIZE)
#define bf192_and_64(lhs, rhs) ((lhs) & (rhs))
#else
ATTR_CONST
static inline bf192_t bf192_and_64(bf192_t lhs, bf64_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] &= rhs;
  }
  return lhs;
}
#endif

#if defined(HAVE_ATTR_VECTOR_SIZE)
#if __has_builtin(__builtin_shufflevector)
#define bf192_shift_right_64(v1) __builtin_shufflevector((v1), bf192_zero(), 4, 0, 1, 5)
#else
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_shift_right_64(bf192_t v1) {
  bf192_t ret;
  BF_VALUE(ret, 0) = 0;
  BF_VALUE(ret, 1) = BF_VALUE(v1, 0);
  BF_VALUE(ret, 2) = BF_VALUE(v1, 1);
  BF_VALUE(ret, 3) = 0;
  return ret;
}
#endif
#endif

ATTR_CONST
static inline bf192_t bf192_shift_left_1(bf192_t value) {
#if defined(HAVE_ATTR_VECTOR_SIZE)
  const bf192_t mask = BF192C(0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff);
  return ((value << 1) | bf192_shift_right_64(value >> 63)) & mask;
#else
  value.values[2] = (value.values[2] << 1) | (value.values[1] >> 63);
  value.values[1] = (value.values[1] << 1) | (value.values[0] >> 63);
  value.values[0] = value.values[0] << 1;
#endif
  return value;
}

ATTR_CONST
static inline uint64_t bf192_bit_to_uint64_mask(bf192_t value, unsigned int bit) {
  const unsigned int byte_idx = bit / 64;
  const unsigned int bit_idx  = bit % 64;

  return -((BF_VALUE(value, byte_idx) >> bit_idx) & 1);
}

bf192_t bf192_mul(bf192_t lhs, bf192_t rhs) {
  bf192_t result = bf192_and_64(lhs, bf192_bit_to_uint64_mask(rhs, 0));
  for (unsigned int idx = 1; idx != 192; ++idx) {
    const uint64_t mask = bf192_bit_to_uint64_mask(lhs, 192 - 1);
    lhs                 = bf192_shift_left_1(lhs);
    BF_VALUE(lhs, 0) ^= (mask & bf192_modulus);

    result = bf192_add(result, bf192_and_64(lhs, bf192_bit_to_uint64_mask(rhs, idx)));
  }
  return result;
}

bf192_t bf192_mul_64(bf192_t lhs, bf64_t rhs) {
  bf192_t result = bf192_and_64(lhs, bf64_bit_to_mask(rhs, 0));
  for (unsigned int idx = 1; idx != 64; ++idx) {
    const uint64_t mask = bf192_bit_to_uint64_mask(lhs, 192 - 1);
    lhs                 = bf192_shift_left_1(lhs);
    BF_VALUE(lhs, 0) ^= (mask & bf192_modulus);

    result = bf192_add(result, bf192_and_64(lhs, bf64_bit_to_mask(rhs, idx)));
  }
  return result;
}

#if !defined(HAVE_ATTR_VECTOR_SIZE)
bf192_t bf192_mul_bit(bf192_t lhs, uint8_t rhs) {
  return bf192_and_64(lhs, -((uint64_t)rhs & 1));
}
#endif

ATTR_CONST static inline bf192_t bf192_dbl(bf192_t lhs) {
  uint64_t mask = bf192_bit_to_uint64_mask(lhs, 192 - 1);
  lhs           = bf192_shift_left_1(lhs);
  BF_VALUE(lhs, 0) ^= (mask & bf192_modulus);

  return lhs;
}

bf192_t bf192_sum_poly(const bf192_t* xs) {
  bf192_t ret = xs[192 - 1];
  for (size_t i = 1; i < 192; ++i) {
    ret = bf192_add(bf192_dbl(ret), xs[192 - 1 - i]);
  }
  return ret;
}

bf192_t bf192_sum_poly_bits(const uint8_t* xs) {
  bf192_t ret = bf192_from_bit(ptr_get_bit(xs, 192 - 1));
  for (size_t i = 1; i < 192; ++i) {
    ret = bf192_add(bf192_dbl(ret), bf192_from_bit(ptr_get_bit(xs, 192 - 1 - i)));
  }
  return ret;
}

// GF(2^256) implementation

static const bf256_t bf256_alpha[7] = {
    BF256C(U64C(0xe7, 0xfe, 0xde, 0x0b, 0x42, 0x88, 0x97, 0x96),
           U64C(0x67, 0x4e, 0x47, 0xa0, 0x38, 0x8d, 0xd6, 0xbe),
           U64C(0x6a, 0xe1, 0xf1, 0xf8, 0x45, 0x98, 0x22, 0xdf),
           U64C(0x33, 0x58, 0xc9, 0x20, 0xcf, 0xa8, 0xc9, 0x04)),
    BF256C(U64C(0xc1, 0x89, 0x22, 0xd5, 0x2a, 0xf5, 0x5a, 0xa9),
           U64C(0x2f, 0x07, 0x42, 0x2c, 0x8d, 0xc4, 0xa5, 0x2b),
           U64C(0xea, 0xb0, 0x00, 0x6c, 0x37, 0x0d, 0x4a, 0xd1),
           U64C(0xf1, 0x4a, 0x5b, 0x9c, 0x69, 0x4d, 0x4e, 0x06)),
    BF256C(U64C(0x1d, 0x9d, 0x80, 0x3f, 0x83, 0xb3, 0xda, 0x55),
           U64C(0x57, 0x0f, 0x3b, 0x53, 0x1e, 0x83, 0x71, 0x17),
           U64C(0x10, 0xac, 0x3f, 0xad, 0x3f, 0x57, 0x96, 0xfb),
           U64C(0x8d, 0xf6, 0x11, 0x70, 0xdb, 0xe3, 0x95, 0x61)),
    BF256C(U64C(0xd5, 0xcd, 0x1b, 0xb0, 0x19, 0x05, 0x01, 0xde),
           U64C(0xf6, 0xe3, 0x30, 0x1a, 0x91, 0x58, 0x27, 0x75),
           U64C(0x3f, 0xa0, 0x9e, 0x48, 0xb6, 0x78, 0x07, 0x2a),
           U64C(0x38, 0x88, 0x76, 0x4f, 0xd6, 0x4f, 0xc2, 0x56)),
    BF256C(U64C(0xb6, 0x30, 0x8a, 0xe9, 0x29, 0xf5, 0xc2, 0x98),
           U64C(0x82, 0x84, 0xf1, 0x40, 0xd4, 0xdb, 0xc4, 0x1b),
           U64C(0x81, 0xa9, 0x49, 0x7d, 0x94, 0x09, 0xbe, 0x2f),
           U64C(0xfc, 0x4f, 0x57, 0x71, 0x6d, 0x0b, 0x27, 0x22)),
    BF256C(U64C(0x0b, 0x67, 0x44, 0xde, 0xb9, 0xaf, 0x75, 0x9e),
           U64C(0xbc, 0xaf, 0xf1, 0x66, 0xc6, 0x66, 0xed, 0xac),
           U64C(0x7e, 0x1f, 0x99, 0xf2, 0x3f, 0x25, 0x01, 0xf0),
           U64C(0xf3, 0x29, 0xfa, 0xd1, 0x2f, 0x37, 0x3d, 0xc0)),
    BF256C(U64C(0x8b, 0xe8, 0x32, 0xb3, 0x98, 0xb6, 0x43, 0xba),
           U64C(0x0d, 0x6f, 0xb8, 0x25, 0xd6, 0xc4, 0x37, 0x52),
           U64C(0x45, 0x15, 0xe8, 0xf4, 0x2a, 0x2b, 0x65, 0x2f),
           U64C(0xb8, 0x7b, 0x6b, 0xd2, 0x09, 0xea, 0x3e, 0x13)),
};

#if defined(FAEST_TESTS)
bf256_t bf256_get_alpha(unsigned int idx) {
  return bf256_alpha[idx];
}
#endif

bf256_t bf256_byte_combine(const bf256_t* x) {
  bf256_t bf_out = x[0];
  for (unsigned int i = 1; i < 8; ++i) {
    bf_out = bf256_add(bf_out, bf256_mul(x[i], bf256_alpha[i - 1]));
  }
  return bf_out;
}

void bf256_sq_bit(bf256_t* out_tag, const bf256_t* in_tag) {
  out_tag[0] = bf256_add(in_tag[0], bf256_add(in_tag[4], in_tag[6]));
  out_tag[1] = bf256_add(in_tag[4], bf256_add(in_tag[6], in_tag[7]));
  out_tag[2] = bf256_add(in_tag[1], in_tag[5]);
  out_tag[3] = bf256_add(bf256_add(in_tag[4], in_tag[5]), bf256_add(in_tag[6], in_tag[7]));
  out_tag[4] = bf256_add(in_tag[2], bf256_add(in_tag[4], in_tag[7]));
  out_tag[5] = bf256_add(in_tag[5], in_tag[6]);
  out_tag[6] = bf256_add(in_tag[3], in_tag[5]);
  out_tag[7] = bf256_add(in_tag[6], in_tag[7]);
}

void bf256_sq_bit_inplace(bf256_t* tag) {
  tag[0]           = bf256_add(tag[0], bf256_add(tag[4], tag[6]));
  const bf256_t i1 = tag[1];
  tag[1]           = bf256_add(tag[4], bf256_add(tag[6], tag[7]));
  const bf256_t i2 = tag[2];
  tag[2]           = bf256_add(i1, tag[5]);
  const bf256_t i3 = tag[3];
  tag[3]           = bf256_add(bf256_add(tag[4], tag[5]), bf256_add(tag[6], tag[7]));
  tag[4]           = bf256_add(i2, bf256_add(tag[4], tag[7]));
  const bf256_t i5 = tag[5];
  tag[5]           = bf256_add(tag[5], tag[6]);
  const bf256_t i6 = tag[6];
  tag[6]           = bf256_add(i3, i5);
  tag[7]           = bf256_add(i6, tag[7]);
}

bf256_t bf256_byte_combine_sq(const bf256_t* x) {
  bf256_t bf_tmp[8];
  bf256_sq_bit(bf_tmp, x);
  return bf256_byte_combine(bf_tmp);
}

bf256_t bf256_byte_combine_bits(uint8_t x) {
#if defined(HAVE_ATTR_VECTOR_SIZE)
  return bf256_from_bit(get_bit(x, 0)) ^ bf256_mul_bit(bf256_alpha[1 - 1], get_bit(x, 1)) ^
         bf256_mul_bit(bf256_alpha[2 - 1], get_bit(x, 2)) ^
         bf256_mul_bit(bf256_alpha[3 - 1], get_bit(x, 3)) ^
         bf256_mul_bit(bf256_alpha[4 - 1], get_bit(x, 4)) ^
         bf256_mul_bit(bf256_alpha[5 - 1], get_bit(x, 5)) ^
         bf256_mul_bit(bf256_alpha[6 - 1], get_bit(x, 6)) ^
         bf256_mul_bit(bf256_alpha[7 - 1], get_bit(x, 7));
#else
  bf256_t bf_out = bf256_from_bit(get_bit(x, 0));
  for (unsigned int i = 1; i < 8; ++i) {
    bf_out = bf256_add(bf_out, bf256_mul_bit(bf256_alpha[i - 1], get_bit(x, i)));
  }
  return bf_out;
#endif
}

bf256_t bf256_byte_combine_bits_sq(uint8_t x) {
  return bf256_byte_combine_bits(bits_sq(x));
}

bf256_t bf256_rand(void) {
  uint8_t buf[BF256_NUM_BYTES];
  rand_bytes(buf, sizeof(buf));
  return bf256_load(buf);
}

#if defined(HAVE_ATTR_VECTOR_SIZE)
#define bf256_and_64(lhs, rhs) ((lhs) & (rhs))
#else
ATTR_CONST
static inline bf256_t bf256_and_64(bf256_t lhs, bf64_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] &= rhs;
  }
  return lhs;
}
#endif

#if defined(HAVE_ATTR_VECTOR_SIZE)
#if __has_builtin(__builtin_shufflevector)
#define bf256_shift_right_64(v1) __builtin_shufflevector((v1), bf256_zero(), 4, 0, 1, 2)
#else
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_shift_right_64(bf256_t v1) {
  bf256_t ret;
  BF_VALUE(ret, 0) = 0;
  BF_VALUE(ret, 1) = BF_VALUE(v1, 0);
  BF_VALUE(ret, 2) = BF_VALUE(v1, 1);
  BF_VALUE(ret, 3) = BF_VALUE(v1, 2);
  return ret;
}
#endif

#define bf256_shift_left_1(value) ((value << 1) | bf256_shift_right_64(value >> 63))
#else
ATTR_CONST
static inline bf256_t bf256_shift_left_1(bf256_t value) {
  value.values[3] = (value.values[3] << 1) | (value.values[2] >> 63);
  value.values[2] = (value.values[2] << 1) | (value.values[1] >> 63);
  value.values[1] = (value.values[1] << 1) | (value.values[0] >> 63);
  value.values[0] = value.values[0] << 1;
  return value;
}
#endif

ATTR_CONST ATTR_ALWAYS_INLINE static inline uint64_t bf256_bit_to_uint64_mask(bf256_t value,
                                                                              unsigned int bit) {
  const unsigned int byte_idx = bit / 64;
  const unsigned int bit_idx  = bit % 64;

  return -((BF_VALUE(value, byte_idx) >> bit_idx) & 1);
}

bf256_t bf256_mul(bf256_t lhs, bf256_t rhs) {
#if defined(HAVE_ATTR_VECTOR_SIZE)
  const bf256_t mod = BF256C(bf256_modulus, 0, 0, 0);
#endif
  bf256_t result = bf256_and_64(lhs, bf256_bit_to_uint64_mask(rhs, 0));
  for (unsigned int idx = 1; idx != 256; ++idx) {
    const uint64_t mask = bf256_bit_to_uint64_mask(lhs, 256 - 1);
    lhs                 = bf256_shift_left_1(lhs);
#if defined(HAVE_ATTR_VECTOR_SIZE)
    lhs ^= bf256_and_64(mod, mask);
#else
    BF_VALUE(lhs, 0) ^= mask & bf256_modulus;
#endif

    result = bf256_add(result, bf256_and_64(lhs, bf256_bit_to_uint64_mask(rhs, idx)));
  }
  return result;
}

bf256_t bf256_mul_64(bf256_t lhs, bf64_t rhs) {
#if defined(HAVE_ATTR_VECTOR_SIZE)
  const bf256_t mod = BF256C(bf256_modulus, 0, 0, 0);
#endif
  bf256_t result = bf256_and_64(lhs, bf64_bit_to_mask(rhs, 0));
  for (unsigned int idx = 1; idx != 64; ++idx) {
    const uint64_t mask = bf256_bit_to_uint64_mask(lhs, 256 - 1);
    lhs                 = bf256_shift_left_1(lhs);
#if defined(HAVE_ATTR_VECTOR_SIZE)
    lhs ^= bf256_and_64(mod, mask);
#else
    BF_VALUE(lhs, 0) ^= mask & bf256_modulus;
#endif

    result = bf256_add(result, bf256_and_64(lhs, bf64_bit_to_mask(rhs, idx)));
  }
  return result;
}

#if !defined(HAVE_ATTR_VECTOR_SIZE)
bf256_t bf256_mul_bit(bf256_t lhs, uint8_t rhs) {
  return bf256_and_64(lhs, -((uint64_t)rhs & 1));
}
#endif

ATTR_CONST static inline bf256_t bf256_dbl(bf256_t lhs) {
  uint64_t mask = bf256_bit_to_uint64_mask(lhs, 256 - 1);
  lhs           = bf256_shift_left_1(lhs);
#if defined(HAVE_ATTR_VECTOR_SIZE)
  const bf256_t mod = BF256C(bf256_modulus, 0, 0, 0);
  return lhs ^ bf256_and_64(mod, mask);
#else
  BF_VALUE(lhs, 0) ^= mask & bf256_modulus;
  return lhs;
#endif
}

bf256_t bf256_sum_poly(const bf256_t* xs) {
  bf256_t ret = xs[256 - 1];
  for (size_t i = 1; i < 256; ++i) {
    ret = bf256_add(bf256_dbl(ret), xs[256 - 1 - i]);
  }
  return ret;
}

bf256_t bf256_sum_poly_bits(const uint8_t* xs) {
  bf256_t ret = bf256_from_bit(ptr_get_bit(xs, 256 - 1));
  for (size_t i = 1; i < 256; ++i) {
    ret = bf256_add(bf256_dbl(ret), bf256_from_bit(ptr_get_bit(xs, 256 - 1 - i)));
  }
  return ret;
}

// GF(2^384)

#if defined(HAVE_ATTR_VECTOR_SIZE)
ATTR_CONST
static inline bf384_t bf384_and_64(bf384_t lhs, bf64_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.inner); ++i) {
    lhs.inner[i] &= rhs;
  }
  return lhs;
}
#else
ATTR_CONST
static inline bf384_t bf384_and_64(bf384_t lhs, bf64_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] &= rhs;
  }
  return lhs;
}
#endif

#if defined(HAVE_ATTR_VECTOR_SIZE)
#if __has_builtin(__builtin_shufflevector)
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf384_t bf384_shift_right_64(bf384_t v1) {
  bf384_t ret;
  ret.inner[0] = __builtin_shufflevector(v1.inner[0], bf128_zero(), 2, 0);
  ret.inner[1] = __builtin_shufflevector(v1.inner[1], bf128_zero(), 2, 0) |
                 __builtin_shufflevector(v1.inner[0], bf128_zero(), 1, 3);
  ret.inner[2] = __builtin_shufflevector(v1.inner[2], bf128_zero(), 2, 0) |
                 __builtin_shufflevector(v1.inner[1], bf128_zero(), 1, 3);
  return ret;
}
#else
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf384_t bf384_shift_right_64(bf384_t v1) {
  bf384_t ret;
  BF_VALUE(ret.inner[0], 0) = 0;
  BF_VALUE(ret.inner[0], 1) = BF_VALUE(v1.inner[0], 0);
  BF_VALUE(ret.inner[1], 0) = BF_VALUE(v1.inner[0], 1);
  BF_VALUE(ret.inner[1], 1) = BF_VALUE(v1.inner[1], 0);
  BF_VALUE(ret.inner[2], 0) = BF_VALUE(v1.inner[1], 1);
  BF_VALUE(ret.inner[2], 1) = BF_VALUE(v1.inner[2], 0);
  return ret;
}
#endif

ATTR_CONST
static inline bf384_t bf384_shift_left_1(bf384_t value) {
  const bf384_t rhs = bf384_shift_right_64(value);
  for (unsigned int i = 0; i != ARRAY_SIZE(value.inner); ++i) {
    value.inner[i] = (value.inner[i] << 1) | (rhs.inner[i] >> 63);
  }
  return value;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline uint64_t bf384_bit_to_uint64_mask(bf384_t value,
                                                                              unsigned int bit) {
  const unsigned int inner_idx = bit / 128;
  const unsigned int inner_bit = bit % 128;
  const unsigned int byte_idx  = inner_bit / 64;
  const unsigned int bit_idx   = inner_bit % 64;

  return -((BF_VALUE(value.inner[inner_idx], byte_idx) >> bit_idx) & 1);
}
#else
ATTR_CONST
static inline bf384_t bf384_shift_left_1(bf384_t value) {
  for (unsigned int i = ARRAY_SIZE(value.values) - 1; i; --i) {
    value.values[i] = (value.values[i] << 1) | (value.values[i - 1] >> 63);
  }
  value.values[0] = value.values[0] << 1;
  return value;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline uint64_t bf384_bit_to_uint64_mask(bf384_t value,
                                                                              unsigned int bit) {
  const unsigned int byte_idx = bit / 64;
  const unsigned int bit_idx  = bit % 64;

  return -((BF_VALUE(value, byte_idx) >> bit_idx) & 1);
}
#endif

bf384_t bf384_mul_128(bf384_t lhs, bf128_t rhs) {
#if defined(HAVE_ATTR_VECTOR_SIZE)
  const bf128_t mod = BF128C(bf384_modulus, 0);
#endif
  bf384_t result = bf384_and_64(lhs, bf128_bit_to_uint64_mask(rhs, 0));
  for (unsigned int idx = 1; idx != 128; ++idx) {
    const uint64_t mask = bf384_bit_to_uint64_mask(lhs, 384 - 1);
    lhs                 = bf384_shift_left_1(lhs);
#if defined(HAVE_ATTR_VECTOR_SIZE)
    lhs.inner[0] ^= bf128_and_64(mod, mask);
#else
    BF_VALUE(lhs, 0) ^= mask & bf384_modulus;
#endif

    result = bf384_add(result, bf384_and_64(lhs, bf128_bit_to_uint64_mask(rhs, idx)));
  }
  return result;
}

// GF(2^576)

#if defined(HAVE_ATTR_VECTOR_SIZE)
ATTR_CONST
static inline bf576_t bf576_and_64(bf576_t lhs, bf64_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.inner); ++i) {
    lhs.inner[i] &= rhs;
  }
  return lhs;
}
#else
ATTR_CONST
static inline bf576_t bf576_and_64(bf576_t lhs, bf64_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] &= rhs;
  }
  return lhs;
}
#endif

#if defined(HAVE_ATTR_VECTOR_SIZE)
#if __has_builtin(__builtin_shufflevector)
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf576_t bf576_shift_right_64(bf576_t v1) {
  bf576_t ret;
  ret.inner[0] = __builtin_shufflevector(v1.inner[0], bf256_zero(), 4, 0, 1, 7);
  ret.inner[1] = __builtin_shufflevector(v1.inner[1], bf256_zero(), 4, 0, 1, 7) |
                 __builtin_shufflevector(v1.inner[0], bf256_zero(), 2, 5, 6, 7);
  ret.inner[2] = __builtin_shufflevector(v1.inner[2], bf256_zero(), 4, 0, 1, 7) |
                 __builtin_shufflevector(v1.inner[1], bf256_zero(), 2, 5, 6, 7);
  return ret;
}
#else
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf576_t bf576_shift_right_64(bf576_t v1) {
  bf576_t ret;
  BF_VALUE(ret.inner[0], 0) = 0;
  BF_VALUE(ret.inner[0], 1) = BF_VALUE(v1.inner[0], 0);
  BF_VALUE(ret.inner[0], 2) = BF_VALUE(v1.inner[0], 1);
  BF_VALUE(ret.inner[1], 0) = BF_VALUE(v1.inner[0], 2);
  BF_VALUE(ret.inner[1], 1) = BF_VALUE(v1.inner[1], 0);
  BF_VALUE(ret.inner[1], 2) = BF_VALUE(v1.inner[1], 1);
  BF_VALUE(ret.inner[2], 0) = BF_VALUE(v1.inner[1], 2);
  BF_VALUE(ret.inner[2], 1) = BF_VALUE(v1.inner[2], 0);
  BF_VALUE(ret.inner[2], 2) = BF_VALUE(v1.inner[2], 1);
  return ret;
}
#endif

ATTR_CONST
static inline bf576_t bf576_shift_left_1(bf576_t value) {
  const bf576_t rhs = bf576_shift_right_64(value);
  for (unsigned int i = 0; i != ARRAY_SIZE(value.inner); ++i) {
    value.inner[i] = (value.inner[i] << 1) | (rhs.inner[i] >> 63);
  }
  return value;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline uint64_t bf576_bit_to_uint64_mask(bf576_t value,
                                                                              unsigned int bit) {
  const unsigned int inner_idx = bit / 192;
  const unsigned int inner_bit = bit % 192;
  const unsigned int byte_idx  = inner_bit / 64;
  const unsigned int bit_idx   = inner_bit % 64;

  return -((BF_VALUE(value.inner[inner_idx], byte_idx) >> bit_idx) & 1);
}
#else
ATTR_CONST
static inline bf576_t bf576_shift_left_1(bf576_t value) {
  for (unsigned int i = ARRAY_SIZE(value.values) - 1; i; --i) {
    value.values[i] = (value.values[i] << 1) | (value.values[i - 1] >> 63);
  }
  value.values[0] = value.values[0] << 1;
  return value;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline uint64_t bf576_bit_to_uint64_mask(bf576_t value,
                                                                              unsigned int bit) {
  const unsigned int byte_idx = bit / 64;
  const unsigned int bit_idx  = bit % 64;

  return -((BF_VALUE(value, byte_idx) >> bit_idx) & 1);
}
#endif

bf576_t bf576_mul_192(bf576_t lhs, bf192_t rhs) {
#if defined(HAVE_ATTR_VECTOR_SIZE)
  const bf192_t mod = BF192C(bf576_modulus, 0, 0);
#endif
  bf576_t result = bf576_and_64(lhs, bf192_bit_to_uint64_mask(rhs, 0));
  for (unsigned int idx = 1; idx != 192; ++idx) {
    const uint64_t mask = bf576_bit_to_uint64_mask(lhs, 576 - 1);
    lhs                 = bf576_shift_left_1(lhs);
#if defined(HAVE_ATTR_VECTOR_SIZE)
    lhs.inner[0] ^= bf192_and_64(mod, mask);
#else
    BF_VALUE(lhs, 0) ^= mask & bf576_modulus;
#endif

    result = bf576_add(result, bf576_and_64(lhs, bf192_bit_to_uint64_mask(rhs, idx)));
  }
  return result;
}

// GF(2^768)

#if defined(HAVE_ATTR_VECTOR_SIZE)
ATTR_CONST
static inline bf768_t bf768_and_64(bf768_t lhs, bf64_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.inner); ++i) {
    lhs.inner[i] &= rhs;
  }
  return lhs;
}
#else
ATTR_CONST
static inline bf768_t bf768_and_64(bf768_t lhs, bf64_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] &= rhs;
  }
  return lhs;
}
#endif

#if defined(HAVE_ATTR_VECTOR_SIZE)
#if __has_builtin(__builtin_shufflevector)
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf768_t bf768_shift_right_64(bf768_t v1) {
  bf768_t ret;
  ret.inner[0] = __builtin_shufflevector(v1.inner[0], bf256_zero(), 4, 0, 1, 2);
  ret.inner[1] = __builtin_shufflevector(v1.inner[1], bf256_zero(), 4, 0, 1, 2) |
                 __builtin_shufflevector(v1.inner[0], bf256_zero(), 3, 5, 6, 7);
  ret.inner[2] = __builtin_shufflevector(v1.inner[2], bf256_zero(), 4, 0, 1, 2) |
                 __builtin_shufflevector(v1.inner[1], bf256_zero(), 3, 5, 6, 7);
  return ret;
}
#else
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf768_t bf768_shift_right_64(bf768_t v1) {
  bf768_t ret;
  BF_VALUE(ret.inner[0], 0) = 0;
  BF_VALUE(ret.inner[0], 1) = BF_VALUE(v1.inner[0], 0);
  BF_VALUE(ret.inner[0], 2) = BF_VALUE(v1.inner[0], 1);
  BF_VALUE(ret.inner[0], 3) = BF_VALUE(v1.inner[0], 2);
  BF_VALUE(ret.inner[1], 0) = BF_VALUE(v1.inner[0], 3);
  BF_VALUE(ret.inner[1], 1) = BF_VALUE(v1.inner[1], 0);
  BF_VALUE(ret.inner[1], 2) = BF_VALUE(v1.inner[1], 1);
  BF_VALUE(ret.inner[1], 3) = BF_VALUE(v1.inner[1], 2);
  BF_VALUE(ret.inner[2], 0) = BF_VALUE(v1.inner[1], 3);
  BF_VALUE(ret.inner[2], 1) = BF_VALUE(v1.inner[2], 0);
  BF_VALUE(ret.inner[2], 2) = BF_VALUE(v1.inner[2], 1);
  BF_VALUE(ret.inner[2], 3) = BF_VALUE(v1.inner[2], 2);
  return ret;
}
#endif

ATTR_CONST
static inline bf768_t bf768_shift_left_1(bf768_t value) {
  bf768_t rhs = bf768_shift_right_64(value);
  for (unsigned int i = 0; i != ARRAY_SIZE(value.inner); ++i) {
    value.inner[i] = (value.inner[i] << 1) | (rhs.inner[i] >> 63);
  }
  return value;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline uint64_t bf768_bit_to_uint64_mask(bf768_t value,
                                                                              unsigned int bit) {
  const unsigned int inner_idx = bit / 256;
  const unsigned int inner_bit = bit % 256;
  const unsigned int byte_idx  = inner_bit / 64;
  const unsigned int bit_idx   = inner_bit % 64;

  return -((BF_VALUE(value.inner[inner_idx], byte_idx) >> bit_idx) & 1);
}
#else
ATTR_CONST
static inline bf768_t bf768_shift_left_1(bf768_t value) {
  for (unsigned int i = ARRAY_SIZE(value.values) - 1; i; --i) {
    value.values[i] = (value.values[i] << 1) | (value.values[i - 1] >> 63);
  }
  value.values[0] = value.values[0] << 1;
  return value;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline uint64_t bf768_bit_to_uint64_mask(bf768_t value,
                                                                              unsigned int bit) {
  const unsigned int byte_idx = bit / 64;
  const unsigned int bit_idx  = bit % 64;

  return -((BF_VALUE(value, byte_idx) >> bit_idx) & 1);
}
#endif

bf768_t bf768_mul_256(bf768_t lhs, bf256_t rhs) {
#if defined(HAVE_ATTR_VECTOR_SIZE)
  const bf256_t mod = BF256C(bf768_modulus, 0, 0, 0);
#endif
  bf768_t result = bf768_and_64(lhs, bf256_bit_to_uint64_mask(rhs, 0));
  for (unsigned int idx = 1; idx != 256; ++idx) {
    const uint64_t mask = bf768_bit_to_uint64_mask(lhs, 768 - 1);
    lhs                 = bf768_shift_left_1(lhs);
#if defined(HAVE_ATTR_VECTOR_SIZE)
    lhs.inner[0] ^= bf256_and_64(mod, mask);
#else
    BF_VALUE(lhs, 0) ^= mask & bf768_modulus;
#endif

    result = bf768_add(result, bf768_and_64(lhs, bf256_bit_to_uint64_mask(rhs, idx)));
  }
  return result;
}
