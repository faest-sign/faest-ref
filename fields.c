/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "fields.h"
#include "endian_compat.h"
#include "randomness.h"

#include <string.h>

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

// GF(2^8) implementation

bf8_t bf8_load(const uint8_t* src) {
  return *src;
}

void bf8_store(uint8_t* dst, bf8_t src) {
  *dst = src;
}

bf8_t bf8_rand() {
  bf8_t ret;
  rand_bytes(&ret, sizeof(ret));
  return ret;
}

bf8_t bf8_zero() {
  return 0;
}

bf8_t bf8_one() {
  return 1;
}

bf8_t bf8_add(bf8_t lhs, bf8_t rhs) {
  return lhs ^ rhs;
}

bf8_t bf8_mul(bf8_t lhs, bf8_t rhs) {
  bf8_t result = 0;
  for (unsigned int idx = 8; idx; --idx, rhs >>= 1) {
    result ^= (-(rhs & 1)) & lhs;
    const uint8_t mask = -((lhs >> 7) & 1);
    lhs                = (lhs << 1) ^ (mask & bf8_modulus);
  }
  return result;
}

bf8_t bf8_inv(bf8_t in) {
  bf8_t t1 = in;
  bf8_t t2 = in;
  for (size_t i = 0; i < 8 - 2; i++) {
    t2 = bf8_mul(t2, t2);
    t1 = bf8_mul(t1, t2);
  }
  return bf8_mul(t1, t1);
}

bf8_t bf8_from_bit(uint8_t bit) {
  return bit & 1;
}

// GF(2^64) implementation

bf64_t bf64_load(const uint8_t* src) {
  bf64_t ret;
  memcpy(&ret, src, sizeof(ret));
#if defined(FAEST_IS_BIG_ENDIAN)
  ret = le64toh(ret);
#endif
  return ret;
}

void bf64_store(uint8_t* dst, bf64_t src) {
#if defined(FAEST_IS_BIG_ENDIAN)
  src = htole64(src);
#endif
  memcpy(dst, &src, sizeof(src));
}

bf64_t bf64_rand() {
  bf64_t ret;
  rand_bytes((uint8_t*)&ret, sizeof(ret));
  return ret;
}

bf64_t bf64_zero() {
  return 0;
}

bf64_t bf64_one() {
  return 1;
}

bf64_t bf64_add(bf64_t lhs, bf64_t rhs) {
  return lhs ^ rhs;
}

bf64_t bf64_mul(bf64_t lhs, bf64_t rhs) {
  bf64_t result = 0;
  for (unsigned int idx = 64; idx; --idx, rhs >>= 1) {
    result ^= (-(rhs & 1)) & lhs;
    const uint64_t mask = -((lhs >> 63) & 1);
    lhs                 = (lhs << 1) ^ (mask & bf64_modulus);
  }
  return result;
}

bf64_t bf64_inv(bf64_t in) {
  bf64_t t1 = in;
  bf64_t t2 = in;
  for (size_t i = 0; i < 64 - 2; i++) {
    t2 = bf64_mul(t2, t2);
    t1 = bf64_mul(t1, t2);
  }
  return bf64_mul(t1, t1);
}

bf64_t bf64_from_bit(uint8_t bit) {
  return bit & 1;
}

// GF(2^128) implementation

static const uint8_t bf128_alpha[16] = {
    0x0d, 0xce, 0x60, 0x55, 0xac, 0xe8, 0x3f, 0xa1, 0x1c, 0x9a, 0x97, 0xa9, 0x55, 0x85, 0x3d, 0x05,
};

bf128_t bf128_byte_combine(const bf128_t* x) {
  bf128_t bf_out        = x[0];
  const bf128_t alpha   = bf128_load(bf128_alpha);
  bf128_t running_alpha = alpha;
  for (unsigned int i = 1; i < 8; ++i, running_alpha = bf128_mul(running_alpha, alpha)) {
    bf_out = bf128_add(bf_out, bf128_mul(x[i], running_alpha));
  }
  return bf_out;
}

bf128_t bf128_byte_combine_bits(uint8_t x) {
  bf128_t bf_out        = bf128_from_bit(x & 1);
  const bf128_t alpha   = bf128_load(bf128_alpha);
  bf128_t running_alpha = alpha;
  for (unsigned int i = 1; i < 8; ++i, running_alpha = bf128_mul(running_alpha, alpha)) {
    const bf128_t xi = bf128_from_bit(x >> i);
    bf_out           = bf128_add(bf_out, bf128_mul(xi, running_alpha));
  }
  return bf_out;
}

bf128_t bf128_load(const uint8_t* src) {
  bf128_t ret;
  for (unsigned int i = 0; i != ARRAY_SIZE(ret.values); ++i, src += sizeof(uint64_t)) {
    memcpy(&ret.values[i], src, sizeof(ret.values[i]));
    ret.values[i] = le64toh(ret.values[i]);
  }
  return ret;
}

bf128_t bf128_from_bf64(bf64_t src) {
  bf128_t ret;
  ret.values[0] = src;
  ret.values[1] = 0;
  return ret;
}

bf128_t bf128_from_bf8(bf8_t src) {
  bf128_t ret;
  ret.values[0] = src;
  ret.values[1] = 0;
  return ret;
}

void bf128_store(uint8_t* dst, bf128_t src) {
  for (unsigned int i = 0; i != ARRAY_SIZE(src.values); ++i, dst += sizeof(uint64_t)) {
    uint64_t tmp = htole64(src.values[i]);
    memcpy(dst, &tmp, sizeof(tmp));
  }
}

bf128_t bf128_rand() {
  bf128_t ret;
  rand_bytes((uint8_t*)&ret, sizeof(ret));
  return ret;
}

bf128_t bf128_zero() {
  bf128_t r = {0};
  return r;
}

bf128_t bf128_one() {
  bf128_t r = {{1, 0}};
  return r;
}

bf128_t bf128_add(bf128_t lhs, bf128_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}

ATTR_CONST
static bf128_t bf128_and(bf128_t lhs, bf128_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] &= rhs.values[i];
  }
  return lhs;
}

ATTR_CONST
static bf128_t bf128_shift_left_1(bf128_t value) {
  value.values[1] = (value.values[1] << 1) | (value.values[0] >> 63);
  value.values[0] = value.values[0] << 1;
  return value;
}

ATTR_CONST
static uint64_t bf128_bit_to_uint64_mask(bf128_t value, unsigned int bit) {
  const unsigned int byte_idx = bit / 64;
  const unsigned int bit_idx  = bit % 64;

  return -((value.values[byte_idx] >> bit_idx) & 1);
}

ATTR_CONST
static bf128_t bf128_bit_to_mask(bf128_t value, unsigned int bit) {
  bf128_t ret;
  ret.values[0] = ret.values[1] = bf128_bit_to_uint64_mask(value, bit);
  return ret;
}

bf128_t bf128_mul(bf128_t lhs, bf128_t rhs) {
  bf128_t result = {0};
  for (unsigned int idx = 0; idx != 128; ++idx) {
    result = bf128_add(result, bf128_and(bf128_bit_to_mask(rhs, idx), lhs));

    const uint64_t mask = bf128_bit_to_uint64_mask(lhs, 127);
    lhs                 = bf128_shift_left_1(lhs);
    lhs.values[0] ^= (mask & bf128_modulus);
  }
  return result;
}

bf128_t bf128_inv(bf128_t in) {
  bf128_t t1 = in;
  bf128_t t2 = in;
  for (size_t i = 0; i < 128 - 2; i++) {
    t2 = bf128_mul(t2, t2);
    t1 = bf128_mul(t1, t2);
  }
  return bf128_mul(t1, t1);
}

bf128_t bf128_from_bit(uint8_t bit) {
  return bf128_from_bf8(bit & 1);
}

bf128_t bf128_sum_poly(const bf128_t* xs) {
  bf128_t ret   = bf128_zero();
  bf128_t alpha = bf128_from_bf64(1);
  for (size_t i = 0; i < 128; ++i, alpha = bf128_shift_left_1(alpha)) {
    ret = bf128_add(ret, bf128_mul(alpha, xs[i]));
  }
  return ret;
}

// GF(2^192) implementation

static const uint8_t bf192_alpha[24] = {
    0x63, 0x97, 0x38, 0x6f, 0xd5, 0xa3, 0xc8, 0xcc, 0xea, 0xbd, 0x6e, 0x96,
    0x6c, 0xd7, 0x65, 0xe6, 0x62, 0x36, 0x6b, 0x0e, 0x14, 0xc8, 0x0b, 0x31,
};

bf192_t bf192_byte_combine(const bf192_t* x) {
  bf192_t bf_out        = x[0];
  const bf192_t alpha   = bf192_load(bf192_alpha);
  bf192_t running_alpha = alpha;
  for (unsigned int i = 1; i < 8; ++i, running_alpha = bf192_mul(running_alpha, alpha)) {
    bf_out = bf192_add(bf_out, bf192_mul(x[i], running_alpha));
  }
  return bf_out;
}

bf192_t bf192_byte_combine_bits(uint8_t x) {
  bf192_t bf_out        = bf192_from_bit(x & 1);
  const bf192_t alpha   = bf192_load(bf192_alpha);
  bf192_t running_alpha = alpha;
  for (unsigned int i = 1; i < 8; ++i, running_alpha = bf192_mul(running_alpha, alpha)) {
    const bf192_t xi = bf192_from_bit(x >> i);
    bf_out           = bf192_add(bf_out, bf192_mul(xi, running_alpha));
  }
  return bf_out;
}

bf192_t bf192_load(const uint8_t* src) {
  bf192_t ret;
  for (unsigned int i = 0; i != ARRAY_SIZE(ret.values); ++i, src += sizeof(uint64_t)) {
    memcpy(&ret.values[i], src, sizeof(ret.values[i]));
    ret.values[i] = le64toh(ret.values[i]);
  }
  return ret;
}

bf192_t bf192_from_bf64(bf64_t src) {
  bf192_t ret;
  ret.values[0] = src;
  ret.values[1] = 0;
  ret.values[2] = 0;
  return ret;
}

bf192_t bf192_from_bf8(bf8_t src) {
  bf192_t ret;
  ret.values[0] = src;
  ret.values[1] = 0;
  ret.values[2] = 0;
  return ret;
}

void bf192_store(uint8_t* dst, bf192_t src) {
  for (unsigned int i = 0; i != ARRAY_SIZE(src.values); ++i, dst += sizeof(uint64_t)) {
    uint64_t tmp = htole64(src.values[i]);
    memcpy(dst, &tmp, sizeof(tmp));
  }
}

bf192_t bf192_rand() {
  bf192_t ret;
  rand_bytes((uint8_t*)&ret, sizeof(ret));
  return ret;
}

bf192_t bf192_zero() {
  bf192_t r = {0};
  return r;
}

bf192_t bf192_one() {
  bf192_t r = {{1, 0, 0}};
  return r;
}

bf192_t bf192_add(bf192_t lhs, bf192_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}

ATTR_CONST
static bf192_t bf192_and(bf192_t lhs, bf192_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] &= rhs.values[i];
  }
  return lhs;
}

ATTR_CONST
static bf192_t bf192_shift_left_1(bf192_t value) {
  value.values[2] = (value.values[2] << 1) | (value.values[1] >> 63);
  value.values[1] = (value.values[1] << 1) | (value.values[0] >> 63);
  value.values[0] = value.values[0] << 1;
  return value;
}

ATTR_CONST
static uint64_t bf192_bit_to_uint64_mask(bf192_t value, unsigned int bit) {
  const unsigned int byte_idx = bit / 64;
  const unsigned int bit_idx  = bit % 64;

  return -((value.values[byte_idx] >> bit_idx) & 1);
}

ATTR_CONST
static bf192_t bf192_bit_to_mask(bf192_t value, unsigned int bit) {
  bf192_t ret;
  ret.values[0] = ret.values[1] = ret.values[2] = bf192_bit_to_uint64_mask(value, bit);
  return ret;
}

bf192_t bf192_mul(bf192_t lhs, bf192_t rhs) {
  bf192_t result = {0};
  for (unsigned int idx = 0; idx != 192; ++idx) {
    result = bf192_add(result, bf192_and(bf192_bit_to_mask(rhs, idx), lhs));

    const uint64_t mask = bf192_bit_to_uint64_mask(lhs, 191);
    lhs                 = bf192_shift_left_1(lhs);
    lhs.values[0] ^= (mask & bf192_modulus);
  }
  return result;
}

bf192_t bf192_inv(bf192_t in) {
  bf192_t t1 = in;
  bf192_t t2 = in;
  for (size_t i = 0; i < 192 - 2; i++) {
    t2 = bf192_mul(t2, t2);
    t1 = bf192_mul(t1, t2);
  }
  return bf192_mul(t1, t1);
}

bf192_t bf192_from_bit(uint8_t bit) {
  return bf192_from_bf8(bit & 1);
}

bf192_t bf192_sum_poly(const bf192_t* xs) {
  bf192_t ret   = bf192_zero();
  bf192_t alpha = bf192_from_bf64(1);
  for (size_t i = 0; i < 192; ++i, alpha = bf192_shift_left_1(alpha)) {
    ret = bf192_add(ret, bf192_mul(alpha, xs[i]));
  }
  return ret;
}

// GF(2^256) implementation

static const uint8_t bf256_alpha[32] = {
    0xe7, 0xfe, 0xde, 0x0b, 0x42, 0x88, 0x97, 0x96, 0x67, 0x4e, 0x47, 0xa0, 0x38, 0x8d, 0xd6, 0xbe,
    0x6a, 0xe1, 0xf1, 0xf8, 0x45, 0x98, 0x22, 0xdf, 0x33, 0x58, 0xc9, 0x20, 0xcf, 0xa8, 0xc9, 0x04,
};

bf256_t bf256_byte_combine(const bf256_t* x) {
  bf256_t bf_out        = x[0];
  const bf256_t alpha   = bf256_load(bf256_alpha);
  bf256_t running_alpha = alpha;
  for (unsigned int i = 1; i < 8; ++i, running_alpha = bf256_mul(running_alpha, alpha)) {
    bf_out = bf256_add(bf_out, bf256_mul(x[i], running_alpha));
  }
  return bf_out;
}

bf256_t bf256_byte_combine_bits(uint8_t x) {
  bf256_t bf_out        = bf256_from_bit(x & 1);
  const bf256_t alpha   = bf256_load(bf256_alpha);
  bf256_t running_alpha = alpha;
  for (unsigned int i = 1; i < 8; ++i, running_alpha = bf256_mul(running_alpha, alpha)) {
    const bf256_t xi = bf256_from_bit(x >> i);
    bf_out           = bf256_add(bf_out, bf256_mul(xi, running_alpha));
  }
  return bf_out;
}

bf256_t bf256_load(const uint8_t* src) {
  bf256_t ret;
  for (unsigned int i = 0; i != ARRAY_SIZE(ret.values); ++i, src += sizeof(uint64_t)) {
    memcpy(&ret.values[i], src, sizeof(ret.values[i]));
    ret.values[i] = le64toh(ret.values[i]);
  }
  return ret;
}

bf256_t bf256_from_bf64(bf64_t src) {
  bf256_t ret;
  ret.values[0] = src;
  ret.values[1] = 0;
  ret.values[2] = 0;
  ret.values[3] = 0;
  return ret;
}

bf256_t bf256_from_bf8(bf8_t src) {
  bf256_t ret;
  ret.values[0] = src;
  ret.values[1] = 0;
  ret.values[2] = 0;
  ret.values[3] = 0;
  return ret;
}

void bf256_store(uint8_t* dst, bf256_t src) {
  for (unsigned int i = 0; i != ARRAY_SIZE(src.values); ++i, dst += sizeof(uint64_t)) {
    uint64_t tmp = htole64(src.values[i]);
    memcpy(dst, &tmp, sizeof(tmp));
  }
}

bf256_t bf256_zero() {
  bf256_t r = {0};
  return r;
}

bf256_t bf256_one() {
  bf256_t r = {{1, 0, 0, 0}};
  return r;
}

bf256_t bf256_rand() {
  bf256_t ret;
  rand_bytes((uint8_t*)&ret, sizeof(ret));
  return ret;
}

bf256_t bf256_add(bf256_t lhs, bf256_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}

ATTR_CONST
static bf256_t bf256_and(bf256_t lhs, bf256_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] &= rhs.values[i];
  }
  return lhs;
}

ATTR_CONST
static bf256_t bf256_shift_left_1(bf256_t value) {
  value.values[3] = (value.values[3] << 1) | (value.values[2] >> 63);
  value.values[2] = (value.values[2] << 1) | (value.values[1] >> 63);
  value.values[1] = (value.values[1] << 1) | (value.values[0] >> 63);
  value.values[0] = value.values[0] << 1;
  return value;
}

ATTR_CONST
static uint64_t bf256_bit_to_uint64_mask(bf256_t value, unsigned int bit) {
  const unsigned int byte_idx = bit / 64;
  const unsigned int bit_idx  = bit % 64;

  return -((value.values[byte_idx] >> bit_idx) & 1);
}

ATTR_CONST
static bf256_t bf256_bit_to_mask(bf256_t value, unsigned int bit) {
  bf256_t ret;
  ret.values[0] = ret.values[1] = ret.values[2] = ret.values[3] =
      bf256_bit_to_uint64_mask(value, bit);
  return ret;
}

bf256_t bf256_mul(bf256_t lhs, bf256_t rhs) {
  bf256_t result = {0};
  for (unsigned int idx = 0; idx != 256; ++idx) {
    result = bf256_add(result, bf256_and(bf256_bit_to_mask(rhs, idx), lhs));

    const uint64_t mask = bf256_bit_to_uint64_mask(lhs, 255);
    lhs                 = bf256_shift_left_1(lhs);
    lhs.values[0] ^= (mask & bf256_modulus);
  }
  return result;
}

bf256_t bf256_inv(bf256_t in) {
  bf256_t t1 = in;
  bf256_t t2 = in;
  for (size_t i = 0; i < 256 - 2; i++) {
    t2 = bf256_mul(t2, t2);
    t1 = bf256_mul(t1, t2);
  }
  return bf256_mul(t1, t1);
}

bf256_t bf256_from_bit(uint8_t bit) {
  return bf256_from_bf8(bit & 1);
}

bf256_t bf256_sum_poly(const bf256_t* xs) {
  bf256_t ret   = bf256_zero();
  bf256_t alpha = bf256_from_bf64(1);
  for (size_t i = 0; i < 256; ++i, alpha = bf256_shift_left_1(alpha)) {
    ret = bf256_add(ret, bf256_mul(alpha, xs[i]));
  }
  return ret;
}
