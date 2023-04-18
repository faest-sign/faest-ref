/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(WITH_CONFIG_H)
#include <config.h>
#endif

#include "fields.h"

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

// GF(2^64) implementation

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

// GF(2^128) implementation

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

// GF(2^192) implementation

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

// GF(2^256) implementation

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
