/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(WITH_CONFIG_H)
#include <config.h>
#endif

#include "fields.h"

// GF(2^8) with X^8 + X^4 + X^3 + X^1 + 1
static const uint8_t bf8_modulus = (1 << 4) | (1 << 3) | (1 << 1) | 1;
// GF(2^64) with X^64 + X^4 + X^3 + X^1 + 1
static const uint64_t bf64_modulus = (1 << 4) | (1 << 3) | (1 << 1) | 1;
// GF(2^128) with X^128 + X^7 + X^2 + X^1 + 1
// static const uint64_t bf128_modulus = (1 << 7) | (1 << 2) | (1 << 1) | 1;
// GF(2^192) with X^192 + X^7 + X^2 + X^1 + 1
// static const uint64_t bf192_modulus = (1 << 7) | (1 << 2) | (1 << 1) | 1;
// GF(2^256) with X^256 + X^10 + X^5 + X^2 + 1
// static const uint64_t bf256_modulus = (1 << 10) | (1 << 5) | (1 << 2) | 1;

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

bf128_t bf128_mul(bf128_t lhs, bf128_t rhs) {
  return lhs;
}

bf192_t bf192_add(bf192_t lhs, bf192_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}

bf192_t bf192_mul(bf192_t lhs, bf192_t rhs) {
  return lhs;
}

bf256_t bf256_add(bf256_t lhs, bf256_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}

bf256_t bf256_mul(bf256_t lhs, bf256_t rhs) {
  return lhs;
}