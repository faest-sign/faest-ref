/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(WITH_CONFIG_H)
#include <config.h>
#endif

#include "fields.h"

static const uint8_t bf8_modulus    = (1 << 4) | (1 << 3) | (1 << 1) | 1;
static const uint64_t bf64_modulus  = (1 << 4) | (1 << 3) | (1 << 1) | 1;
// static const uint64_t bf128_modulus = (1 << 7) | (1 << 2) | (1 << 1) | 1;
// static const uint64_t bf192_modulus = (1 << 7) | (1 << 2) | (1 << 1) | 1;
// static const uint64_t bf256_modulus = (1 << 10) | (1 << 5) | (1 << 2) | 1;


bf8_t bf8_add(bf8_t lhs, bf8_t rhs) {
  return lhs ^ rhs;
}

bf8_t bf8_mul(bf8_t lhs, bf8_t rhs) {
  bf8_t result = 0;
  for (unsigned int idx = 8; idx; --idx, rhs >>= 1) {
    result ^= (-(rhs & 1)) & lhs;
    const uint8_t mask = -((lhs >> 7) & 1);
    lhs = (lhs << 1) ^ (mask & bf8_modulus);
  }
  return result;
}

bf64_t bf64_add(bf64_t lhs, bf64_t rhs) {
  return lhs ^ rhs;
}

bf64_t bf64_mul(bf64_t lhs, bf64_t rhs) {
  bf64_t result = 0;
  for (unsigned int idx = 64; idx; --idx, rhs >>= 1) {
    result ^= (-(rhs & 1)) & lhs;
    const uint64_t mask = -((lhs >> 63) & 1);
    lhs = (lhs << 1) ^ (mask & bf64_modulus);
  }
  return result;
}
