/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_FIELDS_H
#define FAEST_FIELDS_H

#include "macros.h"

#include <stdint.h>

FAEST_BEGIN_C_DECL

typedef uint8_t bf8_t;
typedef uint64_t bf64_t;

typedef struct {
  uint64_t values[2];
} bf128_t;

typedef struct {
  uint64_t values[3];
} bf192_t;

typedef struct {
  uint64_t values[4];
} bf256_t;

bf8_t bf8_load(const uint8_t* src);
void bf8_store(uint8_t* dst, bf8_t src);
bf8_t bf8_rand();
ATTR_CONST bf8_t bf8_zero();
ATTR_CONST bf8_t bf8_one();
ATTR_CONST bf8_t bf8_add(bf8_t lhs, bf8_t rhs);
ATTR_CONST bf8_t bf8_mul(bf8_t lhs, bf8_t rhs);
ATTR_CONST bf8_t bf8_inv(bf8_t lhs);
ATTR_CONST bf8_t bf8_from_bit(uint8_t bit);

bf64_t bf64_load(const uint8_t* src);
void bf64_store(uint8_t* dst, bf64_t src);
bf64_t bf64_rand();
ATTR_CONST bf64_t bf64_zero();
ATTR_CONST bf64_t bf64_one();
ATTR_CONST bf64_t bf64_add(bf64_t lhs, bf64_t rhs);
ATTR_CONST bf64_t bf64_mul(bf64_t lhs, bf64_t rhs);
ATTR_CONST bf64_t bf64_inv(bf64_t lhs);
ATTR_CONST bf64_t bf64_from_bit(uint8_t bit);

bf128_t bf128_byte_combine(const bf128_t* x);
bf128_t bf128_byte_combine_bits(uint8_t x);
bf128_t bf128_load(const uint8_t* src);
ATTR_CONST bf128_t bf128_from_bf64(bf64_t src);
ATTR_CONST bf128_t bf128_from_bf8(bf8_t src);
void bf128_store(uint8_t* dst, bf128_t src);
bf128_t bf128_rand();
ATTR_CONST bf128_t bf128_zero();
ATTR_CONST bf128_t bf128_one();
ATTR_CONST bf128_t bf128_add(bf128_t lhs, bf128_t rhs);
ATTR_CONST bf128_t bf128_mul(bf128_t lhs, bf128_t rhs);
ATTR_CONST bf128_t bf128_inv(bf128_t lhs);
ATTR_CONST bf128_t bf128_from_bit(uint8_t bit);
bf128_t bf128_sum_poly(const bf128_t* xs);

bf192_t bf192_byte_combine(const bf192_t* x);
bf192_t bf192_byte_combine_bits(uint8_t x);
bf192_t bf192_load(const uint8_t* src);
ATTR_CONST bf192_t bf192_from_bf64(bf64_t src);
ATTR_CONST bf192_t bf192_from_bf8(bf8_t src);
void bf192_store(uint8_t* dst, bf192_t src);
bf192_t bf192_rand();
ATTR_CONST bf192_t bf192_zero();
ATTR_CONST bf192_t bf192_one();
ATTR_CONST bf192_t bf192_add(bf192_t lhs, bf192_t rhs);
ATTR_CONST bf192_t bf192_mul(bf192_t lhs, bf192_t rhs);
ATTR_CONST bf192_t bf192_inv(bf192_t lhs);
ATTR_CONST bf192_t bf192_from_bit(uint8_t bit);
bf192_t bf192_sum_poly(const bf192_t* xs);

bf256_t bf256_byte_combine(const bf256_t* x);
bf256_t bf256_byte_combine_bits(uint8_t x);
bf256_t bf256_load(const uint8_t* src);
ATTR_CONST bf256_t bf256_from_bf64(bf64_t src);
ATTR_CONST bf256_t bf256_from_bf8(bf8_t src);
void bf256_store(uint8_t* dst, bf256_t src);
bf256_t bf256_rand();
ATTR_CONST bf256_t bf256_zero();
ATTR_CONST bf256_t bf256_one();
ATTR_CONST bf256_t bf256_add(bf256_t lhs, bf256_t rhs);
ATTR_CONST bf256_t bf256_mul(bf256_t lhs, bf256_t rhs);
ATTR_CONST bf256_t bf256_inv(bf256_t lhs);
ATTR_CONST bf256_t bf256_from_bit(uint8_t bit);
bf256_t bf256_sum_poly(const bf256_t* xs);

FAEST_END_C_DECL

#endif
