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

bf8_t bf8_add(bf8_t lhs, bf8_t rhs);
bf8_t bf8_mul(bf8_t lhs, bf8_t rhs);

bf64_t bf64_add(bf64_t lhs, bf64_t rhs);
bf64_t bf64_mul(bf64_t lhs, bf64_t rhs);

FAEST_END_C_DECL

#endif
