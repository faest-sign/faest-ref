/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_UTILS_H
#define FAEST_UTILS_H

#include <stdlib.h>
#include <stdint.h>

#include "compat.h"
#include "macros.h"

FAEST_BEGIN_C_DECL

static inline void xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out, size_t len) {
  for (size_t i = 0; i < len; i++) {
    out[i] = a[i] ^ b[i];
  }
}

static inline void masked_xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out,
                                       uint8_t mask_bit, size_t len) {
  uint8_t mask = -(mask_bit & 1);
  for (size_t i = 0; i < len; i++) {
    out[i] = a[i] ^ (b[i] & mask);
  }
}

#define get_bit(value, index) (((value) >> (index)) & 1)
#define set_bit(value, index) ((value) << (index))
#define ptr_get_bit(value, index) (((value)[(index) / 8] >> ((index) % 8)) & 1)
#define ptr_set_bit(dst, value, index)                                                             \
  do {                                                                                             \
    (dst)[(index) / 8] |= (value) << ((index) % 8);                                                \
  } while (0)

FAEST_END_C_DECL

#endif
