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

void xorUint8Arr(const uint8_t* a, const uint8_t* b, uint8_t* out, size_t len);
void maskedXorUint8Arr(const uint8_t* a, const uint8_t* b, uint8_t* out, uint8_t mask_bit,
                       size_t len);

ATTR_CONST static inline uint8_t get_bit(uint8_t in, uint8_t index) {
  return (in >> index) & 0x01;
}

ATTR_CONST static inline uint8_t set_bit(uint8_t in, uint8_t index) {
  return (in << index);
}

static inline uint8_t ptr_get_bit(const uint8_t* in, unsigned int index) {
  return (in[index / 8] >> (index % 8)) & 1;
}

static inline void ptr_set_bit(uint8_t* dst, uint8_t in, unsigned int index) {
  dst[index / 8] |= in << (index % 8);
}

FAEST_END_C_DECL

#endif
