/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_UTILS_H
#define FAEST_UTILS_H

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#include "compat.h"
#include "macros.h"
#include "instances.h"

FAEST_BEGIN_C_DECL

#define get_bit(value, index) (((value) >> (index)) & 1)
#define set_bit(value, index) ((value) << (index))
#define ptr_get_bit(value, index) (((value)[(index) / 8] >> ((index) % 8)) & 1)
#define ptr_set_bit(dst, index, value)                                                             \
  do {                                                                                             \
    const unsigned int ptr_set_bit_index_mod_8 = (index) % 8;                                      \
    (dst)[(index) / 8] = ((dst)[(index) / 8] & ~(1 << ptr_set_bit_index_mod_8)) |                  \
                         ((value) << ptr_set_bit_index_mod_8);                                     \
  } while (0)
#define ptr_xor_bit(dst, index, value) (dst)[(index) / 8] ^= ((value) & 1u) << ((index) % 8)

#define ptr_u16_get_bit(value, index) ((uint8_t)(((value)[(index) / 16] >> ((index) % 16)) & 1))
#define ptr_u64_get_bit(value, index) ((uint8_t)(((value)[(index) / 64] >> ((index) % 64)) & 1))

// DecodeAllChall_3
bool decode_all_chall_3(uint16_t* decoded_chall, const uint8_t* chall,
                        const faest_paramset_t* params);

void xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out, size_t len);

void masked_xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out, uint8_t mask_bit,
                         size_t len);

#if 0
void print_u8_array(const char* label, const uint8_t* arr, size_t m);

void print_u8_array_bits(const char* label, const uint8_t* arr, size_t m);
#endif

void transpose_matrix(uint8_t** v, uint8_t** v_row_maj, unsigned int v_row_bits_len,
                      unsigned int v_col_bits_len);

uint8_t** alloc_pointer_array(size_t rows, size_t columns);
void free_pointer_array(uint8_t*** ptr);

FAEST_END_C_DECL

#endif
