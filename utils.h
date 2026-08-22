/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_UTILS_H
#define FAEST_UTILS_H

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "compat.h"
#include "macros.h"
#include "instances.h"
#include <stdio.h>
#include <assert.h>

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

void unpack_uint8(uint8_t* unpacked, const uint8_t* packed, unsigned int row_bits,
                  unsigned int col_bits);

void pack_uint8(uint8_t* packed, const uint8_t* unpacked, unsigned int row_bits,
                unsigned int col_bits);

// DecodeAllChall_3
bool decode_all_chall_3(uint16_t* decoded_chall, const uint8_t* chall,
                        const faest_paramset_t* params);

void xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out, size_t len);

void xor_bit_to_pt(uint8_t* p, size_t i, uint8_t v);

void masked_xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out, uint8_t mask_bit,
                         size_t len);

void print_u8_array(const char* label, const uint8_t* arr, size_t m);

void print_u8_array_bits(const char* label, const uint8_t* arr, size_t m);

void column_to_row_major(uint8_t** v, uint8_t** v_row_maj, unsigned int v_row_bits_len,
                         unsigned int v_col_bits_len);

void gf2_poly_mul_ct(const uint8_t* a, size_t a_bits, const uint8_t* b, size_t b_bits,
                     uint8_t* out);

void gf2_poly_reduce_ct(const uint8_t* a, size_t a_bits, const uint8_t* m, size_t m_bits,
                        uint8_t* out);

FAEST_END_C_DECL

#endif
