/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "utils.h"

#include <assert.h>
#include <string.h>

static inline uint16_t num_rec_2(const uint8_t* v) {
  uint16_t r;
  memcpy(&r, v, sizeof(r));
  return le16toh(r);
}

void unpack_uint8(uint8_t* unpacked, const uint8_t* packed, unsigned int row_bits,
                  unsigned int col_bits) {

  unsigned int unpacked_bit_idx = 0;
  unsigned int unpacked_bytes   = row_bits * ((col_bits + 7) / 8);
  memset(unpacked, 0, unpacked_bytes);

  unsigned int col_bits_bytes = (col_bits + 7) / 8;

  for (unsigned int row_bit_idx = 0; row_bit_idx < row_bits; row_bit_idx++) {
    for (unsigned int col_bit_idx = 0; col_bit_idx < col_bits; col_bit_idx++) {

      ptr_set_bit(unpacked + row_bit_idx * col_bits_bytes, col_bit_idx,
                  ptr_get_bit(packed, unpacked_bit_idx));

      unpacked_bit_idx++;
    }
  }
}

void pack_uint8(uint8_t* packed, const uint8_t* unpacked, unsigned int row_bits,
                unsigned int col_bits) {

  unsigned int packed_bit_idx = 0;
  unsigned int col_bits_bytes = (col_bits + 7) / 8;
  for (unsigned int row_bit_idx = 0; row_bit_idx < row_bits; row_bit_idx++) {
    for (unsigned int col_bit_idx = 0; col_bit_idx < col_bits; col_bit_idx++) {
      ptr_set_bit(packed, packed_bit_idx,
                  ptr_get_bit(unpacked + row_bit_idx * col_bits_bytes, col_bit_idx));

      packed_bit_idx++;
    }
  }
}

// DecodeChall_3
static bool decode_chall_3(uint8_t* decoded_chall, const uint8_t* chall, unsigned int i,
                           const faest_paramset_t* params) {
  if (i >= params->tau) {
    return false;
  }

  const unsigned int t1 = params->tau1;
  const unsigned int k  = params->k;

  unsigned int lo;
  unsigned int hi;
  if (i < t1) {
    lo = i * k;
    hi = (i + 1) * k;
  } else {
    const unsigned int t = i - t1;

    lo = (t1 * k) + (t * (k - 1));
    hi = (t1 * k) + ((t + 1) * (k - 1));
  }

  assert(hi - lo == k || hi - lo == k - 1);
  // TODO: this could be implemented more efficiently using bit shifts
  for (unsigned int j = lo; j < hi; ++j) {
    ptr_set_bit(decoded_chall, j - lo, ptr_get_bit(chall, j));
  }
  return true;
}

bool decode_all_chall_3(uint16_t* decoded_chall, const uint8_t* chall,
                        const faest_paramset_t* params) {
  for (unsigned int i = 0; i != params->tau; ++i) {
    uint8_t tmp[2] = {0};
    if (!decode_chall_3(tmp, chall, i, params)) {
      return false;
    }
    decoded_chall[i] = num_rec_2(tmp);
  }
  return true;
}

void xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out, size_t len) {
  for (size_t i = 0; i < len; i++) {
    out[i] = a[i] ^ b[i];
  }
}

void masked_xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out, uint8_t mask_bit,
                         size_t len) {
  uint8_t mask = -(mask_bit & 1);
  for (size_t i = 0; i < len; i++) {
    out[i] = a[i] ^ (b[i] & mask);
  }
}

void print_u8_array(const char* label, const uint8_t* arr, size_t m) {
  printf("%s: ", label);
  for (size_t i = 0; i < m; i++) {
    printf("%02x ", arr[i]);
  }
  printf("\n");
}

void print_u8_array_bits(const char* label, const uint8_t* arr, size_t m) {
  printf("%s: ", label);
  for (size_t i = 0; i < m; i++) {
    for (int b = 0; b < 8; b++) {
      printf("%u", (arr[i] >> b) & 1u);
    }
    printf(" ");
  }
  printf("\n");
}

void column_to_row_major(uint8_t** v, uint8_t** v_row_maj, unsigned int v_row_bits_len,
                         unsigned int v_col_bits_len) {

  for (unsigned int row = 0; row < v_row_bits_len; row++) {
    for (unsigned int col = 0; col < v_col_bits_len; col++) {
      uint8_t bit = (v[row][col / 8] >> (col % 8)) & 1u;
      v_row_maj[col][row / 8] |= (uint8_t)(bit << (row % 8));
    }
  }

  for (unsigned int row = 0; row < v_row_bits_len; row++) {
    for (unsigned int col = 0; col < v_col_bits_len; col++) {
      assert(ptr_get_bit(v[row], col) == ptr_get_bit(v_row_maj[col], row));
    }
  }
}

void gf2_poly_mul_ct(const uint8_t* a, size_t a_bits, const uint8_t* b, size_t b_bits,
                     uint8_t* out) {
  size_t out_bits  = a_bits + b_bits - 1;
  size_t out_bytes = (out_bits + 7) / 8;
  memset(out, 0, out_bytes);

  for (size_t i = 0; i < a_bits; i++) {
    uint8_t mask = -ptr_get_bit(a, i);
    for (size_t j = 0; j < b_bits; j++) {
      uint8_t bj = ptr_get_bit(b, j);
      ptr_xor_bit(out, i + j, bj & mask);
    }
  }
}

void gf2_poly_reduce_ct(const uint8_t* a, size_t a_bits, const uint8_t* m, size_t m_bits,
                        uint8_t* out) {
  size_t a_bytes   = (a_bits + 7) / 8;
  size_t out_bytes = (m_bits + 7) / 8;

  memset(out, 0, out_bytes);
  memcpy(out, a, a_bytes);

  size_t mdeg = m_bits - 1;
  for (size_t i = a_bits; i-- > mdeg;) {
    uint8_t mask = -ptr_get_bit(out, i);
    size_t shift = i - mdeg;

    for (size_t j = 0; j < m_bits; j++) {
      uint8_t mj = ptr_get_bit(m, j);
      ptr_xor_bit(out, shift + j, mj & mask);
    }

    if (i == 0) {
      break;
    }
  }
}
