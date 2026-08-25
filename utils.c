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

#if 0
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
#endif

void column_to_row_major(uint8_t** v, uint8_t** v_row_maj, unsigned int v_row_bits_len,
                         unsigned int v_col_bits_len) {
  // tranpose 64x64 blocks first
  const unsigned int full_rows_64 = v_row_bits_len & ~63;
  const unsigned int full_cols_64 = v_col_bits_len & ~63;

  for (unsigned int row = 0; row < full_rows_64; row += 64) {
    const unsigned int dst_byte = row / 8;
    for (unsigned int col = 0; col < full_cols_64; col += 64) {
      const unsigned int src_byte = col / 8;
      uint64_t block[64];

      for (unsigned int idx = 0; idx < 64; ++idx) {
        memcpy(&block[idx], v[row + idx] + src_byte, sizeof(block[idx]));
        block[idx] = le64toh(block[idx]);
      }
      transpose_64x64(block);
      for (unsigned int idx = 0; idx < 64; ++idx) {
        block[idx] = htole64(block[idx]);
        memcpy(v_row_maj[col + idx] + dst_byte, &block[idx], sizeof(block[idx]));
      }
    }
  }

  // then do 8x8 blocks
  const unsigned int full_rows_8 = v_row_bits_len & ~7;
  for (unsigned int row = 0; row < full_rows_8; row += 8) {
    const unsigned int dst_byte  = row / 8;
    const unsigned int first_col = row < full_rows_64 ? full_cols_64 : 0;
    for (unsigned int col = first_col; col < v_col_bits_len; col += 8) {
      const unsigned int src_byte = col / 8;
      const unsigned int columns  = MIN(v_col_bits_len - col, 8);
      uint64_t block              = 0;

      for (unsigned int bit = 0; bit < 8; ++bit) {
        block |= (uint64_t)v[row + bit][src_byte] << (bit * 8);
      }
      block = transpose_8x8(block);

      for (unsigned int bit = 0; bit < columns; ++bit) {
        v_row_maj[col + bit][dst_byte] = (uint8_t)(block >> (bit * 8));
      }
    }
  }

  // for the rest, perform bit operations
  const unsigned int remaining_rows = v_row_bits_len - full_rows_8;
  if (remaining_rows != 0) {
    const unsigned int dst_byte = full_rows_8 / 8;
    const uint8_t mask          = (UINT8_C(1) << remaining_rows) - UINT8_C(1);

    for (unsigned int col = 0; col < v_col_bits_len; ++col) {
      uint8_t value = 0;
      for (unsigned int bit = 0; bit < remaining_rows; ++bit) {
        value |= ptr_get_bit(v[full_rows_8 + bit], col) << bit;
      }
      v_row_maj[col][dst_byte] = (v_row_maj[col][dst_byte] & ~mask) | value;
    }
  }

#if !defined(NDEBUG)
  for (unsigned int row = 0; row < v_row_bits_len; row++) {
    for (unsigned int col = 0; col < v_col_bits_len; col++) {
      assert(ptr_get_bit(v[row], col) == ptr_get_bit(v_row_maj[col], row));
    }
  }
#endif
}

uint8_t** alloc_pointer_array(size_t rows, size_t columns) {
  uint8_t** array = malloc(rows * sizeof(uint8_t*));
  assert(array);
  array[0] = calloc(rows, columns);
  assert(array[0]);
  for (unsigned int i = 1; i < rows; ++i) {
    array[i] = array[0] + i * columns;
  }

  return array;
}

void free_pointer_array(uint8_t*** ptr) {
  free((*ptr)[0]);
  free(*ptr);
  *ptr = NULL;
}