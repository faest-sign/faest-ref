/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "utils.h"

#if defined(HAVE_AVX2) || defined(HAVE_SSE2)
#include "cpu.h"
#endif
#if defined(HAVE_AVX2)
#include <immintrin.h>
#elif defined(HAVE_SSE2)
#include <emmintrin.h>
#endif

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

// Based on Henry S. Warren, Jr., Hacker's Delight, "Transposing a Bit Matrix":
// https://github.com/hcs0/Hackers-Delight/blob/master/transpose8.c.txt
static inline uint64_t transpose_8x8(uint64_t value) {
  uint64_t tmp = (value ^ (value >> 7)) & UINT64_C(0x00aa00aa00aa00aa);
  value ^= tmp ^ (tmp << 7);
  tmp = (value ^ (value >> 14)) & UINT64_C(0x0000cccc0000cccc);
  value ^= tmp ^ (tmp << 14);
  tmp = (value ^ (value >> 28)) & UINT64_C(0x00000000f0f0f0f0);
  return value ^ tmp ^ (tmp << 28);
}

// 64-bit generalization of the recursive 32x32 mask-and-swap transpose
// https://github.com/hcs0/Hackers-Delight/blob/master/transpose32.c.txt
static void transpose_64x64(uint64_t* value) {
  uint64_t mask = UINT64_C(0x00000000ffffffff);
  for (unsigned int shift = 32; shift != 0;) {
    for (unsigned int idx = 0; idx < 64; idx = (idx + shift + 1) & ~shift) {
      const uint64_t tmp = ((value[idx] >> shift) ^ value[idx + shift]) & mask;
      value[idx] ^= tmp << shift;
      value[idx + shift] ^= tmp;
    }
    shift >>= 1;
    mask ^= mask << shift;
  }
}

#if defined(HAVE_SSE2)
// Compose 16x8 SSE2 transposes into a 128x128 transpose.
// https://mischasan.wordpress.com/2011/10/03/the-full-sse2-bit-matrix-transpose-routine/
ATTR_TARGET_SSE2 static void transpose_128x128(uint8_t** v, uint8_t** v_row_maj, unsigned int row,
                                               unsigned int col) {
  for (unsigned int row_offset = 0; row_offset < 128; row_offset += 16) {
    const unsigned int dst_byte = (row + row_offset) / 8;
    for (unsigned int col_offset = 0; col_offset < 128; col_offset += 8) {
      const unsigned int src_byte = (col + col_offset) / 8;

      ATTR_ALIGNED(16) uint8_t input[16];
      for (unsigned int idx = 0; idx < 16; ++idx) {
        input[idx] = v[row + row_offset + idx][src_byte];
      }

      __m128i block = _mm_load_si128((const __m128i*)input);
      for (unsigned int bit = 0; bit < 8; ++bit) {
        const uint16_t output = _mm_movemask_epi8(block);
        memcpy(v_row_maj[col + col_offset + 7 - bit] + dst_byte, &output, sizeof(output));
        block = _mm_slli_epi64(block, 1);
      }
    }
  }
}
#endif

#if defined(HAVE_AVX2)
// Apply one stage of a 32x32 mask-and-swap transpose in all eight 32-bit lanes.
#define TRANSPOSE_AVX2_32X32_STAGE(value, shift, mask_value)                                       \
  do {                                                                                             \
    const __m256i mask = _mm256_set1_epi32(mask_value);                                            \
    for (unsigned int idx = 0; idx < 32; idx = (idx + (shift) + 1) & ~(shift)) {                   \
      const __m256i tmp = _mm256_and_si256(                                                        \
          _mm256_xor_si256(_mm256_srli_epi32((value)[idx], shift), (value)[idx + (shift)]), mask); \
      (value)[idx]           = _mm256_xor_si256((value)[idx], _mm256_slli_epi32(tmp, shift));      \
      (value)[idx + (shift)] = _mm256_xor_si256((value)[idx + (shift)], tmp);                      \
    }                                                                                              \
  } while (0)

// Vectorized generalization of the recursive 32x32 transpose
ATTR_TARGET_AVX2 static void transpose_256x256(uint8_t** v, uint8_t** v_row_maj, unsigned int row,
                                               unsigned int col) {
  for (unsigned int row_offset = 0; row_offset < 256; row_offset += 32) {
    const unsigned int dst_byte = (row + row_offset) / 8;
    const unsigned int src_byte = col / 8;
    __m256i block[32];

    for (unsigned int idx = 0; idx < 32; ++idx) {
      block[idx] = _mm256_loadu_si256((const __m256i*)(v[row + row_offset + idx] + src_byte));
    }

    TRANSPOSE_AVX2_32X32_STAGE(block, 16, 0x0000ffff);
    TRANSPOSE_AVX2_32X32_STAGE(block, 8, 0x00ff00ff);
    TRANSPOSE_AVX2_32X32_STAGE(block, 4, 0x0f0f0f0f);
    TRANSPOSE_AVX2_32X32_STAGE(block, 2, 0x33333333);
    TRANSPOSE_AVX2_32X32_STAGE(block, 1, 0x55555555);

    for (unsigned int idx = 0; idx < 32; ++idx) {
      uint32_t output[8];
      _mm256_storeu_si256((__m256i*)output, block[idx]);
      for (unsigned int lane = 0; lane < 8; ++lane) {
        memcpy(v_row_maj[col + 32 * lane + idx] + dst_byte, &output[lane], sizeof(output[lane]));
      }
    }
  }
}
#endif

void column_to_row_major(uint8_t** v, uint8_t** v_row_maj, unsigned int v_row_bits_len,
                         unsigned int v_col_bits_len) {
  // transpose 256x256 blocks first
  unsigned int full_rows_256 = 0;
  unsigned int full_cols_256 = 0;
#if defined(HAVE_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    full_rows_256 = v_row_bits_len & ~255;
    full_cols_256 = v_col_bits_len & ~255;

    for (unsigned int row = 0; row < full_rows_256; row += 256) {
      for (unsigned int col = 0; col < full_cols_256; col += 256) {
        transpose_256x256(v, v_row_maj, row, col);
      }
    }
  }
#endif

  // transpose 128x128 blocks next
  unsigned int full_rows_128 = 0;
  unsigned int full_cols_128 = 0;
#if defined(HAVE_SSE2)
  if (CPU_SUPPORTS_SSE2) {
    full_rows_128 = v_row_bits_len & ~127;
    full_cols_128 = v_col_bits_len & ~127;

    for (unsigned int row = 0; row < full_rows_128; row += 128) {
      const unsigned int first_col = row < full_rows_256 ? full_cols_256 : 0;
      for (unsigned int col = first_col; col < full_cols_128; col += 128) {
        transpose_128x128(v, v_row_maj, row, col);
      }
    }
  }
#endif

  // transpose 64x64 blocks next
  const unsigned int full_rows_64 = v_row_bits_len & ~63;
  const unsigned int full_cols_64 = v_col_bits_len & ~63;

  for (unsigned int row = 0; row < full_rows_64; row += 64) {
    const unsigned int dst_byte  = row / 8;
    const unsigned int first_col = row < full_rows_128 ? full_cols_128 : 0;
    for (unsigned int col = first_col; col < full_cols_64; col += 64) {
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