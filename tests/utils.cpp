/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "utils.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <array>
#include <vector>

BOOST_AUTO_TEST_SUITE(utils)

namespace {
  void test_aligned_alloc(size_t alignment) {
    void* ptr      = faest_aligned_alloc(alignment, 3 * alignment);
    uintptr_t iptr = reinterpret_cast<uintptr_t>(ptr);
    BOOST_TEST((iptr & (alignment - 1)) == 0);
    faest_aligned_free(ptr);
  }

  constexpr size_t alignments[3] = {16, 32, 64};
} // namespace

BOOST_DATA_TEST_CASE(faest_aligned_alloc, alignments, alignment) {
  BOOST_TEST_CONTEXT("alignment: " << alignment) {
    test_aligned_alloc(alignment);
  }
}

BOOST_AUTO_TEST_CASE(test_transpose_matrix) {
  constexpr unsigned int shapes[][2] = {{7, 13}, {64, 64}, {128, 128}, {256, 256}, {461, 459}};
  for (const auto& shape : shapes) {
    const unsigned int rows    = shape[0];
    const unsigned int columns = shape[1];
    BOOST_TEST_CONTEXT("Shape: " << rows << "x" << columns) {
      const size_t input_stride  = (columns + 7) / 8 + 3;
      const size_t output_stride = (rows + 7) / 8 + 5;
      std::vector<uint8_t> input(rows * input_stride + 2);
      std::vector<uint8_t> output(columns * output_stride + 2, 0xa5);
      for (size_t i = 0; i < input.size(); ++i) {
        input[i] = static_cast<uint8_t>(17 * i + 29 * (i / input_stride) + 3);
      }
      const auto original_input = input;
      auto expected             = output;
      const byte_matrix_t src{input.data() + 1, input_stride};
      const byte_matrix_t dst{output.data() + 1, output_stride};
      for (unsigned int row = 0; row < rows; ++row) {
        for (unsigned int col = 0; col < columns; ++col) {
          ptr_set_bit(expected.data() + 1 + col * output_stride, row,
                      ptr_get_bit(matrix_row(&src, row), col));
        }
      }

      transpose_matrix(&src, &dst, rows, columns);

      BOOST_TEST(output == expected);
      BOOST_TEST(input == original_input);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
