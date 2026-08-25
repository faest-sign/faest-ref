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

BOOST_AUTO_TEST_CASE(test_column_to_row_major) {
  constexpr unsigned int rows     = 461;
  constexpr unsigned int columns  = 459;
  constexpr unsigned int row_size = (columns + 7) / 8;
  constexpr unsigned int col_size = (rows + 7) / 8;

  std::array<std::array<uint8_t, row_size>, rows> input{};
  std::array<uint8_t*, rows> input_rows{};
  for (unsigned int row = 0; row < rows; ++row) {
    input_rows[row] = input[row].data();
    for (unsigned int byte = 0; byte < row_size; ++byte) {
      input[row][byte] = static_cast<uint8_t>(17 * row + 29 * byte + 3);
    }
  }

  std::array<std::array<uint8_t, col_size>, columns> output{};
  std::array<uint8_t*, columns> output_rows{};
  for (unsigned int col = 0; col < columns; ++col) {
    output_rows[col] = output[col].data();
  }

  column_to_row_major(input_rows.data(), output_rows.data(), rows, columns);

  for (unsigned int row = 0; row < rows; ++row) {
    for (unsigned int col = 0; col < columns; ++col) {
      BOOST_TEST(ptr_get_bit(input[row].data(), col) == ptr_get_bit(output[col].data(), row));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
