/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "compat.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

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

BOOST_AUTO_TEST_SUITE_END()
