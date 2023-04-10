/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../fields.h"

#include <NTL/GF2X.h>

#define BOOST_TEST_MODULE fields_test
#include <boost/test/included/unit_test.hpp>

using namespace NTL;

namespace {
  const GF2X bf8_residue{(1 << 8) | (1 << 4) | (1 << 3) | (1 << 1) | 1};
}

BOOST_AUTO_TEST_CASE(test_bf8_mul) {
  bf8_t lhs{0xFF}, rhs{0xEE};
  bf8_t result = bf8_mul(lhs, rhs);

  GF2X ntl_lhs{0xFF}, ntl_rhs{0xEE};
  auto ntl_result = MulMod(ntl_lhs, ntl_rhs, bf8_residue);
  BOOST_TEST(ntl_result == GF2X{result});
}

