/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "compat.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(utils)

BOOST_AUTO_TEST_CASE(log2) {
  const unsigned int depth = 7;
  BOOST_TEST(ceil_log2(1 << depth) == depth);
}

BOOST_AUTO_TEST_SUITE_END()
