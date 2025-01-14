/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vc.h"
#include "fields.h"
#include "compat.h"

#include <array>
#include <boost/test/unit_test.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(utils)

BOOST_AUTO_TEST_CASE(log2) {
  const unsigned int depth = 7;
  BOOST_TEST(ceil_log2(1 << depth) == depth);
}

BOOST_AUTO_TEST_SUITE_END()
