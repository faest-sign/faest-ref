/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern "C" {
#include "hash_shake.h"
}

#include <array>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(hash_shake)

BOOST_AUTO_TEST_CASE(shake_256) {
  constexpr std::array<uint8_t, 2> key = {0xab, 0xcd};

  hash_context ctx;
  hash_init(&ctx, 256);
  hash_update(&ctx, key.data(), key.size());
  hash_final(&ctx);

  constexpr std::array<uint8_t, 32> expected_output{0xc3, 0x1f, 0xcb, 0xe0, 0x76, 0x48, 0x87, 0x58,
                                                    0x73, 0x0d, 0xa7, 0xe5, 0xf8, 0x96, 0x4a, 0xfe,
                                                    0xfa, 0x65, 0x9a, 0x6f, 0x52, 0x6b, 0x6f, 0x9c,
                                                    0xa6, 0x7e, 0x59, 0x9f, 0x31, 0x8a, 0x7e, 0xa1};
  std::array<uint8_t, 32> output1, output2;
  hash_squeeze(&ctx, output1.data(), output1.size());
  hash_squeeze(&ctx, output2.data(), output2.size());
  hash_clear(&ctx);

  BOOST_TEST(output1 == expected_output);
  BOOST_TEST(output2 != expected_output);
}

BOOST_AUTO_TEST_SUITE_END()