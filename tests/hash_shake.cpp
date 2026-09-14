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
  constexpr std::array<uint8_t, 2> data1{0xab, 0xcd};
  constexpr uint32_t data2{0x1234};

  hash_context ctx;
  hash_init(&ctx, 256);
  hash_update(&ctx, data1.data(), data1.size());
  hash_update_uint32_le(&ctx, data2);
  hash_final(&ctx);

  constexpr std::array<uint8_t, 32> expected_output{
      0xfd, 0xd0, 0x50, 0x14, 0x9b, 0xaf, 0xcf, 0xfa, 0x03, 0xc4, 0x8b,
      0xa5, 0xc1, 0x5f, 0xfa, 0x86, 0xc2, 0xcc, 0x53, 0x48, 0xe8, 0xc5,
      0xe0, 0x30, 0x93, 0x77, 0xbf, 0x51, 0x4f, 0x56, 0xf7, 0x98,
  };

  std::array<uint8_t, 32> output1{0}, output2{0};
  hash_squeeze(&ctx, output1.data(), output1.size());
  hash_squeeze(&ctx, output2.data(), output2.size());
  hash_clear(&ctx);

  BOOST_TEST(output1 == expected_output);
  BOOST_TEST(output2 != expected_output);
}

BOOST_AUTO_TEST_SUITE_END()