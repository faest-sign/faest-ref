/*
 *  SPDX-License-Identifier: MIT
 */

#include "../aes.h"
#include "prg_tvs.hpp"

#include <boost/test/unit_test.hpp>
#include <array>

namespace {
  typedef std::array<uint8_t, 16> block_t;
} // namespace

using namespace prg_tvs;

BOOST_AUTO_TEST_SUITE(prgs)

BOOST_AUTO_TEST_CASE(test_increment_counter) {
  block_t iv = {
      0x8e, 0xa2, 0xb7, 0xca, 0x51, 0x67, 0x45, 0xbf,
      0xea, 0xfc, 0x49, 0x90, 0x4b, 0x49, 0x6f, 0xff,
  };
  {
    constexpr block_t iv_expected = {
        0x8f, 0xa2, 0xb7, 0xca, 0x51, 0x67, 0x45, 0xbf,
        0xea, 0xfc, 0x49, 0x90, 0x4b, 0x49, 0x6f, 0xff,
    };
    aes_increment_iv(iv.data());
    BOOST_TEST(iv == iv_expected);
  }

  iv = {
      0xff, 0xa2, 0xb7, 0xca, 0x51, 0x67, 0x45, 0xbf,
      0xea, 0xfc, 0x49, 0x90, 0x4b, 0x49, 0x6f, 0xff,
  };
  {
    constexpr block_t iv_expected = {
        0x00, 0xa3, 0xb7, 0xca, 0x51, 0x67, 0x45, 0xbf,
        0xea, 0xfc, 0x49, 0x90, 0x4b, 0x49, 0x6f, 0xff,
    };
    aes_increment_iv(iv.data());
    BOOST_TEST(iv == iv_expected);
  }
}

BOOST_AUTO_TEST_CASE(test_prg_128) {
  std::array<uint8_t, expected_128.size()> output_128;
  prg(key_128.data(), iv_128.data(), tweak_128, output_128.data(), 128, output_128.size());
  BOOST_TEST(output_128 == expected_128);
}

BOOST_AUTO_TEST_CASE(test_prg_192) {
  std::array<uint8_t, expected_192.size()> output_192;
  prg(key_192.data(), iv_192.data(), tweak_192, output_192.data(), 192, output_192.size());
  BOOST_TEST(output_192 == expected_192);
}

BOOST_AUTO_TEST_CASE(test_prg_256) {
  std::array<uint8_t, expected_256.size()> output_256;
  prg(key_256.data(), iv_256.data(), tweak_256, output_256.data(), 256, output_256.size());
  BOOST_TEST(output_256 == expected_256);
}

BOOST_AUTO_TEST_CASE(test_prg_128_2_lambda) {
  std::array<uint8_t, 2 * 128 / 8> expected_output, output;
  prg(key_128.data(), iv_128.data(), tweak_128, expected_output.data(), 128,
      expected_output.size());
  prg_2_lambda(key_128.data(), iv_128.data(), tweak_128, output.data(), 128);
  BOOST_TEST(expected_output == output);
}

BOOST_AUTO_TEST_CASE(test_prg_128_4_lambda) {
  std::array<uint8_t, 4 * 128 / 8> expected_output, output;
  prg(key_128.data(), iv_128.data(), tweak_128, expected_output.data(), 128,
      expected_output.size());
  prg_4_lambda(key_128.data(), iv_128.data(), tweak_128, output.data(), 128);
  BOOST_TEST(expected_output == output);
}

BOOST_AUTO_TEST_CASE(test_prg_192_2_lambda) {
  std::array<uint8_t, 2 * 192 / 8> expected_output, output;
  prg(key_192.data(), iv_192.data(), tweak_192, expected_output.data(), 192,
      expected_output.size());
  prg_2_lambda(key_192.data(), iv_192.data(), tweak_192, output.data(), 192);
  BOOST_TEST(expected_output == output);
}

BOOST_AUTO_TEST_CASE(test_prg_192_4_lambda) {
  std::array<uint8_t, 4 * 192 / 8> expected_output, output;
  prg(key_192.data(), iv_192.data(), tweak_192, expected_output.data(), 192,
      expected_output.size());
  prg_4_lambda(key_192.data(), iv_192.data(), tweak_192, output.data(), 192);
  BOOST_TEST(expected_output == output);
}

BOOST_AUTO_TEST_CASE(test_prg_256_2_lambda) {
  std::array<uint8_t, 2 * 256 / 8> expected_output, output;
  prg(key_256.data(), iv_256.data(), tweak_256, expected_output.data(), 256,
      expected_output.size());
  prg_2_lambda(key_256.data(), iv_256.data(), tweak_256, output.data(), 256);
  BOOST_TEST(expected_output == output);
}

BOOST_AUTO_TEST_CASE(test_prg_256_4_lambda) {
  std::array<uint8_t, 4 * 256 / 8> expected_output, output;
  prg(key_256.data(), iv_256.data(), tweak_256, expected_output.data(), 256,
      expected_output.size());
  prg_4_lambda(key_256.data(), iv_256.data(), tweak_256, output.data(), 256);
  BOOST_TEST(expected_output == output);
}

BOOST_AUTO_TEST_SUITE_END()
