/*
 *  SPDX-License-Identifier: MIT
 */

// Reference - https://nvlpubs.nist.gov/nistpubs/fips/nist.fips.197.pdf
// Tested against Appendix C.1

#include "../aes.h"

#include <boost/test/unit_test.hpp>
#include <array>
#include <cstring>
#include <vector>

#include "faest_tvs.hpp"
#include "utils.hpp"

namespace {
  typedef std::array<uint8_t, 16> block_t;
  typedef std::array<uint8_t, 24> block192_t;
  typedef std::array<uint8_t, 32> block256_t;
} // namespace

BOOST_AUTO_TEST_SUITE(aes)

BOOST_AUTO_TEST_CASE(test_aes128) {
  constexpr uint8_t key_128[16] = {
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
      0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
  };
  constexpr uint8_t plaintext_128[16] = {
      0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
      0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff,
  };
  constexpr block_t expected_128 = {
      0x69, 0xc4, 0xe0, 0xd8, 0x6a, 0x7b, 0x04, 0x30,
      0xd8, 0xcd, 0xb7, 0x80, 0x70, 0xb4, 0xc5, 0x5a,
  };

  aes_round_keys_t ctx;
  aes128_init_round_keys(&ctx, key_128);

  block_t output_128;
  aes128_encrypt_block(&ctx, plaintext_128, output_128.data());

  BOOST_TEST(output_128 == expected_128);
}

BOOST_AUTO_TEST_CASE(test_aes192) {
  constexpr uint8_t key_192[24] = {
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b,
      0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
  };
  constexpr uint8_t plaintext_192[16] = {
      0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
      0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff,
  };
  constexpr block_t expected_192 = {
      0xdd, 0xa9, 0x7c, 0xa4, 0x86, 0x4c, 0xdf, 0xe0,
      0x6e, 0xaf, 0x70, 0xa0, 0xec, 0x0d, 0x71, 0x91,
  };

  aes_round_keys_t ctx;
  aes192_init_round_keys(&ctx, key_192);

  block_t output_192;
  aes192_encrypt_block(&ctx, plaintext_192, output_192.data());

  BOOST_TEST(output_192 == expected_192);
}

BOOST_AUTO_TEST_CASE(test_aes256) {
  constexpr uint8_t key_256[32] = {
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
      0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
      0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  };
  constexpr uint8_t plaintext_256[16] = {
      0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
      0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff,
  };
  constexpr block_t expected_256 = {
      0x8e, 0xa2, 0xb7, 0xca, 0x51, 0x67, 0x45, 0xbf,
      0xea, 0xfc, 0x49, 0x90, 0x4b, 0x49, 0x60, 0x89,
  };

  aes_round_keys_t ctx;
  aes256_init_round_keys(&ctx, key_256);

  block_t output_256;
  aes256_encrypt_block(&ctx, plaintext_256, output_256.data());

  BOOST_TEST(output_256 == expected_256);
}

BOOST_AUTO_TEST_CASE(test_rijndael192) {
  constexpr block192_t key_192       = {0x80, 0x00};
  constexpr block192_t expected_192  = {0x56, 0x4d, 0x36, 0xfd, 0xeb, 0x8b, 0xf7, 0xe2,
                                        0x75, 0xf0, 0x10, 0xb2, 0xf5, 0xee, 0x69, 0xcf,
                                        0xea, 0xe6, 0x7e, 0xa0, 0xe3, 0x7e, 0x32, 0x09};
  constexpr block192_t plaintext_192 = {0};

  aes_round_keys_t ctx;
  rijndael192_init_round_keys(&ctx, key_192.data());

  block192_t output_192;
  rijndael192_encrypt_block(&ctx, plaintext_192.data(), output_192.data());

  BOOST_TEST(output_192 == expected_192);
}

BOOST_AUTO_TEST_CASE(test_rijndael256) {
  constexpr block256_t key_256      = {0x80, 0x00};
  constexpr block256_t expected_256 = {
      0xE6, 0x2A, 0xBC, 0xE0, 0x69, 0x83, 0x7B, 0x65, 0x30, 0x9B, 0xE4,
      0xED, 0xA2, 0xC0, 0xE1, 0x49, 0xFE, 0x56, 0xC0, 0x7B, 0x70, 0x82,
      0xD3, 0x28, 0x7F, 0x59, 0x2C, 0x4A, 0x49, 0x27, 0xA2, 0x77,
  };
  constexpr block256_t plaintext_256 = {0};

  aes_round_keys_t ctx;
  rijndael256_init_round_keys(&ctx, key_256.data());

  block256_t output_256;
  rijndael256_encrypt_block(&ctx, plaintext_256.data(), output_256.data());

  BOOST_TEST(output_256 == expected_256);
}

BOOST_AUTO_TEST_CASE(test_extend_witness_aes128) {
  namespace tv = faest_tvs::faest_128s;
  auto params  = faest_get_paramset(FAEST_128S);

  std::vector<uint8_t> witness;
  witness.resize(params->l / 8);

  aes_extend_witness(witness.data(), tv::sk.data() + params->owf_input_size, tv::sk.data(), params);

  std::vector<uint8_t> expected_witness{tv::witness.begin(), tv::witness.end()};
  BOOST_TEST(expected_witness == witness);
}

BOOST_AUTO_TEST_CASE(test_extend_witness_aes192) {
  namespace tv = faest_tvs::faest_192s;
  auto params  = faest_get_paramset(FAEST_192S);

  std::vector<uint8_t> witness;
  witness.resize(params->l / 8);

  aes_extend_witness(witness.data(), tv::sk.data() + params->owf_input_size, tv::sk.data(), params);

  std::vector<uint8_t> expected_witness{tv::witness.begin(), tv::witness.end()};
  BOOST_TEST(expected_witness == witness);
}

BOOST_AUTO_TEST_CASE(test_extend_witness_aes256) {
  namespace tv = faest_tvs::faest_256s;
  auto params  = faest_get_paramset(FAEST_256S);

  std::vector<uint8_t> witness;
  witness.resize(params->l / 8);

  aes_extend_witness(witness.data(), tv::sk.data() + params->owf_input_size, tv::sk.data(), params);

  std::vector<uint8_t> expected_witness{tv::witness.begin(), tv::witness.end()};
  BOOST_TEST(expected_witness == witness);
}

BOOST_AUTO_TEST_CASE(test_extend_witness_em128) {
  namespace tv = faest_tvs::faest_em_128s;
  auto params  = faest_get_paramset(FAEST_EM_128S);

  std::vector<uint8_t> witness;
  witness.resize(params->l / 8);

  aes_extend_witness(witness.data(), tv::sk.data() + params->owf_input_size, tv::sk.data(), params);

  std::vector<uint8_t> expected_witness{tv::witness.begin(), tv::witness.end()};
  BOOST_TEST(expected_witness == witness);
}

BOOST_AUTO_TEST_CASE(test_extend_witness_em192) {
  namespace tv = faest_tvs::faest_em_192s;
  auto params  = faest_get_paramset(FAEST_EM_192S);

  std::vector<uint8_t> witness;
  witness.resize(params->l / 8);

  aes_extend_witness(witness.data(), tv::sk.data() + params->owf_input_size, tv::sk.data(), params);

  std::vector<uint8_t> expected_witness{tv::witness.begin(), tv::witness.end()};
  BOOST_TEST(expected_witness == witness);
}

BOOST_AUTO_TEST_CASE(test_extend_witness_em256) {
  namespace tv = faest_tvs::faest_em_256s;
  auto params  = faest_get_paramset(FAEST_EM_256S);

  std::vector<uint8_t> witness;
  witness.resize(params->l / 8);

  aes_extend_witness(witness.data(), tv::sk.data() + params->owf_input_size, tv::sk.data(), params);

  std::vector<uint8_t> expected_witness{tv::witness.begin(), tv::witness.end()};
  BOOST_TEST(expected_witness == witness);
}

BOOST_AUTO_TEST_CASE(test_invnorm) {
  BOOST_TEST(invnorm(0x00) == 0);
  BOOST_TEST(invnorm(0x01) == 1);
  BOOST_TEST(invnorm(0x02) == (1 << 3 | 1 << 2 | 1));
  BOOST_TEST(invnorm(0x80) == (1 << 3 | 1 << 2 | 1));
  BOOST_TEST(invnorm(0x88) == (1 << 2 | 1));
}

BOOST_AUTO_TEST_SUITE_END()
