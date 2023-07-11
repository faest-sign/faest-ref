/*
 *  SPDX-License-Identifier: MIT
 */

#include <array>

extern "C" {
#include "api.h"
}

#include <catch_amalgamated.hpp>

TEST_CASE("bench keygen", "[.][bench]") {
  std::array<unsigned char, CRYPTO_SECRETKEYBYTES> sk;
  std::array<unsigned char, CRYPTO_PUBLICKEYBYTES> pk;

  BENCHMARK("keygen") {
    return crypto_sign_keypair(pk.data(), sk.data());
  };
}

TEST_CASE("bench sign", "[.][bench]") {
  std::array<unsigned char, CRYPTO_SECRETKEYBYTES> sk;
  std::array<unsigned char, CRYPTO_PUBLICKEYBYTES> pk;
  crypto_sign_keypair(pk.data(), sk.data());
  const std::string message =
      "This document describes and specifies the FAEST digital signature algorithm.";
  std::vector<unsigned char> signed_message(CRYPTO_BYTES + message.size());
  unsigned long long signed_message_len = CRYPTO_BYTES + message.size();

  BENCHMARK("sign") {
    return crypto_sign(signed_message.data(), &signed_message_len,
                       reinterpret_cast<const unsigned char*>(message.data()), message.size(),
                       sk.data());
  };

  REQUIRE(signed_message_len == signed_message.size());
}

TEST_CASE("bench verify", "[.][bench]") {
  std::array<unsigned char, CRYPTO_SECRETKEYBYTES> sk;
  std::array<unsigned char, CRYPTO_PUBLICKEYBYTES> pk;
  crypto_sign_keypair(pk.data(), sk.data());
  const std::string message =
      "This document describes and specifies the FAEST digital signature algorithm.";
  std::vector<unsigned char> signed_message(CRYPTO_BYTES + message.size());
  unsigned long long signed_message_len = CRYPTO_BYTES + message.size();
  crypto_sign(signed_message.data(), &signed_message_len,
              reinterpret_cast<const unsigned char*>(message.data()), message.size(), sk.data());
  std::vector<unsigned char> opened_message(message.size());
  unsigned long long opened_message_len = message.size();

  BENCHMARK("verify") {
    return crypto_sign_open(opened_message.data(), &opened_message_len, signed_message.data(),
                            signed_message_len, pk.data());
  };

  REQUIRE(opened_message_len == opened_message.size());
  REQUIRE(opened_message ==
          std::vector<unsigned char>(reinterpret_cast<const unsigned char*>(message.c_str()),
                                     reinterpret_cast<const unsigned char*>(message.c_str()) +
                                         message.size()));
}
