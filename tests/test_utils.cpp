/*
 *  SPDX-License-Identifier: MIT
 */

#include "utils.hpp"
#include "../hash_shake.h"

#include <random>
#include <algorithm>

void rand_bytes(std::vector<uint8_t>& v) {
  std::random_device r;
  std::default_random_engine rand_engine(r());
  std::uniform_int_distribution<uint8_t> dist_byte(0, 0xff);

  std::generate(v.begin(), v.end(), [&rand_engine, &dist_byte] { return dist_byte(rand_engine); });
}

std::array<uint8_t, 64> hash_array(const uint8_t* data, size_t len) {
  hash_context ctx;
  hash_init(&ctx, 256);
  hash_update(&ctx, data, len);
  hash_final(&ctx);

  std::array<uint8_t, 64> ret;
  hash_squeeze(&ctx, ret.data(), ret.size());
  return ret;
}

std::array<uint8_t, 64> hash_array(const std::vector<uint8_t>& data) {
  return hash_array(data.data(), data.size());
}