#ifndef TESTS_UTILS_HPP
#define TESTS_UTILS_HPP

#include <array>
#include <boost/io/ios_state.hpp>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../hash_shake.h"

template <class T>
void print_array(const T* data, size_t len) {
  boost::io::ios_flags_saver ifs{std::cout};
  std::cout << "{";
  for (size_t s = 0; s != len; ++s) {
    std::cout << "0x" << std::setw(sizeof(T) * 2) << std::setfill('0') << std::hex
              << static_cast<unsigned int>(data[s]) << ", ";
  }
  std::cout << "};\n";
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

#endif