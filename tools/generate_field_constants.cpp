/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "fields.h"

#include <array>
#include <boost/format.hpp>
#include <iostream>

namespace {
  template <size_t Bytes, class Field>
  void print_field(const Field& value, void (*store)(uint8_t*, const Field*)) {
    std::array<uint8_t, Bytes> bytes;
    store(bytes.data(), &value);
    std::cout << "{";
    for (size_t i = 0; i != bytes.size(); ++i) {
      if (i) {
        std::cout << ", ";
      }
      std::cout << boost::format("0x%02x") % static_cast<unsigned int>(bytes[i]);
    }
    std::cout << "}";
  }

  constexpr uint8_t sbox_x[9] = {0x05, 0x09, 0xf9, 0x25, 0xf4, 0x01, 0xb5, 0x8f, 0x63};

  template <size_t Bytes, class Field>
  void print_constants(unsigned int bits, void (*combine)(Field*, uint8_t),
                       void (*store)(uint8_t*, const Field*)) {
    std::cout << "static const uint8_t bf" << bits << "_c_bytes[9][BF" << bits
              << "_NUM_BYTES] = {\n";
    for (uint8_t x : sbox_x) {
      Field value;
      combine(&value, x);
      std::cout << "    ";
      print_field<Bytes>(value, store);
      std::cout << ",\n";
    }
    std::cout << "};\n";

    Field value;
    combine(&value, 2);
    std::cout << "static const uint8_t bf" << bits << "_bc_2_bytes[BF" << bits << "_NUM_BYTES] = ";
    print_field<Bytes>(value, store);
    std::cout << ";\n";
  }
} // namespace

int main() {
  print_constants<BF128_NUM_BYTES>(128, bf128_byte_combine_bits, bf128_store);
  print_constants<BF192_NUM_BYTES>(192, bf192_byte_combine_bits, bf192_store);
  print_constants<BF256_NUM_BYTES>(256, bf256_byte_combine_bits, bf256_store);
}
