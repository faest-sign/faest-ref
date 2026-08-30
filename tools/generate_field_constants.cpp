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
  void print_u64(uint64_t v) {
    std::cout << boost::format("UINT64_C(0x%016x)") % v;
  }

  void print_bf128(bf128_t v) {
    std::cout << "BF128C(";
    print_u64(BF_VALUE(v, 0));
    std::cout << ", ";
    print_u64(BF_VALUE(v, 1));
    std::cout << ")";
  }

  void print_bf192(bf192_t v) {
    std::cout << "BF192C(";
    print_u64(BF_VALUE(v, 0));
    std::cout << ", ";
    print_u64(BF_VALUE(v, 1));
    std::cout << ", ";
    print_u64(BF_VALUE(v, 2));
    std::cout << ")";
  }

  void print_bf256(bf256_t v) {
    std::cout << "BF256C(";
    print_u64(BF_VALUE(v, 0));
    std::cout << ", ";
    print_u64(BF_VALUE(v, 1));
    std::cout << ", ";
    print_u64(BF_VALUE(v, 2));
    std::cout << ", ";
    print_u64(BF_VALUE(v, 3));
    std::cout << ")";
  }

  constexpr uint8_t sbox_x[9] = {0x05, 0x09, 0xf9, 0x25, 0xf4, 0x01, 0xb5, 0x8f, 0x63};
} // namespace

int main() {
  {
    std::array<bf128_t, 9> c;
    for (unsigned int i = 0; i < 9; i++) {
      bf128_byte_combine_bits(&c[i], sbox_x[i]);
    }

    std::cout << "static const bf128_t bf128_c[9] = {\n";
    for (auto v : c) {
      print_bf128(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
  }

  {
    std::array<bf192_t, 9> c;
    for (unsigned int i = 0; i < 9; i++) {
      bf192_byte_combine_bits(&c[i], sbox_x[i]);
    }

    std::cout << "static const bf192_t bf192_c[9] = {\n";
    for (auto v : c) {
      print_bf192(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
  }

  {
    std::array<bf256_t, 9> c;
    for (unsigned int i = 0; i < 9; i++) {
      bf256_byte_combine_bits(&c[i], sbox_x[i]);
    }

    std::cout << "static const bf256_t bf256_c[9] = {\n";
    for (auto v : c) {
      print_bf256(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
  }

  {
    bf128_t v2;
    bf128_byte_combine_bits(&v2, 2);
    bf128_t v3;
    bf128_byte_combine_bits(&v3, 3);

    std::cout << "static const bf128_t bf128_bc_2 = ";
    print_bf128(v2);
    std::cout << ";\n";

    std::cout << "static const bf128_t bf128_bc_3 = ";
    print_bf128(v3);
    std::cout << ";\n";
  }

  {
    bf192_t v2;
    bf192_byte_combine_bits(&v2, 2);
    bf192_t v3;
    bf192_byte_combine_bits(&v3, 3);

    std::cout << "static const bf192_t bf192_bc_2 = ";
    print_bf192(v2);
    std::cout << ";\n";

    std::cout << "static const bf192_t bf192_bc_3 = ";
    print_bf192(v3);
    std::cout << ";\n";
  }

  {
    bf256_t v2;
    bf256_byte_combine_bits(&v2, 2);
    bf256_t v3;
    bf256_byte_combine_bits(&v3, 3);

    std::cout << "static const bf256_t bf256_bc_2 = ";
    print_bf256(v2);
    std::cout << ";\n";

    std::cout << "static const bf256_t bf256_bc_3 = ";
    print_bf256(v3);
    std::cout << ";\n";
  }
}