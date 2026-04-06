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

  constexpr uint8_t sbox_x[9]    = {0x05, 0x09, 0xf9, 0x25, 0xf4, 0x01, 0xb5, 0x8f, 0x63};
  constexpr uint8_t sbox_x_sq[9] = {0x11, 0x41, 0x07, 0x7d, 0x56, 0x01, 0xfc, 0xcf, 0xc2};
} // namespace

int main() {
  {
    bf128_t beta      = bf128_get_alpha(5);
    const bf128_t tmp = bf128_get_alpha(3);
    bf128_add_inplace(&beta, &tmp);

    std::array<bf128_t, 5> squares;
    std::array<bf128_t, 4> cubes;

    squares[0] = beta;
    for (size_t idx = 1; idx != squares.size(); ++idx) {
      squares[idx] = bf128_mul(squares[idx - 1], squares[idx - 1]);
    }
    cubes[0] = bf128_mul(squares[1], beta);
    for (size_t idx = 1; idx != cubes.size(); ++idx) {
      cubes[idx] = bf128_mul(cubes[idx - 1], cubes[idx - 1]);
    }

    std::cout << "static const bf128_t bf128_beta_squares[5] = {\n";
    for (auto v : squares) {
      print_bf128(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
    std::cout << "static const bf128_t bf128_beta_cubes[4] = {\n";
    for (auto v : cubes) {
      print_bf128(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
  }

  {
    bf192_t beta      = bf192_get_alpha(5);
    const bf192_t tmp = bf192_get_alpha(3);
    bf192_add_inplace(&beta, &tmp);

    std::array<bf192_t, 5> squares;
    std::array<bf192_t, 4> cubes;

    squares[0] = beta;
    for (size_t idx = 1; idx != squares.size(); ++idx) {
      squares[idx] = bf192_mul(squares[idx - 1], squares[idx - 1]);
    }
    cubes[0] = bf192_mul(squares[1], beta);
    for (size_t idx = 1; idx != cubes.size(); ++idx) {
      cubes[idx] = bf192_mul(cubes[idx - 1], cubes[idx - 1]);
    }

    std::cout << "static const bf192_t bf192_beta_squares[5] = {\n";
    for (auto v : squares) {
      print_bf192(v);
      std::cout << ", \n";
    }
    std::cout << "\n};\n";
    std::cout << "static const bf192_t bf192_beta_cubes[4] = {\n";
    for (auto v : cubes) {
      print_bf192(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
  }

  {
    bf256_t beta      = bf256_get_alpha(5);
    const bf256_t tmp = bf256_get_alpha(3);
    bf256_add_inplace(&beta, &tmp);

    std::array<bf256_t, 5> squares;
    std::array<bf256_t, 4> cubes;

    squares[0] = beta;
    for (size_t idx = 1; idx != squares.size(); ++idx) {
      squares[idx] = bf256_mul(squares[idx - 1], squares[idx - 1]);
    }
    cubes[0] = bf256_mul(squares[1], beta);
    for (size_t idx = 1; idx != cubes.size(); ++idx) {
      cubes[idx] = bf256_mul(cubes[idx - 1], cubes[idx - 1]);
    }

    std::cout << "static const bf256_t bf256_beta_squares[5] = {\n";
    for (auto v : squares) {
      print_bf256(v);
      std::cout << ", \n";
    }
    std::cout << "\n};\n";
    std::cout << "static const bf256_t bf256_beta_cubes[4] = {\n";
    for (auto v : cubes) {
      print_bf256(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
  }

  {
    std::array<bf128_t, 9> c_squares, c;
    for (unsigned i = 0; i < 9; i++) {
      bf128_byte_combine_bits(&c[i], sbox_x[i]);
      bf128_byte_combine_bits(&c_squares[i], sbox_x_sq[i]);
    }

    std::cout << "static const bf128_t bf128_c[9] = {\n";
    for (auto v : c) {
      print_bf128(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
    std::cout << "static const bf128_t bf128_c_squares[9] = {\n";
    for (auto v : c_squares) {
      print_bf128(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
  }

  {
    std::array<bf192_t, 9> c_squares, c;
    for (unsigned i = 0; i < 9; i++) {
      bf192_byte_combine_bits(&c[i], sbox_x[i]);
      bf192_byte_combine_bits(&c_squares[i], sbox_x_sq[i]);
    }

    std::cout << "static const bf192_t bf192_c[9] = {\n";
    for (auto v : c) {
      print_bf192(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
    std::cout << "static const bf192_t bf192_c_squares[9] = {\n";
    for (auto v : c_squares) {
      print_bf192(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
  }

  {
    std::array<bf256_t, 9> c_squares, c;
    for (unsigned i = 0; i < 9; i++) {
      bf256_byte_combine_bits(&c[i], sbox_x[i]);
      bf256_byte_combine_bits(&c_squares[i], sbox_x_sq[i]);
    }

    std::cout << "static const bf256_t bf256_c[9] = {\n";
    for (auto v : c) {
      print_bf256(v);
      std::cout << ", \n";
    }
    std::cout << "};\n";
    std::cout << "static const bf256_t bf256_c_squares[9] = {\n";
    for (auto v : c_squares) {
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

    v2 = bf128_mul(v2, v2);
    v3 = bf128_mul(v3, v3);

    std::cout << "static const bf128_t bf128_bc_2_sq = ";
    print_bf128(v2);
    std::cout << ";\n";

    std::cout << "static const bf128_t bf128_bc_3_sq = ";
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

    v2 = bf192_mul(v2, v2);
    v3 = bf192_mul(v3, v3);

    std::cout << "static const bf192_t bf192_bc_2_sq = ";
    print_bf192(v2);
    std::cout << ";\n";

    std::cout << "static const bf192_t bf192_bc_3_sq = ";
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

    v2 = bf256_mul(v2, v2);
    v3 = bf256_mul(v3, v3);

    std::cout << "static const bf256_t bf256_bc_2_sq = ";
    print_bf256(v2);
    std::cout << ";\n";

    std::cout << "static const bf256_t bf256_bc_3_sq = ";
    print_bf256(v3);
    std::cout << ";\n";
  }
}