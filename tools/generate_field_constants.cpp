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
    std::cout << boost::format("UINT64_C(0x%08x)") % v;
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
} // namespace

int main() {
  {
    bf128_t beta = bf128_add(bf128_get_alpha(5), bf128_get_alpha(3));

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
    bf192_t beta = bf192_add(bf192_get_alpha(5), bf192_get_alpha(3));

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
    bf256_t beta = bf256_add(bf256_get_alpha(5), bf256_get_alpha(3));

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
}