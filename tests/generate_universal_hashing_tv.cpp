/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "universal_hashing.h"
#include "utils.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <random>
#include <vector>

int main() {
  std::mt19937_64 rd{0x5eed};
  std::uniform_int_distribution<uint8_t> distrib(0, 255);

  std::cout << "#ifndef TESTS_UNIVERSAL_HASHING_TVS_HPP\n";
  std::cout << "#define TESTS_UNIVERSAL_HASHING_TVS_HPP\n\n";
  std::cout << "#include <array>\n";
  std::cout << "#include <cstdint>\n\n";
  std::cout << "namespace universal_hashing_tvs {\n";
  for (const auto lambda : {128, 192, 256}) {
    std::vector<uint8_t> uhash, x, output;
    uhash.resize(lambda / 8 * 3);
    x.resize(lambda / 8 * 4);
    output.resize(lambda / 8 * 3, 0);

    std::generate(uhash.begin(), uhash.end(), [&rd, &distrib] { return distrib(rd); });
    std::generate(x.begin(), x.end(), [&rd, &distrib] { return distrib(rd); });

    leaf_hash(output.data(), uhash.data(), x.data(), lambda);

    std::cout << "namespace leaf_hash_" << lambda << " {";
    print_named_array("uhash", "uint8_t", uhash);
    print_named_array("x", "uint8_t", x);
    print_named_array("expected_h", "uint8_t", output);
    std::cout << "}\n";
  }
  std::cout << "}\n\n#endif" << std::endl;
}