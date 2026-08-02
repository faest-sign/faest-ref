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
using std::to_string;

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

    if (lambda == 256) {
      leaf_hash_256(output.data(), uhash.data(), x.data());
    }
    else if (lambda == 192) {
      leaf_hash_192(output.data(), uhash.data(), x.data());
    }
    else {
      leaf_hash_128(output.data(), uhash.data(), x.data());
    }

    
    const unsigned int ell = 32;
    const unsigned int d_zk = 3;
    const unsigned int lambda_bytes = lambda / 8;

    std::vector<uint8_t> uhash_sd, uhash_x, uhash_output;
    uhash_sd.resize(lambda_bytes * 6 + 8);
    uhash_x.resize((ell + (d_zk - 1) + 3) * lambda_bytes);
    uhash_output.resize(lambda_bytes * 3, 0);

    std::generate(uhash_sd.begin(), uhash_sd.end(), [&rd, &distrib] { return distrib(rd); });
    std::generate(uhash_x.begin(), uhash_x.end(), [&rd, &distrib] { return distrib(rd); });

    if (lambda == 256) {
      vole_hash_new_256(uhash_output.data(), uhash_sd.data(), uhash_x.data(), ell, d_zk);
    }
    else if (lambda == 192) {
      vole_hash_new_192(uhash_output.data(), uhash_sd.data(), uhash_x.data(), ell, d_zk);
    }
    else {
      vole_hash_new_128(uhash_output.data(), uhash_sd.data(), uhash_x.data(), ell, d_zk);
    }

    std::cout << "namespace leaf_hash_" << lambda << " {";
    print_named_array("uhash", "uint8_t", uhash);
    print_named_array("x", "uint8_t", x);
    print_named_array("expected_h", "uint8_t", output);
    std::cout << "}\n";

    std::cout << "namespace vole_hash_" << lambda << " {";
    print_named_array(("vole_hash_new_" + to_string(lambda) + "_sd").c_str(), "uint8_t", uhash_sd);
    print_named_array(("vole_hash_new_" + to_string(lambda) + "_xs").c_str(), "uint8_t", uhash_x);
    print_named_array(("vole_hash_new_" + to_string(lambda) + "_digest").c_str(), "uint8_t", uhash_output);
    std::cout << "}\n";


    
  }
  std::cout << "}\n\n#endif" << std::endl;
}