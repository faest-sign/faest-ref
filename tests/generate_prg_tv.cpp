/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "aes.h"
#include "utils.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <random>
#include <vector>

int main() {
  std::mt19937_64 rd{0x5eed};
  std::uniform_int_distribution<uint8_t> distrib(0, 255);
  std::uniform_int_distribution<uint32_t> distrib_u32(0, std::numeric_limits<uint32_t>::max());

  std::cout << "#ifndef TESTS_PRG_TVS_HPP\n";
  std::cout << "#define TESTS_PRG_TVS_HPP\n\n";
  std::cout << "#include <array>\n";
  std::cout << "#include <cstdint>\n\n";
  std::cout << "namespace prg_tvs {\n";
  for (const auto lambda : {128, 192, 256}) {
    std::vector<uint8_t> key, output;
    std::array<uint8_t, 16> iv;

    key.resize(lambda / 8);
    output.resize(14 * 16 + 8);

    std::generate(key.begin(), key.end(), [&rd, &distrib] { return distrib(rd); });
    std::generate(iv.begin(), iv.end(), [&rd, &distrib] { return distrib(rd); });

    const uint32_t tweak = distrib_u32(rd);

    prg(key.data(), iv.data(), tweak, output.data(), lambda, output.size());

    std::cout << "constexpr std::array<uint8_t, " << key.size() << "> key_" << lambda;
    print_array(key.data(), key.size());
    std::cout << "constexpr std::array<uint8_t, " << iv.size() << "> iv_" << lambda;
    print_array(iv.data(), iv.size());
    std::cout << "constexpr uint32_t tweak_" << lambda << " = " << tweak << ";\n";
    std::cout << "constexpr std::array<uint8_t, " << output.size() << "> expected_" << lambda;
    print_array(output.data(), output.size());
  }
  std::cout << "}\n\n#endif" << std::endl;
}