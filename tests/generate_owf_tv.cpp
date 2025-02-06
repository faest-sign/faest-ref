/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "owf.h"
#include "utils.hpp"
#include "instances.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <random>
#include <vector>

int main() {
  std::mt19937_64 rd{0x5eed};
  std::uniform_int_distribution<uint8_t> distrib(0, 255);

  std::cout << "#ifndef TESTS_OWF_TVS_HPP\n";
  std::cout << "#define TESTS_OWF_TVS_HPP\n\n";
  std::cout << "#include <array>\n";
  std::cout << "#include <cstdint>\n\n";
  std::cout << "namespace owf_tvs {\n";
  for (const auto param_id : all_parameters) {
    switch (param_id) {
    case FAEST_128S:
    case FAEST_192S:
    case FAEST_256S:
    case FAEST_EM_128S:
    case FAEST_EM_192S:
    case FAEST_EM_256S:
      continue;
    default:
      break;
    }

    const auto params = *faest_get_paramset(param_id);
    std::vector<uint8_t> owf_key, owf_input, owf_output;
    owf_key.resize(params.lambda / 8);
    owf_input.resize(params.owf_input_size);
    owf_output.resize(params.owf_output_size);

    std::generate(owf_key.begin(), owf_key.end(), [&rd, &distrib] { return distrib(rd); });
    std::generate(owf_input.begin(), owf_input.end(), [&rd, &distrib] { return distrib(rd); });

    switch (param_id) {
    case FAEST_128F:
      owf_128(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    case FAEST_192F:
      owf_192(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    case FAEST_256F:
      owf_256(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    case FAEST_EM_128F:
      owf_em_128(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    case FAEST_EM_192F:
      owf_em_192(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    case FAEST_EM_256F:
      owf_em_256(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    default:
      continue;
    }

    std::cout << "namespace " << faest_get_param_name(param_id) << "{\n";
    print_named_array("key", "uint8_t", owf_key);
    print_named_array("input", "uint8_t", owf_input);
    print_named_array("output", "uint8_t", owf_output);
    std::cout << "}\n";
  }
  std::cout << "}\n\n#endif" << std::endl;
}