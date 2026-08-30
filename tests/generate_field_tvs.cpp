/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bavc.h"
#include "fields.h"
#include "instances.hpp"
#include "utils.h"
#include "utils.hpp"

#include <array>
#include <random>
#include <vector>
#include <algorithm>

int main() {
  std::mt19937_64 rd{0x5eed};
  std::uniform_int_distribution<uint8_t> dist(0, 0xff);

  std::cout << "#ifndef TESTS_FIELD_TVS_HPP\n";
  std::cout << "#define TESTS_FIELD_TVS_HPP\n\n";
  std::cout << "#include <array>\n";
  std::cout << "#include <cstdint>\n\n";
  std::cout << "namespace field_tvs {\n";

  for (const auto param_id : all_parameters) {
    const auto& params               = *faest_get_paramset(param_id);
    const auto lambda                = params.lambda;
    const unsigned int lambda_bytes  = lambda / 8;
    const auto w_grind               = params.w_grind;
    const unsigned int w_grind_bytes = (w_grind + 7) / 8;
    const unsigned int deg           = bavc_max_node_depth(0, params.tau1, params.k);

    std::cout << "namespace " << faest_get_param_name(param_id) << "{\n";
    std::vector<uint8_t> lhs(w_grind_bytes, 0);
    std::vector<uint8_t> res(lambda_bytes, 0);

    std::generate(lhs.begin(), lhs.end(), [&rd, &dist] { return dist(rd); });
    for (size_t i = w_grind; i < w_grind_bytes * 8; ++i) {
      ptr_set_bit(lhs.data(), i, 0);
    }

    bf2_poly_mul(res.data(), lhs.data(), w_grind, params.M_TREE, lambda - w_grind + 1);
    print_named_array("lhs", "uint8_t", lhs);
    print_named_array("res", "uint8_t", res);

    std::vector<uint8_t> unreduced((deg * 2 + 7) / 8, 0);
    std::generate(unreduced.begin(), unreduced.end(), [&rd, &dist] { return dist(rd); });
    for (size_t i = deg * 2; i < unreduced.size() * 8; ++i) {
      ptr_set_bit(unreduced.data(), i, 0);
    }

    std::vector<uint8_t> reduced((deg * 2 + 7) / 8, 0);
    bf2_poly_reduce(reduced.data(), unreduced.data(), deg * 2, &params.TREE_MODULI[0], deg + 1);

    std::cout << "constexpr size_t deg = " << deg << ";\n";
    print_named_array("unreduced", "uint8_t", unreduced);
    print_named_array("reduced", "uint8_t", reduced);

    std::cout << "}\n" << std::endl;
  }

  std::cout << "}\n\n#endif" << std::endl;
}