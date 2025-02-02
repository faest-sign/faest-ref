/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bavc.h"
#include "hash_shake.h"
#include "instances.hpp"
#include "utils.hpp"

#include <array>
#include <random>
#include <vector>

namespace {
  constexpr std::array<uint8_t, 16> iv{};
  constexpr std::array<uint8_t, 32> root_key{
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
      0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
      0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  };
} // namespace

int main() {
  std::mt19937_64 rd;

  std::cout << "#ifndef TESTS_BAVC_TVS_HPP\n";
  std::cout << "#define TESTS_BAVC_TVS_HPP\n\n";
  std::cout << "#include <array>\n";
  std::cout << "#include <cstdint>\n\n";
  std::cout << "namespace bavc_tvs {\n";

  for (const auto param_id : all_parameters) {
    const auto params       = *faest_get_paramset(param_id);
    const auto lambda       = params.faest_param.lambda;
    const auto lambda_bytes = lambda / 8;
    const auto com_size     = (faest_is_em(&params) ? 2 : 3) * lambda_bytes;

    bavc_t vc;
    bavc_commit(root_key.data(), iv.data(), &params, &vc);

    auto hashed_k  = hash_array(vc.k, (2 * params.faest_param.L - 1) * lambda_bytes);
    auto hashed_sd = hash_array(vc.sd, params.faest_param.L * lambda_bytes);

    std::cout << "namespace " << faest_get_param_name(param_id) << "{\n";
    print_named_array("h", "uint8_t", vc.h, 2 * lambda_bytes);
    print_named_array("hashed_k", "uint8_t", hashed_k);
    print_named_array("hashed_sd", "uint8_t", hashed_sd);

    std::vector<uint8_t> decom_i;
    std::vector<uint16_t> i_delta;
    i_delta.resize(params.faest_param.tau);

    bool ret = false;
    while (!ret) {
      for (unsigned int i = 0; i != params.faest_param.tau; ++i) {
        std::uniform_int_distribution<> distribution{
            0, static_cast<int>(
                   bavc_max_node_index(i, params.faest_param.tau1, params.faest_param.k)) -
                   1};
        i_delta[i] = distribution(rd);
      }

      decom_i.clear();
      decom_i.resize(com_size * params.faest_param.tau + params.faest_param.T_open * lambda_bytes);

      ret = bavc_open(&vc, i_delta.data(), decom_i.data(), &params);
    }

    auto hashed_decom_i = hash_array(decom_i);

    print_named_array("i_delta", "uint16_t", i_delta);
    print_named_array("hashed_decom_i", "uint8_t", hashed_decom_i);
    bavc_clear(&vc);

    std::vector<uint8_t> rec_h, rec_s;
    rec_h.resize(2 * lambda_bytes);
    rec_s.resize((params.faest_param.L - params.faest_param.tau) * lambda_bytes);

    bavc_rec_t vc_rec;
    vc_rec.h = rec_h.data();
    vc_rec.s = rec_s.data();

    bavc_reconstruct(decom_i.data(), i_delta.data(), iv.data(), &params, &vc_rec);

    auto hashed_rec_sd = hash_array(rec_s);
    print_named_array("hashed_rec_sd", "uint8_t", hashed_rec_sd);
    std::cout << "}\n" << std::endl;
  }

  std::cout << "}\n\n#endif" << std::endl;
}