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
#include <limits>
#include <algorithm>

namespace {
  constexpr std::array<uint8_t, 16> iv{};
  constexpr std::array<uint8_t, 32> root_key{
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
      0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
      0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  };
} // namespace

int main() {
  std::mt19937_64 rd{0x5eed};
  std::uniform_int_distribution<uint8_t> dist(0, 0xff);
  std::uniform_int_distribution<uint32_t> dist_u32(0, std::numeric_limits<uint32_t>::max());

  std::cout << "#ifndef TESTS_BAVC_TVS_HPP\n";
  std::cout << "#define TESTS_BAVC_TVS_HPP\n\n";
  std::cout << "#include <array>\n";
  std::cout << "#include <cstdint>\n\n";
  std::cout << "namespace bavc_tvs {\n";

  for (const auto param_id : all_parameters) {
    const auto params       = *faest_get_paramset(param_id);
    const auto lambda       = params.lambda;
    const auto lambda_bytes = lambda / 8;
    const auto com_size     = (faest_is_em(&params) ? 2 : 3) * lambda_bytes;

    bavc_t vc;
    bavc_commit(&vc, root_key.data(), iv.data(), &params);

    auto hashed_k  = hash_array(vc.k, (2 * params.L - 1) * lambda_bytes);
    auto hashed_sd = hash_array(vc.sd, params.L * lambda_bytes);

    std::cout << "namespace " << faest_get_param_name(param_id) << "{\n";
    print_named_array("h", "uint8_t", vc.h, 2 * lambda_bytes);
    print_named_array("hashed_k", "uint8_t", hashed_k);
    print_named_array("hashed_sd", "uint8_t", hashed_sd);

    std::vector<uint8_t> decom_i;
    std::vector<uint16_t> i_delta;
    i_delta.resize(params.tau);

    bool ret = false;
    while (!ret) {
      for (unsigned int i = 0; i != params.tau; ++i) {
        std::uniform_int_distribution<> distribution{
            0, static_cast<int>(bavc_max_node_index(i, params.tau1, params.k)) - 1};
        i_delta[i] = distribution(rd);
      }

      decom_i.clear();
      decom_i.resize(com_size * params.tau + params.T_open * lambda_bytes);

      ret = bavc_open(decom_i.data(), &vc, i_delta.data(), &params);
    }

    auto hashed_decom_i = hash_array(decom_i);

    print_named_array("i_delta", "uint16_t", i_delta);
    print_named_array("hashed_decom_i", "uint8_t", hashed_decom_i);
    bavc_clear(&vc);

    std::vector<uint8_t> rec_h, rec_s;
    rec_h.resize(2 * lambda_bytes);
    rec_s.resize((params.L - params.tau) * lambda_bytes);

    bavc_rec_t vc_rec;
    vc_rec.h = rec_h.data();
    vc_rec.s = rec_s.data();

    bavc_reconstruct(&vc_rec, decom_i.data(), i_delta.data(), iv.data(), &params);

    auto hashed_rec_sd = hash_array(rec_s);
    print_named_array("hashed_rec_sd", "uint8_t", hashed_rec_sd);

    std::vector<uint8_t> sd, com, key, iv, uhash;
    sd.resize(lambda_bytes);
    com.resize(com_size);
    key.resize(lambda_bytes);
    iv.resize(IV_SIZE);
    if (!faest_is_em(&params)) {
      uhash.resize(3 * lambda_bytes);
    }

    std::generate(key.begin(), key.end(), [&rd, &dist] { return dist(rd); });
    std::generate(iv.begin(), iv.end(), [&rd, &dist] { return dist(rd); });
    std::generate(uhash.begin(), uhash.end(), [&rd, &dist] { return dist(rd); });

    const uint32_t tweak = dist_u32(rd);
    leaf_commit(sd.data(), com.data(), key.data(), iv.data(), tweak, uhash.data(), &params);

    print_named_array("leaf_commit_key", "uint8_t", key);
    print_named_array("leaf_commit_iv", "uint8_t", iv);
    if (!faest_is_em(&params)) {
      print_named_array("leaf_commit_uhash", "uint8_t", uhash);
    }
    std::cout << "constexpr uint32_t tweak = " << tweak << ";\n";
    print_named_array("leaf_commit_expected_sd", "uint8_t", sd);
    print_named_array("leaf_commit_expected_com", "uint8_t", com);

    std::cout << "}\n" << std::endl;
  }

  std::cout << "}\n\n#endif" << std::endl;
}