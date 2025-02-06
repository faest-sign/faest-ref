/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vole.h"
#include "hash_shake.h"
#include "instances.hpp"
#include "utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace {
  constexpr std::array<uint8_t, 16> iv{};
  constexpr std::array<uint8_t, 32> root_key{
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
      0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
      0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  };
} // namespace

int main() {
  std::mt19937_64 mt{0x5eed};
  std::uniform_int_distribution<uint8_t> dist(0, 0xff);

  std::cout << "#ifndef TEST_VOLE_TVS_HPP\n";
  std::cout << "#define TEST_VOLE_TVS_HPP\n\n";
  std::cout << "#include <array>\n";
  std::cout << "#include <cstdint>\n\n";
  std::cout << "namespace vole_tvs {\n";
  for (const auto param_id : all_parameters) {
    const auto params                = *faest_get_paramset(param_id);
    const unsigned int lambda        = params.lambda;
    const unsigned int lambda_bytes  = lambda / 8;
    const unsigned int ell_hat       = params.l + params.lambda * 3 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const auto com_size              = (faest_is_em(&params) ? 2 : 3) * lambda_bytes;

    bavc_t bavc_com;

    std::vector<uint8_t> chal, c, decom_i, u, q_storage, v_storage;
    chal.resize(lambda_bytes);
    c.resize((params.tau - 1) * ell_hat_bytes);
    decom_i.resize(com_size * params.tau + params.T_open * lambda_bytes);
    u.resize(ell_hat_bytes);

    std::vector<uint8_t*> q, v;
    q.resize(lambda);
    v.resize(lambda);

    q_storage.resize(lambda * ell_hat_bytes);
    v_storage.resize(lambda * ell_hat_bytes);

    q[0] = q_storage.data();
    v[0] = v_storage.data();
    for (unsigned int i = 1; i < lambda; ++i) {
      q[i] = q[0] + i * ell_hat_bytes;
      v[i] = v[0] + i * ell_hat_bytes;
    }

    std::cout << "namespace " << faest_get_param_name(param_id) << "{\n";

    vole_commit(root_key.data(), iv.data(), ell_hat, &params, &bavc_com, c.data(), u.data(),
                v.data());
    print_named_array("h", "uint8_t", bavc_com.h, 2 * lambda_bytes);
    print_named_array("hashed_c", "uint8_t", hash_array(c));
    print_named_array("hashed_u", "uint8_t", hash_array(u));
    print_named_array("hashed_v", "uint8_t", hash_array(v_storage));

    while (true) {
      std::generate(chal.begin(), chal.end(), [&mt, &dist] { return dist(mt); });
      for (unsigned int i = lambda - params.w_grind; i != lambda; ++i) {
        ptr_set_bit(chal.data(), i, 0);
      }

      uint16_t i_delta[MAX_TAU];
      decode_all_chall_3(i_delta, chal.data(), &params);
      if (!bavc_open(decom_i.data(), &bavc_com, i_delta, &params)) {
        continue;
      }
      print_named_array("chall", "uint8_t", chal);

      std::vector<uint8_t> hcom_rec;
      hcom_rec.resize(lambda_bytes * 2);
      vole_reconstruct(hcom_rec.data(), q.data(), iv.data(), chal.data(), decom_i.data(), c.data(),
                       ell_hat, &params);
      print_named_array("hashed_q", "uint8_t", hash_array(q_storage));
      break;
    }
    bavc_clear(&bavc_com);
    std::cout << "}\n";
  }

  std::cout << "}\n\n#endif" << std::endl;
}