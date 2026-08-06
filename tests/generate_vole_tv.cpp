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
    const unsigned int ell            = params.ell;
    const unsigned int ell_bytes      = (ell + 7) / 8;
    const unsigned int k              = params.k;
    const unsigned int tau1           = params.tau1;
    const unsigned int d0            = bavc_max_node_depth(0, tau1, k);
    const unsigned int n_mask        = params.n_mask;
    const unsigned int ell_hat       = ell + n_mask * d0;
    const unsigned int ell_hat_bytes = (ell_hat +7) / 8;
    const auto com_size              = (faest_is_em(&params) ? 2 : 3) * lambda_bytes;
    const auto n_mask_bytes           = (n_mask + 7) / 8;
    const auto n_mult                 = params.n_mult;

    bavc_t bavc_com;

    std::vector<uint8_t> chal, c, c_mult, decom_i, u, u_bar, v_bar, q_bar, delta_prime, q_storage, v_storage;
    chal.resize(lambda_bytes);
    c.resize((params.tau - 1) * ell_bytes);
    c_mult.resize(n_mult * n_mask_bytes);
    decom_i.resize(com_size * params.tau + params.T_open * lambda_bytes);
    u.resize(ell_bytes);
    u_bar.resize(n_mask * lambda_bytes);
    v_bar.resize(n_mask * lambda_bytes);
    q_bar.resize(n_mask * lambda_bytes);
    delta_prime.resize(lambda_bytes);

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

    // TODO: Remove
    // NOTE: Just for debugging with a chall TV instead of the random chall from above
    // static const char* chall_hex = "c3810c22de83098321348f823551318d7b22fc5e63c48700";
    // for (unsigned int i = 0; i < lambda_bytes; ++i) {
    //   unsigned int byte;
    //   std::sscanf(chall_hex + 2 * i, "%2x", &byte);
    //   chal[i] = static_cast<uint8_t>(byte);
    // }

    vole_commit(root_key.data(), iv.data(), ell_hat, &params, &bavc_com, c.data(), c_mult.data(),
                u.data(), v.data(), u_bar.data(), v_bar.data());

    print_named_array("h", "uint8_t", bavc_com.h, 2 * lambda_bytes);
    print_named_array("hashed_c", "uint8_t", hash_array(c));
    print_named_array("hashed_u", "uint8_t", hash_array(u));

    // NOTE: --- Claude generated code 
    // This arranges the c_mult into the required ordering similar to pyfaest implementation and its TVs 
    size_t nb        = (n_mult + 7) / 8;
    size_t packed_len = n_mask * nb;
    std::vector<uint8_t> c_mult_packed(packed_len, 0);
    for (size_t i = 0; i < n_mask; ++i) {
        for (size_t j = 0; j < n_mult; ++j) {
          uint16_t word = c_mult[j*2] | (uint16_t)(c_mult[j*2 + 1] << 8);
          uint8_t  bit  = (word >> i) & 1u;
          size_t   pos  = i * nb * 8 + j;
          c_mult_packed[pos >> 3] |= (uint8_t)(bit << (pos & 7));
        }
    }
    // ----
    print_named_array("hashed_c_mult", "uint8_t", hash_array(c_mult_packed));

    // NOTE: --- Claude generated code 
    // This arranges the v_storage into the required ordering similar to pyfaest implementation and its TVs 
    size_t ncols        = ell;
    std::vector<uint8_t> mV_bytes(ncols * lambda_bytes, 0);
    for (size_t ncols_idx = 0; ncols_idx < ncols; ++ncols_idx) {
      for (size_t r = 0; r < lambda; ++r) {
          const uint8_t* srow = v_storage.data() + r * ell_hat_bytes;
          uint8_t bit = (srow[ncols_idx >> 3] >> (ncols_idx & 7)) & 1u;
          size_t  pos = ncols_idx * lambda + r;
          mV_bytes[pos >> 3] |= (uint8_t)(bit << (pos & 7));
      }
    }
    // ---

    print_named_array("hashed_v", "uint8_t", hash_array(mV_bytes));
    print_named_array("hashed_barU", "uint8_t", hash_array(u_bar));
    print_named_array("hashed_barV", "uint8_t", hash_array(v_bar));

    while (true) {
      // TODO: Uncomment to get the random chall
      std::generate(chal.begin(), chal.end(), [&mt, &dist] { return dist(mt); });
      for (unsigned int i = lambda - params.w_grind; i != lambda; ++i) {
        ptr_set_bit(chal.data(), i, 0);
      }

      uint16_t i_delta[MAX_TAU];
      decode_all_chall_3(i_delta, chal.data(), &params);
      if (!bavc_open(decom_i.data(), &bavc_com, i_delta, &params)) {
        continue;
      }

      std::vector<uint8_t> hcom_rec;
      hcom_rec.resize(lambda_bytes * 2);
      vole_reconstruct(hcom_rec.data(), q.data(), iv.data(), chal.data(), decom_i.data(),
                                  c.data(), c_mult.data(), q_bar.data(), delta_prime.data(), ell_hat, &params);

      // NOTE: --- Claude generated code 
      // This arranges the q_storage into the required ordering similar to pyfaest implementation and its TVs
      std::vector<uint8_t> mQ_bytes(ncols * lambda_bytes, 0);
      for (size_t ncols_idx = 0; ncols_idx < ncols; ++ncols_idx) {
        for (size_t r = 0; r < lambda; ++r) {
            const uint8_t* srow = q_storage.data() + r * ell_hat_bytes;
            uint8_t bit = (srow[ncols_idx >> 3] >> (ncols_idx & 7)) & 1u;
            size_t  pos = ncols_idx * lambda + r;
            mQ_bytes[pos >> 3] |= (uint8_t)(bit << (pos & 7));
        }
      }
      // ----
      print_named_array("hashed_q", "uint8_t", hash_array(mQ_bytes));
      print_named_array("hashed_barQ", "uint8_t", hash_array(q_bar));

      break;
    }
    
    print_named_array("chall", "uint8_t", chal);

    bavc_clear(&bavc_com);
    std::cout << "}\n";
  }



  std::cout << "}\n\n#endif" << std::endl;
}