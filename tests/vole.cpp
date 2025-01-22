/*
 *  SPDX-License-Identifier: MIT
 */

#include "vole.h"
#include "instances.hpp"
#include "randomness.h"
#include "universal_hashing.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <cmath>
#include <array>
#include <vector>

namespace {
  constexpr std::array<uint8_t, 32> rootKey{
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
      0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
      0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  };
  constexpr std::array<uint8_t, 16> iv{};

  constexpr unsigned int max_tries = 1000;
} // namespace

BOOST_AUTO_TEST_SUITE(vole)

BOOST_DATA_TEST_CASE(vole_commit_verify, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t params   = faest_get_paramset(param_id);
    const unsigned int lambda       = params.faest_param.lambda;
    const unsigned int lambda_bytes = lambda / 8;
    const unsigned int ell_hat =
        params.faest_param.l + params.faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const auto com_size              = (faest_is_em(&params) ? 2 : 3) * lambda_bytes;

    vec_com_t bavc_com;

    std::vector<uint8_t> chal, c, decom_i, u, q_storage, v_storage;
    chal.resize(lambda_bytes);
    c.resize((params.faest_param.tau - 1) * ell_hat_bytes);
    decom_i.resize(com_size * params.faest_param.tau + params.faest_param.T_open * lambda_bytes);
    u.resize(ell_hat_bytes * params.faest_param.tau);

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

    vole_commit(rootKey.data(), iv.data(), ell_hat, &params, &bavc_com, c.data(), u.data(),
                v.data());

    std::vector<uint8_t> hcom{bavc_com.h, bavc_com.h + lambda_bytes * 2};

    bool tested = false;
    for (unsigned int tries = 0; tries != max_tries; ++tries) {
      rand_bytes(chal.data(), chal.size());
      for (unsigned int i = lambda - params.faest_param.w_grind; i != lambda; ++i) {
        ptr_set_bit(chal.data(), 0, i);
      }

      uint16_t i_delta[MAX_TAU];
      BOOST_TEST(decode_all_chall_3(i_delta, chal.data(), &params));
      if (!bavc_open(&bavc_com, i_delta, decom_i.data(), &params)) {
        continue;
      }

      std::vector<uint8_t> hcom_rec;
      hcom_rec.resize(lambda_bytes * 2);
      BOOST_TEST(vole_reconstruct(hcom_rec.data(), q.data(), iv.data(), chal.data(), decom_i.data(),
                                  c.data(), ell_hat, &params));
      BOOST_TEST(hcom == hcom_rec);
      tested = true;
    }
    BOOST_TEST(tested);
    vec_com_clear(&bavc_com);
  }
}

#if 0
BOOST_DATA_TEST_CASE(convert_to_vole, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t params  = faest_get_paramset(param_id);
    const unsigned int lambda      = params.faest_param.lambda;
    const unsigned int lambdaBytes = lambda / 8;
    const unsigned int ell_hat =
        params.faest_param.l + params.faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const unsigned int max_depth     = std::max(params.faest_param.k0, params.faest_param.k1);
    const unsigned int max_nodes     = 1 << max_depth;

    std::vector<uint8_t> sd, u, v, q, chal_out, chal;
    sd.resize(max_nodes * lambdaBytes);
    rand_bytes(sd.data(), sd.size());
    chal_out.resize(max_depth);
    u.resize(ell_hat_bytes);
    v.resize(ell_hat_bytes * max_depth);
    q.resize(ell_hat_bytes * max_depth);
    chal.resize(lambdaBytes);
    rand_bytes(chal.data(), chal.size());

    for (unsigned int i = 0; i != params.faest_param.tau; ++i) {
      unsigned int depth =
          i < params.faest_param.t0 ? params.faest_param.k0 : params.faest_param.k1;
      unsigned int nodes = 1 << depth;

      ConvertToVole(iv.data(), sd.data(), false, lambda, depth, ell_hat_bytes, u.data(), v.data());

      ChalDec(chal.data(), i, params.faest_param.k0, params.faest_param.t0, params.faest_param.k1,
              params.faest_param.t1, chal_out.data());
      const auto idx = NumRec(depth, chal_out.data());
      std::vector<uint8_t> sdprime;
      sdprime.resize(max_nodes * lambdaBytes, 0);
      for (unsigned int j = 1; j != nodes; ++j) {
        std::copy(&sd[(j ^ idx) * lambdaBytes], &sd[((j ^ idx) + 1) * lambdaBytes],
                  &sdprime[j * lambdaBytes]);
      }

      ConvertToVole(iv.data(), sdprime.data(), true, lambda, depth, ell_hat_bytes, nullptr,
                    q.data());

      for (unsigned int j = 0; j != depth; ++j) {
        BOOST_TEST((chal_out[j] == 0 || chal_out[j] == 1));
        if (chal_out[j]) {
          for (unsigned int inner = 0; inner != ell_hat_bytes; ++inner) {
            q[j * ell_hat_bytes + inner] ^= u[inner];
          }
        }
      }
      BOOST_TEST(q == v);
    }
  }
}
#endif

BOOST_AUTO_TEST_SUITE_END()
