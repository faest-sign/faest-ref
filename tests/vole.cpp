/*
 *  SPDX-License-Identifier: MIT
 */

#include "vole.h"
#include "instances.hpp"
#include "randomness.h"
#include "universal_hashing.h"
#include "utils.hpp"

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
    const auto params                = faest_get_paramset(param_id);
    const unsigned int lambda        = params->lambda;
    const unsigned int lambda_bytes  = lambda / 8;
    const unsigned int ell_hat       = params->l + params->lambda * 3 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const auto com_size              = (faest_is_em(params) ? 2 : 3) * lambda_bytes;

    bavc_t bavc_com;

    std::vector<uint8_t> chal, c, u, q_storage, v_storage;
    chal.resize(lambda_bytes);
    c.resize((params->tau - 1) * ell_hat_bytes);
    u.resize(ell_hat_bytes * params->tau);

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

    vole_commit(rootKey.data(), iv.data(), ell_hat, params, &bavc_com, c.data(), u.data(),
                v.data());

    std::vector<uint8_t> hcom{bavc_com.h, bavc_com.h + lambda_bytes * 2};

    bool tested = false;
    for (unsigned int tries = 0; !tested && tries != max_tries; ++tries) {
      rand_bytes(chal.data(), chal.size());
      for (unsigned int i = lambda - params->w_grind; i != lambda; ++i) {
        ptr_set_bit(chal.data(), i, 0);
      }

      std::vector<uint16_t> i_delta;
      i_delta.resize(params->tau);

      std::vector<uint8_t> decom_i;
      decom_i.resize(com_size * params->tau + params->T_open * lambda_bytes);

      BOOST_TEST(decode_all_chall_3(i_delta.data(), chal.data(), params));
      if (!bavc_open(decom_i.data(), &bavc_com, i_delta.data(), params)) {
        continue;
      }
      tested = true;

      std::vector<uint8_t> hcom_rec;
      hcom_rec.resize(lambda_bytes * 2);
      BOOST_TEST(vole_reconstruct(hcom_rec.data(), q.data(), iv.data(), chal.data(), decom_i.data(),
                                  c.data(), ell_hat, params));
      BOOST_TEST(hcom == hcom_rec);
    }
    BOOST_TEST(tested);
    bavc_clear(&bavc_com);
  }
}

BOOST_DATA_TEST_CASE(test_convert_to_vole, all_parameters, param_id) {
  std::mt19937_64 rd;
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t params    = *faest_get_paramset(param_id);
    const unsigned int lambda        = params.lambda;
    const unsigned int lambda_bytes  = lambda / 8;
    const unsigned int ell_hat_bytes = 16;
    const unsigned int max_depth     = params.k;
    const unsigned int max_nodes     = 1 << max_depth;
    const unsigned int tau           = params.tau;

    std::vector<uint8_t> sd, u, v, q;
    sd.resize(max_nodes * lambda_bytes);
    rand_bytes(sd.data(), sd.size());
    u.resize(ell_hat_bytes);
    v.resize(ell_hat_bytes * max_depth);
    q.resize(ell_hat_bytes * max_depth);

    for (unsigned int i = 0; i != tau; ++i) {
      std::uniform_int_distribution<> distribution{
          0, static_cast<int>(bavc_max_node_index(i, params.tau1, max_depth)) - 1};
      const unsigned int idx = distribution(rd);

      unsigned int depth = bavc_max_node_depth(i, params.tau1, max_depth);
      unsigned int nodes = 1 << depth;

      convert_to_vole(iv.data(), sd.data(), false, i, ell_hat_bytes, u.data(), v.data(), &params);

      std::vector<uint8_t> sdprime;
      sdprime.resize(nodes * lambda_bytes, 0);
      for (unsigned int j = 1; j != nodes; ++j) {
        std::copy(sd.begin() + (j ^ idx) * lambda_bytes,
                  sd.begin() + ((j ^ idx) + 1) * lambda_bytes, &sdprime[j * lambda_bytes]);
      }

      convert_to_vole(iv.data(), sdprime.data(), true, i, ell_hat_bytes, nullptr, q.data(),
                      &params);

      for (unsigned int j = 0; j != depth; ++j) {
        if (idx & (1 << j)) {
          for (unsigned int inner = 0; inner != ell_hat_bytes; ++inner) {
            q[j * ell_hat_bytes + inner] ^= u[inner];
          }
        }
      }
      BOOST_TEST(q == v);
    }
  }
}

namespace {
  template <std::size_t HSize, std::size_t CSize>
  void test_tv(const faest_paramset_t* params, const std::array<uint8_t, CSize>& challenge,
               const std::array<uint8_t, HSize>& expected_h,
               const std::array<uint8_t, 64>& expected_hashed_c,
               const std::array<uint8_t, 64>& expected_hashed_u,
               const std::array<uint8_t, 64>& expected_hashed_v,
               const std::array<uint8_t, 64>& expected_hashed_q) {
    const unsigned int lambda        = params->lambda;
    const unsigned int lambda_bytes  = lambda / 8;
    const unsigned int ell_hat       = params->l + params->lambda * 3 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const auto com_size              = (faest_is_em(params) ? 2 : 3) * lambda_bytes;

    bavc_t bavc_com;

    std::vector<uint8_t> chal, c, decom_i, u, q_storage, v_storage;
    chal.resize(lambda_bytes);
    c.resize((params->tau - 1) * ell_hat_bytes);
    decom_i.resize(com_size * params->tau + params->T_open * lambda_bytes);
    u.resize(ell_hat_bytes * params->tau);

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

    vole_commit(rootKey.data(), iv.data(), ell_hat, params, &bavc_com, c.data(), u.data(),
                v.data());
    std::vector<uint8_t> hcom{bavc_com.h, bavc_com.h + lambda_bytes * 2},
        expected_h_vec{expected_h.begin(), expected_h.end()};
    BOOST_TEST(hcom == expected_h_vec);
    BOOST_TEST(expected_hashed_c == hash_array(c));
    BOOST_TEST(expected_hashed_u == hash_array(u));
    BOOST_TEST(expected_hashed_v == hash_array(v_storage));

    std::vector<uint16_t> i_delta;
    i_delta.resize(params->tau);
    BOOST_TEST(decode_all_chall_3(i_delta.data(), challenge.data(), params));
    BOOST_TEST(bavc_open(decom_i.data(), &bavc_com, i_delta.data(), params));

    std::vector<uint8_t> hcom_rec;
    hcom_rec.resize(lambda_bytes * 2);
    BOOST_TEST(vole_reconstruct(hcom_rec.data(), q.data(), iv.data(), chal.data(), decom_i.data(),
                                c.data(), ell_hat, params));
    BOOST_TEST(hcom_rec == expected_h_vec);
    BOOST_TEST(expected_hashed_q == hash_array(q_storage));
    bavc_clear(&bavc_com);

    for (unsigned int i = 0, running_idx = 0; i < params->tau; ++i) {
      const uint32_t depth = bavc_max_node_depth(i, params->tau1, params->k);
      const auto delta     = i_delta[i];

      for (unsigned int j = 0; j != depth; ++j, ++running_idx) {
        for (unsigned int inner = 0; inner != ell_hat_bytes; ++inner) {
          if ((delta >> j) & 1) {
            // need to correct the vole correlation
            if (i > 0) {
              BOOST_TEST((q[(running_idx)][inner] ^ c[(i - 1) * ell_hat_bytes + inner] ^
                          u[inner]) == v[(running_idx)][inner]);
            } else {
              BOOST_TEST((q[(running_idx)][inner] ^ u[inner]) == v[(running_idx)][inner]);
            }
          } else {
            BOOST_TEST(q[(running_idx)][inner] == v[(running_idx)][inner]);
          }
        }
      }
    }
  }
} // namespace

BOOST_AUTO_TEST_SUITE_END()
