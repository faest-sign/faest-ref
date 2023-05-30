/*
 *  SPDX-License-Identifier: MIT
 */

#include "vole.h"
#include "instances.hpp"
#include "randomness.h"

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
}

BOOST_AUTO_TEST_SUITE(vole)

BOOST_DATA_TEST_CASE(chal_dec, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t params  = faest_get_paramset(param_id);
    const unsigned int lambda      = params.faest_param.lambda;
    const unsigned int lambdaBytes = lambda / 8;

    std::vector<uint8_t> chal;
    chal.resize(lambdaBytes, 0xFF);

    for (unsigned int i = 0; i < params.faest_param.tau; i++) {
      const unsigned int depth =
          (i < params.faest_param.t0) ? params.faest_param.k0 : params.faest_param.k1;
      std::vector<uint8_t> chal_out;
      chal_out.resize(depth, 0);

      ChalDec(chal.data(), i, params.faest_param.k0, params.faest_param.t0, params.faest_param.k1,
              params.faest_param.t1, chal_out.data());
      for (unsigned int j = 0; j != depth; ++j) {
        BOOST_TEST(chal_out[j] == 1);
      }

      BOOST_TEST(NumRec(depth, chal_out.data()) == (1 << depth) - 1);
    }
  }
}

BOOST_DATA_TEST_CASE(vole_commit_verify, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t params  = faest_get_paramset(param_id);
    const unsigned int lambda      = params.faest_param.lambda;
    const unsigned int lambdaBytes = lambda / 8;
    const unsigned int ell_hat =
        params.faest_param.l + params.faest_param.lambda * 2 + params.faest_param.b;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;

    std::vector<uint8_t> hcom, hcomRec, u, b, chal;
    hcom.resize(lambdaBytes * 2);
    hcomRec.resize(lambdaBytes * 2);
    u.resize(ell_hat_bytes);
    b.resize(MAX(params.faest_param.k0, params.faest_param.k1), 0);
    chal.resize(lambdaBytes, 0);

    std::vector<vec_com_t> vec_com;
    vec_com.resize(params.faest_param.tau);

    std::vector<uint8_t*> c, v, q, pdec, com_j;
    c.resize(params.faest_param.tau - 1);
    v.resize(lambda);
    q.resize(lambda);
    pdec.resize(params.faest_param.tau);
    com_j.resize(params.faest_param.tau);

    v[0] = new uint8_t[lambda * ell_hat_bytes];
    q[0] = new uint8_t[lambda * ell_hat_bytes];
    for (unsigned int i = 1; i < lambda; ++i) {
      v[i] = v[0] + i * ell_hat_bytes;
      q[i] = q[0] + i * ell_hat_bytes;
    }

    voleCommit(rootKey.data(), ell_hat, &params, hcom.data(), vec_com.data(), c.data(), u.data(),
               v.data());
    for (unsigned int i = 0; i != params.faest_param.tau - 1; ++i) {
      free(c[i]);
    }

    for (uint32_t i = 0; i < params.faest_param.tau; i++) {
      const uint32_t depth =
          (i < params.faest_param.t0) ? params.faest_param.k0 : params.faest_param.k1;
      const uint32_t numVoleInstances = 1 << depth;

      pdec[i]  = new uint8_t[depth * lambdaBytes];
      com_j[i] = new uint8_t[lambdaBytes * 2];

      vector_open(vec_com[i].k, vec_com[i].com, b.data(), pdec[i], com_j[i], numVoleInstances,
                  lambdaBytes);
      vec_com_clear(&vec_com[i]);
    }

    voleReconstruct(chal.data(), pdec.data(), com_j.data(), hcomRec.data(), q.data(), ell_hat,
                    &params);
    BOOST_TEST(hcom == hcomRec);
    for (unsigned int i = 0; i < lambda; ++i) {
      // they should not be the same
      BOOST_TEST(std::memcmp(v[i], q[i], ell_hat_bytes) != 0);
    }

    for (uint32_t i = 0; i < params.faest_param.tau; i++) {
      delete[] com_j[i];
      delete[] pdec[i];
    }

    delete[] q[0];
    delete[] v[0];
  }
}

BOOST_DATA_TEST_CASE(convert_to_vole, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t params  = faest_get_paramset(param_id);
    const unsigned int lambda      = params.faest_param.lambda;
    const unsigned int lambdaBytes = lambda / 8;
    const unsigned int ell_hat =
        params.faest_param.l + params.faest_param.lambda * 2 + params.faest_param.b;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const unsigned int max_depth     = std::max(params.faest_param.k0, params.faest_param.k1);
    const unsigned int max_nodes     = 1 << max_depth;

    std::vector<uint8_t> sd, sdprime, u, v, q, chal_out, chal;
    sd.resize(max_nodes * lambdaBytes);
    rand_bytes(sd.data(), sd.size());
    sdprime.resize(max_nodes * lambdaBytes, 0);
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

      ConvertToVole(lambda, lambdaBytes, sd.data(), false, nodes, depth, ell_hat_bytes, u.data(),
                    v.data());

      ChalDec(chal.data(), i, params.faest_param.k0, params.faest_param.t0, params.faest_param.k1,
              params.faest_param.t1, chal_out.data());
      const auto idx = NumRec(depth, chal_out.data());
      for (unsigned int j = 1; j != nodes; ++j) {
        std::copy(&sd[(j ^ idx) * lambdaBytes], &sd[((j ^ idx) + 1) * lambdaBytes],
                  &sdprime[j * lambdaBytes]);
      }

      ConvertToVole(lambda, lambdaBytes, sdprime.data(), true, nodes, depth, ell_hat_bytes, nullptr,
                    q.data());

      for (unsigned int j = 0; j != depth; ++j) {
        for (unsigned int inner = 0; inner != ell_hat_bytes; ++inner) {
          q[j * ell_hat_bytes + inner] ^= chal_out[j] * u[inner];
        }
      }
      BOOST_TEST(q == v);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
