/*
 *  SPDX-License-Identifier: MIT
 */

#include "vole.h"
#include "instances.hpp"

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

BOOST_DATA_TEST_CASE(vole_commit_verify, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t params = faest_get_paramset(param_id);

    const uint32_t lambda      = params.faest_param.lambda;
    const uint32_t lambdaBytes = lambda / 8;
    const uint32_t ell_hat =
        params.faest_param.l + params.faest_param.lambda * 2 + params.faest_param.b;
    const uint32_t ell_hat_bytes = (ell_hat + 7) / 8;

    std::vector<uint8_t> hcom, hcomRec, u, b, chal;
    hcom.resize(lambdaBytes * 2);
    hcomRec.resize(lambdaBytes * 2);
    u.resize(ell_hat_bytes);
    b.resize(MAX(params.faest_param.k0, params.faest_param.k1), 0);
    chal.resize(lambdaBytes);

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

BOOST_AUTO_TEST_SUITE_END()