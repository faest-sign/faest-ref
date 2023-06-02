/*
 *  SPDX-License-Identifier: MIT
 */

#include "faest_aes.h"
#include "instances.hpp"
#include "randomness.h"
#include "tvs_aes.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <cmath>
#include <array>
#include <vector>

BOOST_AUTO_TEST_SUITE(test_aes_prove)

BOOST_DATA_TEST_CASE(aes_prove_verify, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t params  = faest_get_paramset(param_id);
    const bool is_em               = params.faest_param.Lke == 0;
    const unsigned int lambda      = params.faest_param.lambda;
    const unsigned int lambdaBytes = lambda / 8;
    const unsigned int ell_hat =
        params.faest_param.l + params.faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const size_t block_size          = is_em ? lambdaBytes : 16;

    // extended witness
    std::vector<uint8_t> w((params.faest_param.l + 7) / 8, 0);
    std::vector<uint8_t> in(params.faest_param.beta * block_size, 0);
    std::vector<uint8_t> out(params.faest_param.beta * block_size, 0);

    if (lambda == 128 && !is_em) {
      BOOST_TEST(in.size() == aes_ctr_128_tv::in.size());
      BOOST_TEST(out.size() == aes_ctr_128_tv::out.size());
      BOOST_TEST(w.size() == aes_ctr_128_tv::expected_extended_witness.size());
      memcpy(in.data(), aes_ctr_128_tv::in.data(), in.size());
      memcpy(out.data(), aes_ctr_128_tv::out.data(), out.size());
      memcpy(w.data(), aes_ctr_128_tv::expected_extended_witness.data(), w.size());
    } else if (lambda == 128 && is_em) {
      BOOST_TEST(in.size() == rijndael_em_128_tv::in.size());
      BOOST_TEST(out.size() == rijndael_em_128_tv::out.size());
      BOOST_TEST(w.size() == rijndael_em_128_tv::expected_extended_witness.size());
      memcpy(in.data(), rijndael_em_128_tv::in.data(), in.size());
      memcpy(out.data(), rijndael_em_128_tv::out.data(), out.size());
      memcpy(w.data(), rijndael_em_128_tv::expected_extended_witness.data(), w.size());
    } else if (lambda == 192 && !is_em) {
      BOOST_TEST(in.size() == aes_ctr_192_tv::in.size());
      BOOST_TEST(out.size() == aes_ctr_192_tv::out.size());
      BOOST_TEST(w.size() == aes_ctr_192_tv::expected_extended_witness.size());
      memcpy(in.data(), aes_ctr_192_tv::in.data(), in.size());
      memcpy(out.data(), aes_ctr_192_tv::out.data(), out.size());
      memcpy(w.data(), aes_ctr_192_tv::expected_extended_witness.data(), w.size());
    } else if (lambda == 192 && is_em) {
      BOOST_TEST(in.size() == rijndael_em_192_tv::in.size());
      BOOST_TEST(out.size() == rijndael_em_192_tv::out.size());
      BOOST_TEST(w.size() == rijndael_em_192_tv::expected_extended_witness.size());
      memcpy(in.data(), rijndael_em_192_tv::in.data(), in.size());
      memcpy(out.data(), rijndael_em_192_tv::out.data(), out.size());
      memcpy(w.data(), rijndael_em_192_tv::expected_extended_witness.data(), w.size());
    } else if (lambda == 256 && !is_em) {
      BOOST_TEST(in.size() == aes_ctr_256_tv::in.size());
      BOOST_TEST(out.size() == aes_ctr_256_tv::out.size());
      BOOST_TEST(w.size() == aes_ctr_256_tv::expected_extended_witness.size());
      memcpy(in.data(), aes_ctr_256_tv::in.data(), in.size());
      memcpy(out.data(), aes_ctr_256_tv::out.data(), out.size());
      memcpy(w.data(), aes_ctr_256_tv::expected_extended_witness.data(), w.size());
    } else if (lambda == 256 && is_em) {
      BOOST_TEST(in.size() == rijndael_em_256_tv::in.size());
      BOOST_TEST(out.size() == rijndael_em_256_tv::out.size());
      BOOST_TEST(w.size() == rijndael_em_256_tv::expected_extended_witness.size());
      memcpy(in.data(), rijndael_em_256_tv::in.data(), in.size());
      memcpy(out.data(), rijndael_em_256_tv::out.data(), out.size());
      memcpy(w.data(), rijndael_em_256_tv::expected_extended_witness.data(), w.size());
    }

    // prepare vole correlation
    std::vector<uint8_t> delta(lambda / 8, 0);
    delta[0] = 42;
    std::vector<uint8_t> u(ell_hat_bytes, 0x13);
    std::vector<uint8_t> vs(ell_hat_bytes * lambda, 0x37);
    std::vector<uint8_t> qs = vs;
    std::vector<uint8_t*> V(lambda, NULL);
    std::vector<uint8_t*> Q(lambda, NULL);

    for (size_t i = 0; i < lambda; ++i) {
      V[i] = vs.data() + i * ell_hat_bytes;
      Q[i] = qs.data() + i * ell_hat_bytes;
      if ((delta[i / 8] >> (i % 8)) & 1) {
        for (size_t j = 0; j < ell_hat_bytes; ++j) {
          Q[i][j] ^= u[j];
        }
      }
    }

    // masked witness d = u ^ w
    std::vector<uint8_t> d(ell_hat_bytes, 0x13);
    for (size_t i = 0; i < w.size(); ++i) {
      d[i] = u[i] ^ w[i];
    }

    std::vector<uint8_t> chall_2((3 * lambda + 64) / 8, 47);

    std::vector<uint8_t> a_tilde(lambda / 8, 0);
    std::vector<uint8_t> b_tilde(lambda / 8, 0);

    aes_prove(w.data(), u.data(), V.data(), in.data(), out.data(), chall_2.data(), a_tilde.data(),
              b_tilde.data(), &params);

    uint8_t* recomputed_b_tilde = aes_verify(d.data(), Q.data(), chall_2.data(), delta.data(),
                                             a_tilde.data(), in.data(), out.data(), &params);

    // check that the proof verifies
    BOOST_TEST(memcmp(recomputed_b_tilde, b_tilde.data(), lambdaBytes) == 0);
    free(recomputed_b_tilde);
  }
}

BOOST_AUTO_TEST_SUITE_END()
