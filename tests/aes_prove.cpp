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

    // extended witness
    std::vector<uint8_t> w;
    std::vector<uint8_t> in;
    std::vector<uint8_t> out;

    if (lambda == 128 && !is_em) {
      std::copy(aes_ctr_128_tv::in.begin(), aes_ctr_128_tv::in.end(),
                std::back_insert_iterator(in));
      std::copy(aes_ctr_128_tv::out.begin(), aes_ctr_128_tv::out.end(),
                std::back_insert_iterator(out));
      std::copy(aes_ctr_128_tv::expected_extended_witness.begin(),
                aes_ctr_128_tv::expected_extended_witness.end(), std::back_insert_iterator(w));
    } else if (lambda == 128 && is_em) {
      std::copy(rijndael_em_128_tv::in.begin(), rijndael_em_128_tv::in.end(),
                std::back_insert_iterator(in));
      std::copy(rijndael_em_128_tv::out.begin(), rijndael_em_128_tv::out.end(),
                std::back_insert_iterator(out));
      std::copy(rijndael_em_128_tv::expected_extended_witness.begin(),
                rijndael_em_128_tv::expected_extended_witness.end(), std::back_insert_iterator(w));
    } else if (lambda == 192 && !is_em) {
      std::copy(aes_ctr_192_tv::in.begin(), aes_ctr_192_tv::in.end(),
                std::back_insert_iterator(in));
      std::copy(aes_ctr_192_tv::out.begin(), aes_ctr_192_tv::out.end(),
                std::back_insert_iterator(out));
      std::copy(aes_ctr_192_tv::expected_extended_witness.begin(),
                aes_ctr_192_tv::expected_extended_witness.end(), std::back_insert_iterator(w));
    } else if (lambda == 192 && is_em) {
      std::copy(rijndael_em_192_tv::in.begin(), rijndael_em_192_tv::in.end(),
                std::back_insert_iterator(in));
      std::copy(rijndael_em_192_tv::out.begin(), rijndael_em_192_tv::out.end(),
                std::back_insert_iterator(out));
      std::copy(rijndael_em_192_tv::expected_extended_witness.begin(),
                rijndael_em_192_tv::expected_extended_witness.end(), std::back_insert_iterator(w));
    } else if (lambda == 256 && !is_em) {
      std::copy(aes_ctr_256_tv::in.begin(), aes_ctr_256_tv::in.end(),
                std::back_insert_iterator(in));
      std::copy(aes_ctr_256_tv::out.begin(), aes_ctr_256_tv::out.end(),
                std::back_insert_iterator(out));
      std::copy(aes_ctr_256_tv::expected_extended_witness.begin(),
                aes_ctr_256_tv::expected_extended_witness.end(), std::back_insert_iterator(w));
    } else if (lambda == 256 && is_em) {
      std::copy(rijndael_em_256_tv::in.begin(), rijndael_em_256_tv::in.end(),
                std::back_insert_iterator(in));
      std::copy(rijndael_em_256_tv::out.begin(), rijndael_em_256_tv::out.end(),
                std::back_insert_iterator(out));
      std::copy(rijndael_em_256_tv::expected_extended_witness.begin(),
                rijndael_em_256_tv::expected_extended_witness.end(), std::back_insert_iterator(w));
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
