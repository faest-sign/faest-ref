/*
 *  SPDX-License-Identifier: MIT
 */

#include "instances.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <cmath>

BOOST_AUTO_TEST_SUITE(instances)

BOOST_DATA_TEST_CASE(instance_parameters, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const auto faest_param = *faest_get_paramset(param_id);

    BOOST_TEST(faest_param.lambda <= MAX_LAMBDA);
    BOOST_TEST(
        (faest_param.lambda == 128 || faest_param.lambda == 192 || faest_param.lambda == 256));

    BOOST_TEST(faest_param.tau1 == (faest_param.lambda - faest_param.w_grind) % faest_param.tau);
    BOOST_TEST(faest_param.tau0 == faest_param.tau - faest_param.tau1);
    BOOST_TEST(faest_param.k == std::floor((faest_param.lambda - faest_param.w_grind) /
                                           (static_cast<double>(faest_param.tau))) +
                                    1);
    BOOST_TEST(faest_param.lambda - faest_param.w_grind ==
               (faest_param.k - 1) * faest_param.tau0 + faest_param.k * faest_param.tau1);
    BOOST_TEST(faest_param.L == faest_param.tau1 * (1 << faest_param.k) +
                                    faest_param.tau0 * (1 << (faest_param.k - 1)));
    BOOST_TEST(faest_param.Lke % 8 == 0);
    BOOST_TEST(faest_param.Lenc % 8 == 0);
    BOOST_TEST(faest_param.l % 8 == 0);
    BOOST_TEST(faest_param.k <= sizeof(uint16_t) * 8);
    BOOST_TEST(faest_param.k <= MAX_DEPTH);
    BOOST_TEST(faest_param.tau <= MAX_TAU);

    const auto ell_bytes    = (faest_param.l + 7) / 8;
    const auto lambda_bytes = faest_param.lambda / 8;
    const auto n_leafcom    = faest_is_em(&faest_param) ? 2 : 3;

    const auto sig_size = faest_param.tau * (ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B) +
                          faest_param.T_open * lambda_bytes +
                          n_leafcom * lambda_bytes * faest_param.tau + lambda_bytes + IV_SIZE +
                          sizeof(uint32_t);
    BOOST_TEST(sig_size == faest_param.sig_size);

    const auto beta = faest_is_em(&faest_param) ? 1 : (faest_param.lambda + 127) / 128;
    BOOST_TEST(faest_param.l == faest_param.Lke + beta * faest_param.Lenc);
  }
}

BOOST_AUTO_TEST_SUITE_END()