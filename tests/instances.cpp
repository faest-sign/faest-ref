/*
 *  SPDX-License-Identifier: MIT
 */

#include "instances.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <cmath>

BOOST_AUTO_TEST_SUITE(instances)

BOOST_DATA_TEST_CASE(test_keys, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const auto param        = faest_get_paramset(param_id);
    const auto& faest_param = param.faest_param;

    BOOST_TEST(faest_param.lambda <= MAX_LAMBDA);
    BOOST_TEST(
        (faest_param.lambda == 128 || faest_param.lambda == 192 || faest_param.lambda == 256));

    BOOST_TEST(faest_param.tau1 == (faest_param.lambda - faest_param.w_grind) % faest_param.tau);
    BOOST_TEST(faest_param.tau0 == faest_param.tau - faest_param.tau1);
    BOOST_TEST(faest_param.k ==
               std::ceil((faest_param.lambda - faest_param.w_grind) / (1.0 * faest_param.tau)));
    BOOST_TEST(faest_param.L == faest_param.tau1 * (1 << faest_param.k) +
                                    faest_param.tau0 * (1 << (faest_param.k - 1)));
    BOOST_TEST(faest_param.Lke % 8 == 0);
    BOOST_TEST(faest_param.Lenc % 8 == 0);
    BOOST_TEST(faest_param.l % 8 == 0);
    BOOST_TEST(faest_param.k <= MAX_DEPTH);
  }
}

BOOST_AUTO_TEST_SUITE_END()