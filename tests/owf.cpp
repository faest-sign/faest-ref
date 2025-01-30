/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "owf.h"
#include "owf_tvs.hpp"
#include "instances.hpp"
#include "utils.hpp"

#include <array>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <random>
#include <vector>

BOOST_AUTO_TEST_SUITE(owf)

namespace {
  template <size_t KeySize, size_t InputSize, size_t OutputSize>
  void test_owf_tv(const faest_paramid_t param_id, const std::array<uint8_t, KeySize>& owf_key,
                   const std::array<uint8_t, InputSize>& owf_input,
                   const std::array<uint8_t, OutputSize>& expected_owf_output) {
    std::array<uint8_t, OutputSize> owf_output;
    switch (param_id) {
    case FAEST_128F:
      owf_128(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    case FAEST_192F:
      owf_192(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    case FAEST_256F:
      owf_256(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    case FAEST_EM_128F:
      owf_em_128(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    case FAEST_EM_192F:
      owf_em_192(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    case FAEST_EM_256F:
      owf_em_256(owf_key.data(), owf_input.data(), owf_output.data());
      break;
    default:
      BOOST_TEST(false);
    }

    BOOST_TEST(owf_output == expected_owf_output);
  }
} // namespace

BOOST_AUTO_TEST_CASE(tv_owf_128) {
  test_owf_tv(FAEST_128F, owf_tvs::FAEST_128F::key, owf_tvs::FAEST_128F::input,
              owf_tvs::FAEST_128F::output);
}

BOOST_AUTO_TEST_CASE(tv_owf_192) {
  test_owf_tv(FAEST_192F, owf_tvs::FAEST_192F::key, owf_tvs::FAEST_192F::input,
              owf_tvs::FAEST_192F::output);
}

BOOST_AUTO_TEST_CASE(tv_owf_256) {
  test_owf_tv(FAEST_256F, owf_tvs::FAEST_256F::key, owf_tvs::FAEST_256F::input,
              owf_tvs::FAEST_256F::output);
}

BOOST_AUTO_TEST_CASE(tv_owf_em_128) {
  test_owf_tv(FAEST_EM_128F, owf_tvs::FAEST_EM_128F::key, owf_tvs::FAEST_EM_128F::input,
              owf_tvs::FAEST_EM_128F::output);
}

BOOST_AUTO_TEST_CASE(tv_owf_em_192) {
  test_owf_tv(FAEST_EM_192F, owf_tvs::FAEST_EM_192F::key, owf_tvs::FAEST_EM_192F::input,
              owf_tvs::FAEST_EM_192F::output);
}

BOOST_AUTO_TEST_CASE(tv_owf_em_256) {
  test_owf_tv(FAEST_EM_256F, owf_tvs::FAEST_EM_256F::key, owf_tvs::FAEST_EM_256F::input,
              owf_tvs::FAEST_EM_256F::output);
}

BOOST_AUTO_TEST_SUITE_END()
