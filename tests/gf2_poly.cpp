/*
 *  SPDX-License-Identifier: MIT
 */

#include "instances.h"
#include "fields.hpp"
#include "crt_vole_commit_tvs.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <array>
#include <vector>

BOOST_AUTO_TEST_SUITE(bf2_matrix)

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_128f) {
  namespace tv = crt_tvs_128f;

  const auto& params = *faest_get_paramset(FAEST_128F);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_128s) {
  namespace tv = crt_tvs_128s;

  const auto& params = *faest_get_paramset(FAEST_128S);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_192f) {
  namespace tv = crt_tvs_192f;

  const auto& params = *faest_get_paramset(FAEST_192F);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_192s) {
  namespace tv = crt_tvs_192s;

  const auto& params = *faest_get_paramset(FAEST_192S);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_256f) {
  namespace tv = crt_tvs_256f;

  const auto& params = *faest_get_paramset(FAEST_256F);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_256s) {
  namespace tv = crt_tvs_256s;

  const auto& params = *faest_get_paramset(FAEST_256S);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_em_128f) {
  namespace tv = crt_tvs_em_128f;

  const auto& params = *faest_get_paramset(FAEST_EM_128F);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_em_128s) {
  namespace tv = crt_tvs_em_128s;

  const auto& params = *faest_get_paramset(FAEST_EM_128S);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_em_192f) {
  namespace tv = crt_tvs_em_192f;

  const auto& params = *faest_get_paramset(FAEST_EM_192F);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_em_192s) {
  namespace tv = crt_tvs_em_192s;

  const auto& params = *faest_get_paramset(FAEST_EM_192S);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_em_256f) {
  namespace tv = crt_tvs_em_256f;

  const auto& params = *faest_get_paramset(FAEST_EM_256F);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_CASE(test_bf2_matrix_mul_tbl_em_256s) {
  namespace tv = crt_tvs_em_256s;

  const auto& params = *faest_get_paramset(FAEST_EM_256S);

  auto delta = decltype(tv::expected_crt_delta_bytes){0};
  bf2_matrix_mul_tbl(delta.data(), tv::delta_bytes.data(), params.W_CRT,
                     params.lambda - params.w_grind, params.lambda, params.w_crt_words);
  BOOST_TEST(delta == tv::expected_crt_delta_bytes);
}

BOOST_AUTO_TEST_SUITE_END()