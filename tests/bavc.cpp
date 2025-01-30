/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bavc.h"
#include "fields.h"
#include "compat.h"
#include "instances.hpp"
#include "bavc_tvs.hpp"
#include "utils.hpp"

#include <array>
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <random>
#include <vector>

namespace {
  constexpr std::array<uint8_t, 16> iv{};
  constexpr std::array<uint8_t, 32> root_key{
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
      0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
      0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  };
} // namespace

BOOST_AUTO_TEST_SUITE(bavc)

BOOST_DATA_TEST_CASE(test_node_indices, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const auto params = faest_get_paramset(param_id);
    const auto tau    = params.faest_param.tau;
    const auto tau_1  = params.faest_param.tau1;
    const auto k      = params.faest_param.k;
    const auto L      = params.faest_param.L;

    unsigned int recomputed_L = 0;
    for (unsigned int i = 0; i != tau; ++i) {
      recomputed_L += bavc_max_node_index(i, tau_1, k);
    }
    BOOST_TEST(L == recomputed_L);
  }
}

BOOST_AUTO_TEST_CASE(test_numrec_bitdec) {
  uint8_t expect_out_1[2] = {0x00, 0x01};
  uint8_t b_1[2];
  BitDec(2, 2, b_1);
  uint64_t idx_1 = NumRec(2, b_1);

  uint8_t expect_out_2[4] = {0x01, 0x01, 0x01, 0x00};
  uint8_t b_2[4];
  BitDec(7, 4, b_2);
  uint64_t idx_2 = NumRec(4, b_2);

  uint8_t expect_out_3[4] = {0x00, 0x01, 0x00, 0x01};
  uint8_t b_3[4];
  BitDec(10, 4, b_3);
  uint64_t idx_3 = NumRec(4, b_3);

  uint8_t expect_out_4[4] = {0x01, 0x00, 0x01, 0x01};
  uint8_t b_4[4];
  BitDec(13, 4, b_4);
  uint64_t idx_4 = NumRec(4, b_4);

  BOOST_TEST(memcmp(b_1, &expect_out_1, 2) == 0);
  BOOST_TEST(idx_1 == 2);
  BOOST_TEST(memcmp(b_2, &expect_out_2, 4) == 0);
  BOOST_TEST(idx_2 == 7);
  BOOST_TEST(memcmp(b_3, &expect_out_3, 4) == 0);
  BOOST_TEST(idx_3 == 10);
  BOOST_TEST(memcmp(b_4, &expect_out_4, 4) == 0);
  BOOST_TEST(idx_4 == 13);
}

namespace {
  template <size_t HSize, size_t IDeltaSize>
  void test_vc_tv(const faest_paramset_t& params, const std::array<uint16_t, IDeltaSize>& i_delta,
                  const std::array<uint8_t, HSize>& expected_h,
                  const std::array<uint8_t, 64>& expected_hashed_k,
                  const std::array<uint8_t, 64>& expected_hashed_sd,
                  const std::array<uint8_t, 64>& expected_hashed_decom_i,
                  const std::array<uint8_t, 64>& expected_hashed_rec_sd) {
    const auto lambda       = params.faest_param.lambda;
    const auto lambda_bytes = lambda / 8;
    const auto com_size     = (faest_is_em(&params) ? 2 : 3) * lambda_bytes;

    bavc_t vc;
    bavc_commit(root_key.data(), iv.data(), &params, &vc);

    const std::vector<uint8_t> h{vc.h, vc.h + HSize},
        expected_h_vec{expected_h.begin(), expected_h.end()};
    BOOST_TEST(h == expected_h_vec);
    // compare hashed k_alpha to reduce size of the TVs
    const auto hashed_k = hash_array(vc.k, (2 * params.faest_param.L - 1) * lambda_bytes);
    BOOST_TEST(hashed_k == expected_hashed_k);
    // compare hashed sd_ij to reduce size of the TVs
    const auto hashed_sd = hash_array(vc.sd, params.faest_param.L * lambda_bytes);
    BOOST_TEST(hashed_sd == expected_hashed_sd);

    std::vector<uint8_t> decom_i;
    decom_i.resize(com_size * params.faest_param.tau + params.faest_param.T_open * lambda_bytes);
    BOOST_TEST(bavc_open(&vc, i_delta.data(), decom_i.data(), &params));

    // compare hashed decom_i to reduce size of the TVs
    const auto hashed_decom_i = hash_array(decom_i);
    BOOST_TEST(hashed_decom_i == expected_hashed_decom_i);

    std::vector<uint8_t> rec_h, rec_s;
    rec_h.resize(2 * lambda_bytes);
    rec_s.resize((params.faest_param.L - params.faest_param.tau) * lambda_bytes);

    bavc_rec_t vc_rec;
    vc_rec.h = rec_h.data();
    vc_rec.s = rec_s.data();

    BOOST_TEST(bavc_reconstruct(decom_i.data(), i_delta.data(), iv.data(), &params, &vc_rec));
    BOOST_TEST(rec_h == expected_h_vec);

    // compare hashed sd_ij to reduce size of the TVs
    const auto hashed_rec_sd = hash_array(rec_s);
    BOOST_TEST(hashed_rec_sd == expected_hashed_rec_sd);
  }
} // namespace

BOOST_AUTO_TEST_CASE(tv_faest_128f) {
  test_vc_tv(faest_get_paramset(FAEST_128F), bavc_tvs::FAEST_128F::i_delta, bavc_tvs::FAEST_128F::h,
             bavc_tvs::FAEST_128F::hashed_k, bavc_tvs::FAEST_128F::hashed_sd,
             bavc_tvs::FAEST_128F::hashed_decom_i, bavc_tvs::FAEST_128F::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_128s) {
  test_vc_tv(faest_get_paramset(FAEST_128S), bavc_tvs::FAEST_128S::i_delta, bavc_tvs::FAEST_128S::h,
             bavc_tvs::FAEST_128S::hashed_k, bavc_tvs::FAEST_128S::hashed_sd,
             bavc_tvs::FAEST_128S::hashed_decom_i, bavc_tvs::FAEST_128S::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_192f) {
  test_vc_tv(faest_get_paramset(FAEST_192F), bavc_tvs::FAEST_192F::i_delta, bavc_tvs::FAEST_192F::h,
             bavc_tvs::FAEST_192F::hashed_k, bavc_tvs::FAEST_192F::hashed_sd,
             bavc_tvs::FAEST_192F::hashed_decom_i, bavc_tvs::FAEST_192F::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_192s) {
  test_vc_tv(faest_get_paramset(FAEST_192S), bavc_tvs::FAEST_192S::i_delta, bavc_tvs::FAEST_192S::h,
             bavc_tvs::FAEST_192S::hashed_k, bavc_tvs::FAEST_192S::hashed_sd,
             bavc_tvs::FAEST_192S::hashed_decom_i, bavc_tvs::FAEST_192S::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_256f) {
  test_vc_tv(faest_get_paramset(FAEST_256F), bavc_tvs::FAEST_256F::i_delta, bavc_tvs::FAEST_256F::h,
             bavc_tvs::FAEST_256F::hashed_k, bavc_tvs::FAEST_256F::hashed_sd,
             bavc_tvs::FAEST_256F::hashed_decom_i, bavc_tvs::FAEST_256F::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_256s) {
  test_vc_tv(faest_get_paramset(FAEST_256S), bavc_tvs::FAEST_256S::i_delta, bavc_tvs::FAEST_256S::h,
             bavc_tvs::FAEST_256S::hashed_k, bavc_tvs::FAEST_256S::hashed_sd,
             bavc_tvs::FAEST_256S::hashed_decom_i, bavc_tvs::FAEST_256S::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_em_128f) {
  test_vc_tv(faest_get_paramset(FAEST_EM_128F), bavc_tvs::FAEST_EM_128F::i_delta,
             bavc_tvs::FAEST_EM_128F::h, bavc_tvs::FAEST_EM_128F::hashed_k,
             bavc_tvs::FAEST_EM_128F::hashed_sd, bavc_tvs::FAEST_EM_128F::hashed_decom_i,
             bavc_tvs::FAEST_EM_128F::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_em_128s) {
  test_vc_tv(faest_get_paramset(FAEST_EM_128S), bavc_tvs::FAEST_EM_128S::i_delta,
             bavc_tvs::FAEST_EM_128S::h, bavc_tvs::FAEST_EM_128S::hashed_k,
             bavc_tvs::FAEST_EM_128S::hashed_sd, bavc_tvs::FAEST_EM_128S::hashed_decom_i,
             bavc_tvs::FAEST_EM_128S::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_em_192f) {
  test_vc_tv(faest_get_paramset(FAEST_EM_192F), bavc_tvs::FAEST_EM_192F::i_delta,
             bavc_tvs::FAEST_EM_192F::h, bavc_tvs::FAEST_EM_192F::hashed_k,
             bavc_tvs::FAEST_EM_192F::hashed_sd, bavc_tvs::FAEST_EM_192F::hashed_decom_i,
             bavc_tvs::FAEST_EM_192F::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_em_192s) {
  test_vc_tv(faest_get_paramset(FAEST_EM_192S), bavc_tvs::FAEST_EM_192S::i_delta,
             bavc_tvs::FAEST_EM_192S::h, bavc_tvs::FAEST_EM_192S::hashed_k,
             bavc_tvs::FAEST_EM_192S::hashed_sd, bavc_tvs::FAEST_EM_192S::hashed_decom_i,
             bavc_tvs::FAEST_EM_192S::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_em_256f) {
  test_vc_tv(faest_get_paramset(FAEST_EM_256F), bavc_tvs::FAEST_EM_256F::i_delta,
             bavc_tvs::FAEST_EM_256F::h, bavc_tvs::FAEST_EM_256F::hashed_k,
             bavc_tvs::FAEST_EM_256F::hashed_sd, bavc_tvs::FAEST_EM_256F::hashed_decom_i,
             bavc_tvs::FAEST_EM_256F::hashed_rec_sd);
}

BOOST_AUTO_TEST_CASE(tv_faest_em_256s) {
  test_vc_tv(faest_get_paramset(FAEST_EM_256S), bavc_tvs::FAEST_EM_256S::i_delta,
             bavc_tvs::FAEST_EM_256S::h, bavc_tvs::FAEST_EM_256S::hashed_k,
             bavc_tvs::FAEST_EM_256S::hashed_sd, bavc_tvs::FAEST_EM_256S::hashed_decom_i,
             bavc_tvs::FAEST_EM_256S::hashed_rec_sd);
}

BOOST_DATA_TEST_CASE(test_keys, all_parameters, param_id) {
  std::random_device rd;

  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const auto params       = faest_get_paramset(param_id);
    const auto lambda       = params.faest_param.lambda;
    const auto lambda_bytes = lambda / 8;

    bavc_t vc;
    bavc_commit(root_key.data(), iv.data(), &params, &vc);

    std::vector<uint8_t> decom_i;
    std::vector<uint16_t> i_delta;
    i_delta.resize(params.faest_param.tau);

    bool ret = false;
    while (!ret) {
      for (unsigned int i = 0; i != params.faest_param.tau; ++i) {
        std::uniform_int_distribution<> distribution{
            0, static_cast<int>(
                   bavc_max_node_index(i, params.faest_param.tau1, params.faest_param.k)) -
                   1};
        i_delta[i] = distribution(rd);
      }

      decom_i.clear();
      decom_i.resize(
          ((faest_is_em(&params) ? 2 : 3) * params.faest_param.tau + params.faest_param.T_open) *
          lambda_bytes);

      ret = bavc_open(&vc, i_delta.data(), decom_i.data(), &params);
    }
    BOOST_TEST(ret);

    std::vector<uint8_t> rec_h, rec_s;
    rec_h.resize(2 * lambda_bytes);
    rec_s.resize((params.faest_param.L - params.faest_param.tau) * lambda_bytes);

    bavc_rec_t vc_rec;
    vc_rec.h = rec_h.data();
    vc_rec.s = rec_s.data();

    BOOST_TEST(bavc_reconstruct(decom_i.data(), i_delta.data(), iv.data(), &params, &vc_rec));
    BOOST_TEST(memcmp(vc.h, vc_rec.h, 2 * lambda_bytes) == 0);

    bavc_clear(&vc);
  }
}

BOOST_AUTO_TEST_SUITE_END()
