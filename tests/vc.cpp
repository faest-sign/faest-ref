/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vc.h"
#include "fields.h"
#include "tree.h"
#include "compat.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(vector_commitments)

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

BOOST_AUTO_TEST_CASE(test_vector_open_128) {
  uint8_t rootKey[16] = {
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
      0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
  };

  faest_paramset_t params = faest_get_paramset(FAEST_128S); // Just using the FAEST-128s
  vec_com_t vecCom;
  uint32_t numVoleInstances = 16;
  uint32_t lambda           = params.faest_param.lambda;
  uint32_t lambdaBytes      = lambda / 8;
  vector_commitment(rootKey, &params, lambda, lambdaBytes, &vecCom, numVoleInstances);

  uint32_t leafIndex       = 7;
  constexpr uint32_t depth = 4;
  uint8_t b[depth];
  BitDec(leafIndex, depth, b);

  std::vector<uint8_t> pdec, com_j;
  pdec.resize(depth * lambdaBytes);
  com_j.resize(lambdaBytes * 2);
  vector_open(vecCom.k, vecCom.com, b, pdec.data(), com_j.data(), numVoleInstances, lambdaBytes);

  BOOST_TEST(memcmp(vecCom.com + (leafIndex * lambdaBytes * 2), com_j.data(), lambdaBytes * 2) ==
             0);
  BOOST_TEST(memcmp(pdec.data(), vecCom.k + (lambdaBytes * getNodeIndex(1, 1)), lambdaBytes) == 0);
  BOOST_TEST(memcmp(pdec.data() + lambdaBytes, vecCom.k + (lambdaBytes * getNodeIndex(2, 0)),
                    lambdaBytes) == 0);
  BOOST_TEST(memcmp(pdec.data() + (lambdaBytes * 2), vecCom.k + (lambdaBytes * getNodeIndex(3, 2)),
                    lambdaBytes) == 0);
  BOOST_TEST(memcmp(pdec.data() + (lambdaBytes * 3), vecCom.k + (lambdaBytes * getNodeIndex(4, 6)),
                    lambdaBytes) == 0);
}

BOOST_AUTO_TEST_CASE(test_vector_open_192) {
  uint8_t rootKey[24] = {
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b,
      0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
  };

  faest_paramset_t params = faest_get_paramset(FAEST_192S); // Just using the FAEST-192s
  vec_com_t vecCom;
  uint32_t numVoleInstances = 16;
  uint32_t lambda           = params.faest_param.lambda;
  uint32_t lambdaBytes      = lambda / 8;
  vector_commitment(rootKey, &params, lambda, lambdaBytes, &vecCom, numVoleInstances);

  uint32_t leafIndex       = 10;
  constexpr uint32_t depth = 4;
  uint8_t b[depth];
  BitDec(leafIndex, depth, b);

  std::vector<uint8_t> pdec, com_j;
  pdec.resize(depth * lambdaBytes);
  com_j.resize(lambdaBytes * 2);
  vector_open(vecCom.k, vecCom.com, b, pdec.data(), com_j.data(), numVoleInstances, lambdaBytes);

  BOOST_TEST(memcmp(vecCom.com + (leafIndex * lambdaBytes * 2), com_j.data(), lambdaBytes * 2) ==
             0);
  BOOST_TEST(memcmp(pdec.data(), vecCom.k + (lambdaBytes * getNodeIndex(1, 0)), lambdaBytes) == 0);
  BOOST_TEST(memcmp(pdec.data() + lambdaBytes, vecCom.k + (lambdaBytes * getNodeIndex(2, 3)),
                    lambdaBytes) == 0);
  BOOST_TEST(memcmp(pdec.data() + (lambdaBytes * 2), vecCom.k + (lambdaBytes * getNodeIndex(3, 4)),
                    lambdaBytes) == 0);
  BOOST_TEST(memcmp(pdec.data() + (lambdaBytes * 3), vecCom.k + (lambdaBytes * getNodeIndex(4, 11)),
                    lambdaBytes) == 0);
}

BOOST_AUTO_TEST_CASE(test_vector_open_256) {
  uint8_t rootKey[32] = {
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
      0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
      0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  };

  faest_paramset_t params = faest_get_paramset(FAEST_256S); // Just using the FAEST-256s
  vec_com_t vecCom;
  uint32_t numVoleInstances = 16;
  uint32_t lambda           = params.faest_param.lambda;
  uint32_t lambdaBytes      = lambda / 8;
  vector_commitment(rootKey, &params, lambda, lambdaBytes, &vecCom, numVoleInstances);

  uint32_t leafIndex       = 7;
  constexpr uint32_t depth = 4;
  uint8_t b[depth];
  BitDec(leafIndex, depth, b);

  std::vector<uint8_t> pdec, com_j;
  pdec.resize(depth * lambdaBytes);
  com_j.resize(lambdaBytes * 2);
  vector_open(vecCom.k, vecCom.com, b, pdec.data(), com_j.data(), numVoleInstances, lambdaBytes);

  BOOST_TEST(memcmp(vecCom.com + (leafIndex * lambdaBytes * 2), com_j.data(), lambdaBytes * 2) ==
             0);
  BOOST_TEST(memcmp(pdec.data(), vecCom.k + (lambdaBytes * getNodeIndex(1, 1)), lambdaBytes) == 0);
  BOOST_TEST(memcmp(pdec.data() + lambdaBytes, vecCom.k + (lambdaBytes * getNodeIndex(2, 0)),
                    lambdaBytes) == 0);
  BOOST_TEST(memcmp(pdec.data() + (lambdaBytes * 2), vecCom.k + (lambdaBytes * getNodeIndex(3, 2)),
                    lambdaBytes) == 0);
  BOOST_TEST(memcmp(pdec.data() + (lambdaBytes * 3), vecCom.k + (lambdaBytes * getNodeIndex(4, 6)),
                    lambdaBytes) == 0);
}

BOOST_AUTO_TEST_CASE(test_vector_reconstruct_and_verify) {
  uint8_t rootKey[16] = {
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
      0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
  };

  faest_paramset_t params = faest_get_paramset(FAEST_128S); // Just using the FAEST-128s
  vec_com_t vecCom;
  uint32_t numVoleInstances = 16;
  uint32_t lambda           = params.faest_param.lambda;
  uint32_t lambdaBytes      = lambda / 8;
  vector_commitment(rootKey, &params, lambda, lambdaBytes, &vecCom, numVoleInstances);

  uint32_t leafIndex       = 7;
  constexpr uint32_t depth = 4;
  uint8_t b[depth];
  BitDec(leafIndex, depth, b);

  std::vector<uint8_t> pdec, com_j;
  pdec.resize(depth * lambdaBytes);
  com_j.resize(lambdaBytes * 2);
  vector_open(vecCom.k, vecCom.com, b, pdec.data(), com_j.data(), numVoleInstances, lambdaBytes);

  BOOST_TEST(vector_verify(pdec.data(), com_j.data(), b, lambda, lambdaBytes, numVoleInstances,
                           nullptr, vecCom.h) == 1);
}

BOOST_AUTO_TEST_SUITE_END()
