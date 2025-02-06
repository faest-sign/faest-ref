/*
 *  SPDX-License-Identifier: MIT
 */

#include "faest_aes.h"
#include "fields.h"
#include "utils.hpp"
#include "utils.h"
#include "instances.hpp"
#include "randomness.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <cmath>
#include <array>
#include <vector>

namespace aes_ctr_128_tv {
  constexpr std::array<uint8_t, 16> key{
      0x42, 0x13, 0x4f, 0x71, 0x34, 0x89, 0x1b, 0x16,
      0x82, 0xa8, 0xab, 0x56, 0x76, 0x27, 0x30, 0x0c,
  };
  constexpr std::array<uint8_t, 16> in{
      0x7b, 0x60, 0x66, 0xd6, 0x5a, 0x73, 0xd6, 0x00,
      0xa0, 0xae, 0xf0, 0x1c, 0x8b, 0x19, 0x17, 0x40,
  };
  constexpr std::array<uint8_t, 16> out{
      0x89, 0x13, 0xaf, 0x07, 0x31, 0xb5, 0x81, 0xdf,
      0x8a, 0x0a, 0xb5, 0x6e, 0x08, 0x3c, 0x6b, 0x7c,
  };
} // namespace aes_ctr_128_tv

namespace rijndael_em_128_tv {
  constexpr std::array<uint8_t, 16> key = aes_ctr_128_tv::key;
  constexpr std::array<uint8_t, 16> in  = aes_ctr_128_tv::in;
  constexpr std::array<uint8_t, 16> out{
      0xb8, 0xcf, 0x49, 0x82, 0x1b, 0x05, 0xe3, 0xfd,
      0x4c, 0xb6, 0x61, 0x5c, 0x5f, 0x79, 0x77, 0x2f,
  };
} // namespace rijndael_em_128_tv

namespace aes_ctr_192_tv {
  constexpr std::array<uint8_t, 24> key{
      0x7f, 0x64, 0xa4, 0x6b, 0xbd, 0x02, 0x67, 0x2c, 0xed, 0x19, 0xfb, 0x73,
      0x5b, 0xf0, 0x46, 0xaf, 0x23, 0x6e, 0x38, 0x79, 0x85, 0x13, 0x79, 0xd3,
  };
  constexpr std::array<uint8_t, 32> in{
      0x2f, 0xae, 0x1b, 0x7c, 0x4a, 0x8f, 0xb6, 0x1c, 0x15, 0x7c, 0x4d,
      0xe2, 0x9d, 0x35, 0x62, 0x33, 0x63, 0x94, 0x75, 0x39, 0x50, 0x2d,
      0x7e, 0xa5, 0xf7, 0x33, 0xd0, 0xca, 0x3c, 0xc2, 0xb5, 0xd0,
  };
  constexpr std::array<uint8_t, 32> out{
      0x1a, 0xb4, 0x2c, 0x3c, 0xde, 0xb6, 0xb5, 0x08, 0xff, 0xc8, 0x3d,
      0x5b, 0x48, 0x9f, 0x62, 0xca, 0xdd, 0x3f, 0x53, 0x92, 0xbb, 0x4b,
      0x0a, 0xe3, 0xed, 0xf0, 0xe9, 0xe7, 0x0c, 0x4d, 0xb4, 0x2c,
  };
} // namespace aes_ctr_192_tv

namespace rijndael_em_192_tv {
  constexpr std::array<uint8_t, 24> key{
      0x24, 0x18, 0x87, 0x72, 0xc5, 0x1f, 0xbe, 0x52, 0xc0, 0xcd, 0x0b, 0xed,
      0xbe, 0x6a, 0x4c, 0x04, 0xb3, 0x75, 0x89, 0x7d, 0x36, 0x9b, 0x7e, 0x62,
  };
  constexpr std::array<uint8_t, 24> in{
      0xc1, 0xa3, 0xc0, 0x22, 0xe7, 0x18, 0x93, 0x5f, 0x46, 0x63, 0x03, 0x86,
      0xaf, 0xa3, 0xd3, 0xf2, 0xc0, 0x72, 0x0b, 0x10, 0xbf, 0x26, 0x6c, 0x19,
  };
  constexpr std::array<uint8_t, 24> out{
      0xbf, 0x71, 0x16, 0xfd, 0x88, 0xb0, 0x75, 0x2f, 0x06, 0x01, 0xc3, 0x14,
      0xdc, 0xbb, 0xa6, 0x25, 0x6b, 0x8e, 0xc4, 0x5b, 0xd2, 0x4e, 0x10, 0x84,
  };

} // namespace rijndael_em_192_tv

namespace aes_ctr_256_tv {
  constexpr std::array<uint8_t, 32> key{
      0xa9, 0x86, 0x63, 0xac, 0x8a, 0x05, 0x78, 0xbe, 0xd7, 0x2c, 0x80,
      0x91, 0x07, 0x67, 0xce, 0x11, 0xf1, 0x79, 0x59, 0xde, 0x6a, 0x99,
      0xbb, 0xdc, 0x75, 0xb2, 0x04, 0x63, 0x6f, 0x1d, 0xd2, 0x5f,
  };
  constexpr std::array<uint8_t, 32> in{
      0x78, 0xff, 0xd6, 0xd1, 0x95, 0x73, 0xdb, 0x9f, 0xa4, 0x08, 0xe8,
      0xcb, 0x29, 0xf8, 0x2e, 0x27, 0x96, 0xe0, 0x8f, 0x0d, 0xf9, 0xf1,
      0x8a, 0x1d, 0x67, 0xd5, 0xec, 0x22, 0x14, 0x92, 0x32, 0x17,
  };
  constexpr std::array<uint8_t, 32> out{
      0xef, 0x34, 0x98, 0xa1, 0x9e, 0xab, 0x0a, 0x9e, 0xad, 0x99, 0xd7,
      0xf2, 0xe1, 0x68, 0xf5, 0xad, 0x18, 0xdf, 0x02, 0x86, 0x6d, 0xf4,
      0x2b, 0x80, 0x24, 0xe3, 0x9e, 0x24, 0xeb, 0x51, 0x83, 0x67,
  };
} // namespace aes_ctr_256_tv

namespace rijndael_em_256_tv {
  constexpr std::array<uint8_t, 32> key{
      0xc0, 0xcd, 0x0b, 0xed, 0xbe, 0x6a, 0x4c, 0x04, 0xb3, 0x75, 0x89,
      0x7d, 0x36, 0x9b, 0x7e, 0x62, 0xaa, 0x6a, 0x6f, 0x17, 0x13, 0xd2,
      0x7a, 0x71, 0xfe, 0x98, 0x9e, 0x93, 0xdc, 0x79, 0xd2, 0x7d,
  };
  constexpr std::array<uint8_t, 32> in{
      0xc1, 0xa3, 0xc0, 0x22, 0xe7, 0x18, 0x93, 0x5f, 0x46, 0x63, 0x03,
      0x86, 0xaf, 0xa3, 0xd3, 0xf2, 0xc0, 0x72, 0x0b, 0x10, 0xbf, 0x26,
      0x6c, 0x19, 0x24, 0x18, 0x87, 0x72, 0xc5, 0x1f, 0xbe, 0x52,
  };
  constexpr std::array<uint8_t, 32> out{
      0xf4, 0x1b, 0x6b, 0xb2, 0x22, 0x5c, 0xaf, 0xfc, 0x82, 0x31, 0xcc,
      0xe0, 0x21, 0x9e, 0x2c, 0xcc, 0xb0, 0xf8, 0x0e, 0x68, 0x0d, 0xf2,
      0xdf, 0xef, 0xca, 0x1b, 0x96, 0x57, 0x18, 0x42, 0x00, 0x25,
  };
} // namespace rijndael_em_256_tv

static bf128_t* column_to_row_major_and_shrink_V_128(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + 2
  // \lambda rows and \lambda columns storing in row-major order
  bf128_t* new_v = (bf128_t*)malloc((ell + FAEST_128F_LAMBDA * 2) * sizeof(bf128_t));
  // faest_aligned_alloc(BF128_ALIGN, (ell + FAEST_128F_LAMBDA*2) * sizeof(bf128_t));
  for (unsigned int row = 0; row != ell + FAEST_128F_LAMBDA * 2; ++row) {
    uint8_t new_row[BF128_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_128F_LAMBDA; ++column) {
      ptr_set_bit(new_row, column, ptr_get_bit(v[column], row));
    }
    new_v[row] = bf128_load(new_row);
  }
  return new_v;
}

static bf192_t* column_to_row_major_and_shrink_V_192(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf192_t* new_v = (bf192_t*)malloc((ell + FAEST_192F_LAMBDA * 2) * sizeof(bf192_t));
  for (unsigned int row = 0; row != ell + FAEST_192F_LAMBDA * 2; ++row) {
    uint8_t new_row[BF192_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_192F_LAMBDA; ++column) {
      ptr_set_bit(new_row, column, ptr_get_bit(v[column], row));
    }
    new_v[row] = bf192_load(new_row);
  }

  return new_v;
}

static bf256_t* column_to_row_major_and_shrink_V_256(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf256_t* new_v = (bf256_t*)malloc((ell + FAEST_256F_LAMBDA * 2) * sizeof(bf256_t));
  for (unsigned int row = 0; row != ell + FAEST_256F_LAMBDA * 2; ++row) {
    uint8_t new_row[BF256_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_256F_LAMBDA; ++column) {
      ptr_set_bit(new_row, column, ptr_get_bit(v[column], row));
    }
    new_v[row] = bf256_load(new_row);
  }

  return new_v;
}

BOOST_AUTO_TEST_SUITE(test_aes_prove)
/*
BOOST_DATA_TEST_CASE(aes_prove_verify, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t* params = faest_get_paramset(param_id);
    const bool is_em               = faest_is_em(params);
    const unsigned int lambda      = params->faest_param.lambda;
    const unsigned int lambdaBytes = lambda / 8;
    const unsigned int ell = params->faest_param.l;
    const unsigned int ell_hat =
        params->faest_param.l + params->faest_param.lambda * 3 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const unsigned int ell_bytes = (params->faest_param.l + 7) / 8;

    // extended witness
    //std::vector<uint8_t> w;
    std::vector<uint8_t> in;
    std::vector<uint8_t> out;

    if (lambda == 256 && !is_em) {
      for (const auto byte : aes_ctr_256_tv::in) {
          for (size_t bit_i = 0; bit_i < ell; bit_i++) {
              in.push_back((byte >> bit_i) & 1);
          }
      }
      for (const auto byte : aes_ctr_256_tv::out) {
          for (size_t bit_i = 0; bit_i < ell; bit_i++) {
              out.push_back((byte >> bit_i) & 1);
          }
      }
    } else if (lambda == 256 && is_em) {
      for (const auto byte : rijndael_em_256_tv::in) {
          for (size_t bit_i = 0; bit_i < ell; bit_i++) {
              in.push_back((byte >> bit_i) & 1);
          }
      }
      for (const auto byte : rijndael_em_256_tv::out) {
          for (size_t bit_i = 0; bit_i < ell; bit_i++) {
              out.push_back((byte >> bit_i) & 1);
          }
      }
    }
    else {
      return;
    }

    uint8_t* w = aes_extend_witness(in.data(), out.data(), params);
    std::vector<uint8_t> w_bits(ell, 0x00);  // 1 bit in per uint8_t
    for (unsigned int bit_i = 0; bit_i < ell; bit_i++) {
      w_bits[bit_i] = (w[bit_i/8] >> bit_i%8) & 1;
    }

    // prepare vole correlation
    std::vector<uint8_t> delta(lambda / 8, 0);
    for (size_t i = 0; i < lambda / 8; ++i) {
      delta[i] = (uint8_t) i;
    }
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

    std::vector<uint8_t> u_bits(2 * lambda, 0x00); // 1 bit in per uint8_t
    for (unsigned int bit_i = 0; bit_i < 2 * lambda; bit_i++) {
      u_bits[bit_i] = (u[(ell + bit_i) / 8] >> (ell + bit_i) % 8) & 1;
    }
    // masked witness d = u ^ w
    std::vector<uint8_t> d(ell_bytes, 0x13);
    for (size_t i = 0; i < ell_bytes; ++i) {
      d[i] = u[i] ^ w[i];
    }
    std::vector<uint8_t> d_bits(ell, 0x00); // 1 bit in per uint8_t
    for (unsigned int bit_i = 0; bit_i < ell; bit_i++) {
      d_bits[bit_i] = (d[bit_i / 8] >> bit_i % 8) & 1;
    }

    std::vector<uint8_t> chall_2((3 * lambda + 64) / 8, 47);

    std::vector<uint8_t> a0_tilde(lambda / 8, 0);
    std::vector<uint8_t> a1_tilde(lambda / 8, 0);
    std::vector<uint8_t> a2_tilde(lambda / 8, 0);

    bf256_t bf_delta = bf256_load(delta.data());
    bf256_t* q = column_to_row_major_and_shrink_V_256(Q.data(), ell);

    bf256_t* w_tag = column_to_row_major_and_shrink_V_256(V.data(), ell);
    bf256_t* bf_u_bits = (bf256_t*) malloc(2*lambda * sizeof(bf256_t));

    for (unsigned int bit_i = 0; bit_i < 2*lambda; bit_i++) {
      u_bits[bit_i] = (u[(ell + bit_i)/8] >> (ell + bit_i)%8) & 1;
      bf_u_bits[bit_i] = bf256_from_bit(u_bits[bit_i]);
    }
    bf256_t q_star_0 = bf256_sum_poly(q + ell);
    bf256_t bf_u_star_0 = bf256_sum_poly(bf_u_bits);
    bf256_t bf_v_star_0 = bf256_sum_poly(w_tag + ell);

    bf256_t test_v0 = bf256_add(q_star_0, bf256_mul(bf_delta, bf_u_star_0));

    BOOST_TEST(memcmp(&test_v0, &bf_v_star_0, lambdaBytes) == 0);

    printf("testing aes_prove\n");

    aes_prove(a0_tilde.data(), a1_tilde.data(), a2_tilde.data(), w_bits.data(), u_bits.data(),
              V.data(), in.data(), out.data(), chall_2.data(), params);

    uint8_t* recomputed_a0_tilde =
        aes_verify(d_bits.data(), Q.data(), chall_2.data(), delta.data(), a1_tilde.data(),
                   a2_tilde.data(), in.data(), out.data(), params);

    // check that the proof verifies
    printf("FAEST - %s\n", faest_get_param_name(param_id));
    for (size_t i = 0; i < 24; i++) {
      printf("%d-%d ", recomputed_a0_tilde[i], a0_tilde.data()[i]);
    }
    printf("\n");

    BOOST_TEST(memcmp(recomputed_a0_tilde, a0_tilde.data(), lambdaBytes) == 0);
    free(recomputed_a0_tilde);
    free(w);
    free(bf_u_bits);
    free(w_tag);
    free(q);
  }
}
 */

/*
BOOST_DATA_TEST_CASE(aes_prove_verify, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t* params = faest_get_paramset(param_id);
    const bool is_em               = faest_is_em(params);
    const unsigned int lambda      = params->faest_param.lambda;
    const unsigned int lambdaBytes = lambda / 8;
    const unsigned int ell = params->faest_param.l;
    const unsigned int ell_hat =
        params->faest_param.l + params->faest_param.lambda * 3 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const unsigned int ell_bytes = (params->faest_param.l + 7) / 8;

    // extended witness
    //std::vector<uint8_t> w;
    std::vector<uint8_t> in;
    std::vector<uint8_t> out;

    if (lambda == 192 && !is_em) {
      for (const auto byte : aes_ctr_192_tv::in) {
          for (size_t bit_i = 0; bit_i < ell; bit_i++) {
              in.push_back((byte >> bit_i) & 1);
          }
      }
      for (const auto byte : aes_ctr_192_tv::out) {
          for (size_t bit_i = 0; bit_i < ell; bit_i++) {
              out.push_back((byte >> bit_i) & 1);
          }
      }
    } else if (lambda == 192 && is_em) {
      for (const auto byte : rijndael_em_192_tv::in) {
          for (size_t bit_i = 0; bit_i < ell; bit_i++) {
              in.push_back((byte >> bit_i) & 1);
          }
      }
      for (const auto byte : rijndael_em_192_tv::out) {
          for (size_t bit_i = 0; bit_i < ell; bit_i++) {
              out.push_back((byte >> bit_i) & 1);
          }
      }
    }
    else {
      return;
    }

    uint8_t* w = aes_extend_witness(in.data(), out.data(), params);
    std::vector<uint8_t> w_bits(ell, 0x00);  // 1 bit in per uint8_t
    for (unsigned int bit_i = 0; bit_i < ell; bit_i++) {
      w_bits[bit_i] = (w[bit_i/8] >> bit_i%8) & 1;
    }

    // prepare vole correlation
    std::vector<uint8_t> delta(lambda / 8, 0);
    for (size_t i = 0; i < lambda / 8; ++i) {
      delta[i] = (uint8_t) i;
    }
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

    std::vector<uint8_t> u_bits(2 * lambda, 0x00); // 1 bit in per uint8_t
    for (unsigned int bit_i = 0; bit_i < 2 * lambda; bit_i++) {
      u_bits[bit_i] = (u[(ell + bit_i) / 8] >> (ell + bit_i) % 8) & 1;
    }
    // masked witness d = u ^ w
    std::vector<uint8_t> d(ell_bytes, 0x13);
    for (size_t i = 0; i < ell_bytes; ++i) {
      d[i] = u[i] ^ w[i];
    }
    std::vector<uint8_t> d_bits(ell, 0x00); // 1 bit in per uint8_t
    for (unsigned int bit_i = 0; bit_i < ell; bit_i++) {
      d_bits[bit_i] = (d[bit_i / 8] >> bit_i % 8) & 1;
    }

    std::vector<uint8_t> chall_2((3 * lambda + 64) / 8, 47);

    std::vector<uint8_t> a0_tilde(lambda / 8, 0);
    std::vector<uint8_t> a1_tilde(lambda / 8, 0);
    std::vector<uint8_t> a2_tilde(lambda / 8, 0);

    bf192_t bf_delta = bf192_load(delta.data());
    bf192_t* q = column_to_row_major_and_shrink_V_192(Q.data(), ell);

    bf192_t* w_tag = column_to_row_major_and_shrink_V_192(V.data(), ell);
    bf192_t* bf_u_bits = (bf192_t*) malloc(2*lambda * sizeof(bf192_t));

    for (unsigned int bit_i = 0; bit_i < 2*lambda; bit_i++) {
      u_bits[bit_i] = (u[(ell + bit_i)/8] >> (ell + bit_i)%8) & 1;
      bf_u_bits[bit_i] = bf192_from_bit(u_bits[bit_i]);
    }
    bf192_t q_star_0 = bf192_sum_poly(q + ell);
    bf192_t bf_u_star_0 = bf192_sum_poly(bf_u_bits);
    bf192_t bf_v_star_0 = bf192_sum_poly(w_tag + ell);

    bf192_t test_v0 = bf192_add(q_star_0, bf192_mul(bf_delta, bf_u_star_0));

    BOOST_TEST(memcmp(&test_v0, &bf_v_star_0, lambdaBytes) == 0);

    printf("testing aes_prove\n");

    aes_prove(a0_tilde.data(), a1_tilde.data(), a2_tilde.data(), w_bits.data(), u_bits.data(),
              V.data(), in.data(), out.data(), chall_2.data(), params);

    uint8_t* recomputed_a0_tilde =
        aes_verify(d_bits.data(), Q.data(), chall_2.data(), delta.data(), a1_tilde.data(),
                   a2_tilde.data(), in.data(), out.data(), params);

    // check that the proof verifies
    printf("FAEST - %s\n", faest_get_param_name(param_id));
    for (size_t i = 0; i < 24; i++) {
      printf("%d-%d ", recomputed_a0_tilde[i], a0_tilde.data()[i]);
    }
    printf("\n");

    BOOST_TEST(memcmp(recomputed_a0_tilde, a0_tilde.data(), lambdaBytes) == 0);
    free(recomputed_a0_tilde);
    free(w);
    free(bf_u_bits);
    free(w_tag);
    free(q);
  }
}
 */

BOOST_DATA_TEST_CASE(aes_prove_verify, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t* params = faest_get_paramset(param_id);
    const bool is_em               = faest_is_em(params);
    const unsigned int lambda      = params->faest_param.lambda;
    const unsigned int lambdaBytes = lambda / 8;
    const unsigned int ell         = params->faest_param.l;
    const unsigned int ell_hat =
        params->faest_param.l + params->faest_param.lambda * 3 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const unsigned int ell_bytes     = (params->faest_param.l + 7) / 8;

    // extended witness
    // std::vector<uint8_t> w;
    std::vector<uint8_t> in;
    std::vector<uint8_t> out;

    if (lambda == 128 && !is_em) {
      for (const auto byte : aes_ctr_128_tv::in) {
        for (size_t bit_i = 0; bit_i < ell; bit_i++) {
          in.push_back((byte >> bit_i) & 1);
        }
      }
      for (const auto byte : aes_ctr_128_tv::out) {
        for (size_t bit_i = 0; bit_i < ell; bit_i++) {
          out.push_back((byte >> bit_i) & 1);
        }
      }
    } else if (lambda == 128 && is_em) {
      for (const auto byte : rijndael_em_128_tv::in) {
        for (size_t bit_i = 0; bit_i < ell; bit_i++) {
          in.push_back((byte >> bit_i) & 1);
        }
      }
      for (const auto byte : rijndael_em_128_tv::out) {
        for (size_t bit_i = 0; bit_i < ell; bit_i++) {
          out.push_back((byte >> bit_i) & 1);
        }
      }
    } else {
      return;
    }

    uint8_t* w = aes_extend_witness(in.data(), out.data(), params);
    std::vector<uint8_t> w_bits(ell, 0x00); // 1 bit in per uint8_t
    for (unsigned int bit_i = 0; bit_i < ell; bit_i++) {
      w_bits[bit_i] = (w[bit_i / 8] >> bit_i % 8) & 1;
    }

    // if (lambda == 128 && !is_em) {
    //   std::copy(aes_ctr_128_tv::in.begin(), aes_ctr_128_tv::in.end(),
    //             std::back_insert_iterator(in));
    //   std::copy(aes_ctr_128_tv::out.begin(), aes_ctr_128_tv::out.end(),
    //             std::back_insert_iterator(out));
    //   std::copy(aes_ctr_128_tv::expected_extended_witness.begin(),
    //             aes_ctr_128_tv::expected_extended_witness.end(), std::back_insert_iterator(w));
    // } else if (lambda == 128 && is_em) {
    //   std::copy(rijndael_em_128_tv::in.begin(), rijndael_em_128_tv::in.end(),
    //             std::back_insert_iterator(in));
    //   std::copy(rijndael_em_128_tv::out.begin(), rijndael_em_128_tv::out.end(),
    //             std::back_insert_iterator(out));
    //   std::copy(rijndael_em_128_tv::expected_extended_witness.begin(),
    //             rijndael_em_128_tv::expected_extended_witness.end(),
    //             std::back_insert_iterator(w));
    // } else if (lambda == 192 && !is_em) {
    //   std::copy(aes_ctr_192_tv::in.begin(), aes_ctr_192_tv::in.end(),
    //             std::back_insert_iterator(in));
    //   std::copy(aes_ctr_192_tv::out.begin(), aes_ctr_192_tv::out.end(),
    //             std::back_insert_iterator(out));
    //   std::copy(aes_ctr_192_tv::expected_extended_witness.begin(),
    //             aes_ctr_192_tv::expected_extended_witness.end(), std::back_insert_iterator(w));
    // } else if (lambda == 192 && is_em) {
    //   std::copy(rijndael_em_192_tv::in.begin(), rijndael_em_192_tv::in.end(),
    //             std::back_insert_iterator(in));
    //   std::copy(rijndael_em_192_tv::out.begin(), rijndael_em_192_tv::out.end(),
    //             std::back_insert_iterator(out));
    //   std::copy(rijndael_em_192_tv::expected_extended_witness.begin(),
    //             rijndael_em_192_tv::expected_extended_witness.end(),
    //             std::back_insert_iterator(w));
    // } else if (lambda == 256 && !is_em) {
    //   std::copy(aes_ctr_256_tv::in.begin(), aes_ctr_256_tv::in.end(),
    //             std::back_insert_iterator(in));
    //   std::copy(aes_ctr_256_tv::out.begin(), aes_ctr_256_tv::out.end(),
    //             std::back_insert_iterator(out));
    //   std::copy(aes_ctr_256_tv::expected_extended_witness.begin(),
    //             aes_ctr_256_tv::expected_extended_witness.end(), std::back_insert_iterator(w));
    // } else if (lambda == 256 && is_em) {
    //   std::copy(rijndael_em_256_tv::in.begin(), rijndael_em_256_tv::in.end(),
    //             std::back_insert_iterator(in));
    //   std::copy(rijndael_em_256_tv::out.begin(), rijndael_em_256_tv::out.end(),
    //             std::back_insert_iterator(out));
    //   std::copy(rijndael_em_256_tv::expected_extended_witness.begin(),
    //             rijndael_em_256_tv::expected_extended_witness.end(),
    //             std::back_insert_iterator(w));
    // }

    // prepare vole correlation
    std::vector<uint8_t> delta(lambda / 8, 0);
    for (size_t i = 0; i < lambda / 8; ++i) {
      delta[i] = (uint8_t)i;
    }
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

    std::vector<uint8_t> u_bits(2 * lambda, 0x00); // 1 bit in per uint8_t
    for (unsigned int bit_i = 0; bit_i < 2 * lambda; bit_i++) {
      u_bits[bit_i] = (u[(ell + bit_i) / 8] >> (ell + bit_i) % 8) & 1;
    }
    // masked witness d = u ^ w
    std::vector<uint8_t> d(ell_bytes, 0x13);
    for (size_t i = 0; i < ell_bytes; ++i) {
      d[i] = u[i] ^ w[i];
    }
    std::vector<uint8_t> d_bits(ell, 0x00); // 1 bit in per uint8_t
    for (unsigned int bit_i = 0; bit_i < ell; bit_i++) {
      d_bits[bit_i] = (d[bit_i / 8] >> bit_i % 8) & 1;
    }

    std::vector<uint8_t> chall_2((3 * lambda + 64) / 8, 47);

    std::vector<uint8_t> a0_tilde(lambda / 8, 0);
    std::vector<uint8_t> a1_tilde(lambda / 8, 0);
    std::vector<uint8_t> a2_tilde(lambda / 8, 0);

    bf128_t bf_delta = bf128_load(delta.data());
    bf128_t* q       = column_to_row_major_and_shrink_V_128(Q.data(), ell);

    bf128_t* w_tag     = column_to_row_major_and_shrink_V_128(V.data(), ell);
    bf128_t* bf_u_bits = (bf128_t*)malloc(2 * lambda * sizeof(bf128_t));

    for (unsigned int bit_i = 0; bit_i < 2 * lambda; bit_i++) {
      u_bits[bit_i]    = (u[(ell + bit_i) / 8] >> (ell + bit_i) % 8) & 1;
      bf_u_bits[bit_i] = bf128_from_bit(u_bits[bit_i]);
    }
    // verifier Delta, q_star
    // printf("delta\n");
    // print_array<uint8_t>((uint8_t*) &bf_delta, lambdaBytes);
    bf128_t q_star_0 = bf128_sum_poly(q + ell);
    // printf("q_star_0\n");
    // print_array<uint8_t>((uint8_t*) &q_star_0, lambdaBytes);

    // prover u,v (want: q = u*Delta + v)
    bf128_t bf_u_star_0 = bf128_sum_poly(bf_u_bits);
    bf128_t bf_v_star_0 = bf128_sum_poly(w_tag + ell);

    // printf("u_star_0\n");
    //  print_array<uint8_t>((uint8_t*) &bf_u_star_0, lambdaBytes);
    // printf("v_star_0\n");
    //  print_array<uint8_t>((uint8_t*) &bf_v_star_0, lambdaBytes);

    bf128_t test_v0 = bf128_add(q_star_0, bf128_mul(bf_delta, bf_u_star_0));
    // printf("q_star_0 + delta * u_star_0\n");
    //  print_array<uint8_t>((uint8_t*) &test_v0, lambdaBytes);

    BOOST_TEST(memcmp(&test_v0, &bf_v_star_0, lambdaBytes) == 0);

    printf("testing aes_prove\n");

    aes_prove(a0_tilde.data(), a1_tilde.data(), a2_tilde.data(), w_bits.data(), u_bits.data(),
              V.data(), in.data(), out.data(), chall_2.data(), params);

    uint8_t* recomputed_a0_tilde =
        aes_verify(d_bits.data(), Q.data(), chall_2.data(), delta.data(), a1_tilde.data(),
                   a2_tilde.data(), in.data(), out.data(), params);

    // check that the proof verifies
    printf("FAEST - %s\n", faest_get_param_name(param_id));
    for (size_t i = 0; i < 16; i++) {
      printf("%d-%d ", recomputed_a0_tilde[i], a0_tilde.data()[i]);
    }
    printf("\n");

    BOOST_TEST(memcmp(recomputed_a0_tilde, a0_tilde.data(), lambdaBytes) == 0);
    free(recomputed_a0_tilde);
    free(w);
    free(bf_u_bits);
    free(w_tag);
    free(q);
  }
}

BOOST_AUTO_TEST_SUITE_END()
