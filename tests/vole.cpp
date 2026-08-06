/*
 *  SPDX-License-Identifier: MIT
 */

#include "vole.h"
#include "instances.hpp"
#include "randomness.h"
#include "universal_hashing.h"
#include "utils.hpp"
#include "vole_tvs.hpp"

#include "tables/tables.h"
// #include "tables/tables_128s.h"
// #include "tables/tables_128f.h"

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
  constexpr std::array<uint8_t, 16> iv{};

  constexpr unsigned int max_tries = 1000;
} // namespace

BOOST_AUTO_TEST_SUITE(vole)

BOOST_DATA_TEST_CASE(vole_commit_verify, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const auto params                 = faest_get_paramset(param_id);
    const unsigned int lambda         = params->lambda;
    const unsigned int lambda_bytes   = lambda / 8;
    const unsigned int ell            = params->ell;
    const unsigned int ell_bytes      = (ell + 7) / 8;
    const unsigned int tau           = params->tau;
    const unsigned int k              = params->k;
    const unsigned int tau1           = params->tau1;
    const unsigned int d0            = bavc_max_node_depth(0, tau1, k);
    const unsigned int n_mask        = params->n_mask;
    const unsigned int ell_hat       = ell + n_mask * d0;
    const unsigned int ell_hat_bytes  = (ell_hat + 7) / 8;
    const auto com_size               = (faest_is_em(params) ? 2 : 3) * lambda_bytes;
    const auto n_mask_bytes           = (n_mask + 7) / 8;
    const auto n_mult                 = params->n_mult;

    bavc_t bavc_com;

    std::vector<uint8_t> chal, c, c_mult, u, u_bar, v_bar, q_bar, delta_prime, q_storage, v_storage;
    chal.resize(lambda_bytes);
    c.resize((params->tau - 1) * ell_bytes);
    c_mult.resize(n_mult * n_mask_bytes);
    u.resize(ell_bytes * params->tau);
    u_bar.resize(n_mask * lambda_bytes);
    v_bar.resize(n_mask * lambda_bytes);
    q_bar.resize(n_mask * lambda_bytes);
    delta_prime.resize(lambda_bytes);

    std::vector<uint8_t*> q, v;
    q.resize(lambda);
    v.resize(lambda);

    q_storage.resize(lambda * ell_hat_bytes);
    v_storage.resize(lambda * ell_hat_bytes);

    q[0] = q_storage.data();
    v[0] = v_storage.data();
    for (unsigned int i = 1; i < lambda; ++i) {
      q[i] = q[0] + i * ell_hat_bytes;
      v[i] = v[0] + i * ell_hat_bytes;
    }

    vole_commit(rootKey.data(), iv.data(), ell_hat, params, &bavc_com, c.data(), c_mult.data(),
                u.data(), v.data(), u_bar.data(), v_bar.data());

    std::vector<uint8_t> hcom{bavc_com.h, bavc_com.h + lambda_bytes * 2};

    bool tested = false;
    for (unsigned int tries = 0; !tested && tries != max_tries; ++tries) {

      rand_bytes(chal.data(), chal.size());
      for (unsigned int i = lambda - params->w_grind; i != lambda; ++i) {
        ptr_set_bit(chal.data(), i, 0);
      }
      // TODO: Remove
      // NOTE: Just for debugging with the faest-128f chall TV instead of the random chall above
      // static const char* chall_hex = "7155b8f0de65bed193b8615bcde68900";
      // for (unsigned int i = 0; i < lambda_bytes; ++i) {
      //   unsigned int byte;
      //   std::sscanf(chall_hex + 2 * i, "%2x", &byte);
      //   chal[i] = static_cast<uint8_t>(byte);
      // }

      std::vector<uint16_t> i_delta;
      i_delta.resize(params->tau);

      std::vector<uint8_t> decom_i;
      decom_i.resize(com_size * params->tau + params->T_open * lambda_bytes);

      BOOST_TEST(decode_all_chall_3(i_delta.data(), chal.data(), params));
      if (!bavc_open(decom_i.data(), &bavc_com, i_delta.data(), params)) {
        continue;
      }
      tested = true;

      std::vector<uint8_t> hcom_rec;
      hcom_rec.resize(lambda_bytes * 2);
      BOOST_TEST(vole_reconstruct(hcom_rec.data(), q.data(), iv.data(), chal.data(), decom_i.data(),
                                  c.data(), c_mult.data(), q_bar.data(), delta_prime.data(), ell_hat, params));
      BOOST_TEST(hcom == hcom_rec);
    }
    BOOST_TEST(tested);
    bavc_clear(&bavc_com);

  }
}

BOOST_DATA_TEST_CASE(test_convert_to_vole, all_parameters, param_id) {
  std::mt19937_64 rd;
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const faest_paramset_t params    = *faest_get_paramset(param_id);
    const unsigned int lambda        = params.lambda;
    const unsigned int lambda_bytes  = lambda / 8;
    const unsigned int tau           = params.tau;
    const unsigned int k              = params.k;
    const unsigned int tau1           = params.tau1;
    const unsigned int d0            = bavc_max_node_depth(0, tau1, k);
    const unsigned int n_mask        = params.n_mask;
    const unsigned int ell           = params.ell;
    const unsigned int ell_hat       = ell + n_mask * d0;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const unsigned int max_depth     = k;
    const unsigned int max_nodes     = 1 << max_depth;

    std::vector<uint8_t> sd, u, v, q;
    sd.resize(max_nodes * lambda_bytes);
    rand_bytes(sd.data(), sd.size());
    u.resize(ell_hat_bytes);
    v.resize(ell_hat_bytes * max_depth);
    q.resize(ell_hat_bytes * max_depth);

    for (unsigned int i = 0; i != tau; ++i) {
      std::uniform_int_distribution<> distribution{
          0, static_cast<int>(bavc_max_node_index(i, params.tau1, max_depth)) - 1};
      const unsigned int idx = distribution(rd);

      unsigned int depth = bavc_max_node_depth(i, params.tau1, max_depth);
      unsigned int nodes = 1 << depth;

      convert_to_vole(iv.data(), sd.data(), false, i, ell_hat_bytes, u.data(), v.data(), &params);

      std::vector<uint8_t> sdprime;
      sdprime.resize(nodes * lambda_bytes, 0);
      for (unsigned int j = 1; j != nodes; ++j) {
        std::copy(sd.begin() + (j ^ idx) * lambda_bytes,
                  sd.begin() + ((j ^ idx) + 1) * lambda_bytes, &sdprime[j * lambda_bytes]);
      }

      convert_to_vole(iv.data(), sdprime.data(), true, i, ell_hat_bytes, nullptr, q.data(),
                      &params);

      for (unsigned int j = 0; j != depth; ++j) {
        if (idx & (1 << j)) {
          for (unsigned int inner = 0; inner != ell_hat_bytes; ++inner) {
            q[j * ell_hat_bytes + inner] ^= u[inner];
          }
        }
      }
      BOOST_TEST(q == v);
    }
  }
}

namespace {
  template <std::size_t HSize, std::size_t CSize>
  void test_tv(const faest_paramset_t* params, 
              const std::array<uint8_t, CSize>& challenge,
              const std::array<uint8_t, HSize>& expected_h,
              const std::array<uint8_t, 64>& expected_hashed_c,
              const std::array<uint8_t, 64>& expected_hashed_u,
              const std::array<uint8_t, 64>& expected_hashed_c_mult,
              const std::array<uint8_t, 64>& expected_hashed_v,
              const std::array<uint8_t, 64>& expected_hashed_barU,
              const std::array<uint8_t, 64>& expected_hashed_barV,
              const std::array<uint8_t, 64>& expected_hashed_q,
              const std::array<uint8_t, 64>& expected_hashed_barQ) {
    const unsigned int lambda        = params->lambda;
    const unsigned int lambda_bytes  = lambda / 8;
    const unsigned int ell            = params->ell;
    const unsigned int ell_bytes      = (ell + 7) / 8;
    const unsigned int tau           = params->tau;
    const unsigned int k              = params->k;
    const unsigned int tau1           = params->tau1;
    const unsigned int d0            = bavc_max_node_depth(0, tau1, k);
    const unsigned int n_mask        = params->n_mask;
    const unsigned int ell_hat       = ell + n_mask * d0;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    const auto com_size              = (faest_is_em(params) ? 2 : 3) * lambda_bytes;
    const auto n_mask_bytes           = (n_mask + 7) / 8;
    const auto n_mult                 = params->n_mult;
    const unsigned int w_grind      = params->w_grind;
    const unsigned int w_grind_bytes = (w_grind + 7) / 8;
    const unsigned int lambda_minus_w_grind = lambda - w_grind;
    const unsigned int lambda_minus_w_grind_bytes = ((lambda_minus_w_grind) + 7 ) / 8;

    const faest_tables_t* T = faest_get_tables(params->id);

    // NOTE: This TV only works for the deg-7 vole commit implementation with the below parameters
    switch (lambda) {
    case 256:
      assert(ell == 2208 || ell == 1792);
      assert(n_mask == 9);
      break;
    case 192:
      assert(ell == 1728 || ell == 1152);
      assert(n_mask == 9);
      break;
    default:
      assert(ell == 960 || ell == 640);
      assert(n_mask == 9);
      break;
    }
    

    bavc_t bavc_com;

    std::vector<uint8_t> chal, c, c_mult, decom_i, u, u_bar, v_bar, q_bar, delta_prime, q_storage, v_storage;
    chal.resize(lambda_bytes);
    c.resize((params->tau - 1) * ell_bytes);
    c_mult.resize(n_mult * n_mask_bytes);
    decom_i.resize(com_size * params->tau + params->T_open * lambda_bytes);
    u.resize(ell_bytes);
    u_bar.resize(n_mask * lambda_bytes);
    v_bar.resize(n_mask * lambda_bytes);
    q_bar.resize(n_mask * lambda_bytes);
    delta_prime.resize(lambda_bytes);

    std::vector<uint8_t*> q, v;
    q.resize(lambda);
    v.resize(lambda);

    q_storage.resize(lambda * ell_hat_bytes);
    v_storage.resize(lambda * ell_hat_bytes);

    q[0] = q_storage.data();
    v[0] = v_storage.data();
    for (unsigned int i = 1; i < lambda; ++i) {
      q[i] = q[0] + i * ell_hat_bytes;
      v[i] = v[0] + i * ell_hat_bytes;
    }

    vole_commit(rootKey.data(), iv.data(), ell_hat, params, &bavc_com, c.data(), c_mult.data(),
                u.data(), v.data(), u_bar.data(), v_bar.data());

    std::vector<uint8_t> hcom{bavc_com.h, bavc_com.h + lambda_bytes * 2},
        expected_h_vec{expected_h.begin(), expected_h.end()};
    BOOST_TEST(hcom == expected_h_vec);

    // print_named_array("challenge", "uint8_t", challenge.data(), challenge.size());

    // print_named_array("c", "uint8_t", c);
    // print_named_array("expected_hashed_c", "uint8_t", expected_hashed_c.data(), expected_hashed_c.size());
    BOOST_TEST(expected_hashed_c == hash_array(c));

    // print_named_array("u", "uint8_t", u);
    // print_named_array("expected_hashed_u", "uint8_t", expected_hashed_u.data(), expected_hashed_u.size());
    BOOST_TEST(expected_hashed_u == hash_array(u));

    BOOST_TEST(expected_hashed_barU == hash_array(u_bar));
    BOOST_TEST(expected_hashed_barV == hash_array(v_bar));

    // NOTE: *** Claude GENERATED CODE
    // This alligns the c_mult to how it is in the python code (swapping the rows and the bits) for the matching TVs
    size_t nb        = (n_mult + 7) / 8;      // 40, matches Python
    size_t packed_len = n_mask * nb;          // 360
    std::vector<uint8_t> c_mult_packed(packed_len, 0);
    for (size_t i = 0; i < n_mask; ++i) {
        for (size_t j = 0; j < n_mult; ++j) {
          uint16_t word = c_mult[j*2] | (uint16_t)(c_mult[j*2 + 1] << 8);
          uint8_t  bit  = (word >> i) & 1u;      // bit i of gate j's word
          size_t   pos  = i * nb * 8 + j;        // dst: row byte-aligned, bit j LSB-first
          c_mult_packed[pos >> 3] |= (uint8_t)(bit << (pos & 7));
        }
    }
    // ***

    // print_named_array("hashed_c_mult", "uint8_t", hash_array(c_mult_packed));
    BOOST_TEST(expected_hashed_c_mult == hash_array(c_mult_packed));

    // NOTE: *** Claude GENERATED CODE
    // This does the row_to_coloumn_major transformation
    size_t ncols        = ell;
    // size_t lambda_bytes = lambda / 8;
    std::vector<uint8_t> mV_bytes(ncols * lambda_bytes, 0);
    for (size_t c = 0; c < ncols; ++c) {
      for (size_t r = 0; r < lambda; ++r) {
          const uint8_t* srow = v_storage.data() + r * ell_hat_bytes;
          uint8_t bit = (srow[c >> 3] >> (c & 7)) & 1u;
          size_t  pos = c * lambda + r;
          mV_bytes[pos >> 3] |= (uint8_t)(bit << (pos & 7));
      }
    }
    // fprintf(stderr, "cc len: %zu\n", mV_bytes.size());
    // for (size_t k = 0; k < 32; ++k) fprintf(stderr, "%02x", mV_bytes[k]);
    // fprintf(stderr, "\n... tail: ");
    // for (size_t k = mV_bytes.size()-32; k < mV_bytes.size(); ++k) fprintf(stderr, "%02x", mV_bytes[k]);
    // fprintf(stderr, "\n");
    // print_named_array("v_storage", "uint8_t", mV_bytes.data(), 32);
    // print_named_array("v_storage", "uint8_t", mV_bytes.data(), 140);
    // ***
    // print_named_array("v_storage", "uint8_t", hash_array(mV_bytes));
    // print_named_array("expected_hashed_v", "uint8_t", expected_hashed_v.data(), expected_hashed_v.size());
    BOOST_TEST(expected_hashed_v == hash_array(mV_bytes));

    std::vector<uint16_t> i_delta;
    i_delta.resize(params->tau);
    BOOST_TEST(decode_all_chall_3(i_delta.data(), challenge.data(), params));
    BOOST_TEST(bavc_open(decom_i.data(), &bavc_com, i_delta.data(), params));

    std::vector<uint8_t> hcom_rec;
    hcom_rec.resize(lambda_bytes * 2);
    BOOST_TEST(vole_reconstruct(hcom_rec.data(), q.data(), iv.data(), challenge.data(), decom_i.data(),
                                  c.data(), c_mult.data(), q_bar.data(), delta_prime.data(), ell_hat, params));
    BOOST_TEST(hcom_rec == expected_h_vec);

    // NOTE: *** Claude GENERATED CODE
    ncols        = ell;
    // size_t lambda_bytes = lambda / 8;
    std::vector<uint8_t> mQ_bytes(ncols * lambda_bytes, 0);
    for (size_t c = 0; c < ncols; ++c) {
      for (size_t r = 0; r < lambda; ++r) {
          const uint8_t* srow = q_storage.data() + r * ell_hat_bytes;
          uint8_t bit = (srow[c >> 3] >> (c & 7)) & 1u;
          size_t  pos = c * lambda + r;
          mQ_bytes[pos >> 3] |= (uint8_t)(bit << (pos & 7));
      }
    }
    // print_named_array("q_storage", "uint8_t", mQ_bytes.data(), 140);
    // ***
    BOOST_TEST(expected_hashed_q == hash_array(mQ_bytes));

    BOOST_TEST(expected_hashed_barQ == hash_array(q_bar));

    bavc_clear(&bavc_com);

    uint8_t* Dp = (uint8_t*)malloc(lambda_minus_w_grind_bytes);
    memset(Dp, 0, lambda_minus_w_grind_bytes);

    unsigned int off = 0;
    for (unsigned int i = 0; i < tau; ++i) {
      const uint32_t deg =
          bavc_max_node_depth(i, tau1, k); // == tree_deg[i]
      const auto delta = i_delta[i];                          // == I[i]
      for (unsigned int t = 0; t < deg; ++t) {
        xor_bit_to_pt(Dp, off + t, (delta >> t) & 1);         // Dp |= ((I[i]>>t)&1) << (off+t)
      }
      off += deg;
    }

    // --- Delta = W_crt * Dp  (crt_lift), same matvec pattern as your v_row_maj example ---
    uint8_t* Delta = (uint8_t*)malloc(lambda_bytes);
    memset(Delta, 0, lambda_bytes);

    for (unsigned int col = 0; col < lambda; ++col) {          // Delta has lambda output bits
      uint8_t acc = 0;
      for (unsigned int row = 0; row < lambda_minus_w_grind; ++row) {   // Dp has n_delta bits
        acc ^= ptr_get_bit(Dp, row)
            // & ptr_get_bit((uint8_t*)&FAEST_128S_W_CRT[col][row / 64], row % 64);
            & ptr_get_bit((uint8_t*)&TBL_U8(T->W_CRT, T->w_crt_words, col)[row / 64], row % 64);
      }
      xor_bit_to_pt(Delta, col, acc);
    }

    assert(off == lambda_minus_w_grind);

    for (unsigned int row = 0; row < ell; ++row) {          // Python: for j in range(ell)
      // rhs = mV[row] ^ (Delta if (u>>row)&1 else 0)
      for (unsigned int b = 0; b < lambda_bytes; ++b) {
        uint8_t rhs_b = mV_bytes[row * lambda_bytes + b];                           // however mV rows are laid out
        if ((ptr_get_bit(u.data(), row)))                        // (u >> row) & 1
          rhs_b ^= Delta[b];
        BOOST_TEST(mQ_bytes[row * lambda_bytes + b] == rhs_b);
      }
    }
  }
} // namespace

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(vole_tv)

BOOST_AUTO_TEST_CASE(vole_tv_128f) {
  vole::test_tv(faest_get_paramset(FAEST_128F),
          vole_tvs::FAEST_128F::chall,
          vole_tvs::FAEST_128F::h,
          vole_tvs::FAEST_128F::hashed_c,
          vole_tvs::FAEST_128F::hashed_u,
          vole_tvs::FAEST_128F::hashed_c_mult,
          vole_tvs::FAEST_128F::hashed_v,
          vole_tvs::FAEST_128F::hashed_barU,
          vole_tvs::FAEST_128F::hashed_barV,
          vole_tvs::FAEST_128F::hashed_q,
          vole_tvs::FAEST_128F::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_128s) {
  vole::test_tv(faest_get_paramset(FAEST_128S),
          vole_tvs::FAEST_128S::chall,
          vole_tvs::FAEST_128S::h,
          vole_tvs::FAEST_128S::hashed_c,
          vole_tvs::FAEST_128S::hashed_u,
          vole_tvs::FAEST_128S::hashed_c_mult,
          vole_tvs::FAEST_128S::hashed_v,
          vole_tvs::FAEST_128S::hashed_barU,
          vole_tvs::FAEST_128S::hashed_barV,
          vole_tvs::FAEST_128S::hashed_q,
          vole_tvs::FAEST_128S::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_em_128f) {
  vole::test_tv(faest_get_paramset(FAEST_EM_128F),
          vole_tvs::FAEST_EM_128F::chall,
          vole_tvs::FAEST_EM_128F::h,
          vole_tvs::FAEST_EM_128F::hashed_c,
          vole_tvs::FAEST_EM_128F::hashed_u,
          vole_tvs::FAEST_EM_128F::hashed_c_mult,
          vole_tvs::FAEST_EM_128F::hashed_v,
          vole_tvs::FAEST_EM_128F::hashed_barU,
          vole_tvs::FAEST_EM_128F::hashed_barV,
          vole_tvs::FAEST_EM_128F::hashed_q,
          vole_tvs::FAEST_EM_128F::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_em_128s) {
  vole::test_tv(faest_get_paramset(FAEST_EM_128S),
          vole_tvs::FAEST_EM_128S::chall,
          vole_tvs::FAEST_EM_128S::h,
          vole_tvs::FAEST_EM_128S::hashed_c,
          vole_tvs::FAEST_EM_128S::hashed_u,
          vole_tvs::FAEST_EM_128S::hashed_c_mult,
          vole_tvs::FAEST_EM_128S::hashed_v,
          vole_tvs::FAEST_EM_128S::hashed_barU,
          vole_tvs::FAEST_EM_128S::hashed_barV,
          vole_tvs::FAEST_EM_128S::hashed_q,
          vole_tvs::FAEST_EM_128S::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_192f) {
  vole::test_tv(faest_get_paramset(FAEST_192F),
          vole_tvs::FAEST_192F::chall,
          vole_tvs::FAEST_192F::h,
          vole_tvs::FAEST_192F::hashed_c,
          vole_tvs::FAEST_192F::hashed_u,
          vole_tvs::FAEST_192F::hashed_c_mult,
          vole_tvs::FAEST_192F::hashed_v,
          vole_tvs::FAEST_192F::hashed_barU,
          vole_tvs::FAEST_192F::hashed_barV,
          vole_tvs::FAEST_192F::hashed_q,
          vole_tvs::FAEST_192F::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_192s) {
  vole::test_tv(faest_get_paramset(FAEST_192S),
          vole_tvs::FAEST_192S::chall,
          vole_tvs::FAEST_192S::h,
          vole_tvs::FAEST_192S::hashed_c,
          vole_tvs::FAEST_192S::hashed_u,
          vole_tvs::FAEST_192S::hashed_c_mult,
          vole_tvs::FAEST_192S::hashed_v,
          vole_tvs::FAEST_192S::hashed_barU,
          vole_tvs::FAEST_192S::hashed_barV,
          vole_tvs::FAEST_192S::hashed_q,
          vole_tvs::FAEST_192S::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_em_192f) {
  vole::test_tv(faest_get_paramset(FAEST_EM_192F),
          vole_tvs::FAEST_EM_192F::chall,
          vole_tvs::FAEST_EM_192F::h,
          vole_tvs::FAEST_EM_192F::hashed_c,
          vole_tvs::FAEST_EM_192F::hashed_u,
          vole_tvs::FAEST_EM_192F::hashed_c_mult,
          vole_tvs::FAEST_EM_192F::hashed_v,
          vole_tvs::FAEST_EM_192F::hashed_barU,
          vole_tvs::FAEST_EM_192F::hashed_barV,
          vole_tvs::FAEST_EM_192F::hashed_q,
          vole_tvs::FAEST_EM_192F::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_em_192s) {
  vole::test_tv(faest_get_paramset(FAEST_EM_192S),
          vole_tvs::FAEST_EM_192S::chall,
          vole_tvs::FAEST_EM_192S::h,
          vole_tvs::FAEST_EM_192S::hashed_c,
          vole_tvs::FAEST_EM_192S::hashed_u,
          vole_tvs::FAEST_EM_192S::hashed_c_mult,
          vole_tvs::FAEST_EM_192S::hashed_v,
          vole_tvs::FAEST_EM_192S::hashed_barU,
          vole_tvs::FAEST_EM_192S::hashed_barV,
          vole_tvs::FAEST_EM_192S::hashed_q,
          vole_tvs::FAEST_EM_192S::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_256f) {
  vole::test_tv(faest_get_paramset(FAEST_256F),
          vole_tvs::FAEST_256F::chall,
          vole_tvs::FAEST_256F::h,
          vole_tvs::FAEST_256F::hashed_c,
          vole_tvs::FAEST_256F::hashed_u,
          vole_tvs::FAEST_256F::hashed_c_mult,
          vole_tvs::FAEST_256F::hashed_v,
          vole_tvs::FAEST_256F::hashed_barU,
          vole_tvs::FAEST_256F::hashed_barV,
          vole_tvs::FAEST_256F::hashed_q,
          vole_tvs::FAEST_256F::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_256s) {
  vole::test_tv(faest_get_paramset(FAEST_256S),
          vole_tvs::FAEST_256S::chall,
          vole_tvs::FAEST_256S::h,
          vole_tvs::FAEST_256S::hashed_c,
          vole_tvs::FAEST_256S::hashed_u,
          vole_tvs::FAEST_256S::hashed_c_mult,
          vole_tvs::FAEST_256S::hashed_v,
          vole_tvs::FAEST_256S::hashed_barU,
          vole_tvs::FAEST_256S::hashed_barV,
          vole_tvs::FAEST_256S::hashed_q,
          vole_tvs::FAEST_256S::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_em_256f) {
  vole::test_tv(faest_get_paramset(FAEST_EM_256F),
          vole_tvs::FAEST_EM_256F::chall,
          vole_tvs::FAEST_EM_256F::h,
          vole_tvs::FAEST_EM_256F::hashed_c,
          vole_tvs::FAEST_EM_256F::hashed_u,
          vole_tvs::FAEST_EM_256F::hashed_c_mult,
          vole_tvs::FAEST_EM_256F::hashed_v,
          vole_tvs::FAEST_EM_256F::hashed_barU,
          vole_tvs::FAEST_EM_256F::hashed_barV,
          vole_tvs::FAEST_EM_256F::hashed_q,
          vole_tvs::FAEST_EM_256F::hashed_barQ);
}

BOOST_AUTO_TEST_CASE(vole_tv_em_256s) {
  vole::test_tv(faest_get_paramset(FAEST_EM_256S),
          vole_tvs::FAEST_EM_256S::chall,
          vole_tvs::FAEST_EM_256S::h,
          vole_tvs::FAEST_EM_256S::hashed_c,
          vole_tvs::FAEST_EM_256S::hashed_u,
          vole_tvs::FAEST_EM_256S::hashed_c_mult,
          vole_tvs::FAEST_EM_256S::hashed_v,
          vole_tvs::FAEST_EM_256S::hashed_barU,
          vole_tvs::FAEST_EM_256S::hashed_barV,
          vole_tvs::FAEST_EM_256S::hashed_q,
          vole_tvs::FAEST_EM_256S::hashed_barQ);
}

BOOST_AUTO_TEST_SUITE_END()
