/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "vole.h"
#include "aes.h"
#include "utils.h"
#include "random_oracle.h"
#include "fields.h"

#include <stdbool.h>
#include <string.h>

static const uint32_t TWEAK_OFFSET = UINT32_C(0x80000000); // 2^31

#define TBL(base, words, i) ((base) + (size_t)(i) * (words))

// to pack bits the way bf2_matrix_mul_tbl expects them to be
static void pack_n_mask_bits(uint8_t* dst, const uint8_t* src, size_t bits, unsigned int mask_idx) {
  memset(dst, 0, (bits + 7) / 8);
  for (size_t bit = 0; bit < bits; ++bit) {
    ptr_set_bit(dst, bit, ptr_get_bit(src + bit * N_MASK_BYTES, mask_idx));
  }
}

#if !defined(FAEST_TESTS)
static
#endif
    unsigned int convert_to_vole(const uint8_t* iv, const uint8_t* sd, bool sd0_bot, unsigned int i,
                                 unsigned int outlen, uint8_t* u, uint8_t* v,
                                 const faest_paramset_t* params) {
  const unsigned int lambda        = params->lambda;
  const unsigned int tau_1         = params->tau1;
  const unsigned int k             = params->k;
  const unsigned int num_instances = bavc_max_node_index(i, tau_1, k);
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int depth         = bavc_max_node_depth(i, tau_1, k);

  // (depth + 1) x num_instances array of outlen; but we only need two rows at a time
  uint8_t* r = calloc(2 * num_instances, outlen);

#define R(row, column) (r + (((row) % 2) * num_instances + (column)) * outlen)
#define V(idx) (v + (idx) * outlen)

  uint32_t tweak = i ^ TWEAK_OFFSET;

  // Step: 2
  if (!sd0_bot) {
    prg(sd, iv, tweak, R(0, 0), lambda, outlen);
  }

  // Step: 3..4
  for (unsigned int j = 1; j < num_instances; ++j) {
    prg(sd + lambda_bytes * j, iv, tweak, R(0, j), lambda, outlen);
  }

  // Step: 5..9
  memset(v, 0, depth * outlen);
  for (unsigned int j = 0; j < depth; j++) {
    unsigned int depthloop = num_instances >> (j + 1);
    for (unsigned int idx = 0; idx < depthloop; idx++) {
      xor_u8_array(V(j), R(j, 2 * idx + 1), V(j), outlen);
      xor_u8_array(R(j, 2 * idx), R(j, 2 * idx + 1), R(j + 1, idx), outlen);
    }
  }
  // Step: 10
  if (!sd0_bot && u != NULL) {
    memcpy(u, R(depth, 0), outlen);
  }
  free(r);
  return depth;
}

#undef R
#undef V

void vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                 const faest_paramset_t* params, bavc_t* bavc, uint8_t* c, uint8_t* c_mult,
                 uint8_t* u, uint8_t** V_dest, uint8_t* u_bar, uint8_t* v_bar) {
  const unsigned int lambda                     = params->lambda;
  const unsigned int lambda_bytes               = lambda / 8;
  const unsigned int ell_hat_bytes              = (ellhat + 7) / 8;
  const unsigned int ell                        = params->ell;
  const unsigned int ell_bytes                  = (ell + 7) / 8;
  const unsigned int tau                        = params->tau;
  const unsigned int tau_1                      = params->tau1;
  const unsigned int k                          = params->k;
  const unsigned int n_mult                     = params->n_mult;
  const unsigned int n_mult_bytes               = (n_mult + 7) / 8;
  const unsigned int w_grind                    = params->w_grind;
  const unsigned int w_grind_bytes              = (w_grind + 7) / 8;
  const unsigned int lambda_minus_w_grind       = lambda - w_grind;
  const unsigned int lambda_minus_w_grind_bytes = ((lambda_minus_w_grind) + 7) / 8;

  // line 1
  bavc_commit(bavc, rootKey, iv, params);

  uint8_t* ui = malloc(tau * ell_hat_bytes);
  assert(ui);

  uint8_t** V = alloc_pointer_array(lambda, ell_hat_bytes);

  // line 3
  unsigned int v_idx = 0;
  uint8_t* sd_i      = bavc->sd;
  for (unsigned int i = 0; i < tau; ++i) {
    v_idx += convert_to_vole(iv, sd_i, false, i, ell_hat_bytes, ui + i * ell_hat_bytes, V[v_idx],
                             params);
    sd_i += lambda_bytes * bavc_max_node_index(i, tau_1, k);
  }

  // line 6
  memcpy(u, ui, ell_bytes);

  // line 8
  for (unsigned int i = 1; i < tau; i++) {
    xor_u8_array(u, ui + i * ell_hat_bytes, c + (i - 1) * ell_bytes, ell_bytes);
  }

  // line 10
  uint8_t* u_hi_all = malloc((N_MASK * w_grind + 7) / 8);
  uint8_t* u_hi     = calloc(N_MASK, w_grind_bytes);
  prg(rootKey, iv, TWEAK_OFFSET - 1, u_hi_all, lambda, (N_MASK * w_grind + 7) / 8);

  for (unsigned n_mask_idx = 0; n_mask_idx < N_MASK; n_mask_idx++) {
    for (unsigned w_grind_idx = 0; w_grind_idx < w_grind; w_grind_idx++) {
      ptr_set_bit(u_hi + n_mask_idx * w_grind_bytes, w_grind_idx,
                  ptr_get_bit(u_hi_all, n_mask_idx * w_grind + w_grind_idx));
    }
  }

  // line 11
  uint8_t* u_low   = calloc(N_MASK, lambda_bytes);
  uint8_t* r_tilde = calloc(N_MASK, lambda_minus_w_grind_bytes);
  for (unsigned m = 0; m < N_MASK; ++m) {
    unsigned int ell_m   = ell + m * bavc_max_node_depth(0, tau_1, k);
    unsigned int bit_idx = 0;

    // line 13
    for (unsigned t = 0; t < tau; t++) {
      const unsigned int depth = bavc_max_node_depth(t, tau_1, k);
      for (unsigned int b = 0; b < depth; ++b, ++bit_idx) {
        ptr_set_bit(u_low, m * lambda + bit_idx, ptr_get_bit(ui + t * ell_hat_bytes, ell_m + b));
      }
    }

    // line 14
    uint8_t first_prod[MAX_LAMBDA_BYTES] = {0};
    bf2_matrix_mul_tbl(first_prod, u_low + m * lambda_bytes, params->W_CRT, lambda_minus_w_grind,
                       lambda, params->w_crt_words);
    uint8_t second_prod[MAX_LAMBDA_BYTES] = {0};
    bf2_poly_mul(second_prod, &u_hi[m * w_grind_bytes], w_grind, params->M_TREE,
                 lambda_minus_w_grind + 1);
    xor_u8_array(first_prod, second_prod, u_bar + m * lambda_bytes, lambda_bytes);

    // line 17
    bit_idx = 0;
    v_idx   = 0;
    for (unsigned int tau_idx = 0; tau_idx < tau; tau_idx++) {
      uint8_t deg = bavc_max_node_depth(tau_idx, tau_1, k);

      uint8_t acc[(MAX_DEPTH * 2 + 7) / 8] = {0};
      for (unsigned int d = 0; d < deg; ++d) {
        for (unsigned int t = 0; t < deg; ++t, ++v_idx) {
          ptr_xor_bit(acc, t + d, ptr_get_bit(V[v_idx], ell_m + d));
        }
        v_idx -= deg;
      }
      v_idx += deg;

      uint8_t reduced_out[(MAX_DEPTH * 2 + 7) / 8] = {0};
      bf2_poly_reduce(reduced_out, acc, deg * 2, &params->TREE_MODULI[tau_idx], deg + 1);

      for (unsigned int b = 0; b < deg; ++b, ++bit_idx) {
        ptr_set_bit(r_tilde + m * lambda_minus_w_grind_bytes, bit_idx, ptr_get_bit(reduced_out, b));
      }
    }
  }

  // line 19
  // a is a copy of V, so directly access V

  // line 21
  uint8_t* v_tilde = malloc(n_mult * N_MASK_BYTES);
  for (unsigned e = 0; e < n_mult; e++) {
    uint8_t L_e_zero[MAX_LAMBDA_BYTES] = {0};
    uint8_t L_e_one[MAX_LAMBDA_BYTES];
    uint8_t h_e_zero[N_MASK_BYTES];
    uint8_t h_e_one[N_MASK_BYTES];

    for (unsigned j = 0; j < lambda_minus_w_grind; j++) {
      if (ptr_u64_get_bit(TBL(params->G, params->g_words, e), j) == 1) {
        xor_u8_array(L_e_zero, V[j], L_e_zero, lambda_bytes);
      }
    }
    xor_u8_array(L_e_zero, u, L_e_one, lambda_bytes);

    // line 23
    H5(h_e_zero, iv, e, L_e_zero, lambda);
    H5(h_e_one, iv, e, L_e_one, lambda);

    // line 26
    for (unsigned m = 0; m < N_MASK; m++) {
      uint8_t dot_product_bit = 0;
      for (unsigned i = 0; i < lambda; i++) {
        uint8_t bit_a = ptr_u64_get_bit(TBL(params->F, params->f_words, e), i);
        uint8_t bit_b = ptr_get_bit(u_bar + m * lambda_bytes, i);

        dot_product_bit ^= bit_a & bit_b;
      }

      ptr_set_bit(v_tilde + e * N_MASK_BYTES, m, ptr_get_bit(h_e_zero, m));
      uint8_t bit_a = ptr_get_bit(v_tilde + e * N_MASK_BYTES, m);
      uint8_t bit_b = ptr_get_bit(h_e_one, m);

      ptr_set_bit(c_mult + m * n_mult_bytes, e, bit_a ^ bit_b ^ dot_product_bit);
    }
  }

  // line 29
  uint8_t* gate_input = malloc(n_mult_bytes);
  for (unsigned int m = 0; m < N_MASK; m++) {
    uint8_t first_prod[MAX_LAMBDA_BYTES] = {0};
    bf2_matrix_mul_tbl(first_prod, r_tilde + m * lambda_minus_w_grind_bytes, params->W_TREE,
                       lambda_minus_w_grind, lambda, params->w_tree_words);

    uint8_t second_prod[MAX_LAMBDA_BYTES] = {0};
    pack_n_mask_bits(gate_input, v_tilde, n_mult, m);
    bf2_matrix_mul_tbl(second_prod, gate_input, params->W_GATE, n_mult, lambda,
                       params->w_gate_words);
    xor_u8_array(first_prod, second_prod, v_bar + m * lambda_bytes, lambda_bytes);
  }

  // line 30
  for (unsigned int ell_idx = 0; ell_idx < ell; ell_idx++) {
    for (unsigned int col = 0; col < lambda_minus_w_grind; col++) {
      uint8_t acc = 0;
      for (unsigned int row = 0; row < lambda_minus_w_grind; row++) {
        acc ^= ptr_get_bit(V[row], ell_idx) &
               ptr_u64_get_bit(TBL(params->W_CRT, params->w_crt_words, col), row);
      }

      ptr_xor_bit(V_dest[col], ell_idx, acc);
    }
  }
  free(gate_input);
  free_pointer_array(&V);

  free(ui);
  free(u_hi_all);
  free(u_hi);
  free(u_low);
  free(r_tilde);
  free(v_tilde);
}

bool vole_reconstruct(uint8_t* com, uint8_t** Q_dest, const uint8_t* iv, const uint8_t* chall_3,
                      const uint8_t* decom_i, const uint8_t* c, const uint8_t* c_mult,
                      uint8_t* q_bar, uint8_t* Delta, unsigned int ellhat,
                      const faest_paramset_t* params) {
  const unsigned int lambda                     = params->lambda;
  const unsigned int lambda_bytes               = lambda / 8;
  const unsigned int ellhat_bytes               = (ellhat + 7) / 8;
  const unsigned int ell                        = params->ell;
  const unsigned int ell_bytes                  = (ell + 7) / 8;
  const unsigned int tau                        = params->tau;
  const unsigned int tau_1                      = params->tau1;
  const unsigned int L                          = params->L;
  const unsigned int k                          = params->k;
  const unsigned int n_mult                     = params->n_mult;
  const unsigned int n_mult_bytes               = (n_mult + 7) / 8;
  const unsigned int w_grind                    = params->w_grind;
  const unsigned int lambda_minus_w_grind       = lambda - w_grind;
  const unsigned int lambda_minus_w_grind_bytes = ((lambda_minus_w_grind) + 7) / 8;

  // line 1
  uint16_t i_delta[MAX_TAU];
  if (!decode_all_chall_3(i_delta, chall_3, params)) {
    return false;
  }

  // line 14
  // Delta = W_crt * Dp
  // Delta has lambda output bits
  bf2_matrix_mul_tbl(Delta, chall_3, params->W_CRT, lambda_minus_w_grind, lambda,
                     params->w_crt_words);
  bavc_rec_t bavc_rec;
  bavc_rec.h = com;
  bavc_rec.s = malloc((L - tau) * lambda_bytes);

  // line 3
  if (!bavc_rec.s || !bavc_reconstruct(&bavc_rec, decom_i, i_delta, iv, params)) {
    free(bavc_rec.s);
    return false;
  }

  uint8_t* sd   = malloc((1 << k) * lambda_bytes);
  uint8_t* qtmp = malloc(MAX_DEPTH * ellhat_bytes);
  assert(sd);
  assert(qtmp);

  uint8_t** Q = alloc_pointer_array(lambda, ellhat_bytes);

  unsigned int q_idx = 0;
  uint8_t* sd_i      = bavc_rec.s;
  for (unsigned int i = 0; i < tau; i++) {
    const unsigned int Ni = bavc_max_node_index(i, tau_1, k);

    for (unsigned int j = 0; j < Ni; j++) {
      if (j < i_delta[i]) {
        memcpy(sd + (j ^ i_delta[i]) * lambda_bytes, sd_i + lambda_bytes * j, lambda_bytes);
      } else if (j > i_delta[i]) {
        memcpy(sd + (j ^ i_delta[i]) * lambda_bytes, sd_i + lambda_bytes * (j - 1), lambda_bytes);
      }
    }

    // line 8
    const unsigned int ki = convert_to_vole(iv, sd, true, i, ellhat_bytes, NULL, qtmp, params);

    if (i == 0) {
      memcpy(Q[q_idx], qtmp, ellhat_bytes * ki);
      q_idx += ki;
    } else {
      for (unsigned int d = 0; d < ki; ++d, ++q_idx) {
        memcpy(Q[q_idx], qtmp + d * ellhat_bytes, ellhat_bytes);

        masked_xor_u8_array(Q[q_idx], c + (i - 1) * ell_bytes, Q[q_idx], (i_delta[i] >> d) & 1,
                            ell_bytes);
      }
    }
    sd_i += lambda_bytes * (Ni - 1);
  }
  free(sd);
  sd = NULL;
  free(bavc_rec.s);
  bavc_rec.s = sd_i = NULL;
  free(qtmp);
  qtmp = NULL;

  // line 13
  // delta_prime can be directly read from chall_3

  // line 19
  uint8_t* r_tilde_prime = calloc(N_MASK, lambda_minus_w_grind_bytes);
  // for scan-build
  assert(lambda_minus_w_grind_bytes > 0);
  for (unsigned int m = 0; m < N_MASK; m++) {
    unsigned int ell_m   = ell + m * bavc_max_node_depth(0, tau_1, k);
    unsigned int bit_idx = 0;
    q_idx                = 0;
    for (unsigned tau_idx = 0; tau_idx < tau; tau_idx++) {
      const uint8_t deg = bavc_max_node_depth(tau_idx, tau_1, k);

      uint8_t acc[(MAX_DEPTH * 2 + 7) / 8] = {0};
      for (unsigned d = 0; d < deg; d++) {
        for (unsigned t = 0; t < deg; t++, ++q_idx) {
          uint8_t bit = ptr_get_bit(Q[q_idx], ell_m + d);
          ptr_xor_bit(acc, t + d, bit);
        }
        q_idx -= deg;
      }
      q_idx += deg;

      uint8_t reduced_out[(MAX_DEPTH * 2 + 7) / 8] = {0};
      bf2_poly_reduce(reduced_out, acc, deg * 2, &params->TREE_MODULI[tau_idx], deg + 1);

      for (unsigned int b = 0; b < deg; ++b, ++bit_idx) {
        ptr_set_bit(r_tilde_prime + m * lambda_minus_w_grind_bytes, bit_idx,
                    ptr_get_bit(reduced_out, b));
      }
    }
  }

  // line 21
  // a_prime is a direct copy of Q, so use Q instead

  uint8_t* q_tilde = calloc(n_mult, N_MASK_BYTES);
  for (unsigned e = 0; e < n_mult; e++) {
    uint8_t L_e_prime[MAX_LAMBDA_BYTES] = {0};
    uint8_t h_e[N_MASK_BYTES];

    // line 24
    uint8_t gamma_e = 0;
    for (unsigned j = 0; j < lambda_minus_w_grind; j++) {
      gamma_e ^= ptr_u64_get_bit(TBL(params->G, params->g_words, e), j) & ptr_get_bit(chall_3, j);

      if (ptr_u64_get_bit(TBL(params->G, params->g_words, e), j) == 1) {
        xor_u8_array(L_e_prime, Q[j], L_e_prime, lambda_bytes);
      }
    }

    // line 25
    H5(h_e, iv, e, L_e_prime, lambda);

    // line 27
    for (unsigned m = 0; m < N_MASK; m++) {
      uint8_t sum = ptr_get_bit(h_e, m) ^ (gamma_e & ptr_get_bit(c_mult + m * n_mult_bytes, e));
      ptr_set_bit(q_tilde + e * N_MASK_BYTES, m, sum);
    }
  }

  // line 30
  uint8_t* gate_input = malloc(n_mult_bytes);
  for (unsigned int m = 0; m < N_MASK; m++) {
    uint8_t first_prod[MAX_LAMBDA_BYTES] = {0};
    bf2_matrix_mul_tbl(first_prod, r_tilde_prime + m * lambda_minus_w_grind_bytes, params->W_TREE,
                       lambda_minus_w_grind, lambda, params->w_tree_words);

    uint8_t second_prod[MAX_LAMBDA_BYTES] = {0};
    pack_n_mask_bits(gate_input, q_tilde, n_mult, m);
    bf2_matrix_mul_tbl(second_prod, gate_input, params->W_GATE, n_mult, lambda,
                       params->w_gate_words);
    xor_u8_array(first_prod, second_prod, q_bar + m * lambda_bytes, lambda_bytes);
  }
  free(gate_input);
  gate_input = NULL;
  free(r_tilde_prime);
  r_tilde_prime = NULL;
  free(q_tilde);
  q_tilde = NULL;

  // line 31
  // TODO refactor this
  for (unsigned int ell_idx = 0; ell_idx < ell; ell_idx++) {
    for (unsigned int col = 0; col < lambda_minus_w_grind; col++) {
      uint8_t acc = 0;
      for (unsigned int row = 0; row < lambda_minus_w_grind; row++) {
        acc ^= ptr_get_bit(Q[row], ell_idx) &
               ptr_u64_get_bit(TBL(params->W_CRT, params->w_crt_words, col), row);
      }

      ptr_xor_bit(Q_dest[col], ell_idx, acc);
    }
  }

  free_pointer_array(&Q);
  return true;
}
