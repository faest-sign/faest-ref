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
#include "tables/tables.h"

#include <stdbool.h>
#include <string.h>

static const uint32_t TWEAK_OFFSET = UINT32_C(0x80000000); // 2^31


#if !defined(FAEST_TESTS)
static
#endif
    unsigned int
    convert_to_vole(const uint8_t* iv, const uint8_t* sd, bool sd0_bot, unsigned int i,
                    unsigned int outlen, uint8_t* u, uint8_t* v, const faest_paramset_t* params) {
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

void vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                 const faest_paramset_t* params, bavc_t* bavc, uint8_t* c, uint8_t* c_mult,
                 uint8_t* u, uint8_t** v, uint8_t* u_bar, uint8_t* v_bar) {
  const unsigned int lambda       = params->lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ell_hat_bytes = (ellhat + 7) / 8;
  const unsigned int ell          = params->ell;
  const unsigned int ell_bytes    = (ell + 7) / 8;
  const unsigned int tau          = params->tau;
  const unsigned int tau_1        = params->tau1;
  const unsigned int k            = params->k;
  const unsigned int n_mask       = params->n_mask;
  const unsigned int n_mask_bytes = (n_mask + 7) / 8;
  const unsigned int n_mult       = params->n_mult;
  const unsigned int w_grind      = params->w_grind;
  const unsigned int w_grind_bytes = (w_grind + 7) / 8;
  const unsigned int lambda_minus_w_grind = lambda - w_grind;
  const unsigned int lambda_minus_w_grind_bytes = ((lambda_minus_w_grind) + 7 ) / 8;

  const faest_tables_t* T = faest_get_tables(params->id);

  bavc_commit(bavc, rootKey, iv, params);

  uint8_t* ui = malloc(tau * ell_hat_bytes);
  assert(ui);

  // line 4
  unsigned int v_idx = 0;
  uint8_t* sd_i      = bavc->sd;
  for (unsigned int i = 0; i < tau; ++i) {

    v_idx +=
        convert_to_vole(iv, sd_i, false, i, ell_hat_bytes, ui + i * ell_hat_bytes, v[v_idx], params);
    sd_i += lambda_bytes * bavc_max_node_index(i, tau_1, k);
  }

  // line 7
  memcpy(u, ui, ell_bytes);

  // line 8
  for (unsigned int i = 1; i < tau; i++) {
    xor_u8_array(u, ui + i * ell_hat_bytes, c + (i - 1) * ell_bytes, ell_bytes);
  }

  uint8_t** v_row_maj = malloc(ellhat * sizeof(uint8_t*));
  for (unsigned i = 0; i < ellhat; i++) {
    v_row_maj[i] = calloc(lambda_minus_w_grind_bytes, 1);
  }
  column_to_row_major(v, v_row_maj, lambda_minus_w_grind, ellhat);

  // line 11
  uint8_t* u_hi_all = malloc((n_mask * w_grind + 7) / 8);
  uint8_t* u_hi = malloc(n_mask * w_grind_bytes);
  memset(u_hi, 0, n_mask * w_grind_bytes);
  prg(rootKey, iv, ((unsigned int)2 << 30) - 1, u_hi_all, lambda, (n_mask * w_grind + 7) / 8);

  for (unsigned n_mask_idx = 0; n_mask_idx < n_mask; n_mask_idx++) {
    for (unsigned w_grind_idx = 0; w_grind_idx < w_grind; w_grind_idx++) {
      ptr_set_bit(u_hi + n_mask_idx * w_grind_bytes, w_grind_idx, ptr_get_bit(u_hi_all, n_mask_idx * w_grind + w_grind_idx));
    }
  }

  // line 14
  uint8_t* u_low = malloc(n_mask * lambda_bytes);
  memset(u_low, 0, n_mask * lambda_bytes);

  uint8_t* r_tilde = malloc(n_mask * lambda_minus_w_grind_bytes);
  memset(r_tilde, 0, n_mask * lambda_minus_w_grind_bytes);
  for (unsigned m = 0; m < n_mask; m++) {

    unsigned int ell_m = ell + m * bavc_max_node_depth(0, tau_1, k);
    unsigned int bit_idx = 0;

    for (unsigned t = 0; t < tau; t++) {
      for (unsigned b = 0; b < bavc_max_node_depth(t, tau_1, k); b++) {
        ptr_set_bit(u_low, m * lambda + bit_idx, ptr_get_bit(ui + t * ell_hat_bytes, ell_m + b));
        bit_idx++;
      }
    }

    // line 15
    uint8_t* first_prod = malloc(lambda_bytes);
    memset(first_prod, 0, lambda_bytes);
    for (unsigned int col = 0; col < lambda; col++) {
      uint8_t sum = 0; 
      for (unsigned int row = 0; row < lambda; row++) {
        uint8_t bit_a = ptr_get_bit((uint8_t*)&TBL_U8(T->W_CRT, T->w_crt_words, col)[row / 64], row%64);
        uint8_t bit_b = ptr_get_bit(u_low, m * lambda + row);
        sum ^= bit_a & bit_b;

      }
      ptr_set_bit(first_prod, col, sum);
    }
    uint8_t* second_prod = malloc(lambda_bytes);
    gf2_poly_mul_ct((uint8_t*)&T->M_TREE[0], lambda_minus_w_grind + 1, (uint8_t*)&u_hi[m * w_grind_bytes], 
                        w_grind, second_prod);
    
    for (unsigned int i = 0; i < lambda_bytes; i++) {
      u_bar[m * lambda_bytes + i] = first_prod[i] ^ second_prod[i];
    }

    // line 18
    bit_idx = 0;
    v_idx = 0;
    for (unsigned tau_idx = 0; tau_idx < tau; tau_idx++) {
      uint8_t deg = bavc_max_node_depth(tau_idx, tau_1, k);

      uint8_t* acc = malloc((deg*2 + 7) / 8);
      memset(acc, 0, (deg*2 + 7) / 8);

      for (unsigned d = 0; d < deg; d++) {
        for (unsigned t = 0; t < deg; t++) {
          uint8_t bit = ptr_get_bit(v_row_maj[ell_m + d], v_idx); 
          xor_bit_to_pt(acc, t + d, bit);
          v_idx++;
        }
        v_idx -= deg;
      }
      v_idx += deg;

      uint8_t* reduced_out = malloc((deg*2 + 7) / 8);
      gf2_poly_reduce_ct(acc, deg*2, (const uint8_t*)&T->TREE_MODULI[tau_idx],
                deg + 1, reduced_out);

      for (unsigned int b = 0; b < deg; b++) {
        ptr_set_bit(r_tilde + m * lambda_minus_w_grind_bytes, bit_idx, ptr_get_bit(reduced_out, b));
        bit_idx++;
      }

      free(acc);
      free(reduced_out);
    }
  
    free(second_prod);
    free(first_prod);
  }

  // line 20
  uint8_t* a = malloc(lambda_minus_w_grind * lambda_bytes);
  for (unsigned int i = 0; i < lambda_minus_w_grind; i++) {
    memcpy(a + i * lambda_bytes, v[i], lambda_bytes);
  }

  // line 22
  uint8_t* v_tilde = malloc(n_mult * n_mask_bytes);
  for (unsigned e = 0; e < n_mult; e++) {

    uint8_t* L_e_zero = malloc(lambda_bytes);
    uint8_t* L_e_one = malloc(lambda_bytes);
    uint8_t* h_e_zero = malloc(n_mask_bytes);
    uint8_t* h_e_one = malloc(n_mask_bytes);

    memset(L_e_zero, 0, lambda_bytes);
    memcpy(L_e_one, u, lambda_bytes);
    for (unsigned j = 0; j < lambda_minus_w_grind; j++) {
      if (ptr_get_bit((const uint8_t*)&TBL_U8(T->G, T->g_words, e)[j/64], j%64) == 1) {
        for (unsigned b = 0; b < lambda_bytes; b++) {
          L_e_zero[b] ^= a[j * lambda_bytes + b];
        }
      }
    }
    for (unsigned b = 0; b < lambda_bytes; b++) {
      L_e_one[b] ^= L_e_zero[b];
    }

    // line 24
    H5(iv, e, L_e_zero, h_e_zero, lambda_bytes, n_mask, lambda);
    H5(iv, e, L_e_one, h_e_one, lambda_bytes, n_mask, lambda);

    // line 27
    for (unsigned m = 0; m < n_mask; m++) {

      uint8_t dot_product_bit = 0;
      for (unsigned i = 0; i < lambda; i++) {
        
        uint8_t bit_a = ptr_get_bit((uint8_t*)&TBL_U8(T->F, T->f_words, e)[i / 64], i%64);
        uint8_t bit_b = ptr_get_bit(u_bar + m * lambda_bytes, i);

        dot_product_bit ^= bit_a & bit_b;
      }

      ptr_set_bit(v_tilde + e * n_mask_bytes, m, ptr_get_bit(h_e_zero, m));
      uint8_t bit_a = ptr_get_bit(v_tilde + e * n_mask_bytes, m);
      uint8_t bit_b = ptr_get_bit(h_e_one, m);

      ptr_set_bit(c_mult + e * n_mask_bytes, m, bit_a ^ bit_b ^ dot_product_bit);
    }

    free(L_e_zero);
    free(L_e_one);
    free(h_e_zero);
    free(h_e_one);
  }

  // line 30
  for (unsigned int m = 0; m < n_mask; m++) {
    
    uint8_t* first_prod = malloc(lambda_bytes);
    memset(first_prod, 0, lambda_bytes);
    for (unsigned int row = 0; row < lambda; row++) {
      for (unsigned int col = 0; col < lambda_minus_w_grind; col++) {

        uint8_t bit_a = ptr_get_bit((uint8_t*)&TBL_U8(T->W_TREE, T->w_tree_words, row)[col / 64], col%64);
        uint8_t bit_b = ptr_get_bit(r_tilde + m * lambda_minus_w_grind_bytes, col);

        xor_bit_to_pt(first_prod, row, bit_a & bit_b);

      }
    }

    uint8_t* second_prod = malloc(lambda_bytes);
    memset(second_prod, 0, lambda_bytes);

    for (unsigned int row = 0; row < lambda; row++) {
      for (unsigned int col = 0; col < n_mult; col++) {   
        
        uint8_t bit_a = ptr_get_bit((uint8_t*)&TBL_U8(T->W_GATE, T->w_gate_words, row)[col / 64], col % 64);
        uint8_t bit_b = ptr_get_bit(v_tilde + col * n_mask_bytes, m);

        xor_bit_to_pt(second_prod, row, bit_a & bit_b);

      }
    }

    for (unsigned int i = 0; i < lambda_bytes; i++) {
      v_bar[m * lambda_bytes + i] = first_prod[i] ^ second_prod[i];
    }

    free(first_prod);
    free(second_prod);

  }

  // line 32
  for (unsigned int i = 0; i < lambda_minus_w_grind; i++) {
    memset(v[i], 0, ell_hat_bytes);
  }
  uint8_t** v_row_maj_new = malloc(ell * sizeof(uint8_t*));
  for (unsigned i = 0; i < ell; i++) {
    v_row_maj_new[i] = malloc(lambda_minus_w_grind_bytes);
    memset(v_row_maj_new[i], 0, lambda_minus_w_grind_bytes);
  }

  for (unsigned int ell_idx = 0; ell_idx < ell; ell_idx++) {

    for (unsigned int col = 0; col < lambda_minus_w_grind; col++) {

      uint8_t acc = 0;
      for (unsigned int row = 0; row < lambda_minus_w_grind; row++) {

        acc ^= ptr_get_bit(v_row_maj[ell_idx], row) 
          & ptr_get_bit((uint8_t*)&TBL_U8(T->W_CRT, T->w_crt_words, col)[row / 64], row % 64);

      }

      xor_bit_to_pt(v_row_maj_new[ell_idx], col, acc);
    }

  }


  column_to_row_major(v_row_maj_new, v, ell, lambda_minus_w_grind); 


  free(ui);
  free(a);
  for (unsigned i = 0; i < ellhat; i++) { free(v_row_maj[i]); }
  free(v_row_maj); 
  free(u_hi_all);
  free(u_hi);
  free(u_low);
  free(r_tilde);
  free(v_tilde);
  for (unsigned i = 0; i < ell; i++) { free(v_row_maj_new[i]); }
  free(v_row_maj_new);  
  
}

bool vole_reconstruct(uint8_t* com, uint8_t** q, const uint8_t* iv, const uint8_t* chall_3,
                      const uint8_t* decom_i, const uint8_t* c, uint8_t* c_mult,
                      uint8_t* q_bar, uint8_t* Delta, unsigned int ellhat, const faest_paramset_t* params) {
  const unsigned int lambda       = params->lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ellhat_bytes = (ellhat + 7) / 8;
  const unsigned int ell          = params->ell;
  const unsigned int ell_bytes    = (ell + 7) / 8;
  const unsigned int tau          = params->tau;
  const unsigned int tau_1        = params->tau1;
  const unsigned int L            = params->L;
  const unsigned int k            = params->k;
  const unsigned int n_mask       = params->n_mask;
  const unsigned int n_mask_bytes = (n_mask + 7) / 8;
  const unsigned int n_mult       = params->n_mult;
  const unsigned int w_grind      = params->w_grind;
  const unsigned int lambda_minus_w_grind = lambda - w_grind;
  const unsigned int lambda_minus_w_grind_bytes = ((lambda_minus_w_grind) + 7 ) / 8;

  const faest_tables_t* T = faest_get_tables(params->id);

  // line 1
  uint16_t i_delta[MAX_TAU];
  if (!decode_all_chall_3(i_delta, chall_3, params)) {
    return false;
  }

  uint8_t* Dp = (uint8_t*)malloc(lambda_minus_w_grind_bytes);
  memset(Dp, 0, lambda_minus_w_grind_bytes);
  unsigned int off = 0;
  for (unsigned int i = 0; i < tau; ++i) {
    const uint32_t deg =
        bavc_max_node_depth(i, tau_1, k); // == tree_deg[i]
    const uint16_t delta = i_delta[i];                          // == I[i]
    for (unsigned int t = 0; t < deg; ++t) {
      xor_bit_to_pt(Dp, off + t, (delta >> t) & 1);         // Dp |= ((I[i]>>t)&1) << (off+t)
    }
    off += deg;
  }
  // --- Delta = W_crt * Dp  (crt_lift), same matvec pattern as your v_row_maj example ---
  // uint8_t* Delta = (uint8_t*)malloc(lambda_bytes);
  // memset(Delta, 0, lambda_bytes);
  for (unsigned int col = 0; col < lambda; ++col) {          // Delta has lambda output bits
    uint8_t acc = 0;
    for (unsigned int row = 0; row < lambda_minus_w_grind; ++row) {   // Dp has n_delta bits
      acc ^= ptr_get_bit(Dp, row)
          // & ptr_get_bit((uint8_t*)&FAEST_128S_W_CRT[col][row / 64], row % 64);
          & ptr_get_bit((uint8_t*)&TBL_U8(T->W_CRT, T->w_crt_words, col)[row / 64], row % 64);
    }
    xor_bit_to_pt(Delta, col, acc);
  }

  bavc_rec_t bavc_rec;
  bavc_rec.h = com;
  bavc_rec.s = malloc((L - tau) * lambda_bytes);
  assert(bavc_rec.s);

  // line 3
  if (!bavc_reconstruct(&bavc_rec, decom_i, i_delta, iv, params)) {
    free(bavc_rec.s);
    return false;
  }

  uint8_t* sd   = malloc((1 << k) * lambda_bytes);
  uint8_t* qtmp = malloc(MAX_DEPTH * ellhat_bytes);
  assert(sd);
  assert(qtmp);

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

    // line 7
    const unsigned int ki = convert_to_vole(iv, sd, true, i, ellhat_bytes, NULL, qtmp, params);

    if (i == 0) {
      memcpy(q[q_idx], qtmp, ellhat_bytes * ki);
      q_idx += ki;
    } else {
      for (unsigned int d = 0; d < ki; ++d, ++q_idx) {

        memcpy(q[q_idx], qtmp + d * ellhat_bytes, ellhat_bytes);

        masked_xor_u8_array(q[q_idx], c + (i - 1) * ell_bytes, q[q_idx],
                            (i_delta[i] >> d) & 1, ell_bytes);
      }
    }
    sd_i += lambda_bytes * (Ni - 1);
  }

  // line 12
  uint8_t* delta_prime = malloc(lambda_bytes);
  memset(delta_prime, 0, lambda_bytes);
  unsigned int delta_prime_idx = 0;
  for (unsigned int tau_idx = 0; tau_idx < tau; tau_idx++) {
    uint8_t deg = bavc_max_node_depth(tau_idx, tau_1, k);
    for (unsigned d = 0; d < deg; d++) {
      ptr_set_bit(delta_prime, delta_prime_idx, ptr_get_bit((uint8_t*)&i_delta[tau_idx], d));
      delta_prime_idx++;
    }
  }

  for (unsigned r = 0; r < lambda_minus_w_grind; ++r) {
    assert(ptr_get_bit(Dp, r) == ptr_get_bit(delta_prime, r));
  }

  // Switching from coloumn major to row major order
  uint8_t** q_row_maj = malloc(ellhat * sizeof(uint8_t*));
  for (unsigned i = 0; i < ellhat; i++) {
    q_row_maj[i] = calloc(lambda_minus_w_grind_bytes, 1);
  }
  column_to_row_major(q, q_row_maj, lambda_minus_w_grind, ellhat);

  // line 14
  uint8_t* r_tilde_prime = malloc(n_mask * lambda_minus_w_grind_bytes);
  memset(r_tilde_prime, 0, n_mask * lambda_minus_w_grind_bytes);
  for (unsigned int m = 0; m < n_mask; m++) {

    unsigned int ell_m = ell + m * bavc_max_node_depth(0, tau_1, k);

    unsigned int bit_idx = 0;
    q_idx = 0;
    for (unsigned tau_idx = 0; tau_idx < tau; tau_idx++) {

      uint8_t deg = bavc_max_node_depth(tau_idx, tau_1, k);

      uint8_t* acc = malloc((deg*2 + 7) / 8);
      memset(acc, 0, (deg*2 + 7) / 8);

      for (unsigned d = 0; d < deg; d++) {
        for (unsigned t = 0; t < deg; t++) {
          uint8_t bit = ptr_get_bit(q_row_maj[ell_m + d], q_idx);
          xor_bit_to_pt(acc, t + d, bit);
          q_idx++;
        }
        q_idx -= deg;
      }
      q_idx += deg;

      // line 17
      uint8_t* reduced_out = malloc((deg*2 + 7) / 8);
      gf2_poly_reduce_ct(acc, deg*2, (const uint8_t*)&T->TREE_MODULI[tau_idx], 
                deg + 1, reduced_out);

      for (unsigned int b = 0; b < deg; b++) {
        ptr_set_bit(r_tilde_prime + m * lambda_minus_w_grind_bytes, bit_idx, ptr_get_bit(reduced_out, b));
        bit_idx++;
      }

      free(acc);
      free(reduced_out);

    }
  }

  // line 19
  uint8_t* a = malloc(lambda_minus_w_grind * lambda_bytes);
  for (unsigned int i = 0; i < lambda_minus_w_grind; i++) {
    memcpy(a + i * lambda_bytes, q[i], lambda_bytes);
  }

  // line 21
  uint8_t* q_tilde = malloc(n_mult * n_mask_bytes);
  memset(q_tilde, 0, n_mult * n_mask_bytes);
  for (unsigned e = 0; e < n_mult; e++) {

    uint8_t* L_e_prime = malloc(lambda_bytes);
    memset(L_e_prime, 0, lambda_bytes);
    uint8_t* h_e = malloc(n_mask_bytes);

    uint8_t gamma_e = 0;
    for (unsigned j = 0; j < lambda_minus_w_grind; j++) {
      gamma_e ^= ptr_get_bit((uint8_t*)&TBL_U8(T->G, T->g_words, e)[j/64], j%64) 
              & ptr_get_bit(delta_prime, j); 

      if (ptr_get_bit((uint8_t*)&TBL_U8(T->G, T->g_words, e)[j/64], j%64) == 1) {
        for (unsigned b = 0; b < lambda_bytes; b++) {
          L_e_prime[b] ^= a[j * lambda_bytes + b];
        }
      }
    }

    // line 23
    H5(iv, e, L_e_prime, h_e, lambda_bytes, n_mask, lambda);
    
    // line 25
    for (unsigned m = 0; m < n_mask; m++) {
        
      uint8_t sum = ptr_get_bit(h_e, m) ^ (gamma_e & ptr_get_bit(c_mult + e * n_mask_bytes, m));
      ptr_set_bit(q_tilde + e * n_mask_bytes, m, sum);
    }

    free(L_e_prime);
    free(h_e);
  }

  // line 28
  for (unsigned int m = 0; m < n_mask; m++) {
    
    uint8_t* first_prod = malloc(lambda_bytes);
    memset(first_prod, 0, lambda_bytes);
    for (unsigned int row = 0; row < lambda; row++) {
      for (unsigned int col = 0; col < lambda_minus_w_grind; col++) {

        uint8_t bit_a = ptr_get_bit((uint8_t*)&TBL_U8(T->W_TREE, T->w_tree_words, row)[col / 64], col%64);
        uint8_t bit_b = ptr_get_bit(r_tilde_prime + m * lambda_minus_w_grind_bytes, col);

        xor_bit_to_pt(first_prod, row, bit_a & bit_b);

      }
    }

    uint8_t* second_prod = malloc(lambda_bytes);
    memset(second_prod, 0, lambda_bytes);

    for (unsigned int row = 0; row < lambda; row++) {
      for (unsigned int col = 0; col < n_mult; col++) {   
        
        uint8_t bit_a = ptr_get_bit((uint8_t*)&TBL_U8(T->W_GATE, T->w_gate_words, row)[col / 64], col % 64);
        uint8_t bit_b = ptr_get_bit(q_tilde + col * n_mask_bytes, m);

        xor_bit_to_pt(second_prod, row, bit_a & bit_b);

      }
    }

    for (unsigned int i = 0; i < lambda_bytes; i++) {
      q_bar[m * lambda_bytes + i] = first_prod[i] ^ second_prod[i];
    }

    free(first_prod);
    free(second_prod);

  }

  for (unsigned int i = 0; i < lambda_minus_w_grind; i++) {
    memset(q[i], 0, ellhat_bytes);
  }
  uint8_t** q_row_maj_new = malloc(ell * sizeof(uint8_t*));
  for (unsigned i = 0; i < ell; i++) {
    q_row_maj_new[i] = malloc(lambda_minus_w_grind_bytes);
    memset(q_row_maj_new[i], 0, lambda_minus_w_grind_bytes);
  }

  // line 29
  for (unsigned int ell_idx = 0; ell_idx < ell; ell_idx++) {

    for (unsigned int col = 0; col < lambda_minus_w_grind; col++) {

      uint8_t acc = 0;
      for (unsigned int row = 0; row < lambda_minus_w_grind; row++) {

        acc ^= ptr_get_bit(q_row_maj[ell_idx], row) 
          & ptr_get_bit((uint8_t*)&TBL_U8(T->W_CRT, T->w_crt_words, col)[row / 64], row % 64);

      }

      xor_bit_to_pt(q_row_maj_new[ell_idx], col, acc);
    }
  }

  column_to_row_major(q_row_maj_new, q, ell, lambda_minus_w_grind); 


  free(qtmp);
  free(sd);
  free(bavc_rec.s);
  for (unsigned i = 0; i < ellhat; i++) { free(q_row_maj[i]); }
  free(q_row_maj);  
  free(r_tilde_prime);
  free(a);
  free(q_tilde);
  for (unsigned i = 0; i < ell; i++) { free(q_row_maj_new[i]); }
  free(q_row_maj_new); 


  return true;
}
