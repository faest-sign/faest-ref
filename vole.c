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
#include "tables/tables_128f.h"

#include <stdbool.h>
#include <string.h>

static const uint32_t TWEAK_OFFSET = UINT32_C(0x80000000); // 2^31

// TODO: For now putting this here
static inline void xor_bit_to_pt(uint8_t *p, size_t i, uint8_t v) {
    size_t byte = i / 8;
    p[byte] ^= (uint8_t)((v & 1u) << (i % 8));
}
void gf2_poly_mul_ct(const uint8_t *a, size_t a_bits,
                     const uint8_t *b, size_t b_bits,
                     uint8_t *out) {
    size_t out_bits  = a_bits + b_bits - 1;
    size_t out_bytes = (out_bits + 7) / 8;
    memset(out, 0, out_bytes);
    for (size_t i = 0; i < a_bits; i++) {
        uint8_t ai = get_bit_from_pt(a, i);
        uint8_t mask = (uint8_t)(0u - ai);
        for (size_t j = 0; j < b_bits; j++) {
            uint8_t bj = get_bit_from_pt(b, j);
            size_t  k  = i + j;
            out[k / 8] ^= (uint8_t)((bj & mask) << (k % 8));
        }
    }
}

void gf2_poly_reduce_ct(const uint8_t *a, size_t a_bits,
                        const uint8_t *m, size_t m_bits,
                        uint8_t *out)
{
  size_t a_bytes = (a_bits + 7) / 8;
  size_t out_bytes = (m_bits + 7) / 8;

  memcpy(out, a, a_bytes);

  size_t top_pad = out_bytes * 8;
  for (size_t i = a_bits; i < top_pad; i++) {
    set_bit_to_pt(out, i, 0);
  }

  size_t mdeg = m_bits - 1;

  for (size_t i = a_bits; i-- > mdeg;) {
      uint8_t lead = get_bit_from_pt(out, i);
      uint8_t mask = (uint8_t)(0u - lead);

      size_t shift = i - mdeg;

      for (size_t j = 0; j < m_bits; j++) {
          uint8_t mj = get_bit_from_pt(m, j);
          size_t  k  = shift + j;
          out[k / 8] ^= (uint8_t)((mj & mask) << (k % 8));
      }

      if (i == 0) break;
  }
}
static inline int highest_set_bit_pos_u64(uint64_t x) {
    int pos = 0;
    while (x) { pos++; x >>= 1; }
    return pos;
}

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
  const unsigned int ellhat_bytes = (ellhat + 7) / 8;
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

  // print_u8_array("rookey", rootKey, 5);
  // print_u8_array("iv", iv, 5);

  bavc_commit(bavc, rootKey, iv, params);

  // print_u8_array("bavc->h", bavc->h, 5);

  uint8_t* ui = malloc(tau * ellhat_bytes);
  assert(ui);

  unsigned int v_idx = 0;
  uint8_t* sd_i      = bavc->sd;
  for (unsigned int i = 0; i < tau; ++i) {

    v_idx +=
        convert_to_vole(iv, sd_i, false, i, ellhat_bytes, ui + i * ellhat_bytes, v[v_idx], params);
    sd_i += lambda_bytes * bavc_max_node_index(i, tau_1, k);
  }

  // line 7
  memcpy(u, ui, ell_bytes);

  // print_u8_array("u", u, 32);

  // line 8-9
  for (unsigned int i = 1; i < tau; i++) {
    xor_u8_array(u, ui + i * ellhat_bytes, c + (i - 1) * ell_bytes, ell_bytes);
  }

  // print_u8_array("c", c, 5);

  // line 11
  // uint8_t* A = malloc(lambda_minus_w_grind * lambda_bytes);
  // for (unsigned int i = 0; i < lambda_minus_w_grind; i++) {
  //   memcpy(A + i * lambda_bytes, v[i], lambda_bytes);
  // }
  // print_u8_array("A", A, (lambda_minus_w_grind) * lambda_bytes);

  // Switching from coloumn major to row major order
  uint8_t** v_row_maj = malloc(ellhat * sizeof(uint8_t*));
  for (unsigned i = 0; i < ellhat; i++) {
    v_row_maj[i] = calloc(lambda_minus_w_grind_bytes, 1);
  }
  column_to_row_major_V(v, v_row_maj, lambda_minus_w_grind, ellhat);

  // print_u8_array("v[0][0]", (uint8_t*)&v[0][0], 2);

  // line 14
  uint8_t* u_hi_all = malloc((n_mask * w_grind + 7) / 8);
  uint8_t* u_hi = malloc(n_mask * w_grind_bytes);
  memset(u_hi, 0, n_mask * w_grind_bytes);
  prg(rootKey, iv, ((unsigned int)2 << 30) - 1, u_hi_all, lambda, (n_mask * w_grind + 7) / 8);

  for (unsigned n_mask_idx = 0; n_mask_idx < n_mask; n_mask_idx++) {
    for (unsigned w_grind_idx = 0; w_grind_idx < w_grind; w_grind_idx++) {
      set_bit_to_pt(u_hi + n_mask_idx * w_grind_bytes, w_grind_idx, get_bit_from_pt(u_hi_all, n_mask_idx * w_grind + w_grind_idx));
    }
  }
  // print_u8_array("u_hi", u_hi, n_mask*w_grind_bytes);

  // line 15
  uint8_t* u_low = malloc(n_mask * lambda_bytes);
  memset(u_low, 0, n_mask * lambda_bytes);

  uint8_t* r_tilde = malloc(n_mask * lambda_minus_w_grind_bytes);
  memset(r_tilde, 0, n_mask * lambda_minus_w_grind_bytes);
  for (unsigned m = 0; m < n_mask; m++) {

    unsigned int ell_m = ell + m * bavc_max_node_depth(0, tau_1, k);
    unsigned int bit_idx = 0;

    for (unsigned t = 0; t < tau; t++) {
      for (unsigned b = 0; b < bavc_max_node_depth(t, tau_1, k); b++) {
        set_bit_to_pt(u_low, m * lambda + bit_idx, get_bit_from_pt(ui, (t * ellhat) + (ell + m*k) + b));
        bit_idx++;
      }
    }

    // print_u8_array("u_low", u_low + m * lambda_bytes, lambda_minus_w_grind_bytes);

    // line 17
    uint8_t* first_prod = malloc(lambda_bytes);
    memset(first_prod, 0, lambda_bytes);
    for (unsigned int col = 0; col < lambda; col++) {
      uint8_t sum = 0; 
      for (unsigned int row = 0; row < lambda; row++) {

        uint8_t bit_a = get_bit_from_pt((uint8_t*)&FAEST_128F_W_CRT[col][row / 64], row%64);
        uint8_t bit_b = get_bit_from_pt(u_low, m * lambda + row);
        sum ^= bit_a & bit_b;

      }
      set_bit_to_pt(first_prod, col, sum);
    }
    uint8_t* second_prod = malloc(lambda_bytes);
    // it has lambda_minus_w_grind degree, thus taking lambda_minus_w_grind + 1 bits
    gf2_poly_mul_ct((uint8_t*)&FAEST_128F_M_TREE, lambda_minus_w_grind + 1, (uint8_t*)&u_hi[m * w_grind_bytes], 
                        w_grind, second_prod);
    
    for (unsigned int i = 0; i < lambda_bytes; i++) {
      u_bar[m * lambda_bytes + i] = first_prod[i] ^ second_prod[i];
    }

    // print_u8_array("first_prod", first_prod, lambda_bytes);
    // print_u8_array("second_prod", second_prod, lambda_bytes);
    // print_u8_array("u_bar", u_bar + m * lambda_bytes, lambda_bytes);

    // line 18
    bit_idx = 0;
    v_idx = 0;
    for (unsigned tau_idx = 0; tau_idx < tau; tau_idx++) {
      uint8_t deg = bavc_max_node_depth(tau_idx, tau_1, k);

      uint8_t* acc = malloc((deg*2 + 7) / 8);
      memset(acc, 0, (deg*2 + 7) / 8);

      for (unsigned d = 0; d < deg; d++) {
        // uint8_t* row = malloc((deg + 7) / 8);
        for (unsigned t = 0; t < deg; t++) {
          uint8_t bit = get_bit_from_pt(v_row_maj[ell_m + d], v_idx); 
          xor_bit_to_pt(acc, t + d, bit);
          v_idx++;
        }
        v_idx -= deg;
      }
      v_idx += deg;
      // print_u8_array("sum_holder", acc, 2);

      uint8_t* reduced_out = malloc((deg*2 + 7) / 8);
      // unsigned int highest_bit_pos = highest_set_bit_pos_u64(FAEST_128F_TREE_MODULI[i]);
      gf2_poly_reduce_ct(acc, deg*2, (uint8_t*)&FAEST_128F_TREE_MODULI[tau_idx], 
                deg + 1, reduced_out);
      // print_u8_array("reduced_out", reduced_out, 1);

      // memcpy(r_tilde + m * lambda_minus_w_grind_bytes + i, reduced_out, (deg + 1 + 7) / 8);

      for (unsigned int b = 0; b < deg; b++) {
        set_bit_to_pt(r_tilde + m * lambda_minus_w_grind_bytes, bit_idx, get_bit_from_pt(reduced_out, b));
        bit_idx++;
      }

      free(acc);
      free(reduced_out);
    }
    // printf("XXXX\n");
    // print_u8_array("r_tilde", r_tilde + m * lambda_minus_w_grind_bytes, tau);
  
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
      if (get_bit_from_pt((uint8_t*)&FAEST_128F_G[e][j/64], j%64) == 1) {
        for (unsigned b = 0; b < lambda_bytes; b++) {
          L_e_zero[b] ^= a[j * lambda_bytes + b];
        }
      }
    }
    // print_u8_array("L_e_zero", L_e_zero, lambda_bytes);
    for (unsigned b = 0; b < lambda_bytes; b++) {
      L_e_one[b] ^= L_e_zero[b];
    }
    // print_u8_array("L_e_one", L_e_one, lambda_bytes);

    // line 25
    H5(iv, e, L_e_zero, h_e_zero, lambda_bytes, n_mask, lambda);
    // print_u8_array("h_e_zero + e * n_mask_bytes", h_e_zero + e * n_mask_bytes, n_mask_bytes);
    H5(iv, e, L_e_one, h_e_one, lambda_bytes, n_mask, lambda);
    // print_u8_array("h_e_one + e * n_mask_bytes",  h_e_one + e * n_mask_bytes, n_mask_bytes);

    // line 26
    for (unsigned m = 0; m < n_mask; m++) {

      uint8_t dot_product_bit = 0;
      for (unsigned i = 0; i < lambda; i++) {
        
        uint8_t bit_a = get_bit_from_pt((uint8_t*)&FAEST_128F_F[e][i / 64], i%64);
        uint8_t bit_b = get_bit_from_pt(u_bar + m * lambda_bytes, i);

        dot_product_bit ^= bit_a & bit_b;
      }

      set_bit_to_pt(v_tilde + e * n_mask_bytes, m, get_bit_from_pt(h_e_zero, m));
      uint8_t bit_a = get_bit_from_pt(v_tilde + e * n_mask_bytes, m);
      uint8_t bit_b = get_bit_from_pt(h_e_one, m);

      set_bit_to_pt(c_mult + e * n_mask_bytes, m, bit_a ^ bit_b ^ dot_product_bit);
    }
    // print_u8_array("c_mult", c_mult + e * n_mask_bytes, n_mask_bytes);
    // print_u8_array_bits("c_mult", c_mult + e * n_mask_bytes, 2);

    free(L_e_zero);
    free(L_e_one);
    free(h_e_zero);
    free(h_e_one);
  }

  // line 29-30
  for (unsigned int m = 0; m < n_mask; m++) {
    
    uint8_t* first_prod = malloc(lambda_bytes);
    memset(first_prod, 0, lambda_bytes);
    for (unsigned int row = 0; row < lambda; row++) {
      for (unsigned int col = 0; col < lambda_minus_w_grind; col++) {

        uint8_t bit_a = get_bit_from_pt((uint8_t*)&FAEST_128F_W_TREE[row][col / 64], col%64);
        uint8_t bit_b = get_bit_from_pt(r_tilde + m * lambda_minus_w_grind_bytes, col);

        xor_bit_to_pt(first_prod, row, bit_a & bit_b);

      }
    }
    // print_u8_array("r_tilde", r_tilde + m * lambda_minus_w_grind_bytes, lambda_minus_w_grind_bytes);
    // print_u8_array("first_prod", first_prod, lambda_bytes);

    uint8_t* second_prod = malloc(lambda_bytes);
    memset(second_prod, 0, lambda_bytes);

    for (unsigned int row = 0; row < lambda; row++) {
      for (unsigned int col = 0; col < n_mult; col++) {   
        
        uint8_t bit_a = get_bit_from_pt((uint8_t*)&FAEST_128F_W_GATE[row][col / 64], col % 64);
        uint8_t bit_b = get_bit_from_pt(v_tilde + col * n_mask_bytes, m);

        xor_bit_to_pt(second_prod, row, bit_a & bit_b);

      }
    }

    for (unsigned int i = 0; i < lambda_bytes; i++) {
      v_bar[m * lambda_bytes + i] = first_prod[i] ^ second_prod[i];
    }
    
    // print_u8_array("second_prod", second_prod, lambda_bytes);
    // print_u8_array("v_bar", v_bar + m * lambda_bytes, lambda_bytes);

    free(first_prod);
    free(second_prod);

  }

  // line 31
  // uint8_t* v_row_maj_combined = malloc(ell * lambda_minus_w_grind_bytes);
  // memset(v_row_maj_combined, 0, ell * lambda_minus_w_grind_bytes);

  // for (unsigned int ell_idx = 0; ell_idx < ell; ell_idx++) {
  //   for (unsigned int lwg_idx = 0; lwg_idx < lambda_minus_w_grind; lwg_idx++) {

  //     set_bit_to_pt(v_row_maj_combined + ell_idx * lambda_minus_w_grind_bytes, lwg_idx, 
  //       get_bit_from_pt(v_row_maj[ell_idx], lwg_idx)
  //     );

  //   }
  // }

  for (unsigned int i = 0; i < lambda_minus_w_grind; i++) {
    memset(v[i], 0, ellhat_bytes);
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

        acc ^= get_bit_from_pt(v_row_maj[ell_idx], row) 
          & get_bit_from_pt((uint8_t*)&FAEST_128F_W_CRT[col][row / 64], row % 64);

      }

      xor_bit_to_pt(v_row_maj_new[ell_idx], col, acc);
      // if (ell_idx == 5) {
      //   break;
      // }
      // print_u8_array("tmp",v_row_maj[ell_idx], lambda_minus_w_grind_bytes);
    }

  }

  for (unsigned int ell_idx = 0; ell_idx < 20; ell_idx++) {
    print_u8_array("v_row_maj_new[col]", v_row_maj_new[ell_idx], lambda_minus_w_grind_bytes);
  }

  column_to_row_major_V(v_row_maj_new, v, ell, lambda_minus_w_grind); 

  for (unsigned int ell_idx = 0; ell_idx < 20; ell_idx++) {
    print_u8_array("v[col]", v[ell_idx], lambda_minus_w_grind_bytes);
  }

  // print_u8_array_bits("v[col]", v[0], 1);
  // print_u8_array_bits("v[col]", v[1], 1);
  // print_u8_array_bits("v[col]", v[2], 1);
  // print_u8_array_bits("v[col]", v[3], 1);
  // print_u8_array_bits("v[col]", v[4], 1);
  // print_u8_array_bits("v[col]", v[5], 1);
  // print_u8_array_bits("v[col]", v[6], 1);
  // print_u8_array_bits("v[col]", v[7], 1);

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
                      uint8_t* q_bar, unsigned int ellhat, const faest_paramset_t* params) {
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
  // const unsigned int n_mult_bytes = (n_mult + 7) / 8;
  const unsigned int w_grind      = params->w_grind;
  // const unsigned int w_grind_bytes = (w_grind + 7) / 8;
  const unsigned int lambda_minus_w_grind = lambda - w_grind;
  const unsigned int lambda_minus_w_grind_bytes = ((lambda_minus_w_grind) + 7 ) / 8;


  uint16_t i_delta[MAX_TAU];
  if (!decode_all_chall_3(i_delta, chall_3, params)) {
    return false;
  }

  bavc_rec_t bavc_rec;
  bavc_rec.h = com;
  bavc_rec.s = malloc((L - tau) * lambda_bytes);
  assert(bavc_rec.s);

  if (!bavc_reconstruct(&bavc_rec, decom_i, i_delta, iv, params)) {
    free(bavc_rec.s);
    return false;
  }

  uint8_t* sd   = malloc((1 << k) * lambda_bytes);
  uint8_t* qtmp = malloc(MAX_DEPTH * ellhat_bytes);
  assert(sd);
  assert(qtmp);

  // Step: 1
  unsigned int q_idx = 0;
  uint8_t* sd_i      = bavc_rec.s;
  for (unsigned int i = 0; i < tau; i++) {
    // Step: 2
    const unsigned int Ni = bavc_max_node_index(i, tau_1, k);

    // Step: 6
    for (unsigned int j = 0; j < Ni; j++) {
      if (j < i_delta[i]) {
        memcpy(sd + (j ^ i_delta[i]) * lambda_bytes, sd_i + lambda_bytes * j, lambda_bytes);
      } else if (j > i_delta[i]) {
        memcpy(sd + (j ^ i_delta[i]) * lambda_bytes, sd_i + lambda_bytes * (j - 1), lambda_bytes);
      }
    }

    // Step: 7..8
    const unsigned int ki = convert_to_vole(iv, sd, true, i, ellhat_bytes, NULL, qtmp, params);

    // print_u8_array("q[q_idx] 1", q[q_idx], 2);

    // Step 11
    if (i == 0) {
      // Step 8
      memcpy(q[q_idx], qtmp, ellhat_bytes * ki);
      q_idx += ki;
    } else {
      // Step 14
      for (unsigned int d = 0; d < ki; ++d, ++q_idx) {

        memcpy(q[q_idx], qtmp + d * ellhat_bytes, ellhat_bytes);

        masked_xor_u8_array(q[q_idx], c + (i - 1) * ell_bytes, q[q_idx],
                            (i_delta[i] >> d) & 1, ell_bytes);

        // print_u8_array("q[q_idx] a c", q[q_idx], 2);
      }
    }
    sd_i += lambda_bytes * (Ni - 1);
  }

  // line 12
  uint8_t* delta_prime = malloc(lambda_minus_w_grind_bytes);
  memset(delta_prime, 0, lambda_minus_w_grind_bytes);
  unsigned int delta_prime_idx = 0;
  for (unsigned int tau_idx = 0; tau_idx < tau; tau_idx++) {
    uint8_t deg = bavc_max_node_depth(tau_idx, tau_1, k);
    for (unsigned d = 0; d < deg; d++) {
      set_bit_to_pt(delta_prime, delta_prime_idx, get_bit_from_pt((uint8_t*)&i_delta[tau_idx], d));
      delta_prime_idx++;
    }
  }

  // Switching from coloumn major to row major order
  uint8_t** q_row_maj = malloc(ellhat * sizeof(uint8_t*));
  for (unsigned i = 0; i < ellhat; i++) {
    q_row_maj[i] = calloc(lambda_minus_w_grind_bytes, 1);
  }
  column_to_row_major_V(q, q_row_maj, lambda_minus_w_grind, ellhat);

  // for (unsigned int lwg = 0; lwg < 3; lwg++) {
  //   print_u8_array("q[lwg]", q[lwg * k], 1);
  // }

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
        // uint8_t* row = malloc((deg + 7) / 8);
        for (unsigned t = 0; t < deg; t++) {
          uint8_t bit = get_bit_from_pt(q_row_maj[ell_m + d], q_idx);
          // uint8_t bit = get_bit_from_pt(q_row_maj[ell_m + d], q_idx + t); 
          // uint8_t bit = get_bit_from_pt(q[q_idx], ell_m + d); 
          // if (m == 0 && tau_idx == 0) {
          //     for (unsigned di = 0; di < deg; di++)
          //         for (unsigned ti = 0; ti < deg; ti++)
          //             printf("C col=%u pos=%u bit=%u\n",
          //                   q_idx + ti, ell_m + di,
          //                   get_bit_from_pt(q_row_maj[ell_m + di], q_idx + ti));
          // }
          xor_bit_to_pt(acc, t + d, bit);
          q_idx++;
        }
        q_idx -= deg;
      }
      q_idx += deg;
      // print_u8_array("sum_holder", acc, 2);

      uint8_t* reduced_out = malloc((deg*2 + 7) / 8);
      // unsigned int highest_bit_pos = highest_set_bit_pos_u64(FAEST_128F_TREE_MODULI[i]);
      gf2_poly_reduce_ct(acc, deg*2, (uint8_t*)&FAEST_128F_TREE_MODULI[tau_idx], 
                deg + 1, reduced_out);
      // print_u8_array("reduced_out", reduced_out, 1);

      // memcpy(r_tilde + m * lambda_minus_w_grind_bytes + i, reduced_out, (deg + 1 + 7) / 8);

      for (unsigned int b = 0; b < deg; b++) {
        set_bit_to_pt(r_tilde_prime + m * lambda_minus_w_grind_bytes, bit_idx, get_bit_from_pt(reduced_out, b));
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
      gamma_e ^= get_bit_from_pt((uint8_t*)&FAEST_128F_G[e][j/64], j%64) 
              & get_bit_from_pt(delta_prime, j); 


      if (get_bit_from_pt((uint8_t*)&FAEST_128F_G[e][j/64], j%64) == 1) {
        for (unsigned b = 0; b < lambda_bytes; b++) {
          L_e_prime[b] ^= a[j * lambda_bytes + b];
        }
      }
    }
    // print_u8_array("L_e_prime", L_e_prime, lambda_bytes);
    // print_u8_array("L_e_one", L_e_one, lambda_bytes);

    // line 25
    H5(iv, e, L_e_prime, h_e, lambda_bytes, n_mask, lambda);
    // print_u8_array_bits("h_e", h_e, 1);
    
    // line 26
    for (unsigned m = 0; m < n_mask; m++) {
        
      uint8_t sum = get_bit_from_pt(h_e, m) ^ (gamma_e & get_bit_from_pt(c_mult + e * n_mask_bytes, m));
      set_bit_to_pt(q_tilde + e * n_mask_bytes, m, sum);
    }
    // print_u8_array("q_tilde", q_tilde + e * n_mask_bytes, 1);
    // print_u8_array_bits("c_mult", c_mult + e * n_mask_bytes, 2);

    free(L_e_prime);
    free(h_e);
  }
  // print_u8_array_bits("q_tilde", q_tilde, 2);
  // print_u8_array_bits("q_tilde", q_tilde + n_mask_bytes, 2);

  // line 29-30
  for (unsigned int m = 0; m < n_mask; m++) {
    
    uint8_t* first_prod = malloc(lambda_bytes);
    memset(first_prod, 0, lambda_bytes);
    for (unsigned int row = 0; row < lambda; row++) {
      for (unsigned int col = 0; col < lambda_minus_w_grind; col++) {

        uint8_t bit_a = get_bit_from_pt((uint8_t*)&FAEST_128F_W_TREE[row][col / 64], col%64);
        uint8_t bit_b = get_bit_from_pt(r_tilde_prime + m * lambda_minus_w_grind_bytes, col);

        xor_bit_to_pt(first_prod, row, bit_a & bit_b);

      }
    }
    // print_u8_array("r_tilde_prime", r_tilde_prime + m * lambda_minus_w_grind_bytes, lambda_minus_w_grind_bytes);
    // print_u8_array("first_prod", first_prod, lambda_bytes);

    uint8_t* second_prod = malloc(lambda_bytes);
    memset(second_prod, 0, lambda_bytes);

    for (unsigned int row = 0; row < lambda; row++) {
      for (unsigned int col = 0; col < n_mult; col++) {   
        
        uint8_t bit_a = get_bit_from_pt((uint8_t*)&FAEST_128F_W_GATE[row][col / 64], col % 64);
        uint8_t bit_b = get_bit_from_pt(q_tilde + col * n_mask_bytes, m);

        xor_bit_to_pt(second_prod, row, bit_a & bit_b);

      }
    }

    for (unsigned int i = 0; i < lambda_bytes; i++) {
      q_bar[m * lambda_bytes + i] = first_prod[i] ^ second_prod[i];
    }
    
    // print_u8_array("first_prod", first_prod, lambda_bytes);
    // print_u8_array("second_prod", second_prod, lambda_bytes);
    // print_u8_array("q_bar", q_bar + m * lambda_bytes, lambda_bytes);

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

  for (unsigned int ell_idx = 0; ell_idx < ell; ell_idx++) {

    for (unsigned int col = 0; col < lambda_minus_w_grind; col++) {

      uint8_t acc = 0;
      for (unsigned int row = 0; row < lambda_minus_w_grind; row++) {

        acc ^= get_bit_from_pt(q_row_maj[ell_idx], row) 
          & get_bit_from_pt((uint8_t*)&FAEST_128F_W_CRT[col][row / 64], row % 64);

      }

      xor_bit_to_pt(q_row_maj_new[ell_idx], col, acc);
    }
    // if (ell_idx == 5) {
    //   break;
    // }
    // print_u8_array("tmp",q_row_maj[ell_idx], lambda_minus_w_grind_bytes);
  }

  for (unsigned int ell_idx = 0; ell_idx < 20; ell_idx++) {
    print_u8_array("q_row_maj_new[col]", q_row_maj_new[ell_idx], lambda_minus_w_grind_bytes);
  }

  column_to_row_major_V(q_row_maj_new, q, ell, lambda_minus_w_grind); 

  for (unsigned int ell_idx = 0; ell_idx < 20; ell_idx++) {
    print_u8_array("q[col]", q[ell_idx], lambda_minus_w_grind_bytes);
  }

  // print_u8_array_bits("v[col]", v[0], 1);
  // print_u8_array_bits("v[col]", v[1], 1);
  // print_u8_array_bits("v[col]", v[2], 1);
  // print_u8_array_bits("v[col]", v[3], 1);
  // print_u8_array_bits("v[col]", v[4], 1);
  // print_u8_array_bits("v[col]", v[5], 1);
  // print_u8_array_bits("v[col]", v[6], 1);
  // print_u8_array_bits("v[col]", v[7], 1);



  free(qtmp);
  free(sd);
  free(bavc_rec.s);
  free(delta_prime);
  for (unsigned i = 0; i < ellhat; i++) { free(q_row_maj[i]); }
  free(q_row_maj);  
  free(r_tilde_prime);
  free(a);
  free(q_tilde);
  for (unsigned i = 0; i < ell; i++) { free(q_row_maj_new[i]); }
  free(q_row_maj_new); 


  return true;
}
