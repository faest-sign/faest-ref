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
static inline uint8_t get_bit_from_pt(const uint8_t *p, size_t i) {
    return (uint8_t)((p[i / 8] >> (i % 8)) & 1u);
}
static inline void set_bit_to_pt(uint8_t *p, size_t i, uint8_t v) {
    size_t byte = i / 8;
    uint8_t bit = (uint8_t)(1u << (i % 8));
    p[byte] = (uint8_t)((p[byte] & ~bit) | (bit & (uint8_t)(0u - v)));
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
    size_t out_bytes = (a_bits + 7) / 8;
    memset(out, 0, out_bytes);
    memcpy(out, a, out_bytes);

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

void column_to_row_major_V(uint8_t** v, uint8_t** v_row_maj, unsigned int v_row_bits_len, unsigned int v_col_bits_len) {

  for (unsigned int row = 0; row < v_row_bits_len; row++) {

        for (unsigned int col = 0; col < v_col_bits_len; col++) {

            uint8_t bit = (v[row][col / 8] >> (col % 8)) & 1u;
            v_row_maj[col][row / 8] |= (uint8_t)(bit << (row % 8));

        }
    }
}

// TODO: Modify vole_commit
void vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                 const faest_paramset_t* params, bavc_t* bavc, uint8_t* c, uint8_t* c_mult,
                 uint8_t* u, uint8_t** v, uint8_t* u_dash, uint8_t* v_dash) {
  const unsigned int lambda       = params->lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ellhat_bytes = (ellhat + 7) / 8;
  const unsigned int ell          = params->ell;
  const unsigned int ell_bytes    = (ell + 7) / 8;
  const unsigned int tau          = params->tau;
  const unsigned int tau_0        = params->tau0;
  const unsigned int tau_1        = params->tau1;
  const unsigned int k            = params->k;
  const unsigned int k_bytes      = (k+7) / 8;
  const unsigned int n_mask       = params->n_mask;
  const unsigned int n_mask_bytes = (n_mask + 7) / 8;
  const unsigned int n_mult       = params->n_mult;
  const unsigned int w_grind      = params->w_grind;
  const unsigned int w_grind_bytes = (w_grind + 7) / 8;
  const unsigned int lambda_minus_w_grind = lambda - w_grind;
  const unsigned int lambda_minus_w_grind_bytes = ((lambda - w_grind) + 7 ) / 8;

  bavc_commit(bavc, rootKey, iv, params);

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

  // line 8-9
  for (unsigned int i = 1; i < tau; i++) {
    xor_u8_array(u, ui + i * ellhat_bytes, c + (i - 1) * ell_bytes, ell_bytes);
  }

  // line 11
  uint8_t* A = malloc((lambda - w_grind) * lambda_bytes);
  for (unsigned int i = 0; i < lambda - w_grind; i++) {
    memcpy(A + i * lambda_bytes, v[i], lambda_bytes);
  }

  // Switching from coloumn major to row major order
  uint8_t** v_row_maj = malloc(ellhat * sizeof(uint8_t*));
  for (unsigned i = 0; i < ellhat; i++) {
    v_row_maj[i] = calloc(lambda_minus_w_grind_bytes, 1);
  }
  column_to_row_major_V(v, v_row_maj, lambda - w_grind, ellhat);

  // line 14
  uint8_t* u_hi_bytes = malloc(n_mask * w_grind_bytes);
  prg(rootKey, iv, ((unsigned int)2 << 30) - 1, u_hi_bytes, lambda, n_mask * w_grind_bytes);

  // line 15
  uint8_t* u_low_bytes = malloc(n_mask * lambda_bytes);
  memset(u_low_bytes, 0, n_mask * lambda_bytes);

  uint8_t* r_tilde = malloc(n_mask * lambda_minus_w_grind_bytes);
  for (unsigned m = 0; m < n_mask; m++) {

    unsigned int bit_idx = 0;
    for (unsigned t = 0; t < tau; t++) {

      // NOTE: I think we have to handle this bitwise since we take only d_i bits from u_0
      for (unsigned b = 0; b < bavc_max_node_depth(t, tau_1, k); b++) {

        set_bit_to_pt(u_low_bytes, m * lambda + bit_idx, get_bit_from_pt(u, t * ellhat + ell + m * k + b));

        bit_idx++;
      }
    }

    // line 17
    uint8_t* first_prod = malloc(lambda_bytes);
    memset(first_prod, 0, lambda_bytes);
    for (unsigned int col = 0; col < lambda; col++) {
      for (unsigned int row = 0; row < lambda; row++) {

        uint8_t bit_a = get_bit_from_pt((uint8_t*)&FAEST_128F_W_CRT[row][col / 64], col%64);
        uint8_t bit_b = get_bit_from_pt(u_low_bytes, m * lambda + row);

        set_bit_to_pt(first_prod, col, bit_a & bit_b);

      }
    }
    uint8_t* second_prod = malloc(lambda_bytes);
    gf2_poly_mul_ct((uint8_t*)&FAEST_128F_M_TREE, lambda_minus_w_grind, (uint8_t*)&u_hi_bytes[m], w_grind, second_prod);
    
    for (unsigned int i = 0; i < lambda_bytes; i++) {
      u_dash[i] = first_prod[i] ^ second_prod[i];
    }

    // line 18
    bit_idx = 0;
    for (unsigned i = 0; i < tau; i++) {

      uint8_t* sum_holder = malloc((bavc_max_node_depth(i, tau_1, k) + 7) / 8);
      for (unsigned t = 0; t < bavc_max_node_depth(i, tau_1, k); t++) {
        uint8_t bit = get_bit_from_pt(v_row_maj[i], m * k + i); 
        set_bit_to_pt(sum_holder, t, bit);
      }

      uint8_t* reduced_out = malloc((bavc_max_node_depth(i, tau_1, k) + 7) / 8);
      unsigned int highest_bit_pos = highest_set_bit_pos_u64(FAEST_128F_TREE_MODULI[i]);
      gf2_poly_reduce_ct(sum_holder, bavc_max_node_depth(i, tau_1, k), (uint8_t*)FAEST_128F_TREE_MODULI[i], 
                highest_bit_pos, reduced_out);
      for (unsigned int b = 0; b < highest_bit_pos - 1; b++) {

        set_bit_to_pt(r_tilde, m * lambda_minus_w_grind + bit_idx, get_bit_from_pt(reduced_out, b));

      }
      bit_idx++;

      free(sum_holder);
      free(reduced_out);

    }
  
    free(second_prod);
    free(first_prod);
  }

  // line 22
  uint8_t* L_e_zero = malloc(lambda_bytes);
  memset(L_e_zero, 0, lambda_bytes);
  uint8_t* L_e_one = malloc(lambda_bytes);
  // line 24
  memcpy(L_e_one, u, lambda_bytes);
  uint8_t* h_e_zero = malloc(n_mult * n_mask_bytes);
  uint8_t* h_e_one = malloc(n_mult * n_mask_bytes);
  for (unsigned e = 0; e < n_mult; e++) {

    for (unsigned j = 0; j < (lambda - params->w_grind); j++) {
      for (unsigned b = 0; b < lambda_bytes; b++) {
        L_e_zero[b] ^= A[j * lambda_bytes + b];
      }
    }
    for (unsigned b = 0; b < lambda_bytes; b++) {
      L_e_one[b] ^= L_e_zero[b];
    }

    // line 25
    H5(iv, e, L_e_zero, h_e_zero + e * n_mask_bytes, lambda_bytes, n_mask_bytes, lambda);
    H5(iv, e, L_e_one, h_e_one + e * n_mask_bytes, lambda_bytes, n_mask_bytes, lambda);

    // line 26
    for (unsigned m = 0; m < n_mask; m++) {

      uint8_t dot_product_bit = 0;
      for (unsigned i = 0; i < lambda; i++) {
        
        uint8_t bit_a = get_bit_from_pt((uint8_t*)&FAEST_128F_F[e][i / 64], i%64);
        uint8_t bit_b = get_bit_from_pt(u_dash, i);

        dot_product_bit ^= bit_a & bit_b;
      }

      uint8_t bit_a = get_bit_from_pt(h_e_zero, e * n_mask_bytes + m);
      uint8_t bit_b = get_bit_from_pt(h_e_one, e * n_mask_bytes + m);

      set_bit_to_pt(c_mult, e * n_mask_bytes + m, bit_a ^ bit_b ^ dot_product_bit);

    }
  }

  // line 29-30
  for (unsigned int m = 0; m < n_mask; m++) {
    
    uint8_t* first_prod = malloc(lambda_bytes);
    memset(first_prod, 0, lambda_bytes);
    for (unsigned int row = 0; row < lambda; row++) {
      for (unsigned int col = 0; col < lambda_minus_w_grind_bytes; col++) {

        uint8_t bit_a = get_bit_from_pt((uint8_t*)&FAEST_128F_W_TREE[row][col / 64], col%64);
        uint8_t bit_b = get_bit_from_pt(r_tilde, m * lambda_minus_w_grind + col);

        set_bit_to_pt(first_prod, row, bit_a & bit_b);

      }
    }

    uint8_t* second_prod = malloc(lambda_bytes);
    memset(second_prod, 0, lambda_bytes);

    for (unsigned int row = 0; row < lambda; row++) {

      for (unsigned int col = 0; col < n_mult; col++) {   
        

        uint8_t bit_a = get_bit_from_pt((uint8_t*)&FAEST_128F_W_GATE[row][col / 64], col % 64);
        uint8_t bit_b = get_bit_from_pt(h_e_zero, col * n_mask_bytes + m);

        set_bit_to_pt(first_prod, row, bit_a & bit_b);

      }
    }

    for (unsigned int i = 0; i < lambda_bytes; i++) {
      v_dash[i] = first_prod[i] ^ second_prod[i];
    }

    free(first_prod);
    free(second_prod);

  }

  // line 31
  uint8_t* v_row_maj_combined = malloc(ell * lambda_minus_w_grind_bytes);
  unsigned int bit_idx = 0;

  for (unsigned int ell_idx = 0; ell_idx < ell; ell_idx++) {
    for (unsigned int lwg_idx = 0; lwg_idx < lambda_minus_w_grind; lwg_idx++) {

      set_bit_to_pt(v_row_maj_combined, bit_idx, 
        get_bit_from_pt(v_row_maj[ell_idx], lwg_idx)
      );
      bit_idx++;

    }
  }

  for (unsigned int i = 0; i < lambda_minus_w_grind; i++) {
    memset(v[i], 0, ell_bytes);
  }

  for (unsigned int col = 0; col < lambda_minus_w_grind; col++) {

    for (unsigned int ell_idx = 0; ell_idx < ell; ell_idx++) {
  
      for (unsigned int row = 0; row < lambda_minus_w_grind; row++) {

        uint8_t tmp = get_bit_from_pt(v_row_maj_combined, ell_idx * lambda_minus_w_grind + row) 
          & get_bit_from_pt((uint8_t*)&FAEST_128F_W_CRT[row][col/64], col%64);

        set_bit_to_pt(v[col], ell_idx, tmp);

      }
    }
  } 

  free(v_row_maj_combined);
  free(L_e_zero);
  free(L_e_one);
  free(h_e_zero);
  free(h_e_one);
  free(r_tilde);
  free(u_low_bytes);
  free(u_hi_bytes);
  for (unsigned i = 0; i < ellhat; i++) { free(v_row_maj[i]); }
  free(v_row_maj);  
  free(A);
  free(ui);
  
}

bool vole_reconstruct(uint8_t* com, uint8_t** q, const uint8_t* iv, const uint8_t* chall_3,
                      const uint8_t* decom_i, const uint8_t* c, unsigned int ellhat,
                      const faest_paramset_t* params) {
  const unsigned int lambda       = params->lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ellhat_bytes = (ellhat + 7) / 8;
  const unsigned int tau          = params->tau;
  const unsigned int tau1         = params->tau1;
  const unsigned int L            = params->L;
  const unsigned int k            = params->k;

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
    const unsigned int Ni = bavc_max_node_index(i, tau1, k);

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

    // Step 11
    if (i == 0) {
      // Step 8
      memcpy(q[q_idx], qtmp, ellhat_bytes * ki);
      q_idx += ki;
    } else {
      // Step 14
      for (unsigned int d = 0; d < ki; ++d, ++q_idx) {
        masked_xor_u8_array(qtmp + d * ellhat_bytes, c + (i - 1) * ellhat_bytes, q[q_idx],
                            (i_delta[i] >> d) & 1, ellhat_bytes);
      }
    }
    sd_i += lambda_bytes * (Ni - 1);
  }

  // ensure 0-padding up to lambda
  for (; q_idx != lambda; ++q_idx) {
    memset(q[q_idx], 0, ellhat_bytes);
  }

  free(qtmp);
  free(sd);
  free(bavc_rec.s);
  return true;
}
