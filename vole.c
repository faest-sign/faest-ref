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

// TODO: Modify vole_commit
void vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                 const faest_paramset_t* params, bavc_t* bavc, uint8_t* c, uint8_t* c_mult,
                 uint8_t* u, uint8_t** v, uint8_t* u_dash_m, uint8_t* v_dash_m) {
  const unsigned int lambda       = params->lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ellhat_bytes = (ellhat + 7) / 8;
  const unsigned int ell_bytes    = (params->ell + 7) / 8;
  const unsigned int tau          = params->tau;
  const unsigned int tau_0        = params->tau0;
  const unsigned int tau_1        = params->tau1;
  const unsigned int k            = params->k;
  const unsigned int k_bytes      = (k+7) / 8;
  const unsigned int n_mask       = params->n_mask;
  const unsigned int n_mask_bytes = (n_mask + 7) / 8;
  const unsigned int n_mult       = params->n_mult;
  const unsigned int w_grind      = params->w_grind;

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

  // TODO: I guess we don't need this?
  // // ensure 0-padding up to lambda
  // for (; v_idx != lambda; ++v_idx) {
  //   memset(v[v_idx], 0, ellhat_bytes);
  // }

  // line 7
  memcpy(u, ui, ell_bytes);

  // line 8-9
  for (unsigned int i = 1; i < tau; i++) {
    xor_u8_array(u, ui + i * ellhat_bytes, c + (i - 1) * ell_bytes, ell_bytes);
  }
  free(ui);

  uint8_t* A = malloc((lambda - params->w_grind) * lambda_bytes);
  // line 11
  for (unsigned int i = 0; i < tau * k; i++) {
    memcpy(A, v[i], lambda_bytes);
  }

  // line 14
  uint8_t* u_hi = malloc(n_mask * (w_grind + 7) / 8);
  prg(rootKey, iv, (2 << 31) - 1, u_hi, lambda, n_mask*w_grind);

  // line 15
  uint8_t* u_low = malloc(n_mask * tau * k_bytes);
  for (unsigned m = 0; m < n_mask; m++) {

    // line 16
    for (unsigned t = 0; t < tau; t++) {
      memcpy(u_low + (m * tau * k_bytes) + (t * k_bytes), u + (t * ellhat_bytes) + (ell_bytes + m * k_bytes), 
              (bavc_max_node_index(t, tau_1, k - 1) + 7) / 8);
    }

    // line 17
    // uint8_t* u_dash_m = malloc(lambda_bytes);
    for (unsigned r = 0; r < lambda; r++) {
      for (unsigned c = 0; c < lambda - w_grind; c++) {
        // TODO:
      }
    }

    // line 18
    uint8_t* r_tilde = malloc(tau * n_mask * lambda_bytes);
    for (unsigned i = 0; i < tau; i++) {
      // TODO:
    }

    // line 22
    uint8_t* L_e_zero = malloc(lambda_bytes);
    memset(L_e_zero, 0, lambda_bytes);
    uint8_t* L_e_one = malloc(lambda_bytes);
    // line 24
    memcpy(L_e_one, u, lambda_bytes);
    uint8_t* h_e_zero = malloc(n_mult * n_mask_bytes);
    uint8_t* h_e_one = malloc(n_mult * n_mask_bytes);
    // uint8_t* c_mult = malloc(n_mult * n_mask_bytes);
    for (unsigned e = 0; e < n_mult; e++) {

      // line 24
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
      for (unsigned m = 0; m < n_mask_bytes; m++) {
        c_mult[e * n_mask_bytes + m] = h_e_zero[e * n_mask_bytes + m] ^ h_e_one[e * n_mask_bytes + m]; // TODO: ^ <F_e, u_dash_,>
      }
    }

    // line 29-30
    // uint8_t* v_dash_m = malloc(lambda_bytes);
    for (unsigned int m = 0; m < n_mask_bytes; m++) {
      // TODO:
    }

    // line 31
    // TODO:

    // free(u_dash_m);
    free(r_tilde);
    free(L_e_zero);
    free(L_e_one);
    free(h_e_zero);
    free(h_e_one);
    // free(c_mult);
    // free(v_dash_m);

  }  

  free(A);
  free(u_hi);
  free(u_low);
  
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
