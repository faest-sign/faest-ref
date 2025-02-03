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

#include <stdbool.h>
#include <string.h>

static const uint32_t TWEAK_OFFSET = UINT32_C(0x80000000); // 2^31

#if !defined(FAEST_TESTS)
static
#endif
    int
    ConvertToVole(const uint8_t* iv, const uint8_t* sd, bool sd0_bot, unsigned int i,
                  unsigned int outLenBytes, uint8_t* u, uint8_t* v,
                  const faest_paramset_t* params) {
  const unsigned int lambda        = params->faest_param.lambda;
  const unsigned int tau_1         = params->faest_param.tau1;
  const unsigned int k             = params->faest_param.k;
  const unsigned int num_instances = bavc_max_node_index(i, tau_1, k);
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int depth         = bavc_max_node_depth(i, tau_1, k);

  // (depth + 1) x num_instances array of outLenBytes; but we only need two rows at a time
  uint8_t* r = calloc(2 * num_instances, outLenBytes);

#define R(row, column) (r + (((row) % 2) * num_instances + (column)) * outLenBytes)
#define V(idx) (v + (idx) * outLenBytes)

  // Step: 2
  if (!sd0_bot) {
    prg(sd, iv, i, R(0, 0), lambda, outLenBytes);
  }

  // Step: 3..4
  for (unsigned int j = 1; j < num_instances; ++j) {
    prg(sd + lambda_bytes * j, iv, i, R(0, j), lambda, outLenBytes);
  }

  // Step: 5..9
  memset(v, 0, depth * outLenBytes);
  for (unsigned int j = 0; j < depth; j++) {
    unsigned int depthloop = num_instances >> (j + 1);
    for (unsigned int idx = 0; idx < depthloop; idx++) {
      xor_u8_array(V(j), R(j, 2 * idx + 1), V(j), outLenBytes);
      xor_u8_array(R(j, 2 * idx), R(j, 2 * idx + 1), R(j + 1, idx), outLenBytes);
    }
  }
  // Step: 10
  if (!sd0_bot && u != NULL) {
    memcpy(u, R(depth, 0), outLenBytes);
  }
  free(r);
  return depth;
}

void vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                 const faest_paramset_t* params, bavc_t* bavc, uint8_t* c, uint8_t* u,
                 uint8_t** v) {
  const unsigned int ellhat_bytes = (ellhat + 7) / 8;
  const unsigned int tau          = params->faest_param.tau;

  bavc_commit(rootKey, iv, params, bavc);

  uint8_t* ui = malloc(tau * ellhat_bytes);

  unsigned int v_idx = 0;
  for (unsigned int i = 0; i < tau; ++i) {
    // Step 6
    v_idx += ConvertToVole(iv, bavc->sd, false, i + TWEAK_OFFSET, ellhat_bytes,
                           ui + i * ellhat_bytes, v[v_idx], params);
  }
  // Step 9
  memcpy(u, ui, ellhat_bytes);
  for (unsigned int i = 1; i < tau; i++) {
    // Step 11
    xor_u8_array(u, ui + i * ellhat_bytes, c + (i - 1) * ellhat_bytes, ellhat_bytes);
  }
  free(ui);
}

bool vole_reconstruct(uint8_t* com, uint8_t** q, const uint8_t* iv, const uint8_t* chall_3,
                      const uint8_t* decom_i, const uint8_t* c, unsigned int ellhat,
                      const faest_paramset_t* params) {
  const unsigned int lambda       = params->faest_param.lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ellhat_bytes = (ellhat + 7) / 8;
  const unsigned int tau          = params->faest_param.tau;
  const unsigned int tau1         = params->faest_param.tau1;
  const unsigned int L            = params->faest_param.L;
  const unsigned int k            = params->faest_param.k;

  uint16_t i_delta[MAX_TAU];
  if (!decode_all_chall_3(i_delta, chall_3, params)) {
    return false;
  }

  bavc_rec_t vec_com_rec;
  vec_com_rec.h = com;
  vec_com_rec.s = malloc((L - tau) * lambda_bytes);

  if (!bavc_reconstruct(decom_i, i_delta, iv, params, &vec_com_rec)) {
    free(vec_com_rec.s);
    return false;
  }

  uint8_t* sd   = malloc((1 << k) * lambda_bytes);
  uint8_t* qtmp = malloc(MAX_DEPTH * ellhat_bytes);

  // Step: 1
  unsigned int q_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    // Step: 2
    const unsigned int Ni = bavc_max_node_index(i, tau1, k);

    // Step: 6
    for (unsigned int j = 1; j < Ni; j++) {
      memcpy(sd + j * lambda_bytes, vec_com_rec.s + lambda_bytes * (j ^ i_delta[i]), lambda_bytes);
    }

    // Step: 7..8
    const unsigned int ki =
        ConvertToVole(iv, sd, true, i + TWEAK_OFFSET, ellhat_bytes, NULL, qtmp, params);

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
  }

  free(qtmp);
  free(sd);
  free(vec_com_rec.s);
  return true;
}
