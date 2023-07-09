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

#if !defined(FAEST_TESTS)
static
#endif
    void
    ConvertToVole(const uint8_t* iv, const uint8_t* sd, bool sd0_bot, unsigned int lambda,
                  unsigned int depth, unsigned int outLenBytes, uint8_t* u, uint8_t* v) {
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;

  // (depth + 1) x num_instances array of outLenBytes; but we only need to rows at a time
  uint8_t* r = calloc(2 * num_instances, outLenBytes);

#define R(row, column) (r + (((row) % 2) * num_instances + (column)) * outLenBytes)
#define V(idx) (v + (idx)*outLenBytes)

  // Step: 2
  if (!sd0_bot) {
    prg(sd, iv, R(0, 0), lambda, outLenBytes);
  }

  // Step: 3..4
  for (unsigned int i = 1; i < num_instances; i++) {
    prg(sd + (lambda_bytes * i), iv, R(0, i), lambda, outLenBytes);
  }

  // Step: 5..9
  memset(v, 0, depth * outLenBytes);
  for (unsigned int j = 0; j < depth; j++) {
    unsigned int depthloop = num_instances >> (j + 1);
    for (unsigned int i = 0; i < depthloop; i++) {
      xor_u8_array(V(j), R(j, 2 * i + 1), V(j), outLenBytes);
      xor_u8_array(R(j, 2 * i), R(j, 2 * i + 1), R(j + 1, i), outLenBytes);
    }
  }
  // Step: 10
  if (!sd0_bot && u != NULL) {
    memcpy(u, R(depth, 0), outLenBytes);
  }
  free(r);
}

int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout) {
  if (i >= t0 + t1) {
    return 0;
  }

  unsigned int lo;
  unsigned int hi;
  if (i < t0) {
    lo = i * k0;
    hi = ((i + 1) * k0);
  } else {
    unsigned int t = i - t0;
    lo             = (t0 * k0) + (t * k1);
    hi             = (t0 * k0) + ((t + 1) * k1);
  }

  assert(hi - lo == k0 || hi - lo == k1);
  for (unsigned int j = lo; j < hi; ++j) {
    // set_bit(chalout, i - lo, get_bit(chal, i));
    chalout[j - lo] = ptr_get_bit(chal, j);
  }
  return 1;
}

void vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                 const faest_paramset_t* params, uint8_t* hcom, vec_com_t* vecCom, uint8_t* c,
                 uint8_t* u, uint8_t** v) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;

  uint8_t* ui = malloc(tau * ellhat_bytes);

  // Step 1
  uint8_t* expanded_keys = malloc(tau * lambda_bytes);
  prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);

  // for Step 12
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);

  unsigned int v_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    // Step 4
    unsigned int depth = i < tau0 ? k0 : k1;

    // Step 5
    vector_commitment(expanded_keys + i * lambda_bytes, iv, params, lambda, &vecCom[i], depth);
    // Step 6
    ConvertToVole(iv, vecCom[i].sd, false, lambda, depth, ellhat_bytes, ui + i * ellhat_bytes,
                  v[v_idx]);
    // Step 7 (and parts of 8)
    v_idx += depth;
    // Step 12 (part)
    H1_update(&h1_ctx, vecCom[i].h, lambda_bytes * 2);
  }
  free(expanded_keys);
  // Step 9
  memcpy(u, ui, ellhat_bytes);
  for (unsigned int i = 1; i < tau; i++) {
    // Step 11
    xor_u8_array(u, ui + i * ellhat_bytes, c + (i - 1) * ellhat_bytes, ellhat_bytes);
  }
  free(ui);

  // Step 12: Generating final commitment from all the com commitments
  H1_final(&h1_ctx, hcom, lambda_bytes * 2);
}

void vole_reconstruct(const uint8_t* iv, const uint8_t* chall, const uint8_t* const* pdec,
                      const uint8_t* const* com_j, uint8_t* hcom, uint8_t** q, unsigned int ellhat,
                      const faest_paramset_t* params) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int tau1         = params->faest_param.t1;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;

  uint8_t* sd = malloc((1 << MAX(k0, k1)) * lambda_bytes);
  memset(sd, 0, lambda_bytes);

  // Step 9
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);

  vec_com_rec_t vecComRec;
  vecComRec.h   = malloc(lambda_bytes * 2);
  vecComRec.k   = calloc(getBinaryTreeNodeCount(MAX(k0, k1)), lambda_bytes);
  vecComRec.com = malloc((1 << MAX(k0, k1)) * lambda_bytes * 2);
  vecComRec.s   = malloc((1 << MAX(k0, k1)) * lambda_bytes);

  // Step: 1
  unsigned int q_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    // Step: 2
    unsigned int depth = i < tau0 ? k0 : k1;
    unsigned int N     = 1 << depth;

    // Step 3
    uint8_t chalout[MAX_DEPTH];
    ChalDec(chall, i, k0, tau0, k1, tau1, chalout);
    // Step 4
    unsigned int idx = NumRec(depth, chalout);

    // Step 5
    vector_reconstruction(iv, pdec[i], com_j[i], chalout, lambda, depth, &vecComRec);

    // Step: 6
    for (unsigned int j = 1; j < N; j++) {
      memcpy(sd + j * lambda_bytes, vecComRec.s + (lambda_bytes * (j ^ idx)), lambda_bytes);
    }

    // Step: 7..8
    ConvertToVole(iv, sd, true, lambda, depth, ellhat_bytes, NULL, q[q_idx]);
    q_idx += depth;

    // Step 9
    H1_update(&h1_ctx, vecComRec.h, lambda_bytes * 2);
  }
  vec_com_rec_clear(&vecComRec);
  free(sd);

  // Step: 9
  H1_final(&h1_ctx, hcom, lambda_bytes * 2);
}
