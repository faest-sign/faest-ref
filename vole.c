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

// TODO: Do not pass lambdaBytes everywhere, compute it in the function....

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

void voleCommit(const uint8_t* rootKey, const uint8_t* iv, uint32_t ellhat,
                const faest_paramset_t* params, uint8_t* hcom, vec_com_t* vecCom, uint8_t** c,
                uint8_t* u, uint8_t** v) {
  uint32_t lambda      = params->faest_param.lambda;
  uint32_t lambdaBytes = lambda / 8;
  uint32_t ellhatBytes = (ellhat + 7) / 8;
  uint32_t tau         = params->faest_param.tau;
  uint32_t tau0        = params->faest_param.t0;
  uint32_t k0          = params->faest_param.k0;
  uint32_t k1          = params->faest_param.k1;

  uint8_t* ui = malloc(tau * ellhatBytes);

  // Step 1
  uint8_t* expanded_keys = malloc(tau * lambdaBytes);
  prg(rootKey, iv, expanded_keys, lambda, lambdaBytes * tau);

  // for Step 12
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);

  uint8_t* tmp_v     = malloc(ellhatBytes * MAX(k0, k1));
  unsigned int v_idx = 0;
  for (uint32_t i = 0; i < tau; i++) {
    // Step 4
    uint32_t depth = i < tau0 ? k0 : k1;
    uint32_t N     = 1 << depth;

    // Step 5
    vector_commitment(expanded_keys + i * lambdaBytes, iv, params, lambda, lambdaBytes, &vecCom[i],
                      N);
    // Step 6
    ConvertToVole(iv, vecCom[i].sd, false, lambda, lambdaBytes, N, depth, ellhatBytes,
                  ui + i * ellhatBytes, tmp_v);
    // Step 7 (and parts of 8)
    for (unsigned int j = 0; j < depth; ++j, ++v_idx) {
      memcpy(v[v_idx], tmp_v + j * ellhatBytes, ellhatBytes);
    }
    // Step 12 (part)
    H1_update(&h1_ctx, vecCom[i].h, lambdaBytes * 2);
  }
  free(tmp_v);
  free(expanded_keys);
  // Step 9
  memcpy(u, ui, ellhatBytes);
  for (uint32_t i = 1; i < tau; i++) {
    // Step 11
    xorUint8Arr(u, ui + i * ellhatBytes, c[i - 1], ellhatBytes);
  }
  free(ui);

  // Step 12: Generating final commitment from all the com commitments
  H1_final(&h1_ctx, hcom, lambdaBytes * 2);
}

void voleReconstruct(const uint8_t* iv, const uint8_t* chall, uint8_t** pdec, uint8_t** com_j,
                     uint8_t* hcom, uint8_t** q, uint32_t ellhat, const faest_paramset_t* params) {
  uint32_t lambda      = params->faest_param.lambda;
  uint32_t lambdaBytes = lambda / 8;
  uint32_t ellhatBytes = (ellhat + 7) / 8;
  uint32_t tau         = params->faest_param.tau;
  uint32_t tau0        = params->faest_param.t0;
  uint32_t tau1        = params->faest_param.t1;
  uint32_t k0          = params->faest_param.k0;
  uint32_t k1          = params->faest_param.k1;

  uint8_t* sd = malloc((1 << MAX(k0, k1)) * lambdaBytes);
  memset(sd, 0, lambdaBytes);

  uint8_t* tmp_q     = malloc(ellhatBytes * MAX(k0, k1));
  unsigned int q_idx = 0;

  // Step 9
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);

  vec_com_rec_t vecComRec;
  vecComRec.h   = malloc(lambdaBytes * 2);
  vecComRec.k   = calloc(getBinaryTreeNodeCount(1 << MAX(k0, k1)), lambdaBytes);
  vecComRec.com = malloc((1 << MAX(k0, k1)) * lambdaBytes * 2);
  vecComRec.s   = malloc((1 << MAX(k0, k1)) * lambdaBytes);

  // Step: 1
  for (uint32_t i = 0; i < tau; i++) {
    // Step: 2
    uint32_t depth = i < tau0 ? k0 : k1;
    uint32_t N     = 1 << depth;

    // Step 3
    uint8_t chalout[MAX_DEPTH];
    ChalDec(chall, i, k0, tau0, k1, tau1, chalout);
    // Step 4
    uint32_t idx = NumRec(depth, chalout);
    assert(idx < N);

    // Step 5
    vector_reconstruction(iv, pdec[i], com_j[i], chalout, lambda, lambdaBytes, N, depth,
                          &vecComRec);

    // Step: 6
    for (uint32_t j = 1; j < N; j++) {
      memcpy(sd + j * lambdaBytes, vecComRec.s + (lambdaBytes * (j ^ idx)), lambdaBytes);
    }

    // Step: 7..8
    ConvertToVole(iv, sd, true, lambda, lambdaBytes, N, depth, ellhatBytes, NULL, tmp_q);
    for (unsigned int j = 0; j < depth; ++j, ++q_idx) {
      // TODO: get rid of this memcpy
      memcpy(q[q_idx], tmp_q + j * ellhatBytes, ellhatBytes);
    }

    // Step 9
    H1_update(&h1_ctx, vecComRec.h, lambdaBytes * 2);
  }
  vec_com_rec_clear(&vecComRec);
  free(tmp_q);
  free(sd);

  // Step: 9
  H1_final(&h1_ctx, hcom, lambdaBytes * 2);
}

void ConvertToVole(const uint8_t* iv, const uint8_t* sd, bool sd0_bot, uint32_t lambda,
                   uint32_t lambdaBytes, uint32_t numVoleInstances, uint32_t depth,
                   uint32_t outLenBytes, uint8_t* u, uint8_t* v) {
  // (depth + 1) x numVoleInstances array of outLenBytes; but we only need to rows at a time
  uint8_t* r = calloc(2 * numVoleInstances, outLenBytes);

#define R(row, column) (r + (((row) % 2) * numVoleInstances + (column)) * outLenBytes)
#define V(idx) (v + (idx)*outLenBytes)

  // Step: 2
  if (!sd0_bot) {
    prg(sd, iv, R(0, 0), lambda, outLenBytes);
  }

  // Step: 3..4
  for (uint32_t i = 1; i < numVoleInstances; i++) {
    prg(sd + (lambdaBytes * i), iv, R(0, i), lambda, outLenBytes);
  }

  // Step: 5..9
  memset(v, 0, depth * outLenBytes);
  for (uint32_t j = 0; j < depth; j++) {
    uint32_t depthloop = (numVoleInstances >> (j + 1));
    for (uint32_t i = 0; i < depthloop; i++) {
      xorUint8Arr(V(j), R(j, 2 * i + 1), V(j), outLenBytes);
      xorUint8Arr(R(j, 2 * i), R(j, 2 * i + 1), R(j + 1, i), outLenBytes);
    }
  }
  // Step: 10
  if (!sd0_bot && u != NULL) {
    memcpy(u, R(depth, 0), outLenBytes);
  }
  free(r);
}
