/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "vole.h"
#include "tree.h"
#include "aes.h"
#include "utils.h"
#include "random_oracle.h"

#include <stdbool.h>
#include <string.h>

// TODO: Do not pass lambdaBytes everywhere, compute it in the function....

static inline uint8_t get_bit(const uint8_t* in, unsigned int index) {
  return in[index / 8] >> (7 - index % 8);
}

#if 0
static inline void set_bit(uint8_t* dst, uint8_t in, unsigned int index) {
  dst[index / 8] |= in << (7 - index % 8);
}
#endif

int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout) {
  if (i >= t0 + t1) {
    return 0;
  }

  unsigned int lo;
  unsigned int hi;
  if (i < t0) {
    lo = i * k0;
    hi = ((i + 1) * k0) - 1;
  } else {
    unsigned int t = i - t0;
    lo             = (t0 * k0) + (t * k1);
    hi             = (t0 * k0) + ((t + 1) * k1) - 1;
  }

  for (unsigned int j = lo; j != hi; ++j) {
    // set_bit(chalout, i - lo, get_bit(chal, i));
    chalout[j] = get_bit(chal, j);
  }
  return 1;
}

void voleCommit(const uint8_t* rootKey, uint32_t ellhat, const faest_paramset_t* params,
                uint8_t* hcom, vec_com_t* vecCom, uint8_t** c, uint8_t* u, uint8_t** v) {

  uint32_t lambda      = params->faest_param.lambda;
  uint32_t lambdaBytes = lambda / 8;
  uint32_t ellhatBytes = (ellhat + 7) / 8;
  uint32_t tau         = params->faest_param.tau;
  uint32_t tau0        = params->faest_param.t0;
  uint32_t k0          = params->faest_param.k0;
  uint32_t k1          = params->faest_param.k1;
  uint8_t** ui         = malloc(params->faest_param.tau * sizeof(uint8_t*));

  // Step 1
  uint8_t iv[16]         = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
  uint8_t* expanded_keys = malloc(tau * lambdaBytes);
  prg(rootKey, iv, expanded_keys, lambda, lambdaBytes * tau);

  // for Step 12
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);

  uint8_t* tmp_v     = malloc(ellhatBytes * MAX(k0, k1));
  unsigned int v_idx = 0;
  for (uint32_t i = 0; i < tau; i++) {
    uint32_t N;
    uint32_t depth;
    // Step 4
    if (i < tau0) {
      N     = 1 << k0;
      depth = k0;
    } else {
      N     = 1 << k1;
      depth = k1;
    }
    ui[i] = malloc(ellhatBytes);

    // Step 5
    vector_commitment(expanded_keys + i * lambdaBytes, params, lambda, lambdaBytes, &vecCom[i], N);
    // Step 6
    ConvertToVole(lambda, lambdaBytes, vecCom[i].sd, N, depth, ellhatBytes, ui[i], tmp_v);
    // Step 7 (and parts of 8)
    for (unsigned int j = 0; j < depth; ++j, ++v_idx) {
      memcpy(v[v_idx], tmp_v + j * ellhatBytes, ellhatBytes);
    }

    // Step 12 (part)
    H1_update(&h1_ctx, vecCom[i].com, lambdaBytes * 2);
  }
  free(tmp_v);
  free(expanded_keys);
  // Step 9
  memcpy(u, ui[0], ellhatBytes);
  free(ui[0]);
  for (uint32_t i = 1; i < tau; i++) {
    // Step 11
    c[i - 1] = malloc(ellhatBytes);
    xorUint8Arr(u, ui[i], c[i - 1], ellhatBytes);
    free(ui[i]);
  }

  // Step 12: Generating final commitment from all the com commitments
  H1_final(&h1_ctx, hcom, lambdaBytes * 2);
}

void voleReconstruct(const uint8_t* chall, uint8_t** pdec, uint8_t** com_j, uint32_t lambda,
                     uint32_t lambdaBytes, uint32_t outlen, uint32_t tau, uint32_t k0, uint32_t k1,
                     uint8_t* hcom, uint8_t** q, vec_com_rec_t* vecComRec) {
  uint32_t t0 = lambda % tau;
  uint32_t t1 = (lambda - (k0 * t0)) / k1;
  uint32_t depth;
  uint32_t N;
  uint8_t** sd = malloc(tau * sizeof(uint8_t*));

  // STep: 1
  for (uint32_t i = 0; i < tau; i++) {
    // Step: 2
    if (i < t0) {
      depth = k0;
      N     = (1 << k0);
    } else {
      depth = k1;
      N     = (1 << k1);
    }

    uint8_t* chalout = malloc(depth);
    ChalDec(chall, i, k0, t0, k1, t1, chalout);
    uint32_t idx = NumRec(depth, chalout);

    vector_reconstruction(pdec[i], com_j[i], chall, lambda, lambdaBytes, N, &vecComRec[i]);
    sd[i] = malloc(N * lambdaBytes);

    // Step: 6
    for (uint32_t j = 1; j < N; j++) {
      memcpy(sd[i] + (j * lambdaBytes), vecComRec[i].k + ((j * lambdaBytes) ^ idx), lambdaBytes);
    }
    // STep: 7..8
    q[i] = malloc(outlen * depth);
    ConvertToVole(lambda, lambdaBytes, sd[i], N, depth, outlen, NULL, q[i]);
    free(sd[i]);
  }
  free(sd);

  // Step: 9
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  for (uint32_t i = 0; i < tau; i++) {
    H1_update(&h1_ctx, vecComRec[i].com, lambdaBytes * 2);
  }
  /* Generating final commitment from all the com commitments */
  H1_final(&h1_ctx, hcom, lambdaBytes * 2);
}

static bool is_all_zeros(const uint8_t* array, size_t len) {
  for (size_t idx = 0; idx != len; ++idx) {
    if (array[idx]) {
      return false;
    }
  }

  return true;
}

void ConvertToVole(uint32_t lambda, uint32_t lambdaBytes, const uint8_t* sd,
                   uint32_t numVoleInstances, uint32_t depth, uint32_t outLenBytes, uint8_t* u,
                   uint8_t* v) {

  uint8_t iv[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
  uint8_t* r     = malloc(getBinaryTreeNodeCount(numVoleInstances) * outLenBytes);

  /* Here r_d,0 is the root and r_0,0 is the first leaf, in this reference implementation we will
  keep it consistent like k_0,0 being the root and k_d,0 being the first leaf by manipulating it
  with getNodeIndex() due to ease of understanding... */

  // Step: 2
  if (is_all_zeros(sd, lambdaBytes)) {
    memset(r + (getNodeIndex(depth, 0) * outLenBytes), 0, outLenBytes);
  } else {
    prg(sd, iv, r + (getNodeIndex(depth, 0) * outLenBytes), lambda, outLenBytes);
  }

  // Step: 3..4
  for (uint32_t i = 1; i < numVoleInstances; i++) {
    uint8_t iv_[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
    prg(sd + (lambdaBytes * i), iv_, r + (outLenBytes * (getNodeIndex(depth, 0) + i)), lambda,
        outLenBytes);
  }

  // Step: 5..9
  memset(v, 0, depth * outLenBytes);
  for (uint32_t d = 0; d < depth; d++) {
    uint32_t depthloop = (numVoleInstances + ((1 << (d + 1)) - 1)) / (1 << (d + 1));
    for (uint32_t i = 0; i < depthloop; i++) {
      for (uint8_t b = 0; b < outLenBytes; b++) {
        *(v + ((depth - 1 - d) * outLenBytes) + b) =
            *(v + ((depth - 1 - d) * outLenBytes) + b) ^
            *(r + (getNodeIndex(depth - d, 2 * i + 1) * outLenBytes) + b);

        *(r + (getNodeIndex(depth - (d + 1), i) * outLenBytes) + b) =
            *(r + (getNodeIndex(depth - d, 2 * i) * outLenBytes) + b) ^
            *(r + (getNodeIndex(depth - d, (2 * i) + 1) * outLenBytes) + b);
      }
    }
  }
  // Step: 10
  if (is_all_zeros(sd, lambdaBytes) == false && u != NULL) {
    memcpy(u, r, outLenBytes);
  }
  free(r);
}