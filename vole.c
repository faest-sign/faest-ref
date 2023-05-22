#include "vole.h"

// TODO: Do not pass lambdaBytes everywhere, compute it in the function....

int ChalDec(const uint8_t* chal, uint32_t i, uint32_t k0, uint32_t t0, uint32_t k1, uint32_t t1,
            uint8_t* chalout) {
  uint32_t lo;
  uint32_t hi;
  uint32_t t;
  if (i >= t0 + t1) {
    return 0;
  }
  if (i < t0) {
    lo = i * k0;
    hi = ((i + 1) * k0) - 1;
  } else {
    t  = i - t0;
    lo = (t0 * k0) + (t * k1);
    hi = (t0 * k0) + ((t + 1) * k1) - 1;
  }
  chalout = malloc(hi - lo);
  // TODO: this should be always 8...
  memcpy(chalout, chal + lo, hi - lo);
  return 1;
}

// TODO: outlen (l) is in bits, change it everywhere
void voleCommit(const uint8_t* rootKey, uint32_t lambda, uint32_t lambdaBytes, uint32_t outlen,
                const faest_paramset_t* params, uint8_t* hcom, vec_com_t** vecCom, uint8_t** c,
                uint8_t* u, uint8_t** v) {

  uint32_t tau = params->faest_param.t;
  uint32_t k0  = params->faest_param.k0;
  uint32_t k1  = params->faest_param.k1;
  uint8_t** ui = malloc(params->faest_param.t * sizeof(uint8_t*));
  uint32_t N;
  uint32_t depth;
  uint8_t** keys = malloc(tau * sizeof(uint8_t*));

  uint8_t iv[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
  uint8_t* out   = malloc(tau * lambdaBytes);
  prg(rootKey, iv, out, lambda, lambdaBytes * tau);
  for (uint32_t i = 0; i < tau; i++) {
    keys[i] = malloc(lambdaBytes);
    memcpy(keys[i], out + (i * lambdaBytes), lambdaBytes);
  }

  for (uint32_t i = 0; i < tau; i++) {
    if (i < (lambda % tau)) {
      N     = 1 << k0;
      depth = k0;
    } else {
      N     = 1 << k1;
      depth = k1;
    }
    ui[i]     = malloc(outlen);
    v[i]      = malloc(depth * outlen);
    vecCom[i] = malloc(sizeof(vec_com_t));
    vector_commitment(keys[i], params, lambda, lambdaBytes, vecCom[i], N);
    ConvertToVoleProver(lambda, lambdaBytes, vecCom[i]->sd, N, depth, outlen, ui[i], v[i]);
  }
  memcpy(u, ui, outlen);
  for (uint32_t i = 1; i < tau; i++) {
    c[i - 1] = malloc(outlen);
    xorUint8Arr(u, ui[i], c[i - 1], outlen);
  }

  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  for (uint32_t i = 0; i < tau; i++) {
    H1_update(&h1_ctx, vecCom[i]->com, lambdaBytes * 2);
  }
  /* Generating final commitment from all the com commitments */
  H1_final(&h1_ctx, hcom, lambdaBytes * 2);
}

void voleVerify(const uint8_t* chal, uint8_t** pdec, uint8_t** com_j, uint32_t lambda,
                uint32_t lambdaBytes, uint32_t outlen, uint32_t tau, uint32_t k0, uint32_t k1,
                uint8_t* hcom, uint8_t** q, vec_com_rec_t** vecComRec) {
  uint32_t t0 = lambda % tau;
  uint32_t t1 = (lambda - (k0 * t0)) / k1;
  uint32_t depth;
  uint32_t N;
  uint8_t** sd = malloc(tau * sizeof(uint8_t*));
  for (uint32_t i = 0; i < tau; i++) {
    if (i < t0) {
      depth = k0;
      N     = (1 << k0);
    } else {
      depth = k1;
      N     = (1 << k1);
    }

    uint8_t* chalout = malloc(depth);
    ChalDec(chal, i, k0, t0, k1, t1, chalout);
    uint32_t idx = NumRec(depth, chalout);
    vecComRec[i] = malloc(sizeof(vec_com_rec_t));

    vector_reconstruction(pdec[i], com_j[i], chal, lambda, lambdaBytes, N, vecComRec[i]);
    sd[i] = malloc(N * lambdaBytes);
    for (uint32_t j = 1; j < N; j++) {
      memcpy(sd[i] + (j * lambdaBytes), vecComRec[i]->k + ((j * lambdaBytes) ^ idx), lambdaBytes);
    }
    q[i] = malloc(outlen * depth);
    ConvertToVoleVerifier(lambda, lambdaBytes, sd[i], N, depth, outlen, q[i]);
  }
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  for (uint32_t i = 0; i < tau; i++) {
    H1_update(&h1_ctx, vecComRec[i]->com, lambdaBytes * 2);
  }
  /* Generating final commitment from all the com commitments */
  H1_final(&h1_ctx, hcom, lambdaBytes * 2);
}

void ConvertToVoleProver(uint32_t lambda, uint32_t lambdaBytes, const uint8_t* sd,
                         uint32_t numVoleInstances, uint32_t depth, uint32_t outLenBytes,
                         uint8_t* u, uint8_t* v) {

  uint8_t iv[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
  uint8_t* r     = malloc(getBinaryTreeNodeCount(numVoleInstances) * outLenBytes);

  /* Here r_d,0 is the root and r_0,0 is the first leaf, in this reference implementation we will
  keep it consistent like k_0,0 being the root and k_d,0 being the first leaf by manipulating it
  with getNodeIndex() due to ease of understanding... */

  uint8_t* allZeros = malloc(lambdaBytes);
  memset(allZeros, 0, lambdaBytes);
  if (memcmp(sd, allZeros, lambdaBytes) == 0) {
    memset(r + (getNodeIndex(depth, 0) * outLenBytes), 0, outLenBytes);
  } else {
    uint8_t* out = malloc(outLenBytes);
    prg(sd, iv, out, lambda, outLenBytes);
    memcpy(r + (getNodeIndex(depth, 0) * outLenBytes), out, outLenBytes);
  }

  for (uint32_t i = 1; i < numVoleInstances; i++) {
    uint8_t iv_[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
    uint8_t* out    = malloc(outLenBytes);
    prg(sd + (lambdaBytes * i), iv_, out, lambda, outLenBytes);
    memcpy(r + (outLenBytes * (getNodeIndex(depth, 0) + i)), out, outLenBytes);
  }

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
  memcpy(u, r, outLenBytes);
}

void ConvertToVoleVerifier(uint32_t lambda, uint32_t lambdaBytes, const uint8_t* sd,
                           uint32_t numVoleInstances, uint32_t depth, uint32_t outLenBytes,
                           uint8_t* v) {
  uint8_t iv[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
  uint8_t* r     = malloc(getBinaryTreeNodeCount(numVoleInstances) * outLenBytes);

  /* Here r_d,0 is the root and r_0,0 is the first leaf, in this reference implementation we will
  keep it consistent like k_0,0 being the root and k_d,0 being the first leaf by manipulating it
  with getNodeIndex() due to ease of understanding... */

  // TODO: stupid way to check,, change it !!
  uint8_t* allZeros = malloc(lambdaBytes);
  memset(allZeros, 0, lambdaBytes);
  if (memcmp(sd, allZeros, lambdaBytes) == 0) {
    memset(r + (getNodeIndex(depth, 0) * outLenBytes), 0, outLenBytes);
  } else {
    uint8_t* out = malloc(outLenBytes);
    prg(sd, iv, out, lambda, outLenBytes);
    memcpy(r + (getNodeIndex(depth, 0) * outLenBytes), out, outLenBytes);
  }

  for (uint32_t i = 1; i < numVoleInstances; i++) {
    uint8_t iv_[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
    uint8_t* out    = malloc(outLenBytes);
    prg(sd + (lambdaBytes * i), iv_, out, lambda, outLenBytes);
    memcpy(r + (outLenBytes * (getNodeIndex(depth, 0) + i)), out, outLenBytes);
  }

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
}