#include "faest.h"

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
  memcpy(chalout, chal + lo, hi - lo);
  return 1;
}

void voleCommit(uint8_t* rootKey, uint32_t lambda, uint32_t lambdaBytes, uint32_t outlen,
                uint32_t tau, uint32_t k0, uint32_t k1, const faest_paramset_t* params,
                uint8_t* hcom, vec_com_t** vecCom, uint8_t** c, uint8_t* u, uint8_t** v) {

  tree_t** tree = malloc(sizeof(tree_t*) * tau);
  uint8_t** ui  = malloc(params->faest_param.t * sizeof(uint8_t*));
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
    c[i] = malloc(outlen);
    xorUint8Arr(u, ui[i], c[i], outlen);
  }

  H1_context_t h1_ctx;
  switch (lambda) {
  case 128:
    H1_init(&h1_ctx, 128);
    break;
  default:
    H1_init(&h1_ctx, 256);
    break;
  }
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
    pdec[i]      = malloc(depth * lambdaBytes);
    com_j[i]     = malloc(lambdaBytes * 2);
    vecComRec[i] = malloc(sizeof(vec_com_rec_t));
    vector_reconstruction(pdec[i], com_j[i], chalout, lambda, lambdaBytes, N, vecComRec[i]);
    sd[i] = malloc(N * lambdaBytes);
    for (uint32_t j = 1; j < N; j++) {
      memcpy(sd[i] + (j * lambdaBytes), vecComRec[i]->k + ((j * lambdaBytes) ^ idx), lambdaBytes);
    }
    q[i] = malloc(outlen * depth);
    ConvertToVoleVerifier(lambda, lambdaBytes, sd[i], N, depth, outlen, q[i]);
  }
  H1_context_t h1_ctx;
  switch (lambda) {
  case 128:
    H1_init(&h1_ctx, 128);
    break;
  default:
    H1_init(&h1_ctx, 256);
    break;
  }
  for (uint32_t i = 0; i < tau; i++) {
    H1_update(&h1_ctx, vecComRec[i]->com, lambdaBytes * 2);
  }
  /* Generating final commitment from all the com commitments */
  H1_final(&h1_ctx, hcom, lambdaBytes * 2);
}