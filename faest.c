/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest.h"
#include "faest_aes.h"
#include "randomness.h"
#include "random_oracle.h"
#include "utils.h"

// TODO: TEST EVERYTHING HERE !!!

// TODO: Do not pass lambdaBytes everywhere, compute it in the function....

// TODO: change q to Q where applicable

void sign(const uint8_t* msg, size_t msglen, const uint8_t* sk, const uint8_t* pk,
          const faest_paramset_t* params, signature_t* signature) {
  const uint32_t l           = params->faest_param.l;
  const uint32_t ell_bytes   = (l + 7) / 8;
  const uint32_t lambda      = params->faest_param.lambda;
  const uint32_t lambdaBytes = lambda / 8;
  const uint32_t tau         = params->faest_param.tau;
  const uint32_t tau0        = params->faest_param.t0;
  const uint32_t ell_hat =
      params->faest_param.l + params->faest_param.lambda * 2 + params->faest_param.b;
  const uint32_t ell_hat_bytes = (ell_hat + 7) / 8;
  const size_t utilde_bytes    = (params->faest_param.lambda + params->faest_param.b + 7) / 8;

  uint8_t* rho = malloc(lambdaBytes);
  rand_bytes(rho, lambdaBytes);

  // Step: 1
  uint8_t* mu = malloc(lambdaBytes * 2);
  {
    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda);
    H1_update(&h1_ctx, pk, params->faest_param.pkSize);
    H1_update(&h1_ctx, msg, msglen);
    H1_final(&h1_ctx, mu, lambdaBytes * 2);
  }

  uint8_t* rootkey = malloc(lambdaBytes);
  {
    H3_context_t h3_ctx;
    H3_init(&h3_ctx, lambda);
    H3_update(&h3_ctx, sk, params->faest_param.skSize);
    H3_update(&h3_ctx, mu, lambdaBytes * 2);
    H3_update(&h3_ctx, rho, lambdaBytes);
    H3_final(&h3_ctx, rootkey, lambdaBytes);
  }

  // Step: 3..4
  vec_com_t* vecCom = calloc(tau, sizeof(vec_com_t));
  uint8_t* u_       = malloc(l);
  uint8_t** v       = malloc(tau * sizeof(uint8_t*));
  voleCommit(rootkey, ell_hat, params, signature->hcom, vecCom, signature->c, u_, v);
  free(rootkey);

  // Step: 5
  uint8_t* chal_1 = malloc((5 * lambdaBytes) + 8);
  {
    H2_context_t h2_ctx;
    H2_init(&h2_ctx, lambda);
    H2_update(&h2_ctx, mu, lambdaBytes * 2);
    H2_update(&h2_ctx, signature->hcom, lambdaBytes * 2);
    for (unsigned int i = 0; i < (tau - 1); ++i) {
      H2_update(&h2_ctx, signature->c[i], ell_hat_bytes);
    }
    H2_final(&h2_ctx, chal_1, (5 * lambdaBytes) + 8);
  }

  // Step: 7
  vole_hash(signature->u_tilde, chal_1, u_, l, lambda);

  // Step: 8 and 9
  {
    H1_context_t h1_ctx_1;
    H1_init(&h1_ctx_1, lambda);

    uint8_t* V_tilde = malloc(lambdaBytes + UNIVERSAL_HASH_B);
    for (unsigned int i = 0; i != lambda; ++i) {
      vole_hash(V_tilde, chal_1, v[i], ell_hat, lambda);
      H1_update(&h1_ctx_1, V_tilde, lambdaBytes + UNIVERSAL_HASH_B);
    }
    free(V_tilde);

    // Step: 9
    H1_final(&h1_ctx_1, signature->h_v, lambdaBytes * 2);
  }
  free(chal_1);

  // Step: 10..11
  const uint8_t* in  = pk;
  const uint8_t* out = pk + params->faest_param.pkSize / 2;
  // Step: 12..13
  // TODO
  uint8_t* w = aes_extend_witness(lambda, l, sk, in);
  xorUint8Arr(w, u_, signature->d, ell_bytes);

  // Step: 14
  uint8_t* chal_2 = malloc(3 * lambdaBytes + 8);
  {
    H2_context_t h2_ctx_1;
    H2_init(&h2_ctx_1, lambda);
    H2_update(&h2_ctx_1, chal_1, lambdaBytes);
    H2_update(&h2_ctx_1, signature->u_tilde, utilde_bytes);
    H2_update(&h2_ctx_1, signature->h_v, 2 * lambdaBytes);
    H2_update(&h2_ctx_1, signature->d, ell_bytes);
    H2_final(&h2_ctx_1, chal_2, (3 * lambdaBytes) + 8);
  }

  // Step: 15..17
  // TODO
  // u = u[0..l+lambda)
  // v = v[0...l+lambda]

  // Step: 18
  // TODO
  unsigned int R, beta, Lke, Lenc, C, Nwd, Ske, Senc;
  aes_prove(w, u_, v, in, out, chal_2, lambda, R, tau, l, beta, Lke, Lenc, C, Nwd, Ske, Senc,
            signature->a_tilde, signature->b_tilde);

  // Step: 19
  uint8_t* chal_3 = malloc(lambdaBytes);
  {
    H2_context_t h2_ctx_2;
    H2_init(&h2_ctx_2, lambda);
    H2_update(&h2_ctx_2, chal_2, 3 * lambdaBytes + 8);
    H2_update(&h2_ctx_2, signature->a_tilde, lambdaBytes);
    H2_update(&h2_ctx_2, signature->b_tilde, lambdaBytes);
    H2_final(&h2_ctx_2, chal_3, lambdaBytes);
  }
  free(chal_2);

  // Step: 20..23
  uint8_t* s_ = malloc(MAX(params->faest_param.k0, params->faest_param.k1));
  for (uint32_t i = 0; i < tau; i++) {
    ChalDec(chal_3, i, params->faest_param.k0, params->faest_param.t0, params->faest_param.k1,
            params->faest_param.t1, s_);
    unsigned int num_vole_instances =
        i < tau0 ? (1 << params->faest_param.k0) : (1 << params->faest_param.k1);
    vector_open(vecCom[i].k, vecCom[i].com, s_, signature->pdec[i], signature->com_j[i],
                num_vole_instances, lambdaBytes);
    vec_com_clear(&vecCom[i]);
  }
  free(s_);
  free(chal_3);
  free(vecCom);
}

// TODO: l is in bits, change it everywhere
int verify(const uint8_t* msg, size_t msglen, const uint8_t* pk, const faest_paramset_t* params,
           const signature_t* signature) {
  const uint32_t l         = params->faest_param.l;
  const uint32_t ell_bytes = (params->faest_param.l + 7) / 8;

  const uint32_t lambda      = params->faest_param.lambda;
  const uint32_t lambdaBytes = lambda / 8;
  const uint32_t tau         = params->faest_param.tau;
  const uint32_t tau0        = params->faest_param.t0;
  const uint32_t ell_hat =
      params->faest_param.l + params->faest_param.lambda * 2 + params->faest_param.b;
  const uint32_t ell_hat_bytes = (ell_hat + 7) / 8;
  const size_t utilde_bytes    = (params->faest_param.lambda + params->faest_param.b + 7) / 8;

  // Step: 2
  const uint8_t* in  = pk;
  const uint8_t* out = pk + params->faest_param.pkSize / 2;

  // Step: 3
  uint8_t* mu = malloc(lambdaBytes * 2);
  {
    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda);
    H1_update(&h1_ctx, pk, params->faest_param.pkSize);
    H1_update(&h1_ctx, msg, msglen);
    H1_final(&h1_ctx, mu, lambdaBytes * 2);
  }

  // Step: 5
  uint8_t* chal_1 = malloc((5 * lambdaBytes) + 8);
  {
    H2_context_t h2_ctx;
    H2_init(&h2_ctx, lambda);
    H2_update(&h2_ctx, mu, lambdaBytes * 2);
    H2_update(&h2_ctx, signature->hcom, lambdaBytes * 2);
    for (unsigned int i = 0; i < (tau - 1); ++i) {
      H2_update(&h2_ctx, signature->c[i], ell_hat_bytes);
    }
    H2_final(&h2_ctx, chal_1, (5 * lambdaBytes) + 8);
  }

  // Step 6
  uint8_t* chal_2 = malloc(3 * lambdaBytes + 8);
  {
    H2_context_t h2_ctx_1;
    H2_init(&h2_ctx_1, lambda);
    H2_update(&h2_ctx_1, chal_1, lambdaBytes);
    H2_update(&h2_ctx_1, signature->u_tilde, utilde_bytes);
    H2_update(&h2_ctx_1, signature->h_v, 2 * lambdaBytes);
    H2_update(&h2_ctx_1, signature->d, ell_bytes);
    H2_final(&h2_ctx_1, chal_2, (3 * lambdaBytes) + 8);
  }

  // Step: 7
  uint8_t* chal_3 = malloc(lambdaBytes);
  {
    H2_context_t h2_ctx_2;
    H2_init(&h2_ctx_2, lambda);
    H2_update(&h2_ctx_2, chal_2, 3 * lambdaBytes + 8);
    H2_update(&h2_ctx_2, signature->a_tilde, lambdaBytes);
    H2_update(&h2_ctx_2, signature->b_tilde, lambdaBytes);
    H2_final(&h2_ctx_2, chal_3, lambdaBytes);
  }
  free(chal_2);

  // Step: 9..11
  uint8_t** q              = malloc(tau * sizeof(uint8_t*));
  vec_com_rec_t* vecComRec = malloc(sizeof(vec_com_rec_t));
  uint8_t* hcomverify      = malloc(lambdaBytes * 2);
  voleVerify(chal_3, signature->pdec, signature->com_j, lambda, lambdaBytes, l, tau,
             params->faest_param.k0, params->faest_param.k1, hcomverify, q, vecComRec);
  if (memcmp(hcomverify, signature->hcom, lambdaBytes * 2) != 0) {
    return 0;
  }

  // Step: 12
  // TODO: no clue what this will take... FIX Q !!
  uint8_t* r0;
  uint8_t** q_tilde = malloc(tau * sizeof(uint8_t*));
  uint32_t depth    = 0;
  size_t ell;
  q_tilde[0] = malloc(l * depth);
  memcpy(q_tilde[0], q[0], sizeof(q[0]));
  for (uint32_t i = 1; i < tau; i++) {
    if (i < (lambda % tau)) { // Computing the num of vole and the depth here
      depth = params->faest_param.k0;
    } else {
      depth = params->faest_param.k1;
    }
    q_tilde[i] = malloc(l * depth);
    for (uint32_t d = 0; d < depth; d++) {
      for (uint32_t k = 0; k < l; k++) {
        *(q_tilde[i] + (d * l) + k) = *(q[i] + (d * l) + k) ^ *(signature->c[i] + (d * l) + k);
      }
    }
  }
  vole_hash(vecComRec->h, chal_1, q_tilde, ell, lambda);

  // Step: 13..16
  uint8_t b         = 0;
  uint8_t** d_tilde = malloc(tau * sizeof(uint8_t*));
  for (uint32_t i = 0; i < tau; i++) {
    if (i < params->faest_param.t0) {
      b = 0;
    } else {
      b = 1;
    }
    // TODO: the bits are represented as uint8,, hopefully its not too bad
    uint8_t* chalout;
    ChalDec(chal_3, i, params->faest_param.k0, params->faest_param.t0, params->faest_param.k1,
            params->faest_param.t1, chalout);
    d_tilde[i] = malloc(8);
    for (uint32_t j = 0; j < 8; j++) {
      d_tilde[i][j] = chalout[j] ^ signature->u_tilde[j];
    }
  }

  // Step: 17..18
  uint8_t* q_tilde_tmp = malloc(tau * 8);
  for (uint32_t i = 0; i < tau; i++) {
    for (uint8_t j = 0; j < 8; j++) {
      *(q_tilde_tmp + (i * tau) + j) = *(q_tilde[i] + j) ^ *(d_tilde[i] + j);
    }
  }
  uint8_t* hv_verify = malloc(lambdaBytes * 2);
  H1_context_t h1_ctx_1;
  H1_init(&h1_ctx_1, lambda);
  H1_update(&h1_ctx_1, q_tilde_tmp, tau * 8);
  H1_final(&h1_ctx_1, hv_verify, lambdaBytes * 2);
  if (memcmp(signature->h_v, hv_verify, lambdaBytes * 2) != 0) {
    return 0;
  }

  unsigned int beta, R, Nwd, Ske, Lke, Lenc, Senc, C;
  bool ret =
      aes_verify(signature->d, q, chal_2, chal_3, signature->a_tilde, signature->b_tilde, in, out,
                 lambda, tau, l, beta, R, Nwd, Ske, Lke, Lenc, Senc, C, params->faest_param.k0,
                 params->faest_param.k1, params->faest_param.t0, params->faest_param.t1);

  if (ret) {
    return 0;
  }
  return 1;
}

int serialize_signature(uint8_t* dst, size_t* len, const signature_t* signature,
                        const faest_paramset_t* params) {
  uint8_t* const old_dst    = dst;
  const unsigned int tau0   = params->faest_param.t0;
  const size_t lambda_bytes = params->faest_param.lambda / 8;
  const size_t ell_bytes    = (params->faest_param.l + 7) / 8;
  const size_t ell_hat =
      params->faest_param.l + params->faest_param.lambda * 2 + params->faest_param.b;
  const size_t ell_hat_bytes = (ell_hat + 7) / 8;
  const size_t utilde_bytes  = (params->faest_param.lambda + params->faest_param.b + 7) / 8;

  // serialize h_com (2 * \lambda)
  memcpy(dst, signature->hcom, lambda_bytes * 2);
  dst += lambda_bytes * 2;

  // serialize c_i
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    memcpy(dst, signature->c[i], ell_hat_bytes);
    dst += ell_hat_bytes;
  }

  // serialize u tilde
  memcpy(dst, signature->u_tilde, utilde_bytes);
  dst += utilde_bytes;

  // serialize h_v
  memcpy(dst, signature->h_v, 2 * lambda_bytes);
  dst += 2 * lambda_bytes;

  // serialize d
  memcpy(dst, signature->d, ell_bytes);
  dst += ell_bytes;

  // serialize a tilde
  memcpy(dst, signature->a_tilde, lambda_bytes);
  dst += lambda_bytes;

  // serialize b tilde
  memcpy(dst, signature->b_tilde, lambda_bytes);
  dst += lambda_bytes;

  // serialize pdec_i, com_i
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;

    memcpy(dst, signature->pdec[i], depth * lambda_bytes);
    dst += depth * lambda_bytes;
    memcpy(dst, signature->com_j[i], lambda_bytes);
    dst += lambda_bytes;
  }

  *len = dst - old_dst;
  return 0;
}

signature_t init_signature(const faest_paramset_t* params) {
  signature_t sig = {NULL};

  const unsigned int tau0   = params->faest_param.t0;
  const size_t lambda_bytes = params->faest_param.lambda / 8;
  const size_t ell_bytes    = (params->faest_param.l + 7) / 8;
  const size_t ell_hat =
      params->faest_param.l + params->faest_param.lambda * 2 + params->faest_param.b;
  const size_t ell_hat_bytes = ell_hat / 8;
  const size_t utilde_bytes  = (params->faest_param.lambda + params->faest_param.b + 7) / 8;

  sig.hcom = malloc(lambda_bytes * 2);
  sig.c    = calloc(params->faest_param.tau, sizeof(uint8_t*));
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    sig.c[i] = malloc(ell_hat_bytes);
  }
  sig.u_tilde = malloc(utilde_bytes);
  sig.h_v     = malloc(2 * lambda_bytes);
  sig.d       = malloc(ell_bytes);
  sig.a_tilde = malloc(lambda_bytes);
  sig.b_tilde = malloc(lambda_bytes);
  sig.pdec    = calloc(params->faest_param.tau, sizeof(uint8_t*));
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;
    sig.pdec[i]        = malloc(depth * lambda_bytes);
  }
  sig.com_j = calloc(params->faest_param.tau, sizeof(uint8_t*));
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    sig.com_j[i] = malloc(lambda_bytes);
  }

  return sig;
}

void free_signature(signature_t sig, const faest_paramset_t* params) {
  if (sig.com_j) {
    for (unsigned int i = params->faest_param.tau; i; --i) {
      free(sig.com_j[i - 1]);
    }
    free(sig.com_j);
  }
  if (sig.pdec) {
    for (unsigned int i = params->faest_param.tau; i; --i) {
      free(sig.pdec[i - 1]);
    }
    free(sig.pdec);
  }
  free(sig.b_tilde);
  free(sig.a_tilde);
  free(sig.d);
  free(sig.h_v);
  free(sig.u_tilde);

  if (sig.c) {
    for (unsigned int i = params->faest_param.tau; i; --i) {
      free(sig.c[i - 1]);
    }
    free(sig.c);
  }
  free(sig.hcom);
}

signature_t deserialize_signature(const uint8_t* src, size_t len, const faest_paramset_t* params) {
  const unsigned int tau0   = params->faest_param.t0;
  const size_t lambda_bytes = params->faest_param.lambda / 8;
  const size_t ell_bytes    = (params->faest_param.l + 7) / 8;
  const size_t ell_hat =
      params->faest_param.l + params->faest_param.lambda * 2 + params->faest_param.b;
  const size_t ell_hat_bytes = (ell_hat + 7) / 8;
  const size_t utilde_bytes  = (params->faest_param.lambda + params->faest_param.b + 7) / 8;

  signature_t sig = init_signature(params);

  // serialize h_com (2 * \lambda)
  memcpy(sig.hcom, src, lambda_bytes * 2);
  src += lambda_bytes * 2;

  // serialize c_i
  for (unsigned int i = 0; i != params->faest_param.tau; ++i, src += ell_hat_bytes) {
    memcpy(sig.c[i], src, ell_hat_bytes);
  }

  // serialize u tilde
  memcpy(sig.u_tilde, src, utilde_bytes);
  src += utilde_bytes;

  // serialize h_v
  memcpy(sig.h_v, src, 2 * lambda_bytes);
  src += 2 * lambda_bytes;

  // serialize d
  memcpy(sig.d, src, ell_bytes);
  src += ell_bytes;

  // serialize a tilde
  memcpy(sig.a_tilde, src, lambda_bytes);
  src += lambda_bytes;

  // serialize b tilde
  memcpy(sig.b_tilde, src, lambda_bytes);
  src += lambda_bytes;

  // serialize pdec_i, com_i
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;
    memcpy(sig.pdec[i], src, depth * lambda_bytes);
    src += depth * lambda_bytes;
    memcpy(sig.com_j[i], src, lambda_bytes);
    src += lambda_bytes;
  }

  return sig;
}
