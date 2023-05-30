/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest.h"
#include "aes.h"
#include "faest_aes.h"
#include "randomness.h"
#include "random_oracle.h"
#include "utils.h"

// TODO: TEST EVERYTHING HERE !!!

// TODO: Do not pass lambdaBytes everywhere, compute it in the function....

// TODO: change q to Q where applicable

static inline uint8_t get_bit(const uint8_t* in, unsigned int index) {
  return (in[index / 8] >> (7 - index % 8)) & 1;
}

static inline void set_bit(uint8_t* dst, uint8_t in, unsigned int index) {
  dst[index / 8] |= in << (7 - index % 8);
}

uint8_t** column_to_row_major_and_shrink_V(uint8_t** v, unsigned int lambda, unsigned int ell) {
  assert(lambda % 8 == 0);
  const unsigned int lambda_bytes = lambda / 8;

  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  uint8_t** new_v = malloc((ell + lambda) * sizeof(uint8_t*));
  new_v[0]        = malloc(lambda * (ell + lambda));
  for (unsigned int i = 1; i < ell + lambda; ++i) {
    new_v[i] = new_v[0] + i * lambda_bytes;
  }

  for (unsigned int row = 0; row != ell + lambda; ++row) {
    for (unsigned int column = 0; column != lambda; ++column) {
      set_bit(new_v[row], column, get_bit(v[column], row));
    }
  }

  return new_v;
}

void sign(const uint8_t* msg, size_t msglen, const uint8_t* sk, const uint8_t* pk,
          const uint8_t* rho, size_t rholen, const faest_paramset_t* params,
          signature_t* signature) {
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
    if (rho && rholen) {
      H3_update(&h3_ctx, rho, rholen);
    }
    H3_final(&h3_ctx, rootkey, lambdaBytes);
  }

  // Step: 3
  uint8_t* hcom     = malloc(lambdaBytes * 2);
  vec_com_t* vecCom = calloc(tau, sizeof(vec_com_t));
  uint8_t* u        = malloc(ell_hat_bytes);
  // v has \hat \ell rows, \lambda columns, storing in column-major order
  uint8_t** v = malloc(lambda * sizeof(uint8_t*));
  v[0]        = malloc(lambda * ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    v[i] = v[0] + i * ell_hat_bytes;
  }
  voleCommit(rootkey, ell_hat, params, hcom, vecCom, signature->c, u, v);
  free(rootkey);
  rootkey = NULL;

  // Step: 4
  uint8_t* chall_1 = malloc((5 * lambdaBytes) + 8);
  {
    H2_context_t h2_ctx;
    H2_init(&h2_ctx, lambda);
    H2_update(&h2_ctx, mu, lambdaBytes * 2);
    H2_update(&h2_ctx, hcom, lambdaBytes * 2);
    for (unsigned int i = 0; i < (tau - 1); ++i) {
      H2_update(&h2_ctx, signature->c[i], ell_hat_bytes);
      // printUint8Arr("sign siganture.c", signature->c[i], 5, 1);
    }
    H2_final(&h2_ctx, chall_1, (5 * lambdaBytes) + 8);
  }
  free(hcom);
  hcom = NULL;
  free(mu);
  mu = NULL;

  // Step: 6
  vole_hash(signature->u_tilde, chall_1, u, l, lambda);

  // Step: 7 and 8
  uint8_t* h_v = malloc(lambdaBytes * 2);
  {
    H1_context_t h1_ctx_1;
    H1_init(&h1_ctx_1, lambda);

    uint8_t* V_tilde = malloc(lambdaBytes + UNIVERSAL_HASH_B);
    for (unsigned int i = 0; i != lambda; ++i) {
      // Step 7
      // printUint8Arr("sign V", v[i], 5, 1);
      vole_hash(V_tilde, chall_1, v[i], l, lambda);
      // Step 8
      // printf("%d ", i);
      // printUint8Arr("sign V_tilde", V_tilde, lambdaBytes + UNIVERSAL_HASH_B, 1);
      H1_update(&h1_ctx_1, V_tilde, lambdaBytes + UNIVERSAL_HASH_B);
    }
    free(V_tilde);

    // Step: 8
    H1_final(&h1_ctx_1, h_v, lambdaBytes * 2);
  }
  // Step: 9
  const uint8_t* in  = pk;
  const uint8_t* out = pk + params->faest_param.pkSize / 2;
  // Step: 10
  uint8_t* w = aes_extend_witness(sk, in, params);
  // Step: 11
  xorUint8Arr(w, u, signature->d, ell_bytes);

  // printUint8Arr("sign chall 1", chall_1, (5 * lambdaBytes) + 8, 1);
  // printUint8Arr("sign signature.u_tilde", signature->u_tilde, utilde_bytes, 1);
  printUint8Arr("sign h_v", h_v, 2 * lambdaBytes, 1);
  // printUint8Arr("sign signature.d", signature->d, ell_bytes, 1);

  // Step: 12
  uint8_t* chall_2 = malloc(3 * lambdaBytes + 8);
  {
    H2_context_t h2_ctx_1;
    H2_init(&h2_ctx_1, lambda);
    H2_update(&h2_ctx_1, chall_1, lambdaBytes);
    H2_update(&h2_ctx_1, signature->u_tilde, utilde_bytes);
    H2_update(&h2_ctx_1, h_v, 2 * lambdaBytes);
    H2_update(&h2_ctx_1, signature->d, ell_bytes);
    H2_final(&h2_ctx_1, chall_2, (3 * lambdaBytes) + 8);
  }
  free(chall_1);
  chall_1 = NULL;
  free(h_v);
  h_v = NULL;

  printUint8Arr("sign chall 2", chall_2, (3 * lambdaBytes) + 8, 1);

  // Step: 14..15
  {
    uint8_t** new_v = column_to_row_major_and_shrink_V(v, lambda, l);
    free(v[0]);
    free(v);
    v = new_v;
  }

  uint8_t* u_ = malloc((l + lambda) / 8);
  memcpy(u_, u, (l + lambda) / 8);

  // Step: 16
  uint8_t* b_tilde = malloc(lambdaBytes);
  aes_prove(w, u_, v, in, out, chall_2, signature->a_tilde, b_tilde, params);
  free(v[0]);
  free(v);
  v = NULL;
  free(u);
  u = NULL;
  free(w);
  w = NULL;
  free(u_);
  u_ = NULL;

  // Step: 17
  {
    H2_context_t h2_ctx_2;
    H2_init(&h2_ctx_2, lambda);
    H2_update(&h2_ctx_2, chall_2, 3 * lambdaBytes + 8);
    H2_update(&h2_ctx_2, signature->a_tilde, lambdaBytes);
    H2_update(&h2_ctx_2, b_tilde, lambdaBytes);
    H2_final(&h2_ctx_2, signature->chall_3, lambdaBytes);
  }
  free(b_tilde);
  b_tilde = NULL;
  free(chall_2);
  chall_2 = NULL;

  // Step: 19..21
  uint8_t* s_ = malloc(MAX(params->faest_param.k0, params->faest_param.k1));
  for (uint32_t i = 0; i < tau; i++) {
    // Step 20
    ChalDec(signature->chall_3, i, params->faest_param.k0, params->faest_param.t0,
            params->faest_param.k1, params->faest_param.t1, s_);
    // printUint8Arr("sign si", s_, MAX(params->faest_param.k0, params->faest_param.k1), 1);
    // Step 21
    const unsigned int num_vole_instances =
        i < tau0 ? (1 << params->faest_param.k0) : (1 << params->faest_param.k1);
    vector_open(vecCom[i].k, vecCom[i].com, s_, signature->pdec[i], signature->com_j[i],
                num_vole_instances, lambdaBytes);
    vec_com_clear(&vecCom[i]);
  }
  free(s_);
  s_ = NULL;
  free(vecCom);
  vecCom = NULL;
}

int verify(const uint8_t* msg, size_t msglen, const uint8_t* pk, const faest_paramset_t* params,
           const signature_t* signature) {
  const unsigned int l           = params->faest_param.l;
  const unsigned int ell_bytes   = (l + 7) / 8;
  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int lambdaBytes = lambda / 8;
  const unsigned int tau         = params->faest_param.tau;
  const unsigned int tau0        = params->faest_param.t0;
  const unsigned int ell_hat =
      params->faest_param.l + params->faest_param.lambda * 2 + params->faest_param.b;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
  const unsigned int utilde_bytes  = (params->faest_param.lambda + params->faest_param.b + 7) / 8;

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
  // q prime is a \hat \ell \times \lambda matrix
  uint8_t** qprime = malloc(lambda * sizeof(uint8_t*));
  qprime[0]        = malloc(lambda * ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    qprime[i] = qprime[0] + i * ell_hat_bytes;
  }
  uint8_t* hcom = malloc(lambdaBytes * 2);
  {
    voleReconstruct(signature->chall_3, signature->pdec, signature->com_j, hcom, qprime, ell_hat,
                    params);
  }

  // Step: 5
  uint8_t* chall_1 = malloc((5 * lambdaBytes) + 8);
  {
    H2_context_t h2_ctx;
    H2_init(&h2_ctx, lambda);
    H2_update(&h2_ctx, mu, lambdaBytes * 2);
    H2_update(&h2_ctx, hcom, lambdaBytes * 2);
    for (unsigned int i = 0; i < (tau - 1); ++i) {
      H2_update(&h2_ctx, signature->c[i], ell_hat_bytes);
    }
    H2_final(&h2_ctx, chall_1, (5 * lambdaBytes) + 8);
  }
  free(mu);
  mu = NULL;

  // Step: 8..14
  uint8_t** q = malloc(lambda * sizeof(uint8_t*));
  // TODO: is each uint8 a bit here ?
  q[0] = malloc(lambda * ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    q[i] = q[0] + i * ell_hat_bytes;
  }

  uint8_t** Dtilde = malloc(tau * sizeof(uint8_t*));
  Dtilde[0]        = calloc(lambda, (lambdaBytes + UNIVERSAL_HASH_B));
  for (unsigned int i = 1; i < lambda; ++i) {
    Dtilde[i] = Dtilde[0] + i * (lambdaBytes + UNIVERSAL_HASH_B);
  }

  const unsigned int max_depth = MAX(params->faest_param.k0, params->faest_param.k0);
  uint8_t* delta               = malloc(max_depth);
  unsigned int Dtilde_idx      = 0;
  for (uint32_t i = 0; i < tau; i++) {
    const unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;

    // Step 11
    ChalDec(signature->chall_3, i, params->faest_param.k0, params->faest_param.t0,
            params->faest_param.k1, params->faest_param.t1, delta);
    // printUint8Arr("verify si", delta, MAX(params->faest_param.k0, params->faest_param.k1), 1);
    // Step 16
    for (unsigned int j = 0; j != depth; ++j, ++Dtilde_idx) {
      if (delta[j]) {
        memcpy(Dtilde[Dtilde_idx], signature->u_tilde, utilde_bytes);
      }
    }

    if (!i) {
      // Step 8
      memcpy(q[i], qprime[i], ell_hat_bytes * depth);
    } else {
      // Step 14
      // printUint8Arr("verify siganture.c", signature->c[i - 1], 5, 1);
      for (uint32_t d = 0; d < depth; d++) {
        if (delta[d]) {
          xorUint8Arr(qprime[i + d], signature->c[i - 1], q[i + d], ell_hat_bytes);
        } else {
          memcpy(q[i + d], qprime[i + d], ell_hat_bytes);
        }
      }
    }
  }

  // Step 15 and 16
  uint8_t* h_v = malloc(lambdaBytes * 2);
  {
    H1_context_t h1_ctx_1;
    H1_init(&h1_ctx_1, lambda);

    uint8_t* Q_tilde = malloc(lambdaBytes + UNIVERSAL_HASH_B);
    for (unsigned int i = 0; i != lambda; ++i) {
      // printUint8Arr("verify Q", q[0] + i * ell_hat_bytes, 5, 1);
      // Step 15
      vole_hash(Q_tilde, chall_1, q[0] + i * ell_hat_bytes, l, lambda);
      // Step 16
      xorUint8Arr(Q_tilde, Dtilde[i], Q_tilde, lambdaBytes + UNIVERSAL_HASH_B);
      H1_update(&h1_ctx_1, Q_tilde, lambdaBytes + UNIVERSAL_HASH_B);
    }
    free(Q_tilde);

    // Step: 16
    H1_final(&h1_ctx_1, h_v, lambdaBytes * 2);
  }

  // printUint8Arr("verify chall 1", chall_1, (5 * lambdaBytes) + 8, 1);
  // printUint8Arr("verify signature.u_tilde", signature->u_tilde, utilde_bytes, 1);
  printUint8Arr("verify h_v", h_v, 2 * lambdaBytes, 1);
  // printUint8Arr("verify signature.d", signature->d, ell_bytes, 1);

  // Step 17
  uint8_t* chall_2 = malloc(3 * lambdaBytes + 8);
  {
    H2_context_t h2_ctx_1;
    H2_init(&h2_ctx_1, lambda);
    H2_update(&h2_ctx_1, chall_1, lambdaBytes);
    H2_update(&h2_ctx_1, signature->u_tilde, utilde_bytes);
    H2_update(&h2_ctx_1, h_v, 2 * lambdaBytes);
    H2_update(&h2_ctx_1, signature->d, ell_bytes);
    H2_final(&h2_ctx_1, chall_2, (3 * lambdaBytes) + 8);
  }
  free(chall_1);
  chall_1 = NULL;

  printUint8Arr("verify chall 2", chall_2, (3 * lambdaBytes) + 8, 1);

  // Step 18
  // TODO: Do we transpose and shorten here before passing to q ?
  uint8_t* b_tilde =
      aes_verify(signature->d, q, chall_2, signature->chall_3, signature->a_tilde, in, out, params);

  // Step: 20
  uint8_t* chall_3 = malloc(lambdaBytes);
  {
    H2_context_t h2_ctx_2;
    H2_init(&h2_ctx_2, lambda);
    H2_update(&h2_ctx_2, chall_2, 3 * lambdaBytes + 8);
    H2_update(&h2_ctx_2, signature->a_tilde, lambdaBytes);
    H2_update(&h2_ctx_2, b_tilde, lambdaBytes);
    H2_final(&h2_ctx_2, chall_3, lambdaBytes);
  }
  free(b_tilde);
  b_tilde = NULL;
  free(chall_2);
  chall_2 = NULL;

  // Step 21
  int ret = memcmp(chall_3, signature->chall_3, lambdaBytes);
  free(chall_3);

  return ret == 0 ? 1 : 0;
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

  sig.c = calloc(params->faest_param.tau - 1, sizeof(uint8_t*));
  for (unsigned int i = 0; i != params->faest_param.tau - 1; ++i) {
    sig.c[i] = malloc(ell_hat_bytes);
  }
  sig.u_tilde = malloc(utilde_bytes);
  sig.d       = malloc(ell_bytes);
  sig.a_tilde = malloc(lambda_bytes);
  sig.pdec    = calloc(params->faest_param.tau, sizeof(uint8_t*));
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;
    sig.pdec[i]        = malloc(depth * lambda_bytes);
  }
  sig.com_j = calloc(params->faest_param.tau, sizeof(uint8_t*));
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    sig.com_j[i] = malloc(lambda_bytes * 2);
  }
  sig.chall_3 = malloc(lambda_bytes);

  return sig;
}

void free_signature(signature_t sig, const faest_paramset_t* params) {
  free(sig.chall_3);
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
  free(sig.a_tilde);
  free(sig.d);
  free(sig.u_tilde);

  if (sig.c) {
    for (unsigned int i = params->faest_param.tau - 1; i; --i) {
      free(sig.c[i - 1]);
    }
    free(sig.c);
  }
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

  // serialize c_i
  for (unsigned int i = 0; i < params->faest_param.tau - 1; ++i) {
    memcpy(dst, signature->c[i], ell_hat_bytes);
    dst += ell_hat_bytes;
  }

  // serialize u tilde
  memcpy(dst, signature->u_tilde, utilde_bytes);
  dst += utilde_bytes;

  // serialize d
  memcpy(dst, signature->d, ell_bytes);
  dst += ell_bytes;

  // serialize a tilde
  memcpy(dst, signature->a_tilde, lambda_bytes);
  dst += lambda_bytes;

  // serialize pdec_i, com_i
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    const unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;
    memcpy(dst, signature->pdec[i], depth * lambda_bytes);
    dst += depth * lambda_bytes;
    memcpy(dst, signature->com_j[i], 2 * lambda_bytes);
    dst += 2 * lambda_bytes;
  }

  // serialize chall_3
  memcpy(dst, signature->chall_3, lambda_bytes);
  dst += lambda_bytes;

  *len = dst - old_dst;
  return 0;
}

signature_t deserialize_signature(const uint8_t* src, const faest_paramset_t* params) {
  const unsigned int tau0   = params->faest_param.t0;
  const size_t lambda_bytes = params->faest_param.lambda / 8;
  const size_t ell_bytes    = (params->faest_param.l + 7) / 8;
  const size_t ell_hat =
      params->faest_param.l + params->faest_param.lambda * 2 + params->faest_param.b;
  const size_t ell_hat_bytes = (ell_hat + 7) / 8;
  const size_t utilde_bytes  = (params->faest_param.lambda + params->faest_param.b + 7) / 8;

  signature_t sig = init_signature(params);

  // serialize c_i
  for (unsigned int i = 0; i != params->faest_param.tau - 1; ++i, src += ell_hat_bytes) {
    memcpy(sig.c[i], src, ell_hat_bytes);
  }

  // serialize u tilde
  memcpy(sig.u_tilde, src, utilde_bytes);
  src += utilde_bytes;

  // serialize d
  memcpy(sig.d, src, ell_bytes);
  src += ell_bytes;

  // serialize a tilde
  memcpy(sig.a_tilde, src, lambda_bytes);
  src += lambda_bytes;

  // serialize pdec_i, com_i
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    const unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;
    memcpy(sig.pdec[i], src, depth * lambda_bytes);
    src += depth * lambda_bytes;
    memcpy(sig.com_j[i], src, 2 * lambda_bytes);
    src += 2 * lambda_bytes;
  }

  // serialize chall_3
  memcpy(sig.chall_3, src, lambda_bytes);
  // src += lambda_bytes;

  return sig;
}
