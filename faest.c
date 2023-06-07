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
#include "vole.h"
#include "universal_hashing.h"

static void hash_mu(uint8_t* mu, const uint8_t* owf_input, const uint8_t* owf_output,
                    size_t owf_size, const uint8_t* msg, size_t msglen, unsigned int lambda) {
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  H1_update(&h1_ctx, owf_input, owf_size);
  H1_update(&h1_ctx, owf_output, owf_size);
  H1_update(&h1_ctx, msg, msglen);
  H1_final(&h1_ctx, mu, 2 * lambda / 8);
}

static void hash_challenge_1(uint8_t* chall_1, const uint8_t* mu, const uint8_t* hcom,
                             const uint8_t* c, const uint8_t* iv, unsigned int lambda,
                             unsigned int ell, unsigned int tau) {
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int ell_hat_bytes = ell / 8 + lambda_bytes * 2 + UNIVERSAL_HASH_B;

  H2_context_t h2_ctx;
  H2_init(&h2_ctx, lambda);
  H2_update(&h2_ctx, mu, lambda_bytes * 2);
  H2_update(&h2_ctx, hcom, lambda_bytes * 2);
  H2_update(&h2_ctx, c, ell_hat_bytes * (tau - 1));
  H2_update(&h2_ctx, iv, IV_SIZE);
  H2_final(&h2_ctx, chall_1, 5 * lambda_bytes + 8);
}

static void hash_challenge_2(uint8_t* chall_2, const uint8_t* chall_1, const uint8_t* u_tilde,
                             const uint8_t* h_v, const uint8_t* d, unsigned int lambda,
                             unsigned int ell) {
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int ell_bytes     = ell / 8;
  const unsigned int u_tilde_bytes = lambda_bytes + UNIVERSAL_HASH_B;

  H2_context_t h2_ctx_1;
  H2_init(&h2_ctx_1, lambda);
  H2_update(&h2_ctx_1, chall_1, 5 * lambda_bytes + 8);
  H2_update(&h2_ctx_1, u_tilde, u_tilde_bytes);
  H2_update(&h2_ctx_1, h_v, 2 * lambda_bytes);
  H2_update(&h2_ctx_1, d, ell_bytes);
  H2_final(&h2_ctx_1, chall_2, 3 * lambda_bytes + 8);
}

static void hash_challenge_3(uint8_t* chall_3, const uint8_t* chall_2, const uint8_t* a_tilde,
                             const uint8_t* b_tilde, unsigned int lambda) {
  const unsigned int lambda_bytes = lambda / 8;

  H2_context_t h2_ctx_2;
  H2_init(&h2_ctx_2, lambda);
  H2_update(&h2_ctx_2, chall_2, 3 * lambda_bytes + 8);
  H2_update(&h2_ctx_2, a_tilde, lambda_bytes);
  H2_update(&h2_ctx_2, b_tilde, lambda_bytes);
  H2_final(&h2_ctx_2, chall_3, lambda_bytes);
}

void sign(const uint8_t* msg, size_t msglen, const uint8_t* owf_key, const uint8_t* owf_input,
          const uint8_t* owf_output, const uint8_t* rho, size_t rholen,
          const faest_paramset_t* params, signature_t* signature) {
  const uint32_t l             = params->faest_param.l;
  const uint32_t ell_bytes     = l / 8;
  const uint32_t lambda        = params->faest_param.lambda;
  const uint32_t lambdaBytes   = lambda / 8;
  const uint32_t tau           = params->faest_param.tau;
  const uint32_t tau0          = params->faest_param.t0;
  const uint32_t ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const uint32_t ell_hat_bytes = ell_hat / 8;

  // Step: 2
  uint8_t mu[MAX_LAMBDA_BYTES * 2];
  hash_mu(mu, owf_input, owf_output, params->faest_param.pkSize / 2, msg, msglen, lambda);

  // Step: 3
  uint8_t rootkey[MAX_LAMBDA_BYTES];
  {
    H3_context_t h3_ctx;
    H3_init(&h3_ctx, lambda);
    H3_update(&h3_ctx, owf_key, lambdaBytes);
    H3_update(&h3_ctx, mu, lambdaBytes * 2);
    if (rho && rholen) {
      H3_update(&h3_ctx, rho, rholen);
    }
    H3_final(&h3_ctx, rootkey, lambdaBytes, signature->iv);
  }

  // Step: 3
  uint8_t hcom[MAX_LAMBDA_BYTES * 2];
  vec_com_t* vecCom = calloc(tau, sizeof(vec_com_t));
  uint8_t* u        = malloc(ell_hat_bytes);
  // v has \hat \ell rows, \lambda columns, storing in column-major order
  uint8_t** V = malloc(lambda * sizeof(uint8_t*));
  V[0]        = calloc(lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    V[i] = V[0] + i * ell_hat_bytes;
  }
  voleCommit(rootkey, signature->iv, ell_hat, params, hcom, vecCom, signature->c, u, V);

  // Step: 4
  uint8_t chall_1[(5 * MAX_LAMBDA_BYTES) + 8];
  hash_challenge_1(chall_1, mu, hcom, signature->c, signature->iv, lambda, l, tau);

  // Step: 6
  vole_hash(signature->u_tilde, chall_1, u, l, lambda);

  // Step: 7 and 8
  uint8_t h_v[MAX_LAMBDA_BYTES * 2];
  {
    H1_context_t h1_ctx_1;
    H1_init(&h1_ctx_1, lambda);

    uint8_t V_tilde[MAX_LAMBDA_BYTES + UNIVERSAL_HASH_B];
    for (unsigned int i = 0; i != lambda; ++i) {
      // Step 7
      vole_hash(V_tilde, chall_1, V[i], l, lambda);
      // Step 8
      H1_update(&h1_ctx_1, V_tilde, lambdaBytes + UNIVERSAL_HASH_B);
    }
    // Step: 8
    H1_final(&h1_ctx_1, h_v, lambdaBytes * 2);
  }
  // Step: 9, 10
  uint8_t* w = aes_extend_witness(owf_key, owf_input, params);
  // Step: 11
  xorUint8Arr(w, u, signature->d, ell_bytes);

  // Step: 12
  uint8_t chall_2[3 * MAX_LAMBDA_BYTES + 8];
  hash_challenge_2(chall_2, chall_1, signature->u_tilde, h_v, signature->d, lambda, l);

  // Step: 14..15
  // transpose is computed in aes_prove

  // Step: 16
  uint8_t b_tilde[MAX_LAMBDA_BYTES];
  aes_prove(w, u, V, owf_input, owf_output, chall_2, signature->a_tilde, b_tilde, params);
  free(V[0]);
  free(V);
  V = NULL;
  free(w);
  w = NULL;
  free(u);
  u = NULL;

  // Step: 17
  hash_challenge_3(signature->chall_3, chall_2, signature->a_tilde, b_tilde, lambda);

  // Step: 19..21
  for (uint32_t i = 0; i < tau; i++) {
    // Step 20
    uint8_t s_[MAX_DEPTH];
    ChalDec(signature->chall_3, i, params->faest_param.k0, params->faest_param.t0,
            params->faest_param.k1, params->faest_param.t1, s_);
    // Step 21
    const unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;
    vector_open(vecCom[i].k, vecCom[i].com, s_, signature->pdec[i], signature->com_j[i], depth,
                lambdaBytes);
    vec_com_clear(&vecCom[i]);
  }
  free(vecCom);
  vecCom = NULL;
}

int verify(const uint8_t* msg, size_t msglen, const uint8_t* owf_input, const uint8_t* owf_output,
           const faest_paramset_t* params, const deserialized_signature_t* signature) {
  const unsigned int l             = params->faest_param.l;
  const unsigned int lambda        = params->faest_param.lambda;
  const unsigned int lambdaBytes   = lambda / 8;
  const unsigned int tau           = params->faest_param.tau;
  const unsigned int tau0          = params->faest_param.t0;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int utilde_bytes  = lambdaBytes + UNIVERSAL_HASH_B;
  const unsigned int k0            = params->faest_param.k0;
  const unsigned int k1            = params->faest_param.k1;

  // Step: 3
  uint8_t mu[MAX_LAMBDA_BYTES * 2];
  hash_mu(mu, owf_input, owf_output, params->faest_param.pkSize / 2, msg, msglen, lambda);

  // Step: 5
  // q prime is a \hat \ell \times \lambda matrix
  uint8_t** qprime = malloc(lambda * sizeof(uint8_t*));
  qprime[0]        = calloc(lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    qprime[i] = qprime[0] + i * ell_hat_bytes;
  }
  uint8_t hcom[MAX_LAMBDA_BYTES * 2];
  voleReconstruct(signature->iv, signature->chall_3, signature->pdec, signature->com_j, hcom,
                  qprime, ell_hat, params);

  // Step: 5
  uint8_t chall_1[(5 * MAX_LAMBDA_BYTES) + 8];
  hash_challenge_1(chall_1, mu, hcom, signature->c, signature->iv, lambda, l, tau);

  // Step: 8..14
  uint8_t** q = malloc(lambda * sizeof(uint8_t*));
  q[0]        = calloc(lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    q[i] = q[0] + i * ell_hat_bytes;
  }

  uint8_t** Dtilde = malloc(lambda * sizeof(uint8_t*));
  Dtilde[0]        = calloc(lambda, (lambdaBytes + UNIVERSAL_HASH_B));
  for (unsigned int i = 1; i < lambda; ++i) {
    Dtilde[i] = Dtilde[0] + i * (lambdaBytes + UNIVERSAL_HASH_B);
  }

  unsigned int Dtilde_idx = 0;
  unsigned int q_idx      = 0;
  for (uint32_t i = 0; i < tau; i++) {
    const unsigned int depth = i < tau0 ? k0 : k1;

    // Step 11
    uint8_t delta[MAX_DEPTH];
    ChalDec(signature->chall_3, i, params->faest_param.k0, params->faest_param.t0,
            params->faest_param.k1, params->faest_param.t1, delta);
    // Step 16
    for (unsigned int j = 0; j != depth; ++j, ++Dtilde_idx) {
      maskedXorUint8Arr(Dtilde[Dtilde_idx], signature->u_tilde, Dtilde[Dtilde_idx], delta[j],
                        utilde_bytes);
    }

    if (i == 0) {
      // Step 8
      memcpy(q[q_idx], qprime[q_idx], ell_hat_bytes * depth);
      q_idx += depth;
    } else {
      // Step 14
      for (unsigned int d = 0; d < depth; ++d, ++q_idx) {
        maskedXorUint8Arr(qprime[q_idx], signature->c + (i - 1) * ell_hat_bytes, q[q_idx], delta[d],
                          ell_hat_bytes);
      }
    }
  }
  free(qprime[0]);
  free(qprime);
  qprime = NULL;

  // Step 15 and 16
  uint8_t h_v[MAX_LAMBDA_BYTES * 2];
  {
    H1_context_t h1_ctx_1;
    H1_init(&h1_ctx_1, lambda);

    uint8_t Q_tilde[MAX_LAMBDA_BYTES + UNIVERSAL_HASH_B];
    for (unsigned int i = 0; i != lambda; ++i) {
      // Step 15
      vole_hash(Q_tilde, chall_1, q[i], l, lambda);
      // Step 16
      xorUint8Arr(Q_tilde, Dtilde[i], Q_tilde, lambdaBytes + UNIVERSAL_HASH_B);
      H1_update(&h1_ctx_1, Q_tilde, lambdaBytes + UNIVERSAL_HASH_B);
    }
    // Step: 16
    H1_final(&h1_ctx_1, h_v, lambdaBytes * 2);
  }
  free(Dtilde[0]);
  free(Dtilde);
  Dtilde = NULL;

  // Step 17
  uint8_t chall_2[3 * MAX_LAMBDA_BYTES + 8];
  hash_challenge_2(chall_2, chall_1, signature->u_tilde, h_v, signature->d, lambda, l);

  // Step 18
  uint8_t* b_tilde = aes_verify(signature->d, q, chall_2, signature->chall_3, signature->a_tilde,
                                owf_input, owf_output, params);
  free(q[0]);
  free(q);
  q = NULL;

  // Step: 20
  uint8_t chall_3[MAX_LAMBDA_BYTES];
  hash_challenge_3(chall_3, chall_2, signature->a_tilde, b_tilde, lambda);
  free(b_tilde);
  b_tilde = NULL;

  // Step 21
  return memcmp(chall_3, signature->chall_3, lambdaBytes) == 0 ? 0 : -1;
}

signature_t init_signature(const faest_paramset_t* params) {
  signature_t sig;
  memset(&sig, 0, sizeof(sig));

  const unsigned int tau0    = params->faest_param.t0;
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  sig.c       = malloc(ell_hat_bytes * (params->faest_param.tau - 1));
  sig.u_tilde = malloc(utilde_bytes);
  sig.d       = malloc(ell_bytes);
  sig.a_tilde = malloc(lambda_bytes);
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;
    sig.pdec[i]        = malloc(depth * lambda_bytes);
  }
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    sig.com_j[i] = malloc(lambda_bytes * 2);
  }
  sig.chall_3 = malloc(lambda_bytes);

  return sig;
}

void free_signature(signature_t sig, const faest_paramset_t* params) {
  free(sig.chall_3);
  for (unsigned int i = params->faest_param.tau; i; --i) {
    free(sig.com_j[i - 1]);
  }
  for (unsigned int i = params->faest_param.tau; i; --i) {
    free(sig.pdec[i - 1]);
  }
  free(sig.a_tilde);
  free(sig.d);
  free(sig.u_tilde);
  free(sig.c);
}

int serialize_signature(uint8_t* dst, const signature_t* signature,
                        const faest_paramset_t* params) {
  const unsigned int tau0    = params->faest_param.t0;
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  // serialize c_i
  memcpy(dst, signature->c, ell_hat_bytes * (params->faest_param.tau - 1));
  dst += ell_hat_bytes * (params->faest_param.tau - 1);

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

  // serialize iv
  memcpy(dst, signature->iv, sizeof(signature->iv));
  dst += sizeof(signature->iv);

  return 0;
}

deserialized_signature_t deserialize_signature(const uint8_t* src, const faest_paramset_t* params) {
  const unsigned int tau0    = params->faest_param.t0;
  const size_t lambda_bytes  = params->faest_param.lambda / 8;
  const size_t ell_bytes     = params->faest_param.l / 8;
  const size_t ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const size_t utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  deserialized_signature_t sig;

  // serialize c_i
  sig.c = src;
  src += ell_hat_bytes * (params->faest_param.tau - 1);

  // serialize u tilde
  sig.u_tilde = src;
  src += utilde_bytes;

  // serialize d
  sig.d = src;
  src += ell_bytes;

  // serialize a tilde
  sig.a_tilde = src;
  src += lambda_bytes;

  // serialize pdec_i, com_i
  for (unsigned int i = 0; i != params->faest_param.tau; ++i) {
    const unsigned int depth = i < tau0 ? params->faest_param.k0 : params->faest_param.k1;

    sig.pdec[i] = src;
    src += depth * lambda_bytes;

    sig.com_j[i] = src;
    src += 2 * lambda_bytes;
  }

  // serialize chall_3
  sig.chall_3 = src;
  src += lambda_bytes;

  // serialize iv
  sig.iv = src;

  return sig;
}
