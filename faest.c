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

#include <string.h>

// helpers to compute position in signature (sign)

ATTR_PURE static inline uint8_t* signature_c(uint8_t* base_ptr, unsigned int index,
                                             const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + index * ell_hat_bytes;
}

ATTR_PURE static inline uint8_t* signature_u_tilde(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes;
}

ATTR_PURE static inline uint8_t* signature_d(uint8_t* base_ptr, const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes;
}

ATTR_PURE static inline uint8_t* signature_a_tilde(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes;
}

ATTR_PURE static inline uint8_t* signature_pdec(uint8_t* base_ptr, unsigned int index,
                                                const faest_paramset_t* params) {
  const unsigned int tau0    = params->faest_param.tau0;
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int k0 = params->faest_param.k;
  const unsigned int k1  = (params->faest_param.tau0 != 0) ? k0 - 1 : k0;

  base_ptr +=
      (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + lambda_bytes;
  if (index < tau0) {
    return base_ptr + index * (k0 + 2) * lambda_bytes;
  } else {
    return base_ptr +
           ((index - tau0) * (k1 + 2) + tau0 * (k0 + 2)) *
               lambda_bytes;
  }
}

ATTR_PURE static inline uint8_t* signature_com(uint8_t* base_ptr, unsigned int index,
                                               const faest_paramset_t* params) {
  const unsigned int tau0    = params->faest_param.tau0;
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int k0 = params->faest_param.k;
  const unsigned int k1  = (params->faest_param.tau0 != 0) ? k0 - 1 : k0;

  base_ptr +=
      (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + lambda_bytes;
  if (index < tau0) {
    return base_ptr +
           (index * (k0 + 2) + k0) * lambda_bytes;
  } else {
    return base_ptr + ((index - tau0) * (k1 + 2) + k1 +
                       tau0 * (k0 + 2)) *
                          lambda_bytes;
  }
}

ATTR_PURE static inline uint8_t* signature_chall_3(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const unsigned int lambda_bytes = params->faest_param.lambda / 8;
  return base_ptr + params->faest_param.sig_size - IV_SIZE - lambda_bytes;
}

ATTR_PURE static inline uint8_t* signature_iv(uint8_t* base_ptr, const faest_paramset_t* params) {
  return base_ptr + params->faest_param.sig_size - IV_SIZE;
}

// helpers to compute position in signature (verify)

ATTR_PURE static inline const uint8_t* dsignature_c(const uint8_t* base_ptr, unsigned int index,
                                                    const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + index * ell_hat_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_u_tilde(const uint8_t* base_ptr,
                                                          const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_d(const uint8_t* base_ptr,
                                                    const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_a_tilde(const uint8_t* base_ptr,
                                                          const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_pdec(const uint8_t* base_ptr, unsigned int index,
                                                       const faest_paramset_t* params) {
  const unsigned int tau0    = params->faest_param.tau0;
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int k0 = params->faest_param.k;
  const unsigned int k1  = (params->faest_param.tau0 != 0) ? k0 - 1 : k0;

  base_ptr +=
      (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + lambda_bytes;
  if (index < tau0) {
    return base_ptr + index * (k0 + 2) * lambda_bytes;
  } else {
    return base_ptr +
           ((index - tau0) * (k1 + 2) + tau0 * (k0 + 2)) *
               lambda_bytes;
  }
}

ATTR_PURE static inline const uint8_t* dsignature_com(const uint8_t* base_ptr, unsigned int index,
                                                      const faest_paramset_t* params) {
  const unsigned int tau0    = params->faest_param.tau0;
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int ell_bytes     = params->faest_param.l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 2 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int k0 = params->faest_param.k;
  const unsigned int k1  = (params->faest_param.tau0 != 0) ? k0 - 1 : k0;

  base_ptr +=
      (params->faest_param.tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + lambda_bytes;
  if (index < tau0) {
    return base_ptr +
           (index * (k0 + 2) + k0) * lambda_bytes;
  } else {
    return base_ptr + ((index - tau0) * (k1 + 2) + k1 +
                       tau0 * (k0 + 2)) *
                          lambda_bytes;
  }
}

ATTR_PURE static inline const uint8_t* dsignature_chall_3(const uint8_t* base_ptr,
                                                          const faest_paramset_t* params) {
  const unsigned int lambda_bytes = params->faest_param.lambda / 8;
  return base_ptr + params->faest_param.sig_size - IV_SIZE - lambda_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_iv(const uint8_t* base_ptr,
                                                     const faest_paramset_t* params) {
  return base_ptr + params->faest_param.sig_size - IV_SIZE;
}

// Called in FAEST.Sign()::3
static void hash_2_0(uint8_t* mu, const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* msg, unsigned int owf_size, unsigned int msg_len, unsigned int lambda) {
  H2_context_t h2_ctx;
  H2_init(&h2_ctx, lambda);
  H2_update(&h2_ctx, owf_in, owf_size);
  H2_update(&h2_ctx, owf_out, owf_size);
  H2_update(&h2_ctx, msg, msg_len);
  H2_final(&h2_ctx, mu, 2 * lambda / 8);
}

// Called in FAEST.Sign()::3
static void hash_2_1(uint8_t* chall_1, const uint8_t* mu, const vec_com_t* vecCom, const uint8_t* c, const uint8_t* iv, unsigned int lambda, const faest_paramset_t* params) {
  H2_context_t h2_ctx;
  H2_init(&h2_ctx, lambda);

  H2_update(&h2_ctx, mu, lambda/8 * 2);
  H2_update(&h2_ctx, vecCom->com, (params->faest_param.L * (lambda/8) * 3)); // L * com_size (where com_size = lambda/8 * 3)  // TODO: check if len is correct
  H2_update(&h2_ctx, c, ((params->faest_param.l/8) + 2 * (lambda/8) + UNIVERSAL_HASH_B) * params->faest_param.tau); // TODO: check if len is correct
  H2_update(&h2_ctx, iv, 16);
  H2_final(&h2_ctx, chall_1, (5 * lambda / 8) + 8);
}

// Called in FAEST.Sign()::15
static void hash_2_2(uint8_t* chall_2, const uint8_t* chall_1, const uint8_t* u_tilde, const uint8_t* h_v, const uint8_t* d, unsigned int lambda, const faest_paramset_t* params) {
  H2_context_t h2_ctx;
  H2_init(&h2_ctx, lambda);

  H2_update(&h2_ctx, chall_1, (5 * lambda/8) + 8);
  H2_update(&h2_ctx, u_tilde, lambda/8 + UNIVERSAL_HASH_B);  // TODO: check if len is correct
  H2_update(&h2_ctx, h_v, lambda/8 * 2);
  H2_update(&h2_ctx, d, params->faest_param.l / 8); // ell_bytes
  H2_final(&h2_ctx, chall_2, (3 * lambda / 8) + 8);
}

// Called in FAEST.Sign()::21
static void hash_2_3(uint8_t* chall_3, const uint8_t* chall_2, const uint8_t* a_tilde_0, const uint8_t* a_tilde_1, const uint8_t* a_tilde_2, unsigned int ctr, unsigned int lambda, const faest_paramset_t* params) {
  H2_context_t h2_ctx;
  H2_init(&h2_ctx, lambda);

  H2_update(&h2_ctx, chall_2, (3 * lambda/8) + 8);
  H2_update(&h2_ctx, a_tilde_0, 999999999999999);  // TODO: FIX a_tilde_0 LEN
  H2_update(&h2_ctx, a_tilde_1, 999999999999999);  // TODO: FIX a_tilde_1 LEN
  H2_update(&h2_ctx, a_tilde_2, 999999999999999);  // TODO: FIX a_tilde_2 LEN
  H2_update(&h2_ctx, ctr, sizeof(unsigned int));
  H2_final(&h2_ctx, chall_3, lambda/8);
}

// Called in FAEST.Sign()::4
static void hash_3(uint8_t* r, uint8_t* sig, const uint8_t* owf_sk, const uint8_t* mu, const uint8_t* rho, unsigned int lambda, const faest_paramset_t* params) {
  H3_context_t h3_ctx;
  H3_init(&h3_ctx, lambda);
  H3_update(&h3_ctx, owf_sk, lambda/8);
  H3_update(&h3_ctx, mu, lambda/8 * 2);
  H3_update(&h3_ctx, rho, lambda/8);
  H3_final(&h3_ctx, r, lambda/8, signature_iv(sig, params));
}

// Called in FAEST.Sign()::5
static void hash_4(uint8_t* iv, const uint8_t* iv_pre, unsigned int lambda) {
  H4_context_t h4_ctx;
  H4_init(&h4_ctx, lambda);
  H4_update(&h4_ctx, iv_pre, 16);
  H4_final(&h4_ctx, iv);
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

// FAEST.Sign()
void faest_sign(uint8_t* sig, const uint8_t* msg, size_t msg_len, const uint8_t* owf_key,
                const uint8_t* owf_input, const uint8_t* owf_output, const uint8_t* rho,
                size_t rholen, const faest_paramset_t* params) {
  const unsigned int l             = params->faest_param.l;
  const unsigned int ell_bytes     = l / 8;
  const unsigned int lambda        = params->faest_param.lambda;
  const unsigned int lambdaBytes   = lambda / 8;
  const unsigned int tau           = params->faest_param.tau;
  const unsigned int tau0          = params->faest_param.tau0;
  const unsigned int tau1          = params->faest_param.tau1;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int k0 = params->faest_param.k;
  const unsigned int k1  = (params->faest_param.tau0 != 0) ? k0 - 1 : k0;
  const unsigned int w_grind = params->faest_param.w_grind;


  // ::1
  unsigned int ctr = 0;

  // ::2
  uint8_t rho[lambda/8]; 
  memset(rho, lambda/8, 0);  // setting it determinstic

  // ::3
  uint8_t mu[MAX_LAMBDA_BYTES * 2];
  unsigned int owf_size = params->faest_param.pk_size / 2;
  hash_2_0(mu, owf_input, owf_output, owf_size, msg, msg_len, lambda);    // hash function for the j-th Fiat-Shamir challenge; modeled as a random oracle

  // ::4
  uint8_t r[MAX_LAMBDA_BYTES];      // root key
  uint8_t iv_pre[128/8];             // pre iv used to generate the iv
  hash_3(r, sig, owf_key, mu, rho, lambda, params); // hash function for secret randomness derivation

  // ::5
  uint8_t iv[128];  // the real iv
  hash_4(iv, iv_pre, lambda); // hash function for IV derivation

  // TODO: For now Line number shifted by 1 due to gray comment

  // ::6-7
  vec_com_t* vecCom = calloc(tau, sizeof(vec_com_t));
  uint8_t* u        = malloc(ell_hat_bytes);
  // v has \hat \ell rows, \lambda columns, storing in column-major order
  uint8_t** V = malloc(lambda * sizeof(uint8_t*));
  V[0]        = calloc(lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    V[i] = V[0] + i * ell_hat_bytes;
  }
  vole_commit(r, signature_iv(sig, params), ell_hat, params, vecCom, signature_c(sig, 0, params), u, V);

  // ::8
  uint8_t chall_1[(5 * MAX_LAMBDA_BYTES) + 8];
  hash_2_1(chall_1, mu, vecCom, signature_c(sig, 0, params), iv, lambda, params);

  // ::9-10 // TODO: I will assume this works fine for now
  vole_hash(signature_u_tilde(sig, params), chall_1, u, l, lambda);

  // TODO: remove this comment
  // ::11-12 (H1 and H5 seems to be same, let's anyway name it H5 for consistency)
  uint8_t h_v[MAX_LAMBDA_BYTES * 2];
  {
    H5_context_t h5_ctx_1;
    H5_init(&h5_ctx_1, lambda);

    uint8_t V_tilde[MAX_LAMBDA_BYTES + UNIVERSAL_HASH_B];
    for (unsigned int i = 0; i != lambda; ++i) {
      vole_hash(V_tilde, chall_1, V[i], l, lambda);
      H5_update(&h5_ctx_1, V_tilde, lambda/8 + UNIVERSAL_HASH_B);
    }

    H5_final(&h5_ctx_1, h_v, lambda/8 * 2);
  }

  // ::13
  uint8_t* w = aes_extend_witness(owf_key, owf_input, params);
  // ::14
  xor_u8_array(w, u, signature_d(sig, params), ell_bytes);

  // ::15
  uint8_t chall_2[3 * MAX_LAMBDA_BYTES + 8];
  hash_2_2(chall_2, chall_1, signature_u_tilde(sig, params), h_v, signature_d(sig, params), lambda, params);

  // TODO: fix this a0, a1, a2 mess T.T
  // TODO: skipping for now
  // ::16-19
  uint8_t b_tilde[MAX_LAMBDA_BYTES];
  aes_prove(w, u, V, owf_input, owf_output, chall_2, signature_a_tilde(sig, params), b_tilde,
            params);
  free(V[0]);
  free(V);
  V = NULL;
  free(w);
  w = NULL;
  free(u);
  u = NULL;

  // ::20-27
  uint8_t chall_3[lambda/8];
  // TODO: For now putting this empty
  uint8_t* a_tilde_0;
  uint8_t* a_tilde_1;
  uint8_t* a_tilde_2;
  uint16_t* decoded_chall_3;
  do {
    hash_2_3(chall_3, chall_2, a_tilde_0, a_tilde_1, a_tilde_2, ctr, lambda, params);

    uint8_t flag = 0;
    for (unsigned int bit_i = lambda - w_grind; bit_i < lambda; bit_i++) {
      unsigned int byte_pos = bit_i/8;
      unsigned int bit_pos = bit_i%8;
      if (((chall_3[byte_pos] >> bit_pos) & 1) != 0) {
        ctr += 1;
        flag = 1;
        break;
      }
    }
    if (flag) { continue; }
    // I guess we do not care about the rest anyway if this fails? This could be fixed in the spec?

    bool flag_i = decode_all_chall_3(decoded_chall_3, chall_3, params);    
    // TODO: update the bavc.open

    break;

  } while (true);


  // // Step: 19..21
  // for (unsigned int i = 0; i < tau; i++) {
  //   // Step 20
  //   uint8_t s_[MAX_DEPTH];
  //   ChalDec(signature_chall_3(sig, params), i, k0, tau0,
  //           k1, tau1, s_);
  //   // Step 21
  //   const unsigned int depth = i < tau0 ? k0 : k1;
  //   bavc_open(vecCom, i, signature_pdec(sig, i, params), params); // TODO: I think this is fine or ?
  //   vec_com_clear(&vecCom[i]);
  // }

  free(vecCom);
  vecCom = NULL;
}

int faest_verify(const uint8_t* msg, size_t msglen, const uint8_t* sig, const uint8_t* owf_input,
                 const uint8_t* owf_output, const faest_paramset_t* params) {
  const unsigned int l             = params->faest_param.l;
  const unsigned int lambda        = params->faest_param.lambda;
  const unsigned int lambdaBytes   = lambda / 8;
  const unsigned int tau           = params->faest_param.tau;
  const unsigned int tau0          = params->faest_param.tau0;
  const unsigned int tau1          = params->faest_param.tau1;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int utilde_bytes  = lambdaBytes + UNIVERSAL_HASH_B;
  const unsigned int k0 = params->faest_param.k;
  const unsigned int k1  = (params->faest_param.tau0 != 0) ? k0 - 1 : k0;

  // Step: 3
  uint8_t mu[MAX_LAMBDA_BYTES * 2];
  hash_mu(mu, owf_input, owf_output, params->faest_param.pk_size / 2, msg, msglen, lambda);

  // Step: 5
  // q prime is a \hat \ell \times \lambda matrix
  uint8_t** qprime = malloc(lambda * sizeof(uint8_t*));
  qprime[0]        = calloc(lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    qprime[i] = qprime[0] + i * ell_hat_bytes;
  }
  uint8_t hcom[MAX_LAMBDA_BYTES * 2];
  {
    // TODO: change vole_reconstruct interface to make this less ugly
    const uint8_t* pdec[MAX_TAU];
    const uint8_t* com[MAX_TAU];
    for (unsigned int i = 0; i < tau; ++i) {
      pdec[i] = dsignature_pdec(sig, i, params);
      com[i]  = dsignature_com(sig, i, params);
    }
    vole_reconstruct(dsignature_iv(sig, params), dsignature_chall_3(sig, params), pdec, com, hcom,
                     qprime, ell_hat, params);
  }

  // Step: 5
  uint8_t chall_1[(5 * MAX_LAMBDA_BYTES) + 8];
  hash_challenge_1(chall_1, mu, hcom, dsignature_c(sig, 0, params), dsignature_iv(sig, params),
                   lambda, l, tau);

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
  for (unsigned int i = 0; i < tau; i++) {
    const unsigned int depth = i < tau0 ? k0 : k1;

    // Step 11
    uint8_t delta[MAX_DEPTH];
    ChalDec(dsignature_chall_3(sig, params), i, k0, tau0,
            k1, tau1, delta);
    // Step 16
    for (unsigned int j = 0; j != depth; ++j, ++Dtilde_idx) {
      // for scan-build
      assert(Dtilde_idx < lambda);
      masked_xor_u8_array(Dtilde[Dtilde_idx], dsignature_u_tilde(sig, params), Dtilde[Dtilde_idx],
                          delta[j], utilde_bytes);
    }

    if (i == 0) {
      // Step 8
      memcpy(q[q_idx], qprime[q_idx], ell_hat_bytes * depth);
      q_idx += depth;
    } else {
      // Step 14
      for (unsigned int d = 0; d < depth; ++d, ++q_idx) {
        masked_xor_u8_array(qprime[q_idx], dsignature_c(sig, i - 1, params), q[q_idx], delta[d],
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
      xor_u8_array(Q_tilde, Dtilde[i], Q_tilde, lambdaBytes + UNIVERSAL_HASH_B);
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
  hash_challenge_2(chall_2, chall_1, dsignature_u_tilde(sig, params), h_v,
                   dsignature_d(sig, params), lambda, l);

  // Step 18
  uint8_t* b_tilde =
      aes_verify(dsignature_d(sig, params), q, chall_2, dsignature_chall_3(sig, params),
                 dsignature_a_tilde(sig, params), owf_input, owf_output, params);
  free(q[0]);
  free(q);
  q = NULL;

  // Step: 20
  uint8_t chall_3[MAX_LAMBDA_BYTES];
  hash_challenge_3(chall_3, chall_2, dsignature_a_tilde(sig, params), b_tilde, lambda);
  free(b_tilde);
  b_tilde = NULL;

  // Step 21
  return memcmp(chall_3, dsignature_chall_3(sig, params), lambdaBytes) == 0 ? 0 : -1;
}

