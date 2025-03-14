/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest_impl.h"
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
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + index * ell_hat_bytes;
}

ATTR_PURE static inline uint8_t* signature_u_tilde(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->tau - 1) * ell_hat_bytes;
}

ATTR_PURE static inline uint8_t* signature_d(uint8_t* base_ptr, const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->tau - 1) * ell_hat_bytes + utilde_bytes;
}

ATTR_PURE static inline uint8_t* signature_a1_tilde(uint8_t* base_ptr,
                                                    const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes;
}

ATTR_PURE static inline uint8_t* signature_a2_tilde(uint8_t* base_ptr,
                                                    const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + lambda_bytes;
}

ATTR_PURE static inline uint8_t* signature_decom_i(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + 2 * lambda_bytes;
}

ATTR_PURE static inline uint8_t* signature_chall_3(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const unsigned int lambda_bytes = params->lambda / 8;
  return base_ptr + params->sig_size - sizeof(uint32_t) - IV_SIZE - lambda_bytes;
}

ATTR_PURE static inline uint8_t* signature_iv_pre(uint8_t* base_ptr,
                                                  const faest_paramset_t* params) {
  return base_ptr + params->sig_size - sizeof(uint32_t) - IV_SIZE;
}

ATTR_PURE static inline uint8_t* signature_ctr(uint8_t* base_ptr, const faest_paramset_t* params) {
  return base_ptr + params->sig_size - sizeof(uint32_t);
}

// helpers to compute position in signature (verify)

ATTR_PURE static inline const uint8_t* dsignature_c(const uint8_t* base_ptr, unsigned int index,
                                                    const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + index * ell_hat_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_u_tilde(const uint8_t* base_ptr,
                                                          const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->tau - 1) * ell_hat_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_d(const uint8_t* base_ptr,
                                                    const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->tau - 1) * ell_hat_bytes + utilde_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_a1_tilde(const uint8_t* base_ptr,
                                                           const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_a2_tilde(const uint8_t* base_ptr,
                                                           const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + lambda_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_decom_i(const uint8_t* base_ptr,
                                                          const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell_bytes     = params->l / 8;
  const unsigned int ell_hat_bytes = ell_bytes + 3 * lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  return base_ptr + (params->tau - 1) * ell_hat_bytes + utilde_bytes + ell_bytes + 2 * lambda_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_chall_3(const uint8_t* base_ptr,
                                                          const faest_paramset_t* params) {
  const unsigned int lambda_bytes = params->lambda / 8;
  return base_ptr + params->sig_size - sizeof(uint32_t) - IV_SIZE - lambda_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_iv_pre(const uint8_t* base_ptr,
                                                         const faest_paramset_t* params) {
  return base_ptr + params->sig_size - sizeof(uint32_t) - IV_SIZE;
}

ATTR_PURE static inline const uint8_t* dsignature_ctr(const uint8_t* base_ptr,
                                                      const faest_paramset_t* params) {
  return base_ptr + params->sig_size - sizeof(uint32_t);
}

// FAEST.Sign: line 3
static void hash_mu(uint8_t* mu, const uint8_t* owf_input, size_t owf_input_size,
                    const uint8_t* owf_output, size_t owf_output_size, const uint8_t* msg,
                    size_t msglen, unsigned int lambda) {
  H2_context_t h1_ctx;
  H2_init(&h1_ctx, lambda);
  H2_update(&h1_ctx, owf_input, owf_input_size);
  H2_update(&h1_ctx, owf_output, owf_output_size);
  H2_update(&h1_ctx, msg, msglen);
  H2_0_final(&h1_ctx, mu, 2 * lambda / 8);
}

static void hash_iv(uint8_t* iv, const uint8_t* iv_pre, unsigned int lambda) {
  H4_context_t h4_ctx;
  H4_init(&h4_ctx, lambda);
  H4_update(&h4_ctx, iv_pre);
  H4_final(&h4_ctx, iv);
}

// FAEST.Sign: line 4 + line 5
static void hash_r_iv(uint8_t* root_key, uint8_t* iv_pre, uint8_t* iv, const uint8_t* owf_key,
                      const uint8_t* mu, const uint8_t* rho, size_t rho_size, unsigned int lambda) {
  const unsigned int lambda_bytes = lambda / 8;

  {
    H3_context_t h3_ctx;
    H3_init(&h3_ctx, lambda);
    H3_update(&h3_ctx, owf_key, lambda_bytes);
    H3_update(&h3_ctx, mu, lambda_bytes * 2);
    if (rho && rho_size) {
      H3_update(&h3_ctx, rho, rho_size);
    }
    H3_final(&h3_ctx, root_key, lambda_bytes, iv_pre);
  }

  hash_iv(iv, iv_pre, lambda);
}

static void hash_challenge_1(uint8_t* chall_1, const uint8_t* mu, const uint8_t* hcom,
                             const uint8_t* c, const uint8_t* iv, unsigned int lambda,
                             unsigned int ell, unsigned int tau) {
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int ell_hat_bytes = ell / 8 + lambda_bytes * 3 + UNIVERSAL_HASH_B;

  H2_context_t h2_ctx;
  H2_init(&h2_ctx, lambda);
  H2_update(&h2_ctx, mu, lambda_bytes * 2);
  H2_update(&h2_ctx, hcom, lambda_bytes * 2);
  H2_update(&h2_ctx, c, ell_hat_bytes * (tau - 1));
  H2_update(&h2_ctx, iv, IV_SIZE);
  H2_1_final(&h2_ctx, chall_1, 5 * lambda_bytes + 8);
}

static void hash_challenge_2_init(H2_context_t* h2_ctx, const uint8_t* chall_1,
                                  const uint8_t* u_tilde, unsigned int lambda) {
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int u_tilde_bytes = lambda_bytes + UNIVERSAL_HASH_B;

  H2_init(h2_ctx, lambda);
  H2_update(h2_ctx, chall_1, 5 * lambda_bytes + 8);
  H2_update(h2_ctx, u_tilde, u_tilde_bytes);
}

static void hash_challenge_2_update_v_tilde(H2_context_t* h2_ctx, const uint8_t* v_tilde,
                                            unsigned int lambda) {
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int v_tilde_bytes = lambda_bytes + UNIVERSAL_HASH_B;

  H2_update(h2_ctx, v_tilde, v_tilde_bytes);
}

static void hash_challenge_2_finalize(uint8_t* chall_2, H2_context_t* h2_ctx, const uint8_t* d,
                                      const unsigned lambda, unsigned int ell) {
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ell_bytes    = ell / 8;

  H2_update(h2_ctx, d, ell_bytes);
  H2_2_final(h2_ctx, chall_2, 3 * lambda_bytes + 8);
}

static void hash_challenge_3_init(H2_context_t* h2_ctx, const uint8_t* chall_2,
                                  const uint8_t* a0_tilde, const uint8_t* a1_tilde,
                                  const uint8_t* a2_tilde, unsigned int lambda) {
  const unsigned int lambda_bytes = lambda / 8;

  H2_init(h2_ctx, lambda);
  H2_update(h2_ctx, chall_2, 3 * lambda_bytes + 8);
  H2_update(h2_ctx, a0_tilde, lambda_bytes);
  H2_update(h2_ctx, a1_tilde, lambda_bytes);
  H2_update(h2_ctx, a2_tilde, lambda_bytes);
}

static void hash_challenge_3_final(uint8_t* chall_3, const H2_context_t* ctx, uint32_t ctr,
                                   unsigned int lambda) {
  const unsigned int lambda_bytes = lambda / 8;

  H2_context_t ctx_copy;
  H2_copy(&ctx_copy, ctx);
  H2_update_u32_le(&ctx_copy, ctr);
  H2_3_final(&ctx_copy, chall_3, lambda_bytes);
}

static void hash_challenge_3(uint8_t* chall_3, const uint8_t* chall_2, const uint8_t* a0_tilde,
                             const uint8_t* a1_tilde, const uint8_t* a2_tilde, const uint8_t* ctr,
                             unsigned int lambda) {
  H2_context_t h2_ctx;
  hash_challenge_3_init(&h2_ctx, chall_2, a0_tilde, a1_tilde, a2_tilde, lambda);
  H2_update(&h2_ctx, ctr, sizeof(uint32_t));
  H2_3_final(&h2_ctx, chall_3, lambda / 8);
}

static bool check_challenge_3(const uint8_t* chall_3, unsigned int start, unsigned int lambda) {
  for (unsigned int bit_i = start; bit_i != lambda; ++bit_i) {
    if (ptr_get_bit(chall_3, bit_i)) {
      return false;
    }
  }
  return true;
}

static void free_pointer_array(uint8_t*** ptr) {
  free((*ptr)[0]);
  free(*ptr);
  *ptr = NULL;
}

// AES(-EM) OWF dispatchers
static inline void aes_prove(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde,
                             const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* owf_in,
                             const uint8_t* owf_out, const uint8_t* chall_2,
                             const faest_paramset_t* params) {
  switch (params->lambda) {
  case 256:
    aes_256_prover(a0_tilde, a1_tilde, a2_tilde, w, u, V, owf_in, owf_out, chall_2, params);
    break;
  case 192:
    aes_192_prover(a0_tilde, a1_tilde, a2_tilde, w, u, V, owf_in, owf_out, chall_2, params);
    break;
  default:
    aes_128_prover(a0_tilde, a1_tilde, a2_tilde, w, u, V, owf_in, owf_out, chall_2, params);
  }
}

static inline void aes_verify(uint8_t* a0_tilde, const uint8_t* d, uint8_t** Q,
                              const uint8_t* chall_2, const uint8_t* chall_3,
                              const uint8_t* a1_tilde, const uint8_t* a2_tilde,
                              const uint8_t* owf_in, const uint8_t* owf_out,
                              const faest_paramset_t* params) {
  switch (params->lambda) {
  case 256:
    aes_256_verifier(a0_tilde, d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
    break;
  case 192:
    aes_192_verifier(a0_tilde, d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
    break;
  default:
    aes_128_verifier(a0_tilde, d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
  }
}

// FAEST.Sign()
void faest_sign(uint8_t* sig, const uint8_t* msg, size_t msg_len, const uint8_t* owf_key,
                const uint8_t* owf_input, const uint8_t* owf_output, const uint8_t* witness,
                const uint8_t* rho, size_t rholen, const faest_paramset_t* params) {
  const unsigned int ell           = params->l;
  const unsigned int ell_bytes     = ell / 8;
  const unsigned int lambda        = params->lambda;
  const unsigned int tau           = params->tau;
  const unsigned int ell_hat       = ell + lambda * 3 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int w_grind       = params->w_grind;

  // ::3
  uint8_t mu[MAX_LAMBDA_BYTES * 2];
  hash_mu(mu, owf_input, params->owf_input_size, owf_output, params->owf_output_size, msg, msg_len,
          lambda);

  // ::4-5
  uint8_t rootkey[MAX_LAMBDA_BYTES], iv[IV_SIZE];
  hash_r_iv(rootkey, signature_iv_pre(sig, params), iv, owf_key, mu, rho, rholen, lambda);

  // ::6-7
  bavc_t bavc;
  uint8_t* u = malloc(ell_hat_bytes);
  assert(u);
  // v has \hat \ell rows, \lambda columns, storing in column-major order
  uint8_t** V = malloc(lambda * sizeof(uint8_t*));
  assert(V);
  V[0] = calloc(lambda, ell_hat_bytes);
  assert(V[0]);
  for (unsigned int i = 1; i < lambda; ++i) {
    V[i] = V[0] + i * ell_hat_bytes;
  }
  vole_commit(rootkey, iv, ell_hat, params, &bavc, signature_c(sig, 0, params), u, V);

  // ::8
  uint8_t chall_1[(5 * MAX_LAMBDA_BYTES) + 8];
  hash_challenge_1(chall_1, mu, bavc.h, signature_c(sig, 0, params), iv, lambda, ell, tau);

  // ::9-10
  vole_hash(signature_u_tilde(sig, params), chall_1, u, ell, lambda);

  // ::11-12
  // To save memory consumption, the chall_2 is computed in an
  // Init-Update-Finalize style as V_tilde is only fed into to the hash and not
  // used elsewhere.
  H2_context_t chall_2_ctx;
  hash_challenge_2_init(&chall_2_ctx, chall_1, signature_u_tilde(sig, params), lambda);
  {
    uint8_t V_tilde[MAX_LAMBDA_BYTES + UNIVERSAL_HASH_B];
    for (unsigned int i = 0; i != lambda; ++i) {
      // Step 11
      vole_hash(V_tilde, chall_1, V[i], ell, lambda);
      // Step 14
      hash_challenge_2_update_v_tilde(&chall_2_ctx, V_tilde, lambda);
    }
  }

  // ::13 witness provided by caller
  // ::14
  xor_u8_array(witness, u, signature_d(sig, params), ell_bytes);

  // :15
  uint8_t chall_2[3 * MAX_LAMBDA_BYTES + 8];
  hash_challenge_2_finalize(chall_2, &chall_2_ctx, signature_d(sig, params), lambda, ell);

  // ::16-20
  uint8_t a0_tilde[MAX_LAMBDA_BYTES];
  aes_prove(a0_tilde, signature_a1_tilde(sig, params), signature_a2_tilde(sig, params), witness,
            u + ell_bytes, V, owf_input, owf_output, chall_2, params);

  free_pointer_array(&V);
  free(u);
  u = NULL;

  // ::21-22
  H2_context_t chall_3_ctx;
  hash_challenge_3_init(&chall_3_ctx, chall_2, a0_tilde, signature_a1_tilde(sig, params),
                        signature_a2_tilde(sig, params), lambda);

  uint32_t ctr = 0;
  for (; true; ++ctr) {
    uint8_t* chall_3 = signature_chall_3(sig, params);
    hash_challenge_3_final(chall_3, &chall_3_ctx, ctr, lambda);
    // declassify chall_3 which is put into the signature
    faest_declassify(chall_3, lambda / 8);

    // ::23
    if (!check_challenge_3(chall_3, lambda - w_grind, lambda)) {
      continue;
    }

    // ::26
    uint16_t decoded_chall_3[MAX_TAU];
    if (!decode_all_chall_3(decoded_chall_3, chall_3, params)) {
      continue;
    }

    // :27
    if (bavc_open(signature_decom_i(sig, params), &bavc, decoded_chall_3, params)) {
      break;
    }
  }
  hash_clear(&chall_3_ctx);
  bavc_clear(&bavc);

  // copy counter to signature
  ctr = htole32(ctr);
  memcpy(signature_ctr(sig, params), &ctr, sizeof(ctr));
}

int faest_verify(const uint8_t* msg, size_t msglen, const uint8_t* sig, const uint8_t* owf_input,
                 const uint8_t* owf_output, const faest_paramset_t* params) {
  const unsigned int ell           = params->l;
  const unsigned int lambda        = params->lambda;
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int tau           = params->tau;
  const unsigned int ell_hat       = ell + lambda * 3 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;

  // ::4-5
  if (!check_challenge_3(dsignature_chall_3(sig, params), lambda - params->w_grind, lambda)) {
    return -1;
  }

  // ::2
  uint8_t mu[MAX_LAMBDA_BYTES * 2];
  hash_mu(mu, owf_input, params->owf_input_size, owf_output, params->owf_output_size, msg, msglen,
          lambda);

  // ::3
  uint8_t iv[IV_SIZE];
  hash_iv(iv, dsignature_iv_pre(sig, params), lambda);

  // Step: 6-7
  // q is a \hat \ell \times \lambda matrix
  uint8_t** q = malloc(lambda * sizeof(uint8_t*));
  assert(q);
  q[0] = calloc(lambda, ell_hat_bytes);
  assert(q[0]);
  for (unsigned int i = 1; i < lambda; ++i) {
    q[i] = q[0] + i * ell_hat_bytes;
  }
  uint8_t hcom[MAX_LAMBDA_BYTES * 2];

  if (!vole_reconstruct(hcom, q, iv, dsignature_chall_3(sig, params),
                        dsignature_decom_i(sig, params), dsignature_c(sig, 0, params), ell_hat,
                        params)) {
    free_pointer_array(&q);
    return -1;
  }

  // ::10
  uint8_t chall_1[5 * MAX_LAMBDA_BYTES + 8];
  hash_challenge_1(chall_1, mu, hcom, dsignature_c(sig, 0, params), iv, lambda, ell, tau);

  // Step 12, 14 and 15
  H2_context_t chall_2_ctx;
  hash_challenge_2_init(&chall_2_ctx, chall_1, dsignature_u_tilde(sig, params), lambda);
  {
    const uint8_t* chall_3 = dsignature_chall_3(sig, params);
    uint8_t Q_tilde[MAX_LAMBDA_BYTES + UNIVERSAL_HASH_B];
    for (unsigned int i = 0; i != lambda; ++i) {
      // Step 12
      vole_hash(Q_tilde, chall_1, q[i], ell, lambda);
      // Step 14
      if (ptr_get_bit(chall_3, i)) {
        xor_u8_array(Q_tilde, dsignature_u_tilde(sig, params), Q_tilde, utilde_bytes);
      }

      // Step 15
      hash_challenge_2_update_v_tilde(&chall_2_ctx, Q_tilde, lambda);
    }
  }

  // Step 15
  uint8_t chall_2[3 * MAX_LAMBDA_BYTES + 8];
  hash_challenge_2_finalize(chall_2, &chall_2_ctx, dsignature_d(sig, params), lambda, ell);

  // Step 18
  const uint8_t* d = dsignature_d(sig, params);
  uint8_t a0_tilde[MAX_LAMBDA_BYTES];
  aes_verify(a0_tilde, d, q, chall_2, dsignature_chall_3(sig, params),
             dsignature_a1_tilde(sig, params), dsignature_a2_tilde(sig, params), owf_input,
             owf_output, params);
  free_pointer_array(&q);

  // Step: 20
  uint8_t chall_3[MAX_LAMBDA_BYTES];
  hash_challenge_3(chall_3, chall_2, a0_tilde, dsignature_a1_tilde(sig, params),
                   dsignature_a2_tilde(sig, params), dsignature_ctr(sig, params), lambda);

  // Step 21
  return memcmp(chall_3, dsignature_chall_3(sig, params), lambda_bytes) == 0 ? 0 : -1;
}
