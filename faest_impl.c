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
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;

  return base_ptr 
          + index * ell_bytes;
}

ATTR_PURE static inline uint8_t* signature_c_mult(uint8_t* base_ptr,
                                             const faest_paramset_t* params) {
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;

  return base_ptr 
          + ((params->tau - 1) * ell_bytes);
}

ATTR_PURE static inline uint8_t* signature_u_tilde(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;

  return base_ptr 
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes;
}

ATTR_PURE static inline uint8_t* signature_d(uint8_t* base_ptr, const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;

  return base_ptr 
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes 
          + utilde_bytes;
}

// We do not need to store a0_tilde in the signature
ATTR_PURE static inline uint8_t* signature_a1toi_tilde(uint8_t* base_ptr,
                                                    const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;

  return base_ptr 
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes 
          + utilde_bytes 
          + ell_bytes;
}

ATTR_PURE static inline uint8_t* signature_decom_i(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;
  const unsigned int d_zk           = params->d_zk;

  return base_ptr 
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes
          + utilde_bytes 
          + ell_bytes 
          + (d_zk - 1) * lambda_bytes;
}

ATTR_PURE static inline uint8_t* signature_chall_3(uint8_t* base_ptr,
                                                   const faest_paramset_t* params) {
  const faest_paramset_t* faest_param = faest_get_paramset(params->id);
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;
  const unsigned int d_zk           = params->d_zk;
  const unsigned int n_leafcom              = faest_is_em(faest_param) ? 2 : 3;
  const unsigned int hcom_size      = lambda_bytes * n_leafcom;
  const unsigned int tau            = params->tau;
  const unsigned int t_open         = params->T_open;

  return base_ptr 
          // + params->sig_size 
          // - sizeof(uint32_t) 
          // - IV_SIZE 
          // - lambda_bytes;
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes
          + utilde_bytes 
          + ell_bytes 
          + (d_zk - 1) * lambda_bytes
          + hcom_size * tau + t_open * lambda_bytes;
}

ATTR_PURE static inline uint8_t* signature_iv_pre(uint8_t* base_ptr,
                                                  const faest_paramset_t* params) {
  const faest_paramset_t* faest_param = faest_get_paramset(params->id);
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;
  const unsigned int d_zk           = params->d_zk;
  const unsigned int n_leafcom              = faest_is_em(faest_param) ? 2 : 3;
  const unsigned int hcom_size      = lambda_bytes * n_leafcom;
  const unsigned int tau            = params->tau;
  const unsigned int t_open         = params->T_open;

  return base_ptr 
          // + params->sig_size 
          // - sizeof(uint32_t) 
          // - IV_SIZE;
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes
          + utilde_bytes 
          + ell_bytes 
          + (d_zk - 1) * lambda_bytes
          + hcom_size * tau + t_open * lambda_bytes
          + lambda_bytes;
}

ATTR_PURE static inline uint8_t* signature_ctr(uint8_t* base_ptr, const faest_paramset_t* params) {

  const faest_paramset_t* faest_param = faest_get_paramset(params->id);
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;
  const unsigned int d_zk           = params->d_zk;
  const unsigned int n_leafcom              = faest_is_em(faest_param) ? 2 : 3;
  const unsigned int hcom_size      = lambda_bytes * n_leafcom;
  const unsigned int tau            = params->tau;
  const unsigned int t_open         = params->T_open;

  return base_ptr 
          // + params->sig_size 
          // - sizeof(uint32_t);
          + ((params->tau - 1) * ell_bytes) //
          + c_mult_bytes          //
          + utilde_bytes //
          + ell_bytes //
          + (d_zk - 1) * lambda_bytes //
          + hcom_size * tau + (t_open * lambda_bytes) //
          + lambda_bytes //
          + IV_SIZE;  // 
          //
}

// helpers to compute position in signature (verify)

ATTR_PURE static inline const uint8_t* dsignature_c(const uint8_t* base_ptr, unsigned int index,
                                                    const faest_paramset_t* params) {
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;

  return base_ptr 
          + index * ell_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_c_mult(const uint8_t* base_ptr,
                                                    const faest_paramset_t* params) {
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;

  return base_ptr 
          + (params->tau - 1) * ell_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_u_tilde(const uint8_t* base_ptr,
                                                          const faest_paramset_t* params) {
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;

  return base_ptr 
          + (params->tau - 1) * ell_bytes
          + c_mult_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_d(const uint8_t* base_ptr,
                                                    const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes  = 4 * lambda_bytes;

  return base_ptr 
          + (params->tau - 1) * ell_bytes 
          + c_mult_bytes 
          + utilde_bytes;
}

// We do not need to store a0_tilde in the signature
ATTR_PURE static inline uint8_t* dsignature_a1toi_tilde(uint8_t* base_ptr,
                                                    const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;

  return base_ptr 
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes 
          + utilde_bytes 
          + ell_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_decom_i(const uint8_t* base_ptr,
                                                          const faest_paramset_t* params) {
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;
  const unsigned int d_zk           = params->d_zk;

  return base_ptr 
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes
          + utilde_bytes 
          + ell_bytes 
          + (d_zk - 1) * lambda_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_chall_3(const uint8_t* base_ptr,
                                                          const faest_paramset_t* params) {
  const faest_paramset_t* faest_param = faest_get_paramset(params->id);
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;
  const unsigned int d_zk           = params->d_zk;
  const unsigned int n_leafcom              = faest_is_em(faest_param) ? 2 : 3;
  const unsigned int hcom_size      = lambda_bytes * n_leafcom;
  const unsigned int tau            = params->tau;
  const unsigned int t_open         = params->T_open;

  return base_ptr 
          // + params->sig_size 
          // - sizeof(uint32_t) 
          // - IV_SIZE 
          // - lambda_bytes;
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes
          + utilde_bytes 
          + ell_bytes 
          + (d_zk - 1) * lambda_bytes
          + hcom_size * tau + t_open * lambda_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_iv_pre(const uint8_t* base_ptr,
                                                         const faest_paramset_t* params) {
  const faest_paramset_t* faest_param = faest_get_paramset(params->id);
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;
  const unsigned int d_zk           = params->d_zk;
  const unsigned int n_leafcom              = faest_is_em(faest_param) ? 2 : 3;
  const unsigned int hcom_size      = lambda_bytes * n_leafcom;
  const unsigned int tau            = params->tau;
  const unsigned int t_open         = params->T_open;

  return base_ptr 
          // + params->sig_size 
          // - sizeof(uint32_t) 
          // - IV_SIZE;
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes
          + utilde_bytes 
          + ell_bytes 
          + (d_zk - 1) * lambda_bytes
          + hcom_size * tau + t_open * lambda_bytes
          + lambda_bytes;
}

ATTR_PURE static inline const uint8_t* dsignature_ctr(const uint8_t* base_ptr,
                                                      const faest_paramset_t* params) {
  const faest_paramset_t* faest_param = faest_get_paramset(params->id);
  const unsigned int lambda_bytes  = params->lambda / 8;
  const unsigned int ell            = params->ell;
  const unsigned int ell_bytes      = (ell + 7) / 8;
  const unsigned n_mask                 = params->n_mask;
  const unsigned n_mult                 = params->n_mult;
  // const unsigned int n_mask_bytes   = (n_mask + 7) / 8;  // NOTE: I think here the c_mult is really stores continously
                                                            // instead of [n_mask_bytes] * n_mult
  const unsigned int c_mult_bytes   = (n_mult * n_mask + 7) / 8;
  const unsigned int utilde_bytes   = lambda_bytes * 4;
  const unsigned int d_zk           = params->d_zk;
  const unsigned int n_leafcom              = faest_is_em(faest_param) ? 2 : 3;
  const unsigned int hcom_size      = lambda_bytes * n_leafcom;
  const unsigned int tau            = params->tau;
  const unsigned int t_open         = params->T_open;

  return base_ptr 
          // + params->sig_size 
          // - sizeof(uint32_t);
          + ((params->tau - 1) * ell_bytes) 
          + c_mult_bytes
          + utilde_bytes 
          + ell_bytes 
          + (d_zk - 1) * lambda_bytes
          + hcom_size * tau + t_open * lambda_bytes
          + lambda_bytes
          + IV_SIZE;
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

ATTR_ALWAYS_INLINE ATTR_ARTIFICIAL static inline void hash_iv(uint8_t* iv, const uint8_t* iv_pre,
                                                              unsigned int lambda) {
  H4(iv, iv_pre, lambda);
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
                             const uint8_t* c, const uint8_t* c_mult_packed, const uint8_t* iv, unsigned int lambda,
                             unsigned int ell, unsigned int tau, unsigned int c_mult_size_bytes) {
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int ell_bytes     = (ell + 7) / 8;

  H2_context_t h2_ctx;
  H2_init(&h2_ctx, lambda);
  H2_update(&h2_ctx, mu, lambda_bytes * 2);
  H2_update(&h2_ctx, hcom, lambda_bytes * 2);
  H2_update(&h2_ctx, c, ell_bytes * (tau - 1));
  H2_update(&h2_ctx, c_mult_packed, c_mult_size_bytes);
  H2_update(&h2_ctx, iv, IV_SIZE);
  H2_1_final(&h2_ctx, chall_1, 8 * lambda_bytes + 8);
}

static void hash_challenge_2_init(H2_context_t* h2_ctx, const uint8_t* chall_1,
                                  const uint8_t* u_tilde, unsigned int lambda) {
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int u_tilde_bytes = lambda_bytes * 4;

  H2_init(h2_ctx, lambda);
  H2_update(h2_ctx, chall_1, 8 * lambda_bytes + 8);
  H2_update(h2_ctx, u_tilde, u_tilde_bytes);
}

static void hash_challenge_2_update_v_tilde(H2_context_t* h2_ctx, const uint8_t* v_tilde,
                                            unsigned int lambda) {
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int v_tilde_bytes = lambda_bytes * 4;

  H2_update(h2_ctx, v_tilde, v_tilde_bytes);
}

static void hash_challenge_2_finalize(uint8_t* chall_2, H2_context_t* h2_ctx, const uint8_t* d,
                                      const unsigned lambda, unsigned int ell) {
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ell_bytes    = (ell + 7) / 8;

  H2_update(h2_ctx, d, ell_bytes);
  H2_2_final(h2_ctx, chall_2, 3 * lambda_bytes + 8);
}

static void hash_challenge_3_init(H2_context_t* h2_ctx, const uint8_t* chall_2,
                                  const uint8_t* a0_tilde, const uint8_t* a1toi_tilde,
                                  unsigned int lambda, unsigned int d_zk) {
  const unsigned int lambda_bytes = lambda / 8;

  H2_init(h2_ctx, lambda);
  H2_update(h2_ctx, chall_2, 3 * lambda_bytes + 8);

  H2_update(h2_ctx, a0_tilde, lambda_bytes);

  for (unsigned int ai_tilde_idx = 0; ai_tilde_idx < d_zk - 1; ai_tilde_idx++) {
    H2_update(h2_ctx, a1toi_tilde + ai_tilde_idx * lambda_bytes, lambda_bytes);
  }
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
                              const uint8_t* a1toi_tilde, const uint8_t* ctr, unsigned int lambda, 
                              unsigned int d_zk) {
  H2_context_t h2_ctx;
  hash_challenge_3_init(&h2_ctx, chall_2, a0_tilde, a1toi_tilde, lambda, d_zk);
  H2_update(&h2_ctx, ctr, sizeof(uint32_t));
  H2_3_final(&h2_ctx, chall_3, lambda / 8);
}

static inline bool check_challenge_3(const uint8_t* chall_3, unsigned int start,
                                     unsigned int lambda) {
  const unsigned int lambda_bytes = lambda / 8;

  uint8_t res = (start % 8) ? (chall_3[start / 8] >> (start % 8)) : 0;
  for (unsigned int bytes = (start + 7) / 8; bytes != lambda_bytes; ++bytes) {
    res |= chall_3[bytes];
  }
  return !res;
}

static inline void free_pointer_array(uint8_t*** ptr) {
  free((*ptr)[0]);
  free(*ptr);
  *ptr = NULL;
}

// AES(-EM) OWF dispatchers
static inline void aes_prove(uint8_t* a0_tilde, uint8_t* a1toi_tilde, const uint8_t* w, uint8_t** V, 
                       const uint8_t* u_bar, const uint8_t* v_bar, const uint8_t* owf_in, 
                       const uint8_t* owf_out, const uint8_t* chall_2, const faest_paramset_t* params) {
  switch (params->lambda) {
  case 256:
    aes_256_prover(a0_tilde, a1toi_tilde, w, V, u_bar, v_bar, owf_in, owf_out, chall_2, params);
    break;
  case 192:
    aes_192_prover(a0_tilde, a1toi_tilde, w, V, u_bar, v_bar, owf_in, owf_out, chall_2, params);
    break;
  default:
    aes_128_prover(a0_tilde, a1toi_tilde, w, V, u_bar, v_bar, owf_in, owf_out, chall_2, params);
  }
}

static inline void aes_verify(uint8_t* a0_tilde, const uint8_t* d, const uint8_t** Q,
                              const uint8_t* q_bar, const uint8_t* owf_in, const uint8_t* owf_out,
                              const uint8_t* chall_2, const uint8_t* chall_3, const uint8_t* a1toi_tilde, const faest_paramset_t* params) {
  switch (params->lambda) {
  case 256:
    aes_256_verifier(a0_tilde, d, Q, q_bar, owf_in, owf_out, chall_2, chall_3, a1toi_tilde, params);
    break;
  case 192:
    aes_192_verifier(a0_tilde, d, Q, q_bar, owf_in, owf_out, chall_2, chall_3, a1toi_tilde, params);
    break;
  default:
    aes_128_verifier(a0_tilde, d, Q, q_bar, owf_in, owf_out, chall_2, chall_3, a1toi_tilde, params);
  }
}

static inline void compute_D(uint8_t* D, const uint8_t* Q_tilde, const uint8_t* delta, uint8_t* u_tilde, 
                              const faest_paramset_t* params) {
  unsigned int lambda_bytes = params->lambda / 8;
  switch (params->lambda) {
  case 256: {
    bf256_t bf_D[4];
    bf256_t bf_Q_tilde[4];
    bf256_t bf_delta;
    bf256_load(&bf_delta, delta);
    bf256_t bf_u_tilde[4];
    for (unsigned int i = 0; i < 4; i++) {
      bf256_load(&bf_Q_tilde[i], Q_tilde + i * lambda_bytes);
      bf256_load(&bf_u_tilde[i], u_tilde + i * lambda_bytes);

      bf256_mul_inplace(&bf_u_tilde[i], &bf_delta);

      bf256_add(&bf_D[i], &bf_Q_tilde[i], &bf_u_tilde[i]);

      bf256_store(D + i * lambda_bytes, &bf_D[i]);
    }
    break;
  }
  case 192: {
    bf192_t bf_D[4];
    bf192_t bf_Q_tilde[4];
    bf192_t bf_delta;
    bf192_load(&bf_delta, delta);
    bf192_t bf_u_tilde[4];
    for (unsigned int i = 0; i < 4; i++) {
      bf192_load(&bf_Q_tilde[i], Q_tilde + i * lambda_bytes);
      bf192_load(&bf_u_tilde[i], u_tilde + i * lambda_bytes);

      bf192_mul_inplace(&bf_u_tilde[i], &bf_delta);

      bf192_add(&bf_D[i], &bf_Q_tilde[i], &bf_u_tilde[i]);

      bf192_store(D + i * lambda_bytes, &bf_D[i]);
    }
    break;
  }
  default:
    bf128_t bf_D[4];
    bf128_t bf_Q_tilde[4];
    bf128_t bf_delta;
    bf128_load(&bf_delta, delta);
    bf128_t bf_u_tilde[4];
    for (unsigned int i = 0; i < 4; i++) {
      bf128_load(&bf_Q_tilde[i], Q_tilde + i * lambda_bytes);
      bf128_load(&bf_u_tilde[i], u_tilde + i * lambda_bytes);

      bf128_mul_inplace(&bf_u_tilde[i], &bf_delta);

      bf128_add(&bf_D[i], &bf_Q_tilde[i], &bf_u_tilde[i]);

      bf128_store(D + i * lambda_bytes, &bf_D[i]);
    }
  }
}

// FAEST.Sign()
void faest_sign(uint8_t* sig, const uint8_t* msg, size_t msg_len, const uint8_t* owf_key,
                const uint8_t* owf_input, const uint8_t* owf_output, const uint8_t* witness,
                const uint8_t* rho, size_t rholen, const faest_paramset_t* params) {
  const unsigned int ell           = params->ell;
  const unsigned int ell_bytes     = (ell + 7) / 8;
  const unsigned int lambda        = params->lambda;
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int d_zk          = params->d_zk;
  const unsigned int tau           = params->tau;
  const unsigned int k              = params->k;
  const unsigned int tau1            = params->tau1;
  const unsigned int d0            = bavc_max_node_depth(0, tau1, k);
  const unsigned int n_mult        = params->n_mult;
  const unsigned int n_mult_bytes        = (n_mult + 7) / 8;
  const unsigned int n_mask        = params->n_mask;
  const unsigned int n_mask_bytes  = (n_mask + 7) / 8;
  const unsigned int ell_hat       = ell + n_mask * d0;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
  const unsigned int w_grind       = params->w_grind;
  const unsigned int c_mult_bytes   = n_mask * n_mult_bytes;

  // line 2
  uint8_t mu[MAX_LAMBDA_BYTES * 2];
  hash_mu(mu, owf_input, params->owf_input_size, owf_output, params->owf_output_size, msg, msg_len,
          lambda);

  // line 3
  uint8_t rootkey[MAX_LAMBDA_BYTES], iv[IV_SIZE];
  hash_r_iv(rootkey, signature_iv_pre(sig, params), iv, owf_key, mu, rho, rholen, lambda);

  // line 6
  bavc_t bavc;
  uint8_t* u = malloc(ell_bytes);
  assert(u);
  // v has \hat \ell rows, \lambda columns, storing in column-major order
  uint8_t** V = malloc(lambda * sizeof(uint8_t*));                        // it is actually lambda - w_grind but keeping it lambda
  assert(V);
  V[0] = calloc(lambda, ell_hat_bytes);
  assert(V[0]);
  for (unsigned int i = 1; i < lambda; ++i) {
    V[i] = V[0] + i * ell_hat_bytes;
  }
  uint8_t* u_bar = malloc(n_mask * lambda_bytes);
  uint8_t* v_bar = malloc(n_mask * lambda_bytes);

  uint8_t* c_mult = malloc(c_mult_bytes);
  memset(c_mult, 0, c_mult_bytes);

  vole_commit(rootkey, iv, ell_hat, params, &bavc, signature_c(sig, 0, params), c_mult, u, V, u_bar, v_bar);

  unsigned int c_mult_idx = 0;
  uint8_t* c_mult_packed = signature_c_mult(sig, params);
  unsigned int c_mult_packed_bytes = (n_mask * n_mult + 7) / 8;
  memset(c_mult_packed, 0, c_mult_packed_bytes);
  pack_uint8(c_mult_packed, c_mult, n_mask, n_mult);

  // print_u8_array("c_mult_ptr", c_mult_ptr, 8);
  // print_u8_array("c_mult_packed", c_mult_packed, 8);

  free(c_mult);

  uint8_t** V_row = malloc(ell * sizeof(uint8_t*));                        // it is actually lambda - w_grind but keeping it lambda
  assert(V_row);
  V_row[0] = calloc(ell, lambda_bytes);
  assert(V_row[0]);
  for (unsigned int i = 1; i < ell; ++i) {
    V_row[i] = V_row[0] + i * lambda_bytes;
  }

  column_to_row_major(V, V_row, lambda, ell);

  
  H2_context_t h2_ctx;
  {
    // line 7
    uint8_t chall_1[(8 * MAX_LAMBDA_BYTES) + 8];
    memset(chall_1, 0, (8 * MAX_LAMBDA_BYTES) + 8);
    hash_challenge_1(chall_1, mu, bavc.h, signature_c(sig, 0, params), c_mult_packed, iv, lambda, ell, tau, c_mult_packed_bytes);

    // print_u8_array("sign chall_1", chall_1, (8 * lambda_bytes) + 8);

    // line 9
    unsigned int uh_0_bytes = lambda_bytes * (ell + d_zk - 1);
    uint8_t* uh_0 = malloc(uh_0_bytes);
    memset(uh_0, 0, uh_0_bytes);
    for (unsigned int bit_idx = 0; bit_idx < ell; bit_idx++) {
      ptr_set_bit(uh_0, bit_idx * lambda, ptr_get_bit(u, bit_idx));   // [u[0],..0 lambda times, u[1], 0 lambda times, ...]
    }

    for (unsigned int u_bar_idx = 0; u_bar_idx <= d_zk - 2; u_bar_idx++) {
      memcpy(uh_0 + (lambda_bytes * ell) + (u_bar_idx * lambda_bytes), u_bar + (u_bar_idx * lambda_bytes), lambda_bytes);
    }

    // line 10
    unsigned int uh_1_bytes = lambda_bytes * 4;
    uint8_t* uh_1 = malloc(uh_1_bytes);
    memset(uh_1, 0, uh_1_bytes);
    memcpy(uh_1, u_bar + ((d_zk - 1) * lambda_bytes), uh_1_bytes);

    // line 11
    unsigned int vh_0_bytes = lambda_bytes * (ell + d_zk - 1);
    uint8_t* vh_0 = malloc(vh_0_bytes);
    memset(vh_0, 0, vh_0_bytes);
    for (unsigned int v_idx = 0; v_idx < ell; v_idx++) {
      memcpy(vh_0 + v_idx * lambda_bytes, V_row[v_idx], lambda_bytes);
    }
    for (unsigned int v_idx = 0; v_idx <= d_zk - 2; v_idx++) {
      memcpy(vh_0 + (ell * lambda_bytes) + (v_idx * lambda_bytes), v_bar + v_idx * lambda_bytes, lambda_bytes);
    }
    // line 12
    unsigned int vh_1_bytes = lambda_bytes * 4;
    uint8_t* vh_1 = malloc(vh_1_bytes);
    memset(vh_1, 0, vh_1_bytes);
    memcpy(vh_1, v_bar + ((d_zk - 1) * lambda_bytes), vh_1_bytes);

    // line 13
    uint8_t* uh = malloc(uh_0_bytes + uh_1_bytes);
    memcpy(uh, uh_0, uh_0_bytes);
    memcpy(uh + uh_0_bytes, uh_1, uh_1_bytes);
    vole_hash_new(signature_u_tilde(sig, params), chall_1, uh, ell, d_zk, lambda);

    // To save memory consumption, the chall_2 is computed in an
    // Init-Update-Finalize style as V_tilde is only fed into to the hash and not
    // used elsewhere.
    // print_u8_array("sign signature_u_tilde", signature_u_tilde(sig, params), 4 * lambda_bytes);
    
    hash_challenge_2_init(&h2_ctx, chall_1, signature_u_tilde(sig, params), lambda);
    {
      // line 14
      uint8_t V_tilde[MAX_LAMBDA_BYTES * 4];
      memset(V_tilde, 0, MAX_LAMBDA_BYTES * 4);
      uint8_t* vh = malloc(vh_0_bytes + vh_1_bytes);
      memcpy(vh, vh_0, vh_0_bytes);
      memcpy(vh + vh_0_bytes, vh_1, vh_1_bytes);
      vole_hash_new(V_tilde, chall_1, vh, ell, d_zk, lambda);

      // print_u8_array("sign V_tilde", V_tilde, lambda_bytes * 4);

      // line 20
      hash_challenge_2_update_v_tilde(&h2_ctx, V_tilde, lambda);

      free(vh);
    }

    // free(c_mult_packed);

    free(V_row[0]);
    free(V_row);  

    free(uh);
    free(vh_0);
    free(vh_1);
    free(uh_0);
    free(uh_1);
  }

  // line 19
  xor_u8_array(witness, u, signature_d(sig, params), ell_bytes);

  {
    // print_u8_array("sign signature_d", signature_d(sig, params), ell_bytes);

    // line 20
    uint8_t chall_2[(3 * MAX_LAMBDA_BYTES) + 8];
    hash_challenge_2_finalize(chall_2, &h2_ctx, signature_d(sig, params), lambda, ell);

    // print_u8_array("sign chall_2", chall_2, 3 * lambda_bytes + 8);

    // line 22
    uint8_t* a0_tilde = malloc(lambda_bytes);
    aes_prove(a0_tilde, signature_a1toi_tilde(sig, params), witness, V, u_bar, v_bar, owf_input, owf_output, chall_2, params);

    free_pointer_array(&V);
    free(u);
    u = NULL;

    // print_u8_array("sign a0_tilde", a0_tilde, lambda_bytes);
    // print_u8_array("sign a1_tilde", signature_a1toi_tilde(sig, params), lambda_bytes);
    // print_u8_array("sign a2_tilde", signature_a1toi_tilde(sig, params) + lambda_bytes, lambda_bytes);

    // line 26
    hash_challenge_3_init(&h2_ctx, chall_2, a0_tilde, signature_a1toi_tilde(sig, params), lambda, d_zk);

    free(a0_tilde);
  }

  uint32_t ctr = 0;
  for (; true; ++ctr) {
    uint8_t* chall_3 = signature_chall_3(sig, params);
    hash_challenge_3_final(chall_3, &h2_ctx, ctr, lambda);
    
    
    // declassify chall_3 which is put into the signature
    faest_declassify(chall_3, lambda / 8);

    // line 27
    if (!check_challenge_3(chall_3, lambda - w_grind, lambda)) {
      continue;
    }

    // line 29
    uint16_t decoded_chall_3[MAX_TAU];
    if (!decode_all_chall_3(decoded_chall_3, chall_3, params)) {
      continue;
    }

    // line 30
    if (bavc_open(signature_decom_i(sig, params), &bavc, decoded_chall_3, params)) {
      break;
    }
  }

  // size_t fwd_end = (signature_decom_i(sig,params) - sig)
  //         + (lambda_bytes * 2)*tau + params->T_open*lambda_bytes;
  // size_t bwd_start = params->sig_size - sizeof(uint32_t) - IV_SIZE - lambda_bytes;
  // fprintf(stderr, "decom_i end=%zu, chall_3 start=%zu, sig_size=%u\n",
  //         fwd_end, bwd_start, params->sig_size);

  hash_clear(&h2_ctx);
  bavc_clear(&bavc);

  // print_u8_array("sign chall_3", signature_chall_3(sig, params), lambda_bytes);

  // copy counter to signature
  ctr = htole32(ctr);
  memcpy(signature_ctr(sig, params), &ctr, sizeof(ctr));

  free(u_bar);
  free(v_bar);

  const unsigned int n_leafcom              = faest_is_em(params) ? 2 : 3;
  const unsigned int t_open         = params->T_open;
  const unsigned int hcom_size      = lambda_bytes * n_leafcom;

  // // NOTE: Claude debug code ----
  // const size_t needed = 
  //       ((params->tau - 1) * ell_bytes)
  //       + c_mult_bytes
  //       + lambda_bytes * 4
  //       + ell_bytes
  //       + (d_zk - 1) * lambda_bytes
  //       + hcom_size * tau + (t_open * lambda_bytes)
  //       + lambda_bytes
  //       + IV_SIZE
  //       + sizeof(uint32_t);
  // printf("\nsig_size=%u needed=%zu diff=%td\n",
  //      params->sig_size, needed, (ptrdiff_t)params->sig_size - (ptrdiff_t)needed);

  // const size_t decom_i_end = (size_t)(signature_decom_i(sig, params) - sig)
  //                        + 2u * tau * lambda_bytes + params->T_open * lambda_bytes;
  // const size_t tail_start  = (size_t)(signature_chall_3(sig, params) - sig);
  // assert(decom_i_end <= tail_start);

  // printf("hcom_size=%u  2*lambda_bytes=%u  tau=%u t_open=%u\n",
  //      hcom_size, 2*lambda_bytes, params->tau, t_open);
  // // ----

  

  // assert(sig - (((params->tau - 1) * ell_bytes)
  //         + c_mult_bytes
  //         + lambda_bytes * 4
  //         + ell_bytes
  //         + (d_zk - 1) * lambda_bytes
  //         + hcom_size * tau + (t_open * lambda_bytes)
  //         + lambda_bytes
  //         + IV_SIZE
  //         + sizeof(uint32_t)) == 0);

}

int faest_verify(const uint8_t* msg, size_t msglen, const uint8_t* sig, const uint8_t* owf_input,
                 const uint8_t* owf_output, const faest_paramset_t* params) {
  const unsigned int ell           = params->ell;
  const unsigned int ell_bytes           = (ell + 7) / 8;
  const unsigned int lambda        = params->lambda;
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int d_zk          = params->d_zk;
  const unsigned int tau           = params->tau;
  const unsigned int k              = params->k;
  const unsigned int tau1           = params->tau1;
  const unsigned int d0            = bavc_max_node_depth(0, tau1, k);
  const unsigned int n_mult        = params->n_mult;
  const unsigned int n_mult_bytes        = (n_mult + 7) / 8;
  const unsigned int n_mask        = params->n_mask;
  const unsigned int n_mask_bytes  = (n_mask + 7) / 8;
  const unsigned int ell_hat       = ell + n_mask * d0;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
  const unsigned int utilde_bytes  = lambda_bytes * 4;
  const unsigned int c_mult_bytes   = n_mask * n_mult_bytes;
  const unsigned int w_grind      = params->w_grind;
  const unsigned int lambda_minus_w_grind = lambda - w_grind;
  const unsigned int lambda_minus_w_grind_bytes = ((lambda_minus_w_grind) + 7 ) / 8;

  // line 3
  if (!check_challenge_3(dsignature_chall_3(sig, params), lambda - params->w_grind, lambda)) {
    return -1;
  }

  // line 1
  uint8_t mu[MAX_LAMBDA_BYTES * 2];
  hash_mu(mu, owf_input, params->owf_input_size, owf_output, params->owf_output_size, msg, msglen,
          lambda);

  // line 2
  uint8_t iv[IV_SIZE];
  hash_iv(iv, dsignature_iv_pre(sig, params), lambda);

  uint8_t* c_mult = malloc(c_mult_bytes);
  unsigned int c_mult_packed_bytes = (n_mask * n_mult + 7) / 8;
  uint8_t* c_mult_packed = dsignature_c_mult(sig, params);
  unpack_uint8(c_mult, c_mult_packed, n_mask, n_mult);

  // unsigned int c_mult_idx = 0;
  // const uint8_t* c_mult_ptr = dsignature_c_mult(sig, params);
  // for (unsigned int n_mask_idx = 0; n_mask_idx < n_mask; n_mask_idx++) {
  //   for (unsigned int n_mult_idx = 0; n_mult_idx < n_mult; n_mult_idx++) {
  //     ptr_set_bit(c_mult_packed + n_mask_idx * n_mult_bytes, 
  //                 n_mult_idx, 
  //                 ptr_get_bit(c_mult_ptr, c_mult_idx));

  //     c_mult_idx++;
  //   }
  // }

  // print_u8_array("c_mult_ptr", c_mult_ptr, 8);
  // print_u8_array("c_mult_as_n_mask_bytes", c_mult_as_n_mask_bytes, 8);

  // line 6
  // q is a \hat \ell \times \lambda matrix
  uint8_t** Q = malloc(lambda * sizeof(uint8_t*));
  assert(Q);
  Q[0] = calloc(lambda, ell_hat_bytes);
  assert(Q[0]);
  for (unsigned int i = 1; i < lambda; ++i) {
    Q[i] = Q[0] + i * ell_hat_bytes;
  }
  uint8_t hcom[MAX_LAMBDA_BYTES * 2];
  uint8_t* q_bar = malloc(n_mask * lambda_bytes);
  uint8_t* Delta = malloc(lambda_bytes);
  memset(Delta, 0, lambda_bytes);
  if (!vole_reconstruct(hcom, Q, iv, dsignature_chall_3(sig, params),
                        dsignature_decom_i(sig, params), dsignature_c(sig, 0, params), c_mult, 
                        q_bar, Delta, ell_hat, params)) {
    free_pointer_array(&Q);
    return -1;
  }

  free(c_mult);

  uint8_t** Q_row = malloc(ell * sizeof(uint8_t*));                        // it is actually lambda - w_grind but keeping it lambda
  assert(Q_row);
  Q_row[0] = calloc(ell, lambda_bytes);
  assert(Q_row[0]);
  for (unsigned int i = 1; i < ell; ++i) {
    Q_row[i] = Q_row[0] + i * lambda_bytes;
  }

  column_to_row_major(Q, Q_row, lambda, ell);

  // line 9
  uint8_t chall_1[(8 * MAX_LAMBDA_BYTES) + 8];
  memset(chall_1, 0, (8 * MAX_LAMBDA_BYTES) + 8);
  hash_challenge_1(chall_1, mu, hcom, dsignature_c(sig, 0, params), c_mult_packed, iv, lambda, ell, tau, c_mult_packed_bytes);

  // print_u8_array("verify chall_1", chall_1, 8 * lambda_bytes + 8);

  // line 10
  // print_u8_array("verify dsignature_u_tilde", dsignature_u_tilde(sig, params), 4 * lambda_bytes);
  H2_context_t chall_2_ctx;
  hash_challenge_2_init(&chall_2_ctx, chall_1, dsignature_u_tilde(sig, params), lambda);
  {

    unsigned int qh_0_bytes = lambda_bytes * (ell + d_zk - 1);
    uint8_t* qh_0 = malloc(qh_0_bytes);
    memset(qh_0, 0, qh_0_bytes);
    for (unsigned int q_idx = 0; q_idx < ell; q_idx++) {
      memcpy(qh_0 + q_idx * lambda_bytes, Q_row[q_idx], lambda_bytes);
    }
    for (unsigned int v_idx = 0; v_idx <= d_zk - 2; v_idx++) {
      memcpy(qh_0 + (ell * lambda_bytes) + (v_idx * lambda_bytes), q_bar + v_idx * lambda_bytes, lambda_bytes);
    }
    unsigned int qh_1_bytes = lambda_bytes * 4;
    uint8_t* qh_1 = malloc(qh_1_bytes);
    memset(qh_1, 0, qh_1_bytes);
    memcpy(qh_1, q_bar + ((d_zk - 1) * lambda_bytes), qh_1_bytes);
    

    // line 13
    uint8_t Q_tilde[MAX_LAMBDA_BYTES * 4];
    memset(Q_tilde, 0, MAX_LAMBDA_BYTES * 4);
    uint8_t* qh = malloc(qh_0_bytes + qh_1_bytes);
    memset(qh, 0, qh_0_bytes + qh_1_bytes);
    memcpy(qh, qh_0, qh_0_bytes);
    memcpy(qh + qh_0_bytes, qh_1, qh_1_bytes);

    vole_hash_new(Q_tilde, chall_1, qh, ell, d_zk, lambda);

    // print_u8_array("delta_prime", Delta, 16);

    compute_D(Q_tilde, Q_tilde, Delta, dsignature_u_tilde(sig, params), params);

    // print_u8_array("verify Q_tilde", Q_tilde, lambda_bytes * 4);

    // line 15
    hash_challenge_2_update_v_tilde(&chall_2_ctx, Q_tilde, lambda);

    free(qh);
    free(qh_0);
    free(qh_1);
  }

  // print_u8_array("verify dsignature_d", dsignature_d(sig, params), ell_bytes);

  // line 16
  uint8_t chall_2[3 * MAX_LAMBDA_BYTES + 8];
  hash_challenge_2_finalize(chall_2, &chall_2_ctx, dsignature_d(sig, params), lambda, ell);

  // print_u8_array("verify chall_2", chall_2, 3 * lambda_bytes + 8);

  // line 18
  uint8_t a0_tilde[MAX_LAMBDA_BYTES];
  aes_verify(a0_tilde, dsignature_d(sig, params), (const uint8_t**)Q, q_bar, owf_input, owf_output, chall_2, Delta, dsignature_a1toi_tilde(sig, params), params);
  free_pointer_array(&Q);

  // print_u8_array("verify a0_tilde", a0_tilde, lambda_bytes);
  // print_u8_array("verify a1_tilde", dsignature_a1toi_tilde(sig, params), lambda_bytes);
  // print_u8_array("verify a2_tilde", dsignature_a1toi_tilde(sig, params) + lambda_bytes, lambda_bytes);

  // line 19
  uint8_t chall_3_prime[MAX_LAMBDA_BYTES];
  hash_challenge_3(chall_3_prime, chall_2, a0_tilde, dsignature_a1toi_tilde(sig, params), dsignature_ctr(sig, params), lambda, d_zk);

  // print_u8_array("verify chall_3_prime", chall_3_prime, lambda_bytes);

  free(Q_row[0]);
  free(Q_row); 
  free(q_bar);
  // free(c_mult_packed);
  free(Delta);

  // Step 21
  return memcmp(chall_3_prime, dsignature_chall_3(sig, params), lambda_bytes) == 0 ? 0 : -1;
}
