#include <stddef.h>
#include <stdio.h>

#include "macros.h"
#include "vbb.h"
#include "vole.h"
#include "vole_stream.h"
#include "vc_stream.h"
#include "vc.h"
#include "instances.h"
#include "faest.h"
#include "faest_aes.h"
#include "fields.h"
#include "parameters.h"
#include "utils.h"

static void recompute_hash(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int ellhat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;

  unsigned int amount       = MIN(len, vbb->params->faest_param.lambda - start);
  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  partial_vole_commit_cmo(vbb->root_key, vbb->iv, ellhat, vbb->params, sVecCom, vbb->vole_V_cache,
                          start, amount, NULL, NULL, NULL);
  free(sVecCom);

  vbb->cache_idx = start;
}

static void recompute_prove(vbb_t* vbb, unsigned int start, unsigned int len) {
  if (len >= vbb->params->faest_param.l + vbb->params->faest_param.lambda) {
    start = 0;
  } else if (start + len > vbb->params->faest_param.l + vbb->params->faest_param.lambda) {
    start = vbb->params->faest_param.l + vbb->params->faest_param.lambda - len;
  }

  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  partial_vole_commit_rmo(vbb->root_key, vbb->iv, start, len, vbb->params, sVecCom,
                          vbb->vole_V_cache);
  free(sVecCom);

  vbb->cache_idx = start;
}

// len is the number of OLE v's that is allowed to be stored in memory.
// Hence we store (at most) len*lambda in memory.
void init_vbb(vbb_t* vbb, unsigned int len, const uint8_t* root_key, const uint8_t* iv, uint8_t* c,
              const faest_paramset_t* params) {
  vbb->iv       = iv;
  vbb->com_hash = calloc(MAX_LAMBDA_BYTES * 2, sizeof(uint8_t));
  vbb->params   = params;
  vbb->root_key = root_key;

  unsigned int lambda       = vbb->params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  const unsigned int ellhat = vbb->params->faest_param.l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  vbb->full_size            = len >= ellhat;
  vbb->vole_U               = malloc(ellhat_bytes);

  // PROVE cache
  // FIXME - would MAX(len, vbb->params->faest_param.Lenc) be correct for all variants?
  // FIXME - do not allocate too much memory...
  unsigned int row_count = len;
  vbb->row_count         = row_count;
  vbb->vole_V_cache      = calloc(row_count, lambda_bytes);

  // Setup u hcom c.
  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  if (vbb->full_size) {
    vbb->v_buf = malloc(lambda_bytes);
    partial_vole_commit_cmo(vbb->root_key, vbb->iv, ellhat, vbb->params, sVecCom, vbb->vole_V_cache,
                            0, ellhat, vbb->vole_U, vbb->com_hash, c);
  } else {
    vole_commit_u_hcom_c(vbb->root_key, vbb->iv, ellhat, vbb->params, vbb->com_hash, sVecCom, c,
                         vbb->vole_U);
  }
  free(sVecCom);
  // HASH cache
  unsigned int column_count = (size_t)vbb->row_count * (size_t)lambda_bytes / (size_t)ellhat_bytes;
  assert(column_count >= 1);
  vbb->column_count = column_count;
}

void clean_vbb(vbb_t* vbb) {
  free(vbb->com_hash);
  free(vbb->vole_U);
  free(vbb->vole_V_cache);
  if (vbb->full_size) {
    free(vbb->v_buf);
  }
}

void prepare_hash(vbb_t* vbb) {
  vbb->cache_idx = 0;
  if (vbb->full_size) {
    return;
  }
  recompute_hash(vbb, 0, vbb->column_count);
}

void prepare_prove(vbb_t* vbb) {
  vbb->cache_idx = 0;
  if (vbb->full_size) {
    return;
  }
  unsigned int len = vbb->row_count;
  recompute_prove(vbb, 0, len);
}

inline uint8_t* get_vole_v_hash(vbb_t* vbb, unsigned int idx) {
  const unsigned int ellhat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;

  assert(idx < vbb->params->faest_param.lambda);
  if (!(idx >= vbb->cache_idx && idx < vbb->cache_idx + vbb->column_count)) {
    recompute_hash(vbb, idx, vbb->column_count);
  }

  // FIXME - should it be ell + lambda instead? (also change offset in partial_vole_commit_cmo
  // then!)
  unsigned int offset = idx - vbb->cache_idx;
  return vbb->vole_V_cache + offset * ellhat_bytes;
}

static inline uint8_t* get_vole_v_prove(vbb_t* vbb, unsigned int idx) {
  unsigned int lambda       = vbb->params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat       = vbb->params->faest_param.l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;

  // FIXME: verifier 
  if (vbb->full_size) {
    memset(vbb->v_buf, 0, lambda_bytes);
    // Transpose on the fly into v_buf
    for (unsigned int column = 0; column != lambda; ++column) {
      ptr_set_bit(vbb->v_buf, ptr_get_bit(vbb->vole_V_cache + column * ellhat_bytes, idx), column);
    }
    return vbb->v_buf;
  }

  if (!(idx >= vbb->cache_idx && idx < vbb->cache_idx + vbb->row_count)) {
    recompute_prove(vbb, idx, vbb->row_count);
  }
  unsigned int offset = (idx - vbb->cache_idx) * (vbb->params->faest_param.lambda / 8);
  return vbb->vole_V_cache + offset;
}

bf256_t* get_vole_v_prove_256(vbb_t* vbb, unsigned int idx) {
  return (bf256_t*)get_vole_v_prove(vbb, idx);
}

bf192_t* get_vole_v_prove_192(vbb_t* vbb, unsigned int idx) {
  return (bf192_t*)get_vole_v_prove(vbb, idx);
}

bf128_t* get_vole_v_prove_128(vbb_t* vbb, unsigned int idx) {
  return (bf128_t*)get_vole_v_prove(vbb, idx);
}

uint8_t* get_vole_u(vbb_t* vbb) {
  return vbb->vole_U;
}

uint8_t* get_com_hash(vbb_t* vbb) {
  return vbb->com_hash;
}

// TODO - refactor this stuff
void vector_open_ondemand(vbb_t* vbb, unsigned int idx, const uint8_t* s_, uint8_t* sig_pdec,
                          uint8_t* sig_com, unsigned int depth) {
  unsigned int lambda       = vbb->params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int tau          = vbb->params->faest_param.tau;
  uint8_t* expanded_keys    = malloc(tau * lambda_bytes);
  prg(vbb->root_key, vbb->iv, expanded_keys, lambda, lambda_bytes * tau);

  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  stream_vector_commitment(expanded_keys + lambda_bytes * idx, lambda, &sVecCom[idx], depth);
  stream_vector_open(&sVecCom[idx], s_, sig_pdec, sig_com, depth, vbb->iv, lambda);
  free(sVecCom);
  free(expanded_keys);
}

// QBB Implementation
void init_qbb(vbb_t* qbb, unsigned int len, const uint8_t* iv, uint8_t* c, uint8_t* pdec_sig,
              uint8_t* com_sig, uint8_t* chall3, uint8_t* u_tilde, const faest_paramset_t* params,
              const uint8_t* sig) {
  qbb->iv        = iv;
  qbb->row_count = len;
  qbb->params    = params;
  qbb->iv        = iv;
  qbb->c         = c;
  qbb->com_hash  = calloc(MAX_LAMBDA_BYTES * 2, sizeof(uint8_t));
  qbb->full_size = false;

  const unsigned int lambda        = params->faest_param.lambda;
  const unsigned int l             = params->faest_param.l;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int tau           = params->faest_param.tau;
  const unsigned int tau0          = params->faest_param.t0;
  const size_t lambda_bytes        = params->faest_param.lambda / 8;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int k0            = params->faest_param.k0;
  const unsigned int k1            = params->faest_param.k1;

  uint8_t** qprime = malloc(lambda * sizeof(uint8_t*));
  qprime[0]        = calloc(lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    qprime[i] = qprime[0] + i * ell_hat_bytes;
  }
  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  for (unsigned int i = 0; i < tau; ++i) {
    pdec[i] = dsignature_pdec(sig, i, params);
    com[i]  = dsignature_com(sig, i, params);
  }

  vole_reconstruct(iv, chall3, pdec, com, qbb->com_hash, qprime, ell_hat, params);

  uint8_t** q = malloc(lambda * sizeof(uint8_t*));
  q[0]        = calloc(lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    q[i] = q[0] + i * ell_hat_bytes;
  }

  qbb->Dtilde    = malloc(lambda * sizeof(uint8_t*));
  qbb->Dtilde[0] = calloc(lambda, (lambda_bytes + UNIVERSAL_HASH_B));
  for (unsigned int i = 1; i < lambda; ++i) {
    qbb->Dtilde[i] = qbb->Dtilde[0] + i * (lambda_bytes + UNIVERSAL_HASH_B);
  }

  unsigned int Dtilde_idx = 0;
  unsigned int q_idx      = 0;
  for (unsigned int i = 0; i < tau; i++) {
    const unsigned int depth = i < tau0 ? k0 : k1;

    // Step 11
    uint8_t delta[MAX_DEPTH];
    ChalDec(chall3, i, params->faest_param.k0, params->faest_param.t0, params->faest_param.k1,
            params->faest_param.t1, delta);
    // Step 16
    for (unsigned int j = 0; j != depth; ++j, ++Dtilde_idx) {
      // for scan-build
      assert(Dtilde_idx < lambda);
      masked_xor_u8_array(qbb->Dtilde[Dtilde_idx], u_tilde, qbb->Dtilde[Dtilde_idx], delta[j],
                          utilde_bytes);
    }

    if (i == 0) {
      // Step 8
      memcpy(q[q_idx], qprime[q_idx], ell_hat_bytes * depth);
      q_idx += depth;
    } else {
      // Step 14
      for (unsigned int d = 0; d < depth; ++d, ++q_idx) {
        masked_xor_u8_array(qprime[q_idx], c + (i - 1) * ell_hat_bytes, q[q_idx], delta[d],
                            ell_hat_bytes);
      }
    }
  }
  free(qprime[0]);
  free(qprime);
  qprime                  = NULL;
  qbb->vole_Q_cache       = q[0];
  qbb->vole_Q_cache_index = q;
}

uint8_t* get_vole_q_hash(vbb_t* qbb, unsigned int idx) {
  qbb->cache_idx = 0;
  const unsigned int ellhat =
      qbb->params->faest_param.l + qbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int offset       = idx - qbb->cache_idx;
  return qbb->vole_Q_cache + offset * ellhat_bytes;
}

void prepare_verify_qbb(vbb_t* qbb, const uint8_t* sig_d, const uint8_t* sig_chall_3) {
  qbb->cache_idx                   = 0;
  const unsigned int lambda        = qbb->params->faest_param.lambda;
  const unsigned int l             = qbb->params->faest_param.l;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int tau = qbb->params->faest_param.tau;
  const unsigned int t0  = qbb->params->faest_param.t0;
  const unsigned int k0  = qbb->params->faest_param.k0;
  const unsigned int t1  = qbb->params->faest_param.t1;
  const unsigned int k1  = qbb->params->faest_param.k1;

  // TODO: Actually EM use Lenc, but Lenc == L for all EM..
  unsigned int size = qbb->params->faest_param.l;

  for (unsigned int i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t decoded_challenge[MAX_DEPTH];
    ChalDec(sig_chall_3, i, k0, t0, k1, t1, decoded_challenge);
    for (unsigned int j = 0; j < depth; j++, ++col) {
      if (decoded_challenge[j] == 1) {
        xor_u8_array(sig_d, qbb->vole_Q_cache_index[col], qbb->vole_Q_cache_index[col], (size + 7) / 8);
      }
    }
  }

  switch (lambda) {
  case 256:
    qbb->vole_Q_cache_RMO =
        (uint8_t*)column_to_row_major_and_shrink_V_256(qbb->vole_Q_cache_index, FAEST_256F_L);
    break;
  case 192:
    qbb->vole_Q_cache_RMO =
        (uint8_t*)column_to_row_major_and_shrink_V_192(qbb->vole_Q_cache_index, FAEST_192F_L);
    break;
  default:
    qbb->vole_Q_cache_RMO =
        (uint8_t*)column_to_row_major_and_shrink_V_128(qbb->vole_Q_cache_index, FAEST_128F_L);
    break;
  }
  qbb->vole_V_cache = qbb->vole_Q_cache_RMO;
  //memcpy(qbb->vole_V_cache, qbb->vole_Q_cache_RMO, (l + lambda) * lambda / 8);
}

uint8_t* get_vole_q_verify(vbb_t* qbb, unsigned int idx) {
  unsigned int offset = (idx - qbb->cache_idx) * (qbb->params->faest_param.lambda / 8);
  return qbb->vole_Q_cache_RMO + offset;
}

bf128_t* get_vole_q_verify_128(vbb_t* qbb, unsigned int idx) {
  return (bf128_t*)get_vole_q_verify(qbb, idx);
}

bf192_t* get_vole_q_verify_192(vbb_t* qbb, unsigned int idx) {
  return (bf192_t*)get_vole_q_verify(qbb, idx);
}

bf256_t* get_vole_q_verify_256(vbb_t* qbb, unsigned int idx) {
  return (bf256_t*)get_vole_q_verify(qbb, idx);
}
