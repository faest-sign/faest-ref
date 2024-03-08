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

static void recompute_hash_prove(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int ellhat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;

  unsigned int amount       = MIN(len, vbb->params->faest_param.lambda - start);
  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  partial_vole_commit_cmo(vbb->root_key, vbb->iv, ellhat, vbb->params, sVecCom, vbb->vole_cache,
                          start, amount, NULL, NULL, NULL);
  free(sVecCom);

  vbb->cache_idx = start;
}

static void recompute_aes_prove(vbb_t* vbb, unsigned int start, unsigned int len) {
  unsigned int lambda = vbb->params->faest_param.lambda;
  unsigned int ell    = vbb->params->faest_param.l;
  unsigned int tau    = vbb->params->faest_param.tau;

  if (len >= ell + lambda) {
    start = 0;
  } else if (start + len > ell + lambda) {
    start = ell + lambda - len;
  }

  stream_vec_com_t* sVecCom = calloc(tau, sizeof(stream_vec_com_t));
  partial_vole_commit_rmo(vbb->root_key, vbb->iv, start, len, vbb->params, sVecCom,
                          vbb->vole_cache);
  free(sVecCom);

  vbb->cache_idx = start;
}

// len is the number of OLE v's that is allowed to be stored in memory.
// Hence we store (at most) len*lambda in memory.
void init_vbb_prove(vbb_t* vbb, unsigned int len, const uint8_t* root_key, const uint8_t* iv,
                    uint8_t* c, const faest_paramset_t* params) {
  vbb->party    = SIGNER;
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

  unsigned int row_count = len;
  vbb->row_count         = row_count;
  vbb->vole_cache        = calloc(row_count, lambda_bytes);

  // Setup u, hcom, and c.
  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  if (vbb->full_size) {
    vbb->v_buf = malloc(lambda_bytes);
    partial_vole_commit_cmo(vbb->root_key, vbb->iv, ellhat, vbb->params, sVecCom, vbb->vole_cache,
                            0, lambda, vbb->vole_U, vbb->com_hash, c);
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
  free(vbb->vole_cache);
  if (vbb->full_size) {
    free(vbb->v_buf);
  }
}

void prepare_hash_prove(vbb_t* vbb) {
  if (vbb->full_size) {
    return;
  }
  recompute_hash_prove(vbb, 0, vbb->column_count);
}

void prepare_aes_prove(vbb_t* vbb) {
  if (vbb->full_size) {
    return;
  }
  unsigned int len = vbb->row_count;
  recompute_aes_prove(vbb, 0, len);
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

// Verifier implementation
static void recompute_hash_verify(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda        = vbb->params->faest_param.lambda;
  const unsigned int l             = vbb->params->faest_param.l;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
  const unsigned int tau           = vbb->params->faest_param.tau;
  const unsigned int tau0          = vbb->params->faest_param.t0;
  const unsigned int k0            = vbb->params->faest_param.k0;
  const unsigned int k1            = vbb->params->faest_param.k1;
  const uint8_t* chall3            = dsignature_chall_3(vbb->sig, vbb->params);

  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  for (unsigned int i = 0; i < tau; ++i) {
    pdec[i] = dsignature_pdec(vbb->sig, i, vbb->params);
    com[i]  = dsignature_com(vbb->sig, i, vbb->params);
  }

  unsigned int amount = MIN(len, lambda - start);

  partial_vole_reconstruct_cmo(vbb->iv, chall3, pdec, com, NULL, vbb->vole_cache, ell_hat,
                               vbb->params, start, amount);

  // Apply correction values
  const uint8_t* c     = dsignature_c(vbb->sig, 0, vbb->params);
  unsigned int col_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    const unsigned int depth = i < tau0 ? k0 : k1;
    if (start >= col_idx + depth) {
      col_idx += depth;
      continue;
    }
    if (i == 0) { // FIXME remove this first iteration
      col_idx += depth;
      continue;
    }

    uint8_t delta[MAX_DEPTH];
    ChalDec(chall3, i, vbb->params->faest_param.k0, vbb->params->faest_param.t0,
            vbb->params->faest_param.k1, vbb->params->faest_param.t1, delta);

    
    // Construct q inplace of qprime
    for (unsigned int d = 0; d < depth; d++, col_idx++) {
      if (start > col_idx) {
        continue;
      }
      masked_xor_u8_array(
          vbb->vole_cache + (col_idx - start) * ell_hat_bytes, c + (i - 1) * ell_hat_bytes,
          vbb->vole_cache + (col_idx - start) * ell_hat_bytes, delta[d], ell_hat_bytes);

      if (col_idx + 1 >= start + amount) {
        col_idx++;
        break;
      }
    }
    if (col_idx >= start + amount) {
      break;
    }
  }
  vbb->cache_idx = start;
}

static void recompute_aes_verify(vbb_t* vbb, unsigned int start, unsigned int len) {
  unsigned int lambda        = vbb->params->faest_param.lambda;
  unsigned int lambda_bytes  = lambda / 8;
  unsigned int ell           = vbb->params->faest_param.l;
  unsigned int tau           = vbb->params->faest_param.tau;
  unsigned int ell_hat       = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ell_hat_bytes = ell_hat / 8;
  unsigned int tau0          = vbb->params->faest_param.t0;
  unsigned int tau1          = vbb->params->faest_param.t1;
  unsigned int k0            = vbb->params->faest_param.k0;
  unsigned int k1            = vbb->params->faest_param.k1;

  if (len >= ell + lambda) {
    start = 0;
  } else if (start + len > ell + lambda) {
    start = ell + lambda - len;
  }

  const uint8_t* chall3 = dsignature_chall_3(vbb->sig, vbb->params);
  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  for (unsigned int i = 0; i < tau; ++i) {
    pdec[i] = dsignature_pdec(vbb->sig, i, vbb->params);
    com[i]  = dsignature_com(vbb->sig, i, vbb->params);
  }
  partial_vole_reconstruct_rmo(vbb->iv, chall3, pdec, com, vbb->vole_cache, ell_hat, vbb->params,
                               start, len);

  // Apply transposed correction values
  uint8_t buf[MAX_LAMBDA_BYTES] = {0};
  const uint8_t* c              = dsignature_c(vbb->sig, 0, vbb->params);
  for (unsigned int row_idx = 0; row_idx < len; row_idx++) {
    memset(buf, 0, sizeof(buf));
    // uint8_t* buf_prt = buf + 1;
    uint8_t byte             = 0;
    unsigned int bit_counter = k0;
    // Pack all c_i into buf lambda
    for (unsigned int t = 1; t < tau; t++) {
      unsigned int depth = t < tau0 ? k0 : k1;

      const uint8_t* c_idx = c + (t - 1) * ell_hat_bytes;
      unsigned int c_byte  = (row_idx + start) / 8;
      unsigned int c_bit   = (row_idx + start) % 8;
      uint8_t bit          = c_idx[c_byte] >> c_bit & 1;
      // pack bit into buf
      for (unsigned int i = 0; i < depth; i++) {
        byte = byte >> 1;
        byte |= (bit << 7);
        bit_counter++;
        if (bit_counter % 8 == 0) {
          buf[bit_counter / 8 - 1] = byte;
          byte                     = 0;
        }
      }
    }

    // AND buf with Delta
    for (unsigned int i = 0; i < lambda_bytes; i++) {
      buf[i] &= chall3[i];
      vbb->vole_cache[row_idx * lambda_bytes + i] ^= buf[i];
    }
  }

  // Apply transposed witness values
  const uint8_t* d          = dsignature_d(vbb->sig, vbb->params);
  unsigned int full_col_idx = 0;
  unsigned int end_row_idx  = MIN(len, ell - start);
  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;
    uint8_t delta[MAX_DEPTH];
    ChalDec(dsignature_chall_3(vbb->sig, vbb->params), i, k0, tau0, k1, tau1, delta);
    for (unsigned int col_idx = 0; col_idx < depth; col_idx++) {
      uint8_t delta_i = delta[col_idx];
      if (delta_i == 0) {
        continue;
      }
      for (unsigned int row_idx = 0; row_idx < end_row_idx; row_idx++) {
        if (row_idx + start > ell) {
          break;
        }
        unsigned int d_byte = (row_idx + start) / 8;
        unsigned int d_bit  = (row_idx + start) % 8;
        uint8_t bit         = d[d_byte] >> d_bit & 1;
        if (bit == 0) {
          continue;
        }
        unsigned int q_byte = (full_col_idx + col_idx) / 8;

        vbb->vole_cache[row_idx * lambda_bytes + q_byte] ^= bit << (full_col_idx + col_idx) % 8;
      }
    }
    full_col_idx += depth;
  }

  vbb->cache_idx = start;
}

void init_vbb_verify(vbb_t* vbb, unsigned int len, const faest_paramset_t* params,
                     const uint8_t* sig) {
  vbb->party                       = VERIFIER;
  vbb->params                      = params;
  const unsigned int lambda        = params->faest_param.lambda;
  const size_t lambda_bytes        = params->faest_param.lambda / 8;
  const unsigned int l             = params->faest_param.l;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
  const unsigned int tau           = vbb->params->faest_param.tau;

  vbb->iv        = dsignature_iv(sig, params);
  vbb->com_hash  = calloc(MAX_LAMBDA_BYTES * 2, sizeof(uint8_t));
  vbb->full_size = len >= ell_hat;
  vbb->sig       = sig;

  unsigned int row_count = len;
  vbb->row_count         = row_count;
  vbb->vole_cache        = calloc(row_count, lambda_bytes);

  // HASH cache
  unsigned int column_count = (size_t)vbb->row_count * (size_t)lambda_bytes / (size_t)ell_hat_bytes;
  assert(column_count >= 1);
  vbb->column_count = column_count;

  // Compute hcom
  const uint8_t* chall3 = dsignature_chall_3(vbb->sig, vbb->params);
  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  for (unsigned int i = 0; i < tau; ++i) {
    pdec[i] = dsignature_pdec(vbb->sig, i, vbb->params);
    com[i]  = dsignature_com(vbb->sig, i, vbb->params);
  }
  if (vbb->full_size) {
    vbb->v_buf = malloc(lambda_bytes);
    partial_vole_reconstruct_cmo(vbb->iv, chall3, pdec, com, vbb->com_hash, vbb->vole_cache,
                                 ell_hat, vbb->params, 0, lambda);
  } else {
    vole_reconstruct_hcom(vbb->iv, chall3, pdec, com, vbb->com_hash, ell_hat, vbb->params);
  }

  vbb->Dtilde_buf = malloc(lambda_bytes + UNIVERSAL_HASH_B);
}

void prepare_hash_verify(vbb_t* vbb) {
  if (vbb->full_size) {
    return;
  }
  recompute_hash_verify(vbb, 0, vbb->column_count);
}

const uint8_t* get_dtilde(vbb_t* vbb, unsigned int idx) {
  const unsigned int tau0         = vbb->params->faest_param.t0;
  const unsigned int tau1         = vbb->params->faest_param.t1;
  const size_t lambda_bytes       = vbb->params->faest_param.lambda / 8;
  const unsigned int utilde_bytes = lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int k0           = vbb->params->faest_param.k0;
  const unsigned int k1           = vbb->params->faest_param.k1;

  unsigned int i;
  unsigned int j;
  if (idx < k0 * tau0) {
    i = idx / k0;
    j = idx % k0;
  } else {
    i = tau0 + (idx - k0 * tau0) / k1;
    j = (idx - k0 * tau0) % k1;
  }

  uint8_t delta[MAX_DEPTH];
  ChalDec(dsignature_chall_3(vbb->sig, vbb->params), i, k0, tau0, k1, tau1, delta);

  memset(vbb->Dtilde_buf, 0, utilde_bytes);
  masked_xor_u8_array(vbb->Dtilde_buf, dsignature_u_tilde(vbb->sig, vbb->params), vbb->Dtilde_buf,
                      delta[j], utilde_bytes);

  return vbb->Dtilde_buf;
}

void prepare_aes_verify(vbb_t* vbb) {
  if (vbb->full_size) {
    const unsigned int tau           = vbb->params->faest_param.tau;
    const unsigned int t0            = vbb->params->faest_param.t0;
    const unsigned int k0            = vbb->params->faest_param.k0;
    const unsigned int t1            = vbb->params->faest_param.t1;
    const unsigned int k1            = vbb->params->faest_param.k1;
    const unsigned int lambda        = vbb->params->faest_param.lambda;
    const unsigned int l             = vbb->params->faest_param.l;
    const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;

    // Apply withness to CMO cache
    unsigned int size = vbb->params->faest_param.l;
    for (unsigned int i = 0, col = 0; i < tau; i++) {
      unsigned int depth = i < t0 ? k0 : k1;
      uint8_t decoded_challenge[MAX_DEPTH];
      ChalDec(dsignature_chall_3(vbb->sig, vbb->params), i, k0, t0, k1, t1, decoded_challenge);
      for (unsigned int j = 0; j < depth; j++, ++col) {
        if (decoded_challenge[j] == 1) {
          xor_u8_array(dsignature_d(vbb->sig, vbb->params), vbb->vole_cache + col * ell_hat_bytes,
                       vbb->vole_cache + col * ell_hat_bytes, (size + 7) / 8);
        }
      }
    }
    return;
  }

  // Compute first chunk with withness applied
  recompute_aes_verify(vbb, 0, vbb->row_count);
}

// Get voles for hashing
const uint8_t* get_vole_v_hash(vbb_t* vbb, unsigned int idx) {
  const unsigned int ellhat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;

  assert(idx < vbb->params->faest_param.lambda);
  if (!(idx >= vbb->cache_idx && idx < vbb->cache_idx + vbb->column_count)) {
    recompute_hash_prove(vbb, idx, vbb->column_count);
  }

  unsigned int offset = idx - vbb->cache_idx;
  return vbb->vole_cache + offset * ellhat_bytes;
}

const uint8_t* get_vole_q_hash(vbb_t* vbb, unsigned int idx) {
  unsigned int lambda = vbb->params->faest_param.lambda;
  const unsigned int ellhat =
      vbb->params->faest_param.l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;

  assert(idx < lambda);
  if (!(idx >= vbb->cache_idx && idx < vbb->cache_idx + vbb->column_count)) {
    recompute_hash_verify(vbb, idx, vbb->column_count);
  }

  unsigned int offset = idx - vbb->cache_idx;
  return vbb->vole_cache + offset * ellhat_bytes;
}

// Get voles for AES
static inline uint8_t* get_vole_aes(vbb_t* vbb, unsigned int idx) {
  unsigned int lambda       = vbb->params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat       = vbb->params->faest_param.l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;

  if (vbb->full_size) {
    memset(vbb->v_buf, 0, lambda_bytes);
    // Transpose on the fly into v_buf
    for (unsigned int column = 0; column != lambda; ++column) {
      ptr_set_bit(vbb->v_buf, ptr_get_bit(vbb->vole_cache + column * ellhat_bytes, idx), column);
    }
    return vbb->v_buf;
  }

  if (!(idx >= vbb->cache_idx && idx < vbb->cache_idx + vbb->row_count)) {
    if (vbb->party == VERIFIER) {
      recompute_aes_verify(vbb, idx, vbb->row_count);
    } else {
      recompute_aes_prove(vbb, idx, vbb->row_count);
    }
  }

  unsigned int offset = (idx - vbb->cache_idx) * lambda_bytes;
  return vbb->vole_cache + offset;
}

const bf256_t* get_vole_aes_256(vbb_t* vbb, unsigned int idx) {
  return (bf256_t*)get_vole_aes(vbb, idx);
}

const bf192_t* get_vole_aes_192(vbb_t* vbb, unsigned int idx) {
  return (bf192_t*)get_vole_aes(vbb, idx);
}

const bf128_t* get_vole_aes_128(vbb_t* vbb, unsigned int idx) {
  return (bf128_t*)get_vole_aes(vbb, idx);
}

const uint8_t* get_vole_u(vbb_t* vbb) {
  return vbb->vole_U;
}

const uint8_t* get_com_hash(vbb_t* vbb) {
  return vbb->com_hash;
}
