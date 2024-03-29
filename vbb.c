#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>

#include "macros.h"
#include "vbb.h"
#include "vole.h"
#include "vc.h"
#include "instances.h"
#include "faest.h"
#include "faest_aes.h"
#include "fields.h"
#include "parameters.h"

#define ACCESS_PATTERN_TEST 0

static void setup_vk_cache(vbb_t* vbb);

ATTR_CONST ATTR_ALWAYS_INLINE static inline bool is_em_variant(faest_paramid_t id) {
  return id > 6;
}

ATTR_CONST ATTR_ALWAYS_INLINE static bool is_column_cached(vbb_t* vbb, unsigned int index) {
  bool above_cache_start = index >= vbb->cache_idx;
  bool below_cache_end = index < vbb->cache_idx + vbb->column_count;
  return above_cache_start && below_cache_end;
}

ATTR_CONST ATTR_ALWAYS_INLINE static bool is_row_cached(vbb_t* vbb, unsigned int index) {
  bool above_cache_start = index >= vbb->cache_idx;
  bool below_cache_end = index < vbb->cache_idx + vbb->row_count;
  return above_cache_start && below_cache_end;
}


static void recompute_hash_sign(vbb_t* vbb, unsigned int start, unsigned int end) {
  const unsigned int lambda = vbb->params->faest_param.lambda;
  const unsigned int ell    = vbb->params->faest_param.l;
  const unsigned int ellhat = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int capped_end   = MIN(end, lambda);

  partial_vole_commit_cmo(vbb->root_key, vbb->iv, ellhat, 
                          start, capped_end,
                          vole_mode_v(vbb->vole_cache),
                          vbb->params);
  vbb->cache_idx = start;
}

static void recompute_aes_sign(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda = vbb->params->faest_param.lambda;
  const unsigned int ell    = vbb->params->faest_param.l;

  if(start >= ell){
    start = start - len + 1;
  }
  if (len >= ell + lambda) {
    start = 0;
  } else if (start + len > ell + lambda) {
    start = ell + lambda - len;
  }
  partial_vole_commit_rmo(vbb->root_key, vbb->iv, start, len, vbb->params, vbb->vole_cache);
  vbb->cache_idx = start;
}

// len is the number of OLE v's that is allowed to be stored in memory.
// Hence we store (at most) len*lambda in memory.
void init_vbb_sign(vbb_t* vbb, unsigned int len, const uint8_t* root_key, const uint8_t* iv,
                   uint8_t* c, const faest_paramset_t* params) {
  const unsigned int lambda       = params->faest_param.lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ellhat       = params->faest_param.l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ellhat_bytes = (ellhat + 7) / 8;
  const unsigned int row_count    = len;
  const unsigned int column_count = (size_t)row_count * (size_t)lambda_bytes / (size_t)ellhat_bytes;
  assert(column_count >= 1);

  vbb->party        = SIGNER;
  vbb->iv           = iv;
  vbb->com_hash     = calloc(MAX_LAMBDA_BYTES * 2, sizeof(uint8_t));
  vbb->params       = params;
  vbb->root_key     = root_key;
  vbb->full_size    = len >= ellhat;
  vbb->vole_U       = malloc(ellhat_bytes);
  vbb->row_count    = row_count;
  vbb->vole_cache   = calloc(row_count, lambda_bytes);
  vbb->v_buf        = malloc(lambda_bytes);
  vbb->column_count = column_count;

  // Setup vk_buf if we are not in an EM variant
  if (!is_em_variant(vbb->params->faest_paramid)) {
    vbb->vk_buf = malloc(lambda_bytes);
  }
  
  sign_vole_mode_ctx_t mode = vbb->full_size ? vole_mode_all_sign(vbb->vole_cache, vbb->vole_U, vbb->com_hash, c)
                                    : vole_mode_u_hcom_c(vbb->vole_U, vbb->com_hash, c);

  partial_vole_commit_cmo(vbb->root_key, vbb->iv, ellhat, 0, lambda, mode, vbb->params);
}

void prepare_hash_sign(vbb_t* vbb) {
  if (vbb->full_size) {
    vbb->cache_idx = 0;
    return;
  }
  recompute_hash_sign(vbb, 0, vbb->column_count);
}

void prepare_aes_sign(vbb_t* vbb) {
  if (vbb->full_size) {
    vbb->cache_idx = 0;
  } else {
    recompute_aes_sign(vbb, 0, vbb->row_count);
  }
  setup_vk_cache(vbb);
}

void vector_open_ondemand(vbb_t* vbb, unsigned int idx, const uint8_t* s_, uint8_t* sig_pdec,
                          uint8_t* sig_com, unsigned int depth) {
  const unsigned int lambda       = vbb->params->faest_param.lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int tau          = vbb->params->faest_param.tau;
  uint8_t* expanded_keys          = malloc(tau * lambda_bytes);
  prg(vbb->root_key, vbb->iv, expanded_keys, lambda, lambda_bytes * tau);

  vec_com_t vec_com;
  vector_commitment(expanded_keys + lambda_bytes * idx, lambda, depth, NULL, &vec_com);
  vector_open(&vec_com, s_, sig_pdec, sig_com, depth, vbb->iv, lambda);
  free(expanded_keys);
}

static inline void apply_correction_values_cmo(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda        = vbb->params->faest_param.lambda;
  const unsigned int ell           = vbb->params->faest_param.l;
  const unsigned int ell_hat       = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
  const unsigned int tau           = vbb->params->faest_param.tau;
  const unsigned int tau0          = vbb->params->faest_param.t0;
  const unsigned int k0            = vbb->params->faest_param.k0;
  const unsigned int k1            = vbb->params->faest_param.k1;

  const uint8_t* chall3 = dsignature_chall_3(vbb->sig, vbb->params);
  const uint8_t* c      = dsignature_c(vbb->sig, 0, vbb->params);
  unsigned int col_idx  = k0;
  for (unsigned int i = 1; i < tau; i++) {
    const unsigned int depth = i < tau0 ? k0 : k1;
    if (start >= col_idx + depth) {
      col_idx += depth;
      continue;
    }
    uint8_t delta[MAX_DEPTH];
    ChalDec(chall3, i, vbb->params->faest_param.k0, vbb->params->faest_param.t0,
            vbb->params->faest_param.k1, vbb->params->faest_param.t1, delta);

    for (unsigned int d = 0; d < depth; d++, col_idx++) {
      if (start > col_idx) {
        continue;
      }
      masked_xor_u8_array(
          vbb->vole_cache + (col_idx - start) * ell_hat_bytes, c + (i - 1) * ell_hat_bytes,
          vbb->vole_cache + (col_idx - start) * ell_hat_bytes, delta[d], ell_hat_bytes);

      if (col_idx + 1 >= start + len) {
        col_idx++;
        break;
      }
    }
    if (col_idx >= start + len) {
      break;
    }
  }
  vbb->cache_idx = start;
}

static void setup_pdec_com(vbb_t* vbb, const uint8_t** pdec, const uint8_t** com) {
  const unsigned int tau = vbb->params->faest_param.tau;
  for (unsigned int i = 0; i < tau; ++i) {
    pdec[i] = dsignature_pdec(vbb->sig, i, vbb->params);
    com[i]  = dsignature_com(vbb->sig, i, vbb->params);
  }
}

// Verifier implementation
static void recompute_hash_verify(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda  = vbb->params->faest_param.lambda;
  const unsigned int ell     = vbb->params->faest_param.l;
  const unsigned int ell_hat = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int amount  = MIN(len, lambda - start);
  const uint8_t* chall3      = dsignature_chall_3(vbb->sig, vbb->params);

  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  setup_pdec_com(vbb, pdec, com);

  partial_vole_reconstruct_cmo(vbb->iv, chall3, pdec, com, ell_hat,
                               start, amount,
                               vole_mode_q(vbb->vole_cache),
                               vbb->params);
  apply_correction_values_cmo(vbb, start, amount);
  vbb->cache_idx = start;
}

static void apply_correction_values_rmo(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda        = vbb->params->faest_param.lambda;
  const unsigned int lambda_bytes  = lambda / 8;
  const unsigned int ell           = vbb->params->faest_param.l;
  const unsigned int tau           = vbb->params->faest_param.tau;
  const unsigned int ell_hat       = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int tau0          = vbb->params->faest_param.t0;
  const unsigned int k0            = vbb->params->faest_param.k0;
  const unsigned int k1            = vbb->params->faest_param.k1;
  const uint8_t* chall3            = dsignature_chall_3(vbb->sig, vbb->params);
  const uint8_t* c                 = dsignature_c(vbb->sig, 0, vbb->params);
  uint8_t buf[MAX_LAMBDA_BYTES]    = {0};
  for (unsigned int row_idx = 0; row_idx < len; row_idx++) {
    memset(buf, 0, sizeof(buf));
    uint8_t packed_byte      = 0;
    unsigned int bit_counter = k0;

    for (unsigned int t = 1; t < tau; t++) {
      unsigned int depth   = t < tau0 ? k0 : k1;
      const uint8_t* c_idx = c + (t - 1) * ell_hat_bytes;

      unsigned int abs_idx = row_idx + start;
      unsigned int c_byte  = abs_idx / 8;
      unsigned int c_bit   = abs_idx % 8;
      uint8_t bit          = c_idx[c_byte] >> c_bit & 1;

      for (unsigned int i = 0; i < depth; i++) {
        packed_byte = packed_byte >> 1;
        packed_byte |= (bit << 7);
        bit_counter++;
        if (bit_counter % 8 == 0) {
          buf[bit_counter / 8 - 1] = packed_byte;
          packed_byte              = 0;
        }
      }
    }

    // AND buf with Delta
    for (unsigned int i = 0; i < lambda_bytes; i++) {
      buf[i] &= chall3[i];
      vbb->vole_cache[row_idx * lambda_bytes + i] ^= buf[i];
    }
  }
}

static void apply_witness_values_rmo(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda       = vbb->params->faest_param.lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int ell          = vbb->params->faest_param.l;
  const unsigned int tau          = vbb->params->faest_param.tau;
  const unsigned int tau0         = vbb->params->faest_param.t0;
  const unsigned int tau1         = vbb->params->faest_param.t1;
  const unsigned int k0           = vbb->params->faest_param.k0;
  const unsigned int k1           = vbb->params->faest_param.k1;
  const uint8_t* d                = dsignature_d(vbb->sig, vbb->params);
  unsigned int full_col_idx       = 0;
  unsigned int end_row_idx        = MIN(len, ell - start);
  uint8_t delta[MAX_DEPTH];

  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;
    ChalDec(dsignature_chall_3(vbb->sig, vbb->params), i, k0, tau0, k1, tau1, delta);
    for (unsigned int col_idx = 0; col_idx < depth; col_idx++) {
      uint8_t delta_i = delta[col_idx];
      if (delta_i == 0) {
        continue;
      }
      for (unsigned int row_idx = 0; row_idx < end_row_idx; row_idx++) {
        unsigned int q_byte = (full_col_idx + col_idx) / 8;
        unsigned int d_byte = (row_idx + start) / 8;
        unsigned int d_bit  = (row_idx + start) % 8;
        uint8_t bit         = d[d_byte] >> d_bit & 1;
        if (bit == 0) {
          continue;
        }
        vbb->vole_cache[row_idx * lambda_bytes + q_byte] ^= bit << (full_col_idx + col_idx) % 8;
      }
    }
    full_col_idx += depth;
  }
}

static void apply_witness_values_cmo(vbb_t* vbb) {
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
}

static void recompute_aes_verify(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int lambda  = vbb->params->faest_param.lambda;
  const unsigned int ell     = vbb->params->faest_param.l;
  const unsigned int ell_hat = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;

  if(start >= ell){
    start = start - len + 1;
  }
  if (len >= ell + lambda) {
    start = 0;
  } else if (start + len > ell + lambda) {
    start = ell + lambda - len;
  }

  const uint8_t* chall3 = dsignature_chall_3(vbb->sig, vbb->params);
  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  setup_pdec_com(vbb, pdec, com);
  partial_vole_reconstruct_rmo(vbb->iv, chall3, pdec, com, vbb->vole_cache, ell_hat, vbb->params,
                               start, len);
  apply_correction_values_rmo(vbb, start, len);
  apply_witness_values_rmo(vbb, start, len);
  vbb->cache_idx = start;
}

void init_vbb_verify(vbb_t* vbb, unsigned int len, const faest_paramset_t* params,
                     const uint8_t* sig) {
  const unsigned int lambda        = params->faest_param.lambda;
  const unsigned int lambda_bytes  = params->faest_param.lambda / 8;
  const unsigned int l             = params->faest_param.l;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
  const unsigned int row_count     = len;
  const unsigned int column_count =
      (size_t)row_count * (size_t)lambda_bytes / (size_t)ell_hat_bytes;
  assert(column_count >= 1);

  vbb->party        = VERIFIER;
  vbb->params       = params;
  vbb->iv           = dsignature_iv(sig, params);
  vbb->com_hash     = calloc(MAX_LAMBDA_BYTES * 2, sizeof(uint8_t));
  vbb->full_size    = len >= ell_hat;
  vbb->sig          = sig;
  vbb->row_count    = row_count;
  vbb->vole_cache   = calloc(row_count, lambda_bytes);
  vbb->column_count = column_count;
  vbb->Dtilde_buf   = malloc(lambda_bytes + UNIVERSAL_HASH_B);
  vbb->v_buf        = malloc(lambda_bytes);

  // Setup vk_buf if we are not in an EM variant
  if (!is_em_variant(vbb->params->faest_paramid)) {
    vbb->vk_buf = malloc(lambda_bytes);
  }

  const uint8_t* chall3 = dsignature_chall_3(vbb->sig, vbb->params);
  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  setup_pdec_com(vbb, pdec, com);

  verify_vole_mode_ctx_t vole_mode = (vbb->full_size) ? vole_mode_all_verify(vbb->vole_cache, vbb->com_hash)
                                                      : vole_mode_hcom(vbb->com_hash);
  partial_vole_reconstruct_cmo(vbb->iv, chall3, pdec, com, ell_hat, 0, lambda, vole_mode, vbb->params);
  if (vbb->full_size) {
    apply_correction_values_cmo(vbb, 0, lambda);
  }
}

void prepare_hash_verify(vbb_t* vbb) {
  if (vbb->full_size) {
    vbb->cache_idx = 0;
    return;
  }
  recompute_hash_verify(vbb, 0, vbb->column_count);
}

const uint8_t* get_dtilde(vbb_t* vbb, unsigned int idx) {
  const unsigned int tau0         = vbb->params->faest_param.t0;
  const unsigned int tau1         = vbb->params->faest_param.t1;
  const unsigned int lambda       = vbb->params->faest_param.lambda;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int utilde_bytes = lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int k0           = vbb->params->faest_param.k0;
  const unsigned int k1           = vbb->params->faest_param.k1;

  unsigned int i = 0;
  unsigned int j = 0;
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
    apply_witness_values_cmo(vbb);
    vbb->cache_idx = 0;
  } else {
    recompute_aes_verify(vbb, 0, vbb->row_count);
  }
  setup_vk_cache(vbb);
}

// Get voles for hashing
const uint8_t* get_vole_v_hash(vbb_t* vbb, unsigned int idx) {
  const unsigned int lambda        = vbb->params->faest_param.lambda;
  const unsigned int ell           = vbb->params->faest_param.l;
  const unsigned int ell_hat       = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;

  if (!is_column_cached(vbb, idx)) {
    unsigned int cmo_budget = vbb->column_count;
    recompute_hash_sign(vbb, idx, idx + cmo_budget);
  }
  const unsigned int offset = idx - vbb->cache_idx;

  return vbb->vole_cache + offset * ell_hat_bytes;
}

const uint8_t* get_vole_q_hash(vbb_t* vbb, unsigned int idx) {
  const unsigned int lambda        = vbb->params->faest_param.lambda;
  const unsigned int ell           = vbb->params->faest_param.l;
  const unsigned int ell_hat       = ell + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;


  if (!is_column_cached(vbb, idx)) {
    unsigned int cmo_budget = vbb->column_count;
    recompute_hash_verify(vbb, idx, cmo_budget);
  }

  unsigned int offset = idx - vbb->cache_idx;
  return vbb->vole_cache + offset * ell_hat_bytes;
}

// Get voles for AES
static inline uint8_t* get_vole_aes(vbb_t* vbb, unsigned int idx) {
  unsigned int lambda       = vbb->params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat       = vbb->params->faest_param.l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;

#if ACCESS_PATTERN_TEST
  // Store the idx value in a file
  FILE* file;
  if (vbb->party == VERIFIER) {
    file = fopen("access_pattern_verifier.txt", "a");
  } else {
    file = fopen("access_pattern_prover.txt", "a");
  }
  fprintf(file, "%d\n", idx);
  fclose(file);
#endif

  if (vbb->full_size) {
    memset(vbb->v_buf, 0, lambda_bytes);
    // Transpose on the fly into v_buf
    for (unsigned int column = 0; column != lambda; ++column) {
      ptr_set_bit(vbb->v_buf, ptr_get_bit(vbb->vole_cache + column * ellhat_bytes, idx), column);
    }
    return vbb->v_buf;
  }

  if (!is_row_cached(vbb, idx)) {
    unsigned int rmo_budget = vbb->row_count;
    if (vbb->party == VERIFIER) {
      recompute_aes_verify(vbb, idx, rmo_budget);
    } else {
      recompute_aes_sign(vbb, idx, rmo_budget);
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

void clean_vbb(vbb_t* vbb) {
  free(vbb->vole_cache);
  free(vbb->com_hash);

  if (vbb->full_size) {
    free(vbb->v_buf);
  }

  if (vbb->party == VERIFIER) {
    free(vbb->Dtilde_buf);
  } else {
    free(vbb->vole_U);
  }

  // V_k cache
  if (!is_em_variant(vbb->params->faest_paramid)) {
    free(vbb->vk_buf);
    if (!vbb->full_size) {
      free(vbb->vk_cache);
    }
  }
}

// V_k cache

static void setup_vk_cache(vbb_t* vbb) {
  unsigned int lambda_bytes = vbb->params->faest_param.lambda / 8;
  if (is_em_variant(vbb->params->faest_paramid)) {
    return;
  }

  vbb->vk_cache = calloc(vbb->params->faest_param.Lke, lambda_bytes);

  for (unsigned int i = 0; i < vbb->params->faest_param.Lke; i++) {
    unsigned int offset = i * lambda_bytes;
    memcpy(vbb->vk_cache + offset, get_vole_aes(vbb, i), lambda_bytes);
  }
}

static inline uint8_t* get_vk(vbb_t* vbb, unsigned int idx) {
  assert(idx < vbb->params->faest_param.Lke);
  unsigned int offset = idx * (vbb->params->faest_param.lambda / 8);
  return (vbb->vk_cache + offset);
}

const bf128_t* get_vk_128(vbb_t* vbb, unsigned int idx) {
  if (idx < FAEST_128F_LAMBDA) {
    const bf128_t* vk = (bf128_t*)get_vk(vbb, idx);
    memcpy(vbb->vk_buf, vk, sizeof(bf128_t));
    return (bf128_t*)vbb->vk_buf;
  }

  unsigned int j = idx / 32 + FAEST_128F_Nwd;
  if ((j % FAEST_128F_Nwd) == 0 || (FAEST_128F_Nwd > 6 && (j % FAEST_128F_Nwd) == 4)) {
    unsigned int i_wd       = FAEST_128F_LAMBDA;
    unsigned int factor_128 = (idx / 128) - 1;
    unsigned int offset_128 = idx % 128;
    unsigned int index      = i_wd + factor_128 * 32 + offset_128;
    const bf128_t* vk       = (bf128_t*)get_vk(vbb, index);
    memcpy(vbb->vk_buf, vk, sizeof(bf128_t));
    return (bf128_t*)vbb->vk_buf;
  }

  // Lhs recursive call
  const bf128_t* lhs_ptr = get_vk_128(vbb, idx - FAEST_128F_Nwd * 32);
  bf128_t lhs            = *lhs_ptr;
  // Rhs recursive call
  const bf128_t* rhs_ptr = get_vk_128(vbb, idx - 32);
  bf128_t rhs            = *rhs_ptr;

  bf128_t vk = bf128_add(lhs, rhs);
  memcpy(vbb->vk_buf, &vk, sizeof(bf128_t));
  return (bf128_t*)vbb->vk_buf;
}

const bf192_t* get_vk_192(vbb_t* vbb, unsigned int idx) {
  if (idx < FAEST_192F_LAMBDA) {
    const bf192_t* vk = (bf192_t*)get_vk(vbb, idx);
    memcpy(vbb->vk_buf, vk, sizeof(bf192_t));
    return (bf192_t*)vbb->vk_buf;
  }

  unsigned int j = idx / 32 + FAEST_192F_Nwd;
  if ((j % FAEST_192F_Nwd) == 0 || (FAEST_192F_Nwd > 6 && (j % FAEST_192F_Nwd) == 4)) {
    unsigned int i_wd       = FAEST_192F_LAMBDA;
    unsigned int factor_192 = (idx / 192) - 1;
    unsigned int offset_192 = idx % 192;
    unsigned int index      = i_wd + factor_192 * 32 + offset_192;
    const bf192_t* vk       = (bf192_t*)get_vk(vbb, index);
    memcpy(vbb->vk_buf, vk, sizeof(bf192_t));
    return (bf192_t*)vbb->vk_buf;
  }

  // Lhs recursive call
  const bf192_t* lhs_ptr = get_vk_192(vbb, idx - FAEST_192F_Nwd * 32);
  bf192_t lhs            = *lhs_ptr;
  // Rhs recursive call
  const bf192_t* rhs_ptr = get_vk_192(vbb, idx - 32);
  bf192_t rhs            = *rhs_ptr;

  bf192_t vk = bf192_add(lhs, rhs);
  memcpy(vbb->vk_buf, &vk, 3 * sizeof(uint64_t));
  return (bf192_t*)vbb->vk_buf;
}

const bf256_t* get_vk_256(vbb_t* vbb, unsigned int idx) {
  if (idx < FAEST_256F_LAMBDA) {
    const bf256_t* vk = (bf256_t*)get_vk(vbb, idx);
    memcpy(vbb->vk_buf, vk, sizeof(bf256_t));
    return (bf256_t*)vbb->vk_buf;
  }

  // We go from j=N_k to j=4(R+1)
  // In our case this is j=8 to j=4(14+1)=60
  // Based on each j we have 32 tags in vk
  unsigned int j = idx / 32 + FAEST_256F_Nwd;
  if ((j % FAEST_256F_Nwd) == 0 || (FAEST_256F_Nwd > 6 && (j % FAEST_256F_Nwd) == 4)) {
    unsigned int i_wd       = FAEST_256F_LAMBDA;
    unsigned int factor_128 = (idx / 128) - 2;
    unsigned int offset_128 = idx % 128;
    unsigned int index      = i_wd + factor_128 * 32 + offset_128;
    const bf256_t* vk       = (bf256_t*)get_vk(vbb, index);
    memcpy(vbb->vk_buf, vk, sizeof(bf256_t));
    return (bf256_t*)vbb->vk_buf;
  }

  // Lhs recursive call
  const bf256_t* lhs_ptr = get_vk_256(vbb, idx - FAEST_256F_Nwd * 32);
  bf256_t lhs            = *lhs_ptr;
  // Rhs recursive call
  const bf256_t* rhs_ptr = get_vk_256(vbb, idx - 32);
  bf256_t rhs            = *rhs_ptr;

  bf256_t vk = bf256_add(lhs, rhs);
  memcpy(vbb->vk_buf, &vk, sizeof(bf256_t));
  return (bf256_t*)vbb->vk_buf;
}
