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

static void recompute_hash(vbb_t* vbb, unsigned int start, unsigned int len) {
  const unsigned int ellhat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;

  unsigned int amount       = MIN(len, vbb->params->faest_param.lambda - start);
  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  partial_vole_commit_cmo(vbb->root_key, vbb->iv, ellhat, vbb->params, sVecCom, vbb->vole_V_cache,
                          start, amount);
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
  vbb->vole_U               = malloc(ellhat_bytes);

  // PROVE cache
  // FIXME - would MAX(len, vbb->params->faest_param.Lenc) be correct for all variants?
  unsigned int row_count = len;
  vbb->row_count         = row_count;
  vbb->vole_V_cache      = calloc(row_count, lambda_bytes);

  // Setup u hcom c.
  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  vole_commit_u_hcom_c(vbb->root_key, vbb->iv, ellhat, vbb->params, vbb->com_hash, sVecCom, c,
                       vbb->vole_U);
  free(sVecCom);
  // HASH cache
  unsigned int column_count = (size_t)vbb->row_count * (size_t)lambda_bytes / (size_t)ellhat_bytes;
  printf("%lu\t%u\n", column_count, lambda);
  assert(column_count >= 1);
  vbb->column_count = column_count;

  // TODO: Make cleanup to free all malloc
}

void prepare_hash(vbb_t* vbb) {
  vbb->cache_idx = 0;
  recompute_hash(vbb, 0, vbb->column_count);
}

void prepare_prove(vbb_t* vbb) {
  unsigned int len = vbb->row_count;

  vbb->cache_idx = 0;
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
  if (!(idx >= vbb->cache_idx && idx < vbb->cache_idx + vbb->row_count)) {
    recompute_prove(vbb, idx, vbb->row_count);
  }

  unsigned int offset = (idx - vbb->cache_idx) * (vbb->params->faest_param.lambda / 8);
  return vbb->vole_V_cache + offset;
}

inline bf256_t* get_vole_v_prove_256(vbb_t* vbb, unsigned int idx) {
  return (bf256_t*)get_vole_v_prove(vbb, idx);
}

inline bf192_t* get_vole_v_prove_192(vbb_t* vbb, unsigned int idx) {
  return (bf192_t*)get_vole_v_prove(vbb, idx);
}

inline bf128_t* get_vole_v_prove_128(vbb_t* vbb, unsigned int idx) {
  return (bf128_t*)get_vole_v_prove(vbb, idx);
}

inline uint8_t* get_vole_u(vbb_t* vbb) {
  return vbb->vole_U;
}

inline uint8_t* get_com_hash(vbb_t* vbb) {
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
