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
  const unsigned int ell_hat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;

  unsigned int amount       = MIN(len, vbb->params->faest_param.lambda - start);
  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  partial_vole_commit_cmo(vbb->root_key, vbb->iv, ell_hat, vbb->params, sVecCom,
                          vbb->vole_V_cache_hash, start, amount);
  free(sVecCom);

  vbb->start_idx_hash = start;
}

static void recompute_prove(vbb_t* vbb, unsigned int start, unsigned int len) {
  if (len >= vbb->params->faest_param.l + vbb->params->faest_param.lambda) {
    start = 0;
  } else if (start + len > vbb->params->faest_param.l + vbb->params->faest_param.lambda) {
    start = vbb->params->faest_param.l + vbb->params->faest_param.lambda - len;
  }
  memset(vbb->vole_V_cache_prove, 0, ((size_t)len) * (size_t)vbb->params->faest_param.lambda / 8);

  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  partial_vole_commit_rmo(vbb->root_key, vbb->iv, start, len, vbb->params, sVecCom,
                          vbb->vole_V_cache_prove);
  free(sVecCom);

  vbb->start_idx_prove = start;
}

// len is the number of OLE v's that is allowed to be stored in memory.
// Hence we store len*lambda in memory.
void init_vbb(vbb_t* vbb, unsigned int len, const uint8_t* root_key, const uint8_t* iv, uint8_t* c,
              const faest_paramset_t* params) {
  vbb->len = len;
  vbb->iv  = iv;
  vbb->c   = c;
  // uint8_t hcom[MAX_LAMBDA_BYTES * 2]
  vbb->com_hash = calloc(MAX_LAMBDA_BYTES * 2, sizeof(uint8_t));
  vbb->params   = params;
  vbb->root_key = root_key;
  vbb->vecCom   = calloc(params->faest_param.tau, sizeof(vec_com_t));

  const unsigned int ell_hat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  vbb->vole_U                      = malloc(ell_hat_bytes);

  // FIXME: FAKE IT - make all variant work pre change
  vbb->vole_V_cache    = malloc(params->faest_param.lambda * sizeof(uint8_t*));
  vbb->vole_V_cache[0] = calloc(params->faest_param.lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < params->faest_param.lambda; ++i) {
    vbb->vole_V_cache[i] = vbb->vole_V_cache[0] + i * ell_hat_bytes;
  }
  vole_commit(vbb->root_key, vbb->iv, ell_hat, vbb->params, vbb->com_hash, vbb->vecCom, vbb->c,
              vbb->vole_U, vbb->vole_V_cache);
  // PROVE cache
  vbb->vole_V_cache_prove = calloc(len, vbb->params->faest_param.lambda / 8);
  vbb->start_idx_prove    = 0;
  recompute_prove(vbb, 0, len);

  // Setup u hcom c.
  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  vole_commit_u_hcom_c(vbb->root_key, vbb->iv, ell_hat, vbb->params, vbb->com_hash, sVecCom, vbb->c,
                       vbb->vole_U);
  free(sVecCom);

  // HASH cache
  unsigned int long_len     = ((uint32_t)vbb->len * (uint32_t)params->faest_param.lambda) / ell_hat;
  vbb->long_len             = long_len;
  vbb->vole_V_cache_hash    = malloc(long_len * sizeof(uint8_t*));
  vbb->vole_V_cache_hash[0] = calloc(long_len, ell_hat_bytes);
  for (unsigned int i = 1; i < long_len; ++i) {
    vbb->vole_V_cache_hash[i] = vbb->vole_V_cache_hash[0] + i * ell_hat_bytes;
  }
  vbb->start_idx_hash = 0;
  recompute_hash(vbb, 0, long_len);

  // TODO: Make cleanup to free all malloc
}

uint8_t* get_vole_v_hash(vbb_t* vbb, unsigned int idx) {
  assert(idx < vbb->params->faest_param.lambda);
  if (!(idx >= vbb->start_idx_hash && idx < vbb->start_idx_hash + vbb->long_len)) {
    recompute_hash(vbb, idx, vbb->long_len);
  }

  unsigned int offset = idx - vbb->start_idx_hash;
  return vbb->vole_V_cache_hash[offset];
}

bf256_t* get_vole_v_prove_256(vbb_t* vbb, unsigned int idx) {
  assert(idx < vbb->params->faest_param.l + vbb->params->faest_param.lambda);
  if (!(idx >= vbb->start_idx_prove && idx < vbb->start_idx_prove + vbb->len)) {
    recompute_prove(vbb, idx, vbb->len);
  }

  unsigned int offset = idx - vbb->start_idx_prove;
  return ((bf256_t*)(vbb->vole_V_cache_prove)) + offset;
}

bf192_t* get_vole_v_prove_192(vbb_t* vbb, unsigned int idx) {
  assert(idx < vbb->params->faest_param.l + vbb->params->faest_param.lambda);
  if (!(idx >= vbb->start_idx_prove && idx < vbb->start_idx_prove + vbb->len)) {
    recompute_prove(vbb, idx, vbb->len);
  }

  unsigned int offset = idx - vbb->start_idx_prove;
  return (bf192_t*)(vbb->vole_V_cache_prove + offset*vbb->params->faest_param.lambda/8);
}

uint8_t* get_vole_u(vbb_t* vbb) {
  return vbb->vole_U;
}

uint8_t* get_com_hash(vbb_t* vbb) {
  return vbb->com_hash;
}

// TODO - refactor this stuff
void vector_open_ondemand(vbb_t* vbb, unsigned int idx, const uint8_t* s_, uint8_t* sig_pdec, uint8_t* sig_com, unsigned int depth) {
  unsigned int lambda       = vbb->params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int tau          = vbb->params->faest_param.tau;
  uint8_t* expanded_keys = malloc(tau * lambda_bytes);
  prg(vbb->root_key, vbb->iv, expanded_keys, lambda, lambda_bytes * tau);

  stream_vec_com_t* sVecCom = calloc(vbb->params->faest_param.tau, sizeof(stream_vec_com_t));
  stream_vector_commitment(expanded_keys + lambda_bytes * idx, lambda, &sVecCom[idx], depth);
  stream_vector_open(&sVecCom[idx], s_, sig_pdec, sig_com, depth, vbb->iv, lambda);
  free(sVecCom);
  free(expanded_keys);
}
