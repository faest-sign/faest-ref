#include <stddef.h>
#include <stdio.h>

#include "macros.h"
#include "vbb.h"
#include "vole.h"
#include "vc.h"
#include "instances.h"
#include "faest.h"

static void recompute_hash(vbb_t* vbb, int start, int len) {
  const unsigned int ell_hat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;

  // TODO: FAKE IT
  uint8_t** V = malloc(vbb->params->faest_param.lambda * sizeof(uint8_t*));
  V[0]        = calloc(vbb->params->faest_param.lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < vbb->params->faest_param.lambda; ++i) {
    V[i] = V[0] + i * ell_hat_bytes;
  }

  // TODO: Modify vole_commit to ouput specified size
  vole_commit(vbb->root_key, vbb->iv, ell_hat, vbb->params, vbb->com_hash, vbb->vecCom, vbb->c,
              vbb->vole_U, V);

  size_t amount = MIN(len, vbb->params->faest_param.lambda - start);
  memcpy(vbb->vole_V_cache_hash[0], V[start], amount * ell_hat_bytes);

  free(V[0]);
  free(V);
  vbb->start_idx = start;
}

// len is the number of OLE v's that is allowed to be stored in memory.
// Hence we store len*lambda in memory.
void init_vbb(vbb_t* vbb, int len, const uint8_t* root_key, const uint8_t* iv, const uint8_t* c,
              const faest_paramset_t* params) {
  vbb->len = len;
  vbb->iv  = iv;
  vbb->c   = c;
  // uint8_t hcom[MAX_LAMBDA_BYTES * 2]
  vbb->com_hash    = calloc(MAX_LAMBDA_BYTES * 2, sizeof(uint8_t));
  vbb->initialized = 0;
  vbb->params      = params;
  vbb->root_key    = root_key;
  vbb->start_idx   = 0;
  vbb->vecCom      = calloc(params->faest_param.tau, sizeof(vec_com_t));

  const unsigned int ell_hat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  vbb->vole_U                      = malloc(ell_hat_bytes);

  // FIXME: correct sizes
  // FIXME: Transposed matrix
  vbb->vole_V_cache    = malloc(params->faest_param.lambda * sizeof(uint8_t*));
  vbb->vole_V_cache[0] = calloc(params->faest_param.lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < params->faest_param.lambda; ++i) {
    vbb->vole_V_cache[i] = vbb->vole_V_cache[0] + i * ell_hat_bytes;
  }

  vole_commit(vbb->root_key, vbb->iv, ell_hat, vbb->params, vbb->com_hash, vbb->vecCom, vbb->c,
              vbb->vole_U, vbb->vole_V_cache);

  // FIXME: Use len to limit size of vole_V_cache_hash
  // Len = mem / \lambda 
  // This case: mem / \ell_hat_bytes
  // ell_hat * lambda / (ell_hat/8) = ell_hat * lambda * (8/ell_hat) = lambda * 8 = 1024 = 128*8
  int long_len = (vbb->len * params->faest_param.lambda) / ell_hat;
  vbb->long_len = long_len;
  vbb->vole_V_cache_hash    = malloc(long_len * sizeof(uint8_t*));
  vbb->vole_V_cache_hash[0] = calloc(long_len, ell_hat_bytes);
  for (unsigned int i = 1; i < long_len; ++i) {
    vbb->vole_V_cache_hash[i] = vbb->vole_V_cache_hash[0] + i * ell_hat_bytes;
  }

  recompute_hash(vbb, 0, long_len);
  // TODO: Make cleanup to free all malloc
}

uint8_t* get_vole_v_hash(vbb_t* vbb, int idx) {
  if (!(idx >= vbb->start_idx && idx < vbb->start_idx + vbb->long_len)) {
    recompute_hash(vbb, idx, vbb->long_len);
  }

  int offset = idx - vbb->start_idx;
  printf("offset %d", offset);
  return vbb->vole_V_cache_hash[offset];
}

uint8_t* get_vole_v_prove(vbb_t* vbb, int idx) {
  // TODO: FAKE IT
  int offset = idx - vbb->start_idx;
  return vbb->vole_V_cache[offset];
}

uint8_t* get_vole_u(vbb_t* vbb) {
  return vbb->vole_U;
}

uint8_t* get_com_hash(vbb_t* vbb) {
  return vbb->com_hash;
}