#include <stddef.h>

#include "vbb.h"
#include "vole.h"
#include "vc.h"
#include "instances.h"
#include "faest.h"

static void recompute(vbb_t* vbb, int start) { // FIXME: use start
  const unsigned int ell_hat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;

  vole_commit(vbb->root_key, vbb->iv, ell_hat, vbb->params, vbb->com_hash, vbb->vecCom, vbb->c,
              vbb->vole_U, vbb->vole_V_cache);
}

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
  vbb->vole_V_cache    = malloc(params->faest_param.lambda * sizeof(uint8_t*));
  vbb->vole_V_cache[0] = calloc(params->faest_param.lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < params->faest_param.lambda; ++i) {
    vbb->vole_V_cache[i] = vbb->vole_V_cache[0] + i * ell_hat_bytes;
  }

  recompute(vbb, 0);
}

uint8_t* get_vole_v(vbb_t* vbb, int idx) {
  if (!(idx >= vbb->start_idx && idx < vbb->start_idx + vbb->len)) {
    recompute(vbb, idx);
  }

  const unsigned int ell_hat =
      vbb->params->faest_param.l + vbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;

  int offset = idx - vbb->start_idx;
  return vbb->vole_V_cache[offset];
}

uint8_t* get_vole_u(vbb_t* vbb) {
  return vbb->vole_U;
}

uint8_t* get_com_hash(vbb_t* vbb) {
  return vbb->com_hash;
}