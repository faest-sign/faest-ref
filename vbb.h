#ifndef VBB_H
#define VBB_H

#include <stdint.h>
#include <stdbool.h>

#include "vole.h"

typedef struct {
  int start_idx;
  int len; // TODO: Compile time
  vec_com_t* vecCom;
  uint8_t** vole_V_cache;
  uint8_t* vole_U;
  uint8_t* com_hash;
  const uint8_t* root_key;
  bool initialized; // TODO: remove if recompute when initialize
  const faest_paramset_t* params;
  const uint8_t* iv;
  const uint8_t* c;
} vbb_t;

void init_vbb(vbb_t* vbb, int len, const uint8_t* root_key, const uint8_t* iv, const uint8_t* c,
              const faest_paramset_t* params);
uint8_t* get_vole_v(vbb_t* vbb, int idx);
uint8_t* get_vole_u(vbb_t* vbb);
uint8_t* get_com_hash(vbb_t* vbb);

#endif
