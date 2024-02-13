#ifndef VBB_H
#define VBB_H

#include <stdint.h>
#include <stdbool.h>

#include "vole.h"

typedef struct vbb_t vbb_t;
#include "fields.h"
#include "faest_aes.h"

struct vbb_t {
  unsigned int start_idx_hash;
  unsigned int start_idx_prove;
  unsigned int len; // TODO: Compile time
  unsigned int long_len;
  vec_com_t* vecCom;
  uint8_t** vole_V_cache;
  uint8_t** vole_V_cache_hash;
  bf256_t* vole_V_cache_prove;
  uint8_t* vole_U;
  uint8_t* com_hash;
  const uint8_t* root_key;
  const faest_paramset_t* params;
  const uint8_t* iv;
  uint8_t* c;
};

void init_vbb(vbb_t* vbb, unsigned int len, const uint8_t* root_key, const uint8_t* iv, uint8_t* c,
              const faest_paramset_t* params);
uint8_t* get_vole_v_hash(vbb_t* vbb, unsigned int idx);
bf256_t* get_vole_v_prove(vbb_t* vbb, unsigned int idx);
uint8_t* get_vole_u(vbb_t* vbb);
uint8_t* get_com_hash(vbb_t* vbb);

#endif
