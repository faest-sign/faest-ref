#ifndef VBB_H
#define VBB_H

#include <stdint.h>
#include <stdbool.h>

#include "vole.h"

typedef struct vbb_t vbb_t;
#include "fields.h"
#include "faest_aes.h"

typedef enum party_t {
  VERIFIER,
  SIGNER,
} party_t;

struct vbb_t {
  // Signer
  uint8_t* vole_U;
  // Verifier
  const uint8_t* sig;
  uint8_t* Dtilde_buf;
  // Common
  party_t party;
  unsigned int row_count;    // Budget of RMO computation
  unsigned int column_count; // Budget of CMO computation
  unsigned int cache_idx;
  const uint8_t* root_key;
  const faest_paramset_t* params;
  const uint8_t* iv;
  uint8_t* vole_cache;
  uint8_t* com_hash;
  // Optimizing parameters
  bool full_size;
  uint8_t* v_buf;
  // Vk_cache
  uint8_t* vk_buf;
  uint8_t* vk_cache;
};

void init_vbb_sign(vbb_t* vbb, unsigned int len, const uint8_t* root_key, const uint8_t* iv,
                   uint8_t* c, const faest_paramset_t* params);
void clean_vbb(vbb_t* vbb);
void prepare_hash_sign(vbb_t* vbb);
void prepare_aes_sign(vbb_t* vbb);
const uint8_t* get_vole_v_hash(vbb_t* vbb, unsigned int idx);
const bf256_t* get_vole_aes_256(vbb_t* vbb, unsigned int idx);
const bf192_t* get_vole_aes_192(vbb_t* vbb, unsigned int idx);
const bf128_t* get_vole_aes_128(vbb_t* vbb, unsigned int idx);
const uint8_t* get_vole_u(vbb_t* vbb);
const uint8_t* get_com_hash(vbb_t* vbb);
void vector_open_ondemand(vbb_t* vbb, unsigned int idx, const uint8_t* s_, uint8_t* sig_pdec,
                          uint8_t* sig_com, unsigned int depth);

// Verifier
void init_vbb_verify(vbb_t* vbb, unsigned int len, const faest_paramset_t* params,
                     const uint8_t* sig);
void prepare_hash_verify(vbb_t* vbb);
const uint8_t* get_vole_q_hash(vbb_t* vbb, unsigned int idx);
void prepare_aes_verify(vbb_t* vbb);
const uint8_t* get_dtilde(vbb_t* vbb, unsigned int idx);

// Vk_box
const bf128_t* get_vk_128(vbb_t* vbb, unsigned int idx);
const bf192_t* get_vk_192(vbb_t* vbb, unsigned int idx);
const bf256_t* get_vk_256(vbb_t* vbb, unsigned int idx);

#endif
