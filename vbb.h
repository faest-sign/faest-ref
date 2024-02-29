#ifndef VBB_H
#define VBB_H

#include <stdint.h>
#include <stdbool.h>

#include "vole.h"

typedef struct vbb_t vbb_t;
#include "fields.h"
#include "faest_aes.h"

struct vbb_t {
  unsigned int row_count;
  unsigned int column_count;
  unsigned int cache_idx;
  uint8_t* vole_V_cache;
  uint8_t* vole_U;
  uint8_t* com_hash;
  const uint8_t* root_key;
  const faest_paramset_t* params;
  const uint8_t* iv;
  // Optimized parameters
  bool full_size;
  uint8_t* v_buf;
  // QBB fields
  const uint8_t* c;
  uint8_t** Dtilde;
  uint8_t* vole_Q_cache;
  uint8_t** vole_Q_cache_index;
};

void init_vbb(vbb_t* vbb, unsigned int len, const uint8_t* root_key, const uint8_t* iv, uint8_t* c,
              const faest_paramset_t* params);
void clean_vbb(vbb_t* vbb);
void prepare_hash(vbb_t* vbb);
void prepare_prove(vbb_t* vbb);
uint8_t* get_vole_v_hash(vbb_t* vbb, unsigned int idx);
bf256_t* get_vole_v_prove_256(vbb_t* vbb, unsigned int idx);
bf192_t* get_vole_v_prove_192(vbb_t* vbb, unsigned int idx);
bf128_t* get_vole_v_prove_128(vbb_t* vbb, unsigned int idx);
uint8_t* get_vole_u(vbb_t* vbb);
uint8_t* get_com_hash(vbb_t* vbb);
void vector_open_ondemand(vbb_t* vbb, unsigned int idx, const uint8_t* s_, uint8_t* sig_pdec,
                          uint8_t* sig_com, unsigned int depth);

// QBB signatures
void init_qbb(vbb_t* vbb, unsigned int len, const uint8_t* iv, uint8_t* c, uint8_t* pdec_sig,
              uint8_t* com_sig, uint8_t* chall3, uint8_t* u_tilde, const faest_paramset_t* params,
              const uint8_t* sig);
uint8_t* get_vole_q_hash(vbb_t* vbb, unsigned int idx);
void prepare_verify(vbb_t* vbb, const uint8_t* sig_d, const uint8_t* sig_chall_3);
/*
bf128_t* get_vole_q_verify_128(vbb_t* vbb, unsigned int idx);
bf192_t* get_vole_q_verify_192(vbb_t* vbb, unsigned int idx);
bf256_t* get_vole_q_prove_256(vbb_t* vbb, unsigned int idx);
*/

#endif
