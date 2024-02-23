#ifndef QBB_H
#define QBB_H

#include <stdint.h>
#include <stdbool.h>

typedef struct qbb_t qbb_t;
#include "fields.h"
#include "faest_aes.h"

struct qbb_t {
  unsigned int row_count;
  unsigned int column_count;
  unsigned int cache_idx;
  uint8_t* vole_Q_cache;
  const faest_paramset_t* params;
  const uint8_t* iv;
  const uint8_t* c;
  uint8_t* com_hash;
};

void init_qbb(qbb_t* qbb, unsigned int len, const uint8_t* iv, uint8_t* c,
              uint8_t* pdec_sig, uint8_t* com_sig, uint8_t* chall3, uint8_t* u_tilde, const faest_paramset_t* params);
uint8_t* get_vole_q_hash(qbb_t* qbb, unsigned int idx);
bf256_t* get_vole_q_prove_256(qbb_t* qbb, unsigned int idx);

#endif