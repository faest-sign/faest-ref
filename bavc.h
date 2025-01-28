#ifndef FAEST_VC_H
#define FAEST_VC_H

#include <stdint.h>

#include "instances.h"
#include "utils.h"

FAEST_BEGIN_C_DECL

typedef struct bavc_t {
  uint8_t* h;
  uint8_t* k;
  uint8_t* com;
  uint8_t* sd;
} bavc_t;

typedef struct bavc_rec_t {
  uint8_t* h;
  uint8_t* s;
} bavc_rec_t;

int BitDec(unsigned int leafIndex, unsigned int depth, uint8_t* out);
unsigned int NumRec(unsigned int depth, const uint8_t* bi);

static inline ATTR_CONST unsigned int bavc_max_node_depth(unsigned int i, unsigned int tau_1,
                                                          unsigned int k) {
  return (i < tau_1) ? k : (k - 1);
}

static inline ATTR_CONST unsigned int bavc_max_node_index(unsigned int i, unsigned int tau_1,
                                                          unsigned int k) {
  return 1 << bavc_max_node_depth(i, tau_1, k);
}

void bavc_commit(const uint8_t* rootKey, const uint8_t* iv, const faest_paramset_t* params,
                 bavc_t* vecCom);

bool bavc_open(const bavc_t* vc, const uint16_t* i_delta, uint8_t* decom_i,
               const faest_paramset_t* params);

bool bavc_reconstruct(const uint8_t* decom_i, const uint16_t* i_delta, const uint8_t* iv,
                      const faest_paramset_t* params, bavc_rec_t* vecComRec);

void bavc_clear(bavc_t* com);
void bavc_rec_clear(bavc_rec_t* rec);

FAEST_END_C_DECL

#endif
