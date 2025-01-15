#ifndef FAEST_VC_H
#define FAEST_VC_H

#include <stdint.h>

#include "instances.h"
#include "utils.h"

FAEST_BEGIN_C_DECL

typedef struct vec_com_t {
  uint8_t* h;
  uint8_t* k;
  uint8_t* com;
  uint8_t* sd;
} vec_com_t;

typedef struct vec_com_rec_t {
  uint8_t* h;
  uint8_t* s;
} vec_com_rec_t;

int BitDec(unsigned int leafIndex, unsigned int depth, uint8_t* out);
unsigned int NumRec(unsigned int depth, const uint8_t* bi);

static inline ATTR_CONST unsigned int bavc_max_node_index(unsigned int i, unsigned int tau_1,
                                                          unsigned int k) {
  return (i < tau_1) ? (1 << k) : (1 << (k - 1));
}

void vector_commitment(const uint8_t* rootKey, const uint8_t* iv, const faest_paramset_t* params,
                       vec_com_t* vecCom);

bool vector_open(const vec_com_t* vc, const uint16_t* i_delta, uint8_t* decom_i,
                 const faest_paramset_t* params);

bool vector_reconstruction(const uint8_t* decom_i, const uint16_t* i_delta, const uint8_t* iv,
                           const faest_paramset_t* params, vec_com_rec_t* vecComRec);
#if defined(FAEST_TESTS)
int vector_verify(const uint8_t* decom_i, const uint16_t* i_delta, const uint8_t* iv,
                  const uint8_t* vecComH, const faest_paramset_t* params);
#endif

void vec_com_clear(vec_com_t* com);
void vec_com_rec_clear(vec_com_rec_t* rec);

FAEST_END_C_DECL

#endif
