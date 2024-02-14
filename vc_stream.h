#ifndef FAEST_VC_STREAM_H
#define FAEST_VC_STREAM_H

#include <stdint.h>

#include "instances.h"
#include "utils.h"

FAEST_BEGIN_C_DECL

typedef struct stream_vec_com_t {
  uint8_t rootKey[MAX_LAMBDA_BYTES];
  unsigned int depth;
  unsigned int index;
  uint8_t *path;
} stream_vec_com_t;

//unsigned int NumRec(unsigned int depth, const uint8_t* bi);

void get_sd_com(stream_vec_com_t* sVecCom, const uint8_t* iv, uint32_t lambda, unsigned int index, uint8_t* sd, uint8_t* com);

void stream_vector_commitment(const uint8_t* rootKey, uint32_t lambda, stream_vec_com_t* sVecCom, uint32_t depth);

FAEST_END_C_DECL

#endif