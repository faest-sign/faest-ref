#ifndef FAEST_VC_H
#define FAEST_VC_H

#include <stdint.h>
#include <stdbool.h>

#include "instances.h"
#include "utils.h"

FAEST_BEGIN_C_DECL

typedef struct path_t {
  bool empty;
  unsigned int index;
  uint8_t* nodes;
} path_t;

typedef struct vec_com_t {
  unsigned int depth;
  uint8_t rootKey[MAX_LAMBDA_BYTES];
  path_t path;
} vec_com_t;

typedef struct vec_com_rec_t {
  unsigned int depth;
  uint8_t* nodes;
  uint8_t* b;
  uint8_t* com_j;
  path_t path;
} vec_com_rec_t;

unsigned int NumRec(unsigned int depth, const uint8_t* bi);

void get_sd_com(vec_com_t* vec_com, const uint8_t* iv, uint32_t lambda, unsigned int index,
                uint8_t* sd, uint8_t* com);
void get_sd_com_rec(vec_com_rec_t* vec_com_rec, const uint8_t* iv, uint32_t lambda,
                    unsigned int index, uint8_t* sd, uint8_t* com);

void vector_commitment(const uint8_t* rootKey, uint32_t lambda, uint32_t depth, uint8_t* path_nodes,
                       vec_com_t* vec_com_rec);
void vector_open(vec_com_t* vec_com_rec, const uint8_t* b, uint8_t* cop, uint8_t* com_j,
                 uint32_t depth, const uint8_t* iv, uint32_t lambda);
void vector_reconstruction(const uint8_t* cop, const uint8_t* com_j, const uint8_t* b,
                           uint32_t lambda, uint32_t depth, uint8_t* tree_nodes, vec_com_rec_t* vec_com_rec);

FAEST_END_C_DECL

#endif
