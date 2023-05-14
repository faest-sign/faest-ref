#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include "tree.h"
#include "aes.h"
#include "random_oracle.h"
#include "instances.h"
#include "utils.h"

typedef struct vec_com_t {
  uint8_t* h;
  uint8_t* k;
  uint32_t k_uint_size;
  uint8_t* com;
  uint32_t com_unit_size;
  uint8_t* sd;
  uint32_t sd_uint_size;
} vec_com_t;

typedef struct vec_com_rec_t {
  uint8_t* h;
  uint8_t* k; // contains vole_instances-1 elements
  uint32_t k_uint_size;
  uint8_t* com; // contains vole_instances-1 elements
  uint32_t com_unit_size;
  uint8_t* com_j;
  uint8_t* m; // contains vole_instances-1 elements
  uint32_t m_uint_size;
} vec_com_rec_t;

uint64_t NumRec(uint32_t depth, const uint8_t* bi);

void vector_commitment(const uint8_t* rootKey, const faest_paramset_t* params, vec_com_t* vecCom,
                       tree_t* tree, uint32_t voleInstances);

void vector_open(const uint8_t* k, const uint8_t* com, const uint8_t* b, uint8_t* pdec,
                 uint8_t* com_j, uint32_t numVoleInstances, uint32_t lambdabits,
                 vec_com_rec_t* vecComRec, const vec_com_t* vecCom);

void vector_reconstruction(const faest_paramset_t* params, const uint8_t* pdec,
                           const uint8_t* com_j, const uint8_t* b, uint32_t lambdaBits,
                           uint32_t NumVoleInstances, const vec_com_t* VecCom,
                           vec_com_rec_t* VecComRec);

int vector_verify(const faest_paramset_t* params, const uint8_t* pdec, const uint8_t* com_j,
                  const uint8_t* b, uint32_t lambdaBits, uint32_t numVoleInstances,
                  const vec_com_t* VecCom, vec_com_rec_t* VecComRec);