#ifndef FAEST_VC_H
#define FAEST_VC_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include "instances.h"

FAEST_BEGIN_C_DECL

// TODO: Remove the sizese,, looks reduddant !!

typedef struct vec_com_t {
  uint8_t* h;
  uint8_t* k;
  uint8_t* com;
  uint8_t* sd;
} vec_com_t;

typedef struct vec_com_rec_t {
  uint8_t* h;
  uint8_t* k;
  uint8_t* com;
  uint8_t* com_j;
  uint8_t* m;
} vec_com_rec_t;

int BitDec(uint32_t leafIndex, uint32_t depth, uint8_t* out);

uint64_t NumRec(uint32_t depth, const uint8_t* bi);

void vector_commitment(const uint8_t* rootKey, const faest_paramset_t* params, uint32_t lambda,
                       uint32_t lambdaBytes, vec_com_t* vecCom, uint32_t numVoleInstances);

void vector_open(const uint8_t* k, const uint8_t* com, const uint8_t* b, uint8_t* pdec,
                 uint8_t* com_j, uint32_t numVoleInstances, uint32_t lambdaBytes);

void vector_reconstruction(const uint8_t* pdec, const uint8_t* com_j, const uint8_t* b,
                           uint32_t lambda, uint32_t lambdaBytes, uint32_t numVoleInstances,
                           vec_com_rec_t* vecComRec);

int vector_verify(const uint8_t* pdec, const uint8_t* com_j, const uint8_t* b, uint32_t lambda,
                  uint32_t lambdaBytes, uint32_t numVoleInstances, vec_com_rec_t* vecComRec,
                  uint8_t* vecComH);

FAEST_END_C_DECL

#endif