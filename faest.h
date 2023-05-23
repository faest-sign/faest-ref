#ifndef FAEST_H
#define FAEST_H

#include <stdlib.h>

#include "vole.h"
#include "universal_hashing.h"
#include "owf.h"
#include "faest_aes.h"

typedef struct signature_t {
  uint8_t* hcom;
  uint8_t** c;
  uint8_t* u_tilde;
  uint8_t* h_v;
  uint8_t* d;
  uint8_t* a_tilde;
  uint8_t* b_tilde;
  uint8_t** pdec;
  uint8_t** com_j;
} signature_t;

void keyGen(uint32_t lambda, uint32_t lambdaBytes, uint8_t* sk, uint8_t* pk);

void sign(const uint8_t* msg, const uint8_t* sk, const uint8_t* pk, const faest_paramset_t* params,
          uint32_t l, signature_t* signature);

int verify(const uint8_t* msg, const uint8_t* pk, const faest_paramset_t* params, uint32_t l,
           const signature_t* signature);

#endif