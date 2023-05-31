#ifndef FAEST_H
#define FAEST_H

#include <stdlib.h>

#include "vole.h"
#include "universal_hashing.h"
#include "owf.h"
#include "faest_aes.h"

typedef struct signature_t {
  uint8_t** c;
  uint8_t* u_tilde;
  uint8_t* d;
  uint8_t* a_tilde;
  uint8_t** pdec;
  uint8_t** com_j;
  uint8_t* chall_3;
  uint8_t iv[16];
} signature_t;

signature_t init_signature(const faest_paramset_t* params);
void free_signature(signature_t sig, const faest_paramset_t* params);

void sign(const uint8_t* msg, size_t msglen, const uint8_t* sk, const uint8_t* pk,
          const uint8_t* rho, size_t rholen, const faest_paramset_t* params,
          signature_t* signature);

int verify(const uint8_t* msg, size_t msglen, const uint8_t* pk, const faest_paramset_t* params,
           const signature_t* signature);

int serialize_signature(uint8_t* dest, size_t* len, const signature_t* signature,
                        const faest_paramset_t* params);
signature_t deserialize_signature(const uint8_t* src, const faest_paramset_t* params);

#endif
