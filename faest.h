/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_H
#define FAEST_H

#include <stdint.h>
#include <stddef.h>

#include "instances.h"

typedef struct signature_s {
  uint8_t* c[MAX_TAU - 1];
  uint8_t* u_tilde;
  uint8_t* d;
  uint8_t* a_tilde;
  uint8_t* pdec[MAX_TAU];
  uint8_t* com_j[MAX_TAU];
  uint8_t* chall_3;
  uint8_t iv[IV_SIZE];
} signature_t;

typedef struct deserialized_signature_s {
  const uint8_t* c[MAX_TAU - 1];
  const uint8_t* u_tilde;
  const uint8_t* d;
  const uint8_t* a_tilde;
  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com_j[MAX_TAU];
  const uint8_t* chall_3;
  const uint8_t* iv;
} deserialized_signature_t;

signature_t init_signature(const faest_paramset_t* params);
void free_signature(signature_t sig, const faest_paramset_t* params);

void sign(const uint8_t* msg, size_t msglen, const uint8_t* owf_key, const uint8_t* owf_input,
          const uint8_t* owf_output, const uint8_t* rho, size_t rholen,
          const faest_paramset_t* params, signature_t* signature);

int verify(const uint8_t* msg, size_t msglen, const uint8_t* owf_input, const uint8_t* owf_output,
           const faest_paramset_t* params, const deserialized_signature_t* signature);

int serialize_signature(uint8_t* dest, size_t* len, const signature_t* signature,
                        const faest_paramset_t* params);
deserialized_signature_t deserialize_signature(const uint8_t* src, const faest_paramset_t* params);

#endif
