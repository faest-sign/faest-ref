/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_FAEST_AES_H
#define FAEST_FAEST_AES_H

#include <stdint.h>
#include <assert.h>

#include "instances.h"
#include "aes.h"

FAEST_BEGIN_C_DECL

void aes_prove(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
               const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params);

uint8_t* aes_verify(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2, const uint8_t* chall_3,
                    const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out,
                    const faest_paramset_t* params);

FAEST_END_C_DECL

#endif
