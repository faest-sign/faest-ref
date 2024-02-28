/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_FAEST_AES_H
#define FAEST_FAEST_AES_H

#include <stdint.h>
#include <assert.h>

#include "instances.h"
#include "aes.h"
#include "vbb.h"
#include "qbb.h"

FAEST_BEGIN_C_DECL

void aes_prove(const uint8_t* w, vbb_t* vbb, const uint8_t* in,
               const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params);

uint8_t* aes_verify(const uint8_t* d, qbb_t* qbb, const uint8_t* chall_2, const uint8_t* chall_3,
                    const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out,
                    const faest_paramset_t* params);

bf256_t* column_to_row_major_and_shrink_V_256(uint8_t** v, unsigned int ell);
bf192_t* column_to_row_major_and_shrink_V_192(uint8_t** v, unsigned int ell);
bf128_t* column_to_row_major_and_shrink_V_128(uint8_t** v, unsigned int ell);

FAEST_END_C_DECL

#endif
