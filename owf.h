/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef OWF_H
#define OWF_H

#include "macros.h"

#include <stddef.h>
#include <stdint.h>

FAEST_BEGIN_C_DECL

void owf_128(const uint8_t* key, const uint8_t* input, uint8_t* output);
void owf_192(const uint8_t* key, const uint8_t* input, uint8_t* output);
void owf_256(const uint8_t* key, const uint8_t* input, uint8_t* output);

void owf_em_128(const uint8_t* key, const uint8_t* input, uint8_t* output);
void owf_em_192(const uint8_t* key, const uint8_t* input, uint8_t* output);
void owf_em_256(const uint8_t* key, const uint8_t* input, uint8_t* output);

FAEST_END_C_DECL

#endif