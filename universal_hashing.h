/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_UNIVERSAL_HASHING_H
#define FAEST_UNIVERSAL_HASHING_H

#include <stdint.h>

#include "macros.h"
#include "fields.h"

FAEST_BEGIN_C_DECL

void vole_hash_128(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell);
void vole_hash_192(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell);
void vole_hash_256(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell);
void vole_hash(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell, uint32_t lambda);

void zk_hash_128(uint8_t* h, const uint8_t* sd, const bf128_t* x, unsigned int ell);
void zk_hash_192(uint8_t* h, const uint8_t* sd, const bf192_t* x, unsigned int ell);
void zk_hash_256(uint8_t* h, const uint8_t* sd, const bf256_t* x, unsigned int ell);

FAEST_END_C_DECL

#endif
