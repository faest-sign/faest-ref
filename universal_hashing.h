/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_UNIVERSAL_HASHING_H
#define FAEST_UNIVERSAL_HASHING_H

#include <stdint.h>
#include <stddef.h>

#define UNIVERSAL_HASH_B 16

void vole_hash_128(uint8_t* h, const uint8_t* r0, const uint8_t* r1, const uint8_t* s,
                   const uint8_t* t, const uint8_t* x, size_t ell);
void vole_hash_192(uint8_t* h, const uint8_t* r0, const uint8_t* r1, const uint8_t* s,
                   const uint8_t* t, const uint8_t* x, size_t ell);
void vole_hash_256(uint8_t* h, const uint8_t* r0, const uint8_t* r1, const uint8_t* s,
                   const uint8_t* t, const uint8_t* x, size_t ell);

void zk_hash_128(uint8_t* h, const uint8_t* r, const uint8_t* s, const uint8_t* t, const bf128_t* x,
                 size_t ell);
void zk_hash_192(uint8_t* h, const uint8_t* r, const uint8_t* s, const uint8_t* t, const bf192_t* x,
                 size_t ell);
void zk_hash_256(uint8_t* h, const uint8_t* r, const uint8_t* s, const uint8_t* t, const bf256_t* x,
                 size_t ell);

#endif