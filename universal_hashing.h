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

#if defined(FAEST_TESTS)
void zk_hash_128(uint8_t* h, const uint8_t* sd, const bf128_t* x, unsigned int ell);
void zk_hash_192(uint8_t* h, const uint8_t* sd, const bf192_t* x, unsigned int ell);
void zk_hash_256(uint8_t* h, const uint8_t* sd, const bf256_t* x, unsigned int ell);
#endif

typedef struct {
  bf128_t h0;
  bf128_t h1;
  bf128_t s;
  bf64_t t;
  const uint8_t* sd;
} zk_hash_128_ctx;

void zk_hash_128_init(zk_hash_128_ctx* ctx, const uint8_t* sd);
void zk_hash_128_update(zk_hash_128_ctx* ctx, bf128_t v);
void zk_hash_128_finalize(uint8_t* h, zk_hash_128_ctx* ctx, bf128_t x1);

typedef struct {
  bf192_t h0;
  bf192_t h1;
  bf192_t s;
  bf64_t t;
  const uint8_t* sd;
} zk_hash_192_ctx;

void zk_hash_192_init(zk_hash_192_ctx* ctx, const uint8_t* sd);
void zk_hash_192_update(zk_hash_192_ctx* ctx, bf192_t v);
void zk_hash_192_finalize(uint8_t* h, zk_hash_192_ctx* ctx, bf192_t x1);

typedef struct {
  bf256_t h0;
  bf256_t h1;
  bf256_t s;
  bf64_t t;
  const uint8_t* sd;
} zk_hash_256_ctx;

void zk_hash_256_init(zk_hash_256_ctx* ctx, const uint8_t* sd);
void zk_hash_256_update(zk_hash_256_ctx* ctx, bf256_t v);
void zk_hash_256_finalize(uint8_t* h, zk_hash_256_ctx* ctx, bf256_t x1);

FAEST_END_C_DECL

#endif
