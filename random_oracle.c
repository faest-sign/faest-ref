/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "random_oracle.h"

static const uint8_t domain_sep_H0 = 0;
static const uint8_t domain_sep_H1 = 1;
static const uint8_t domain_sep_H2 = 2;
static const uint8_t domain_sep_H3 = 3;

// H_0
void H0_init(H0_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H0_update(H0_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H0_final(H0_context_t* ctx, uint8_t* seed, size_t seed_len, uint8_t* commitment,
              size_t commitment_len) {
  hash_update(ctx, &domain_sep_H0, sizeof(domain_sep_H0));
  hash_final(ctx);
  hash_squeeze(ctx, seed, seed_len);
  hash_squeeze(ctx, commitment, commitment_len);
}

// H_1
void H1_init(H1_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H1_update(H1_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H1_final(H1_context_t* ctx, uint8_t* digest, size_t len) {
  hash_update(ctx, &domain_sep_H1, sizeof(domain_sep_H1));
  hash_final(ctx);
  hash_squeeze(ctx, digest, len);
}

// H_2
void H2_init(H2_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H2_update(H2_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H2_final(H2_context_t* ctx, uint8_t* digest, size_t len) {
  hash_update(ctx, &domain_sep_H2, sizeof(domain_sep_H2));
  hash_final(ctx);
  hash_squeeze(ctx, digest, len);
}

// H_3
void H3_init(H3_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H3_update(H3_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H3_final(H3_context_t* ctx, uint8_t* digest, size_t len, uint8_t* iv) {
  hash_update(ctx, &domain_sep_H3, sizeof(domain_sep_H3));
  hash_final(ctx);
  hash_squeeze(ctx, digest, len);
  hash_squeeze(ctx, iv, 16);
}
