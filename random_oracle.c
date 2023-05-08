/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(WITH_CONFIG_H)
#include <config.h>
#endif

#include "random_oracle.h"

static const uint8_t domain_sep_H0 = 0;
static const uint8_t domain_sep_H1 = 1;
static const uint8_t domain_sep_H2 = 2;

void H0_init(H0_t* ro, unsigned int security_param) {
  hash_init(ro, security_param == 128 ? 128 : 256);
}

void H0_update(H0_t* ro, const uint8_t* src, size_t len) {
  hash_update(ro, src, len);
}

void H0_final(H0_t* ro, uint8_t* seed, size_t seed_len, uint8_t* commitment,
              size_t commitment_len) {
  hash_update(ro, &domain_sep_H0, sizeof(domain_sep_H0));
  hash_final(ro);
  hash_squeeze(ro, seed, seed_len);
  hash_squeeze(ro, commitment, commitment_len);
}

void H1_init(H1_t* ro, unsigned int security_param) {
  hash_init(ro, security_param == 128 ? 128 : 256);
}

void H1_update(H1_t* ro, const uint8_t* src, size_t len) {
  hash_update(ro, src, len);
}

void H1_final(H1_t* ro, uint8_t* digest, size_t len) {
  hash_update(ro, &domain_sep_H1, sizeof(domain_sep_H1));
  hash_final(ro);
  hash_squeeze(ro, digest, len);
}

void random_oracle_init(random_oracle_t* ro, unsigned int security_param) {
  hash_init(ro, security_param == 128 ? 128 : 256);
}

void random_oracle_update(random_oracle_t* ro, const uint8_t* src, size_t len) {
  hash_update(ro, src, len);
}

void random_oracle_final(random_oracle_t* ro, uint8_t* digest, size_t len) {
  hash_update(ro, &domain_sep_H2, sizeof(domain_sep_H2));
  hash_final(ro);
  hash_squeeze(ro, digest, len);
}
