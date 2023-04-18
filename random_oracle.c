/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(WITH_CONFIG_H)
#include <config.h>
#endif

#include "random_oracle.h"

void random_oracle_shake128_init(random_oracle_shake128_t* ro) {
  hash_init(ro, 128);
}

void random_oracle_shake128_update(random_oracle_shake128_t* ro, const uint8_t* src, size_t len) {
  hash_update(ro, src, len);
}

void random_oracle_shake128_final(random_oracle_shake128_t* ro, uint8_t* digest, size_t len) {
  hash_final(ro);
  hash_squeeze(ro, digest, len);
}

void random_oracle_shake256_init(random_oracle_shake256_t* ro) {
  hash_init(ro, 256);
}

void random_oracle_shake256_update(random_oracle_shake256_t* ro, const uint8_t* src, size_t len) {
  hash_update(ro, src, len);
}

void random_oracle_shake256_final(random_oracle_shake256_t* ro, uint8_t* digest, size_t len) {
  hash_final(ro);
  hash_squeeze(ro, digest, len);
}
