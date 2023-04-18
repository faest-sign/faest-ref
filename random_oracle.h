/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_RANDOM_ORACLE_H
#define FAEST_RANDOM_ORACLE_H

#include "hash_shake.h"

typedef hash_context random_oracle_shake128_t;
typedef hash_context random_oracle_shake256_t;

void random_oracle_shake128_init(random_oracle_shake128_t* ro);
void random_oracle_shake128_update(random_oracle_shake128_t* ro, const uint8_t* src, size_t len);
void random_oracle_shake128_final(random_oracle_shake128_t* ro, uint8_t* digest, size_t len);

void random_oracle_shake256_init(random_oracle_shake256_t* ro);
void random_oracle_shake256_update(random_oracle_shake256_t* ro, const uint8_t* src, size_t len);
void random_oracle_shake256_final(random_oracle_shake256_t* ro, uint8_t* digest, size_t len);

#endif