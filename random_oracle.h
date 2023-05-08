/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_RANDOM_ORACLE_H
#define FAEST_RANDOM_ORACLE_H

#include "hash_shake.h"

// implementation of H_0

typedef hash_context H0_t;

void H0_init(H0_t* ro, unsigned int security_param);
void H0_update(H0_t* ro, const uint8_t* src, size_t len);
void H0_final(H0_t* ro, uint8_t* seed, size_t seed_len, uint8_t* commitment, size_t commitment_len);

// implementation of H_1

typedef hash_context H1_t;

void H1_init(H1_t* ro, unsigned int security_param);
void H1_update(H1_t* ro, const uint8_t* src, size_t len);
void H1_final(H1_t* ro, uint8_t* digest, size_t len);

// implementation of H_2

typedef hash_context random_oracle_t;

void random_oracle_init(random_oracle_t* ro, unsigned int security_param);
void random_oracle_update(random_oracle_t* ro, const uint8_t* src, size_t len);
void random_oracle_final(random_oracle_t* ro, uint8_t* digest, size_t len);

#endif