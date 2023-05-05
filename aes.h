/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_AES_H
#define FAEST_AES_H

#include <stdint.h>

#include "macros.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

FAEST_BEGIN_C_DECL

typedef uint8_t aes_round_key_t[4][4];

typedef struct {
  aes_round_key_t keys[10 + 1];
} aes128_round_keys_t;

typedef struct {
  aes_round_key_t keys[12 + 1];
} aes192_round_keys_t;

typedef struct {
  aes_round_key_t keys[14 + 1];
} aes256_round_keys_t;

void aes128_init_round_keys(aes128_round_keys_t* round_key, const uint8_t* key);
void aes192_init_round_keys(aes192_round_keys_t* round_key, const uint8_t* key);
void aes256_init_round_keys(aes256_round_keys_t* round_key, const uint8_t* key);

void aes128_ctr_encrypt(const aes128_round_keys_t* round_key, const uint8_t* iv,
                        const uint8_t* plaintext, uint8_t* ciphertext);
void aes192_ctr_encrypt(const aes192_round_keys_t* round_key, const uint8_t* iv,
                        const uint8_t* plaintext, uint8_t* ciphertext);
void aes256_ctr_encrypt(const aes256_round_keys_t* round_key, const uint8_t* iv,
                        const uint8_t* plaintext, uint8_t* ciphertext);

void aes_increment_ctr(uint8_t* ctr);

FAEST_END_C_DECL

#endif