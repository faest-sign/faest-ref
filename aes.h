/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_AES_H
#define FAEST_AES_H

#include "macros.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

FAEST_BEGIN_C_DECL

#define AES_MAX_ROUNDS 14

typedef uint8_t aes_word_t[4];
// round key with 4 (AES) up to 8 (Rijndael-256) units
typedef aes_word_t aes_round_key_t[8];

typedef struct {
  aes_round_key_t round_keys[AES_MAX_ROUNDS + 1];
} aes_round_keys_t;

bool aes128_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);
bool aes192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);
bool aes256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);
bool rijndael192_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);
bool rijndael256_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);

void aes128_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                          uint8_t* ciphertext);
void aes192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                          uint8_t* ciphertext);
void aes256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                          uint8_t* ciphertext);
void rijndael192_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                               uint8_t* ciphertext);
void rijndael256_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                               uint8_t* ciphertext);

void aes128_ctr_encrypt(const aes_round_keys_t* key, const uint8_t* iv, const uint8_t* plaintext,
                        uint8_t* ciphertext);
void aes192_ctr_encrypt(const aes_round_keys_t* key, const uint8_t* iv, const uint8_t* plaintext,
                        uint8_t* ciphertext);
void aes256_ctr_encrypt(const aes_round_keys_t* key, const uint8_t* iv, const uint8_t* plaintext,
                        uint8_t* ciphertext);

void aes_increment_iv(uint8_t* iv);

void aes_prg(const uint8_t* key, uint8_t* iv, uint8_t* out, uint16_t seclvl, uint64_t outSizeBits);

FAEST_END_C_DECL

#endif
