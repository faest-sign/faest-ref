/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "owf.h"
#include "aes.h"

void owf_128(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  aes128_init_round_keys(&round_keys, key);
  aes128_encrypt_block(&round_keys, input, output);
}

void owf_192(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  aes192_init_round_keys(&round_keys, key);
  aes192_encrypt_block(&round_keys, input, output);
  aes192_encrypt_block(&round_keys, input + 16, output + 16);
}

void owf_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  aes256_init_round_keys(&round_keys, key);
  aes256_encrypt_block(&round_keys, input, output);
  aes256_encrypt_block(&round_keys, input + 16, output + 16);
}

void owf_em_128(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  aes128_init_round_keys(&round_keys, key);
  aes128_encrypt_block(&round_keys, input, output);
  for (unsigned int i = 0; i != 16; ++i) {
    output[i] ^= key[i];
  }
}

void owf_em_192(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  rijndael192_init_round_keys(&round_keys, key);
  rijndael192_encrypt_block(&round_keys, input, output);
  for (unsigned int i = 0; i != 24; ++i) {
    output[i] ^= key[i];
  }
}

void owf_em_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  rijndael256_init_round_keys(&round_keys, key);
  rijndael256_encrypt_block(&round_keys, input, output);
  for (unsigned int i = 0; i != 32; ++i) {
    output[i] ^= key[i];
  }
}
