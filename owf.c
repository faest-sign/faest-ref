/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "owf.h"
#include "aes.h"

bool owf_128(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  int ret = aes128_init_round_keys(&round_keys, key);
  ret |= aes128_encrypt_block(&round_keys, input, output);
  return ret == 0;
}

bool owf_192(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  int ret = aes192_init_round_keys(&round_keys, key);
  ret |= aes192_encrypt_block(&round_keys, input, output);
  ret |= aes192_encrypt_block(&round_keys, input + 16, output + 16);
  return ret == 0;
}

bool owf_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  int ret = aes256_init_round_keys(&round_keys, key);
  ret |= aes256_encrypt_block(&round_keys, input, output);
  ret |= aes256_encrypt_block(&round_keys, input + 16, output + 16);
  return ret == 0;
}

bool owf_em_128(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  aes128_init_round_keys(&round_keys, input);
  int ret = aes128_encrypt_block(&round_keys, key, output);
  for (unsigned int i = 0; i != 16; ++i) {
    output[i] ^= key[i];
  }
  return ret == 0;
}

bool owf_em_192(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  rijndael192_init_round_keys(&round_keys, input);
  int ret = rijndael192_encrypt_block(&round_keys, key, output);
  for (unsigned int i = 0; i != 24; ++i) {
    output[i] ^= key[i];
  }
  return ret == 0;
}

bool owf_em_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  rijndael256_init_round_keys(&round_keys, input);
  int ret = rijndael256_encrypt_block(&round_keys, key, output);
  for (unsigned int i = 0; i != 32; ++i) {
    output[i] ^= key[i];
  }
  return ret == 0;
}
