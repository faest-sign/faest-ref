/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "owf.h"
#include "aes.h"

#include <string.h>

#if defined(HAVE_OPENSSL)
#include <openssl/evp.h>
#include <assert.h>
#endif

void owf_128(const uint8_t* key, const uint8_t* input, uint8_t* output) {
#if defined(HAVE_OPENSSL)
  const EVP_CIPHER* cipher = EVP_aes_128_ecb();
  assert(cipher);
  EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
  assert(ctx);
  EVP_CIPHER_CTX_set_padding(ctx, 0);

  EVP_EncryptInit_ex(ctx, cipher, NULL, key, NULL);
  int len = 0;
  EVP_EncryptUpdate(ctx, output, &len, input, IV_SIZE);
  EVP_EncryptFinal_ex(ctx, output + IV_SIZE, &len);
  EVP_CIPHER_CTX_free(ctx);
#else
  aes_round_keys_t round_keys;
  aes128_init_round_keys(&round_keys, key);
  aes128_encrypt_block(&round_keys, input, output);
#endif
}

void owf_192(const uint8_t* key, const uint8_t* input, uint8_t* output) {
#if defined(HAVE_OPENSSL)
  const EVP_CIPHER* cipher = EVP_aes_192_ecb();
  assert(cipher);
  EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
  assert(ctx);
  EVP_CIPHER_CTX_set_padding(ctx, 0);

  EVP_EncryptInit_ex(ctx, cipher, NULL, key, NULL);
  int len = 0;
  EVP_EncryptUpdate(ctx, output, &len, input, IV_SIZE);
  uint8_t buf[16];
  memcpy(buf, input, sizeof(buf));
  buf[0] ^= 0x1;
  EVP_EncryptUpdate(ctx, output + IV_SIZE, &len, buf, IV_SIZE);
  EVP_EncryptFinal_ex(ctx, output + 2 * IV_SIZE, &len);
  EVP_CIPHER_CTX_free(ctx);
#else
  aes_round_keys_t round_keys;
  aes192_init_round_keys(&round_keys, key);
  aes192_encrypt_block(&round_keys, input, output);

  uint8_t buf[16];
  memcpy(buf, input, sizeof(buf));
  buf[0] ^= 0x1;
  aes192_encrypt_block(&round_keys, buf, output + 16);
#endif
}

void owf_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
#if defined(HAVE_OPENSSL)
  const EVP_CIPHER* cipher = EVP_aes_256_ecb();
  assert(cipher);
  EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
  assert(ctx);
  EVP_CIPHER_CTX_set_padding(ctx, 0);

  EVP_EncryptInit_ex(ctx, cipher, NULL, key, NULL);
  int len = 0;
  EVP_EncryptUpdate(ctx, output, &len, input, IV_SIZE);
  uint8_t buf[16];
  memcpy(buf, input, sizeof(buf));
  buf[0] ^= 0x1;
  EVP_EncryptUpdate(ctx, output + IV_SIZE, &len, buf, IV_SIZE);
  EVP_EncryptFinal_ex(ctx, output + 2 * IV_SIZE, &len);
  EVP_CIPHER_CTX_free(ctx);
#else
  aes_round_keys_t round_keys;
  aes256_init_round_keys(&round_keys, key);
  aes256_encrypt_block(&round_keys, input, output);

  uint8_t buf[16];
  memcpy(buf, input, sizeof(buf));
  buf[0] ^= 0x1;
  aes256_encrypt_block(&round_keys, buf, output + 16);
#endif
}

void owf_em_128(const uint8_t* key, const uint8_t* input, uint8_t* output) {
#if defined(HAVE_OPENSSL)
  const EVP_CIPHER* cipher = EVP_aes_128_ecb();
  assert(cipher);
  EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
  assert(ctx);
  EVP_CIPHER_CTX_set_padding(ctx, 0);

  EVP_EncryptInit_ex(ctx, cipher, NULL, input, NULL);
  int len = 0;
  EVP_EncryptUpdate(ctx, output, &len, key, IV_SIZE);
  EVP_EncryptFinal_ex(ctx, output + IV_SIZE, &len);
  EVP_CIPHER_CTX_free(ctx);
  for (unsigned int i = 0; i != 16; ++i) {
    output[i] ^= key[i];
  }
#else
  aes_round_keys_t round_keys;
  aes128_init_round_keys(&round_keys, input);
  aes128_encrypt_block(&round_keys, key, output);
  for (unsigned int i = 0; i != 16; ++i) {
    output[i] ^= key[i];
  }
#endif
}

void owf_em_192(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  rijndael192_init_round_keys(&round_keys, input);
  rijndael192_encrypt_block(&round_keys, key, output);
  for (unsigned int i = 0; i != 24; ++i) {
    output[i] ^= key[i];
  }
}

void owf_em_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  rijndael256_init_round_keys(&round_keys, input);
  rijndael256_encrypt_block(&round_keys, key, output);
  for (unsigned int i = 0; i != 32; ++i) {
    output[i] ^= key[i];
  }
}
