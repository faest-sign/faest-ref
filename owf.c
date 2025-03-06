/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "owf.h"
#include "aes.h"
#include "utils.h"

#include <string.h>

#if defined(HAVE_OPENSSL)
#include <openssl/evp.h>
#include <assert.h>
#endif

#if defined(HAVE_AESNI)
#include "cpu.h"
#include "aesni.h"

ATTR_TARGET_AESNI static void owf_128_aesni(const uint8_t* key, const uint8_t* input,
                                            uint8_t* output) {
  __m128i rk[AES_ROUNDS_128 + 1];
  aes128_expand_key_aesni(rk, key);

  __m128i m = _mm_xor_si128(_mm_loadu_si128((const __m128i_u*)input), rk[0]);
  for (unsigned int round = 1; round != AES_ROUNDS_128; ++round) {
    m = _mm_aesenc_si128(m, rk[round]);
  }
  m = _mm_aesenclast_si128(m, rk[AES_ROUNDS_128]);
  _mm_storeu_si128((__m128i_u*)output, m);
}

ATTR_TARGET_AESNI static void owf_192_aesni(const uint8_t* key, const uint8_t* input,
                                            uint8_t* output) {
  __m128i rk[AES_ROUNDS_192 + 1];
  aes192_expand_key_aesni(rk, key);

  __m128i temp[2];
  temp[1] = _mm_loadu_si128((const __m128i_u*)input);
  temp[0] = _mm_xor_si128(temp[1], rk[0]);
  temp[1] = _mm_xor_si128(temp[0], _mm_setr_epi32(1, 0, 0, 0));
  for (unsigned int round = 1; round != AES_ROUNDS_192; ++round) {
    temp[0] = _mm_aesenc_si128(temp[0], rk[round]);
    temp[1] = _mm_aesenc_si128(temp[1], rk[round]);
  }
  temp[0] = _mm_aesenclast_si128(temp[0], rk[AES_ROUNDS_192]);
  temp[1] = _mm_aesenclast_si128(temp[1], rk[AES_ROUNDS_192]);
  _mm_storeu_si128((__m128i_u*)output, temp[0]);
  _mm_storeu_si128((__m128i_u*)(output + IV_SIZE), temp[1]);
}

ATTR_TARGET_AESNI static void owf_256_aesni(const uint8_t* key, const uint8_t* input,
                                            uint8_t* output) {
  __m128i rk[AES_ROUNDS_256 + 1];
  aes256_expand_key_aesni(rk, key);

  __m128i temp[2];
  temp[1] = _mm_loadu_si128((const __m128i_u*)input);
  temp[0] = _mm_xor_si128(temp[1], rk[0]);
  temp[1] = _mm_xor_si128(temp[0], _mm_setr_epi32(1, 0, 0, 0));
  for (unsigned int round = 1; round != AES_ROUNDS_256; ++round) {
    temp[0] = _mm_aesenc_si128(temp[0], rk[round]);
    temp[1] = _mm_aesenc_si128(temp[1], rk[round]);
  }
  temp[0] = _mm_aesenclast_si128(temp[0], rk[AES_ROUNDS_256]);
  temp[1] = _mm_aesenclast_si128(temp[1], rk[AES_ROUNDS_256]);
  _mm_storeu_si128((__m128i_u*)output, temp[0]);
  _mm_storeu_si128((__m128i_u*)(output + IV_SIZE), temp[1]);
}

ATTR_TARGET_AESNI static void owf_em_128_aesni(const uint8_t* key, const uint8_t* input,
                                               uint8_t* output) {
  __m128i rk[AES_ROUNDS_128 + 1];
  aes128_expand_key_aesni(rk, input);

  __m128i mkey = _mm_loadu_si128((const __m128i_u*)key);
  __m128i m    = _mm_xor_si128(mkey, rk[0]);
  for (unsigned int round = 1; round != AES_ROUNDS_128; ++round) {
    m = _mm_aesenc_si128(m, rk[round]);
  }
  m = _mm_aesenclast_si128(m, rk[AES_ROUNDS_128]);
  m = _mm_xor_si128(m, mkey);
  _mm_storeu_si128((__m128i_u*)output, m);
}

#if defined(HAVE_AVX2)
ATTR_TARGET_AESNI_AVX2 static void owf_128_aesni_avx2(const uint8_t* key, const uint8_t* input,
                                                      uint8_t* output) {
  __m128i rk[AES_ROUNDS_128 + 1];
  aes128_expand_key_aesni_avx2(rk, key);

  __m128i m = _mm_xor_si128(_mm_loadu_si128((const __m128i_u*)input), rk[0]);
  for (unsigned int round = 1; round != AES_ROUNDS_128; ++round) {
    m = _mm_aesenc_si128(m, rk[round]);
  }
  m = _mm_aesenclast_si128(m, rk[AES_ROUNDS_128]);
  _mm_storeu_si128((__m128i_u*)output, m);
}

ATTR_TARGET_AESNI_AVX2 static void owf_192_aesni_avx2(const uint8_t* key, const uint8_t* input,
                                                      uint8_t* output) {
  __m128i rk[AES_ROUNDS_192 + 1];
  aes192_expand_key_aesni_avx2(rk, key);

  __m128i temp[2];
  temp[1] = _mm_loadu_si128((const __m128i_u*)input);
  temp[0] = _mm_xor_si128(temp[1], rk[0]);
  temp[1] = _mm_xor_si128(temp[0], _mm_setr_epi32(1, 0, 0, 0));
  for (unsigned int round = 1; round != AES_ROUNDS_192; ++round) {
    temp[0] = _mm_aesenc_si128(temp[0], rk[round]);
    temp[1] = _mm_aesenc_si128(temp[1], rk[round]);
  }
  temp[0] = _mm_aesenclast_si128(temp[0], rk[AES_ROUNDS_192]);
  temp[1] = _mm_aesenclast_si128(temp[1], rk[AES_ROUNDS_192]);
  _mm_storeu_si128((__m128i_u*)output, temp[0]);
  _mm_storeu_si128((__m128i_u*)(output + IV_SIZE), temp[1]);
}

ATTR_TARGET_AESNI_AVX2 static void owf_256_aesni_avx2(const uint8_t* key, const uint8_t* input,
                                                      uint8_t* output) {
  __m128i rk[AES_ROUNDS_256 + 1];
  aes256_expand_key_aesni_avx2(rk, key);

  __m128i temp[2];
  temp[1] = _mm_loadu_si128((const __m128i_u*)input);
  temp[0] = _mm_xor_si128(temp[1], rk[0]);
  temp[1] = _mm_xor_si128(temp[0], _mm_setr_epi32(1, 0, 0, 0));
  for (unsigned int round = 1; round != AES_ROUNDS_256; ++round) {
    temp[0] = _mm_aesenc_si128(temp[0], rk[round]);
    temp[1] = _mm_aesenc_si128(temp[1], rk[round]);
  }
  temp[0] = _mm_aesenclast_si128(temp[0], rk[AES_ROUNDS_256]);
  temp[1] = _mm_aesenclast_si128(temp[1], rk[AES_ROUNDS_256]);
  _mm_storeu_si128((__m128i_u*)output, temp[0]);
  _mm_storeu_si128((__m128i_u*)(output + IV_SIZE), temp[1]);
}

ATTR_TARGET_AESNI_AVX2 static void owf_em_128_aesni_avx2(const uint8_t* key, const uint8_t* input,
                                                         uint8_t* output) {
  __m128i rk[AES_ROUNDS_128 + 1];
  aes128_expand_key_aesni_avx2(rk, input);

  __m128i mkey = _mm_loadu_si128((const __m128i_u*)key);
  __m128i m    = _mm_xor_si128(mkey, rk[0]);
  for (unsigned int round = 1; round != AES_ROUNDS_128; ++round) {
    m = _mm_aesenc_si128(m, rk[round]);
  }
  m = _mm_aesenclast_si128(m, rk[AES_ROUNDS_128]);
  m = _mm_xor_si128(m, mkey);
  _mm_storeu_si128((__m128i_u*)output, m);
}
#endif
#endif

void owf_128(const uint8_t* key, const uint8_t* input, uint8_t* output) {
#if defined(HAVE_AESNI)
#if defined(HAVE_AVX2)
  if (CPU_SUPPORTS_AESNI_AVX2) {
    owf_128_aesni_avx2(key, input, output);
    return;
  }
#endif
  if (CPU_SUPPORTS_AESNI) {
    owf_128_aesni(key, input, output);
    return;
  }
#endif

#if defined(HAVE_OPENSSL)
  const EVP_CIPHER* cipher = EVP_aes_128_ecb();
  assert(cipher);
  EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
  assert(ctx);

  EVP_EncryptInit_ex(ctx, cipher, NULL, key, NULL);
  int len = 0;
  EVP_EncryptUpdate(ctx, output, &len, input, IV_SIZE);
  assert((unsigned int)len == IV_SIZE);
  EVP_CIPHER_CTX_free(ctx);
#else
  aes_round_keys_t round_keys;
  aes128_init_round_keys(&round_keys, key);
  aes128_encrypt_block(&round_keys, input, output);
#endif
}

void owf_192(const uint8_t* key, const uint8_t* input, uint8_t* output) {
#if defined(HAVE_AESNI)
#if defined(HAVE_AVX2)
  if (CPU_SUPPORTS_AESNI_AVX2) {
    owf_192_aesni_avx2(key, input, output);
    return;
  }
#endif
  if (CPU_SUPPORTS_AESNI) {
    owf_192_aesni(key, input, output);
    return;
  }
#endif

#if defined(HAVE_OPENSSL)
  const EVP_CIPHER* cipher = EVP_aes_192_ecb();
  assert(cipher);
  EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
  assert(ctx);

  EVP_EncryptInit_ex(ctx, cipher, NULL, key, NULL);
  int len = 0;
  EVP_EncryptUpdate(ctx, output, &len, input, IV_SIZE);
  assert((unsigned int)len == IV_SIZE);
  uint8_t buf[16];
  memcpy(buf, input, sizeof(buf));
  buf[0] ^= 0x1;
  EVP_EncryptUpdate(ctx, output + IV_SIZE, &len, buf, IV_SIZE);
  assert((unsigned int)len == IV_SIZE);
  EVP_CIPHER_CTX_free(ctx);
#else
  aes_round_keys_t round_keys;
  aes192_init_round_keys(&round_keys, key);

  // first block
  aes192_encrypt_block(&round_keys, input, output);
  // second block
  uint8_t buf[16];
  memcpy(buf, input, sizeof(buf));
  buf[0] ^= 0x1;
  aes192_encrypt_block(&round_keys, buf, output + 16);
#endif
}

void owf_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
#if defined(HAVE_AESNI)
#if defined(HAVE_AVX2)
  if (CPU_SUPPORTS_AESNI_AVX2) {
    owf_256_aesni_avx2(key, input, output);
    return;
  }
#endif
  if (CPU_SUPPORTS_AESNI) {
    owf_256_aesni(key, input, output);
    return;
  }
#endif

#if defined(HAVE_OPENSSL)
  const EVP_CIPHER* cipher = EVP_aes_256_ecb();
  assert(cipher);
  EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
  assert(ctx);

  EVP_EncryptInit_ex(ctx, cipher, NULL, key, NULL);
  int len = 0;
  // first block
  EVP_EncryptUpdate(ctx, output, &len, input, IV_SIZE);
  assert((unsigned int)len == IV_SIZE);
  // second block
  uint8_t buf[16];
  memcpy(buf, input, sizeof(buf));
  buf[0] ^= 0x1;
  EVP_EncryptUpdate(ctx, output + IV_SIZE, &len, buf, IV_SIZE);
  assert((unsigned int)len == IV_SIZE);
  EVP_CIPHER_CTX_free(ctx);
#else
  aes_round_keys_t round_keys;
  aes256_init_round_keys(&round_keys, key);

  // first block
  aes256_encrypt_block(&round_keys, input, output);
  // second block
  uint8_t buf[16];
  memcpy(buf, input, sizeof(buf));
  buf[0] ^= 0x1;
  aes256_encrypt_block(&round_keys, buf, output + 16);
#endif
}

void owf_em_128(const uint8_t* key, const uint8_t* input, uint8_t* output) {
#if defined(HAVE_AESNI)
#if defined(HAVE_AVX2)
  if (CPU_SUPPORTS_AESNI_AVX2) {
    owf_em_128_aesni_avx2(key, input, output);
    return;
  }
#endif
  if (CPU_SUPPORTS_AESNI) {
    owf_em_128_aesni(key, input, output);
    return;
  }
#endif

  // same as owf_128 with swapped keys and the additional xor
  owf_128(input, key, output);
  xor_u8_array(output, key, output, 16);
}

void owf_em_192(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  rijndael192_init_round_keys(&round_keys, input);
  rijndael192_encrypt_block(&round_keys, key, output);
  xor_u8_array(output, key, output, 24);
}

void owf_em_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  rijndael256_init_round_keys(&round_keys, input);
  rijndael256_encrypt_block(&round_keys, key, output);
  xor_u8_array(output, key, output, 32);
}
