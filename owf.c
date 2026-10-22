/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include <assert.h>
#include <string.h>

#include "owf.h"
#include "aes.h"
#include "utils.h"
#include "cpu.h"

#if defined(HAVE_ARM_AES)
#include "aes_arm.h"
#endif

#if defined(HAVE_AESNI)
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

#if defined(HAVE_ARM_AES)
ATTR_TARGET_ARM_AES static void owf_arm_aes(const uint8_t* key, const uint8_t* input,
                                            uint8_t* output, unsigned int key_words) {
  const unsigned int rounds = key_words + 6;
  const unsigned int blocks = key_words == 4 ? 1 : 2;
  uint8x16_t rk[AES_MAX_ROUNDS + 1];
  arm_aes_expand_key(rk, key, key_words, rounds);

  uint8x16_t states[2];
  states[0] = vld1q_u8(input);
  if (blocks == 2) {
    const uint8_t flip[IV_SIZE] = {1};
    states[1]                   = veorq_u8(states[0], vld1q_u8(flip));
  }
  arm_aes_encrypt_blocks(states, blocks, rk, rounds);
  for (unsigned int i = 0; i < blocks; ++i) {
    vst1q_u8(output + i * IV_SIZE, states[i]);
  }
}

ATTR_TARGET_ARM_AES static void owf_em_128_arm_aes(const uint8_t* key, const uint8_t* input,
                                                   uint8_t* output) {
  uint8x16_t rk[AES_ROUNDS_128 + 1];
  arm_aes_expand_key(rk, input, 4, AES_ROUNDS_128);
  const uint8x16_t mkey = vld1q_u8(key);
  uint8x16_t state      = mkey;
  arm_aes_encrypt_blocks(&state, 1, rk, AES_ROUNDS_128);
  vst1q_u8(output, veorq_u8(state, mkey));
}
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

#if defined(HAVE_ARM_AES)
  if (CPU_SUPPORTS_ARM_AES) {
    owf_arm_aes(key, input, output, 4);
    return;
  }
#endif

  generic_aes_ecb_t ctx;
  int ret = generic_aes_ecb_new(&ctx, key, 128);
  assert(ret == 0);
  (void)ret;

  ret = generic_aes_ecb_encrypt(&ctx, output, input, 1);
  assert(ret == 0);
  (void)ret;

  generic_aes_ecb_free(&ctx);
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

#if defined(HAVE_ARM_AES)
  if (CPU_SUPPORTS_ARM_AES) {
    owf_arm_aes(key, input, output, 6);
    return;
  }
#endif

  generic_aes_ecb_t ctx;
  int ret = generic_aes_ecb_new(&ctx, key, 192);
  assert(ret == 0);
  (void)ret;

  uint8_t buf[2 * IV_SIZE];
  memcpy(buf, input, IV_SIZE);
  memcpy(buf + IV_SIZE, input, IV_SIZE);
  buf[IV_SIZE] ^= 0x1;

  ret = generic_aes_ecb_encrypt(&ctx, output, buf, 2);
  assert(ret == 0);
  (void)ret;

  generic_aes_ecb_free(&ctx);
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

#if defined(HAVE_ARM_AES)
  if (CPU_SUPPORTS_ARM_AES) {
    owf_arm_aes(key, input, output, 8);
    return;
  }
#endif

  generic_aes_ecb_t ctx;
  int ret = generic_aes_ecb_new(&ctx, key, 256);
  assert(ret == 0);
  (void)ret;

  uint8_t buf[2 * IV_SIZE];
  memcpy(buf, input, IV_SIZE);
  memcpy(buf + IV_SIZE, input, IV_SIZE);
  buf[IV_SIZE] ^= 0x1;

  ret = generic_aes_ecb_encrypt(&ctx, output, buf, 2);
  assert(ret == 0);
  (void)ret;

  generic_aes_ecb_free(&ctx);
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

#if defined(HAVE_ARM_AES)
  if (CPU_SUPPORTS_ARM_AES) {
    owf_em_128_arm_aes(key, input, output);
    return;
  }
#endif

  // same as owf_128 with swapped keys and the additional xor
  owf_128(input, key, output);
  xor_u8_array_inplace(output, key, 16);
}

void owf_em_192(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  rijndael192_init_round_keys(&round_keys, input);
  rijndael192_encrypt_block(&round_keys, key, output);
  xor_u8_array_inplace(output, key, 24);
}

void owf_em_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  aes_round_keys_t round_keys;
  rijndael256_init_round_keys(&round_keys, input);
  rijndael256_encrypt_block(&round_keys, key, output);
  xor_u8_array_inplace(output, key, 32);
}
