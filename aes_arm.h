/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_AES_ARM_H
#define FAEST_AES_ARM_H

#include "aes.h"
#include "utils.h"
#include "endian_compat.h"

#if defined(_MSC_VER) && !defined(__clang__) && (defined(_M_ARM64) || defined(_M_ARM64EC))
#include <arm64_neon.h>
#else
#include <arm_neon.h>
#endif

// Repeat the final word across all four columns for SubBytes without ShiftRows.
ATTR_TARGET_ARM_AES static inline uint8x16_t arm_aes_last_word(uint8x16_t key) {
  return vreinterpretq_u8_u32(vdupq_lane_u32(vget_high_u32(vreinterpretq_u32_u8(key)), 1));
}

// Prefix XOR of four words: (a, a^b, a^b^c, a^b^c^d).
ATTR_TARGET_ARM_AES static inline uint8x16_t arm_aes_key_prefix(uint8x16_t key) {
  const uint8x16_t zero = vdupq_n_u8(0);
  key                   = veorq_u8(key, vextq_u8(zero, key, 12));
  return veorq_u8(key, vextq_u8(zero, key, 8));
}

ATTR_TARGET_ARM_AES static inline uint8x16_t arm_aes_key_core(uint8x16_t key, uint8_t rcon) {
  uint8x16_t word = arm_aes_last_word(key);
  word            = vaeseq_u8(vextq_u8(word, word, 1), vdupq_n_u8(0));
  return veorq_u8(word, vreinterpretq_u8_u32(vdupq_n_u32(htole32(rcon))));
}

ATTR_TARGET_ARM_AES static inline void arm_aes128_expand_key(uint8x16_t* rk, const uint8_t* key) {
  uint8x16_t value = rk[0] = vld1q_u8(key);
  for (unsigned int r = 1; r <= AES_ROUNDS_128; ++r) {
    value =
        veorq_u8(arm_aes_key_prefix(value), arm_aes_key_core(value, aes_round_constants[r - 1]));
    rk[r] = value;
  }
}

ATTR_TARGET_ARM_AES static inline void arm_aes192_expand_key(uint8x16_t* rk, const uint8_t* key) {
  uint8x16_t first = rk[0] = vld1q_u8(key);
  // Repeat the last two words so that the final lane holds the last valid word.
  const uint8x8_t last = vld1_u8(key + IV_SIZE);
  uint8x16_t tail      = vcombine_u8(last, last);
  unsigned int out     = 1;
  for (unsigned int r = 0; r < 8; ++r) {
    const uint8x16_t previous_tail = tail;
    first = veorq_u8(arm_aes_key_prefix(first), arm_aes_key_core(tail, aes_round_constants[r]));
    tail  = veorq_u8(arm_aes_key_prefix(tail), arm_aes_last_word(first));
    tail  = vcombine_u8(vget_low_u8(tail), vget_low_u8(tail));
    // Two six-word expansion groups supply three four-word round keys.
    if (r & 1u) {
      rk[out++] = first;
    } else {
      rk[out++] = vcombine_u8(vget_low_u8(previous_tail), vget_low_u8(first));
      rk[out++] = vcombine_u8(vget_high_u8(first), vget_low_u8(tail));
    }
  }
}

ATTR_TARGET_ARM_AES static inline void arm_aes256_expand_key(uint8x16_t* rk, const uint8_t* key) {
  uint8x16_t first = rk[0] = vld1q_u8(key);
  uint8x16_t second = rk[1] = vld1q_u8(key + IV_SIZE);
  for (unsigned int r = 2; r <= AES_ROUNDS_256; r += 2) {
    first = veorq_u8(arm_aes_key_prefix(first),
                     arm_aes_key_core(second, aes_round_constants[r / 2 - 1]));
    rk[r] = first;
    if (r < AES_ROUNDS_256) {
      const uint8x16_t subword = vaeseq_u8(arm_aes_last_word(first), vdupq_n_u8(0));
      second                   = veorq_u8(arm_aes_key_prefix(second), subword);
      rk[r + 1]                = second;
    }
  }
}

ATTR_TARGET_ARM_AES static inline void arm_aes_expand_key(uint8x16_t* rk, const uint8_t* key,
                                                          unsigned int key_words,
                                                          unsigned int rounds) {
  assert(rounds == key_words + 6);
  (void)rounds;
  switch (key_words) {
  case 4:
    arm_aes128_expand_key(rk, key);
    break;
  case 6:
    arm_aes192_expand_key(rk, key);
    break;
  default:
    assert(key_words == 8);
    arm_aes256_expand_key(rk, key);
    break;
  }
}

ATTR_TARGET_ARM_AES static inline void arm_aes_encrypt_blocks(uint8x16_t* states,
                                                              unsigned int blocks,
                                                              const uint8x16_t* rk,
                                                              unsigned int rounds) {
  // AESE adds the round key before SubBytes/ShiftRows. AESMC supplies MixColumns.
  for (unsigned int r = 0; r < rounds - 1; ++r) {
    for (unsigned int i = 0; i < blocks; ++i) {
      states[i] = vaesmcq_u8(vaeseq_u8(states[i], rk[r]));
    }
  }
  for (unsigned int i = 0; i < blocks; ++i) {
    states[i] = veorq_u8(vaeseq_u8(states[i], rk[rounds - 1]), rk[rounds]);
  }
}

#endif
