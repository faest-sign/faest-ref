/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest.h"
#include "faest_aes.h"
#include "fields.h"
#include "vole.h"
#include "universal_hashing.h"
#include "utils.h"
#include "parameters.h"

#include <string.h>
#include <stdlib.h>

// TODO: Make it central somewhere
#define ALLOW_ZERO_SBOX

static_assert(FAEST_128F_ELL == FAEST_128S_ELL, "Invalid parameters");
static_assert(FAEST_128F_LAMBDA == FAEST_128S_LAMBDA, "Invalid parameters");
static_assert(FAEST_128F_Lke == FAEST_128S_Lke, "Invalid parameters");
static_assert(FAEST_128F_Nwd == FAEST_128S_Nwd, "Invalid parameters");
static_assert(FAEST_128F_R == FAEST_128S_R, "Invalid parameters");
static_assert(FAEST_128F_Senc == FAEST_128S_Senc, "Invalid parameters");
static_assert(FAEST_128F_Ske == FAEST_128S_Ske, "Invalid parameters");

static_assert(FAEST_192F_ELL == FAEST_192S_ELL, "Invalid parameters");
static_assert(FAEST_192F_LAMBDA == FAEST_192S_LAMBDA, "Invalid parameters");
static_assert(FAEST_192F_Lke == FAEST_192S_Lke, "Invalid parameters");
static_assert(FAEST_192F_Nwd == FAEST_192S_Nwd, "Invalid parameters");
static_assert(FAEST_192F_R == FAEST_192S_R, "Invalid parameters");
static_assert(FAEST_192F_Senc == FAEST_192S_Senc, "Invalid parameters");
static_assert(FAEST_192F_Ske == FAEST_192S_Ske, "Invalid parameters");

static_assert(FAEST_256F_ELL == FAEST_256S_ELL, "Invalid parameters");
static_assert(FAEST_256F_LAMBDA == FAEST_256S_LAMBDA, "Invalid parameters");
static_assert(FAEST_256F_Lke == FAEST_256S_Lke, "Invalid parameters");
static_assert(FAEST_256F_Nwd == FAEST_256S_Nwd, "Invalid parameters");
static_assert(FAEST_256F_R == FAEST_256S_R, "Invalid parameters");
static_assert(FAEST_256F_Senc == FAEST_256S_Senc, "Invalid parameters");
static_assert(FAEST_256F_Ske == FAEST_256S_Ske, "Invalid parameters");

static_assert(FAEST_EM_128F_LAMBDA == FAEST_EM_128S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_128F_Lenc == FAEST_EM_128S_Lenc, "Invalid parameters");
static_assert(FAEST_EM_128F_Nwd == FAEST_EM_128S_Nwd, "Invalid parameters");
static_assert(FAEST_EM_128F_R == FAEST_EM_128S_R, "Invalid parameters");
static_assert(FAEST_EM_128F_Senc == FAEST_EM_128S_Senc, "Invalid parameters");
// for scan-build
static_assert(FAEST_EM_128F_LAMBDA * (FAEST_EM_128F_R + 1) / 8 ==
                  sizeof(aes_word_t) * FAEST_EM_128F_Nwd * (FAEST_EM_128F_R + 1),
              "Invalid parameters");

static_assert(FAEST_EM_192F_LAMBDA == FAEST_EM_192S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_192F_Lenc == FAEST_EM_192S_Lenc, "Invalid parameters");
static_assert(FAEST_EM_192F_Nwd == FAEST_EM_192S_Nwd, "Invalid parameters");
static_assert(FAEST_EM_192F_R == FAEST_EM_192S_R, "Invalid parameters");
static_assert(FAEST_EM_192F_Senc == FAEST_EM_192S_Senc, "Invalid parameters");
// for scan-build
static_assert(FAEST_EM_192F_LAMBDA * (FAEST_EM_192F_R + 1) / 8 ==
                  sizeof(aes_word_t) * FAEST_EM_192F_Nwd * (FAEST_EM_192F_R + 1),
              "Invalid parameters");

static_assert(FAEST_EM_256F_LAMBDA == FAEST_EM_256S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_256F_Lenc == FAEST_EM_256S_Lenc, "Invalid parameters");
static_assert(FAEST_EM_256F_Nwd == FAEST_EM_256S_Nwd, "Invalid parameters");
static_assert(FAEST_EM_256F_R == FAEST_EM_256S_R, "Invalid parameters");
static_assert(FAEST_EM_256F_Senc == FAEST_EM_256S_Senc, "Invalid parameters");
// for scan-build
static_assert(FAEST_EM_256F_LAMBDA * (FAEST_EM_256F_R + 1) / 8 ==
                  sizeof(aes_word_t) * FAEST_EM_256F_Nwd * (FAEST_EM_256F_R + 1),
              "Invalid parameters");

static const bf8_t Rcon[30] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a,
    0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91,
};

// ###########################################################################################################################################
// ##################################           LAMBDA = 128            ######################################################################
// ###########################################################################################################################################

static bf128_t* column_to_row_major_and_shrink_V_128(const uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf128_t* new_v = faest_aligned_alloc(BF128_ALIGN, (ell + FAEST_128F_LAMBDA*2) * sizeof(bf128_t));
  for (unsigned int row = 0; row != ell + FAEST_128F_LAMBDA*2; ++row) {
    uint8_t new_row[BF128_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_128F_LAMBDA; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf128_load(new_row);
  }
  return new_v;
}

#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_forward_128_1(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  bf128_t* bf_y, bf128_t* bf_y_sq)
#else
static void aes_enc_forward_128_1(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  bf128_t* bf_y)
#endif
{
  // called only with Mtag == Mkey == 0

  // Step: 2
  for (unsigned int i = 0; i < 16; i++) {
    // Step: 3, 4 (bit spliced)
    // -((1 ^ Mtag) & (1 ^ Mkey)) == 0xFF
    const uint8_t xin = in[i];
    // Step: 5
    bf_y[i] = bf128_add(bf128_byte_combine_bits(xin), bf128_byte_combine_bits(xk[i]));
    #if defined(ALLOW_ZERO_SBOX)
    uint8_t tmp;
    for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
      tmp ^= (  ( (xin >> bit_j) ^ (xk[i] >> bit_j)  ) & 1 ) << bit_j;
    }
    bf_y_sq[i] = bf128_byte_combine_bits_sq(tmp);
    #endif
  }

  const bf128_t bf_two   = bf128_byte_combine_bits(2);
  const bf128_t bf_three = bf128_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_128F_R; j++) {

    for (unsigned int c = 0; c < 4; c++) {

      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_xk_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      uint8_t x_bits[4*8];
      uint8_t x_bits_by_two[8];
      uint8_t output_bits[8];
      #endif
      

      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          x_bits[r*8 + bit_l] = (x[(ix + 8 * r) / 8] >> bit_l) & 1;
        }
        // mul by two
        x_bits_by_two[0] = x_bits[r*8 + 7];
        x_bits_by_two[1] = x_bits[r*8 + 7] ^ x_bits[r*8 + 0];
        x_bits_by_two[2] = x_bits[r*8 + 1];
        x_bits_by_two[3] = x_bits[r*8 + 7] ^ x_bits[r*8 + 2];
        x_bits_by_two[4] = x_bits[r*8 + 7] ^ x_bits[r*8 + 3];
        x_bits_by_two[5] = x_bits[r*8 + 4];
        x_bits_by_two[6] = x_bits[r*8 + 5];
        x_bits_by_two[7] = x_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_x_hat[r]  = bf128_byte_combine_bits(x[(ix + 8 * r) / 8]);
        bf_xk_hat[r] = bf128_byte_combine_bits(xk[(ik + 8 * r) / 8]);
        #endif

      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        uint8_t output_bits = 0;
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits ^= (x_bits_by_two[r*8 + bit_l] ^ 
                            (x_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (x_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (x_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (x_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ xk[(ik + 8 * r) / 8])
                                )
                              )
                            )
                          ) << bit_l;
        }
        bf_y[iy + r] = bf128_byte_combine_bits(output_bits);
        bf_y_sq[iy + r] = bf128_byte_combine_bits_sq(output_bits);
      }
      #else
      // Step : 14
      bf_y[iy + 0] = bf128_add(bf_xk_hat[0], bf128_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[2]);
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[3]);

      // Step: 15
      bf_y[iy + 1] = bf128_add(bf_xk_hat[1], bf_x_hat[0]);
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[3]);

      // Step: 16
      bf_y[iy + 2] = bf128_add(bf_xk_hat[2], bf_x_hat[0]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_x_hat[1]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf128_add(bf_xk_hat[3], bf128_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[1]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[2]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[3], bf_two));
      #endif

    }
  }
}


#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_forward_128(const bf128_t* bf_x, const bf128_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf128_t* bf_y, bf128_t* bf_y_sq)
#else
static void aes_enc_forward_128(const bf128_t* bf_x, const bf128_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf128_t* bf_y)
#endif
{
  const bf128_t bf_delta  = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t bf_factor = bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey));

  // Step: 2..4
  for (unsigned int i = 0; i < 16; i++) {
    bf128_t bf_xin[8];
    for (unsigned int j = 0; j < 8; j++) {
      bf_xin[j] = bf128_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
    }
    // Step: 5
    bf_y[i] = bf128_add(bf128_byte_combine(bf_xin), bf128_byte_combine(bf_xk + (8 * i)));
    #if defined(ALLOW_ZERO_SBOX)
    bf128_t tmp[8];
    for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
      tmp[bit_j] = bf_xin[i] ^ bf_xk[i * 8 + bit_j];
    }
    bf_y_sq[i] = bf128_byte_combine_sq(tmp);
    #endif
  }

  const bf128_t bf_two   = bf128_byte_combine_bits(2);
  const bf128_t bf_three = bf128_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_128F_R; j++) {

    for (unsigned int c = 0; c < 4; c++) {

      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_xk_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      bf128_t bf_x_bits[4*8];
      bf128_t bf_x_bits_by_two[8];
      bf128_t output_bits[8];
      #endif

      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          bf_x_bits[r*8 + bit_l] = bf_x[(ix + 8 * r) + bit_l];
        }
        // mul by two
        bf_x_bits_by_two[0] = bf_x_bits[r*8 + 7];
        bf_x_bits_by_two[1] = bf_x_bits[r*8 + 7] ^ bf_x_bits[r*8 + 0];
        bf_x_bits_by_two[2] = bf_x_bits[r*8 + 1];
        bf_x_bits_by_two[3] = bf_x_bits[r*8 + 7] ^ bf_x_bits[r*8 + 2];
        bf_x_bits_by_two[4] = bf_x_bits[r*8 + 7] ^ bf_x_bits[r*8 + 3];
        bf_x_bits_by_two[5] = bf_x_bits[r*8 + 4];
        bf_x_bits_by_two[6] = bf_x_bits[r*8 + 5];
        bf_x_bits_by_two[7] = bf_x_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_x_hat[r]  = bf128_byte_combine(bf_x + (ix + 8 * r));
        bf_xk_hat[r] = bf128_byte_combine(bf_xk + (ik + 8 * r));
        #endif

      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        bf128_t output_bits[8];
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits[bit_l] = (bf_x_bits_by_two[r*8 + bit_l] ^ 
                            (bf_x_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (bf_x_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (bf_x_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (bf_x_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ bf_xk[(ik + 8 * r) / 8])
                                )
                              )
                            )
                          );
        }
        bf_y[iy + r] = bf128_byte_combine(output_bits);
        bf_y_sq[iy + r] = bf128_byte_combine_sq(output_bits);
      }
      #else
      // Step : 14
      bf_y[iy + 0] = bf128_add(bf_xk_hat[0], bf128_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[2]);
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[3]);

      // Step: 15
      bf_y[iy + 1] = bf128_add(bf_xk_hat[1], bf_x_hat[0]);
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[3]);

      // Step: 16
      bf_y[iy + 2] = bf128_add(bf_xk_hat[2], bf_x_hat[0]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_x_hat[1]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf128_add(bf_xk_hat[3], bf128_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[1]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[2]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[3], bf_two));
      #endif
    }
  }
}


#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_backward_128_1(const uint8_t* x, const uint8_t* xk, const uint8_t* out,
                                   bf128_t* y_out, bf128_t* y_out_sq)
#else
static void aes_enc_backward_128_1(const uint8_t* x, const uint8_t* xk, const uint8_t* out,
                                   bf128_t* y_out)
#endif
{
  // called only with Mtag == Mkey == 0

  uint8_t xtilde;
  // Step:2..4
  for (unsigned int j = 0; j < FAEST_128F_R; j++) {

    for (unsigned int c = 0; c < 4; c++) {

      for (unsigned int r = 0; r < 4; r++) {

        // Step: 5..6
        unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
        if (j < (FAEST_128F_R - 1)) {
          // Step: 7
          xtilde = x[ird / 8];
        } else {
          // Step: 9..11 (bit spliced)
          // -((1 ^ Mtag) & (1 ^ Mkey)) == 0xff
          const uint8_t xout = out[(ird - 128 * (FAEST_128F_R - 1)) / 8];
          xtilde             = xout ^ xk[(128 + ird) / 8];
        }

        // Step: 12..17 (bit spliced)
        // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
        const uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^ rotr8(xtilde, 2) ^ 0x5;

        // Step: 18
        y_out[16 * j + 4 * c + r] = bf128_byte_combine_bits(ytilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[16 * j + 4 * c + r] = bf128_byte_combine_bits_sq(ytilde);
        #endif
      }
    }
  }
}


#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_backward_128(const bf128_t* bf_x, const bf128_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf128_t* y_out, bf128_t* y_out_sq)
#else
static void aes_enc_backward_128(const bf128_t* bf_x, const bf128_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf128_t* y_out)
#endif
{
  // Step: 1
  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t factor =
      bf128_mul_bit(bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  // Step: 2..4
  for (unsigned int j = 0; j < FAEST_128F_R; j++) {

    for (unsigned int c = 0; c < 4; c++) {

      for (unsigned int r = 0; r < 4; r++) {

        bf128_t bf_x_tilde[8];
        // Step: 5
        unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
        // Step: 6
        if (j < (FAEST_128F_R - 1)) {
          // Step: 7
          memcpy(bf_x_tilde, bf_x + ird, sizeof(bf_x_tilde));
        } else {
          // Step: 10
          for (unsigned int i = 0; i < 8; ++i) {
            // Step: 11
            bf128_t bf_xout =
                bf128_mul_bit(factor, get_bit(out[(ird - 128 * (FAEST_128F_R - 1)) / 8], i));
            // Step: 12
            bf_x_tilde[i] = bf128_add(bf_xout, bf_xk[128 + ird + i]);
          }
        }
        // Step: 13..17
        bf128_t bf_y_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf128_add(bf128_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                    bf_x_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

        // Step: 18
        y_out[16 * j + 4 * c + r] = bf128_byte_combine(bf_y_tilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[16 * j + 4 * c + r] = bf128_byte_combine_sq(bf_y_tilde);
        #endif
      }
    }
  }
}


// Mkey = 0, this is for the prover
static void aes_enc_constraints_Mkey_0_128(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                           const bf128_t* v, const uint8_t* k, const bf128_t* vk,
                                           zk_hash_128_ctx* a0_ctx, zk_hash_128_ctx* a1_ctx) {
  
  bf128_t s[FAEST_128F_Senc]; // inverse input
  bf128_t vs[FAEST_128F_Senc]; // inverse input tag
  #if defined(ALLOW_ZERO_SBOX)
  bf128_t s_sq[FAEST_128F_Senc]; // inverse input sq
  bf128_t vs_sq[FAEST_128F_Senc]; // input input tag sq
  aes_enc_forward_128_1(w, k, in, s, s_sq);
  aes_enc_forward_128(v, vk, in, 1, 0, NULL, vs, vs_sq);
  #else
  aes_enc_forward_128_1(w, k, in, s);
  aes_enc_forward_128(v, vk, in, 1, 0, NULL, vs);
  #endif

  bf128_t s_dash[FAEST_128F_Senc]; // inverse output
  bf128_t vs_dash[FAEST_128F_Senc]; // inverse output tag
  #if defined(ALLOW_ZERO_SBOX)
  bf128_t s_dash_sq[FAEST_128F_Senc]; // inverse output sq
  bf128_t vs_dash_sq[FAEST_128F_Senc]; // inverse output tag sq
  aes_enc_backward_128_1(w, k, out, s_dash, s_dash_sq);
  aes_enc_backward_128(v, vk, 1, 0, NULL, out, vs_dash, vs_dash_sq);
  #else
  aes_enc_backward_128_1(w, k, out, s_dash);
  aes_enc_backward_128(v, vk, 1, 0, NULL, out, vs_dash);
  #endif
  
  for (unsigned int j = 0; j < FAEST_128F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Psuedoinverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf128_t mul_tag_1 = bf128_mul(vs_sq[j], vs_dash[j]);   // mul1.mac0
    const bf128_t mul_tag_2 = bf128_mul(vs[j], vs_dash_sq[j]);   // mul2.mac0
    // add the tags with the values and multiply the results
    const bf128_t mul_val_1 = bf128_mul(                        // mul1.mac1
                        bf128_add(s_sq[j], vs_sq[j]), 
                        bf128_add(s_dash[j], vs_dash[j])
                        );
    const bf128_t mul_val_2 = bf128_mul(                        // mul2.mac1
                        bf128_add(s[j], vs[j]),
                        bf128_add(s_dash_sq[j], vs_dash_sq[j])
                        );

    // the constnat term
    zk_hash_128_update(a0_ctx, mul_tag_1);
    // the linear term
    zk_hash_128_update(a1_ctx, bf128_add(
                                        bf128_add(mul_tag_1, s[j]),   // mul1.mac0 + x.val  (because here we check x^2*y = x)
                                        bf128_add(mul_val_1, vs[j])   // mul1.mac1 + x.mac
                                        )
                      );
    // the constant term
    zk_hash_128_update(a0_ctx, mul_tag_2);
    // the linear term
    zk_hash_128_update(a1_ctx, bf128_add(
                                        bf128_add(mul_tag_2, s_dash[j]),  // mul2.mac0 + y.val  (because here we check x*y^2 = y)
                                        bf128_add(mul_val_2, vs_dash[j])  // mul2.mac1 + y.mac
                                        )
                      );    
    #else
    // Inverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf128_t mul_tag = bf128_mul(vs[j], vs_dash[j]);   // mac0
    // add the tags with the values and multiply the results
    const bf128_t mul_val = bf128_mul(bf128_add(s[j], vs[j]), bf128_add(s_dash[j], vs_dash[j]));    // mac1

    // this is the constant term
    zk_hash_128_update(a0_ctx, mul_tag);    
    // this is the linear term
    zk_hash_128_update(a1_ctx, bf128_add(
                                        bf128_add(mul_tag, bf128_one()),    // mac0 + 1     (because here we check x*y = 1)
                                        mul_val                             // (mac0 + 1) + mac1
                                        )
                      );
    #endif
  }
}

// Mkey = 1, this is for the verifier
static void aes_enc_constraints_Mkey_1_128(const uint8_t* in, const uint8_t* out, const bf128_t* q,
                                           const bf128_t* qk, const uint8_t* delta,
                                           zk_hash_128_ctx* b0_ctx) {

  // Step: 11..12
  bf128_t qs[FAEST_128F_Senc];
  bf128_t qs_dash[FAEST_128F_Senc];
  #if defined(ALLOW_ZERO_SBOX)
  bf128_t qs_sq[FAEST_128F_Senc];
  bf128_t qs_dash_sq[FAEST_128F_Senc];
  aes_enc_forward_128(q, qk, in, 0, 1, delta, qs, qs_sq);
  aes_enc_backward_128(q, qk, 0, 1, delta, out, qs_dash, qs_dash_sq);
  #else
  aes_enc_forward_128(q, qk, in, 0, 1, delta, qs);
  aes_enc_backward_128(q, qk, 0, 1, delta, out, qs_dash);
  #endif
  

  // Step: 13..14
  bf128_t delta_squared = bf128_mul(bf128_load(delta), bf128_load(delta));

  for (unsigned int j = 0; j < FAEST_128F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Quicksilver multiplication
    // multiply the tags
    const bf128_t mul_tag_1 = bf128_mul(qs_sq[j], qs_dash[j]);   // mul1.mac0
    const bf128_t mul_tag_2 = bf128_mul(qs[j], qs_dash_sq[j]);   // mul2.mac0

    // the constnat term
    zk_hash_128_update(b0_ctx, bf128_add(
                                        mul_tag_1,
                                        bf128_mul(delta_squared, qs[j])
                                        )
    );
    // the constnat term
    zk_hash_128_update(b0_ctx, bf128_add(
                                        mul_tag_2,
                                        bf128_mul(delta_squared, qs_dash[j])
                                        )
    );

    #else
    // the constnat term
    zk_hash_128_update(b0_ctx, bf128_add(
                                        bf128_mul(qs[j], qs_dash[j]),
                                        delta_squared));
    #endif
    
  }
}




static void aes_128_inverse_affine(uint8_t* y, bf128_t* y_tag, uint8_t x, bf128_t* x_tag, bool isprover, bf128_t delta) {
  
  if (isprover) {
    y[0] = (rotr8(x, 7) ^ rotr8(x, 5) ^ rotr8(x, 2)) ^ 0x05; // the compressed form
  }
  for (unsigned int i = 0; i < 8; ++i) {
    y_tag[i] = bf128_add(bf128_add(x_tag[(i + 7) % 8], x_tag[(i + 5) % 8]),
                                  x_tag[(i + 2) % 8]);
    
    // TODO: The below part is missing in the InverseAffine(), 
    // this was present in Round 1, re-check later
    if (!isprover) {
      y_tag[0] = bf128_add(y_tag[0], delta);
      y_tag[2] = bf128_add(y_tag[2], delta);
    }
  }
}

static void aes_128_keyexp_backward(uint8_t* y, bf128_t* y_tag, const uint8_t* x, const bf128_t* x_tag, uint8_t* key, bf128_t* key_tag,
                                const faest_paramset_t* params, bool isprover, bf128_t delta) {

  const unsigned int lambda = params->faest_param.lambda;
  const unsigned int Ske    = params->faest_param.Ske;

  // ::2
  uint8_t x_tilde;  // contains all (F_2)^8 bits
  bf128_t x_tilde_tag[8];
  // ::3
  unsigned int iwd   = 0;
  // ::4
  bool rmvRcon       = true;
  unsigned int ircon = 0;
  // ::5-6
  for (unsigned int j = 0; j < Ske; j++) {

    // ::7
    if (isprover) {
      x_tilde = x[j] ^ key[iwd/8 + j%4];  // for the witness
    }
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_tilde_tag[bit_i] = bf128_add(x_tag[j*8 + bit_i], key_tag[iwd + (j%4)*8 + bit_i]); // for the tags of each witness bit
    }

    // ::8-10
    if (rmvRcon == true && j % 4 == 0) {
      if (isprover) {
        // adding round constant to the witness
        x_tilde = x_tilde ^ Rcon[ircon];
      }
      // adding round constant to the tags
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // for prover, no multiplication with delta
        bf128_t rcon_tag;
        if (isprover) {
          rcon_tag = bf128_from_bit(get_bit(Rcon[ircon], bit_i));
        } 
        // for verifier, multiplication with delta
        else {
          rcon_tag = bf128_mul_bit(delta, get_bit(Rcon[ircon], bit_i));
        }
        x_tilde_tag[bit_i] = bf128_add(x_tilde_tag[bit_i], rcon_tag);

      }
      ++ircon;
    }

    // ::11
    aes_inverse_affine(y, y_tag, x_tilde, x_tilde_tag, isprover, delta);

    // ::12-16 lines only relavant for aes-128
    if (j%4 == 4) {
      iwd += 128;
    }
  }
}

static void aes_128_keyexp_forward(uint8_t* y, bf128_t* y_tag, const uint8_t* w, const bf128_t* w_tag, 
                                const faest_paramset_t* params, bool isprover, bf128_t delta) {

  unsigned int lambda = params->faest_param.lambda;
  unsigned int Nk = lambda/8;
  unsigned int R = params->faest_param.R;

  // ::1-3 copying the Nk words
  if (isprover) { // witness exists only for prover
    for (unsigned int i = 0; i < Nk; i++) {
      y[i] = w[i];  // unit8 contains 8 bits
    }
  }
  // used for tag (prover) or key (verifier) 
  for (unsigned int i = 0; i < lambda; i++) {
    y_tag[i] = w_tag[i];  // each bf128 is tag for one bit
  }
  unsigned int i_wd = Nk;

  // ::4-10
  for (unsigned int j = Nk; j < 4*(R + 1); j++) {

    // ::5
    if (j % Nk == 0) {

      // ::6
      if (isprover) {
        for (unsigned int word_idx = i_wd; word_idx < 4; word_idx++) {
          y[4*j + word_idx] = w[i_wd + word_idx];   // storing byte by byte
        }
      }
      for (unsigned int word_bit_idx = i_wd*8; word_bit_idx < 4*8; word_bit_idx++) {
        y_tag[4*8*j + word_bit_idx] = w_tag[i_wd*8 + word_bit_idx]; // storing bit tags, "bit by bit"
      }
      // ::7
      i_wd += 4;    // 32 bits -> 4 words
    // ::8
    } else {
      // ::9-10
      if (isprover) {
        for (unsigned int word_idx = 0; word_idx < 4; word_idx++) {
          y[4*j + word_idx] = y[4*(j - Nk) + word_idx] ^ y[4*(j - 1) + word_idx]; // adding bitwise
        }
      }
      for (unsigned int word_bit_idx = 0; word_bit_idx < 32; word_bit_idx++) {
        y_tag[4*8*j + word_bit_idx] = bf128_add(y_tag[4*8*(j - Nk) + word_bit_idx], y_tag[4*8*(j -1) + word_bit_idx]);  // adding the tags in F2lambda "bitwwise"
      }
    }
  }
}

static void aes_128_expkey_constraints(bf128_t* z0, bf128_t* z1, uint8_t* k, bf128_t* k_tag, const uint8_t* w, const bf128_t* w_tag, 
                                    const faest_paramset_t* params, bool isprover, bf128_t delta) {


  unsigned int Ske = params->faest_param.Ske;
  unsigned int R = params->faest_param.R;
  unsigned int lambda = params->faest_param.lambda;
  unsigned int Nk = lambda/32;

  // ::1
  aes_128_keyexp_forward(k, k_tag, w, w_tag, params, isprover, delta);
  // ::2
  uint8_t w_flat[Ske];
  bf128_t w_flat_tag[8*Ske];
  aes_128_keyexp_backward(w_flat, w_flat_tag, w, w_tag, k, k_tag, params, isprover, delta);

  // ::3-5
  unsigned int iwd = 32*(Nk - 1);  // as 1 unit8 has 8 bits
  // ::6 Used only on AES-256
  // ::7
  for (unsigned int j = 0; j < FAEST_128F_Ske / 4; j++) {
    // ::8
    bf128_t k_hat[4];     // expnaded key witness
    bf128_t w_hat[4];     // inverse output
    bf128_t k_hat_sq[4];  // expanded key witness sq
    bf128_t w_hat_sq[4];      // inverse output sq

    bf128_t k_hat_tag[4]; // expanded key witness tag
    bf128_t w_hat_tag[4]; // inverse output tag
    bf128_t k_hat_tag_sq[4];  // expanded key tag sq
    bf128_t w_hat_tag_sq[4];  // inverser output tag sq
    // ::9
    for (unsigned int r = 0; r < 4; r++) {
      // ::10
      unsigned int r_prime = (r + 3)%4;
      // ::11 used only for AES-256
      // ::12-15
      if(isprover) { // only prover has the witness
        k_hat[r_prime] = bf128_byte_combine_bits(k[(iwd + 8 * r) / 8]); // lifted key witness        
        k_hat_sq[r_prime] = bf128_byte_combine_bits_sq(k[(iwd + 8 * r) / 8]); // lifted key witness sq

        w_hat[r] = bf128_byte_combine_bits(w_flat[(32 * j + 8 * r) / 8]); // lifted output
        w_hat_sq[r] = bf128_byte_combine_bits_sq(w_flat[(32 * j + 8 * r) / 8]);  // lifted output sq
      }
      // done by both prover and verifier
      k_hat_tag[r_prime] = bf128_byte_combine(k_tag + (iwd + 8 * r)); // lifted key tag
      k_hat_tag_sq[r_prime] = bf128_byte_combine_sq(k_tag + (iwd + 8 * r)); // lifted key tag sq

      w_hat_tag[r] = bf128_byte_combine(w_flat_tag + ((32 * j + 8 * r))); // lifted output tag
      w_hat_tag_sq[r] = bf128_byte_combine_sq(w_flat_tag + (32 * j + 8 * r)); // lifted output tag sq
    }
    // ::16 used only for AES-256
    // ::17
    for (unsigned int r = 0; r < 4; r++) {
      // ::18-19
      if (isprover) { // only prover has the witness vals
        z0[8*j + 2*r] = bf128_add(bf128_mul(k_hat_sq[r], w_hat[r]), k_hat[r]);
        z0[8*j + 2*r + 1] = bf128_add(bf128_mul(k_hat[r], w_hat_sq[r]), w_hat[r]);
      }
      // for both prover and verifier
      z1[8*j + 2*r] = bf128_add(bf128_mul(k_hat_tag_sq[r], w_hat_tag[r]), k_hat_tag[r]);
      z1[8*j + 2*r + 1] = bf128_add(bf128_mul(k_hat_tag[r], w_hat_tag_sq[r]), k_hat_tag[r]);
    }
    iwd = iwd + 128;
  }

}

static aes_128_enc_constraints(bf128_t* z0, bf128_t* z1, uint8_t* owf_in, 
                                uint8_t* owf_out, uint8_t* w, 
                                bf128_t* w_tag, uint8_t* k, bf128_t* k_tag, 
                                const faest_paramset_t* params, bool isprover, bf128_t delta) {

    unsigned int Nst = 4;
    unsigned int Nstbits = 32 * Nst;
    unsigned int R = params->faest_param.R;
    unsigned int Nstbytes = Nstbits/8;

    /// ::1 AddRoundKey
    bf128_t state_bits[Nstbytes];
    bf128_t state_bits_tag[Nstbytes];
    // TODO: Unsure of the squaring part come here???, spec not defined
    for (unsigned int i = 0; i < Nstbytes; i++) {
      if (isprover) {
        state_bits[i] = bf128_add(
                        bf128_byte_combine_bits(owf_in[i]), bf128_byte_combine_bits(k[i]));
        // TODO: do we square, I guess yes ???
        // uint8_t tmp;
        // for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
        //   tmp ^= (  ( (xin >> bit_j) ^ (xk[i] >> bit_j)  ) & 1 ) << bit_j;
        // }
        // bf_y_sq[i] = bf128_byte_combine_bits_sq(tmp);
        bf128_t bf_tmp[8];
        for (unsigned int j = 0; j < 8; j++) {
          bf_tmp[j] = bf128_from_bit(get_bit(owf_in[i], j));
        }
        state_bits_tag[i] = bf128_add(    // For prover we do not multiply with delta 
                        bf128_byte_combine(bf_tmp), bf128_byte_combine(k_tag + i*8));
        // TODO: do we square, I guess yes ???
        // bf128_t tmp[8];
        // for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
        //   tmp[bit_j] = bf_xin[i] ^ bf_xk[i * 8 + bit_j];
        // }
        // bf_y_sq[i] = bf128_byte_combine_sq(tmp);       
      } else {
        bf128_t bf_tmp[8];
        for (unsigned int j = 0; j < 8; j++) {  // For verifier we multiply in with delta
          bf_tmp[j] = bf128_mul_bit(delta, get_bit(owf_in[i], j));
        }
        state_bits_tag[i] = bf128_add(
                        bf128_byte_combine(bf_tmp), bf128_byte_combine(k_tag + i*8));
      }
    }

    // TODO: here need to be careful with the 0000||4bits invnorm witness values!!

    // ::2
    for (unsigned int r = 0; r < R/2; r++) {

      // ::3-4
      bf128_t state_conj[Nstbytes];
      bf128_t state_conj_tag[Nstbytes];
      // TODO: not implemented function below
      aes_128_conjugates(state_conj, state_conj_tag, state_bits, state_bits_tag, 
                          isprover, delta);

      // ::5-6
      uint8_t n[Nstbits/2*4]; // 1 uint8 contains 4 bits of the Invnorm 0000||4bits
      bf128_t n_tag[Nstbits/2]; // tag for each bit
      for (unsigned int i = 0; i < Nstbits/2; i++) {
        if (isprover) {
          // only prover has the witness
          n[i] = w[(((3*Nstbits)/2)*r)/4];
        }
        // verifier has the tags (keys)
        n_tag[i] = w_tag[(((3*Nstbits)/2)*r)];
      }

      // ::7
      bf128_t st[Nstbytes*8];
      bf128_t st_tag[Nstbytes*8];
      for (unsigned int i = 0; i < Nstbytes; i++) {
        // ::8-9
        uint8_t y[4];
        bf128_t y_tag[4];
        // TODO: not implemented function below
        aes_128_inv_norm_to_conjugates(y, y_tag, n, n_tag, isprover, delta);

        // ::10-11
        // TODO: not implemented function below
        aes_128_inv_norm_constraints(z0 + (3*r*Nstbytes), z1 + (3*r*Nstbytes), state_conj, 
                                      state_conj_tag, y, y_tag, isprover, delta);

        // ::12
        for (unsigned int j = 0; j < 8; j++) {
          // ::13-14
          // TODO: unsure what happens here, help! T.T
        }
      }

      // ::15-16
      // TODO: not implemented function below
      bf8_t k_0[16];
      aes_128_state_to_bytes(k_0, k + (Nstbits*(2*r +1)));
      // ::17
      // TODO:

      // ::18
      for (unsigned int b = 0; b < 2; b++) {
        // TODO:
        aes_128_sbox_affine();
        aes_128_shiftrows();
        aes_128_mix_coloumns();
        aes_128_add_round_key();
      }

      if (r == R/2 - 1) {
        aes_128_add_round_key();
      } else {

      }


static void aes_128_inv_norm_to_conjugates_1(bf128_t* y, uint8_t x) {
  bf128_t beta_4        = bf128_add(bf128_get_alpha(5), bf128_get_alpha(3));
  bf128_t beta_square   = beta_4;
  bf128_t beta_square_1 = bf128_mul(beta_4, beta_4);
  bf128_t beta_cube     = bf128_mul(beta_square_1, beta_4);
  for (unsigned int i = 0; i != 4; ++i) {
    y[i]          = bf128_add(bf128_add(bf128_mul_bit(bf128_one(), get_bit(x, 0)),
                                        bf128_mul_bit(beta_square, get_bit(x, 1)))
                                  bf128_add(bf128_mul_bit(beta_square_1, get_bit(x, 2)),
                                            bf128_mul_bit(beta_cube, get_bit(x, 3))));
    beta_square   = bf128_mul(beta_square, beta_square);
    beta_square_1 = bf128_mul(beta_square_1, beta_square_1);
    beta_cube     = bf128_mul(beta_cube, beta_cube);
  }
}

static void aes_128_inv_norm_to_conjugates_lambda(bf128_t* y, const bf128_t* x) {
  bf128_t beta_4        = bf128_add(bf128_get_alpha(5), bf128_get_alpha(3));
  bf128_t beta_square   = beta_4;
  bf128_t beta_square_1 = bf128_mul(beta_4, beta_4);
  bf128_t beta_cube     = bf128_mul(beta_square_1, beta_4);
  for (unsigned int i = 0; i != 4; ++i) {
    y[i] = bf128_add(
        bf128_add(bf128_mul(bf128_one(), x[0]), bf128_mul_bit(beta_square, x[1]))
            bf128_add(bf128_mul_bit(beta_square_1, x[2]), bf128_mul_bit(beta_cube, x[3])));
    beta_square   = bf128_mul(beta_square, beta_square);
    beta_square_1 = bf128_mul(beta_square_1, beta_square_1);
    beta_cube     = bf128_mul(beta_cube, beta_cube);
  }
}

static void aes_192_inv_norm_to_conjugates_1(bf192_t* y, uint8_t x) {
  bf192_t beta_4        = bf192_add(bf192_get_alpha(5), bf192_get_alpha(3));
  bf192_t beta_square   = beta_4;
  bf192_t beta_square_1 = bf192_mul(beta_4, beta_4);
  bf192_t beta_cube     = bf192_mul(beta_square_1, beta_4);
  for (unsigned int i = 0; i != 4; ++i) {
    y[i]          = bf192_add(bf192_add(bf192_mul_bit(bf192_one(), get_bit(x, 0)),
                                        bf192_mul_bit(beta_square, get_bit(x, 1)))
                                  bf192_add(bf192_mul_bit(beta_square_1, get_bit(x, 2)),
                                            bf192_mul_bit(beta_cube, get_bit(x, 3))));
    beta_square   = bf192_mul(beta_square, beta_square);
    beta_square_1 = bf192_mul(beta_square_1, beta_square_1);
    beta_cube     = bf192_mul(beta_cube, beta_cube);
  }
}

static void aes_192_inv_norm_to_conjugates_lambda(bf192_t* y, const bf192_t* x) {
  bf192_t beta_4        = bf192_add(bf192_get_alpha(5), bf192_get_alpha(3));
  bf192_t beta_square   = beta_4;
  bf192_t beta_square_1 = bf192_mul(beta_4, beta_4);
  bf192_t beta_cube     = bf192_mul(beta_square_1, beta_4);
  for (unsigned int i = 0; i != 4; ++i) {
    y[i] = bf192_add(
        bf192_add(bf192_mul(bf192_one(), x[0]), bf192_mul_bit(beta_square, x[1]))
            bf192_add(bf192_mul_bit(beta_square_1, x[2]), bf192_mul_bit(beta_cube, x[3])));
    beta_square   = bf192_mul(beta_square, beta_square);
    beta_square_1 = bf192_mul(beta_square_1, beta_square_1);
    beta_cube     = bf192_mul(beta_cube, beta_cube);
  }
}

static void aes_256_inv_norm_to_conjugates_1(bf256_t* y, uint8_t x) {
  bf256_t beta_4        = bf256_add(bf256_get_alpha(5), bf256_get_alpha(3));
  bf256_t beta_square   = beta_4;
  bf256_t beta_square_1 = bf256_mul(beta_4, beta_4);
  bf256_t beta_cube     = bf256_mul(beta_square_1, beta_4);
  for (unsigned int i = 0; i != 4; ++i) {
    y[i]          = bf256_add(bf256_add(bf256_mul_bit(bf256_one(), get_bit(x, 0)),
                                        bf256_mul_bit(beta_square, get_bit(x, 1)))
                                  bf256_add(bf256_mul_bit(beta_square_1, get_bit(x, 2)),
                                            bf256_mul_bit(beta_cube, get_bit(x, 3))));
    beta_square   = bf256_mul(beta_square, beta_square);
    beta_square_1 = bf256_mul(beta_square_1, beta_square_1);
    beta_cube     = bf256_mul(beta_cube, beta_cube);
  }
}

static void aes_256_inv_norm_to_conjugates_lambda(bf256_t* y, const bf256_t* x) {
  bf256_t beta_4        = bf256_add(bf256_get_alpha(5), bf256_get_alpha(3));
  bf256_t beta_square   = beta_4;
  bf256_t beta_square_1 = bf256_mul(beta_4, beta_4);
  bf256_t beta_cube     = bf256_mul(beta_square_1, beta_4);
  for (unsigned int i = 0; i != 4; ++i) {
    y[i] = bf256_add(
        bf256_add(bf256_mul(bf256_one(), x[0]), bf256_mul_bit(beta_square, x[1]))
            bf256_add(bf256_mul_bit(beta_square_1, x[2]), bf256_mul_bit(beta_cube, x[3])));
    beta_square   = bf256_mul(beta_square, beta_square);
    beta_square_1 = bf256_mul(beta_square_1, beta_square_1);
    beta_cube     = bf256_mul(beta_cube, beta_cube);
  }
}

static void aes_128_f256_f2_conjugates_1(bf128_t* y, const uint8_t* state) {
  for (unsigned int i = 0; i != 4; ++i) {
    uint8_t x0 = state[i];
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf128_byte_combine_bits(x0);
      x0           = bits_square(x0);
    }
    y[i * 8 + 7] = bf128_byte_combine_bits_sq(x0);
  }
}

static void aes_128_f256_f2_conjugates_128(bf128_t* y, const bf128_t* state) {
  for (unsigned int i = 0; i != 4; ++i) {
    bf128_t x[8];
    memcpy(x, state[i * 8], sizeof(x));
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf128_byte_combine(x);
      bf128_t tmp[8];
      memcpy(tmp, x, sizeof(x));
      bf128_square_bit(x, tmp);
    }
    y[i * 8 + 7] = bf128_byte_combine_sq(x);
  }
}

static void aes_192_f256_f2_conjugates_1(bf192_t* y, const uint8_t* state) {
  for (unsigned int i = 0; i != 4; ++i) {
    uint8_t x0 = state[i];
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf192_byte_combine_bits(x0);
      x0           = bits_square(x0);
    }
    y[i * 8 + 7] = bf192_byte_combine_bits_sq(x0);
  }
}

static void aes_192_f256_f2_conjugates_192(bf192_t* y, const bf192_t* state) {
  for (unsigned int i = 0; i != 4; ++i) {
    bf192_t x[8];
    memcpy(x, state[i * 8], sizeof(x));
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf192_byte_combine(x);
      bf192_t tmp[8];
      memcpy(tmp, x, sizeof(x));
      bf192_square_bit(x, tmp);
    }
    y[i * 8 + 7] = bf192_byte_combine_sq(x);
  }
}

static void aes_256_f256_f2_conjugates_1(bf256_t* y, const uint8_t* state) {
  for (unsigned int i = 0; i != 4; ++i) {
    uint8_t x0 = state[i];
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf256_byte_combine_bits(x0);
      x0           = bits_square(x0);
    }
    y[i * 8 + 7] = bf256_byte_combine_bits_sq(x0);
  }
}

static void aes_256_f256_f2_conjugates_256(bf256_t* y, const bf256_t* state) {
  for (unsigned int i = 0; i != 4; ++i) {
    bf256_t x[8];
    memcpy(x, state[i * 8], sizeof(x));
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf256_byte_combine(x);
      bf256_t tmp[8];
      memcpy(tmp, x, sizeof(x));
      bf256_square_bit(x, tmp);
    }
    y[i * 8 + 7] = bf256_byte_combine_sq(x);
  }
}

static void aes_em_192_f256_f2_conjugates_1(bf192_t* y, const uint8_t* state) {
  for (unsigned int i = 0; i != 6; ++i) {
    uint8_t x0 = state[i];
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf192_byte_combine_bits(x0);
      x0           = bits_square(x0);
    }
    y[i * 8 + 7] = bf192_byte_combine_bits_sq(x0);
  }
}

static void aes_192_f256_f2_conjugates_192(bf192_t* y, const bf192_t* state) {
  for (unsigned int i = 0; i != 6; ++i) {
    bf192_t x[8];
    memcpy(x, state[i * 8], sizeof(x));
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf192_byte_combine(x);
      bf192_t tmp[8];
      memcpy(tmp, x, sizeof(x));
      bf192_square_bit(x, tmp);
    }
    y[i * 8 + 7] = bf192_byte_combine_sq(x);
  }
}

static void aes_em_256_f256_f2_conjugates_1(bf256_t* y, const uint8_t* state) {
  for (unsigned int i = 0; i != 8; ++i) {
    uint8_t x0 = state[i];
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf256_byte_combine_bits(x0);
      x0           = bits_square(x0);
    }
    y[i * 8 + 7] = bf256_byte_combine_bits_sq(x0);
  }
}

static void aes_em_256_f256_f2_conjugates_256(bf256_t* y, const bf256_t* state) {
  for (unsigned int i = 0; i != 8; ++i) {
    bf256_t x[8];
    memcpy(x, state[i * 8], sizeof(x));
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf256_byte_combine(x);
      bf256_t tmp[8];
      memcpy(tmp, x, sizeof(x));
      bf256_square_bit(x, tmp);
    }
    y[i * 8 + 7] = bf256_byte_combine_sq(x);
  }
}

// TODO:
void aes_128_inv_norm_constraints(bf128_t* z0, bf128_t* z1, const bf128_t* state_bits,
                                  const bf128_t* state_bits_tags, const uint8_t* y,
                                  const bf128_t* y_tag, uint8_t isprover) {}

void aes_128_state_to_bytes(bf128_t* out, bf128_t* out_tag, const uint8_t* state,
                            const bf128_t* state_tag, bool isprover) {

  for (unsigned int i = 0; i < 16; i++) {
    if (isprover) {
      out[i] = bf128_byte_combine_bits(state[i]);
    }
    out_tag[i] = bf128_byte_combine(state_tag + i * 8);
  }
}


static void aes_128_sbox_affine(bf128_t* out, bf128_t* out_tag, const uint8_t* in, const bf128_t* in_tag, bool dosq, const faest_paramset_t* params, bool isprover) {

  unsigned int Nst_bytes = params->faest_param.lambda/8;

  unsigned int s, t;
  if (dosq) {
    s = 2; t = 1;
  } else {
    s = 1; t = 0;
  }

  bf128_t C[9];
  uint8_t x[9] = {0x05,0x09,0xf9,0x25,0xf4,0x01,0xb5,0x8f,0x63};

  // ::5-6
  // TODO: we sq the bits and them do bytecombine?
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf128_byte_combine_bits_sq(x[i]);
    }
  } else {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf128_byte_combine_bits(x[i]);
    }
  }

  // TODO: this looks ugly :(
  if (dosq) {
    for (unsigned int i = 0; i < Nst_bytes; i++) {
      for (unsigned int Cidx = 0; Cidx < 9; Cidx++) {
        if (isprover) {
          // TODO: re-check for possible indexing error
          out[i] += bf128_mul(C[Cidx], bf128_byte_combine_bits_sq(in[i]));  // TODO: I think this is wrong, ask Peter
        }
        // TODO: re-check for possible indexing error
        out_tag[i] += bf128_mul(C[Cidx], bf128_byte_combine_sq(in_tag + i*8));  // TODO: I think this is wrong, ask Peter
      }
    }
  } else {
    for (unsigned int i = 0; i < Nst_bytes; i++) {
      for (unsigned int Cidx = 0; Cidx < 9; Cidx++) {
        if (isprover) {
          // TODO: re-check for possible indexing error
          out[i] += bf128_mul(C[Cidx], bf128_byte_combine_bits(in[i]));  // TODO: I think this is wrong, ask Peter
        }
        // TODO: re-check for possible indexing error
          out_tag[i] += bf128_mul(C[Cidx], bf128_byte_combine(in_tag + i*8));  // TODO: I think this is wrong, ask Peter
      }
    }
  }

}


static void aes_128_inverse_shiftrows(uint8_t* out, bf128_t* out_tag, const uint8_t* in, const bf128_t* in_tag, const faest_paramset_t* params, bool isprover) {
  unsigned int Nst = 4;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      unsigned int i;
      if (r <= 1) {
        i = 4*((c-r)%4) + r;
      } 
      else {
        i = 4*((c-r-1) % 4) + r;
      }

      if (isprover) {
        for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
          out[8*(4*c + r) + byte_idx] = in[8*i + byte_idx];
        }
      }
      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
          out_tag[8*(4*c + r) + byte_idx] = in_tag[8*i + byte_idx];
        }
    }
  }
}


static void aes_128_shiftrows(uint8_t* out, bf128_t* out_tag, const uint8_t* in, const bf128_t* in_tag, const faest_paramset_t* params, bool isprover) {
  unsigned int Nst = 4;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if (r <= 1) {
        if (isprover) {
          out[4*c + r] = in[4*((c + r) % 4) + r];
        }
        out_tag[4*c + r] = in_tag[(4*((c + r) % 4) + r)];
      } 
      else {
        if (isprover) {
          out[4*c + r] = in[4*((c + r) % 4) + r];
        }
        out_tag[4*c + r] = in_tag[(4*((c + r + 1) % 4) + r)];
      }
    }
  }
}


static void aes_128_mix_coloumns_prover(bf128_t* y, bf128_t* y_tag, const uint8_t* in, const uint8_t* in_tag, bool dosq, const faest_paramset_t* params) {
  
  unsigned int Nst = 3;
  
  //  ::2-4
  bf128_t v1 = bf128_byte_combine_bits(1);
  bf128_t v2 = bf128_byte_combine_bits(2);
  bf128_t v3 = bf128_byte_combine_bits(3);
  if (dosq) {
    v1 = bf128_mul(v1, v1);
    v2 = bf128_mul(v2, v2);
    v3 = bf128_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {

    unsigned int i0 = 4*c;
    unsigned int i1 = 4*c + 1;
    unsigned int i2 = 4*c + 2;
    unsigned int i3 = 4*c + 3;

    // ::7
    bf128_t tmp1_val = bf128_mul(bf128_byte_combine_bits_sq(in[i0]), v2);
    bf128_t tmp2_val = bf128_mul(bf128_byte_combine_bits_sq(in[i1]), v3);
    bf128_t tmp3_val = bf128_mul(bf128_byte_combine_bits_sq(in[i2]), v1);
    bf128_t tmp4_val = bf128_mul(bf128_byte_combine_bits_sq(in[i3]), v1);
    y[i0] = bf128_add(bf128_add(tmp1_val, tmp2_val), bf128_add(tmp3_val, tmp4_val));

    bf128_t tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i0]), v2);
    bf128_t tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i1]), v3);
    bf128_t tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i2]), v1);
    bf128_t tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i3]), v1);
    y_tag[i0] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    // ::8
    tmp1_val = bf128_mul(bf128_byte_combine_bits_sq(in[i0]), v1);
    tmp2_val = bf128_mul(bf128_byte_combine_bits_sq(in[i1]), v2);
    tmp3_val = bf128_mul(bf128_byte_combine_bits_sq(in[i2]), v3);
    tmp4_val = bf128_mul(bf128_byte_combine_bits_sq(in[i3]), v1);
    y[i1] = bf128_add(bf128_add(tmp1_val, tmp2_val), bf128_add(tmp3_val, tmp4_val));

    tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i1]), v2);
    tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i2]), v3);
    tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i3]), v1);
    y_tag[i1] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    // ::9
    tmp1_val = bf128_mul(bf128_byte_combine_bits_sq(in[i0]), v1);
    tmp2_val = bf128_mul(bf128_byte_combine_bits_sq(in[i1]), v1);
    tmp3_val = bf128_mul(bf128_byte_combine_bits_sq(in[i2]), v2);
    tmp4_val = bf128_mul(bf128_byte_combine_bits_sq(in[i3]), v3);
    y[i2] = bf128_add(bf128_add(tmp1_val, tmp2_val), bf128_add(tmp3_val, tmp4_val));

    tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i1]), v1);
    tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i2]), v2);
    tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i3]), v3);
    y_tag[i2] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    // ::10
    tmp1_val = bf128_mul(bf128_byte_combine_bits_sq(in[i0]), v3);
    tmp2_val = bf128_mul(bf128_byte_combine_bits_sq(in[i1]), v1);
    tmp3_val = bf128_mul(bf128_byte_combine_bits_sq(in[i2]), v1);
    tmp4_val = bf128_mul(bf128_byte_combine_bits_sq(in[i3]), v2);
    y[i3] = bf128_add(bf128_add(tmp1_val, tmp2_val), bf128_add(tmp3_val, tmp4_val));

    tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i1]), v1);
    tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i2]), v2);
    tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i3]), v3);
    y_tag[i3] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));
  
  }

}


static void aes_128_mix_coloumns_prover(bf128_t* y_key, const uint8_t* in_key, bool dosq, const faest_paramset_t* params) {
  
  unsigned int Nst = 3;
  
  //  ::2-4
  bf128_t v1 = bf128_byte_combine_bits(1);
  bf128_t v2 = bf128_byte_combine_bits(2);
  bf128_t v3 = bf128_byte_combine_bits(3);
  if (dosq) {
    v1 = bf128_mul(v1, v1);
    v2 = bf128_mul(v2, v2);
    v3 = bf128_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {

    unsigned int i0 = 4*c;
    unsigned int i1 = 4*c + 1;
    unsigned int i2 = 4*c + 2;
    unsigned int i3 = 4*c + 3;

    // ::7
      bf128_t tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i0]), v2);
    bf128_t tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i1]), v3);
    bf128_t tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i2]), v1);
    bf128_t tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i3]), v1);
    y_key[i0] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    // ::8
    tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i1]), v2);
    tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i2]), v3);
    tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i3]), v1);
    y_key[i1] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    // ::9
    tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i1]), v1);
    tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i2]), v2);
    tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i3]), v3);
    y_key[i2] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    // ::10
    tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i1]), v1);
    tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i2]), v2);
    tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_tag[i3]), v3);
    y_key[i3] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));
  
  }

}



// TODO:
static void aes_128_add_round_key() {

}


// TODO:
static void aes_128_conjugates(uint8_t* state_conj, bf128_t* state_conj_tag, const uint8_t* state, const bf128_t* state_tag, bool isprover, bf128_t delta) {

}

static void aes_128_deg2to3(bf128_t* z0, bf128_t* z1, bf128_t val, bf128_t tag, bool isprover, bf128_t delta) {
  if(isprover) {
    z0[0] = val;
    z1[0] = tag;
  } else {
    // verifier does not have tag
    z1[0] = bf128_mul(tag, delta);
  }
}

static void constant_to_vole_128(bf128_t* tag, const uint8_t* val, bool isprover, bf128_t delta) {

  if (isprover) {
    // the val stay the same as the val is a pub const!
    for (unsigned int i = 0; i < 128; i++) {
      tag[i] = bf128_zero();  // for constant values the tag is zero
    }   
  } else {
    for (unsigned int i = 0; i < 128; i++) {
      tag[i] = bf128_mul(bf128_from_bit(get_bit(val[i/8], i%8)), delta);  // multiply delta with each bit
    }  
  }
}

static void aes_constraints_128(bf128_t* z0, bf128_t* z1, const uint8_t* w, const bf128_t* w_tag, const uint8_t* owf_in, const uint8_t* owf_out, const faest_paramset_t* params, bool isprover, bf128_t delta) {

  unsigned int lambda = params->faest_param.lambda;
  unsigned int R = params->faest_param.R;
  unsigned int Ske = params->faest_param.Ske;
  unsigned int Lke = lambda + 8*Ske;
  unsigned int Lenc = params->faest_param.Lenc;
  unsigned int Senc = params->faest_param.Senc;
  // ::1-3 owf_in, owf_out, z and z_tag

  // ::4-5
  aes_deg2to3(&z0[0], &z1[0], (w[0]&1) & ((w[0]>>1)&1), bf128_mul(w_tag[0], w_tag[1]), isprover, delta);   // multiplying the first and second bit of the witness and the tags

  // jump to ::13 for AES
  bf128_t owf_in_tag[128];  // don't really use it anywhere, just sticking close to the spec! // TODO: can also remove this
  constant_to_vole_128(owf_in_tag, owf_in, true, bf128_one());

  // ::14
  bf128_t owf_out_tag[128];  // don't really use it anywhere, just sticking close to the spec! // TODO: can also remove this
  constant_to_vole_128(owf_out_tag, owf_out, true, bf128_one());

  // ::15 skiped as B = 1
  // ::16
  bf128_t z_tilde_0_expkey[FAEST_128F_Ske / 4 + 2*4];
  bf128_t z_tilde_1_expkey[FAEST_128F_Ske / 4 + 2*4];
  uint8_t k[(R+1)*lambda/8];
  bf128_t k_tag[(R+1)*lambda];
   // if isprover == true, z_tilde_0 will be empty after aes_128_expkey_constraints() returns
  aes_128_expkey_constraints(z_tilde_0_expkey, z_tilde_1_expkey, k, k_tag, w, w_tag, params, isprover, delta);
  
  // ::17
  for (unsigned int i = 0; i < FAEST_128F_Ske / 4 + 2*4; i++) {
    aes_128_deg2to3(z0 + 1+i, z1 + 1+i, z_tilde_0_expkey[i], z_tilde_1_expkey[i], isprover, delta);
  }

  // ::18 b = 0
  // ::19
  uint8_t w_tilde[Lenc/8];
  bf128_t w_tilde_tag[Lenc];
  for (unsigned int i = 0; i < Lenc/8; i++) {   // can also do a memcpy
    w_tilde[i] = w[Lke/8 + i];  // copying 8 bits at a time
  }
  for (unsigned int i = 0; i < Lenc; i++) {
    w_tilde_tag[i] = w_tag[Lke + i];  // copying 1 bit's tag at a time
  }
  // ::20 not needed for aes128
  // ::21
  bf128_t z_tilde_0[Senc];
  bf128_t z_tilde_1[Senc];

  aes_128_enc_constraints(z_tilde_0, z_tilde_1, owf_in, owf_out, w, w_tag, k, k_tag, params, isprover, delta);

  for (unsigned int i = 0; i < Senc; i++) {
    z0[1+(FAEST_128F_Ske / 4 + 2*4) + i] = z_tilde_0[i];
    z1[1+(FAEST_128F_Ske / 4 + 2*4) + i] = z_tilde_1[i];
  }

}

static void aes_prove_128(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde, const uint8_t* w, const uint8_t* u, 
                          const uint8_t** V, const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* chall_2, const faest_paramset_t* params) {

  unsigned int lambda = params->faest_param.lambda;
  unsigned int ske = params->faest_param.Ske;
  unsigned int senc = params->faest_param.Senc;
  // TODO: CHANGE THIS FOR OTHER SETTING WHEN COPY PASTING!!!!!
  unsigned int beta = 1;
  unsigned int c = 2*ske + (3/2)*senc + 1;

  // ::1-5
  // also includes the lifting of V at ::5
  bf128_t* bf_v = column_to_row_major_and_shrink_V_128(V, FAEST_128F_ELL); // This is the tag for w
  // we have w in its f2 form

  // ::6-9 embed VOLE masks
  bf128_t bf_u_star_0 = bf128_load(u);
  bf128_t bf_u_star_1 = bf128_load(u + lambda);
  bf128_t bf_v_star_0 = bf128_sum_poly(bf_v);
  bf128_t bf_v_star_1 = bf128_sum_poly(bf_v + lambda);

  // ::10-12 <z>^3_\lambda -> z \in F^{<=c}_{2^lambda}
  // (0 , z'_0 , z'_1), z'_i \in F_{2^lambda}, z'_0 takes the val, z'_1 takes the tag
  bf128_t bf_z0[c];
  bf128_t bf_z1[c];
  aes_constraints_128(bf_z0, bf_z1, w, bf_v, owf_in, owf_out, params, false, bf128_one());  // delta set to 1 for prover

  // ::13
  bf128_t a0[c*3];
  bf128_t a1[c*3];
  bf128_t a2[c*3];
  // TODO: the magical parsing, asked Peter already

  // Step: 16..18
  zk_hash_128_ctx a0_ctx;
  zk_hash_128_ctx a1_ctx;
  zk_hash_128_ctx a2_ctx;
  zk_hash_128_init(&a0_ctx, chall_2);
  zk_hash_128_init(&a1_ctx, chall_2);
  zk_hash_128_init(&a2_ctx, chall_2);

  // TODO: ugly zk update here
  for (unsigned int i = 0; i < c*3; i++) {
    zk_hash_128_update(&a0_ctx, a0[i]);
    zk_hash_128_update(&a1_ctx, a1[i]);
    zk_hash_128_update(&a2_ctx, a2[i]);
  }

  zk_hash_128_finalize(a0_tilde, &a0_ctx, bf_u_star_0);
  zk_hash_128_finalize(a1_tilde, &a1_ctx, bf_v_star_0 + bf_u_star_1);
  zk_hash_128_finalize(a2_tilde, &a2_ctx, bf_v_star_1);

}


static uint8_t* aes_verify_128(const uint8_t* d, const uint8_t** Q, const uint8_t* owf_in, const uint8_t* owf_out,
                                 const uint8_t* chall_2, const uint8_t* chall_3,  const uint8_t* a1_tilde, const uint8_t* a2_tilde, const faest_paramset_t* params) {

  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.tau0;
  const unsigned int k0  = params->faest_param.k;
  const unsigned int t1  = params->faest_param.tau1;
  const unsigned int k1  = (t0 != 0) ? k0 - 1 : k0;
  unsigned int lambda = params->faest_param.lambda;
  unsigned int ske = params->faest_param.Ske;
  unsigned int senc = params->faest_param.Senc;
  // TODO: CHANGE THIS FOR OTHER SETTING WHEN COPY PASTING!!!!!
  unsigned int beta = 1;
  unsigned int c = 2*ske + (3/2)*senc + 1;

  // ::2
  bf128_t bf_delta = bf128_load(chall_3);
  bf128_t bf_delta_sq = bf128_mul(bf_delta, bf_delta);

  // ::3-7
  bf128_t* bf_q = column_to_row_major_and_shrink_V_128(Q, FAEST_128F_ELL); // This is the vole key for masked extended witness d
  // w is "d" in F2

  // ::8-10
  bf128_t bf_q_star_0 = bf128_sum_poly(bf_q);
  bf128_t bf_q_star_1 = bf128_sum_poly(bf_q + lambda);

  // ::11
  bf128_t bf_q_star = bf128_add(bf_q_star_0, bf128_mul(bf_delta, bf_q_star_1));

  // ::12-13
  bf128_t bf_z[c*3];
  bf128_t bf_z_tag[c*3];
  aes_constraints(bf_z, bf_z_tag, d, bf_q, owf_in, owf_out, params);

  // ::14
  bf128_t b[c];
  // TODO: the magical parsing, asked Peter already

  // ::15
  zk_hash_128_ctx b_ctx;
  zk_hash_128_init(&b_ctx, chall_2);
  // TODO: ugly zk update here
  for (unsigned int i = 0; i < c; i++) {
    zk_hash_128_update(&b_ctx, b[i]);
  }
  uint8_t q_tilde[lambda/8*3];
  zk_hash_128_finalize(q_tilde, &b_ctx, bf_q_star);

  // ::16-17
  bf128_t tmp1 = bf128_mul(bf128_load(a1_tilde), bf_delta);
  bf128_t tmp2 = bf128_mul(bf128_load(a2_tilde), bf_delta_sq);
  bf128_t tmp3 = bf128_add(tmp1, tmp2);
  bf128_t ret = bf128_add(bf128_load(q_tilde), tmp3);
  uint8_t* a0_tilde = malloc(lambda/8*3);   // freed in faest_verify
  bf128_store(a0_tilde, ret);
  return a0_tilde;

}

// ###########################################################################################################################################
// ##################################           LAMBDA = 192            ######################################################################
// ###########################################################################################################################################

static bf192_t* column_to_row_major_and_shrink_V_192(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf192_t* new_v = faest_aligned_alloc(BF192_ALIGN, (ell + FAEST_192F_LAMBDA) * sizeof(bf192_t));
  for (unsigned int row = 0; row != ell + FAEST_192F_LAMBDA; ++row) {
    uint8_t new_row[BF192_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_192F_LAMBDA; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf192_load(new_row);
  }

  return new_v;
}

static void aes_key_schedule_forward_192(const bf192_t* v, bf192_t* bf_out) {
  // Step: 1 sanity check (skipped)

  memcpy(bf_out, v, FAEST_192F_LAMBDA * sizeof(bf192_t));

  // Step: 4
  unsigned int i_wd = FAEST_192F_LAMBDA;
  // Step: 5..10
  for (unsigned int j = FAEST_192F_Nwd; j < 4 * (FAEST_192F_R + 1); j++) {
    if ((j % FAEST_192F_Nwd) == 0 || (FAEST_192F_Nwd > 6 && (j % FAEST_192F_Nwd) == 4)) {
      memcpy(bf_out + j * 32, v + i_wd, sizeof(bf192_t) * 32);
      i_wd += 32;
    } else {
      for (unsigned int i = 0; i < 32; i++) {
        bf_out[(32 * j) + i] =
            bf192_add(bf_out[32 * (j - FAEST_192F_Nwd) + i], bf_out[32 * (j - 1) + i]);
      }
    }
  }
}

static void aes_key_schedule_backward_192(const bf192_t* v, const bf192_t* Vk, uint8_t Mtag,
                                          uint8_t Mkey, const uint8_t* delta, bf192_t* bf_out) {
  // Step: 1
  assert(!((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)));

  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  unsigned int iwd       = 0;
  unsigned int c         = 0;
  unsigned int ircon     = 0;

  bf192_t bf_minus_mkey       = bf192_from_bit(1 ^ Mkey);
  uint8_t minus_mtag          = 1 ^ Mtag;
  bf192_t bf_mkey_times_delta = bf192_mul_bit(bf_delta, Mkey);
  bf_mkey_times_delta         = bf192_add(bf_mkey_times_delta, bf_minus_mkey);

  for (unsigned int j = 0; j < FAEST_192F_Ske; j++) {
    // Step 7
    bf192_t bf_x_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      bf_x_tilde[i] = bf192_add(v[8 * j + i], Vk[iwd + 8 * c + i]);
    }

    if (Mtag == 0 && c == 0) {
      // Step 9
      uint8_t r = Rcon[ircon];
      ircon     = ircon + 1;

      bf192_t bf_r[8];
      for (unsigned int i = 0; i < 8; i++) {
        // Step 12
        bf_r[i] = bf192_mul_bit(bf_mkey_times_delta, get_bit(r, i));
        // Step 13
        bf_x_tilde[i] = bf192_add(bf_x_tilde[i], bf_r[i]);
      }
    }

    for (unsigned int i = 0; i < 8; ++i) {
      bf_out[i + 8 * j] = bf192_add(bf192_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                    bf_x_tilde[(i + 2) % 8]);
    }
    bf_out[0 + 8 * j] =
        bf192_add(bf_out[0 + 8 * j], bf192_mul_bit(bf_mkey_times_delta, minus_mtag));
    bf_out[2 + 8 * j] =
        bf192_add(bf_out[2 + 8 * j], bf192_mul_bit(bf_mkey_times_delta, minus_mtag));
    c = c + 1;

    if (c == 4) {
      c = 0;
      iwd += 192;
    }
  }
}

// Mkey = 0, this is for the prover
static void aes_key_schedule_constraints_Mkey_0_192(const uint8_t* w, const bf192_t* v,
                                                    zk_hash_192_ctx* a0_ctx,
                                                    zk_hash_192_ctx* a1_ctx, uint8_t* k,
                                                    bf192_t* vk, const faest_paramset_t* params) {
  // for scan-build
  assert(FAEST_192F_Ske == params->faest_param.Ske);

  // Step: 2
  aes_key_schedule_forward_1(w, k, params);

  // Step: 3
  aes_key_schedule_forward_192(v, vk);

  // Step: 4, this is for the inverse witness
  uint8_t w_dash[FAEST_192F_Ske];
  aes_key_schedule_backward_1(w + FAEST_192F_LAMBDA / 8, k, w_dash, params);

  // Step: 5
  bf192_t v_w_dash[FAEST_192F_Ske * 8];
  aes_key_schedule_backward_192(v + FAEST_192F_LAMBDA, vk, 1, 0, NULL, v_w_dash);

  // Step: 6..8
  unsigned int iwd = 32 * (FAEST_192F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_192F_Ske / 4; j++) {
    bf192_t bf_k_hat[4];
    bf192_t bf_v_k_hat[4];
    bf192_t bf_w_hat[4];
    bf192_t bf_v_w_hat[4];
    #if defined(ALLOW_ZERO_SBOX)
    bf192_t bf_k_hat_sq[4];
    bf192_t bf_v_k_hat_sq[4];
    bf192_t bf_w_hat_sq[4];
    bf192_t bf_v_w_hat_sq[4];
    #endif

    for (unsigned int r = 0; r < 4; r++) {
      // Step: 10..16
      bf_k_hat[(r + 3) % 4]   = bf192_byte_combine_bits(k[(iwd + 8 * r) / 8]); // lifted inverse input
      bf_v_k_hat[(r + 3) % 4] = bf192_byte_combine(vk + (iwd + 8 * r));         // lifted inverse input tag
      bf_w_hat[r]        = bf192_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]); // lifted inverse output
      bf_v_w_hat[r]      = bf192_byte_combine(v_w_dash + (32 * j + 8 * r));       // lifted inverse output tag

      // squaring the bits and macs and lifting
      #if defined(ALLOW_ZERO_SBOX)
      bf_k_hat_sq[(r + 3) % 4] = bf192_byte_combine_bits_sq(k[(iwd + 8 * r) / 8]);    // lifted inverse input sq
      bf_v_k_hat_sq[(r + 3) % 4] = bf192_byte_combine_sq(vk + (iwd + 8 * r));         // lifted inverse input tag sq
      bf_w_hat_sq[r] = bf192_byte_combine_bits_sq(w_dash[(32 * j + 8 * r) / 8]);  // lifted inverse output sq
      bf_v_w_hat_sq[r] = bf192_byte_combine_sq(v_w_dash + (32 * j + 8 * r));      // lifted inverse output tag sq
      #endif

    }
    // Step: 13..17
    for (unsigned int r = 0; r < 4; r++) {

      #if defined(ALLOW_ZERO_SBOX)
      // Psuedoinverse constraints
      // QS multiplication
      const bf192_t mul_tag_1 = bf192_mul(bf_v_k_hat_sq[r], bf_v_w_hat[r]); // mul1.mac0
      const bf192_t mul_tag_2 = bf192_mul(bf_v_k_hat[r], bf_v_w_hat_sq[r]); // mul2.mac0
      // multiply the tags with the values and multiply the results
      const bf192_t mul_val_1 = bf192_mul(                      // mul1.mac1
                                          bf192_add(bf_k_hat_sq[r], bf_v_k_hat_sq[r]),
                                          bf192_add(bf_w_hat[r], bf_v_w_hat[r])
                                          );
      const bf192_t mul_val_2 = bf192_mul(                      // mul2.mac1
                                          bf192_add(bf_k_hat[r], bf_v_k_hat[r]),
                                          bf192_add(bf_w_hat_sq[r], bf_v_w_hat_sq[r])
                                          );

      // the constant term
      zk_hash_192_update(a0_ctx, mul_tag_1);
      // the linear term
      zk_hash_192_update(a1_ctx, bf192_add(
                                            bf192_add(mul_tag_1, bf_k_hat[r]), // mul1.mac0 + x.val  (because here we check x^2*y = x)
                                            bf192_add(mul_val_1, bf_v_k_hat[r]) // mul1.mac1 + x.mac
                                          )
                        );

      // the constant term
      zk_hash_192_update(a0_ctx, mul_tag_2);
      // the linear term
      zk_hash_192_update(a1_ctx, bf192_add(
                                            bf192_add(mul_tag_2, bf_w_hat[r]), // mul2.mac0 + y.val  (because here we check x*y^2 = y)
                                            bf192_add(mul_val_2, bf_v_w_hat[r]) // mul2.mac1 + y.mac
                                          )
                        );
      
      #else
      // Inverse constraint
      const bf192_t tmp_a0 = bf192_mul(bf_v_k_hat[r], bf_v_w_hat[r]);
      zk_hash_192_update(a0_ctx, tmp_a0);
      const bf192_t tmp_a1 = bf192_add(
                                      bf192_add(
                                                bf192_mul(
                                                          bf192_add(
                                                                    bf_k_hat[r], 
                                                                    bf_v_k_hat[r]),
                                                          bf192_add(bf_w_hat[r], 
                                                                    bf_v_w_hat[r])),
                                                bf192_one()),
                                      tmp_a0);
      zk_hash_192_update(a1_ctx, tmp_a1);
      #endif
    }

    iwd = iwd + 192;
  }
}

// Mkey = 1, this is for the verifier
static void aes_key_schedule_constraints_Mkey_1_192(const bf192_t* q, const uint8_t* delta,
                                                    zk_hash_192_ctx* b0_ctx, bf192_t* qk) {
  // Step: 19..20
  aes_key_schedule_forward_192(q, qk);
  bf192_t q_w_dash[FAEST_192F_Ske * 8];
  aes_key_schedule_backward_192(&q[FAEST_192F_LAMBDA], qk, 0, 1, delta, q_w_dash);

  const bf192_t bf_delta      = bf192_load(delta);
  const bf192_t delta_squared = bf192_mul(bf_delta, bf_delta);

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_192F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_192F_Ske / 4; j++) {

    bf192_t bf_q_hat_k[4];
    bf192_t bf_q_hat_w_dash[4];
    #if defined(ALLOW_ZERO_SBOX)
    bf192_t bf_q_hat_k_sq[4];
    bf192_t bf_q_hat_w_dash_sq[4];
    #endif

    for (unsigned int r = 0; r < 4; r++) {

      // Step: 25..26
      bf_q_hat_k[(r + 3) % 4] = bf192_byte_combine(qk + ((iwd + 8 * r)));
      bf_q_hat_w_dash[r]      = bf192_byte_combine(q_w_dash + ((32 * j + 8 * r)));
      #if defined(ALLOW_ZERO_SBOX)
      bf_q_hat_k_sq[(r + 3) % 4] = bf192_byte_combine_sq(qk + ((iwd + 8 * r)));
      bf_q_hat_w_dash_sq[r]      = bf192_byte_combine_sq(q_w_dash + ((32 * j + 8 * r)));
      #endif
    }
    // Step: 38
    for (unsigned int r = 0; r < 4; r++) {

      #if defined(ALLOW_ZERO_SBOX)
      // Quicksilver multiplication
      // multiply the tags
      const bf192_t mul_tag_1 = bf192_mul(bf_q_hat_k_sq[r], bf_q_hat_w_dash[r]);  // mul1.mac0
      const bf192_t mul_tag_2 = bf192_mul(bf_q_hat_k[r], bf_q_hat_w_dash_sq[r]);  // mul2.mac0
      
      // the constnat term
      zk_hash_192_update(b0_ctx, bf192_add(
                                          mul_tag_1,
                                          bf192_mul(delta_squared, bf_q_hat_k[r])
                                          )
      );
      // the constnat term
      zk_hash_192_update(b0_ctx, bf192_add(
                                          mul_tag_2,
                                          bf192_mul(delta_squared, bf_q_hat_w_dash[r])
                                          )
      );
      #else
      // the constant term
      zk_hash_192_update(b0_ctx, bf192_add(
                                          bf192_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]),
                                          delta_squared));
      #endif
    }
    iwd = iwd + 192;
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_forward_192_1(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  bf192_t* bf_y, bf192_t* bf_y_sq)
#else
static void aes_enc_forward_192_1(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  bf192_t* bf_y)
#endif
{
  // called only with Mtag == Mkey == 0

  // Step: 2
  for (unsigned int i = 0; i < 16; i++) {
    // Step: 3,4 (bit spliced)
    // -((1 ^ Mtag) & (1 ^ Mkey)) == 0xFF
    const uint8_t xin = in[i];
    // Step: 5
    bf_y[i] = bf192_add(bf192_byte_combine_bits(xin), bf192_byte_combine_bits(xk[i]));
    #if defined(ALLOW_ZERO_SBOX)
    uint8_t tmp;
    for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
      tmp ^= (  ( (xin >> bit_j) ^ (xk[i] >> bit_j)  ) & 1 ) << bit_j;
    }
    bf_y_sq[i] = bf192_byte_combine_bits_sq(tmp);
    #endif
  }

  const bf192_t bf_two   = bf192_byte_combine_bits(2);
  const bf192_t bf_three = bf192_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_192F_R; j++) {

    for (unsigned int c = 0; c < 4; c++) {

      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_xk_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      uint8_t x_bits[4*8];
      uint8_t x_bits_by_two[8];
      uint8_t output_bits[8];
      #endif


      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          x_bits[r*8 + bit_l] = (x[(ix + 8 * r) / 8] >> bit_l) & 1;
        }
        // mul by two
        x_bits_by_two[0] = x_bits[r*8 + 7];
        x_bits_by_two[1] = x_bits[r*8 + 7] ^ x_bits[r*8 + 0];
        x_bits_by_two[2] = x_bits[r*8 + 1];
        x_bits_by_two[3] = x_bits[r*8 + 7] ^ x_bits[r*8 + 2];
        x_bits_by_two[4] = x_bits[r*8 + 7] ^ x_bits[r*8 + 3];
        x_bits_by_two[5] = x_bits[r*8 + 4];
        x_bits_by_two[6] = x_bits[r*8 + 5];
        x_bits_by_two[7] = x_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_x_hat[r]  = bf192_byte_combine_bits(x[(ix + 8 * r) / 8]);
        bf_xk_hat[r] = bf192_byte_combine_bits(xk[(ik + 8 * r) / 8]);
        #endif
        
      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        uint8_t output_bits = 0;
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits ^= (x_bits_by_two[r*8 + bit_l] ^ 
                            (x_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (x_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (x_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (x_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ xk[(ik + 8 * r) / 8])
                                )
                              )
                            )
                          ) << bit_l;
        }
        bf_y[iy + r] = bf192_byte_combine_bits(output_bits);
        bf_y_sq[iy + r] = bf192_byte_combine_bits_sq(output_bits);
      }
      #else
      // Step : 14
      bf_y[iy + 0] = bf192_add(bf_xk_hat[0], bf192_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf192_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_x_hat[2]);
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_x_hat[3]);

      // Step: 15
      bf_y[iy + 1] = bf192_add(bf_xk_hat[1], bf_x_hat[0]);
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf_x_hat[3]);

      // Step: 16
      bf_y[iy + 2] = bf192_add(bf_xk_hat[2], bf_x_hat[0]);
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf_x_hat[1]);
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf192_add(bf_xk_hat[3], bf192_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_x_hat[1]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_x_hat[2]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_x_hat[3], bf_two));
      #endif

    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_forward_192(const bf192_t* bf_x, const bf192_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf192_t* bf_y, bf192_t* bf_y_sq)
#else
static void aes_enc_forward_192(const bf192_t* bf_x, const bf192_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf192_t* bf_y)
#endif
{
  const bf192_t bf_delta  = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t bf_factor = bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey));

  // Step: 2..4
  for (unsigned int i = 0; i < 16; i++) {
    bf192_t bf_xin[8];
    for (unsigned int j = 0; j < 8; j++) {
      bf_xin[j] = bf192_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
    }
    // Step: 5
    bf_y[i] = bf192_add(bf192_byte_combine(bf_xin), bf192_byte_combine(bf_xk + (8 * i)));
    #if defined(ALLOW_ZERO_SBOX)
    bf192_t tmp[8];
    for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
      tmp[bit_j] = bf_xin[i] ^ bf_xk[i * 8 + bit_j];
    }
    bf_y_sq[i] = bf192_byte_combine_sq(tmp);
    #endif
  }

  const bf192_t bf_two   = bf192_byte_combine_bits(2);
  const bf192_t bf_three = bf192_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_192F_R; j++) {

    for (unsigned int c = 0; c < 4; c++) {

      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_xk_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      bf192_t bf_x_bits[4*8];
      bf192_t bf_x_bits_by_two[8];
      bf192_t output_bits[8];
      #endif

      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          bf_x_bits[r*8 + bit_l] = bf_x[(ix + 8 * r) + bit_l];
        }
        // mul by two
        bf_x_bits_by_two[0] = bf_x_bits[r*8 + 7];
        bf_x_bits_by_two[1] = bf_x_bits[r*8 + 7] ^ bf_x_bits[r*8 + 0];
        bf_x_bits_by_two[2] = bf_x_bits[r*8 + 1];
        bf_x_bits_by_two[3] = bf_x_bits[r*8 + 7] ^ bf_x_bits[r*8 + 2];
        bf_x_bits_by_two[4] = bf_x_bits[r*8 + 7] ^ bf_x_bits[r*8 + 3];
        bf_x_bits_by_two[5] = bf_x_bits[r*8 + 4];
        bf_x_bits_by_two[6] = bf_x_bits[r*8 + 5];
        bf_x_bits_by_two[7] = bf_x_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_x_hat[r]  = bf192_byte_combine(bf_x + (ix + 8 * r));
        bf_xk_hat[r] = bf192_byte_combine(bf_xk + (ik + 8 * r));
        #endif

      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        bf192_t output_bits[8];
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits[bit_l] = (bf_x_bits_by_two[r*8 + bit_l] ^ 
                            (bf_x_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (bf_x_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (bf_x_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (bf_x_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ bf_xk[(ik + 8 * r) / 8])
                                )
                              )
                            )
                          );
        }
        bf_y[iy + r] = bf192_byte_combine(output_bits);
        bf_y_sq[iy + r] = bf192_byte_combine_sq(output_bits);
      }
      #else
      // Step : 14
      bf_y[iy + 0] = bf192_add(bf_xk_hat[0], bf192_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf192_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_x_hat[2]);
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_x_hat[3]);

      // Step: 15
      bf_y[iy + 1] = bf192_add(bf_xk_hat[1], bf_x_hat[0]);
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf_x_hat[3]);

      // Step: 16
      bf_y[iy + 2] = bf192_add(bf_xk_hat[2], bf_x_hat[0]);
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf_x_hat[1]);
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf192_add(bf_xk_hat[3], bf192_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_x_hat[1]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_x_hat[2]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_x_hat[3], bf_two));
      #endif
    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_backward_192_1(const uint8_t* x, const uint8_t* xk, const uint8_t* out, 
                                  bf192_t* y_out, bf192_t* y_out_sq)
#else
static void aes_enc_backward_192_1(const uint8_t* x, const uint8_t* xk, const uint8_t* out, 
                                  bf192_t* y_out)
#endif
{
  // called only with Mtag == Mkey == 0

  uint8_t xtilde;
  // Step:2..4
  for (unsigned int j = 0; j < FAEST_192F_R; j++) {
    for (unsigned int c = 0; c < 4; c++) {
      for (unsigned int r = 0; r < 4; r++) {
        // Step: 5..6
        unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
        if (j < (FAEST_192F_R - 1)) {
          // Step: 7
          xtilde = x[ird / 8];
        } else {
          // Step: 9..11 (bit spliced)
          // -((1 ^ Mtag) & (1 ^ Mkey)) == 0xff
          uint8_t xout = out[(ird - 128 * (FAEST_192F_R - 1)) / 8];
          xtilde       = xout ^ xk[(128 + ird) / 8];
        }

        // Step: 12..17 (bit spliced)
        // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
        uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^ rotr8(xtilde, 2) ^ 0x5;

        // Step: 18
        y_out[16 * j + 4 * c + r] = bf192_byte_combine_bits(ytilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[16 * j + 4 * c + r] = bf192_byte_combine_bits_sq(ytilde);
        #endif
      }
    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_backward_192(const bf192_t* bf_x, const bf192_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf192_t* y_out, bf192_t* y_out_sq)
#else
static void aes_enc_backward_192(const bf192_t* bf_x, const bf192_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf192_t* y_out)
#endif
{
  // Step: 1
  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t factor =
      bf192_mul_bit(bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  // Step: 2..4
  for (unsigned int j = 0; j < FAEST_192F_R; j++) {
    for (unsigned int c = 0; c < 4; c++) {
      for (unsigned int r = 0; r < 4; r++) {
        bf192_t bf_x_tilde[8];
        // Step: 5
        unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
        // Step: 6
        if (j < (FAEST_192F_R - 1)) {
          // Step: 7
          memcpy(bf_x_tilde, bf_x + ird, sizeof(bf_x_tilde));
        } else {
          // Step: 10
          for (unsigned int i = 0; i < 8; ++i) {
            // Step: 11
            bf192_t bf_xout =
                bf192_mul_bit(factor, get_bit(out[(ird - 128 * (FAEST_192F_R - 1)) / 8], i));
            // Step: 12
            bf_x_tilde[i] = bf192_add(bf_xout, bf_xk[128 + ird + i]);
          }
        }
        // Step: 13..17
        bf192_t bf_y_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf192_add(bf192_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                    bf_x_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);

        // Step: 18
        y_out[16 * j + 4 * c + r] = bf192_byte_combine(bf_y_tilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[16 * j + 4 * c + r] = bf192_byte_combine_sq(bf_y_tilde);
        #endif
      }
    }
  }
}

// MKey = 0, this is for the prover
static void aes_enc_constraints_Mkey_0_192(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                           const bf192_t* v, const uint8_t* k, const bf192_t* vk,
                                           zk_hash_192_ctx* a0_ctx, zk_hash_192_ctx* a1_ctx) {
  bf192_t s[FAEST_192F_Senc]; // inverse input
  bf192_t vs[FAEST_192F_Senc]; // inverse input tag
  #if defined(ALLOW_ZERO_SBOX)
  bf192_t s_sq[FAEST_192F_Senc]; // inverse input sq
  bf192_t vs_sq[FAEST_192F_Senc]; // input input tag sq
  aes_enc_forward_192_1(w, k, in, s, s_sq);
  aes_enc_forward_192(v, vk, in, 1, 0, NULL, vs, vs_sq);
  #else
  aes_enc_forward_192_1(w, k, in, s);
  aes_enc_forward_192(v, vk, in, 1, 0, NULL, vs);
  #endif

  bf192_t s_dash[FAEST_192F_Senc]; // inverse output
  bf192_t vs_dash[FAEST_192F_Senc]; // inverse output tag
  #if defined(ALLOW_ZERO_SBOX)
  bf192_t s_dash_sq[FAEST_192F_Senc]; // inverse output sq
  bf192_t vs_dash_sq[FAEST_192F_Senc]; // inverse output tag sq
  aes_enc_backward_192_1(w, k, out, s_dash, s_dash_sq);
  aes_enc_backward_192(v, vk, 1, 0, NULL, out, vs_dash, vs_dash_sq);
  #else
  aes_enc_backward_192_1(w, k, out, s_dash);
  aes_enc_backward_192(v, vk, 1, 0, NULL, out, vs_dash);
  #endif

  for (unsigned int j = 0; j < FAEST_192F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Psuedoinverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf192_t mul_tag_1 = bf192_mul(vs_sq[j], vs_dash[j]);   // mul1.mac0
    const bf192_t mul_tag_2 = bf192_mul(vs[j], vs_dash_sq[j]);   // mul2.mac0
    // add the tags with the values and multiply the results
    const bf192_t mul_val_1 = bf192_mul(                        // mul1.mac1
                        bf192_add(s_sq[j], vs_sq[j]), 
                        bf192_add(s_dash[j], vs_dash[j])
                        );
    const bf192_t mul_val_2 = bf192_mul(                        // mul2.mac1
                        bf192_add(s[j], vs[j]),
                        bf192_add(s_dash_sq[j], vs_dash_sq[j])
                        );

    // the constnat term
    zk_hash_192_update(a0_ctx, mul_tag_1);
    // the linear term
    zk_hash_192_update(a1_ctx, bf192_add(
                                        bf192_add(mul_tag_1, s[j]),   // mul1.mac0 + x.val  (because here we check x^2*y = x)
                                        bf192_add(mul_val_1, vs[j])   // mul1.mac1 + x.mac
                                        )
                      );
    // the constant term
    zk_hash_192_update(a0_ctx, mul_tag_2);
    // the linear term
    zk_hash_192_update(a1_ctx, bf192_add(
                                        bf192_add(mul_tag_2, s_dash[j]),  // mul2.mac0 + y.val  (because here we check x*y^2 = y)
                                        bf192_add(mul_val_2, vs_dash[j])  // mul2.mac1 + y.mac
                                        )
                      );    

    #else

    // Inverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf192_t mul_tag = bf192_mul(vs[j], vs_dash[j]);   // mac0
    // add the tags with the values and multiply the results
    const bf192_t mul_val = bf192_mul(bf192_add(s[j], vs[j]), bf192_add(s_dash[j], vs_dash[j]));    // mac1

    // this is the constant term
    zk_hash_192_update(a0_ctx, mul_tag);    
    // this is the linear term
    zk_hash_192_update(a1_ctx, bf192_add(
                                        bf192_add(mul_tag, bf192_one()),    // mac0 + 1     (because here we check x*y = 1)
                                        mul_val                             // (mac0 + 1) + mac1
                                        )
                      );
    #endif
    // // instead of storing in A0, A1, hash it
    // const bf192_t tmp = bf192_mul(vs[j], vs_dash[j]);
    // zk_hash_192_update(a0_ctx, tmp);
    // zk_hash_192_update(a1_ctx, bf192_add(bf192_add(bf192_mul(bf192_add(s[j], vs[j]),
    //                                                          bf192_add(s_dash[j], vs_dash[j])),
    //                                                tmp),
    //                                      bf192_one()));
  }
}

// Mkey = 1, this is for the verifier
static void aes_enc_constraints_Mkey_1_192(const uint8_t* in, const uint8_t* out, const bf192_t* q,
                                           const bf192_t* qk, const uint8_t* delta,
                                           zk_hash_192_ctx* b0_ctx) {
  // Step: 11..12
  bf192_t qs[FAEST_192F_Senc];
  bf192_t qs_dash[FAEST_192F_Senc];
  #if defined(ALLOW_ZERO_SBOX)
  bf192_t qs_sq[FAEST_192F_Senc];
  bf192_t qs_dash_sq[FAEST_192F_Senc];
  aes_enc_forward_192(q, qk, in, 0, 1, delta, qs, qs_sq);
  aes_enc_backward_192(q, qk, 0, 1, delta, out, qs_dash, qs_dash_sq);
  #else
  aes_enc_forward_192(q, qk, in, 0, 1, delta, qs);
  aes_enc_backward_192(q, qk, 0, 1, delta, out, qs_dash);
  #endif
  

  // Step: 13..14
  bf192_t delta_squared = bf192_mul(bf192_load(delta), bf192_load(delta));
  
  for (unsigned int j = 0; j < FAEST_192F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Quicksilver multiplication
    // multiply the tags
    const bf192_t mul_tag_1 = bf192_mul(qs_sq[j], qs_dash[j]);   // mul1.mac0
    const bf192_t mul_tag_2 = bf192_mul(qs[j], qs_dash_sq[j]);   // mul2.mac0

    // the constnat term
    zk_hash_192_update(b0_ctx, bf192_add(
                                        mul_tag_1,
                                        bf192_mul(delta_squared, qs[j])
                                        )
    );
    // the constnat term
    zk_hash_192_update(b0_ctx, bf192_add(
                                        mul_tag_2,
                                        bf192_mul(delta_squared, qs_dash[j])
                                        )
    );
    #else
    // the constnat term
    zk_hash_192_update(b0_ctx, bf192_add(
                                        bf192_mul(qs[j], qs_dash[j]),
                                        delta_squared));
    #endif
    // // instead of storing it, hash it
    // zk_hash_192_update(b0_ctx, bf192_add(bf192_mul(qs[j], qs_dash[j]), delta_squared));
  }
}

/* static void aes_prove_192(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* owf_in,
                          const uint8_t* owf_out, const uint8_t* chall_2, uint8_t* a0_tilde,
                          uint8_t* a12_tilde, const faest_paramset_t* params) {
  // Step: 1..2
  bf192_t* bf_v = column_to_row_major_and_shrink_V_192(V, FAEST_192F_ELL);

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7 + 18
  uint8_t* k  = malloc((FAEST_192F_R + 1) * 128 / 8);
  bf192_t* vk = faest_aligned_alloc(BF192_ALIGN, sizeof(bf192_t) * ((FAEST_192F_R + 1) * 128));
  zk_hash_192_ctx a0_ctx;
  zk_hash_192_ctx a1_ctx;

  zk_hash_192_init(&a0_ctx, chall_2);
  zk_hash_192_init(&a1_ctx, chall_2);
  aes_key_schedule_constraints_Mkey_0_192(w, bf_v, &a0_ctx, &a1_ctx, k, vk, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  aes_enc_constraints_Mkey_0_192(owf_in, owf_out, w + FAEST_192F_Lke / 8, bf_v + FAEST_192F_Lke, k, vk,
                                 &a0_ctx, &a1_ctx);
  // Step: 12-15
  aes_enc_constraints_Mkey_0_192(owf_in + 16, owf_out + 16, w + (FAEST_192F_Lke + FAEST_192F_Lenc) / 8,
                                 bf_v + FAEST_192F_Lke + FAEST_192F_Lenc, k, vk, &a0_ctx, &a1_ctx);
  faest_aligned_free(vk);
  free(k);

  // Step: 16..18
  zk_hash_192_finalize(a_tilde, &a1_ctx, bf192_load(u + FAEST_192F_ELL / 8));
  zk_hash_192_finalize(b_tilde, &a0_ctx, bf192_sum_poly(bf_v + FAEST_192F_ELL));
  faest_aligned_free(bf_v);
}

static uint8_t* aes_verify_192(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                               const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                               const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.tau0;
  const unsigned int k0  = params->faest_param.k;
  const unsigned int t1  = params->faest_param.tau1;
  const unsigned int k1  = (t0 != 0) ? k0 - 1 : k0;

  // Step: 1
  const uint8_t* delta = chall_3;
  // Step: 2,3
  // do nothing

  // Step: 4..10
  for (unsigned int i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t decoded_challenge[MAX_DEPTH];
    ChalDec(chall_3, i, k0, t0, k1, t1, decoded_challenge);
    for (unsigned int j = 0; j < depth; j++, ++col) {
      if (decoded_challenge[j] == 1) {
        xor_u8_array(d, Q[col], Q[col], (FAEST_192F_ELL + 7) / 8);
      }
    }
  }

  // Step: 11..12
  bf192_t* bf_q = column_to_row_major_and_shrink_V_192(Q, FAEST_192F_ELL);

  // Step: 13 + 21
  bf192_t* qk = faest_aligned_alloc(BF192_ALIGN, sizeof(bf192_t) * ((FAEST_192F_R + 1) * 128));
  // instead of storing B_0 in an array, we process the values with zk_hash_128
  zk_hash_192_ctx b0_ctx;
  zk_hash_192_init(&b0_ctx, chall_2);
  aes_key_schedule_constraints_Mkey_1_192(bf_q, delta, &b0_ctx, qk);

  // Step: 14
  aes_enc_constraints_Mkey_1_192(in, out, bf_q + FAEST_192F_Lke, qk, delta, &b0_ctx);

  // Step: 18
  aes_enc_constraints_Mkey_1_192(in + 16, out + 16, bf_q + FAEST_192F_Lke + FAEST_192F_Lenc, qk,
                                 delta, &b0_ctx);
  faest_aligned_free(qk);

  // Step: 20+21
  uint8_t* q_tilde = malloc(FAEST_192F_LAMBDA / 8);
  zk_hash_192_finalize(q_tilde, &b0_ctx, bf192_sum_poly(bf_q + FAEST_192F_ELL));
  faest_aligned_free(bf_q);

  bf192_t bf_qtilde = bf192_load(q_tilde);
  bf192_store(q_tilde, bf192_add(bf_qtilde, bf192_mul(bf192_load(a_tilde), bf192_load(delta))));

  return q_tilde;
}

 */


// ###########################################################################################################################################
// ##################################           LAMBDA = 256            ######################################################################
// ###########################################################################################################################################

static bf256_t* column_to_row_major_and_shrink_V_256(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf256_t* new_v = faest_aligned_alloc(BF256_ALIGN, (ell + FAEST_256F_LAMBDA) * sizeof(bf256_t));
  for (unsigned int row = 0; row != ell + FAEST_256F_LAMBDA; ++row) {
    uint8_t new_row[BF256_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_256F_LAMBDA; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf256_load(new_row);
  }

  return new_v;
}

static void aes_key_schedule_forward_256(const bf256_t* v, bf256_t* bf_out) {
  // Step: 1 sanity check (skipped)

  memcpy(bf_out, v, sizeof(bf256_t) * FAEST_256F_LAMBDA);

  // Step: 4
  unsigned int i_wd = FAEST_256F_LAMBDA;
  // Step: 5..10
  for (unsigned int j = FAEST_256F_Nwd; j < 4 * (FAEST_256F_R + 1); j++) {
    if ((j % FAEST_256F_Nwd) == 0 || (FAEST_256F_Nwd > 6 && (j % FAEST_256F_Nwd) == 4)) {
      memcpy(bf_out + j * 32, v + i_wd, sizeof(bf256_t) * 32);
      i_wd += 32;
    } else {
      for (unsigned int i = 0; i < 32; i++) {
        bf_out[(32 * j) + i] =
            bf256_add(bf_out[32 * (j - FAEST_256F_Nwd) + i], bf_out[32 * (j - 1) + i]);
      }
    }
  }
}

static void aes_key_schedule_backward_256(const bf256_t* v, const bf256_t* Vk, uint8_t Mtag,
                                          uint8_t Mkey, const uint8_t* delta, bf256_t* bf_out) {
  // Step: 1
  assert(!((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)));

  unsigned int iwd   = 0;
  unsigned int c     = 0;
  bool rmvRcon       = true;
  unsigned int ircon = 0;

  const bf256_t bf_delta      = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t bf_minus_mkey = bf256_from_bit(1 ^ Mkey);
  const uint8_t minus_mtag    = 1 ^ Mtag;
  bf256_t bf_mkey_times_delta = bf256_mul_bit(bf_delta, Mkey);
  bf_mkey_times_delta         = bf256_add(bf_mkey_times_delta, bf_minus_mkey);

  for (unsigned int j = 0; j < FAEST_256F_Ske; j++) {
    // Step 7
    bf256_t bf_x_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      bf_x_tilde[i] = bf256_add(v[8 * j + i], Vk[iwd + 8 * c + i]);
    }

    if (Mtag == 0 && rmvRcon == true && c == 0) {
      // Step 9
      uint8_t r = Rcon[ircon];
      ircon     = ircon + 1;

      bf256_t bf_r[8];
      for (unsigned int i = 0; i < 8; i++) {
        // Step 12
        bf_r[i] = bf256_mul_bit(bf_mkey_times_delta, get_bit(r, i));
        // Step 13
        bf_x_tilde[i] = bf256_add(bf_x_tilde[i], bf_r[i]);
      }
    }

    for (unsigned int i = 0; i < 8; ++i) {
      bf_out[i + 8 * j] = bf256_add(bf256_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                    bf_x_tilde[(i + 2) % 8]);
    }
    bf_out[0 + 8 * j] =
        bf256_add(bf_out[0 + 8 * j], bf256_mul_bit(bf_mkey_times_delta, minus_mtag));
    bf_out[2 + 8 * j] =
        bf256_add(bf_out[2 + 8 * j], bf256_mul_bit(bf_mkey_times_delta, minus_mtag));
    c = c + 1;

    if (c == 4) {
      c = 0;
      iwd += 128;
      rmvRcon = !rmvRcon;
    }
  }
}

// Mkey = 0, this is for the prover
static void aes_key_schedule_constraints_Mkey_0_256(const uint8_t* w, const bf256_t* v,
                                                    zk_hash_256_ctx* a0_ctx,
                                                    zk_hash_256_ctx* a1_ctx, uint8_t* k,
                                                    bf256_t* vk, const faest_paramset_t* params) {
  // for scan-build
  assert(FAEST_256F_Ske == params->faest_param.Ske);

  // Step: 2
  aes_key_schedule_forward_1(w, k, params);

  // Step: 3
  aes_key_schedule_forward_256(v, vk);

  // Step: 4
  uint8_t w_dash[FAEST_256F_Ske];
  aes_key_schedule_backward_1(w + FAEST_256F_LAMBDA / 8, k, w_dash, params);

  // Step: 5
  bf256_t v_w_dash[FAEST_256F_Ske * 8];
  aes_key_schedule_backward_256(v + FAEST_256F_LAMBDA, vk, 1, 0, NULL, v_w_dash);

  // Step: 6..8
  bool rotate_word = true;
  unsigned int iwd = 32 * (FAEST_256F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_256F_Ske / 4; j++) {
    bf256_t bf_k_hat[4];
    bf256_t bf_v_k_hat[4];
    bf256_t bf_w_hat[4];
    bf256_t bf_v_w_hat[4];
    #if defined(ALLOW_ZERO_SBOX)
    bf256_t bf_k_hat_sq[4];
    bf256_t bf_v_k_hat_sq[4];
    bf256_t bf_w_hat_sq[4];
    bf256_t bf_v_w_hat_sq[4];
    #endif

    for (unsigned int r = 0; r < 4; r++) {
      // Step: 10..16
      if (rotate_word) {
        bf_k_hat[(r + 3) % 4]   = bf256_byte_combine_bits(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat[(r + 3) % 4] = bf256_byte_combine(vk + (iwd + 8 * r));
        bf_w_hat[r]        = bf256_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
        bf_v_w_hat[r]      = bf256_byte_combine(v_w_dash + (32 * j + 8 * r));
      } else {
        bf_k_hat[r]        = bf256_byte_combine_bits(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat[r]      = bf256_byte_combine(vk + (iwd + 8 * r));
        bf_w_hat[r]   = bf256_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
        bf_v_w_hat[r] = bf256_byte_combine(v_w_dash + (32 * j + 8 * r));
      }

      // squaring the bits and macs and lifting
      #if defined(ALLOW_ZERO_SBOX)
      if (rotate_word) {
        bf_k_hat_sq[(r + 3) % 4]   = bf256_byte_combine_bits_sq(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat_sq[(r + 3) % 4] = bf256_byte_combine_sq(vk + (iwd + 8 * r));
        bf_w_hat_sq[r]        = bf256_byte_combine_bits_sq(w_dash[(32 * j + 8 * r) / 8]);
        bf_v_w_hat_sq[r]      = bf256_byte_combine_sq(v_w_dash + (32 * j + 8 * r));
      } else {
        bf_k_hat_sq[r]        = bf256_byte_combine_bits_sq(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat_sq[r]      = bf256_byte_combine_sq(vk + (iwd + 8 * r));
        bf_w_hat_sq[r]   = bf256_byte_combine_bits_sq(w_dash[(32 * j + 8 * r) / 8]);
        bf_v_w_hat_sq[r] = bf256_byte_combine_sq(v_w_dash + (32 * j + 8 * r));
      }
      #endif

    }
    // Step: 13..17
    for (unsigned int r = 0; r < 4; r++) {

      #if defined(ALLOW_ZERO_SBOX)
      // Psuedoinverse constraints
      // QS multiplication
      const bf256_t mul_tag_1 = bf256_mul(bf_v_k_hat_sq[r], bf_v_w_hat[r]); // mul1.mac0
      const bf256_t mul_tag_2 = bf256_mul(bf_v_k_hat[r], bf_v_w_hat_sq[r]); // mul2.mac0
      // multiply the tags with the values and multiply the results
      const bf256_t mul_val_1 = bf256_mul(                      // mul1.mac1
                                          bf256_add(bf_k_hat_sq[r], bf_v_k_hat_sq[r]),
                                          bf256_add(bf_w_hat[r], bf_v_w_hat[r])
                                          );
      const bf256_t mul_val_2 = bf256_mul(                      // mul2.mac1
                                          bf256_add(bf_k_hat[r], bf_v_k_hat[r]),
                                          bf256_add(bf_w_hat_sq[r], bf_v_w_hat_sq[r])
                                          );

      // the constant term
      zk_hash_256_update(a0_ctx, mul_tag_1);
      // the linear term
      zk_hash_256_update(a1_ctx, bf256_add(
                                            bf256_add(mul_tag_1, bf_k_hat[r]), // mul1.mac0 + x.val  (because here we check x^2*y = x)
                                            bf256_add(mul_val_1, bf_v_k_hat[r]) // mul1.mac1 + x.mac
                                          )
                        );

      // the constant term
      zk_hash_256_update(a0_ctx, mul_tag_2);
      // the linear term
      zk_hash_256_update(a1_ctx, bf256_add(
                                            bf256_add(mul_tag_2, bf_w_hat[r]), // mul2.mac0 + y.val  (because here we check x*y^2 = y)
                                            bf256_add(mul_val_2, bf_v_w_hat[r]) // mul2.mac1 + y.mac
                                          )
                        );
      #else
      // Inverse constraint
      const bf256_t tmp_a0 = bf256_mul(bf_v_k_hat[r], bf_v_w_hat[r]);
      zk_hash_256_update(a0_ctx, tmp_a0);
      const bf256_t tmp_a1 = bf256_add(
                                      bf256_add(
                                                bf256_mul(
                                                          bf256_add(
                                                                    bf_k_hat[r], 
                                                                    bf_v_k_hat[r]),
                                                          bf256_add(bf_w_hat[r], 
                                                                    bf_v_w_hat[r])),
                                                bf256_one()),
                                      tmp_a0);
      zk_hash_256_update(a1_ctx, tmp_a1);
      #endif

      // const bf256_t tmp = bf256_mul(bf_v_k_hat[r], bf_v_w_hat[r]);
      // zk_hash_256_update(a0_ctx, tmp);
      // zk_hash_256_update(
      //     a1_ctx, bf256_add(bf256_add(bf256_mul(bf256_add(bf_k_hat[r], bf_v_k_hat[r]),
      //                                           bf256_add(bf_w_hat[r], bf_v_w_hat[r])),
      //                                 bf256_one()),
      //                       tmp));


    }
    iwd         = iwd + 128;
    rotate_word = !rotate_word;
  }
}

// Mkey = 1, this is for the verifier
static void aes_key_schedule_constraints_Mkey_1_256(const bf256_t* q, const uint8_t* delta,
                                                    zk_hash_256_ctx* b0_ctx, bf256_t* qk) {
  // Step: 19..20
  aes_key_schedule_forward_256(q, qk);
  bf256_t q_w_dash[FAEST_256F_Ske * 8];
  aes_key_schedule_backward_256(&q[FAEST_256F_LAMBDA], qk, 0, 1, delta, q_w_dash);

  const bf256_t bf_delta         = bf256_load(delta);
  const bf256_t delta_squared = bf256_mul(bf_delta, bf_delta);
  bool rotate_word               = true;

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_256F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_256F_Ske / 4; j++) {

    bf256_t bf_q_hat_k[4];
    bf256_t bf_q_hat_w_dash[4];
    #if defined(ALLOW_ZERO_SBOX)
    bf256_t bf_q_hat_k_sq[4];
    bf256_t bf_q_hat_w_dash_sq[4];
    #endif

    for (unsigned int r = 0; r < 4; r++) {

      // Step: 25..26
      if (rotate_word) {
        bf_q_hat_k[(r + 3) % 4] = bf256_byte_combine(qk + ((iwd + 8 * r)));
        bf_q_hat_w_dash[r]      = bf256_byte_combine(q_w_dash + ((32 * j + 8 * r)));
      } else {
        bf_q_hat_k[r]      = bf256_byte_combine(qk + ((iwd + 8 * r)));
        bf_q_hat_w_dash[r] = bf256_byte_combine(q_w_dash + ((32 * j + 8 * r)));
      }

      #if defined(ALLOW_ZERO_SBOX)
      if (rotate_word) {
        bf_q_hat_k_sq[(r + 3) % 4] = bf256_byte_combine_sq(qk + ((iwd + 8 * r)));
        bf_q_hat_w_dash_sq[r]      = bf256_byte_combine_sq(q_w_dash + ((32 * j + 8 * r)));
      } else {
        bf_q_hat_k_sq[r]      = bf256_byte_combine_sq(qk + ((iwd + 8 * r)));
        bf_q_hat_w_dash_sq[r] = bf256_byte_combine_sq(q_w_dash + ((32 * j + 8 * r)));
      }
      #endif
    }
    // Step: 38
    for (unsigned int r = 0; r < 4; r++) {

      #if defined(ALLOW_ZERO_SBOX)
      // Quicksilver multiplication
      // multiply the tags
      const bf256_t mul_tag_1 = bf256_mul(bf_q_hat_k_sq[r], bf_q_hat_w_dash[r]);  // mul1.mac0
      const bf256_t mul_tag_2 = bf256_mul(bf_q_hat_k[r], bf_q_hat_w_dash_sq[r]);  // mul2.mac0
      
      // the constnat term
      zk_hash_256_update(b0_ctx, bf256_add(
                                          mul_tag_1,
                                          bf256_mul(delta_squared, bf_q_hat_k[r])
                                          )
      );
      // the constnat term
      zk_hash_256_update(b0_ctx, bf256_add(
                                          mul_tag_2,
                                          bf256_mul(delta_squared, bf_q_hat_w_dash[r])
                                          )
      );
      #else
      // the constant term
      zk_hash_256_update(b0_ctx, bf256_add(
                                          bf256_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]),
                                          delta_squared));
      #endif
      // bf256_t bf_tmp = bf256_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      // zk_hash_256_update(b0_ctx, bf256_add(bf_tmp, bf_delta_squared));
    }
    iwd         = iwd + 128;
    rotate_word = !rotate_word;
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_forward_256_1(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  bf256_t* bf_y, bf256_t* bf_y_sq) 
#else
static void aes_enc_forward_256_1(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  bf256_t* bf_y)
#endif                       
{
  // called only with Mtag == Mkey == 0

  // Step: 2
  for (unsigned int i = 0; i < 16; i++) {
    // Step: 3,4 (bit spliced)
    const uint8_t xin = in[i];
    // Step: 5
    bf_y[i] = bf256_add(bf256_byte_combine_bits(xin), bf256_byte_combine_bits(xk[i]));
    #if defined(ALLOW_ZERO_SBOX)
    uint8_t tmp;
    for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
      tmp ^= (  ( (xin >> bit_j) ^ (xk[i] >> bit_j)  ) & 1 ) << bit_j;
    }
    bf_y_sq[i] = bf256_byte_combine_bits_sq(tmp);
    #endif
  }

  const bf256_t bf_two   = bf256_byte_combine_bits(2);
  const bf256_t bf_three = bf256_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_256F_R; j++) {

    for (unsigned int c = 0; c < 4; c++) {

      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_xk_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      uint8_t x_bits[4*8];
      uint8_t x_bits_by_two[8];
      uint8_t output_bits[8];
      #endif


      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          x_bits[r*8 + bit_l] = (x[(ix + 8 * r) / 8] >> bit_l) & 1;
        }
        // mul by two
        x_bits_by_two[0] = x_bits[r*8 + 7];
        x_bits_by_two[1] = x_bits[r*8 + 7] ^ x_bits[r*8 + 0];
        x_bits_by_two[2] = x_bits[r*8 + 1];
        x_bits_by_two[3] = x_bits[r*8 + 7] ^ x_bits[r*8 + 2];
        x_bits_by_two[4] = x_bits[r*8 + 7] ^ x_bits[r*8 + 3];
        x_bits_by_two[5] = x_bits[r*8 + 4];
        x_bits_by_two[6] = x_bits[r*8 + 5];
        x_bits_by_two[7] = x_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_x_hat[r]  = bf256_byte_combine_bits(x[(ix + 8 * r) / 8]);
        bf_xk_hat[r] = bf256_byte_combine_bits(xk[(ik + 8 * r) / 8]);
        #endif

      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        uint8_t output_bits = 0;
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits ^= (x_bits_by_two[r*8 + bit_l] ^ 
                            (x_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (x_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (x_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (x_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ xk[(ik + 8 * r) / 8])
                                )
                              )
                            )
                          ) << bit_l;
        }
        bf_y[iy + r] = bf256_byte_combine_bits(output_bits);
        bf_y_sq[iy + r] = bf256_byte_combine_bits_sq(output_bits);
      }
      #else
      // Step : 14
      bf_y[iy + 0] = bf256_add(bf_xk_hat[0], bf256_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf256_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_x_hat[2]);
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_x_hat[3]);

      // Step: 15
      bf_y[iy + 1] = bf256_add(bf_xk_hat[1], bf_x_hat[0]);
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf_x_hat[3]);

      // Step: 16
      bf_y[iy + 2] = bf256_add(bf_xk_hat[2], bf_x_hat[0]);
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf_x_hat[1]);
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf256_add(bf_xk_hat[3], bf256_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_x_hat[1]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_x_hat[2]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_x_hat[3], bf_two));
      #endif

    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_forward_256(const bf256_t* bf_x, const bf256_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf256_t* bf_y, bf256_t* bf_y_sq)
#else
static void aes_enc_forward_256(const bf256_t* bf_x, const bf256_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf256_t* bf_y)
#endif

{
  const bf256_t bf_delta  = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t bf_factor = bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey));

  // Step: 2..4
  for (unsigned int i = 0; i < 16; i++) {
    bf256_t bf_xin[8];
    for (unsigned int j = 0; j < 8; j++) {
      bf_xin[j] = bf256_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
    }
    // Step: 5
    bf_y[i] = bf256_add(bf256_byte_combine(bf_xin), bf256_byte_combine(bf_xk + (8 * i)));
    #if defined(ALLOW_ZERO_SBOX)
    bf256_t tmp[8];
    for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
      tmp[bit_j] = bf_xin[i] ^ bf_xk[i * 8 + bit_j];
    }
    bf_y_sq[i] = bf256_byte_combine_sq(tmp);
    #endif
  }

  const bf256_t bf_two   = bf256_byte_combine_bits(2);
  const bf256_t bf_three = bf256_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_256S_R; j++) {

    for (unsigned int c = 0; c < 4; c++) {

      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_xk_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      bf256_t bf_x_bits[4*8];
      bf256_t bf_x_bits_by_two[8];
      bf256_t output_bits[8];
      #endif

      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          bf_x_bits[r*8 + bit_l] = bf_x[(ix + 8 * r) + bit_l];
        }
        // mul by two
        bf_x_bits_by_two[0] = bf_x_bits[r*8 + 7];
        bf_x_bits_by_two[1] = bf_x_bits[r*8 + 7] ^ bf_x_bits[r*8 + 0];
        bf_x_bits_by_two[2] = bf_x_bits[r*8 + 1];
        bf_x_bits_by_two[3] = bf_x_bits[r*8 + 7] ^ bf_x_bits[r*8 + 2];
        bf_x_bits_by_two[4] = bf_x_bits[r*8 + 7] ^ bf_x_bits[r*8 + 3];
        bf_x_bits_by_two[5] = bf_x_bits[r*8 + 4];
        bf_x_bits_by_two[6] = bf_x_bits[r*8 + 5];
        bf_x_bits_by_two[7] = bf_x_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_x_hat[r]  = bf256_byte_combine(bf_x + (ix + 8 * r));
        bf_xk_hat[r] = bf256_byte_combine(bf_xk + (ik + 8 * r));
        #endif

      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        bf256_t output_bits[8];
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits[bit_l] = (bf_x_bits_by_two[r*8 + bit_l] ^ 
                            (bf_x_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (bf_x_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (bf_x_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (bf_x_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ bf_xk[(ik + 8 * r) / 8])
                                )
                              )
                            )
                          );
        }
        bf_y[iy + r] = bf256_byte_combine(output_bits);
        bf_y_sq[iy + r] = bf256_byte_combine_sq(output_bits);
      }
      #else
      // Step : 14
      bf_y[iy + 0] = bf256_add(bf_xk_hat[0], bf256_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf256_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_x_hat[2]);
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_x_hat[3]);

      // Step: 15
      bf_y[iy + 1] = bf256_add(bf_xk_hat[1], bf_x_hat[0]);
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf_x_hat[3]);

      // Step: 16
      bf_y[iy + 2] = bf256_add(bf_xk_hat[2], bf_x_hat[0]);
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf_x_hat[1]);
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf256_add(bf_xk_hat[3], bf256_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_x_hat[1]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_x_hat[2]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_x_hat[3], bf_two));
      #endif
      
    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_backward_256_1(const uint8_t* x, const uint8_t* xk, const uint8_t* out, 
                                  bf256_t* y_out, bf256_t* y_out_sq) 
#else
static void aes_enc_backward_256_1(const uint8_t* x, const uint8_t* xk, const uint8_t* out, 
                                  bf256_t* y_out) 
#endif                                  
{
  uint8_t xtilde;
  // Step:2..4
  for (unsigned int j = 0; j < FAEST_256F_R; j++) {
    for (unsigned int c = 0; c < 4; c++) {
      for (unsigned int r = 0; r < 4; r++) {
        // Step: 5..6
        unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
        if (j < (FAEST_256F_R - 1)) {
          // Step: 7
          xtilde = x[ird / 8];
        } else {
          // Step: 9..11 (bit spliced)
          // -((1 ^ Mtag) & (1 ^ Mkey)) == 0xff
          uint8_t xout = out[(ird - 128 * (FAEST_256F_R - 1)) / 8];
          xtilde       = xout ^ xk[(128 + ird) / 8];
        }

        // Step: 12..17 (bit spliced)
        // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
        uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^ rotr8(xtilde, 2) ^ 0x5;

        // Step: 18
        y_out[16 * j + 4 * c + r] = bf256_byte_combine_bits(ytilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[16 * j + 4 * c + r] = bf256_byte_combine_bits_sq(ytilde);
        #endif
      }
    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void aes_enc_backward_256(const bf256_t* bf_x, const bf256_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf256_t* y_out, bf256_t* y_out_sq) {
#else
static void aes_enc_backward_256(const bf256_t* bf_x, const bf256_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf256_t* y_out)
#endif
  // Step: 1
  const bf256_t bf_delta = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t factor =
      bf256_mul_bit(bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  // Step: 2..4
  for (unsigned int j = 0; j < FAEST_256F_R; j++) {
    for (unsigned int c = 0; c < 4; c++) {
      for (unsigned int r = 0; r < 4; r++) {
        bf256_t bf_x_tilde[8];
        // Step: 5
        unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
        // Step: 6
        if (j < (FAEST_256F_R - 1)) {
          // Step: 7
          memcpy(bf_x_tilde, bf_x + ird, sizeof(bf_x_tilde));
        } else {
          // Step: 10
          for (unsigned int i = 0; i < 8; ++i) {
            // Step: 11
            bf256_t bf_xout =
                bf256_mul_bit(factor, get_bit(out[(ird - 128 * (FAEST_256F_R - 1)) / 8], i));
            // Step: 12
            bf_x_tilde[i] = bf256_add(bf_xout, bf_xk[128 + ird + i]);
          }
        }
        // Step: 13..17
        bf256_t bf_y_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf256_add(bf256_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                    bf_x_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);

        // Step: 18
        y_out[16 * j + 4 * c + r] = bf256_byte_combine(bf_y_tilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[16 * j + 4 * c + r] = bf256_byte_combine_sq(bf_y_tilde);
        #endif
      }
    }
  }
}

// MKey = 0, this is for the prover
static void aes_enc_constraints_Mkey_0_256(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                           const bf256_t* v, const uint8_t* k, const bf256_t* vk,
                                           zk_hash_256_ctx* a0_ctx, zk_hash_256_ctx* a1_ctx) {
  bf256_t s[FAEST_256F_Senc];
  bf256_t vs[FAEST_256F_Senc];
  #if defined(ALLOW_ZERO_SBOX)
  bf256_t s_sq[FAEST_256F_Senc]; // inverse input sq
  bf256_t vs_sq[FAEST_256F_Senc]; // input input tag sq
  aes_enc_forward_256_1(w, k, in, s, s_sq);
  aes_enc_forward_256(v, vk, in, 1, 0, NULL, vs, vs_sq);
  #else
  aes_enc_forward_256_1(w, k, in, s);
  aes_enc_forward_256(v, vk, in, 1, 0, NULL, vs);
  #endif

  bf256_t s_dash[FAEST_256F_Senc]; // inverse output
  bf256_t vs_dash[FAEST_256F_Senc]; // inverse output tag
  #if defined(ALLOW_ZERO_SBOX)
  bf256_t s_dash_sq[FAEST_256F_Senc]; // inverse output sq
  bf256_t vs_dash_sq[FAEST_256F_Senc]; // inverse output tag sq
  aes_enc_backward_256_1(w, k, out, s_dash, s_dash_sq);
  aes_enc_backward_256(v, vk, 1, 0, NULL, out, vs_dash, vs_dash_sq);
  #else
  aes_enc_backward_256_1(w, k, out, s_dash);
  aes_enc_backward_256(v, vk, 1, 0, NULL, out, vs_dash);
  #endif

  for (unsigned int j = 0; j < FAEST_256F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Psuedoinverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf256_t mul_tag_1 = bf256_mul(vs_sq[j], vs_dash[j]);   // mul1.mac0
    const bf256_t mul_tag_2 = bf256_mul(vs[j], vs_dash_sq[j]);   // mul2.mac0
    // add the tags with the values and multiply the results
    const bf256_t mul_val_1 = bf256_mul(                        // mul1.mac1
                        bf256_add(s_sq[j], vs_sq[j]), 
                        bf256_add(s_dash[j], vs_dash[j])
                        );
    const bf256_t mul_val_2 = bf256_mul(                        // mul2.mac1
                        bf256_add(s[j], vs[j]),
                        bf256_add(s_dash_sq[j], vs_dash_sq[j])
                        );

    // the constnat term
    zk_hash_256_update(a0_ctx, mul_tag_1);
    // the linear term
    zk_hash_256_update(a1_ctx, bf256_add(
                                        bf256_add(mul_tag_1, s[j]),   // mul1.mac0 + x.val  (because here we check x^2*y = x)
                                        bf256_add(mul_val_1, vs[j])   // mul1.mac1 + x.mac
                                        )
                      );
    // the constant term
    zk_hash_256_update(a0_ctx, mul_tag_2);
    // the linear term
    zk_hash_256_update(a1_ctx, bf256_add(
                                        bf256_add(mul_tag_2, s_dash[j]),  // mul2.mac0 + y.val  (because here we check x*y^2 = y)
                                        bf256_add(mul_val_2, vs_dash[j])  // mul2.mac1 + y.mac
                                        )
                      );    

    #else

    // Inverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf256_t mul_tag = bf256_mul(vs[j], vs_dash[j]);   // mac0
    // add the tags with the values and multiply the results
    const bf256_t mul_val = bf256_mul(bf256_add(s[j], vs[j]), bf256_add(s_dash[j], vs_dash[j]));    // mac1

    // this is the constant term
    zk_hash_256_update(a0_ctx, mul_tag);    
    // this is the linear term
    zk_hash_256_update(a1_ctx, bf256_add(
                                        bf256_add(mul_tag, bf256_one()),    // mac0 + 1     (because here we check x*y = 1)
                                        mul_val                             // (mac0 + 1) + mac1
                                        )
                      );
    #endif
    // const bf256_t tmp = bf256_mul(vs[j], vs_dash[j]);
    // zk_hash_256_update(a0_ctx, tmp);
    // zk_hash_256_update(a1_ctx, bf256_add(bf256_add(bf256_mul(bf256_add(s[j], vs[j]),
    //                                                          bf256_add(s_dash[j], vs_dash[j])),
    //                                                tmp),
    //                                      bf256_one()));
  }
}

// Mkey = 1, this is for the verifier
static void aes_enc_constraints_Mkey_1_256(const uint8_t* in, const uint8_t* out, const bf256_t* q,
                                           const bf256_t* qk, const uint8_t* delta,
                                           zk_hash_256_ctx* b0_ctx) {
  // Step: 11..12
  bf256_t qs[FAEST_256F_Senc];
  bf256_t qs_dash[FAEST_256F_Senc];
  #if defined(ALLOW_ZERO_SBOX)
  bf256_t qs_sq[FAEST_256F_Senc];
  bf256_t qs_dash_sq[FAEST_256F_Senc];
  aes_enc_forward_256(q, qk, in, 0, 1, delta, qs, qs_sq);
  aes_enc_backward_256(q, qk, 0, 1, delta, out, qs_dash, qs_dash_sq);
  #else
  aes_enc_forward_256(q, qk, in, 0, 1, delta, qs);
  aes_enc_backward_256(q, qk, 0, 1, delta, out, qs_dash);
  #endif

  // Step: 13..14
  bf256_t delta_squared = bf256_mul(bf256_load(delta), bf256_load(delta));
  
  for (unsigned int j = 0; j < FAEST_256F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Quicksilver multiplication
    // multiply the tags
    const bf256_t mul_tag_1 = bf256_mul(qs_sq[j], qs_dash[j]);   // mul1.mac0
    const bf256_t mul_tag_2 = bf256_mul(qs[j], qs_dash_sq[j]);   // mul2.mac0

    // the constnat term
    zk_hash_256_update(b0_ctx, bf256_add(
                                        mul_tag_1,
                                        bf256_mul(delta_squared, qs[j])
                                        )
    );
    // the constnat term
    zk_hash_256_update(b0_ctx, bf256_add(
                                        mul_tag_2,
                                        bf256_mul(delta_squared, qs_dash[j])
                                        )
    );
    #else
    // the constnat term
    zk_hash_256_update(b0_ctx, bf256_add(
                                        bf256_mul(qs[j], qs_dash[j]),
                                        delta_squared));
    #endif
    // zk_hash_256_update(b0_ctx, bf256_add(bf256_mul(qs[j], qs_dash[j]), delta_squared));
  }
}

/* 
static void aes_prove_256(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* owf_in,
                          const uint8_t* owf_out, const uint8_t* chall_2, uint8_t* a0_tilde,
                          uint8_t* a12_tilde, const faest_paramset_t* params) {
  // Step: 1..2
  bf256_t* bf_v = column_to_row_major_and_shrink_V_256(V, FAEST_256F_ELL);

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7 + 18
  uint8_t* k  = malloc((FAEST_256F_R + 1) * 128 / 8);
  bf256_t* vk = faest_aligned_alloc(BF256_ALIGN, sizeof(bf256_t) * ((FAEST_256F_R + 1) * 128));
  zk_hash_256_ctx a0_ctx;
  zk_hash_256_ctx a1_ctx;

  zk_hash_256_init(&a0_ctx, chall);
  zk_hash_256_init(&a1_ctx, chall);
  aes_key_schedule_constraints_Mkey_0_256(w, bf_v, &a0_ctx, &a1_ctx, k, vk, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  aes_enc_constraints_Mkey_0_256(in, out, w + FAEST_256F_Lke / 8, bf_v + FAEST_256F_Lke, k, vk,
                                 &a0_ctx, &a1_ctx);
  // Step: 12-15
  aes_enc_constraints_Mkey_0_256(in + 16, out + 16, w + (FAEST_256F_Lke + FAEST_256F_Lenc) / 8,
                                 bf_v + FAEST_256F_Lke + FAEST_256F_Lenc, k, vk, &a0_ctx, &a1_ctx);
  faest_aligned_free(vk);
  free(k);

  // Step: 16..18
  zk_hash_256_finalize(a_tilde, &a1_ctx, bf256_load(u + FAEST_256F_ELL / 8));
  zk_hash_256_finalize(b_tilde, &a0_ctx, bf256_sum_poly(bf_v + FAEST_256F_ELL));
  faest_aligned_free(bf_v);
}

static uint8_t* aes_verify_256(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                               const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                               const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.tau0;
  const unsigned int k0  = params->faest_param.k;
  const unsigned int t1  = params->faest_param.tau1;
  const unsigned int k1  = (t0 != 0) ? k0 - 1 : k0;

  // Step: 1
  const uint8_t* delta = chall_3;
  // Step: 2,3
  // do nothing

  // Step: 4..10
  for (unsigned int i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t decoded_challenge[MAX_DEPTH];
    ChalDec(chall_3, i, k0, t0, k1, t1, decoded_challenge);
    for (unsigned int j = 0; j < depth; j++, ++col) {
      if (decoded_challenge[j] == 1) {
        xor_u8_array(d, Q[col], Q[col], (FAEST_256F_ELL + 7) / 8);
      }
    }
  }

  // Step: 11..12
  bf256_t* bf_q = column_to_row_major_and_shrink_V_256(Q, FAEST_256F_ELL);

  // Step: 13, 21
  bf256_t* qk = faest_aligned_alloc(BF256_ALIGN, sizeof(bf256_t) * ((FAEST_256F_R + 1) * 128));
  zk_hash_256_ctx b0_ctx;

  zk_hash_256_init(&b0_ctx, chall_2);
  aes_key_schedule_constraints_Mkey_1_256(bf_q, delta, &b0_ctx, qk);

  // Step: 14
  aes_enc_constraints_Mkey_1_256(in, out, bf_q + FAEST_256F_Lke, qk, delta, &b0_ctx);
  // Step: 18
  aes_enc_constraints_Mkey_1_256(in + 16, out + 16, bf_q + FAEST_256F_Lke + FAEST_256F_Lenc, qk,
                                 delta, &b0_ctx);
  faest_aligned_free(qk);

  // Step: 20, 21
  uint8_t* q_tilde = malloc(FAEST_256F_LAMBDA / 8);
  zk_hash_256_finalize(q_tilde, &b0_ctx, bf256_sum_poly(bf_q + FAEST_256F_ELL));
  faest_aligned_free(bf_q);

  bf256_t bf_qtilde = bf256_load(q_tilde);
  bf256_store(q_tilde, bf256_add(bf_qtilde, bf256_mul(bf256_load(a_tilde), bf256_load(delta))));

  return q_tilde;
}
 */

// ###########################################################################################################################################
// ##################################           LAMBDA = EM-128            ###################################################################
// ###########################################################################################################################################


#if defined(ALLOW_ZERO_SBOX)
static void em_enc_forward_128_1(const uint8_t* z, const uint8_t* x, 
                                bf128_t* bf_y, bf128_t* bf_y_sq)
#else
static void em_enc_forward_128_1(const uint8_t* z, const uint8_t* x, bf128_t* bf_y)
#endif
{
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_128F_Nwd; j++) {

    bf_y[j] = bf128_add(bf128_byte_combine_bits(z[j]), bf128_byte_combine_bits(x[j]));
    #if defined(ALLOW_ZERO_SBOX)
    uint8_t tmp;
    for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
      tmp ^= (  ( (z[j] >> bit_j) ^ (x[j] >> bit_j)  ) & 1 ) << bit_j;
    }
    bf_y_sq[j] = bf128_byte_combine_bits_sq(tmp);
    #endif
  }

  const bf128_t bf_two   = bf128_byte_combine_bits(2);
  const bf128_t bf_three = bf128_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_128F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {

      const unsigned int i  = 32 * FAEST_EM_128F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_128F_Nwd * j + 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_z_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      uint8_t z_bits[4*8];
      uint8_t z_bits_by_two[8];
      #endif

      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          z_bits[r*8 + bit_l] = (z[(i + 8 * r) / 8] >> bit_l) & 1;
        }
        // mul by two
        z_bits_by_two[0] = z_bits[r*8 + 7];
        z_bits_by_two[1] = z_bits[r*8 + 7] ^ z_bits[r*8 + 0];
        z_bits_by_two[2] = z_bits[r*8 + 1];
        z_bits_by_two[3] = z_bits[r*8 + 7] ^ z_bits[r*8 + 2];
        z_bits_by_two[4] = z_bits[r*8 + 7] ^ z_bits[r*8 + 3];
        z_bits_by_two[5] = z_bits[r*8 + 4];
        z_bits_by_two[6] = z_bits[r*8 + 5];
        z_bits_by_two[7] = z_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_z_hat[r] = bf128_byte_combine_bits(z[(i + 8 * r) / 8]);
        bf_x_hat[r] = bf128_byte_combine_bits(x[(i + 8 * r) / 8]);
        #endif
      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        uint8_t output_bits = 0;
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits ^= (z_bits_by_two[r*8 + bit_l] ^ 
                            (z_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (z_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (z_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (z_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ x[(i + 8 * r) / 8])
                                )
                              )
                            )
                          ) << bit_l;
        }
        bf_y[iy + r] = bf128_byte_combine_bits(output_bits);
        bf_y_sq[iy + r] = bf128_byte_combine_bits_sq(output_bits);
      }
      #else
      bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_two), bf128_mul(bf_z_hat[1], bf_three));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_z_hat[2]);
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_z_hat[3]);
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[0]);

      bf_y[iy + 1] = bf128_add(bf_z_hat[0], bf128_mul(bf_z_hat[1], bf_two));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_z_hat[2], bf_three));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_z_hat[3]);
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[1]);

      bf_y[iy + 2] = bf128_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_z_hat[2], bf_two));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_z_hat[3], bf_three));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_x_hat[2]);

      bf_y[iy + 3] = bf128_add(bf128_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_z_hat[2]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_z_hat[3], bf_two));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[3]);
      #endif

    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void em_enc_forward_128(const bf128_t* bf_z, const bf128_t* bf_x, 
                                bf128_t* bf_y, bf128_t* bf_y_sq)
#else
static void em_enc_forward_128(const bf128_t* bf_z, const bf128_t* bf_x, 
                                bf128_t* bf_y)
#endif
{

  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_128F_Nwd; j++) {

    bf_y[j] = bf128_byte_combine(bf_z + 8 * j);
    if (bf_x) {
      bf_y[j] = bf128_add(bf_y[j], bf128_byte_combine(bf_x + 8 * j));
    }

    #if defined(ALLOW_ZERO_SBOX)
    if(bf_x) {
      bf128_t tmp[8];
      for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
        tmp[bit_j] = bf_z[j * 8 + bit_j] ^ bf_x[j * 8 + bit_j];
      }
      bf_y_sq[j] = bf128_byte_combine_sq(tmp);
    } 
    else {
      bf_y_sq[j] = bf128_byte_combine_sq(bf_z + 8 * j);
    }
    #endif
  }

  const bf128_t bf_two   = bf128_byte_combine_bits(2);
  const bf128_t bf_three = bf128_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_128F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {

      const unsigned int i  = 32 * FAEST_EM_128F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_128F_Nwd * j + 4 * c;

      bf128_t bf_z_hat[4];
      bf128_t bf_x_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      bf128_t bf_z_bits[4*8];
      bf128_t bf_z_bits_by_two[8];
      bf128_t output_bits[8];
      #endif
      
      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          bf_z_bits[r*8 + bit_l] = bf_z[(i + 8 * r) + bit_l];
        }
        // mul by two
        bf_z_bits_by_two[0] = bf_z_bits[r*8 + 7];
        bf_z_bits_by_two[1] = bf_z_bits[r*8 + 7] ^ bf_z_bits[r*8 + 0];
        bf_z_bits_by_two[2] = bf_z_bits[r*8 + 1];
        bf_z_bits_by_two[3] = bf_z_bits[r*8 + 7] ^ bf_z_bits[r*8 + 2];
        bf_z_bits_by_two[4] = bf_z_bits[r*8 + 7] ^ bf_z_bits[r*8 + 3];
        bf_z_bits_by_two[5] = bf_z_bits[r*8 + 4];
        bf_z_bits_by_two[6] = bf_z_bits[r*8 + 5];
        bf_z_bits_by_two[7] = bf_z_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_z_hat[r] = bf128_byte_combine(bf_z + (i + 8 * r));
        if (bf_x) {
          bf_x_hat[r] = bf128_byte_combine(bf_x + (i + 8 * r));
        } else {
          bf_x_hat[r] = bf128_zero();
        }
        #endif
      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        bf128_t output_bits[8];
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits[bit_l] = (bf_z_bits_by_two[r*8 + bit_l] ^ 
                            (bf_z_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (bf_z_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (bf_z_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (bf_z_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ bf_x[(i + 8 * r) / 8])
                                )
                              )
                            )
                          );
        }
        bf_y[iy + r] = bf128_byte_combine(output_bits);
        bf_y_sq[iy + r] = bf128_byte_combine_sq(output_bits);
      }
      #else
      bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_two), bf128_mul(bf_z_hat[1], bf_three));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_z_hat[2]);
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_z_hat[3]);
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[0]);

      bf_y[iy + 1] = bf128_add(bf_z_hat[0], bf128_mul(bf_z_hat[1], bf_two));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_z_hat[2], bf_three));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_z_hat[3]);
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[1]);

      bf_y[iy + 2] = bf128_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_z_hat[2], bf_two));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_z_hat[3], bf_three));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_x_hat[2]);

      bf_y[iy + 3] = bf128_add(bf128_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_z_hat[2]);
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_z_hat[3], bf_two));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[3]);
      #endif
    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void em_enc_backward_128_1(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf128_t* y_out, bf128_t* y_out_sq)
#else
static void em_enc_backward_128_1(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf128_t* y_out)
#endif
{
  // only called with Mtag == Mkey == 0

  for (unsigned int j = 0; j < FAEST_EM_128F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {

      for (unsigned int r = 0; r < 4; r++) {

        const unsigned int icol = (c - r + FAEST_EM_128F_Nwd) % FAEST_EM_128F_Nwd;
        const unsigned int ird =
            FAEST_EM_128F_LAMBDA + 32 * FAEST_EM_128F_Nwd * j + 32 * icol + 8 * r;
        uint8_t z_tilde = 0;
        if (j < (FAEST_EM_128F_R - 1)) {
          z_tilde = z[ird / 8];
        } else {
          z_tilde = z_out[(ird - 32 * FAEST_EM_128F_Nwd * (j + 1)) / 8] ^ x[ird / 8];
        }

        // (bit spliced)
        // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
        const uint8_t y_tilde = rotr8(z_tilde, 7) ^ rotr8(z_tilde, 5) ^ rotr8(z_tilde, 2) ^ 0x5;

        // Step: 18
        y_out[4 * FAEST_EM_128F_Nwd * j + 4 * c + r] = bf128_byte_combine_bits(y_tilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[4 * FAEST_EM_128F_Nwd * j + 4 * c + r] = bf128_byte_combine_bits_sq(y_tilde);
        #endif
      }
    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void em_enc_backward_128(const bf128_t* bf_z, const bf128_t* bf_x, const bf128_t* bf_z_out,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, 
                                bf128_t* y_out, bf128_t* y_out_sq)
#else
static void em_enc_backward_128(const bf128_t* bf_z, const bf128_t* bf_x, const bf128_t* bf_z_out,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, 
                                bf128_t* y_out)
#endif
{
  // Step: 1
  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t factor =
      bf128_mul_bit(bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int j = 0; j < FAEST_EM_128F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {

      for (unsigned int r = 0; r < 4; r++) {

        bf128_t bf_z_tilde[8];
        const unsigned int icol = (c - r + FAEST_EM_128F_Nwd) % FAEST_EM_128F_Nwd;
        const unsigned int ird =
            FAEST_EM_128F_LAMBDA + 32 * FAEST_EM_128F_Nwd * j + 32 * icol + 8 * r;

        if (j < (FAEST_EM_128F_R - 1)) {
          memcpy(bf_z_tilde, bf_z + ird, sizeof(bf_z_tilde));
        } else {
          for (unsigned int i = 0; i < 8; ++i) {
            // Step: 12
            bf_z_tilde[i] = bf_z_out[ird - 32 * FAEST_EM_128F_Nwd * (j + 1) + i];
            if (bf_x) {
              bf_z_tilde[i] = bf128_add(bf_z_tilde[i], bf_x[ird + i]);
            }
          }
        }

        bf128_t bf_y_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf128_add(bf128_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                    bf_z_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

        // Step: 18
        y_out[4 * FAEST_EM_128F_Nwd * j + 4 * c + r] = bf128_byte_combine(bf_y_tilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[4 * FAEST_EM_128F_Nwd * j + 4 * c + r] = bf128_byte_combine_sq(bf_y_tilde);
        #endif
      }
    }
  }
}

// Mkey = 0, this is for the prover
static void em_enc_constraints_Mkey_0_128(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                          const bf128_t* bf_v, zk_hash_128_ctx* a0_ctx,
                                          zk_hash_128_ctx* a1_ctx) {
  // Step 6
  uint8_t w_out[FAEST_EM_128F_LAMBDA / 8];
  xor_u8_array(out, w, w_out, sizeof(w_out));

  bf128_t bf_s[FAEST_EM_128F_Senc]; // inverse input
  bf128_t bf_vs[FAEST_EM_128F_Senc]; // inverse input tag
  #if defined(ALLOW_ZERO_SBOX)
  bf128_t bf_s_sq[FAEST_EM_128F_Senc]; // inverse input sq
  bf128_t bf_vs_sq[FAEST_EM_128F_Senc]; // input input tag sq
  em_enc_forward_128_1(w, x, bf_s, bf_s_sq);
  em_enc_forward_128(bf_v, NULL, bf_vs, bf_vs_sq);
  #else
  em_enc_forward_128_1(w, x, bf_s);
  em_enc_forward_128(bf_v, NULL, bf_vs);
  #endif

  bf128_t bf_s_dash[FAEST_EM_128F_Senc]; // inverse output
  bf128_t bf_vs_dash[FAEST_EM_128F_Senc]; // inverse output tag
  #if defined(ALLOW_ZERO_SBOX)
  bf128_t bf_s_dash_sq[FAEST_EM_128F_Senc]; // inverse output sq
  bf128_t bf_vs_dash_sq[FAEST_EM_128F_Senc]; // inverse output tag sq
  em_enc_backward_128_1(w, x, w_out, bf_s_dash, bf_s_dash_sq);
  em_enc_backward_128(bf_v, NULL, bf_v, 1, 0, NULL, bf_vs_dash, bf_vs_dash_sq);
  #else
  em_enc_backward_128_1(w, x, w_out, bf_s_dash);
  em_enc_backward_128(bf_v, NULL, bf_v, 1, 0, NULL, bf_vs_dash);
  #endif  

  for (unsigned int j = 0; j < FAEST_EM_128F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Psuedoinverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf128_t mul_tag_1 = bf128_mul(bf_vs_sq[j], bf_vs_dash[j]);   // mul1.mac0
    const bf128_t mul_tag_2 = bf128_mul(bf_vs[j], bf_vs_dash_sq[j]);   // mul2.mac0
    // add the tags with the values and multiply the results
    const bf128_t mul_val_1 = bf128_mul(                        // mul1.mac1
                        bf128_add(bf_s_sq[j], bf_vs_sq[j]), 
                        bf128_add(bf_s_dash[j], bf_vs_dash[j])
                        );
    const bf128_t mul_val_2 = bf128_mul(                        // mul2.mac1
                        bf128_add(bf_s[j], bf_vs[j]),
                        bf128_add(bf_s_dash_sq[j], bf_vs_dash_sq[j])
                        );

    // the constnat term
    zk_hash_128_update(a0_ctx, mul_tag_1);
    // the linear term
    zk_hash_128_update(a1_ctx, bf128_add(
                                        bf128_add(mul_tag_1, bf_s[j]),   // mul1.mac0 + x.val  (because here we check x^2*y = x)
                                        bf128_add(mul_val_1, bf_vs[j])   // mul1.mac1 + x.mac
                                        )
                      );
    // the constant term
    zk_hash_128_update(a0_ctx, mul_tag_2);
    // the linear term
    zk_hash_128_update(a1_ctx, bf128_add(
                                        bf128_add(mul_tag_2, bf_s_dash[j]),  // mul2.mac0 + y.val  (because here we check x*y^2 = y)
                                        bf128_add(mul_val_2, bf_vs_dash[j])  // mul2.mac1 + y.mac
                                        )
                      );    
    #else
    // Inverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf128_t mul_tag = bf128_mul(bf_vs[j], bf_vs_dash[j]); // mac0
    // add the tags with the values and multiply the results
    const bf128_t mul_val = bf128_mul(bf128_add(bf_s[j], bf_vs[j]), bf128_add(bf_s_dash[j], bf_vs_dash[j])); // mac1
    
    // this is the constant term
    zk_hash_128_update(a0_ctx, mul_tag);
    // this is the linear term
    zk_hash_128_update(a1_ctx, bf128_add(
                                        bf128_add(mul_tag, bf128_one),    // mac0 + 1     (because here we check x*y = 1)
                                        mul_val                             // (mac0 + 1) + mac1
                                        )
                      );
    #endif

    // const bf128_t tmp = bf128_mul(bf_vs[j], bf_vs_dash[j]);
    // zk_hash_128_update(a0_ctx, tmp);
    // zk_hash_128_update(a1_ctx,
    //                    bf128_add(bf128_add(bf128_mul(bf128_add(bf_s[j], bf_vs[j]),
    //                                                  bf128_add(bf_s_dash[j], bf_vs_dash[j])),
    //                                        tmp),
    //                              bf128_one()));
  }
}

// Mkey = 1, this is for the verifier
static void em_enc_constraints_Mkey_1_128(const uint8_t* out, const uint8_t* x, const bf128_t* bf_q,
                                          const uint8_t* delta, zk_hash_128_ctx* b0_ctx) {
  
  // Step: 18, 19
  const bf128_t bf_delta = bf128_load(delta);
  bf128_t* bf_x = faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * 128 * (FAEST_EM_128F_R + 1));
  for (unsigned int i = 0; i < 128 * (FAEST_EM_128F_R + 1); i++) {
    bf_x[i] = bf128_mul_bit(bf_delta, ptr_get_bit(x, i));
  }

  // Step 21
  bf128_t* bf_q_out = faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * FAEST_EM_128F_LAMBDA);
  for (unsigned int i = 0; i < FAEST_EM_128F_LAMBDA; i++) {
    bf_q_out[i] = bf128_add(bf128_mul_bit(bf_delta, ptr_get_bit(out, i)), bf_q[i]);
  }

  bf128_t bf_qs[FAEST_EM_128F_Senc];
  bf128_t bf_qs_dash[FAEST_EM_128F_Senc];
  #if defined(ALLOW_ZERO_SBOX)
  bf128_t bf_qs_sq[FAEST_EM_128F_Senc];
  bf128_t bf_qs_dash_sq[FAEST_EM_128F_Senc];
  em_enc_forward_128(bf_q, bf_x, bf_qs, bf_qs_sq);
  em_enc_backward_128(bf_q, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash, bf_qs_dash_sq);
  #else
  em_enc_forward_128(bf_q, bf_x, bf_qs);
  em_enc_backward_128(bf_q, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash);
  #endif
  faest_aligned_free(bf_q_out);
  faest_aligned_free(bf_x);

  // Step: 13..14
  bf128_t delta_squared = bf128_mul(bf_delta, bf_delta);

  for (unsigned int j = 0; j < FAEST_EM_128F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Quicksilver multiplication
    // multiply the tags
    const bf128_t mul_tag_1 = bf128_mul(bf_qs_sq[j], bf_qs_dash[j]);   // mul1.mac0
    const bf128_t mul_tag_2 = bf128_mul(bf_qs[j], bf_qs_dash_sq[j]);   // mul2.mac0

    // the constnat term
    zk_hash_128_update(b0_ctx, bf128_add(
                                        mul_tag_1,
                                        bf128_mul(delta_squared, bf_qs[j])
                                        )
    );
    // the constnat term
    zk_hash_128_update(b0_ctx, bf128_add(
                                        mul_tag_2,
                                        bf128_mul(delta_squared, bf_qs_dash[j])
                                        )
    );
    #else
    // the constant term
    zk_hash_128_update(b0_ctx, bf128_add(
                                        bf128_mul(bf_qs[j], bf_qs_dash[j]), 
                                        delta_squared));
    #endif
  }
}

/* 
static void em_prove_128(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* owf_in,
                          const uint8_t* owf_out, const uint8_t* chall_2, uint8_t* a0_tilde,
                          uint8_t* a12_tilde, const faest_paramset_t* params) {
  // copy expanded key in to an array
  uint8_t x[FAEST_EM_128F_LAMBDA * (FAEST_EM_128F_R + 1)];
  {
    aes_round_keys_t round_keys;
    aes128_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_128F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_128F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  bf128_t* bf_v = column_to_row_major_and_shrink_V_128(V, FAEST_EM_128F_Lenc);
  zk_hash_128_ctx a0_ctx;
  zk_hash_128_ctx a1_ctx;

  zk_hash_128_init(&a0_ctx, chall);
  zk_hash_128_init(&a1_ctx, chall);
  em_enc_constraints_Mkey_0_128(out, x, w, bf_v, &a0_ctx, &a1_ctx);

  zk_hash_128_finalize(a_tilde, &a1_ctx, bf128_load(u + FAEST_EM_128F_Lenc / 8));
  zk_hash_128_finalize(b_tilde, &a0_ctx, bf128_sum_poly(bf_v + FAEST_EM_128F_Lenc));

  faest_aligned_free(bf_v);
}

static uint8_t* em_verify_128(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                              const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                              const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.tau0;
  const unsigned int k0  = params->faest_param.k;
  const unsigned int t1  = params->faest_param.tau1;
  const unsigned int k1  = (t0 != 0) ? k0 - 1 : k0;

  const uint8_t* delta = chall_3;

  for (unsigned int i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t decoded_challenge[MAX_DEPTH];
    ChalDec(chall_3, i, k0, t0, k1, t1, decoded_challenge);
    for (unsigned int j = 0; j < depth; j++, ++col) {
      if (decoded_challenge[j] == 1) {
        xor_u8_array(d, Q[col], Q[col], (FAEST_EM_128F_Lenc + 7) / 8);
      }
    }
  }

  bf128_t* bf_q = column_to_row_major_and_shrink_V_128(Q, FAEST_EM_128F_Lenc);

  // copy expanded key in to an array
  uint8_t x[FAEST_EM_128F_LAMBDA * (FAEST_EM_128F_R + 1) / 8];
  {
    aes_round_keys_t round_keys;
    aes128_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_128F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_128F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  zk_hash_128_ctx b0_ctx;
  zk_hash_128_init(&b0_ctx, chall_2);
  em_enc_constraints_Mkey_1_128(out, x, bf_q, delta, &b0_ctx);

  uint8_t* q_tilde = malloc(FAEST_EM_128F_LAMBDA / 8);
  zk_hash_128_finalize(q_tilde, &b0_ctx, bf128_sum_poly(bf_q + FAEST_EM_128F_Lenc));
  faest_aligned_free(bf_q);

  bf128_t bf_qtilde = bf128_load(q_tilde);
  bf128_store(q_tilde, bf128_add(bf_qtilde, bf128_mul(bf128_load(a_tilde), bf128_load(delta))));

  return q_tilde;
}

 */

// ###########################################################################################################################################
// ##################################           LAMBDA = EM-192            ###################################################################
// ###########################################################################################################################################

#if defined(ALLOW_ZERO_SBOX)
static void em_enc_forward_192_1(const uint8_t* z, const uint8_t* x, 
                                  bf192_t* bf_y, bf192_t* bf_y_sq)
#else
static void em_enc_forward_192_1(const uint8_t* z, const uint8_t* x, bf192_t* bf_y)
#endif
{
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_192F_Nwd; j++) {

    bf_y[j] = bf192_add(bf192_byte_combine_bits(z[j]), bf192_byte_combine_bits(x[j]));
    #if defined(ALLOW_ZERO_SBOX)
    uint8_t tmp;
    for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
      tmp ^= (  ( (z[j] >> bit_j) ^ (x[j] >> bit_j)  ) & 1 ) << bit_j;
    }
    bf_y_sq[j] = bf192_byte_combine_bits_sq(tmp);
    #endif
  }

  const bf192_t bf_two   = bf192_byte_combine_bits(2);
  const bf192_t bf_three = bf192_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_192F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {

      const unsigned int i  = 32 * FAEST_EM_192F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_192F_Nwd * j + 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_z_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      uint8_t z_bits[4*8];
      uint8_t z_bits_by_two[8];
      #endif

      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          z_bits[r*8 + bit_l] = (z[(i + 8 * r) / 8] >> bit_l) & 1;
        }
        // mul by two
        z_bits_by_two[0] = z_bits[r*8 + 7];
        z_bits_by_two[1] = z_bits[r*8 + 7] ^ z_bits[r*8 + 0];
        z_bits_by_two[2] = z_bits[r*8 + 1];
        z_bits_by_two[3] = z_bits[r*8 + 7] ^ z_bits[r*8 + 2];
        z_bits_by_two[4] = z_bits[r*8 + 7] ^ z_bits[r*8 + 3];
        z_bits_by_two[5] = z_bits[r*8 + 4];
        z_bits_by_two[6] = z_bits[r*8 + 5];
        z_bits_by_two[7] = z_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_z_hat[r] = bf192_byte_combine_bits(z[(i + 8 * r) / 8]);
        bf_x_hat[r] = bf192_byte_combine_bits(x[(i + 8 * r) / 8]);
        #endif
      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        uint8_t output_bits = 0;
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits ^= (z_bits_by_two[r*8 + bit_l] ^ 
                            (z_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (z_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (z_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (z_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ x[(i + 8 * r) / 8])
                                )
                              )
                            )
                          ) << bit_l;
        }
        bf_y[iy + r] = bf192_byte_combine_bits(output_bits);
        bf_y_sq[iy + r] = bf192_byte_combine_bits_sq(output_bits);
      }
      #else
      bf_y[iy + 0] = bf192_add(bf192_mul(bf_z_hat[0], bf_two), bf192_mul(bf_z_hat[1], bf_three));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_z_hat[2]);
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_z_hat[3]);
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_x_hat[0]);

      bf_y[iy + 1] = bf192_add(bf_z_hat[0], bf192_mul(bf_z_hat[1], bf_two));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_z_hat[2], bf_three));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf_z_hat[3]);
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf_x_hat[1]);

      bf_y[iy + 2] = bf192_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_z_hat[2], bf_two));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_z_hat[3], bf_three));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf_x_hat[2]);

      bf_y[iy + 3] = bf192_add(bf192_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_z_hat[2]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_z_hat[3], bf_two));
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_x_hat[3]);
      #endif
      
    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void em_enc_forward_192(const bf192_t* bf_z, const bf192_t* bf_x, 
                              bf192_t* bf_y, bf192_t* bf_y_sq)
#else
static void em_enc_forward_192(const bf192_t* bf_z, const bf192_t* bf_x, 
                                bf192_t* bf_y)
#endif
{

  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_192F_Nwd; j++) {

    bf_y[j] = bf192_byte_combine(bf_z + 8 * j);
    if (bf_x) {
      bf_y[j] = bf192_add(bf_y[j], bf192_byte_combine(bf_x + 8 * j));
    }

    #if defined(ALLOW_ZERO_SBOX)
    if(bf_x) {
      bf192_t tmp[8];
      for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
        tmp[bit_j] = bf_z[j * 8 + bit_j] ^ bf_x[j * 8 + bit_j];
      }
      bf_y_sq[j] = bf192_byte_combine_sq(tmp);
    } 
    else {
      bf_y_sq[j] = bf192_byte_combine_sq(bf_z + 8 * j);
    }
    #endif
  }

  const bf192_t bf_two   = bf192_byte_combine_bits(2);
  const bf192_t bf_three = bf192_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_192F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {

      const unsigned int i  = 32 * FAEST_EM_192F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_192F_Nwd * j + 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_z_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      bf192_t bf_z_bits[4*8];
      bf192_t bf_z_bits_by_two[8];
      bf192_t output_bits[8];
      #endif

      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          bf_z_bits[r*8 + bit_l] = bf_z[(i + 8 * r) + bit_l];
        }
        // mul by two
        bf_z_bits_by_two[0] = bf_z_bits[r*8 + 7];
        bf_z_bits_by_two[1] = bf_z_bits[r*8 + 7] ^ bf_z_bits[r*8 + 0];
        bf_z_bits_by_two[2] = bf_z_bits[r*8 + 1];
        bf_z_bits_by_two[3] = bf_z_bits[r*8 + 7] ^ bf_z_bits[r*8 + 2];
        bf_z_bits_by_two[4] = bf_z_bits[r*8 + 7] ^ bf_z_bits[r*8 + 3];
        bf_z_bits_by_two[5] = bf_z_bits[r*8 + 4];
        bf_z_bits_by_two[6] = bf_z_bits[r*8 + 5];
        bf_z_bits_by_two[7] = bf_z_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_z_hat[r] = bf192_byte_combine(bf_z + (i + 8 * r));
        if (bf_x) {
          bf_x_hat[r] = bf192_byte_combine(bf_x + (i + 8 * r));
        } else {
          bf_x_hat[r] = bf192_zero();
        }
        #endif
      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        bf192_t output_bits[8];
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits[bit_l] = (bf_z_bits_by_two[r*8 + bit_l] ^ 
                            (bf_z_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (bf_z_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (bf_z_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (bf_z_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ bf_x[(i + 8 * r) / 8])
                                )
                              )
                            )
                          );
        }
        bf_y[iy + r] = bf192_byte_combine(output_bits);
        bf_y_sq[iy + r] = bf192_byte_combine_sq(output_bits);
      }
      #else
      bf_y[iy + 0] = bf192_add(bf192_mul(bf_z_hat[0], bf_two), bf192_mul(bf_z_hat[1], bf_three));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_z_hat[2]);
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_z_hat[3]);
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf_x_hat[0]);

      bf_y[iy + 1] = bf192_add(bf_z_hat[0], bf192_mul(bf_z_hat[1], bf_two));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_z_hat[2], bf_three));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf_z_hat[3]);
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf_x_hat[1]);

      bf_y[iy + 2] = bf192_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_z_hat[2], bf_two));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_z_hat[3], bf_three));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf_x_hat[2]);

      bf_y[iy + 3] = bf192_add(bf192_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_z_hat[2]);
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_z_hat[3], bf_two));
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf_x_hat[3]);
      #endif
    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void em_enc_backward_192_1(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf192_t* y_out, bf192_t* y_out_sq)
#else
static void em_enc_backward_192_1(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf192_t* y_out)
#endif
{
  // only called with Mtag == Mkey == 0

  for (unsigned int j = 0; j < FAEST_EM_192F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {

      for (unsigned int r = 0; r < 4; r++) {

        const unsigned int icol = (c - r + FAEST_EM_192F_Nwd) % FAEST_EM_192F_Nwd;
        const unsigned int ird =
            FAEST_EM_192F_LAMBDA + 32 * FAEST_EM_192F_Nwd * j + 32 * icol + 8 * r;
        uint8_t z_tilde = 0;
        if (j < (FAEST_EM_192F_R - 1)) {
          z_tilde = z[ird / 8];
        } else {
          z_tilde = z_out[(ird - 32 * FAEST_EM_192F_Nwd * (j + 1)) / 8] ^ x[ird / 8];
        }

        // (bit spliced)
        // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
        const uint8_t y_tilde = rotr8(z_tilde, 7) ^ rotr8(z_tilde, 5) ^ rotr8(z_tilde, 2) ^ 0x5;

        // Step: 18
        y_out[4 * FAEST_EM_192F_Nwd * j + 4 * c + r] = bf192_byte_combine_bits(y_tilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[4 * FAEST_EM_192F_Nwd * j + 4 * c + r] = bf192_byte_combine_bits_sq(y_tilde);
        #endif
      }
    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void em_enc_backward_192(const bf192_t* bf_z, const bf192_t* bf_x, const bf192_t* bf_z_out,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, 
                                bf192_t* y_out, bf192_t* y_out_sq)
#else
static void em_enc_backward_192(const bf192_t* bf_z, const bf192_t* bf_x, const bf192_t* bf_z_out,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, 
                                bf192_t* y_out, bf192_t* y_out_sq)
#endif
{
  // Step: 1
  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t factor =
      bf192_mul_bit(bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int j = 0; j < FAEST_EM_192F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {

      for (unsigned int r = 0; r < 4; r++) {

        bf192_t bf_z_tilde[8];
        const unsigned int icol = (c - r + FAEST_EM_192F_Nwd) % FAEST_EM_192F_Nwd;
        const unsigned int ird =
            FAEST_EM_192F_LAMBDA + 32 * FAEST_EM_192F_Nwd * j + 32 * icol + 8 * r;

        if (j < (FAEST_EM_192F_R - 1)) {
          memcpy(bf_z_tilde, bf_z + ird, sizeof(bf_z_tilde));
        } else {
          for (unsigned int i = 0; i < 8; ++i) {
            // Step: 12
            bf_z_tilde[i] = bf_z_out[ird - 32 * FAEST_EM_192F_Nwd * (j + 1) + i];
            if (bf_x) {
              bf_z_tilde[i] = bf192_add(bf_z_tilde[i], bf_x[ird + i]);
            }
          }
        }

        bf192_t bf_y_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf192_add(bf192_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                    bf_z_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);

        // Step: 18
        y_out[4 * FAEST_EM_192F_Nwd * j + 4 * c + r] = bf192_byte_combine(bf_y_tilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[4 * FAEST_EM_192F_Nwd * j + 4 * c + r] = bf192_byte_combine_sq(bf_y_tilde);
        #endif
      }
    }
  }
}

// Mkey = 0, this is for the prover
static void em_enc_constraints_Mkey_0_192(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                          const bf192_t* bf_v, zk_hash_192_ctx* a0_ctx,
                                          zk_hash_192_ctx* a1_ctx) {
  // Step 6
  uint8_t w_out[FAEST_EM_192F_LAMBDA / 8];
  xor_u8_array(out, w, w_out, sizeof(w_out));

  bf192_t bf_s[FAEST_EM_192F_Senc]; // inverse input
  bf192_t bf_vs[FAEST_EM_192F_Senc]; // inverse input tag
  #if defined(ALLOW_ZERO_SBOX)
  bf192_t bf_s_sq[FAEST_EM_192F_Senc]; // inverse input sq
  bf192_t bf_vs_sq[FAEST_EM_192F_Senc]; // input input tag sq
  em_enc_forward_192_1(w, x, bf_s, bf_s_sq);
  em_enc_forward_192(bf_v, NULL, bf_vs, bf_vs_sq);
  #else
  em_enc_forward_192_1(w, x, bf_s);
  em_enc_forward_192(bf_v, NULL, bf_vs);
  #endif


  bf192_t bf_s_dash[FAEST_EM_192F_Senc]; // inverse output
  bf192_t bf_vs_dash[FAEST_EM_192F_Senc]; // inverse output tag
  #if defined(ALLOW_ZERO_SBOX)
  bf192_t bf_s_dash_sq[FAEST_EM_192F_Senc]; // inverse output sq
  bf192_t bf_vs_dash_sq[FAEST_EM_192F_Senc]; // inverse output tag sq
  em_enc_backward_192_1(w, x, w_out, bf_s_dash, bf_s_dash_sq);
  em_enc_backward_192(bf_v, NULL, bf_v, 1, 0, NULL, bf_vs_dash, bf_vs_dash_sq);
  #else
  em_enc_backward_192_1(w, x, w_out, bf_s_dash);
  em_enc_backward_192(bf_v, NULL, bf_v, 1, 0, NULL, bf_vs_dash);
  #endif  
  
  for (unsigned int j = 0; j < FAEST_EM_192F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Psuedoinverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf192_t mul_tag_1 = bf192_mul(bf_vs_sq[j], bf_vs_dash[j]);   // mul1.mac0
    const bf192_t mul_tag_2 = bf192_mul(bf_vs[j], bf_vs_dash_sq[j]);   // mul2.mac0
    // add the tags with the values and multiply the results
    const bf192_t mul_val_1 = bf192_mul(                        // mul1.mac1
                        bf192_add(bf_s_sq[j], bf_vs_sq[j]), 
                        bf192_add(bf_s_dash[j], bf_vs_dash[j])
                        );
    const bf192_t mul_val_2 = bf192_mul(                        // mul2.mac1
                        bf192_add(bf_s[j], bf_vs[j]),
                        bf192_add(bf_s_dash_sq[j], bf_vs_dash_sq[j])
                        );

    // the constnat term
    zk_hash_192_update(a0_ctx, mul_tag_1);
    // the linear term
    zk_hash_192_update(a1_ctx, bf192_add(
                                        bf192_add(mul_tag_1, bf_s[j]),   // mul1.mac0 + x.val  (because here we check x^2*y = x)
                                        bf192_add(mul_val_1, bf_vs[j])   // mul1.mac1 + x.mac
                                        )
                      );
    // the constant term
    zk_hash_192_update(a0_ctx, mul_tag_2);
    // the linear term
    zk_hash_192_update(a1_ctx, bf192_add(
                                        bf192_add(mul_tag_2, bf_s_dash[j]),  // mul2.mac0 + y.val  (because here we check x*y^2 = y)
                                        bf192_add(mul_val_2, bf_vs_dash[j])  // mul2.mac1 + y.mac
                                        )
                      );    
    #else
    // Inverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf192_t mul_tag = bf192_mul(bf_vs[j], bf_vs_dash[j]); // mac0
    // add the tags with the values and multiply the results
    const bf192_t mul_val = bf192_mul(bf192_add(bf_s[j], bf_vs[j]), bf192_add(bf_s_dash[j], bf_vs_dash[j])); // mac1
    
    // this is the constant term
    zk_hash_192_update(a0_ctx, mul_tag);
    // this is the linear term
    zk_hash_192_update(a1_ctx, bf192_add(
                                        bf192_add(mul_tag, bf192_one),    // mac0 + 1     (because here we check x*y = 1)
                                        mul_val                             // (mac0 + 1) + mac1
                                        )
                      );
    #endif
  }
}

// Mkey = 1, this is for the verifier
static void em_enc_constraints_Mkey_1_192(const uint8_t* out, const uint8_t* x, const bf192_t* bf_q,
                                          const uint8_t* delta, zk_hash_192_ctx* b0_ctx) {
  // Step: 18, 19
  const bf192_t bf_delta = bf192_load(delta);
  bf192_t* bf_x = faest_aligned_alloc(BF192_ALIGN, sizeof(bf192_t) * 192 * (FAEST_EM_192F_R + 1));
  for (unsigned int i = 0; i < 192 * (FAEST_EM_192F_R + 1); i++) {
    bf_x[i] = bf192_mul_bit(bf_delta, ptr_get_bit(x, i));
  }

  // Step 21
  bf192_t* bf_q_out = faest_aligned_alloc(BF192_ALIGN, sizeof(bf192_t) * FAEST_EM_192F_LAMBDA);
  for (unsigned int i = 0; i < FAEST_EM_192F_LAMBDA; i++) {
    bf_q_out[i] = bf192_add(bf192_mul_bit(bf_delta, ptr_get_bit(out, i)), bf_q[i]);
  }

  bf192_t bf_qs[FAEST_EM_192F_Senc];
  bf192_t bf_qs_dash[FAEST_EM_192F_Senc];
  #if defined(ALLOW_ZERO_SBOX)
  bf192_t bf_qs_sq[FAEST_EM_192F_Senc];
  bf192_t bf_qs_dash_sq[FAEST_EM_192F_Senc];
  em_enc_forward_192(bf_q, bf_x, bf_qs, bf_qs_sq);
  em_enc_backward_192(bf_q, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash, bf_qs_dash_sq);
  #else
  em_enc_forward_192(bf_q, bf_x, bf_qs);
  em_enc_backward_192(bf_q, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash);
  #endif
  faest_aligned_free(bf_q_out);
  faest_aligned_free(bf_x);

  // Step: 13..14
  bf192_t delta_squared = bf192_mul(bf_delta, bf_delta);

  for (unsigned int j = 0; j < FAEST_EM_192F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Quicksilver multiplication
    // multiply the tags
    const bf192_t mul_tag_1 = bf192_mul(bf_qs_sq[j], bf_qs_dash[j]);   // mul1.mac0
    const bf192_t mul_tag_2 = bf192_mul(bf_qs[j], bf_qs_dash_sq[j]);   // mul2.mac0

    // the constnat term
    zk_hash_192_update(b0_ctx, bf192_add(
                                        mul_tag_1,
                                        bf192_mul(delta_squared, bf_qs[j])
                                        )
    );
    // the constnat term
    zk_hash_192_update(b0_ctx, bf192_add(
                                        mul_tag_2,
                                        bf192_mul(delta_squared, bf_qs_dash[j])
                                        )
    );
    #else
    // the constant term
    zk_hash_192_update(b0_ctx, bf192_add(
                                        bf192_mul(bf_qs[j], bf_qs_dash[j]), 
                                        delta_squared));
    #endif
  }
}

/* 

static void em_prove_192(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* owf_in,
                          const uint8_t* owf_out, const uint8_t* chall_2, uint8_t* a0_tilde,
                          uint8_t* a12_tilde, const faest_paramset_t* params) {
  // copy expanded key in to an array
  uint8_t x[FAEST_EM_192F_LAMBDA * (FAEST_EM_192F_R + 1) / 8];
  {
    aes_round_keys_t round_keys;
    rijndael192_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_192F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_192F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  bf192_t* bf_v = column_to_row_major_and_shrink_V_192(V, FAEST_EM_192F_Lenc);
  zk_hash_192_ctx a0_ctx;
  zk_hash_192_ctx a1_ctx;

  zk_hash_192_init(&a0_ctx, chall);
  zk_hash_192_init(&a1_ctx, chall);
  em_enc_constraints_Mkey_0_192(out, x, w, bf_v, &a0_ctx, &a1_ctx);

  zk_hash_192_finalize(a_tilde, &a1_ctx, bf192_load(u + FAEST_EM_192F_Lenc / 8));
  zk_hash_192_finalize(b_tilde, &a0_ctx, bf192_sum_poly(bf_v + FAEST_EM_192F_Lenc));

  faest_aligned_free(bf_v);
}

static uint8_t* em_verify_192(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                              const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                              const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.tau0;
  const unsigned int k0  = params->faest_param.k;
  const unsigned int t1  = params->faest_param.tau1;
  const unsigned int k1  = (t0 != 0) ? k0 - 1 : k0;

  const uint8_t* delta = chall_3;

  for (unsigned int i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t decoded_challenge[MAX_DEPTH];
    ChalDec(chall_3, i, k0, t0, k1, t1, decoded_challenge);
    for (unsigned int j = 0; j < depth; j++, ++col) {
      if (decoded_challenge[j] == 1) {
        xor_u8_array(d, Q[col], Q[col], (FAEST_EM_192F_Lenc + 7) / 8);
      }
    }
  }

  bf192_t* bf_q = column_to_row_major_and_shrink_V_192(Q, FAEST_EM_192F_Lenc);

  // copy expanded key in to an array
  uint8_t x[FAEST_EM_192F_LAMBDA * (FAEST_EM_192F_R + 1) / 8];
  {
    aes_round_keys_t round_keys;
    rijndael192_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_192F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_192F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  zk_hash_192_ctx b0_ctx;
  zk_hash_192_init(&b0_ctx, chall_2);
  em_enc_constraints_Mkey_1_192(out, x, bf_q, delta, &b0_ctx);

  uint8_t* q_tilde = malloc(FAEST_EM_192F_LAMBDA / 8);
  zk_hash_192_finalize(q_tilde, &b0_ctx, bf192_sum_poly(bf_q + FAEST_EM_192F_Lenc));
  faest_aligned_free(bf_q);

  bf192_t bf_qtilde = bf192_load(q_tilde);
  bf192_store(q_tilde, bf192_add(bf_qtilde, bf192_mul(bf192_load(a_tilde), bf192_load(delta))));

  return q_tilde;
}

 */

// ###########################################################################################################################################
// ##################################           LAMBDA = EM-256            ###################################################################
// ###########################################################################################################################################

#if defined(ALLOW_ZERO_SBOX)
static void em_enc_forward_256_1(const uint8_t* z, const uint8_t* x, 
                                  bf256_t* bf_y, bf256_t* bf_y_sq)
#else
static void em_enc_forward_256_1(const uint8_t* z, const uint8_t* x, bf256_t* bf_y)
#endif
{
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_256F_Nwd; j++) {

    bf_y[j] = bf256_add(bf256_byte_combine_bits(z[j]), bf256_byte_combine_bits(x[j]));
    #if defined(ALLOW_ZERO_SBOX)
    uint8_t tmp;
    for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
      tmp ^= (  ( (z[j] >> bit_j) ^ (x[j] >> bit_j)  ) & 1 ) << bit_j;
    }
    bf_y_sq[j] = bf256_byte_combine_bits_sq(tmp);
    #endif
  }

  const bf256_t bf_two   = bf256_byte_combine_bits(2);
  const bf256_t bf_three = bf256_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_256F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {

      const unsigned int i  = 32 * FAEST_EM_256F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_256F_Nwd * j + 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_z_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      uint8_t z_bits[4*8];
      uint8_t z_bits_by_two[8];
      #endif

      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          z_bits[r*8 + bit_l] = (z[(i + 8 * r) / 8] >> bit_l) & 1;
        }
        // mul by two
        z_bits_by_two[0] = z_bits[r*8 + 7];
        z_bits_by_two[1] = z_bits[r*8 + 7] ^ z_bits[r*8 + 0];
        z_bits_by_two[2] = z_bits[r*8 + 1];
        z_bits_by_two[3] = z_bits[r*8 + 7] ^ z_bits[r*8 + 2];
        z_bits_by_two[4] = z_bits[r*8 + 7] ^ z_bits[r*8 + 3];
        z_bits_by_two[5] = z_bits[r*8 + 4];
        z_bits_by_two[6] = z_bits[r*8 + 5];
        z_bits_by_two[7] = z_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_z_hat[r] = bf256_byte_combine_bits(z[(i + 8 * r) / 8]);
        bf_x_hat[r] = bf256_byte_combine_bits(x[(i + 8 * r) / 8]);
        #endif
      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        uint8_t output_bits = 0;
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits ^= (z_bits_by_two[r*8 + bit_l] ^ 
                            (z_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (z_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (z_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (z_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ x[(i + 8 * r) / 8])
                                )
                              )
                            )
                          ) << bit_l;
        }
        bf_y[iy + r] = bf256_byte_combine_bits(output_bits);
        bf_y_sq[iy + r] = bf256_byte_combine_bits_sq(output_bits);
      }
      #else
      bf_y[iy + 0] = bf256_add(bf256_mul(bf_z_hat[0], bf_two), bf256_mul(bf_z_hat[1], bf_three));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_z_hat[2]);
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_z_hat[3]);
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_x_hat[0]);

      bf_y[iy + 1] = bf256_add(bf_z_hat[0], bf256_mul(bf_z_hat[1], bf_two));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_z_hat[2], bf_three));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf_z_hat[3]);
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf_x_hat[1]);

      bf_y[iy + 2] = bf256_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_z_hat[2], bf_two));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_z_hat[3], bf_three));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf_x_hat[2]);

      bf_y[iy + 3] = bf256_add(bf256_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_z_hat[2]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_z_hat[3], bf_two));
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_x_hat[3]);
      #endif

    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void em_enc_forward_256(const bf256_t* bf_z, const bf256_t* bf_x, 
                                bf256_t* bf_y, bf256_t* bf_y_sq)
#else
static void em_enc_forward_256(const bf256_t* bf_z, const bf256_t* bf_x, 
                                bf256_t* bf_y)
#endif
{

  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_256F_Nwd; j++) {

    bf_y[j] = bf256_byte_combine(bf_z + 8 * j);
    if (bf_x) {
      bf_y[j] = bf256_add(bf_y[j], bf256_byte_combine(bf_x + 8 * j));
    }

    #if defined(ALLOW_ZERO_SBOX)
    if(bf_x) {
      bf256_t tmp[8];
      for (unsigned int bit_j = 0; bit_j < 8; ++bit_j) {
        tmp[bit_j] = bf_z[j * 8 + bit_j] ^ bf_x[j * 8 + bit_j];
      }
      bf_y_sq[j] = bf256_byte_combine_sq(tmp);
    } 
    else {
      bf_y_sq[j] = bf256_byte_combine_sq(bf_z + 8 * j);
    }
    #endif
  }

  const bf256_t bf_two   = bf256_byte_combine_bits(2);
  const bf256_t bf_three = bf256_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_256F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {

      const unsigned int i  = 32 * FAEST_EM_256F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_256F_Nwd * j + 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_z_hat[4];
      #if defined(ALLOW_ZERO_SBOX)
      bf256_t bf_z_bits[4*8];
      bf256_t bf_z_bits_by_two[8];
      bf256_t output_bits[8];
      #endif
      
      for (unsigned int r = 0; r < 4; r++) {

        #if defined(ALLOW_ZERO_SBOX)
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          bf_z_bits[r*8 + bit_l] = bf_z[(i + 8 * r) + bit_l];
        }
        // mul by two
        bf_z_bits_by_two[0] = bf_z_bits[r*8 + 7];
        bf_z_bits_by_two[1] = bf_z_bits[r*8 + 7] ^ bf_z_bits[r*8 + 0];
        bf_z_bits_by_two[2] = bf_z_bits[r*8 + 1];
        bf_z_bits_by_two[3] = bf_z_bits[r*8 + 7] ^ bf_z_bits[r*8 + 2];
        bf_z_bits_by_two[4] = bf_z_bits[r*8 + 7] ^ bf_z_bits[r*8 + 3];
        bf_z_bits_by_two[5] = bf_z_bits[r*8 + 4];
        bf_z_bits_by_two[6] = bf_z_bits[r*8 + 5];
        bf_z_bits_by_two[7] = bf_z_bits[r*8 + 6];
        #else
        // Step: 12..13
        bf_z_hat[r] = bf256_byte_combine(bf_z + (i + 8 * r));
        if (bf_x) {
          bf_x_hat[r] = bf256_byte_combine(bf_x + (i + 8 * r));
        } else {
          bf_x_hat[r] = bf256_zero();
        }
        #endif
      }

      #if defined(ALLOW_ZERO_SBOX)
      for (unsigned int r = 0; r < 4; ++r) {
        bf256_t output_bits[8];
        // Similar to the else
        for (unsigned int bit_l = 0; bit_l < 8; ++bit_l) {
          output_bits[bit_l] = (bf_z_bits_by_two[r*8 + bit_l] ^ 
                            (bf_z_bits[((r + 2) % 4) * 8 + bit_l]  ^ 
                              (bf_z_bits[((r + 3) % 4) * 8 + bit_l] ^
                                (bf_z_bits[((r + 1) % 4) * 8 + bit_l] ^
                                  (bf_z_bits_by_two[((r + 1) % 4) * 8 + bit_l] ^ bf_x[(i + 8 * r) / 8])
                                )
                              )
                            )
                          );
        }
        bf_y[iy + r] = bf256_byte_combine(output_bits);
        bf_y_sq[iy + r] = bf256_byte_combine_sq(output_bits);
      }
      #else
      bf_y[iy + 0] = bf256_add(bf256_mul(bf_z_hat[0], bf_two), bf256_mul(bf_z_hat[1], bf_three));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_z_hat[2]);
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_z_hat[3]);
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf_x_hat[0]);

      bf_y[iy + 1] = bf256_add(bf_z_hat[0], bf256_mul(bf_z_hat[1], bf_two));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_z_hat[2], bf_three));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf_z_hat[3]);
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf_x_hat[1]);

      bf_y[iy + 2] = bf256_add(bf_z_hat[0], bf_z_hat[1]);
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_z_hat[2], bf_two));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_z_hat[3], bf_three));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf_x_hat[2]);

      bf_y[iy + 3] = bf256_add(bf256_mul(bf_z_hat[0], bf_three), bf_z_hat[1]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_z_hat[2]);
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_z_hat[3], bf_two));
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf_x_hat[3]);
      #endif
    }
  }
}

#if defined(ALLOW_ZERO_SBOX)
static void em_enc_backward_256_1(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf256_t* y_out, bf256_t* y_out_sq)
#else
static void em_enc_backward_256_1(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf256_t* y_out)
#endif
{
  // only called with Mtag == Mkey == 0

  for (unsigned int j = 0; j < FAEST_EM_256F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {

      for (unsigned int r = 0; r < 4; r++) {

        unsigned int icol = (c - r + FAEST_EM_256F_Nwd) % FAEST_EM_256F_Nwd;
        if (FAEST_EM_256F_Nwd == 8 && r >= 2) {
          icol = (icol - 1 + FAEST_EM_256F_Nwd) % FAEST_EM_256F_Nwd;
        }
        unsigned int ird = FAEST_EM_256F_LAMBDA + 32 * FAEST_EM_256F_Nwd * j + 32 * icol + 8 * r;
        uint8_t z_tilde  = 0;
        if (j < (FAEST_EM_256F_R - 1)) {
          z_tilde = z[ird / 8];
        } else {
          z_tilde = z_out[(ird - 32 * FAEST_EM_256F_Nwd * (j + 1)) / 8] ^ x[ird / 8];
        }

        // (bit spliced)
        // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
        const uint8_t y_tilde = rotr8(z_tilde, 7) ^ rotr8(z_tilde, 5) ^ rotr8(z_tilde, 2) ^ 0x5;

        // Step: 18
        y_out[4 * FAEST_EM_256F_Nwd * j + 4 * c + r] = bf256_byte_combine_bits(y_tilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[4 * FAEST_EM_256F_Nwd * j + 4 * c + r] = bf256_byte_combine_bits_sq(y_tilde);
        #endif
      }
    }
  }
}


#if defined(ALLOW_ZERO_SBOX)
static void em_enc_backward_256(const bf256_t* bf_z, const bf256_t* bf_x, const bf256_t* bf_z_out,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, 
                                bf256_t* y_out, bf256_t* y_out_sq)
#else
static void em_enc_backward_256(const bf256_t* bf_z, const bf256_t* bf_x, const bf256_t* bf_z_out,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, 
                                bf256_t* y_out)
#endif
{
  // Step: 1
  const bf256_t bf_delta = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t factor =
      bf256_mul_bit(bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int j = 0; j < FAEST_EM_256F_R; j++) {

    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {

      for (unsigned int r = 0; r <= 3; r++) {

        bf256_t bf_z_tilde[8];
        unsigned int icol = (c - r + FAEST_EM_256F_Nwd) % FAEST_EM_256F_Nwd;
        if (FAEST_EM_256F_Nwd == 8 && r >= 2) {
          icol = (icol - 1 + FAEST_EM_256F_Nwd) % FAEST_EM_256F_Nwd;
        }
        const unsigned int ird =
            FAEST_EM_256F_LAMBDA + 32 * FAEST_EM_256F_Nwd * j + 32 * icol + 8 * r;

        if (j < (FAEST_EM_256F_R - 1)) {
          memcpy(bf_z_tilde, bf_z + ird, sizeof(bf_z_tilde));
        } else {
          for (unsigned int i = 0; i < 8; ++i) {
            // Step: 12
            bf_z_tilde[i] = bf_z_out[ird - 32 * FAEST_EM_256F_Nwd * (j + 1) + i];
            if (bf_x) {
              bf_z_tilde[i] = bf256_add(bf_z_tilde[i], bf_x[ird + i]);
            }
          }
        }

        bf256_t bf_y_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf256_add(bf256_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                    bf_z_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);

        // Step: 18
        y_out[4 * FAEST_EM_256F_Nwd * j + 4 * c + r] = bf256_byte_combine(bf_y_tilde);
        #if defined(ALLOW_ZERO_SBOX)
        y_out_sq[4 * FAEST_EM_256F_Nwd * j + 4 * c + r] = bf256_byte_combine_sq(bf_y_tilde);
        #endif
      }
    }
  }
}

// Mkey = 0, this is for the prover
static void em_enc_constraints_Mkey_0_256(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                          const bf256_t* bf_v, zk_hash_256_ctx* a0_ctx,
                                          zk_hash_256_ctx* a1_ctx) {
  // Step 6
  uint8_t w_out[FAEST_EM_256F_LAMBDA / 8];
  xor_u8_array(out, w, w_out, sizeof(w_out));

  bf256_t bf_s[FAEST_EM_256F_Senc]; // inverse input
  bf256_t bf_vs[FAEST_EM_256F_Senc]; // inverse input tag
  #if defined(ALLOW_ZERO_SBOX)
  bf256_t bf_s_sq[FAEST_EM_256F_Senc]; // inverse input sq
  bf256_t bf_vs_sq[FAEST_EM_256F_Senc]; // input input tag sq
  em_enc_forward_256_1(w, x, bf_s, bf_s_sq);
  em_enc_forward_256(bf_v, NULL, bf_vs, bf_vs_sq);
  #else
  em_enc_forward_256_1(w, x, bf_s);
  em_enc_forward_256(bf_v, NULL, bf_vs);
  #endif

  bf256_t bf_s_dash[FAEST_EM_256F_Senc]; // inverse output
  bf256_t bf_vs_dash[FAEST_EM_256F_Senc]; // inverse output tag
  #if defined(ALLOW_ZERO_SBOX)
  bf256_t bf_s_dash_sq[FAEST_EM_256F_Senc]; // inverse output sq
  bf256_t bf_vs_dash_sq[FAEST_EM_256F_Senc]; // inverse output tag sq
  em_enc_backward_256_1(w, x, w_out, bf_s_dash, bf_s_dash_sq);
  em_enc_backward_256(bf_v, NULL, bf_v, 1, 0, NULL, bf_vs_dash, bf_vs_dash_sq);
  #else
  em_enc_backward_256_1(w, x, w_out, bf_s_dash);
  em_enc_backward_256(bf_v, NULL, bf_v, 1, 0, NULL, bf_vs_dash);
  #endif  

  for (unsigned int j = 0; j < FAEST_EM_256F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Psuedoinverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf256_t mul_tag_1 = bf256_mul(bf_vs_sq[j], bf_vs_dash[j]);   // mul1.mac0
    const bf256_t mul_tag_2 = bf256_mul(bf_vs[j], bf_vs_dash_sq[j]);   // mul2.mac0
    // add the tags with the values and multiply the results
    const bf256_t mul_val_1 = bf256_mul(                        // mul1.mac1
                        bf256_add(bf_s_sq[j], bf_vs_sq[j]), 
                        bf256_add(bf_s_dash[j], bf_vs_dash[j])
                        );
    const bf256_t mul_val_2 = bf256_mul(                        // mul2.mac1
                        bf256_add(bf_s[j], bf_vs[j]),
                        bf256_add(bf_s_dash_sq[j], bf_vs_dash_sq[j])
                        );

    // the constnat term
    zk_hash_256_update(a0_ctx, mul_tag_1);
    // the linear term
    zk_hash_256_update(a1_ctx, bf256_add(
                                        bf256_add(mul_tag_1, bf_s[j]),   // mul1.mac0 + x.val  (because here we check x^2*y = x)
                                        bf256_add(mul_val_1, bf_vs[j])   // mul1.mac1 + x.mac
                                        )
                      );
    // the constant term
    zk_hash_256_update(a0_ctx, mul_tag_2);
    // the linear term
    zk_hash_256_update(a1_ctx, bf256_add(
                                        bf256_add(mul_tag_2, bf_s_dash[j]),  // mul2.mac0 + y.val  (because here we check x*y^2 = y)
                                        bf256_add(mul_val_2, bf_vs_dash[j])  // mul2.mac1 + y.mac
                                        )
                      ); 
    #else
    // Inverse contraint
    // Quicksilver multiplication
    // multiply the tags
    const bf256_t mul_tag = bf256_mul(bf_vs[j], bf_vs_dash[j]); // mac0
    // add the tags with the values and multiply the results
    const bf256_t mul_val = bf256_mul(bf256_add(bf_s[j], bf_vs[j]), bf256_add(bf_s_dash[j], bf_vs_dash[j])); // mac1
    
    // this is the constant term
    zk_hash_256_update(a0_ctx, mul_tag);
    // this is the linear term
    zk_hash_256_update(a1_ctx, bf256_add(
                                        bf256_add(mul_tag, bf256_one),    // mac0 + 1     (because here we check x*y = 1)
                                        mul_val                             // (mac0 + 1) + mac1
                                        )
                      );
    #endif
  }
}

// Mkey = 1, this is for the verifier
static void em_enc_constraints_Mkey_1_256(const uint8_t* out, const uint8_t* x, const bf256_t* bf_q,
                                          const uint8_t* delta, zk_hash_256_ctx* b0_ctx) {
  // Step: 18, 19
  const bf256_t bf_delta = bf256_load(delta);
  bf256_t* bf_x = faest_aligned_alloc(BF256_ALIGN, sizeof(bf256_t) * 256 * (FAEST_EM_256F_R + 1));
  for (unsigned int i = 0; i < 256 * (FAEST_EM_256F_R + 1); i++) {
    bf_x[i] = bf256_mul_bit(bf_delta, ptr_get_bit(x, i));
  }

  // Step 21
  bf256_t* bf_q_out = faest_aligned_alloc(BF256_ALIGN, sizeof(bf256_t) * FAEST_EM_256F_LAMBDA);
  for (unsigned int i = 0; i < FAEST_EM_256F_LAMBDA; i++) {
    bf_q_out[i] = bf256_add(bf256_mul_bit(bf_delta, ptr_get_bit(out, i)), bf_q[i]);
  }

  bf256_t bf_qs[FAEST_EM_256F_Senc];
  bf256_t bf_qs_dash[FAEST_EM_256F_Senc];
  #if defined(ALLOW_ZERO_SBOX)
  bf256_t bf_qs_sq[FAEST_EM_256F_Senc];
  bf256_t bf_qs_dash_sq[FAEST_EM_256F_Senc];
  em_enc_forward_256(bf_q, bf_x, bf_qs, bf_qs_sq);
  em_enc_backward_256(bf_q, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash, bf_qs_dash_sq);
  #else
  em_enc_forward_256(bf_q, bf_x, bf_qs);
  em_enc_backward_256(bf_q, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash);
  #endif
  faest_aligned_free(bf_q_out);
  faest_aligned_free(bf_x);

  // Step: 13..14
  bf256_t delta_squared = bf256_mul(bf_delta, bf_delta);

  for (unsigned int j = 0; j < FAEST_EM_256F_Senc; j++) {

    #if defined(ALLOW_ZERO_SBOX)
    // Quicksilver multiplication
    // multiply the tags
    const bf256_t mul_tag_1 = bf256_mul(bf_qs_sq[j], bf_qs_dash[j]);   // mul1.mac0
    const bf256_t mul_tag_2 = bf256_mul(bf_qs[j], bf_qs_dash_sq[j]);   // mul2.mac0

    // the constnat term
    zk_hash_256_update(b0_ctx, bf256_add(
                                        mul_tag_1,
                                        bf256_mul(delta_squared, bf_qs[j])
                                        )
    );
    // the constnat term
    zk_hash_256_update(b0_ctx, bf256_add(
                                        mul_tag_2,
                                        bf256_mul(delta_squared, bf_qs_dash[j])
                                        )
    );
    #else
    // the constant term
    zk_hash_256_update(b0_ctx, bf256_add(
                                        bf256_mul(bf_qs[j], bf_qs_dash[j]), 
                                        delta_squared));
    #endif
  }
}

/* 
static void em_prove_256(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* owf_in,
                          const uint8_t* owf_out, const uint8_t* chall_2, uint8_t* a0_tilde,
                          uint8_t* a12_tilde, const faest_paramset_t* params) {
  // copy expanded key in to an array
  uint8_t x[FAEST_EM_256F_LAMBDA * (FAEST_EM_256F_R + 1) / 8];
  {
    aes_round_keys_t round_keys;
    rijndael256_init_round_keys(&round_keys, owf_in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_256F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_256F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  bf256_t* bf_v = column_to_row_major_and_shrink_V_256(V, FAEST_EM_256F_Lenc);
  zk_hash_256_ctx a0_ctx;
  zk_hash_256_ctx a1_ctx;

  zk_hash_256_init(&a0_ctx, chall_2);
  zk_hash_256_init(&a1_ctx, chall_2);
  em_enc_constraints_Mkey_0_256(owf_out, x, w, bf_v, &a0_ctx, &a1_ctx);

  zk_hash_256_finalize(a_tilde, &a1_ctx, bf256_load(u + FAEST_EM_256F_Lenc / 8));
  zk_hash_256_finalize(b_tilde, &a0_ctx, bf256_sum_poly(bf_v + FAEST_EM_256F_Lenc));

  faest_aligned_free(bf_v);
}

static uint8_t* em_verify_256(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                              const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                              const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.tau0;
  const unsigned int k0  = params->faest_param.k;
  const unsigned int t1  = params->faest_param.tau1;
  const unsigned int k1  = (t0 != 0) ? k0 - 1 : k0;

  const uint8_t* delta = chall_3;

  for (unsigned int i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t decoded_challenge[MAX_DEPTH];
    ChalDec(chall_3, i, k0, t0, k1, t1, decoded_challenge);
    for (unsigned int j = 0; j < depth; j++, ++col) {
      if (decoded_challenge[j] == 1) {
        xor_u8_array(d, Q[col], Q[col], (FAEST_EM_256F_Lenc + 7) / 8);
      }
    }
  }

  bf256_t* bf_q = column_to_row_major_and_shrink_V_256(Q, FAEST_EM_256F_Lenc);

  // copy expanded key in to an array
  uint8_t x[FAEST_EM_256F_LAMBDA * (FAEST_EM_256F_R + 1) / 8];
  {
    aes_round_keys_t round_keys;
    rijndael256_init_round_keys(&round_keys, in);
    uint8_t* tmp_x = x;
    for (unsigned int r = 0; r != FAEST_EM_256F_R + 1; ++r) {
      for (unsigned int i = 0; i != FAEST_EM_256F_Nwd; ++i) {
        memcpy(tmp_x, round_keys.round_keys[r][i], sizeof(aes_word_t));
        tmp_x += sizeof(aes_word_t);
      }
    }
  }

  zk_hash_256_ctx b0_ctx;
  zk_hash_256_init(&b0_ctx, chall_2);
  em_enc_constraints_Mkey_1_256(out, x, bf_q, delta, &b0_ctx);

  uint8_t* q_tilde = malloc(FAEST_EM_256F_LAMBDA / 8);
  zk_hash_256_finalize(q_tilde, &b0_ctx, bf256_sum_poly(bf_q + FAEST_EM_256F_Lenc));
  faest_aligned_free(bf_q);

  bf256_t bf_qtilde = bf256_load(q_tilde);
  bf256_store(q_tilde, bf256_add(bf_qtilde, bf256_mul(bf256_load(a_tilde), bf256_load(delta))));

  return q_tilde;
}

 */

// dispatchers

void aes_prove(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde, const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* owf_in,
               const uint8_t* owf_out, const uint8_t* chall_2, const faest_paramset_t* params) {
  switch (params->faest_param.lambda) {
  case 256:
    if (params->faest_param.Lke) {
      // aes_prove_256(w, u, V, owf_in, owf_out, chall_2, a0_tilde, a12_tilde, params);
    } else {
      // em_prove_256(w, u, V, owf_in, owf_out, chall_2, a0_tilde, a12_tilde);
    }
    break;
  case 192:
    if (params->faest_param.Lke) {
      // aes_prove_192(w, u, V, owf_in, owf_out, chall_2, a0_tilde, a12_tilde, params);
    } else {
      // em_prove_192(w, u, V, owf_in, owf_out, chall_2, a0_tilde, a12_tilde);
    }
    break;
  default:
    if (params->faest_param.Lke) {
      aes_prove_128(a0_tilde, a1_tilde, a2_tilde, w, u, V, owf_in, owf_out, chall_2, params);
    } else {
      // em_prove_128(w, u, V, owf_in, owf_out, chall_2, a0_tilde, a12_tilde);
    }
  }
}

uint8_t* aes_verify(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2, const uint8_t* chall_3,
                    const uint8_t* a1_tilde, const uint8_t* a2_tilde, const uint8_t* owf_in, const uint8_t* owf_out,
                    const faest_paramset_t* params) {
  switch (params->faest_param.lambda) {
  case 256:
    if (params->faest_param.Lke) {
      // return aes_verify_256(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    } else {
      // return em_verify_256(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    }
  case 192:
    if (params->faest_param.Lke) {
      // return aes_verify_192(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    } else {
      // return em_verify_192(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    }
  default:
    if (params->faest_param.Lke) {
      return aes_verify_128(d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
    } else {
      // return em_verify_128(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    }
  }
}
