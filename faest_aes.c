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

static_assert(FAEST_128F_L == FAEST_128S_L, "Invalid parameters");
static_assert(FAEST_128F_LAMBDA == FAEST_128S_LAMBDA, "Invalid parameters");
static_assert(FAEST_128F_Lke == FAEST_128S_Lke, "Invalid parameters");
static_assert(FAEST_128F_Nwd == FAEST_128S_Nwd, "Invalid parameters");
static_assert(FAEST_128F_R == FAEST_128S_R, "Invalid parameters");
static_assert(FAEST_128F_Senc == FAEST_128S_Senc, "Invalid parameters");
static_assert(FAEST_128F_Ske == FAEST_128S_Ske, "Invalid parameters");

static_assert(FAEST_192F_L == FAEST_192S_L, "Invalid parameters");
static_assert(FAEST_192F_LAMBDA == FAEST_192S_LAMBDA, "Invalid parameters");
static_assert(FAEST_192F_Lke == FAEST_192S_Lke, "Invalid parameters");
static_assert(FAEST_192F_Nwd == FAEST_192S_Nwd, "Invalid parameters");
static_assert(FAEST_192F_R == FAEST_192S_R, "Invalid parameters");
static_assert(FAEST_192F_Senc == FAEST_192S_Senc, "Invalid parameters");
static_assert(FAEST_192F_Ske == FAEST_192S_Ske, "Invalid parameters");

static_assert(FAEST_256F_L == FAEST_256S_L, "Invalid parameters");
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

static bf128_t* column_to_row_major_and_shrink_V_128(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf128_t* new_v = malloc((ell + sizeof(bf128_t) * 8) * sizeof(bf128_t));
  for (unsigned int row = 0; row != ell + sizeof(bf128_t) * 8; ++row) {
    uint8_t new_row[sizeof(bf128_t)] = {0};
    for (unsigned int column = 0; column != sizeof(bf128_t) * 8; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf128_load(new_row);
  }

  return new_v;
}

static bf192_t* column_to_row_major_and_shrink_V_192(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf192_t* new_v = malloc((ell + sizeof(bf192_t) * 8) * sizeof(bf192_t));
  for (unsigned int row = 0; row != ell + sizeof(bf192_t) * 8; ++row) {
    uint8_t new_row[sizeof(bf192_t)] = {0};
    for (unsigned int column = 0; column != sizeof(bf192_t) * 8; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf192_load(new_row);
  }

  return new_v;
}

static bf256_t* column_to_row_major_and_shrink_V_256(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf256_t* new_v = malloc((ell + sizeof(bf256_t) * 8) * sizeof(bf256_t));
  for (unsigned int row = 0; row != ell + sizeof(bf256_t) * 8; ++row) {
    uint8_t new_row[sizeof(bf256_t)] = {0};
    for (unsigned int column = 0; column != sizeof(bf256_t) * 8; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf256_load(new_row);
  }

  return new_v;
}

// m == 1 implementations

static void aes_key_schedule_forward_1(const uint8_t* x, uint8_t* out,
                                       const faest_paramset_t* params) {
  // Step: 1 skipped (sanity check)

  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int R           = params->faest_param.R;
  const unsigned int Nwd         = params->faest_param.Nwd;
  const unsigned int lambdaBytes = lambda / 8;

  const unsigned int out_len = (R + 1) * 128 / 8;
  // Step 3
  memcpy(out, x, lambdaBytes);
  memset(out + lambdaBytes, 0, out_len - lambdaBytes);

  // Step: 4
  unsigned int i_wd = lambda;
  // Step: 5..10
  for (unsigned int j = Nwd; j < 4 * (R + 1); j++) {
    if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {
      memcpy(out + 32 * j / 8, x + i_wd / 8, 4);
      i_wd += 32;
    } else {
      for (unsigned int i = 0; i < 32; i += 8) {
        // bit spliced
        out[(32 * j + i) / 8] |= out[(32 * (j - Nwd) + i) / 8] ^ out[(32 * (j - 1) + i) / 8];
      }
    }
  }
}

static void aes_key_schedule_backward_1(const uint8_t* x, const uint8_t* xk, uint8_t* out,
                                        const faest_paramset_t* params) {
  // Step: 1 skipped (sanity check)

  const unsigned int lambda = params->faest_param.lambda;
  const unsigned int Ske    = params->faest_param.Ske;

  // Step: 2
  unsigned int iwd   = 0;
  unsigned int c     = 0;
  bool rmvRcon       = true;
  unsigned int ircon = 0;

  for (unsigned int j = 0; j < Ske; j++) {
    // Step 7 (bit sliced)
    uint8_t x_tilde = x[j] ^ xk[iwd + c];

    // Step 8
    // this function is only called with Mtag == Mkey == 0
    if (/* Mtag == 0 && */ rmvRcon == true && c == 0) {
      // Steps 12 and 13, bitsliced; delta is always 0
      x_tilde ^= Rcon[ircon];
      ++ircon;
    }

    // Step: 15..19 (bit spliced)
    const uint8_t y_tilde = rotr8(x_tilde, 7) ^ rotr8(x_tilde, 5) ^ rotr8(x_tilde, 2);
    // this function is only called with Mtag == Mkey == 0
    // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
    out[j] = y_tilde ^ 0x5;

    // Step: 20
    ++c;
    if (c == 4) {
      c = 0;
      if (lambda == 192) {
        iwd += 192 / 8;
      } else {
        iwd += 128 / 8;
        if (lambda == 256) {
          rmvRcon = !rmvRcon;
        }
      }
    }
  }
}

// lambda == 128 implementation

static void aes_key_schedule_forward_128(const bf128_t* v, bf128_t* bf_out) {
  // Step: 1 sanity check (skipped)

  memcpy(bf_out, v, FAEST_128F_LAMBDA * sizeof(bf128_t));

  // Step: 4
  unsigned int i_wd = FAEST_128F_LAMBDA;
  // Step: 5..10
  for (unsigned int j = FAEST_128F_Nwd; j < 4 * (FAEST_128F_R + 1); j++) {
    if ((j % FAEST_128F_Nwd) == 0 || (FAEST_128F_Nwd > 6 && (j % FAEST_128F_Nwd) == 4)) {
      // copy all at once
      memcpy(bf_out + j * 32, v + i_wd, sizeof(bf128_t) * 32);
      i_wd += 32;
    } else {
      for (unsigned int i = 0; i < 32; i++) {
        bf_out[(32 * j) + i] =
            bf128_add(bf_out[32 * (j - FAEST_128F_Nwd) + i], bf_out[32 * (j - 1) + i]);
      }
    }
  }
}

static void aes_key_schedule_backward_128(const bf128_t* v, const bf128_t* Vk, uint8_t Mtag,
                                          uint8_t Mkey, const uint8_t* delta, bf128_t* bf_out) {
  // Step: 1
  assert(!((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)));

  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();

  unsigned int iwd   = 0;
  unsigned int c     = 0;
  unsigned int ircon = 0;

  bf128_t bf_minus_mkey       = bf128_from_bit(1 ^ Mkey);
  uint8_t minus_mtag          = 1 ^ Mtag;
  bf128_t bf_mkey_times_delta = bf128_mul_bit(bf_delta, Mkey);
  bf_mkey_times_delta         = bf128_add(bf_mkey_times_delta, bf_minus_mkey);

  for (unsigned int j = 0; j < FAEST_128F_Ske; j++) {
    // Step 7
    bf128_t bf_x_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      bf_x_tilde[i] = bf128_add(v[8 * j + i], Vk[iwd + 8 * c + i]);
    }

    if (Mtag == 0 && c == 0) {
      // Step 9
      uint8_t r = Rcon[ircon];
      ircon     = ircon + 1;

      bf128_t bf_r[8];
      for (unsigned int i = 0; i < 8; i++) {
        // Step 12
        bf_r[i] = bf128_mul_bit(bf_mkey_times_delta, get_bit(r, i));
        // Step 13
        bf_x_tilde[i] = bf128_add(bf_x_tilde[i], bf_r[i]);
      }
    }

    for (unsigned int i = 0; i < 8; ++i) {
      bf_out[i + 8 * j] = bf128_add(bf128_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                    bf_x_tilde[(i + 2) % 8]);
    }
    bf_out[0 + 8 * j] =
        bf128_add(bf_out[0 + 8 * j], bf128_mul_bit(bf_mkey_times_delta, minus_mtag));
    bf_out[2 + 8 * j] =
        bf128_add(bf_out[2 + 8 * j], bf128_mul_bit(bf_mkey_times_delta, minus_mtag));
    c = c + 1;

    if (c == 4) {
      c = 0;
      iwd += 128;
    }
  }
}

static void aes_key_schedule_constraints_128(const uint8_t* w, const bf128_t* v, const uint8_t Mkey,
                                             const bf128_t* q, const uint8_t* delta, bf128_t* A0,
                                             bf128_t* A1, uint8_t* k, bf128_t* vk, bf128_t* B,
                                             bf128_t* qk, const faest_paramset_t* params) {
  if (Mkey == 0) {
    // Step: 2
    aes_key_schedule_forward_1(w, k, params);

    // Step: 3
    aes_key_schedule_forward_128(v, vk);

    // Step: 4
    uint8_t w_dash[FAEST_128F_Ske];
    aes_key_schedule_backward_1(w + FAEST_128F_LAMBDA / 8, k, w_dash, params);

    // Step: 5
    bf128_t v_w_dash[FAEST_128F_Ske * 8];
    aes_key_schedule_backward_128(v + FAEST_128F_LAMBDA, vk, 1, 0, NULL, v_w_dash);

    // Step: 6..8
    unsigned int iwd = 32 * (FAEST_128F_Nwd - 1);
    for (unsigned int j = 0; j < FAEST_128F_Ske / 4; j++) {
      bf128_t bf_k_hat[4];
      bf128_t bf_v_k_hat[4];
      bf128_t bf_w_dash_hat[4];
      bf128_t bf_v_w_dash_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 10..11
        bf_k_hat[(r + 3) % 4]   = bf128_byte_combine_bits(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat[(r + 3) % 4] = bf128_byte_combine(vk + (iwd + 8 * r));
        bf_w_dash_hat[r]        = bf128_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
        bf_v_w_dash_hat[r]      = bf128_byte_combine(v_w_dash + (32 * j + 8 * r));
      }
      // Step: 13..17
      for (unsigned int r = 0; r <= 3; r++) {
        A0[4 * j + r] = bf128_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
        A1[4 * j + r] =
            bf128_add(bf128_add(bf128_mul(bf128_add(bf_k_hat[r], bf_v_k_hat[r]),
                                          bf128_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                bf128_from_bf8(bf8_one())),
                      A0[4 * j + r]);
      }
      if (FAEST_128F_LAMBDA == 192) {
        iwd = iwd + 192;
      } else {
        iwd = iwd + 128;
      }
    }
    return;
  }

  // Step: 19..20
  aes_key_schedule_forward_128(q, qk);
  bf128_t q_w_dash[FAEST_128F_Ske * 8];
  aes_key_schedule_backward_128(&q[FAEST_128F_LAMBDA], qk, 0, 1, delta, q_w_dash);

  const bf128_t bf_delta = bf128_load(delta);

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_128F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_128F_Ske / 4; j++) {
    bf128_t bf_q_hat_k[4];
    bf128_t bf_q_hat_w_dash[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 25..26
      bf_q_hat_k[(r + 3) % 4] = bf128_byte_combine(qk + ((iwd + 8 * r)));
      bf_q_hat_w_dash[r]      = bf128_byte_combine(q_w_dash + ((32 * j + 8 * r)));
    }
    // Step: 27
    for (unsigned int r = 0; r <= 3; r++) {
      bf128_t bf_tmp = bf128_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      B[4 * j + r]   = bf128_add(bf_tmp, bf128_mul(bf_delta, bf_delta));
    }
    if (FAEST_128F_LAMBDA == 192) {
      iwd = iwd + 192;
    } else {
      iwd = iwd + 128;
    }
  }
}

static void aes_enc_forward_128_1(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  bf128_t* bf_y) {
  // called only with Mtag == Mkey == 0

  // Step: 2
  for (unsigned int i = 0; i < 16; i++) {
    // Step: 3, 4 (bit spliced)
    // -((1 ^ Mtag) & (1 ^ Mkey)) == 0xFF
    const uint8_t xin = in[i];
    // Step: 5
    bf_y[i] = bf128_add(bf128_byte_combine_bits(xin), bf128_byte_combine_bits(xk[i]));
  }

  const bf128_t bf_two   = bf128_byte_combine_bits(2);
  const bf128_t bf_three = bf128_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_128F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_xk_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf128_byte_combine_bits(x[(ix + 8 * r) / 8]);
        bf_xk_hat[r] = bf128_byte_combine_bits(xk[(ik + 8 * r) / 8]);
      }

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
    }
  }
}

static void aes_enc_forward_128(const bf128_t* bf_x, const bf128_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf128_t* bf_y) {
  const bf128_t bf_delta      = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t bf_minus_mtag = bf128_from_bit(1 ^ Mtag);
  const bf128_t bf_minus_mkey = bf128_from_bit(1 ^ Mkey);
  const bf128_t bf_mkey       = bf128_from_bit(Mkey);

  // Step: 2..4
  for (unsigned int i = 0; i < 16; i++) {
    bf128_t bf_xin[8];
    for (unsigned int j = 0; j < 8; j++) {
      bf_xin[j] = bf128_mul(bf128_mul_bit(bf_minus_mtag, get_bit(in[i], j)),
                            bf128_add(bf128_mul(bf_mkey, bf_delta), bf_minus_mkey));
    }
    // Step: 5
    bf_y[i] = bf128_add(bf128_byte_combine(bf_xin), bf128_byte_combine(bf_xk + (8 * i)));
  }

  const bf128_t bf_two   = bf128_byte_combine_bits(2);
  const bf128_t bf_three = bf128_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_128F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_xk_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf128_byte_combine(bf_x + (ix + 8 * r));
        bf_xk_hat[r] = bf128_byte_combine(bf_xk + (ik + 8 * r));
      }

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
    }
  }
}

static void aes_enc_backward_128_1(const uint8_t* x, const uint8_t* xk, const uint8_t* out,
                                   bf128_t* y_out) {
  // called only with Mtag == Mkey == 0

  uint8_t xtilde;
  // Step:2..4
  for (unsigned int j = 0; j < FAEST_128F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
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
      }
    }
  }
}

static void aes_enc_backward_128(const bf128_t* bf_x, const bf128_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf128_t* y_out) {
  // Step: 1
  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t factor =
      bf128_mul_bit(bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  // Step: 2..4
  for (unsigned int j = 0; j < FAEST_128F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
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
      }
    }
  }
}

static void aes_enc_constraints_128(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                    const bf128_t* v, const uint8_t* k, const bf128_t* vk,
                                    uint8_t Mkey, const bf128_t* q, const bf128_t* qk,
                                    const uint8_t* delta, bf128_t* A0, bf128_t* A1, bf128_t* B) {
  if (Mkey == 0) {
    bf128_t s[FAEST_128F_Senc];
    bf128_t vs[FAEST_128F_Senc];
    bf128_t s_dash[FAEST_128F_Senc];
    bf128_t vs_dash[FAEST_128F_Senc];
    aes_enc_forward_128_1(w, k, in, s);
    aes_enc_forward_128(v, vk, in, 1, 0, NULL, vs);
    aes_enc_backward_128_1(w, k, out, s_dash);
    aes_enc_backward_128(v, vk, 1, 0, NULL, out, vs_dash);

    for (unsigned int j = 0; j < FAEST_128F_Senc; j++) {
      A0[j] = bf128_mul(vs[j], vs_dash[j]);
      A1[j] = bf128_add(
          bf128_add(bf128_mul(bf128_add(s[j], vs[j]), bf128_add(s_dash[j], vs_dash[j])), A0[j]),
          bf128_one());
    }
  } else {
    // Step: 11..12
    bf128_t qs[FAEST_128F_Senc];
    bf128_t qs_dash[FAEST_128F_Senc];
    aes_enc_forward_128(q, qk, in, 0, 1, delta, qs);
    aes_enc_backward_128(q, qk, 0, 1, delta, out, qs_dash);

    // Step: 13..14
    bf128_t minus_part = bf128_mul(bf128_load(delta), bf128_load(delta));
    for (unsigned int j = 0; j < FAEST_128F_Senc; j++) {
      B[j] = bf128_add(bf128_mul(qs[j], qs_dash[j]), minus_part);
    }
  }
}

static void aes_prove_128(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
                          const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde,
                          uint8_t* b_tilde, const faest_paramset_t* params) {
  // Step: 1..2
  bf128_t* bf_v = column_to_row_major_and_shrink_V_128(V, FAEST_128F_L);

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7
  bf128_t* A0 = malloc(sizeof(bf128_t) * (FAEST_128F_Ske + FAEST_128F_Senc + 1));
  bf128_t* A1 = malloc(sizeof(bf128_t) * (FAEST_128F_Ske + FAEST_128F_Senc + 1));
  uint8_t* k  = malloc((FAEST_128F_R + 1) * 128 / 8);
  bf128_t* vk = malloc(sizeof(bf128_t) * ((FAEST_128F_R + 1) * 128));
  bf128_t* qk = malloc(sizeof(bf128_t) * ((FAEST_128F_R + 1) * 128));
  aes_key_schedule_constraints_128(w, bf_v, 0, NULL, NULL, A0, A1, k, vk, NULL, qk, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  aes_enc_constraints_128(in, out, w + FAEST_128F_Lke / 8, bf_v + FAEST_128F_Lke, k, vk, 0, NULL,
                          NULL, NULL, A0 + FAEST_128F_Ske, A1 + FAEST_128F_Ske, NULL);
  // Step: 12 (beta == 1)
  free(qk);
  free(vk);
  free(k);

  // Step: 16..18
  A1[FAEST_128F_Ske + FAEST_128F_Senc] = bf128_load(u + FAEST_128F_L / 8);
  A0[FAEST_128F_Ske + FAEST_128F_Senc] = bf128_sum_poly(bf_v + FAEST_128F_L);
  free(bf_v);

  zk_hash_128(a_tilde, chall, A1, FAEST_128F_Ske + FAEST_128F_Senc);
  zk_hash_128(b_tilde, chall, A0, FAEST_128F_Ske + FAEST_128F_Senc);

  free(A0);
  free(A1);
}

static uint8_t* aes_verify_128(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                               const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                               const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.t0;
  const unsigned int k0  = params->faest_param.k0;
  const unsigned int t1  = params->faest_param.t1;
  const unsigned int k1  = params->faest_param.k1;

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
        xor_u8_array(d, Q[col], Q[col], (FAEST_128F_L + 7) / 8);
      }
    }
  }

  // Step: 11..12
  bf128_t* bf_q = column_to_row_major_and_shrink_V_128(Q, FAEST_128F_L);

  // Step: 13
  bf128_t* qk  = malloc(sizeof(bf128_t) * ((FAEST_128F_R + 1) * 128));
  bf128_t* B_0 = malloc(sizeof(bf128_t) * (FAEST_128F_Ske + FAEST_128F_Senc + 1));
  aes_key_schedule_constraints_128(NULL, NULL, 1, bf_q, delta, NULL, NULL, NULL, NULL, B_0, qk,
                                   params);

  // Step: 14
  bf128_t* B_1 = B_0 + FAEST_128F_Ske;
  aes_enc_constraints_128(in, out, NULL, NULL, NULL, NULL, 1, bf_q + FAEST_128F_Lke, qk, delta,
                          NULL, NULL, B_1);
  // Step: 18 (beta == 1)
  free(qk);

  // Step: 20
  B_0[FAEST_128F_Ske + FAEST_128F_Senc] = bf128_sum_poly(bf_q + FAEST_128F_L);
  free(bf_q);

  // Step 21
  uint8_t* q_tilde = malloc(FAEST_128F_LAMBDA / 8);
  zk_hash_128(q_tilde, chall_2, B_0, FAEST_128F_Ske + FAEST_128F_Senc);
  free(B_0);

  bf128_t bf_qtilde = bf128_load(q_tilde);
  bf128_store(q_tilde, bf128_add(bf_qtilde, bf128_mul(bf128_load(a_tilde), bf128_load(delta))));

  return q_tilde;
}

// lambda == 192 implementation

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

static void aes_key_schedule_constraints_192(const uint8_t* w, const bf192_t* v, const uint8_t Mkey,
                                             const bf192_t* q, const uint8_t* delta, bf192_t* A0,
                                             bf192_t* A1, uint8_t* k, bf192_t* vk, bf192_t* B,
                                             bf192_t* qk, const faest_paramset_t* params) {
  if (Mkey == 0) {
    // Step: 2
    aes_key_schedule_forward_1(w, k, params);

    // Step: 3
    aes_key_schedule_forward_192(v, vk);

    // Step: 4
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
      bf192_t bf_w_dash_hat[4];
      bf192_t bf_v_w_dash_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 10..11
        bf_k_hat[(r + 3) % 4]   = bf192_byte_combine_bits(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat[(r + 3) % 4] = bf192_byte_combine(vk + (iwd + 8 * r));
        bf_w_dash_hat[r]        = bf192_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
        bf_v_w_dash_hat[r]      = bf192_byte_combine(v_w_dash + (32 * j + 8 * r));
      }
      // Step: 13..17
      for (unsigned int r = 0; r <= 3; r++) {
        A0[4 * j + r] = bf192_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
        A1[4 * j + r] =
            bf192_add(bf192_add(bf192_mul(bf192_add(bf_k_hat[r], bf_v_k_hat[r]),
                                          bf192_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                bf192_from_bf8(bf8_one())),
                      A0[4 * j + r]);
      }
      if (FAEST_192F_LAMBDA == 192) {
        iwd = iwd + 192;
      } else {
        iwd = iwd + 128;
      }
    }
    return;
  }

  // Step: 19..20
  aes_key_schedule_forward_192(q, qk);
  bf192_t q_w_dash[FAEST_192F_Ske * 8];
  aes_key_schedule_backward_192(&q[FAEST_192F_LAMBDA], qk, 0, 1, delta, q_w_dash);

  const bf192_t bf_delta = bf192_load(delta);

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_192F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_192F_Ske / 4; j++) {
    bf192_t bf_q_hat_k[4];
    bf192_t bf_q_hat_w_dash[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 25..26
      bf_q_hat_k[(r + 3) % 4] = bf192_byte_combine(qk + ((iwd + 8 * r)));
      bf_q_hat_w_dash[r]      = bf192_byte_combine(q_w_dash + ((32 * j + 8 * r)));
    }
    // Step: 27
    for (unsigned int r = 0; r <= 3; r++) {
      bf192_t bf_tmp = bf192_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      B[4 * j + r]   = bf192_add(bf_tmp, bf192_mul(bf_delta, bf_delta));
    }
    if (FAEST_192F_LAMBDA == 192) {
      iwd = iwd + 192;
    } else {
      iwd = iwd + 128;
    }
  }
}

static void aes_enc_forward_192_1(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  uint8_t Mtag, uint8_t Mkey, bf192_t* bf_y) {
  // Step: 2
  for (unsigned int i = 0; i < 16; i++) {
    // Step: 3,4 (bit spliced)
    const uint8_t xin = in[i] & -((1 ^ Mtag) & (1 ^ Mkey));
    // Step: 5
    bf_y[i] = bf192_add(bf192_byte_combine_bits(xin), bf192_byte_combine_bits(xk[i]));
  }

  const bf192_t bf_two   = bf192_byte_combine_bits(2);
  const bf192_t bf_three = bf192_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_192F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_xk_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf192_byte_combine_bits(x[(ix + 8 * r) / 8]);
        bf_xk_hat[r] = bf192_byte_combine_bits(xk[(ik + 8 * r) / 8]);
      }

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
    }
  }
  return;
}

static void aes_enc_forward_192(const bf192_t* bf_x, const bf192_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf192_t* bf_y) {
  const bf192_t bf_delta      = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t bf_minus_mtag = bf192_from_bit(1 ^ Mtag);
  const bf192_t bf_minus_mkey = bf192_from_bit(1 ^ Mkey);

  // Step: 2..4
  for (unsigned int i = 0; i < 16; i++) {
    bf192_t bf_xin[8];
    for (unsigned int j = 0; j < 8; j++) {
      bf_xin[j] = bf192_mul(bf192_mul_bit(bf_minus_mtag, get_bit(in[i], j)),
                            bf192_add(bf192_mul_bit(bf_delta, Mkey), bf_minus_mkey));
    }
    // Step: 5
    bf_y[i] = bf192_add(bf192_byte_combine(bf_xin), bf192_byte_combine(bf_xk + (8 * i)));
  }

  const bf192_t bf_two   = bf192_byte_combine_bits(2);
  const bf192_t bf_three = bf192_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_192F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_xk_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf192_byte_combine(bf_x + (ix + 8 * r));
        bf_xk_hat[r] = bf192_byte_combine(bf_xk + (ik + 8 * r));
      }

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
    }
  }
}

static void aes_enc_backward_192_1(const uint8_t* x, const uint8_t* xk, uint8_t Mtag, uint8_t Mkey,
                                   const uint8_t* out, bf192_t* y_out) {
  uint8_t xtilde;
  // Step:2..4
  for (unsigned int j = 0; j < FAEST_192F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 5..6
        unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
        if (j < (FAEST_192F_R - 1)) {
          // Step: 7
          xtilde = x[ird / 8];
        } else {
          // Step: 9..11 (bit spliced)
          uint8_t xout = out[(ird - 128 * (FAEST_192F_R - 1)) / 8] & -((1 ^ Mtag) & (1 ^ Mkey));
          xtilde       = xout ^ xk[(128 + ird) / 8];
        }

        // Step: 12..17 (bit spliced)
        uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^ rotr8(xtilde, 2);
        ytilde ^= set_bit((1 ^ Mtag) & (1 ^ Mkey), 0);
        ytilde ^= set_bit((1 ^ Mtag) & (1 ^ Mkey), 2);

        // Step: 18
        y_out[16 * j + 4 * c + r] = bf192_byte_combine_bits(ytilde);
      }
    }
  }
  return;
}

static void aes_enc_backward_192(const bf192_t* bf_x, const bf192_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf192_t* y_out) {
  // Step: 1
  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t factor =
      bf192_mul_bit(bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  // Step: 2..4
  for (unsigned int j = 0; j < FAEST_192F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
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
      }
    }
  }
}

static void aes_enc_constraints_192(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                    const bf192_t* v, const uint8_t* k, const bf192_t* vk,
                                    uint8_t Mkey, const bf192_t* q, const bf192_t* qk,
                                    const uint8_t* delta, bf192_t* A0, bf192_t* A1, bf192_t* B) {
  if (Mkey == 0) {
    bf192_t s[FAEST_192F_Senc];
    bf192_t vs[FAEST_192F_Senc];
    bf192_t s_dash[FAEST_192F_Senc];
    bf192_t vs_dash[FAEST_192F_Senc];
    aes_enc_forward_192_1(w, k, in, 0, 0, s);
    aes_enc_forward_192(v, vk, in, 1, 0, NULL, vs);
    aes_enc_backward_192_1(w, k, 0, 0, out, s_dash);
    aes_enc_backward_192(v, vk, 1, 0, NULL, out, vs_dash);

    for (unsigned int j = 0; j < FAEST_192F_Senc; j++) {
      A0[j] = bf192_mul(vs[j], vs_dash[j]);
      A1[j] = bf192_add(
          bf192_add(bf192_mul(bf192_add(s[j], vs[j]), bf192_add(s_dash[j], vs_dash[j])), A0[j]),
          bf192_one());
    }
  } else {
    // Step: 11..12
    bf192_t qs[FAEST_192F_Senc];
    bf192_t qs_dash[FAEST_192F_Senc];
    aes_enc_forward_192(q, qk, in, 0, 1, delta, qs);
    aes_enc_backward_192(q, qk, 0, 1, delta, out, qs_dash);

    // Step: 13..14
    bf192_t minus_part = bf192_mul(bf192_load(delta), bf192_load(delta));
    for (unsigned int j = 0; j < FAEST_192F_Senc; j++) {
      B[j] = bf192_add(bf192_mul(qs[j], qs_dash[j]), minus_part);
    }
  }
}

static void aes_prove_192(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
                          const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde,
                          uint8_t* b_tilde, const faest_paramset_t* params) {
  // Step: 1..2
  bf192_t* bf_v = column_to_row_major_and_shrink_V_192(V, FAEST_192F_L);

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7
  bf192_t* A0 = malloc(sizeof(bf192_t) * (FAEST_192F_Ske + 2 * FAEST_192F_Senc + 1));
  bf192_t* A1 = malloc(sizeof(bf192_t) * (FAEST_192F_Ske + 2 * FAEST_192F_Senc + 1));
  uint8_t* k  = malloc((FAEST_192F_R + 1) * 128 / 8);
  bf192_t* vk = malloc(sizeof(bf192_t) * ((FAEST_192F_R + 1) * 128));
  bf192_t* qk = malloc(sizeof(bf192_t) * ((FAEST_192F_R + 1) * 128));
  aes_key_schedule_constraints_192(w, bf_v, 0, NULL, NULL, A0, A1, k, vk, NULL, qk, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  aes_enc_constraints_192(in, out, w + FAEST_192F_Lke / 8, bf_v + FAEST_192F_Lke, k, vk, 0, NULL,
                          NULL, NULL, A0 + FAEST_192F_Ske, A1 + FAEST_192F_Ske, NULL);
  // Step: 12-15
  aes_enc_constraints_192(in + 16, out + 16, w + (FAEST_192F_Lke + FAEST_192F_Lenc) / 8,
                          bf_v + FAEST_192F_Lke + FAEST_192F_Lenc, k, vk, 0, NULL, NULL, NULL,
                          A0 + (FAEST_192F_Ske + FAEST_192F_Senc),
                          A1 + (FAEST_192F_Ske + FAEST_192F_Senc), NULL);
  free(qk);
  free(vk);
  free(k);

  // Step: 16..18
  A1[FAEST_192F_Ske + 2 * FAEST_192F_Senc] = bf192_load(u + FAEST_192F_L / 8);
  A0[FAEST_192F_Ske + 2 * FAEST_192F_Senc] = bf192_sum_poly(bf_v + FAEST_192F_L);
  free(bf_v);

  zk_hash_192(a_tilde, chall, A1, FAEST_192F_Ske + 2 * FAEST_192F_Senc);
  zk_hash_192(b_tilde, chall, A0, FAEST_192F_Ske + 2 * FAEST_192F_Senc);

  free(A0);
  free(A1);
}

static uint8_t* aes_verify_192(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                               const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                               const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.t0;
  const unsigned int k0  = params->faest_param.k0;
  const unsigned int t1  = params->faest_param.t1;
  const unsigned int k1  = params->faest_param.k1;

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
        xor_u8_array(d, Q[col], Q[col], (FAEST_192F_L + 7) / 8);
      }
    }
  }

  // Step: 11..12
  bf192_t* bf_q = column_to_row_major_and_shrink_V_192(Q, FAEST_192F_L);

  // Step: 13
  uint8_t* k   = malloc((FAEST_192F_R + 1) * 128);
  bf192_t* vk  = malloc(sizeof(bf192_t) * ((FAEST_192F_R + 1) * 128));
  bf192_t* qk  = malloc(sizeof(bf192_t) * ((FAEST_192F_R + 1) * 128));
  bf192_t* B_0 = malloc(sizeof(bf192_t) * (FAEST_192F_Ske + 2 * FAEST_192F_Senc + 1));
  aes_key_schedule_constraints_192(NULL, NULL, 1, bf_q, delta, NULL, NULL, k, vk, B_0, qk, params);

  // Step: 14
  bf192_t* B_1 = B_0 + FAEST_192F_Ske;
  aes_enc_constraints_192(in, out, NULL, NULL, NULL, NULL, 1, bf_q + FAEST_192F_Lke, qk, delta,
                          NULL, NULL, B_1);

  // Step: 18
  bf192_t* B_2 = B_0 + FAEST_192F_Ske + FAEST_192F_Senc;
  aes_enc_constraints_192(in + 16, out + 16, NULL, NULL, NULL, NULL, 1,
                          bf_q + FAEST_192F_Lke + FAEST_192F_Lenc, qk, delta, NULL, NULL, B_2);
  free(qk);
  free(vk);
  free(k);

  // Step: 20
  B_0[FAEST_192F_Ske + 2 * FAEST_192F_Senc] = bf192_sum_poly(bf_q + FAEST_192F_L);
  free(bf_q);

  // Step 21
  uint8_t* q_tilde = malloc(FAEST_192F_LAMBDA / 8);
  zk_hash_192(q_tilde, chall_2, B_0, FAEST_192F_Ske + 2 * FAEST_192F_Senc);
  free(B_0);

  bf192_t bf_qtilde = bf192_load(q_tilde);
  bf192_store(q_tilde, bf192_add(bf_qtilde, bf192_mul(bf192_load(a_tilde), bf192_load(delta))));

  return q_tilde;
}

// lambda == 256 implementation

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

static void aes_key_schedule_constraints_256(const uint8_t* w, const bf256_t* v, const uint8_t Mkey,
                                             const bf256_t* q, const uint8_t* delta, bf256_t* A0,
                                             bf256_t* A1, uint8_t* k, bf256_t* vk, bf256_t* B,
                                             bf256_t* qk, const faest_paramset_t* params) {
  bool rotate_word = true;

  if (Mkey == 0) {
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
    unsigned int iwd = 32 * (FAEST_256F_Nwd - 1);
    for (unsigned int j = 0; j < FAEST_256F_Ske / 4; j++) {
      bf256_t bf_k_hat[4];
      bf256_t bf_v_k_hat[4];
      bf256_t bf_w_dash_hat[4];
      bf256_t bf_v_w_dash_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 10..11
        if (rotate_word) {
          bf_k_hat[(r + 3) % 4]   = bf256_byte_combine_bits(k[(iwd + 8 * r) / 8]);
          bf_v_k_hat[(r + 3) % 4] = bf256_byte_combine(vk + (iwd + 8 * r));
          bf_w_dash_hat[r]        = bf256_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
          bf_v_w_dash_hat[r]      = bf256_byte_combine(v_w_dash + (32 * j + 8 * r));
        } else {
          bf_k_hat[r]        = bf256_byte_combine_bits(k[(iwd + 8 * r) / 8]);
          bf_v_k_hat[r]      = bf256_byte_combine(vk + (iwd + 8 * r));
          bf_w_dash_hat[r]   = bf256_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
          bf_v_w_dash_hat[r] = bf256_byte_combine(v_w_dash + (32 * j + 8 * r));
        }
      }
      // Step: 13..17
      for (unsigned int r = 0; r <= 3; r++) {
        A0[4 * j + r] = bf256_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
        A1[4 * j + r] =
            bf256_add(bf256_add(bf256_mul(bf256_add(bf_k_hat[r], bf_v_k_hat[r]),
                                          bf256_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                bf256_from_bf8(bf8_one())),
                      A0[4 * j + r]);
      }
      if (FAEST_256F_LAMBDA == 192) {
        iwd = iwd + 192;
      } else {
        iwd = iwd + 128;
        if (FAEST_256F_LAMBDA == 256) {
          rotate_word = !rotate_word;
        }
      }
    }
    return;
  }

  // Step: 19..20
  aes_key_schedule_forward_256(q, qk);
  bf256_t* q_w_dash = malloc(FAEST_256F_Ske * 8 * sizeof(bf256_t));
  aes_key_schedule_backward_256(&q[FAEST_256F_LAMBDA], qk, 0, 1, delta, q_w_dash);

  const bf256_t bf_delta = bf256_load(delta);

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_256F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_256F_Ske / 4; j++) {
    bf256_t bf_q_hat_k[4];
    bf256_t bf_q_hat_w_dash[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 25..26
      if (rotate_word) {
        bf_q_hat_k[(r + 3) % 4] = bf256_byte_combine(qk + ((iwd + 8 * r)));
        bf_q_hat_w_dash[r]      = bf256_byte_combine(q_w_dash + ((32 * j + 8 * r)));
      } else {
        bf_q_hat_k[r]      = bf256_byte_combine(qk + ((iwd + 8 * r)));
        bf_q_hat_w_dash[r] = bf256_byte_combine(q_w_dash + ((32 * j + 8 * r)));
      }
    }
    // Step: 27
    for (unsigned int r = 0; r <= 3; r++) {
      bf256_t bf_tmp = bf256_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      B[4 * j + r]   = bf256_add(bf_tmp, bf256_mul(bf_delta, bf_delta));
    }
    if (FAEST_256F_LAMBDA == 192) {
      iwd = iwd + 192;
    } else {
      iwd = iwd + 128;
      if (FAEST_256F_LAMBDA == 256) {
        rotate_word = !rotate_word;
      }
    }
  }
  free(q_w_dash);
}

static void aes_enc_forward_256_1(const uint8_t* x, const uint8_t* xk, const uint8_t* in,
                                  uint8_t Mtag, uint8_t Mkey, bf256_t* bf_y) {
  // Step: 2
  for (unsigned int i = 0; i < 16; i++) {
    // Step: 3,4 (bit spliced)
    const uint8_t xin = in[i] & -((1 ^ Mtag) & (1 ^ Mkey));
    // Step: 5
    bf_y[i] = bf256_add(bf256_byte_combine_bits(xin), bf256_byte_combine_bits(xk[i]));
  }

  const bf256_t bf_two   = bf256_byte_combine_bits(2);
  const bf256_t bf_three = bf256_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_256F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_xk_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf256_byte_combine_bits(x[(ix + 8 * r) / 8]);
        bf_xk_hat[r] = bf256_byte_combine_bits(xk[(ik + 8 * r) / 8]);
      }

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
    }
  }
}

static void aes_enc_forward_256(const bf256_t* bf_x, const bf256_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf256_t* bf_y) {
  const bf256_t bf_delta      = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t bf_minus_mtag = bf256_from_bit(1 ^ Mtag);
  const bf256_t bf_minus_mkey = bf256_from_bit(1 ^ Mkey);

  // Step: 2..4
  for (unsigned int i = 0; i < 16; i++) {
    bf256_t bf_xin[8];
    for (unsigned int j = 0; j < 8; j++) {
      bf_xin[j] = bf256_mul(bf256_mul_bit(bf_minus_mtag, get_bit(in[i], j)),
                            bf256_add(bf256_mul_bit(bf_delta, Mkey), bf_minus_mkey));
    }
    // Step: 5
    bf_y[i] = bf256_add(bf256_byte_combine(bf_xin), bf256_byte_combine(bf_xk + (8 * i)));
  }

  const bf256_t bf_two   = bf256_byte_combine_bits(2);
  const bf256_t bf_three = bf256_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_256S_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      const unsigned int ix = 128 * (j - 1) + 32 * c;
      const unsigned int ik = 128 * j + 32 * c;
      const unsigned int iy = 16 * j + 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_xk_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf256_byte_combine(bf_x + (ix + 8 * r));
        bf_xk_hat[r] = bf256_byte_combine(bf_xk + (ik + 8 * r));
      }

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
    }
  }
}

static void aes_enc_backward_256_1(const uint8_t* x, const uint8_t* xk, uint8_t Mtag, uint8_t Mkey,
                                   const uint8_t* out, bf256_t* y_out) {
  uint8_t xtilde;
  // Step:2..4
  for (unsigned int j = 0; j < FAEST_256F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 5..6
        unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
        if (j < (FAEST_256F_R - 1)) {
          // Step: 7
          xtilde = x[ird / 8];
        } else {
          // Step: 9..11 (bit spliced)
          uint8_t xout = out[(ird - 128 * (FAEST_256F_R - 1)) / 8] & -((1 ^ Mtag) & (1 ^ Mkey));
          xtilde       = xout ^ xk[(128 + ird) / 8];
        }

        // Step: 12..17 (bit spliced)
        uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^ rotr8(xtilde, 2);
        ytilde ^= set_bit((1 ^ Mtag) & (1 ^ Mkey), 0);
        ytilde ^= set_bit((1 ^ Mtag) & (1 ^ Mkey), 2);

        // Step: 18
        y_out[16 * j + 4 * c + r] = bf256_byte_combine_bits(ytilde);
      }
    }
  }
}

static void aes_enc_backward_256(const bf256_t* bf_x, const bf256_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf256_t* y_out) {
  // Step: 1
  const bf256_t bf_delta = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t factor =
      bf256_mul_bit(bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  // Step: 2..4
  for (unsigned int j = 0; j < FAEST_256F_R; j++) {
    for (unsigned int c = 0; c <= 3; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
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
      }
    }
  }
}

static void aes_enc_constraints_256(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                    const bf256_t* v, const uint8_t* k, const bf256_t* vk,
                                    uint8_t Mkey, const bf256_t* q, const bf256_t* qk,
                                    const uint8_t* delta, bf256_t* A0, bf256_t* A1, bf256_t* B) {
  if (Mkey == 0) {
    bf256_t s[FAEST_256F_Senc];
    bf256_t vs[FAEST_256F_Senc];
    bf256_t s_dash[FAEST_256F_Senc];
    bf256_t vs_dash[FAEST_256F_Senc];
    aes_enc_forward_256_1(w, k, in, 0, 0, s);
    aes_enc_forward_256(v, vk, in, 1, 0, NULL, vs);
    aes_enc_backward_256_1(w, k, 0, 0, out, s_dash);
    aes_enc_backward_256(v, vk, 1, 0, NULL, out, vs_dash);

    for (unsigned int j = 0; j < FAEST_256F_Senc; j++) {
      A0[j] = bf256_mul(vs[j], vs_dash[j]);
      A1[j] = bf256_add(
          bf256_add(bf256_mul(bf256_add(s[j], vs[j]), bf256_add(s_dash[j], vs_dash[j])), A0[j]),
          bf256_one());
    }
  } else {
    // Step: 11..12
    bf256_t qs[FAEST_256F_Senc];
    bf256_t qs_dash[FAEST_256F_Senc];
    aes_enc_forward_256(q, qk, in, 0, 1, delta, qs);
    aes_enc_backward_256(q, qk, 0, 1, delta, out, qs_dash);

    // Step: 13..14
    bf256_t minus_part = bf256_mul(bf256_load(delta), bf256_load(delta));
    for (unsigned int j = 0; j < FAEST_256F_Senc; j++) {
      B[j] = bf256_add(bf256_mul(qs[j], qs_dash[j]), minus_part);
    }
  }
}

static void aes_prove_256(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
                          const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde,
                          uint8_t* b_tilde, const faest_paramset_t* params) {
  // Step: 1..2
  bf256_t* bf_v = column_to_row_major_and_shrink_V_256(V, FAEST_256F_L);

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7
  bf256_t* A0 = malloc(sizeof(bf256_t) * (FAEST_256F_Ske + 2 * FAEST_256F_Senc + 1));
  bf256_t* A1 = malloc(sizeof(bf256_t) * (FAEST_256F_Ske + 2 * FAEST_256F_Senc + 1));
  uint8_t* k  = malloc((FAEST_256F_R + 1) * 128 / 8);
  bf256_t* vk = malloc(sizeof(bf256_t) * ((FAEST_256F_R + 1) * 128));
  bf256_t* qk = malloc(sizeof(bf256_t) * ((FAEST_256F_R + 1) * 128));
  aes_key_schedule_constraints_256(w, bf_v, 0, NULL, NULL, A0, A1, k, vk, NULL, qk, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  aes_enc_constraints_256(in, out, w + FAEST_256F_Lke / 8, bf_v + FAEST_256F_Lke, k, vk, 0, NULL,
                          NULL, NULL, A0 + FAEST_256F_Ske, A1 + FAEST_256F_Ske, NULL);
  // Step: 12-15
  aes_enc_constraints_256(in + 16, out + 16, w + (FAEST_256F_Lke + FAEST_256F_Lenc) / 8,
                          bf_v + FAEST_256F_Lke + FAEST_256F_Lenc, k, vk, 0, NULL, NULL, NULL,
                          A0 + FAEST_256F_Ske + FAEST_256F_Senc,
                          A1 + FAEST_256F_Ske + FAEST_256F_Senc, NULL);
  free(qk);
  free(vk);
  free(k);

  // Step: 16..18
  A1[FAEST_256F_Ske + 2 * FAEST_256F_Senc] = bf256_load(u + FAEST_256F_L / 8);
  A0[FAEST_256F_Ske + 2 * FAEST_256F_Senc] = bf256_sum_poly(bf_v + FAEST_256F_L);
  free(bf_v);

  zk_hash_256(a_tilde, chall, A1, FAEST_256F_Ske + 2 * FAEST_256F_Senc);
  zk_hash_256(b_tilde, chall, A0, FAEST_256F_Ske + 2 * FAEST_256F_Senc);

  free(A0);
  free(A1);
}

static uint8_t* aes_verify_256(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                               const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                               const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.t0;
  const unsigned int k0  = params->faest_param.k0;
  const unsigned int t1  = params->faest_param.t1;
  const unsigned int k1  = params->faest_param.k1;

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
        xor_u8_array(d, Q[col], Q[col], (FAEST_256F_L + 7) / 8);
      }
    }
  }

  // Step: 11..12
  bf256_t* bf_q = column_to_row_major_and_shrink_V_256(Q, FAEST_256F_L);

  // Step: 13
  uint8_t* k   = malloc((FAEST_256F_R + 1) * 128);
  bf256_t* vk  = malloc(sizeof(bf256_t) * ((FAEST_256F_R + 1) * 128));
  bf256_t* qk  = malloc(sizeof(bf256_t) * ((FAEST_256F_R + 1) * 128));
  bf256_t* B_0 = malloc(sizeof(bf256_t) * (FAEST_256F_Ske + 2 * FAEST_256F_Senc + 1));
  aes_key_schedule_constraints_256(NULL, NULL, 1, bf_q, delta, NULL, NULL, k, vk, B_0, qk, params);

  // Step: 14
  bf256_t* B_1 = B_0 + FAEST_256F_Ske;
  aes_enc_constraints_256(in, out, NULL, NULL, NULL, NULL, 1, bf_q + FAEST_256F_Lke, qk, delta,
                          NULL, NULL, B_1);

  // Step: 18
  bf256_t* B_2 = B_0 + (FAEST_256F_Ske + FAEST_256F_Senc);
  aes_enc_constraints_256(in + 16, out + 16, NULL, NULL, NULL, NULL, 1,
                          bf_q + (FAEST_256F_Lke + FAEST_256F_Lenc), qk, delta, NULL, NULL, B_2);
  free(qk);
  free(vk);
  free(k);

  // Step: 20
  B_0[FAEST_256F_Ske + 2 * FAEST_256F_Senc] = bf256_sum_poly(bf_q + FAEST_256F_L);
  free(bf_q);

  // Step 21
  uint8_t* q_tilde = malloc(FAEST_256F_LAMBDA / 8);
  zk_hash_256(q_tilde, chall_2, B_0, FAEST_256F_Ske + 2 * FAEST_256F_Senc);
  free(B_0);

  bf256_t bf_qtilde = bf256_load(q_tilde);
  bf256_store(q_tilde, bf256_add(bf_qtilde, bf256_mul(bf256_load(a_tilde), bf256_load(delta))));

  return q_tilde;
}

// EM-128

static void em_enc_forward_128_1(const uint8_t* z, const uint8_t* x, bf128_t* bf_y) { // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_128F_Nwd; j++) {
    bf_y[j] = bf128_add(bf128_byte_combine_bits(z[j]), bf128_byte_combine_bits(x[j]));
  }

  const bf128_t bf_two   = bf128_byte_combine_bits(2);
  const bf128_t bf_three = bf128_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_128F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_128F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_128F_Nwd * j + 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_z_hat[r] = bf128_byte_combine_bits(z[(i + 8 * r) / 8]);
        bf_x_hat[r] = bf128_byte_combine_bits(x[(i + 8 * r) / 8]);
      }

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
    }
  }
}

static void em_enc_forward_128(const bf128_t* bf_z, const bf128_t* bf_x, bf128_t* bf_y) {
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_128F_Nwd; j++) {
    bf_y[j] = bf128_byte_combine(bf_z + 8 * j);
    if (bf_x) {
      bf_y[j] = bf128_add(bf_y[j], bf128_byte_combine(bf_x + 8 * j));
    }
  }

  const bf128_t bf_two   = bf128_byte_combine_bits(2);
  const bf128_t bf_three = bf128_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_128F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_128F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_128F_Nwd * j + 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_z_hat[r] = bf128_byte_combine(bf_z + (i + 8 * r));
        if (bf_x) {
          bf_x_hat[r] = bf128_byte_combine(bf_x + (i + 8 * r));
        } else {
          bf_x_hat[r] = bf128_zero();
        }
      }

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
    }
  }
}

static void em_enc_backward_128_1(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf128_t* y_out) {
  // only called with Mtag == Mkey == 0

  for (unsigned int j = 0; j < FAEST_EM_128F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
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
        // delta is always bot
        // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
        const uint8_t y_tilde = rotr8(z_tilde, 7) ^ rotr8(z_tilde, 5) ^ rotr8(z_tilde, 2) ^ 0x5;

        // Step: 18
        y_out[4 * FAEST_EM_128F_Nwd * j + 4 * c + r] = bf128_byte_combine_bits(y_tilde);
      }
    }
  }
}

static void em_enc_backward_128(const bf128_t* bf_z, const bf128_t* bf_x, const bf128_t* bf_z_out,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf128_t* y_out) {
  // Step: 1
  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t factor =
      bf128_mul_bit(bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int j = 0; j < FAEST_EM_128F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
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
      }
    }
  }
}

static void em_enc_constraints_128(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                   const bf128_t* bf_v, uint8_t Mkey, const bf128_t* bf_q,
                                   const uint8_t* delta, bf128_t* A0, bf128_t* A1, bf128_t* B) {
  if (Mkey == 0) {
    // Step 6
    uint8_t w_out[FAEST_EM_128F_LAMBDA / 8];
    xor_u8_array(out, w, w_out, sizeof(w_out));

    bf128_t bf_s[FAEST_EM_128F_Senc];
    bf128_t bf_vs[FAEST_EM_128F_Senc];
    bf128_t bf_s_dash[FAEST_EM_128F_Senc];
    bf128_t bf_vs_dash[FAEST_EM_128F_Senc];
    em_enc_forward_128_1(w, x, bf_s);
    em_enc_forward_128(bf_v, NULL, bf_vs);
    em_enc_backward_128_1(w, x, w_out, bf_s_dash);
    em_enc_backward_128(bf_v, NULL, bf_v, 1, 0, NULL, bf_vs_dash);

    for (unsigned int j = 0; j < FAEST_EM_128F_Senc; j++) {
      A0[j] = bf128_mul(bf_vs[j], bf_vs_dash[j]);
      A1[j] = bf128_add(
          bf128_add(bf128_mul(bf128_add(bf_s[j], bf_vs[j]), bf128_add(bf_s_dash[j], bf_vs_dash[j])),
                    A0[j]),
          bf128_one());
    }
  } else {
    // Step: 18, 19
    // TODO: compute these on demand in em_enc_backward_128
    const bf128_t bf_delta = bf128_load(delta);
    bf128_t* bf_x          = malloc(sizeof(bf128_t) * 128 * (FAEST_EM_128F_R + 1));
    for (unsigned int i = 0; i < 128 * (FAEST_EM_128F_R + 1); i++) {
      bf_x[i] = bf128_mul_bit(bf_delta, ptr_get_bit(x, i));
    }

    // Step 21
    bf128_t* bf_q_out = malloc(sizeof(bf128_t) * FAEST_EM_128F_LAMBDA);
    for (unsigned int i = 0; i < FAEST_EM_128F_LAMBDA; i++) {
      bf_q_out[i] = bf128_add(bf128_mul_bit(bf_delta, ptr_get_bit(out, i)), bf_q[i]);
    }

    bf128_t bf_qs[FAEST_EM_128F_Senc];
    bf128_t bf_qs_dash[FAEST_EM_128F_Senc];
    em_enc_forward_128(bf_q, bf_x, bf_qs);
    em_enc_backward_128(bf_q, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash);
    free(bf_q_out);

    // Step: 13..14
    bf128_t minus_part = bf128_mul(bf_delta, bf_delta);
    for (unsigned int j = 0; j < FAEST_EM_128F_Senc; j++) {
      B[j] = bf128_add(bf128_mul(bf_qs[j], bf_qs_dash[j]), minus_part);
    }
    free(bf_x);
  }
}

static void em_prove_128(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
                         const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde,
                         uint8_t* b_tilde) {
  // copy expanded key in to an array
  uint8_t* x = malloc(FAEST_EM_128F_LAMBDA * (FAEST_EM_128F_R + 1) / 8);
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

  bf128_t* A0 = malloc(sizeof(bf128_t) * (FAEST_EM_128F_Senc + 1));
  bf128_t* A1 = malloc(sizeof(bf128_t) * (FAEST_EM_128F_Senc + 1));
  em_enc_constraints_128(out, x, w, bf_v, 0, NULL, NULL, A0, A1, NULL);
  free(x);

  A1[FAEST_EM_128F_Senc] = bf128_load(u + FAEST_EM_128F_Lenc / 8);
  A0[FAEST_EM_128F_Senc] = bf128_sum_poly(bf_v + FAEST_EM_128F_Lenc);
  free(bf_v);

  zk_hash_128(a_tilde, chall, A1, FAEST_EM_128F_Senc);
  zk_hash_128(b_tilde, chall, A0, FAEST_EM_128F_Senc);

  free(A0);
  free(A1);
}

static uint8_t* em_verify_128(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                              const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                              const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.t0;
  const unsigned int k0  = params->faest_param.k0;
  const unsigned int t1  = params->faest_param.t1;
  const unsigned int k1  = params->faest_param.k1;

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
  uint8_t* x = malloc(FAEST_EM_128F_LAMBDA * (FAEST_EM_128F_R + 1) / 8);
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

  bf128_t* B = malloc(sizeof(bf128_t) * (FAEST_EM_128F_Senc + 1));

  em_enc_constraints_128(out, x, NULL, NULL, 1, bf_q, delta, NULL, NULL, B);
  free(x);

  B[FAEST_EM_128F_Senc] = bf128_sum_poly(bf_q + FAEST_EM_128F_Lenc);
  free(bf_q);

  uint8_t* q_tilde = malloc(FAEST_EM_128F_LAMBDA / 8);
  zk_hash_128(q_tilde, chall_2, B, FAEST_EM_128F_Senc);
  free(B);

  bf128_t bf_qtilde = bf128_load(q_tilde);
  bf128_store(q_tilde, bf128_add(bf_qtilde, bf128_mul(bf128_load(a_tilde), bf128_load(delta))));

  return q_tilde;
}

// EM-192

static void em_enc_forward_192_1(const uint8_t* z, const uint8_t* x, bf192_t* bf_y) {
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_192F_Nwd; j++) {
    bf_y[j] = bf192_add(bf192_byte_combine_bits(z[j]), bf192_byte_combine_bits(x[j]));
  }

  const bf192_t bf_two   = bf192_byte_combine_bits(2);
  const bf192_t bf_three = bf192_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_192F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_192F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_192F_Nwd * j + 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_z_hat[r] = bf192_byte_combine_bits(z[(i + 8 * r) / 8]);
        bf_x_hat[r] = bf192_byte_combine_bits(x[(i + 8 * r) / 8]);
      }

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
    }
  }
}

static void em_enc_forward_192(const bf192_t* bf_z, const bf192_t* bf_x, bf192_t* bf_y) {
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_192F_Nwd; j++) {
    bf_y[j] = bf192_byte_combine(bf_z + 8 * j);
    if (bf_x) {
      bf_y[j] = bf192_add(bf_y[j], bf192_byte_combine(bf_x + 8 * j));
    }
  }

  const bf192_t bf_two   = bf192_byte_combine_bits(2);
  const bf192_t bf_three = bf192_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_192F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_192F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_192F_Nwd * j + 4 * c;

      bf192_t bf_x_hat[4];
      bf192_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_z_hat[r] = bf192_byte_combine(bf_z + (i + 8 * r));
        if (bf_x) {
          bf_x_hat[r] = bf192_byte_combine(bf_x + (i + 8 * r));
        } else {
          bf_x_hat[r] = bf192_zero();
        }
      }

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
    }
  }
}

static void em_enc_backward_192_1(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf192_t* y_out) {
  // only called with Mtag == Mkey == 0

  for (unsigned int j = 0; j < FAEST_EM_192F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
        const unsigned int icol = (c - r + FAEST_EM_192F_Nwd) % FAEST_EM_192F_Nwd;
        const unsigned int ird =
            FAEST_EM_192F_LAMBDA + 32 * FAEST_EM_192F_Nwd * j + 32 * icol + 8 * r;
        uint8_t z_tilde = 0;
        if (j < (FAEST_EM_192F_R - 1)) {
          z_tilde = z[ird / 8];
        } else {
          z_tilde = z_out[(ird - 32 * FAEST_EM_192F_Nwd * (j + 1)) / 8] ^ x[ird / 8];
        }

        // delta is always bot
        // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
        const uint8_t y_tilde = rotr8(z_tilde, 7) ^ rotr8(z_tilde, 5) ^ rotr8(z_tilde, 2) ^ 0x5;

        // Step: 18
        y_out[4 * FAEST_EM_192F_Nwd * j + 4 * c + r] = bf192_byte_combine_bits(y_tilde);
      }
    }
  }
}

static void em_enc_backward_192(const bf192_t* bf_z, const bf192_t* bf_x, const bf192_t* bf_z_out,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf192_t* y_out) {
  // Step: 1
  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t factor =
      bf192_mul_bit(bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int j = 0; j < FAEST_EM_192F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
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
      }
    }
  }
}

static void em_enc_constraints_192(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                   const bf192_t* bf_v, uint8_t Mkey, const bf192_t* bf_q,
                                   const uint8_t* delta, bf192_t* A0, bf192_t* A1, bf192_t* B) {
  if (Mkey == 0) {
    // Step 6
    uint8_t w_out[FAEST_EM_192F_LAMBDA / 8];
    xor_u8_array(out, w, w_out, sizeof(w_out));

    bf192_t bf_s[FAEST_EM_192F_Senc];
    bf192_t bf_vs[FAEST_EM_192F_Senc];
    bf192_t bf_s_dash[FAEST_EM_192F_Senc];
    bf192_t bf_vs_dash[FAEST_EM_192F_Senc];
    em_enc_forward_192_1(w, x, bf_s);
    em_enc_forward_192(bf_v, NULL, bf_vs);
    em_enc_backward_192_1(w, x, w_out, bf_s_dash);
    em_enc_backward_192(bf_v, NULL, bf_v, 1, 0, NULL, bf_vs_dash);

    for (unsigned int j = 0; j < FAEST_EM_192F_Senc; j++) {
      A0[j] = bf192_mul(bf_vs[j], bf_vs_dash[j]);
      A1[j] = bf192_add(
          bf192_add(bf192_mul(bf192_add(bf_s[j], bf_vs[j]), bf192_add(bf_s_dash[j], bf_vs_dash[j])),
                    A0[j]),
          bf192_one());
    }
  } else {
    // Step: 18, 19
    // TODO: compute these on demand in em_enc_backward_192
    const bf192_t bf_delta = bf192_load(delta);
    bf192_t* bf_x          = malloc(sizeof(bf192_t) * 192 * (FAEST_EM_192F_R + 1));
    for (unsigned int i = 0; i < 192 * (FAEST_EM_192F_R + 1); i++) {
      bf_x[i] = bf192_mul_bit(bf_delta, ptr_get_bit(x, i));
    }

    // Step 21
    bf192_t* bf_q_out = malloc(sizeof(bf192_t) * FAEST_EM_192F_LAMBDA);
    for (unsigned int i = 0; i < FAEST_EM_192F_LAMBDA; i++) {
      bf_q_out[i] = bf192_add(bf192_mul_bit(bf_delta, ptr_get_bit(out, i)), bf_q[i]);
    }

    bf192_t bf_qs[FAEST_EM_192F_Senc];
    bf192_t bf_qs_dash[FAEST_EM_192F_Senc];
    em_enc_forward_192(bf_q, bf_x, bf_qs);
    em_enc_backward_192(bf_q, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash);
    free(bf_q_out);

    // Step: 13..14
    bf192_t minus_part = bf192_mul(bf_delta, bf_delta);
    for (unsigned int j = 0; j < FAEST_EM_192F_Senc; j++) {
      B[j] = bf192_add(bf192_mul(bf_qs[j], bf_qs_dash[j]), minus_part);
    }
    free(bf_x);
  }
}

static void em_prove_192(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
                         const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde,
                         uint8_t* b_tilde) {
  // copy expanded key in to an array
  uint8_t* x = malloc(FAEST_EM_192F_LAMBDA * (FAEST_EM_192F_R + 1) / 8);
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

  bf192_t* A0 = malloc(sizeof(bf192_t) * (FAEST_EM_192F_Senc + 1));
  bf192_t* A1 = malloc(sizeof(bf192_t) * (FAEST_EM_192F_Senc + 1));
  em_enc_constraints_192(out, x, w, bf_v, 0, NULL, NULL, A0, A1, NULL);
  free(x);

  A1[FAEST_EM_192F_Senc] = bf192_load(u + FAEST_EM_192F_Lenc / 8);
  A0[FAEST_EM_192F_Senc] = bf192_sum_poly(bf_v + FAEST_EM_192F_Lenc);
  free(bf_v);

  zk_hash_192(a_tilde, chall, A1, FAEST_EM_192F_Senc);
  zk_hash_192(b_tilde, chall, A0, FAEST_EM_192F_Senc);

  free(A0);
  free(A1);
}

static uint8_t* em_verify_192(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                              const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                              const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.t0;
  const unsigned int k0  = params->faest_param.k0;
  const unsigned int t1  = params->faest_param.t1;
  const unsigned int k1  = params->faest_param.k1;

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
  uint8_t* x = malloc(FAEST_EM_192F_LAMBDA * (FAEST_EM_192F_R + 1) / 8);
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

  bf192_t* B = malloc(sizeof(bf192_t) * (FAEST_EM_192F_Senc + 1));

  em_enc_constraints_192(out, x, NULL, NULL, 1, bf_q, delta, NULL, NULL, B);
  free(x);

  B[FAEST_EM_192F_Senc] = bf192_sum_poly(bf_q + FAEST_EM_192F_Lenc);
  free(bf_q);

  uint8_t* q_tilde = malloc(FAEST_EM_192F_LAMBDA / 8);
  zk_hash_192(q_tilde, chall_2, B, FAEST_EM_192F_Senc);
  free(B);

  bf192_t bf_qtilde = bf192_load(q_tilde);
  bf192_store(q_tilde, bf192_add(bf_qtilde, bf192_mul(bf192_load(a_tilde), bf192_load(delta))));

  return q_tilde;
}

// EM-256

static void em_enc_forward_256_1(const uint8_t* z, const uint8_t* x, bf256_t* bf_y) {
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_256F_Nwd; j++) {
    bf_y[j] = bf256_add(bf256_byte_combine_bits(z[j]), bf256_byte_combine_bits(x[j]));
  }

  const bf256_t bf_two   = bf256_byte_combine_bits(2);
  const bf256_t bf_three = bf256_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_256F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_256F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_256F_Nwd * j + 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_z_hat[r] = bf256_byte_combine_bits(z[(i + 8 * r) / 8]);
        bf_x_hat[r] = bf256_byte_combine_bits(x[(i + 8 * r) / 8]);
      }

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
    }
  }
}

static void em_enc_forward_256(const bf256_t* bf_z, const bf256_t* bf_x, bf256_t* bf_y) {
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_256F_Nwd; j++) {
    bf_y[j] = bf256_byte_combine(bf_z + 8 * j);
    if (bf_x) {
      bf_y[j] = bf256_add(bf_y[j], bf256_byte_combine(bf_x + 8 * j));
    }
  }

  const bf256_t bf_two   = bf256_byte_combine_bits(2);
  const bf256_t bf_three = bf256_byte_combine_bits(3);

  for (unsigned int j = 1; j < FAEST_EM_256F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {
      const unsigned int i  = 32 * FAEST_EM_256F_Nwd * j + 32 * c;
      const unsigned int iy = 4 * FAEST_EM_256F_Nwd * j + 4 * c;

      bf256_t bf_x_hat[4];
      bf256_t bf_z_hat[4];
      for (unsigned int r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_z_hat[r] = bf256_byte_combine(bf_z + (i + 8 * r));
        if (bf_x) {
          bf_x_hat[r] = bf256_byte_combine(bf_x + (i + 8 * r));
        } else {
          bf_x_hat[r] = bf256_zero();
        }
      }

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
    }
  }
}

static void em_enc_backward_256_1(const uint8_t* z, const uint8_t* x, const uint8_t* z_out,
                                  bf256_t* y_out) {
  // only called with Mtag == Mkey == 0

  for (unsigned int j = 0; j < FAEST_EM_256F_R; j++) {
    for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {
      for (unsigned int r = 0; r <= 3; r++) {
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
        // delta is always bot
        // set_bit((1 ^ Mtag) & (1 ^ Mkey), 0) ^ set_bit((1 ^ Mtag) & (1 ^ Mkey), 2) == 0x5
        const uint8_t y_tilde = rotr8(z_tilde, 7) ^ rotr8(z_tilde, 5) ^ rotr8(z_tilde, 2) ^ 0x5;

        // Step: 18
        y_out[4 * FAEST_EM_256F_Nwd * j + 4 * c + r] = bf256_byte_combine_bits(y_tilde);
      }
    }
  }
}

static void em_enc_backward_256(const bf256_t* bf_z, const bf256_t* bf_x, const bf256_t* bf_z_out,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf256_t* y_out) {
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
      }
    }
  }
}

static void em_enc_constraints_256(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                   const bf256_t* bf_v, uint8_t Mkey, const bf256_t* bf_q,
                                   const uint8_t* delta, bf256_t* A0, bf256_t* A1, bf256_t* B) {
  if (Mkey == 0) {
    // Step 6
    uint8_t w_out[FAEST_EM_256F_LAMBDA / 8];
    xor_u8_array(out, w, w_out, sizeof(w_out));

    bf256_t bf_s[FAEST_EM_256F_Senc];
    bf256_t bf_vs[FAEST_EM_256F_Senc];
    bf256_t bf_s_dash[FAEST_EM_256F_Senc];
    bf256_t bf_vs_dash[FAEST_EM_256F_Senc];
    em_enc_forward_256_1(w, x, bf_s);
    em_enc_forward_256(bf_v, NULL, bf_vs);
    em_enc_backward_256_1(w, x, w_out, bf_s_dash);
    em_enc_backward_256(bf_v, NULL, bf_v, 1, 0, NULL, bf_vs_dash);

    for (unsigned int j = 0; j < FAEST_EM_256F_Senc; j++) {
      A0[j] = bf256_mul(bf_vs[j], bf_vs_dash[j]);
      A1[j] = bf256_add(
          bf256_add(bf256_mul(bf256_add(bf_s[j], bf_vs[j]), bf256_add(bf_s_dash[j], bf_vs_dash[j])),
                    A0[j]),
          bf256_one());
    }
  } else {
    // Step: 18, 19
    // TODO: compute these on demand in em_enc_backward_256
    const bf256_t bf_delta = bf256_load(delta);
    bf256_t* bf_x          = malloc(sizeof(bf256_t) * 256 * (FAEST_EM_256F_R + 1));
    for (unsigned int i = 0; i < 256 * (FAEST_EM_256F_R + 1); i++) {
      bf_x[i] = bf256_mul_bit(bf_delta, ptr_get_bit(x, i));
    }

    // Step 21
    bf256_t* bf_q_out = malloc(sizeof(bf256_t) * FAEST_EM_256F_LAMBDA);
    for (unsigned int i = 0; i < FAEST_EM_256F_LAMBDA; i++) {
      bf_q_out[i] = bf256_add(bf256_mul_bit(bf_delta, ptr_get_bit(out, i)), bf_q[i]);
    }

    bf256_t bf_qs[FAEST_EM_256F_Senc];
    bf256_t bf_qs_dash[FAEST_EM_256F_Senc];
    em_enc_forward_256(bf_q, bf_x, bf_qs);
    em_enc_backward_256(bf_q, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash);
    free(bf_q_out);

    // Step: 13..14
    bf256_t minus_part = bf256_mul(bf_delta, bf_delta);
    for (unsigned int j = 0; j < FAEST_EM_256F_Senc; j++) {
      B[j] = bf256_add(bf256_mul(bf_qs[j], bf_qs_dash[j]), minus_part);
    }
    free(bf_x);
  }
}

static void em_prove_256(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
                         const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde,
                         uint8_t* b_tilde) {
  // copy expanded key in to an array
  uint8_t* x = malloc(FAEST_EM_256F_LAMBDA * (FAEST_EM_256F_R + 1) / 8);
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

  bf256_t* bf_v = column_to_row_major_and_shrink_V_256(V, FAEST_EM_256F_Lenc);

  bf256_t* A0 = malloc(sizeof(bf256_t) * (FAEST_EM_256F_Senc + 1));
  bf256_t* A1 = malloc(sizeof(bf256_t) * (FAEST_EM_256F_Senc + 1));
  em_enc_constraints_256(out, x, w, bf_v, 0, NULL, NULL, A0, A1, NULL);
  free(x);

  A1[FAEST_EM_256F_Senc] = bf256_load(u + FAEST_EM_256F_Lenc / 8);
  A0[FAEST_EM_256F_Senc] = bf256_sum_poly(bf_v + FAEST_EM_256F_Lenc);
  free(bf_v);

  zk_hash_256(a_tilde, chall, A1, FAEST_EM_256F_Senc);
  zk_hash_256(b_tilde, chall, A0, FAEST_EM_256F_Senc);

  free(A0);
  free(A1);
}

static uint8_t* em_verify_256(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                              const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                              const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int tau = params->faest_param.tau;
  const unsigned int t0  = params->faest_param.t0;
  const unsigned int k0  = params->faest_param.k0;
  const unsigned int t1  = params->faest_param.t1;
  const unsigned int k1  = params->faest_param.k1;

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
  uint8_t* x = malloc(FAEST_EM_256F_LAMBDA * (FAEST_EM_256F_R + 1) / 8);
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

  bf256_t* B = malloc(sizeof(bf256_t) * (FAEST_EM_256F_Senc + 1));

  em_enc_constraints_256(out, x, NULL, NULL, 1, bf_q, delta, NULL, NULL, B);
  free(x);

  B[FAEST_EM_256F_Senc] = bf256_sum_poly(bf_q + FAEST_EM_256F_Lenc);
  free(bf_q);

  uint8_t* q_tilde = malloc(FAEST_EM_256F_LAMBDA / 8);
  zk_hash_256(q_tilde, chall_2, B, FAEST_EM_256F_Senc);
  free(B);

  bf256_t bf_qtilde = bf256_load(q_tilde);
  bf256_store(q_tilde, bf256_add(bf_qtilde, bf256_mul(bf256_load(a_tilde), bf256_load(delta))));

  return q_tilde;
}

// dispatchers

void aes_prove(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
               const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params) {
  switch (params->faest_param.lambda) {
  case 256:
    if (params->faest_param.Lke) {
      aes_prove_256(w, u, V, in, out, chall, a_tilde, b_tilde, params);
    } else {
      em_prove_256(w, u, V, in, out, chall, a_tilde, b_tilde);
    }
    break;
  case 192:
    if (params->faest_param.Lke) {
      aes_prove_192(w, u, V, in, out, chall, a_tilde, b_tilde, params);
    } else {
      em_prove_192(w, u, V, in, out, chall, a_tilde, b_tilde);
    }
    break;
  default:
    if (params->faest_param.Lke) {
      aes_prove_128(w, u, V, in, out, chall, a_tilde, b_tilde, params);
    } else {
      em_prove_128(w, u, V, in, out, chall, a_tilde, b_tilde);
    }
  }
}

uint8_t* aes_verify(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2, const uint8_t* chall_3,
                    const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out,
                    const faest_paramset_t* params) {
  switch (params->faest_param.lambda) {
  case 256:
    if (params->faest_param.Lke) {
      return aes_verify_256(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    } else {
      return em_verify_256(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    }
  case 192:
    if (params->faest_param.Lke) {
      return aes_verify_192(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    } else {
      return em_verify_192(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    }
  default:
    if (params->faest_param.Lke) {
      return aes_verify_128(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    } else {
      return em_verify_128(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    }
  }
}
