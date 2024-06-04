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

bf128_t* column_to_row_major_and_shrink_V_128(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf128_t* new_v = faest_aligned_alloc(BF128_ALIGN, (ell + FAEST_128F_LAMBDA) * sizeof(bf128_t));
  for (unsigned int row = 0; row != ell + FAEST_128F_LAMBDA; ++row) {
    uint8_t new_row[BF128_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_128F_LAMBDA; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf128_load(new_row);
  }

  return new_v;
}

bf192_t* column_to_row_major_and_shrink_V_192(uint8_t** v, unsigned int ell) {
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

bf256_t* column_to_row_major_and_shrink_V_256(uint8_t** v, unsigned int ell) {
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

static void aes_key_schedule_backward_128_vbb_vk(vbb_t* vbb, uint8_t Mtag, uint8_t Mkey,
                                                 const uint8_t* delta, bf128_t* bf_out) {
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
      bf_x_tilde[i] = bf128_add(*get_vole_aes_128(vbb, (8 * j + i) + FAEST_128F_LAMBDA),
                                *get_vk_128(vbb, iwd + 8 * c + i));
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

static void aes_key_schedule_constraints_Mkey_0_128(const uint8_t* w, vbb_t* vbb,
                                                    zk_hash_128_ctx* a0_ctx,
                                                    zk_hash_128_ctx* a1_ctx, uint8_t* k,
                                                    const faest_paramset_t* params) {
  // for scan-build
  assert(FAEST_128F_Ske == params->faest_param.Ske);

  // Step: 2
  aes_key_schedule_forward_1(w, k, params);

  // Step: 3
  // aes_key_schedule_forward_128_vbb(vbb, vk);

  // Step: 4
  uint8_t w_dash[FAEST_128F_Ske];
  aes_key_schedule_backward_1(w + FAEST_128F_LAMBDA / 8, k, w_dash, params);

  // Step: 5
  bf128_t v_w_dash[FAEST_128F_Ske * 8];
  aes_key_schedule_backward_128_vbb_vk(vbb, 1, 0, NULL, v_w_dash);

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
      bf_v_k_hat[(r + 3) % 4] = bf128_byte_combine_vk(vbb, (iwd + 8 * r));
      bf_w_dash_hat[r]        = bf128_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
      bf_v_w_dash_hat[r]      = bf128_byte_combine(v_w_dash + (32 * j + 8 * r));
    }
    // Step: 13..17
    for (unsigned int r = 0; r <= 3; r++) {
      // instead of storing in A0, A1, hash it
      const bf128_t tmp = bf128_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
      zk_hash_128_update(a0_ctx, tmp);
      zk_hash_128_update(
          a1_ctx, bf128_add(bf128_add(bf128_mul(bf128_add(bf_k_hat[r], bf_v_k_hat[r]),
                                                bf128_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                      bf128_one()),
                            tmp));
    }
    iwd = iwd + 128;
  }
}

static void aes_key_schedule_constraints_Mkey_1_128(vbb_t* vbb, const uint8_t* delta,
                                                    zk_hash_128_ctx* b0_ctx) {
  // Step: 19..20
  // aes_key_schedule_forward_128_vbb(vbb, qk);
  bf128_t q_w_dash[FAEST_128F_Ske * 8];
  aes_key_schedule_backward_128_vbb_vk(vbb, 0, 1, delta, q_w_dash);

  const bf128_t bf_delta      = bf128_load(delta);
  const bf128_t delta_squared = bf128_mul(bf_delta, bf_delta);

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_128F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_128F_Ske / 4; j++) {
    bf128_t bf_q_hat_k[4];
    bf128_t bf_q_hat_w_dash[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 25..26
      bf_q_hat_k[(r + 3) % 4] = bf128_byte_combine_vk(vbb, ((iwd + 8 * r)));
      bf_q_hat_w_dash[r]      = bf128_byte_combine(q_w_dash + ((32 * j + 8 * r)));
    }
    // Step: 27
    for (unsigned int r = 0; r <= 3; r++) {
      bf128_t bf_tmp = bf128_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      // instead of storing B, hash it
      zk_hash_128_update(b0_ctx, bf128_add(bf_tmp, delta_squared));
    }
    iwd = iwd + 128;
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

static void aes_enc_forward_128_vbb_vk(vbb_t* vbb, unsigned int offset, const uint8_t* in,
                                       uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                       bf128_t* bf_y) {
  const bf128_t bf_delta  = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t bf_factor = bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey));

  // Step: 2..4
  for (unsigned int i = 0; i < 16; i++) {
    bf128_t bf_xin[8];
    for (unsigned int j = 0; j < 8; j++) {
      bf_xin[j] = bf128_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
    }
    // Step: 5
    bf_y[i] = bf128_add(bf128_byte_combine(bf_xin), bf128_byte_combine_vk(vbb, (8 * i)));
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
        bf_x_hat[r]  = bf128_byte_combine_vbb(vbb, offset + ix + 8 * r);
        bf_xk_hat[r] = bf128_byte_combine_vk(vbb, (ik + 8 * r));
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

static void aes_enc_backward_128_vbb_linear_access(vbb_t* vbb, unsigned int offset, uint8_t Mtag,
                                                   uint8_t Mkey, const uint8_t* delta,
                                                   const uint8_t* out, bf128_t* y_out) {
  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t factor =
      bf128_mul_bit(bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)), 1 ^ Mtag);
  // Access pattern friendly for OLEs
  bf128_t bf_x_tilde[8];
  for (unsigned int i = 0; i < (FAEST_128F_R - 1) * 8 * 16; i++) {
    unsigned int j = i / 128;
    unsigned int r = ((i - (128 * j)) % 32) / 8;
    unsigned int c = ((i - 128 * j - 8 * r) / 32 + r) % 4;

    memcpy(bf_x_tilde + (i % 8), get_vole_aes_128(vbb, i + offset), sizeof(bf128_t));

    if (i % 8 == 7) {
      bf128_t bf_y_tilde[8];
      for (unsigned int k = 0; k < 8; ++k) {
        bf_y_tilde[k] = bf128_add(bf128_add(bf_x_tilde[(k + 7) % 8], bf_x_tilde[(k + 5) % 8]),
                                  bf_x_tilde[(k + 2) % 8]);
      }
      bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
      bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

      // Step: 18
      y_out[16 * j + 4 * c + r] = bf128_byte_combine(bf_y_tilde);
    }
  }

  unsigned int j = FAEST_128F_R - 1;
  for (unsigned int c = 0; c <= 3; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      memset(bf_x_tilde, 0, sizeof(bf_x_tilde));
      // Step: 5
      unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);

      // Step: 10
      for (unsigned int i = 0; i < 8; ++i) {
        // Step: 11
        bf128_t bf_xout =
            bf128_mul_bit(factor, get_bit(out[(ird - 128 * (FAEST_128F_R - 1)) / 8], i));
        // Step: 12
        bf_x_tilde[i] = bf128_add(bf_xout, *get_vk_128(vbb, 128 + ird + i));
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

static void aes_enc_constraints_Mkey_0_128(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                           vbb_t* vbb, unsigned int offset, const uint8_t* k,
                                           zk_hash_128_ctx* a0_ctx, zk_hash_128_ctx* a1_ctx) {
  unsigned int w_offset = offset / 8;
  w += w_offset;
  bf128_t s[FAEST_128F_Senc];
  bf128_t vs[FAEST_128F_Senc];
  bf128_t s_dash[FAEST_128F_Senc];
  bf128_t vs_dash[FAEST_128F_Senc];
  aes_enc_forward_128_1(w, k, in, s);
  aes_enc_forward_128_vbb_vk(vbb, offset, in, 1, 0, NULL, vs);
  aes_enc_backward_128_1(w, k, out, s_dash);
  aes_enc_backward_128_vbb_linear_access(vbb, offset, 1, 0, NULL, out, vs_dash);

  for (unsigned int j = 0; j < FAEST_128F_Senc; j++) {
    // instead of storing in A0, A!, hash it
    const bf128_t tmp = bf128_mul(vs[j], vs_dash[j]);
    zk_hash_128_update(a0_ctx, tmp);
    zk_hash_128_update(a1_ctx, bf128_add(bf128_add(bf128_mul(bf128_add(s[j], vs[j]),
                                                             bf128_add(s_dash[j], vs_dash[j])),
                                                   tmp),
                                         bf128_one()));
  }
}

static void aes_enc_constraints_Mkey_1_128(const uint8_t* in, const uint8_t* out, vbb_t* vbb,
                                           unsigned int offset, const uint8_t* delta,
                                           zk_hash_128_ctx* b0_ctx) {

  // Step: 11..12
  bf128_t qs[FAEST_128F_Senc];
  bf128_t qs_dash[FAEST_128F_Senc];
  aes_enc_forward_128_vbb_vk(vbb, offset, in, 0, 1, delta, qs);
  aes_enc_backward_128_vbb_linear_access(vbb, offset, 0, 1, delta, out, qs_dash);

  // Step: 13..14
  bf128_t minus_part = bf128_mul(bf128_load(delta), bf128_load(delta));
  for (unsigned int j = 0; j < FAEST_128F_Senc; j++) {
    // instead of storing it, hash it
    zk_hash_128_update(b0_ctx, bf128_add(bf128_mul(qs[j], qs_dash[j]), minus_part));
  }
}

static void aes_prove_128(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                          const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
                          const faest_paramset_t* params) {
  // Step: 1..2
  // compute on the fly

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7 + 18
  uint8_t* k = malloc((FAEST_128F_R + 1) * 128 / 8);
  // bf128_t* vk = faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * ((FAEST_128F_R + 1) * 128));
  zk_hash_128_ctx a0_ctx;
  zk_hash_128_ctx a1_ctx;

  zk_hash_128_init(&a0_ctx, chall);
  zk_hash_128_init(&a1_ctx, chall);
  aes_key_schedule_constraints_Mkey_0_128(w, vbb, &a0_ctx, &a1_ctx, k, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  unsigned int offset = FAEST_128F_Lke;
  aes_enc_constraints_Mkey_0_128(in, out, w, vbb, offset, k, &a0_ctx, &a1_ctx);
  // Step: 12 (beta == 1)
  // faest_aligned_free(vk);
  free(k);

  // Step: 16..18
  zk_hash_128_finalize(a_tilde, &a1_ctx, bf128_load(get_vole_u(vbb) + FAEST_128F_L / 8));
  zk_hash_128_finalize(b_tilde, &a0_ctx, bf128_sum_poly_vbb(vbb, FAEST_128F_L));
}

static uint8_t* aes_verify_128(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                               const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out) {
  // Step: 1
  const uint8_t* delta = chall_3;

  // Step: 13 + 21
  // bf128_t* qk = faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * ((FAEST_128F_R + 1) * 128));
  // instead of storing B_0 in an array, we process the values with zk_hash_128
  zk_hash_128_ctx b0_ctx;
  zk_hash_128_init(&b0_ctx, chall_2);
  aes_key_schedule_constraints_Mkey_1_128(vbb, delta, &b0_ctx);

  // Step: 14
  aes_enc_constraints_Mkey_1_128(in, out, vbb, FAEST_128F_Lke, delta, &b0_ctx);
  // Step: 18 (beta == 1)
  // faest_aligned_free(qk);

  // Step: 20+21
  uint8_t* q_tilde = malloc(FAEST_128F_LAMBDA / 8);
  zk_hash_128_finalize(q_tilde, &b0_ctx, bf128_sum_poly_vbb(vbb, FAEST_128F_L));

  bf128_t bf_qtilde = bf128_load(q_tilde);
  bf128_store(q_tilde, bf128_add(bf_qtilde, bf128_mul(bf128_load(a_tilde), bf128_load(delta))));

  return q_tilde;
}

// lambda == 192 implementation

static void aes_key_schedule_backward_192_vbb_vk(vbb_t* vbb, uint8_t Mtag, uint8_t Mkey,
                                                 const uint8_t* delta, bf192_t* bf_out) {
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
      bf_x_tilde[i] = bf192_add(*get_vole_aes_192(vbb, (8 * j + i) + FAEST_192F_LAMBDA),
                                *get_vk_192(vbb, iwd + 8 * c + i));
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

static void aes_key_schedule_constraints_Mkey_0_192(const uint8_t* w, vbb_t* vbb,
                                                    zk_hash_192_ctx* a0_ctx,
                                                    zk_hash_192_ctx* a1_ctx, uint8_t* k,
                                                    const faest_paramset_t* params) {
  // for scan-build
  assert(FAEST_192F_Ske == params->faest_param.Ske);

  // Step: 2
  aes_key_schedule_forward_1(w, k, params);

  // Step: 3
  // aes_key_schedule_forward_192_vbb(vbb, vk);

  // Step: 4
  uint8_t w_dash[FAEST_192F_Ske];
  aes_key_schedule_backward_1(w + FAEST_192F_LAMBDA / 8, k, w_dash, params);

  // Step: 5
  bf192_t v_w_dash[FAEST_192F_Ske * 8];
  aes_key_schedule_backward_192_vbb_vk(vbb, 1, 0, NULL, v_w_dash);

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
      bf_v_k_hat[(r + 3) % 4] = bf192_byte_combine_vk(vbb, (iwd + 8 * r));
      bf_w_dash_hat[r]        = bf192_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
      bf_v_w_dash_hat[r]      = bf192_byte_combine(v_w_dash + (32 * j + 8 * r));
    }
    // Step: 13..17
    for (unsigned int r = 0; r <= 3; r++) {
      // instead of storing in A0, A1, hash it
      const bf192_t tmp = bf192_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
      zk_hash_192_update(a0_ctx, tmp);
      zk_hash_192_update(
          a1_ctx, bf192_add(bf192_add(bf192_mul(bf192_add(bf_k_hat[r], bf_v_k_hat[r]),
                                                bf192_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                      bf192_one()),
                            tmp));
    }
    iwd = iwd + 192;
  }
}

static void aes_key_schedule_constraints_Mkey_1_192(vbb_t* vbb, const uint8_t* delta,
                                                    zk_hash_192_ctx* b0_ctx) {
  // Step: 19..20
  // aes_key_schedule_forward_192_vbb(vbb, qk);
  bf192_t q_w_dash[FAEST_192F_Ske * 8];
  aes_key_schedule_backward_192_vbb_vk(vbb, 0, 1, delta, q_w_dash);

  const bf192_t bf_delta      = bf192_load(delta);
  const bf192_t delta_squared = bf192_mul(bf_delta, bf_delta);

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_192F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_192F_Ske / 4; j++) {
    bf192_t bf_q_hat_k[4];
    bf192_t bf_q_hat_w_dash[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 25..26
      bf_q_hat_k[(r + 3) % 4] = bf192_byte_combine_vk(vbb, ((iwd + 8 * r)));
      bf_q_hat_w_dash[r]      = bf192_byte_combine(q_w_dash + ((32 * j + 8 * r)));
    }
    // Step: 27
    for (unsigned int r = 0; r <= 3; r++) {
      bf192_t bf_tmp = bf192_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      // instead of storing B, hash it
      zk_hash_192_update(b0_ctx, bf192_add(bf_tmp, delta_squared));
    }
    iwd = iwd + 192;
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

static void aes_enc_forward_192_vbb_vk(vbb_t* vbb, unsigned int offset, const uint8_t* in,
                                       uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                       bf192_t* bf_y) {
  const bf192_t bf_delta  = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t bf_factor = bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey));

  // Step: 2..4
  for (unsigned int i = 0; i < 16; i++) {
    bf192_t bf_xin[8];
    for (unsigned int j = 0; j < 8; j++) {
      bf_xin[j] = bf192_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
    }
    // Step: 5
    bf_y[i] = bf192_add(bf192_byte_combine(bf_xin), bf192_byte_combine_vk(vbb, (8 * i)));
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
        bf_x_hat[r]  = bf192_byte_combine_vbb(vbb, offset + (ix + 8 * r));
        bf_xk_hat[r] = bf192_byte_combine_vk(vbb, (ik + 8 * r));
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

static void aes_enc_backward_192_vbb_linear_access(vbb_t* vbb, unsigned int offset, uint8_t Mtag,
                                                   uint8_t Mkey, const uint8_t* delta,
                                                   const uint8_t* out, bf192_t* y_out) {
  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t factor =
      bf192_mul_bit(bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey)), 1 ^ Mtag);
  // Access pattern friendly for OLEs
  bf192_t bf_x_tilde[8];
  for (unsigned int i = 0; i < (FAEST_192F_R - 1) * 8 * 16; i++) {
    unsigned int j = i / 128;
    unsigned int r = ((i - (128 * j)) % 32) / 8;
    unsigned int c = ((i - 128 * j - 8 * r) / 32 + r) % 4;

    memcpy(bf_x_tilde + (i % 8), get_vole_aes_192(vbb, i + offset), sizeof(bf192_t));

    if (i % 8 == 7) {
      bf192_t bf_y_tilde[8];
      for (unsigned int k = 0; k < 8; ++k) {
        bf_y_tilde[k] = bf192_add(bf192_add(bf_x_tilde[(k + 7) % 8], bf_x_tilde[(k + 5) % 8]),
                                  bf_x_tilde[(k + 2) % 8]);
      }
      bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
      bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);

      // Step: 18
      y_out[16 * j + 4 * c + r] = bf192_byte_combine(bf_y_tilde);
    }
  }

  unsigned int j = FAEST_192F_R - 1;
  for (unsigned int c = 0; c <= 3; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      memset(bf_x_tilde, 0, sizeof(bf_x_tilde));
      // Step: 5
      unsigned int ird = (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);

      // Step: 10
      for (unsigned int i = 0; i < 8; ++i) {
        // Step: 11
        bf192_t bf_xout =
            bf192_mul_bit(factor, get_bit(out[(ird - 128 * (FAEST_192F_R - 1)) / 8], i));
        // Step: 12
        bf_x_tilde[i] = bf192_add(bf_xout, *get_vk_192(vbb, 128 + ird + i));
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

static void aes_enc_constraints_Mkey_0_192(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                           vbb_t* vbb, unsigned int offset, const uint8_t* k,
                                           zk_hash_192_ctx* a0_ctx, zk_hash_192_ctx* a1_ctx) {
  unsigned int w_offset = offset / 8;
  w += w_offset;
  bf192_t s[FAEST_192F_Senc];
  bf192_t vs[FAEST_192F_Senc];
  bf192_t s_dash[FAEST_192F_Senc];
  bf192_t vs_dash[FAEST_192F_Senc];
  aes_enc_forward_192_1(w, k, in, 0, 0, s);
  aes_enc_forward_192_vbb_vk(vbb, offset, in, 1, 0, NULL, vs);
  aes_enc_backward_192_1(w, k, 0, 0, out, s_dash);
  aes_enc_backward_192_vbb_linear_access(vbb, offset, 1, 0, NULL, out, vs_dash);

  for (unsigned int j = 0; j < FAEST_192F_Senc; j++) {
    // instead of storing in A0, A!, hash it
    const bf192_t tmp = bf192_mul(vs[j], vs_dash[j]);
    zk_hash_192_update(a0_ctx, tmp);
    zk_hash_192_update(a1_ctx, bf192_add(bf192_add(bf192_mul(bf192_add(s[j], vs[j]),
                                                             bf192_add(s_dash[j], vs_dash[j])),
                                                   tmp),
                                         bf192_one()));
  }
}

static void aes_enc_constraints_Mkey_1_192(const uint8_t* in, const uint8_t* out, vbb_t* vbb,
                                           unsigned int offset, const uint8_t* delta,
                                           zk_hash_192_ctx* b0_ctx) {
  // Step: 11..12
  bf192_t qs[FAEST_192F_Senc];
  bf192_t qs_dash[FAEST_192F_Senc];
  aes_enc_forward_192_vbb_vk(vbb, offset, in, 0, 1, delta, qs);
  aes_enc_backward_192_vbb_linear_access(vbb, offset, 0, 1, delta, out, qs_dash);

  // Step: 13..14
  bf192_t minus_part = bf192_mul(bf192_load(delta), bf192_load(delta));
  for (unsigned int j = 0; j < FAEST_192F_Senc; j++) {
    // instead of storing it, hash it
    zk_hash_192_update(b0_ctx, bf192_add(bf192_mul(qs[j], qs_dash[j]), minus_part));
  }
}

static void aes_prove_192(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                          const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
                          const faest_paramset_t* params) {
  // Step: 1..2
  // compute on the fly

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7 + 18
  uint8_t* k = malloc((FAEST_192F_R + 1) * 128 / 8);
  // bf192_t* vk = faest_aligned_alloc(BF192_ALIGN, sizeof(bf192_t) * ((FAEST_192F_R + 1) * 128));
  zk_hash_192_ctx a0_ctx;
  zk_hash_192_ctx a1_ctx;

  zk_hash_192_init(&a0_ctx, chall);
  zk_hash_192_init(&a1_ctx, chall);
  aes_key_schedule_constraints_Mkey_0_192(w, vbb, &a0_ctx, &a1_ctx, k, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  unsigned int offset = FAEST_192F_Lke;
  aes_enc_constraints_Mkey_0_192(in, out, w, vbb, offset, k, &a0_ctx, &a1_ctx);
  // Step: 12-15
  offset = FAEST_192F_Lke + FAEST_192F_Lenc;
  aes_enc_constraints_Mkey_0_192(in + 16, out + 16, w, vbb, offset, k, &a0_ctx, &a1_ctx);
  // faest_aligned_free(vk);
  free(k);

  // Step: 16..18
  zk_hash_192_finalize(a_tilde, &a1_ctx, bf192_load(get_vole_u(vbb) + FAEST_192F_L / 8));
  zk_hash_192_finalize(b_tilde, &a0_ctx, bf192_sum_poly_vbb(vbb, FAEST_192F_L));
}

static uint8_t* aes_verify_192(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                               const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out) {
  // Step: 1
  const uint8_t* delta = chall_3;
  // Step: 2,3
  // do nothing

  // Step: 13 + 21
  // bf192_t* qk = faest_aligned_alloc(BF192_ALIGN, sizeof(bf192_t) * ((FAEST_192F_R + 1) * 128));
  // instead of storing B_0 in an array, we process the values with zk_hash_128
  zk_hash_192_ctx b0_ctx;
  zk_hash_192_init(&b0_ctx, chall_2);
  aes_key_schedule_constraints_Mkey_1_192(vbb, delta, &b0_ctx);

  // Step: 14
  aes_enc_constraints_Mkey_1_192(in, out, vbb, FAEST_192F_Lke, delta, &b0_ctx);

  // Step: 18
  aes_enc_constraints_Mkey_1_192(in + 16, out + 16, vbb, FAEST_192F_Lke + FAEST_192F_Lenc, delta,
                                 &b0_ctx);
  // faest_aligned_free(qk);

  // Step: 20+21
  uint8_t* q_tilde = malloc(FAEST_192F_LAMBDA / 8);
  zk_hash_192_finalize(q_tilde, &b0_ctx, bf192_sum_poly_vbb(vbb, FAEST_192F_L));

  bf192_t bf_qtilde = bf192_load(q_tilde);
  bf192_store(q_tilde, bf192_add(bf_qtilde, bf192_mul(bf192_load(a_tilde), bf192_load(delta))));

  return q_tilde;
}

// lambda == 256 implementation

static void aes_key_schedule_backward_256_vbb_vk(vbb_t* vbb, uint8_t Mtag, uint8_t Mkey,
                                                 const uint8_t* delta, bf256_t* bf_out) {
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
      bf_x_tilde[i] = bf256_add(*get_vole_aes_256(vbb, (8 * j + i) + FAEST_256F_LAMBDA),
                                *get_vk_256(vbb, iwd + 8 * c + i)); // Vk[iwd + 8 * c + i]);
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

static void aes_key_schedule_constraints_Mkey_0_256(const uint8_t* w, vbb_t* vbb,
                                                    zk_hash_256_ctx* a0_ctx,
                                                    zk_hash_256_ctx* a1_ctx, uint8_t* k,
                                                    const faest_paramset_t* params) {
  // for scan-build
  assert(FAEST_256F_Ske == params->faest_param.Ske);

  // Step: 2
  aes_key_schedule_forward_1(w, k, params);

  // Step: 3
  // aes_key_schedule_forward_256_vbb(vbb, vk);

  // Step: 4
  uint8_t w_dash[FAEST_256F_Ske];
  aes_key_schedule_backward_1(w + FAEST_256F_LAMBDA / 8, k, w_dash, params);

  // Step: 5
  bf256_t v_w_dash[FAEST_256F_Ske * 8];
  aes_key_schedule_backward_256_vbb_vk(vbb, 1, 0, NULL, v_w_dash);

  // Step: 6..8
  bool rotate_word = true;
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
        bf_v_k_hat[(r + 3) % 4] = bf256_byte_combine_vk(vbb, (iwd + 8 * r));
        bf_w_dash_hat[r]        = bf256_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
        bf_v_w_dash_hat[r]      = bf256_byte_combine(v_w_dash + (32 * j + 8 * r));
      } else {
        bf_k_hat[r]        = bf256_byte_combine_bits(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat[r]      = bf256_byte_combine_vk(vbb, (iwd + 8 * r));
        bf_w_dash_hat[r]   = bf256_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
        bf_v_w_dash_hat[r] = bf256_byte_combine(v_w_dash + (32 * j + 8 * r));
      }
    }
    // Step: 13..17
    for (unsigned int r = 0; r <= 3; r++) {
      const bf256_t tmp = bf256_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
      zk_hash_256_update(a0_ctx, tmp);
      zk_hash_256_update(
          a1_ctx, bf256_add(bf256_add(bf256_mul(bf256_add(bf_k_hat[r], bf_v_k_hat[r]),
                                                bf256_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                      bf256_one()),
                            tmp));
    }
    iwd         = iwd + 128;
    rotate_word = !rotate_word;
  }
}

static void aes_key_schedule_constraints_Mkey_1_256(vbb_t* vbb, const uint8_t* delta,
                                                    zk_hash_256_ctx* b0_ctx) {
  // Step: 19..20
  // aes_key_schedule_forward_256_vbb(vbb, qk);
  bf256_t q_w_dash[FAEST_256F_Ske * 8];
  aes_key_schedule_backward_256_vbb_vk(vbb, 0, 1, delta, q_w_dash);

  const bf256_t bf_delta         = bf256_load(delta);
  const bf256_t bf_delta_squared = bf256_mul(bf_delta, bf_delta);
  bool rotate_word               = true;

  // Step 23..24
  unsigned int iwd = 32 * (FAEST_256F_Nwd - 1);
  for (unsigned int j = 0; j < FAEST_256F_Ske / 4; j++) {
    bf256_t bf_q_hat_k[4];
    bf256_t bf_q_hat_w_dash[4];
    for (unsigned int r = 0; r <= 3; r++) {
      // Step: 25..26
      if (rotate_word) {
        bf_q_hat_k[(r + 3) % 4] = bf256_byte_combine_vk(vbb, ((iwd + 8 * r)));
        bf_q_hat_w_dash[r]      = bf256_byte_combine(q_w_dash + ((32 * j + 8 * r)));
      } else {
        bf_q_hat_k[r]      = bf256_byte_combine_vk(vbb, ((iwd + 8 * r)));
        bf_q_hat_w_dash[r] = bf256_byte_combine(q_w_dash + ((32 * j + 8 * r)));
      }
    }
    // Step: 27
    for (unsigned int r = 0; r <= 3; r++) {
      bf256_t bf_tmp = bf256_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      zk_hash_256_update(b0_ctx, bf256_add(bf_tmp, bf_delta_squared));
    }
    iwd         = iwd + 128;
    rotate_word = !rotate_word;
  }
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

static void aes_enc_forward_256_vbb_vk(vbb_t* vbb, unsigned int offset, const uint8_t* in,
                                       uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                       bf256_t* bf_y) {
  const bf256_t bf_delta  = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t bf_factor = bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey));

  // Step: 2..4
  for (unsigned int i = 0; i < 16; i++) {
    bf256_t bf_xin[8];
    for (unsigned int j = 0; j < 8; j++) {
      bf_xin[j] = bf256_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
    }
    // Step: 5
    bf_y[i] = bf256_add(bf256_byte_combine(bf_xin), bf256_byte_combine_vk(vbb, (8 * i)));
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
        bf_x_hat[r]  = bf256_byte_combine_vbb(vbb, offset + ix + 8 * r);
        bf_xk_hat[r] = bf256_byte_combine_vk(vbb, (ik + 8 * r));
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

static void aes_enc_backward_256_vbb_linear_access(vbb_t* vbb, unsigned int offset, uint8_t Mtag,
                                                   uint8_t Mkey, const uint8_t* delta,
                                                   const uint8_t* out, bf256_t* y_out) {
  const bf256_t bf_delta = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t factor =
      bf256_mul_bit(bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey)), 1 ^ Mtag);
  // Access pattern friendly for OLEs
  bf256_t bf_x_tilde[8];
  for (unsigned int ird = 0; ird < FAEST_256F_R*128; ird++) {
    unsigned int j = ird / 128;
    unsigned int r = ((ird - (128 * j)) % 32) / 8;
    unsigned int c = ((ird - 128 * j - 8 * r) / 32 + r) % 4;

    if (j != FAEST_256F_R - 1) {
      memcpy(bf_x_tilde + (ird % 8), get_vole_aes_256(vbb, ird + offset), sizeof(bf256_t));
    } else {
      memset(bf_x_tilde, 0, sizeof(bf256_t));
      unsigned int z = (ird/8)*8;

      for (unsigned int i = 0; i < 8; ++i) {
        bf256_t bf_xout =
            bf256_mul_bit(factor, get_bit(out[(z - 128*(FAEST_256F_R-1)) / 8], i));
        bf_x_tilde[i] = bf256_add(bf_xout, *get_vk_256(vbb, 128+z + i));
      }
    }

    if (ird % 8 == 7) {
        bf256_t bf_y_tilde[8];
        for (unsigned int i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf256_add(bf256_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                    bf_x_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);

        y_out[16 * j + 4 * c + r] = bf256_byte_combine(bf_y_tilde);
    }
  }
}

static void aes_enc_constraints_Mkey_0_256(const uint8_t* in, const uint8_t* out, const uint8_t* w,
                                           vbb_t* vbb, unsigned int offset, const uint8_t* k,
                                           zk_hash_256_ctx* a0_ctx, zk_hash_256_ctx* a1_ctx) {
  unsigned int w_offset = offset / 8;
  w += w_offset;
  bf256_t s[FAEST_256F_Senc];
  bf256_t vs[FAEST_256F_Senc];
  bf256_t s_dash[FAEST_256F_Senc];
  bf256_t vs_dash[FAEST_256F_Senc];
  aes_enc_forward_256_1(w, k, in, 0, 0, s);
  aes_enc_forward_256_vbb_vk(vbb, offset, in, 1, 0, NULL, vs);
  aes_enc_backward_256_1(w, k, 0, 0, out, s_dash);
  aes_enc_backward_256_vbb_linear_access(vbb, offset, 1, 0, NULL, out, vs_dash);

  for (unsigned int j = 0; j < FAEST_256F_Senc; j++) {
    const bf256_t tmp = bf256_mul(vs[j], vs_dash[j]);
    zk_hash_256_update(a0_ctx, tmp);
    zk_hash_256_update(a1_ctx, bf256_add(bf256_add(bf256_mul(bf256_add(s[j], vs[j]),
                                                             bf256_add(s_dash[j], vs_dash[j])),
                                                   tmp),
                                         bf256_one()));
  }
}

static void aes_enc_constraints_Mkey_1_256(const uint8_t* in, const uint8_t* out, vbb_t* vbb,
                                           unsigned int offset, const uint8_t* delta,
                                           zk_hash_256_ctx* b0_ctx) {
  // Step: 11..12
  bf256_t qs[FAEST_256F_Senc];
  bf256_t qs_dash[FAEST_256F_Senc];
  aes_enc_forward_256_vbb_vk(vbb, offset, in, 0, 1, delta, qs);
  aes_enc_backward_256_vbb_linear_access(vbb, offset, 0, 1, delta, out, qs_dash);

  // Step: 13..14
  bf256_t minus_part = bf256_mul(bf256_load(delta), bf256_load(delta));
  for (unsigned int j = 0; j < FAEST_256F_Senc; j++) {
    zk_hash_256_update(b0_ctx, bf256_add(bf256_mul(qs[j], qs_dash[j]), minus_part));
  }
}

static void aes_prove_256(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                          const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
                          const faest_paramset_t* params) {
  // Step: 1..2
  // compute on the fly

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7 + 18
  uint8_t* k = malloc((FAEST_256F_R + 1) * 128 / 8);
  // bf256_t* vk = faest_aligned_alloc(BF256_ALIGN, sizeof(bf256_t) * ((FAEST_256F_R + 1) * 128));
  zk_hash_256_ctx a0_ctx;
  zk_hash_256_ctx a1_ctx;

  zk_hash_256_init(&a0_ctx, chall);
  zk_hash_256_init(&a1_ctx, chall);
  aes_key_schedule_constraints_Mkey_0_256(w, vbb, &a0_ctx, &a1_ctx, k, params);

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  unsigned int offset = FAEST_256F_Lke;
  aes_enc_constraints_Mkey_0_256(in, out, w, vbb, offset, k, &a0_ctx, &a1_ctx);
  // Step: 12-15
  offset = FAEST_256F_Lke + FAEST_256F_Lenc;
  aes_enc_constraints_Mkey_0_256(in + 16, out + 16, w, vbb, offset, k, &a0_ctx, &a1_ctx);
  // faest_aligned_free(vk);
  free(k);

  // Step: 16..18
  zk_hash_256_finalize(a_tilde, &a1_ctx, bf256_load(get_vole_u(vbb) + FAEST_256F_L / 8));
  zk_hash_256_finalize(b_tilde, &a0_ctx, bf256_sum_poly_vbb(vbb, FAEST_256F_L, 0));
}

static uint8_t* aes_verify_256(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                               const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out) {
  // Step: 1
  const uint8_t* delta = chall_3;
  // Step: 2,3
  // do nothing

  // Step: 13, 21
  // bf256_t* qk = faest_aligned_alloc(BF256_ALIGN, sizeof(bf256_t) * ((FAEST_256F_R + 1) * 128));
  zk_hash_256_ctx b0_ctx;

  zk_hash_256_init(&b0_ctx, chall_2);
  aes_key_schedule_constraints_Mkey_1_256(vbb, delta, &b0_ctx);

  // Step: 14
  aes_enc_constraints_Mkey_1_256(in, out, vbb, FAEST_256F_Lke, delta, &b0_ctx);
  // Step: 18
  aes_enc_constraints_Mkey_1_256(in + 16, out + 16, vbb, FAEST_256F_Lke + FAEST_256F_Lenc, delta,
                                 &b0_ctx);
  // faest_aligned_free(qk);

  // Step: 20, 21
  uint8_t* q_tilde = malloc(FAEST_256F_LAMBDA / 8);
  zk_hash_256_finalize(q_tilde, &b0_ctx, bf256_sum_poly_vbb(vbb, FAEST_256F_L, 0));

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

static void em_enc_forward_128_vbb(vbb_t* vbb, const bf128_t* bf_x, bf128_t* bf_y) {
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_128F_Nwd; j++) {
    bf_y[j] = bf128_byte_combine_vbb(vbb, 8 * j);
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
        bf_z_hat[r] = bf128_byte_combine_vbb(vbb, (i + 8 * r));
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

static void em_enc_backward_128_linear_access_verify(vbb_t* vbb, const bf128_t* bf_x,
                                                     const bf128_t* bf_z_out, uint8_t Mtag,
                                                     uint8_t Mkey, const uint8_t* delta,
                                                     bf128_t* y_out) {
  // Step: 1
  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t factor =
      bf128_mul_bit(bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int c = 0; c < FAEST_EM_128F_Nwd; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      unsigned int j          = FAEST_EM_128F_R - 1;
      const unsigned int icol = (c - r + FAEST_EM_128F_Nwd) % FAEST_EM_128F_Nwd;
      const unsigned int ird =
          FAEST_EM_128F_LAMBDA + 32 * FAEST_EM_128F_Nwd * j + 32 * icol + 8 * r;

      bf128_t bf_z_tilde[8];
      for (unsigned int i = 0; i < 8; ++i) {
        // Step: 12
        bf_z_tilde[i] = bf_z_out[ird - 32 * FAEST_EM_128F_Nwd * (j + 1) + i];
        if (bf_x) {
          bf_z_tilde[i] = bf128_add(bf_z_tilde[i], bf_x[ird + i]);
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

  for (unsigned int chunk_idx = FAEST_EM_128F_LAMBDA;
       chunk_idx < FAEST_EM_128F_R * FAEST_EM_128F_Nwd * 4 * 8; chunk_idx += 8) {
    unsigned int ird = chunk_idx;
    unsigned int j   = (ird - FAEST_EM_128F_LAMBDA) / (32 * FAEST_EM_128F_Nwd);
    unsigned int r   = (ird - FAEST_EM_128F_LAMBDA - j * 32 * FAEST_EM_128F_Nwd) % 32 / 8;
    unsigned int c = ((ird - FAEST_EM_128F_LAMBDA - j * 32 * FAEST_EM_128F_Nwd - 8 * r) / 32 + r) %
                     FAEST_EM_128F_Nwd;

    bf128_t bf_z_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      memcpy(bf_z_tilde + i, get_vole_aes_128(vbb, chunk_idx + i), sizeof(bf128_t));
    }

    bf128_t bf_y_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf128_add(bf128_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                bf_z_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
    bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

    y_out[4 * FAEST_EM_128F_Nwd * j + 4 * c + r] = bf128_byte_combine(bf_y_tilde);
  }
}

static void em_enc_backward_128_vbb_linear_access(vbb_t* vbb, const bf128_t* bf_x, vbb_t* vbb_out,
                                                  uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                                  bf128_t* y_out) {
  // Step: 1
  const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
  const bf128_t factor =
      bf128_mul_bit(bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int chunk_idx = 0; chunk_idx < FAEST_EM_128F_LAMBDA; chunk_idx += 8) {
    unsigned int j   = FAEST_EM_128F_R - 1;
    unsigned int ird = chunk_idx + 32 * FAEST_EM_128F_Nwd * (j + 1);
    unsigned int r   = chunk_idx % 32 / 8;
    unsigned int c   = ((chunk_idx - 8 * r) / 32 + r) % FAEST_EM_128F_Nwd;

    bf128_t bf_z_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_z_tilde[i] = *get_vole_aes_128(vbb_out, chunk_idx + i);
      if (bf_x) {
        bf_z_tilde[i] = bf128_add(bf_z_tilde[i], bf_x[ird + i]);
      }
    }

    bf128_t bf_y_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf128_add(bf128_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                bf_z_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
    bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

    y_out[4 * FAEST_EM_128F_Nwd * j + 4 * c + r] = bf128_byte_combine(bf_y_tilde);
  }

  for (unsigned int chunk_idx = FAEST_EM_128F_LAMBDA;
       chunk_idx < FAEST_EM_128F_R * FAEST_EM_128F_Nwd * 4 * 8; chunk_idx += 8) {
    unsigned int ird = chunk_idx;
    unsigned int j   = (ird - FAEST_EM_128F_LAMBDA) / (32 * FAEST_EM_128F_Nwd);
    unsigned int r   = (ird - FAEST_EM_128F_LAMBDA - j * 32 * FAEST_EM_128F_Nwd) % 32 / 8;
    unsigned int c = ((ird - FAEST_EM_128F_LAMBDA - j * 32 * FAEST_EM_128F_Nwd - 8 * r) / 32 + r) %
                     FAEST_EM_128F_Nwd;

    bf128_t bf_z_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      memcpy(bf_z_tilde + i, get_vole_aes_128(vbb, chunk_idx + i), sizeof(bf128_t));
    }

    bf128_t bf_y_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf128_add(bf128_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                bf_z_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
    bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

    y_out[4 * FAEST_EM_128F_Nwd * j + 4 * c + r] = bf128_byte_combine(bf_y_tilde);
  }
}

static void em_enc_constraints_Mkey_0_128(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                          vbb_t* vbb, zk_hash_128_ctx* a0_ctx,
                                          zk_hash_128_ctx* a1_ctx) {
  // Step 6
  uint8_t w_out[FAEST_EM_128F_LAMBDA / 8];
  xor_u8_array(out, w, w_out, sizeof(w_out));

  bf128_t bf_s[FAEST_EM_128F_Senc];
  bf128_t bf_vs[FAEST_EM_128F_Senc];
  bf128_t bf_s_dash[FAEST_EM_128F_Senc];
  bf128_t bf_vs_dash[FAEST_EM_128F_Senc];
  em_enc_forward_128_1(w, x, bf_s);
  em_enc_forward_128_vbb(vbb, NULL, bf_vs);
  em_enc_backward_128_1(w, x, w_out, bf_s_dash);
  em_enc_backward_128_vbb_linear_access(vbb, NULL, vbb, 1, 0, NULL, bf_vs_dash);

  for (unsigned int j = 0; j < FAEST_EM_128F_Senc; j++) {
    const bf128_t tmp = bf128_mul(bf_vs[j], bf_vs_dash[j]);
    zk_hash_128_update(a0_ctx, tmp);
    zk_hash_128_update(a1_ctx,
                       bf128_add(bf128_add(bf128_mul(bf128_add(bf_s[j], bf_vs[j]),
                                                     bf128_add(bf_s_dash[j], bf_vs_dash[j])),
                                           tmp),
                                 bf128_one()));
  }
}

static void em_enc_constraints_Mkey_1_128(const uint8_t* out, const uint8_t* x, vbb_t* vbb,
                                          const uint8_t* delta, zk_hash_128_ctx* b0_ctx) {
  // Step: 18, 19
  // TODO: compute these on demand in em_enc_backward_128
  const bf128_t bf_delta = bf128_load(delta);
  bf128_t* bf_x = faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * 128 * (FAEST_EM_128F_R + 1));
  for (unsigned int i = 0; i < 128 * (FAEST_EM_128F_R + 1); i++) {
    bf_x[i] = bf128_mul_bit(bf_delta, ptr_get_bit(x, i));
  }

  // Step 21
  bf128_t* bf_q_out = faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * FAEST_EM_128F_LAMBDA);
  for (unsigned int i = 0; i < FAEST_EM_128F_LAMBDA; i++) {
    bf_q_out[i] =
        bf128_add(bf128_mul_bit(bf_delta, ptr_get_bit(out, i)), *get_vole_aes_128(vbb, i));
  }

  bf128_t bf_qs[FAEST_EM_128F_Senc];
  bf128_t bf_qs_dash[FAEST_EM_128F_Senc];
  em_enc_forward_128_vbb(vbb, bf_x, bf_qs);
  em_enc_backward_128_linear_access_verify(vbb, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash);
  faest_aligned_free(bf_q_out);
  faest_aligned_free(bf_x);

  // Step: 13..14
  bf128_t minus_part = bf128_mul(bf_delta, bf_delta);
  for (unsigned int j = 0; j < FAEST_EM_128F_Senc; j++) {
    zk_hash_128_update(b0_ctx, bf128_add(bf128_mul(bf_qs[j], bf_qs_dash[j]), minus_part));
  }
}

static void em_prove_128(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                         const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde) {
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

  zk_hash_128_ctx a0_ctx;
  zk_hash_128_ctx a1_ctx;

  zk_hash_128_init(&a0_ctx, chall);
  zk_hash_128_init(&a1_ctx, chall);
  em_enc_constraints_Mkey_0_128(out, x, w, vbb, &a0_ctx, &a1_ctx);

  zk_hash_128_finalize(a_tilde, &a1_ctx, bf128_load(get_vole_u(vbb) + FAEST_EM_128F_Lenc / 8));
  zk_hash_128_finalize(b_tilde, &a0_ctx, bf128_sum_poly_vbb(vbb, FAEST_EM_128F_Lenc));
}

static uint8_t* em_verify_128(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                              const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out) {
  const uint8_t* delta = chall_3;

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
  em_enc_constraints_Mkey_1_128(out, x, vbb, delta, &b0_ctx);

  uint8_t* q_tilde = malloc(FAEST_EM_128F_LAMBDA / 8);
  zk_hash_128_finalize(q_tilde, &b0_ctx, bf128_sum_poly_vbb(vbb, FAEST_EM_128F_Lenc));

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

static void em_enc_forward_192_vbb(vbb_t* vbb, const bf192_t* bf_x, bf192_t* bf_y) {
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_192F_Nwd; j++) {
    bf_y[j] = bf192_byte_combine_vbb(vbb, 8 * j);
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
        bf_z_hat[r] = bf192_byte_combine_vbb(vbb, i + 8 * r);
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

static void em_enc_backward_192_linear_access_verify(vbb_t* vbb, const bf192_t* bf_x,
                                                     const bf192_t* bf_z_out, uint8_t Mtag,
                                                     uint8_t Mkey, const uint8_t* delta,
                                                     bf192_t* y_out) {
  // Step: 1
  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t factor =
      bf192_mul_bit(bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int c = 0; c < FAEST_EM_192F_Nwd; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      unsigned int j          = FAEST_EM_192F_R - 1;
      const unsigned int icol = (c - r + FAEST_EM_192F_Nwd) % FAEST_EM_192F_Nwd;
      const unsigned int ird =
          FAEST_EM_192F_LAMBDA + 32 * FAEST_EM_192F_Nwd * j + 32 * icol + 8 * r;

      bf192_t bf_z_tilde[8];
      for (unsigned int i = 0; i < 8; ++i) {
        // Step: 12
        bf_z_tilde[i] = bf_z_out[ird - 32 * FAEST_EM_192F_Nwd * (j + 1) + i];
        if (bf_x) {
          bf_z_tilde[i] = bf192_add(bf_z_tilde[i], bf_x[ird + i]);
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

  for (unsigned int chunk_idx = FAEST_EM_192F_LAMBDA;
       chunk_idx < FAEST_EM_192F_R * FAEST_EM_192F_Nwd * 4 * 8; chunk_idx += 8) {
    unsigned int ird = chunk_idx;
    unsigned int j   = (ird - FAEST_EM_192F_LAMBDA) / (32 * FAEST_EM_192F_Nwd);
    unsigned int r   = (ird - FAEST_EM_192F_LAMBDA - j * 32 * FAEST_EM_192F_Nwd) % 32 / 8;
    unsigned int c = ((ird - FAEST_EM_192F_LAMBDA - j * 32 * FAEST_EM_192F_Nwd - 8 * r) / 32 + r) %
                     FAEST_EM_192F_Nwd;

    bf192_t bf_z_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      memcpy(bf_z_tilde + i, get_vole_aes_192(vbb, chunk_idx + i), sizeof(bf192_t));
    }

    bf192_t bf_y_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf192_add(bf192_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                bf_z_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
    bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);

    y_out[4 * FAEST_EM_192F_Nwd * j + 4 * c + r] = bf192_byte_combine(bf_y_tilde);
  }
}

static void em_enc_backward_192_vbb_linear_access(vbb_t* vbb, const bf192_t* bf_x, vbb_t* vbb_out,
                                                  uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                                  bf192_t* y_out) {
  // Step: 1
  const bf192_t bf_delta = delta ? bf192_load(delta) : bf192_zero();
  const bf192_t factor =
      bf192_mul_bit(bf192_add(bf192_mul_bit(bf_delta, Mkey), bf192_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int chunk_idx = 0; chunk_idx < FAEST_EM_192F_LAMBDA; chunk_idx += 8) {
    unsigned int j   = FAEST_EM_192F_R - 1;
    unsigned int ird = chunk_idx + 32 * FAEST_EM_192F_Nwd * (j + 1);
    unsigned int r   = chunk_idx % 32 / 8;
    unsigned int c   = ((chunk_idx - 8 * r) / 32 + r) % FAEST_EM_192F_Nwd;

    bf192_t bf_z_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_z_tilde[i] = *get_vole_aes_192(vbb_out, chunk_idx + i);
      if (bf_x) {
        bf_z_tilde[i] = bf192_add(bf_z_tilde[i], bf_x[ird + i]);
      }
    }

    bf192_t bf_y_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf192_add(bf192_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                bf_z_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
    bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);

    y_out[4 * FAEST_EM_192F_Nwd * j + 4 * c + r] = bf192_byte_combine(bf_y_tilde);
  }

  for (unsigned int chunk_idx = FAEST_EM_192F_LAMBDA;
       chunk_idx < FAEST_EM_192F_R * FAEST_EM_192F_Nwd * 4 * 8; chunk_idx += 8) {
    unsigned int ird = chunk_idx;
    unsigned int j   = (ird - FAEST_EM_192F_LAMBDA) / (32 * FAEST_EM_192F_Nwd);
    unsigned int r   = (ird - FAEST_EM_192F_LAMBDA - j * 32 * FAEST_EM_192F_Nwd) % 32 / 8;
    unsigned int c = ((ird - FAEST_EM_192F_LAMBDA - j * 32 * FAEST_EM_192F_Nwd - 8 * r) / 32 + r) %
                     FAEST_EM_192F_Nwd;

    bf192_t bf_z_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      memcpy(bf_z_tilde + i, get_vole_aes_192(vbb, chunk_idx + i), sizeof(bf192_t));
    }

    bf192_t bf_y_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf192_add(bf192_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                bf_z_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf192_add(bf_y_tilde[0], factor);
    bf_y_tilde[2] = bf192_add(bf_y_tilde[2], factor);

    y_out[4 * FAEST_EM_192F_Nwd * j + 4 * c + r] = bf192_byte_combine(bf_y_tilde);
  }
}

static void em_enc_constraints_Mkey_0_192(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                          vbb_t* vbb, zk_hash_192_ctx* a0_ctx,
                                          zk_hash_192_ctx* a1_ctx) {
  // Step 6
  uint8_t w_out[FAEST_EM_192F_LAMBDA / 8];
  xor_u8_array(out, w, w_out, sizeof(w_out));

  bf192_t bf_s[FAEST_EM_192F_Senc];
  bf192_t bf_vs[FAEST_EM_192F_Senc];
  bf192_t bf_s_dash[FAEST_EM_192F_Senc];
  bf192_t bf_vs_dash[FAEST_EM_192F_Senc];
  em_enc_forward_192_1(w, x, bf_s);
  em_enc_forward_192_vbb(vbb, NULL, bf_vs);
  em_enc_backward_192_1(w, x, w_out, bf_s_dash);
  em_enc_backward_192_vbb_linear_access(vbb, NULL, vbb, 1, 0, NULL, bf_vs_dash);

  for (unsigned int j = 0; j < FAEST_EM_192F_Senc; j++) {
    const bf192_t tmp = bf192_mul(bf_vs[j], bf_vs_dash[j]);
    zk_hash_192_update(a0_ctx, tmp);
    zk_hash_192_update(a1_ctx,
                       bf192_add(bf192_add(bf192_mul(bf192_add(bf_s[j], bf_vs[j]),
                                                     bf192_add(bf_s_dash[j], bf_vs_dash[j])),
                                           tmp),
                                 bf192_one()));
  }
}

static void em_enc_constraints_Mkey_1_192(const uint8_t* out, const uint8_t* x, vbb_t* vbb,
                                          const uint8_t* delta, zk_hash_192_ctx* b0_ctx) {
  // Step: 18, 19
  // TODO: compute these on demand in em_enc_backward_192
  const bf192_t bf_delta = bf192_load(delta);
  bf192_t* bf_x = faest_aligned_alloc(BF192_ALIGN, sizeof(bf192_t) * 192 * (FAEST_EM_192F_R + 1));
  for (unsigned int i = 0; i < 192 * (FAEST_EM_192F_R + 1); i++) {
    bf_x[i] = bf192_mul_bit(bf_delta, ptr_get_bit(x, i));
  }

  // Step 21
  bf192_t* bf_q_out = faest_aligned_alloc(BF192_ALIGN, sizeof(bf192_t) * FAEST_EM_192F_LAMBDA);
  for (unsigned int i = 0; i < FAEST_EM_192F_LAMBDA; i++) {
    bf_q_out[i] =
        bf192_add(bf192_mul_bit(bf_delta, ptr_get_bit(out, i)), *get_vole_aes_192(vbb, i));
  }

  bf192_t bf_qs[FAEST_EM_192F_Senc];
  bf192_t bf_qs_dash[FAEST_EM_192F_Senc];
  em_enc_forward_192_vbb(vbb, bf_x, bf_qs);
  em_enc_backward_192_linear_access_verify(vbb, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash);
  faest_aligned_free(bf_q_out);
  faest_aligned_free(bf_x);

  // Step: 13..14
  bf192_t minus_part = bf192_mul(bf_delta, bf_delta);
  for (unsigned int j = 0; j < FAEST_EM_192F_Senc; j++) {
    zk_hash_192_update(b0_ctx, bf192_add(bf192_mul(bf_qs[j], bf_qs_dash[j]), minus_part));
  }
}

static void em_prove_192(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                         const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde) {
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

  zk_hash_192_ctx a0_ctx;
  zk_hash_192_ctx a1_ctx;

  zk_hash_192_init(&a0_ctx, chall);
  zk_hash_192_init(&a1_ctx, chall);
  em_enc_constraints_Mkey_0_192(out, x, w, vbb, &a0_ctx, &a1_ctx);

  zk_hash_192_finalize(a_tilde, &a1_ctx, bf192_load(get_vole_u(vbb) + FAEST_EM_192F_Lenc / 8));
  zk_hash_192_finalize(b_tilde, &a0_ctx, bf192_sum_poly_vbb(vbb, FAEST_EM_192F_Lenc));
}

static uint8_t* em_verify_192(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                              const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out) {
  const uint8_t* delta = chall_3;

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
  em_enc_constraints_Mkey_1_192(out, x, vbb, delta, &b0_ctx);

  uint8_t* q_tilde = malloc(FAEST_EM_192F_LAMBDA / 8);
  zk_hash_192_finalize(q_tilde, &b0_ctx, bf192_sum_poly_vbb(vbb, FAEST_EM_192F_Lenc));

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

static void em_enc_forward_256_vbb(vbb_t* vbb, const bf256_t* bf_x, bf256_t* bf_y) {
  // Step: 2
  for (unsigned int j = 0; j < 4 * FAEST_EM_256F_Nwd; j++) {
    bf_y[j] = bf256_byte_combine_vbb(vbb, 8 * j);
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
        bf_z_hat[r] = bf256_byte_combine_vbb(vbb, i + 8 * r);
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

static void em_enc_backward_256_linear_access_verify(vbb_t* vbb, const bf256_t* bf_x, const bf256_t* bf_z_out,
                                       uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                       bf256_t* y_out) {
  // Step: 1
  const bf256_t bf_delta = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t factor =
      bf256_mul_bit(bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {
    for (unsigned int r = 0; r <= 3; r++) {
      unsigned int j    = FAEST_EM_256F_R - 1;
      unsigned int icol = (c - r + FAEST_EM_256F_Nwd) % FAEST_EM_256F_Nwd;
      if (FAEST_EM_256F_Nwd == 8 && r >= 2) {
        icol = (icol - 1 + FAEST_EM_256F_Nwd) % FAEST_EM_256F_Nwd;
      }
      const unsigned int ird =
          FAEST_EM_256F_LAMBDA + 32 * FAEST_EM_256F_Nwd * j + 32 * icol + 8 * r;

      bf256_t bf_z_tilde[8];
      for (unsigned int i = 0; i < 8; ++i) {
        // Step: 12
        bf_z_tilde[i] = bf_z_out[ird - 32 * FAEST_EM_256F_Nwd * (j + 1) + i];
        if (bf_x) {
          bf_z_tilde[i] = bf256_add(bf_z_tilde[i], bf_x[ird + i]);
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

  for (unsigned int chunk_idx = FAEST_EM_256F_LAMBDA;
       chunk_idx < FAEST_EM_256F_R * FAEST_EM_256F_Nwd * 4 * 8; chunk_idx += 8) {
    unsigned int ird = chunk_idx;
    unsigned int j   = (ird - FAEST_EM_256F_LAMBDA) / (32 * FAEST_EM_256F_Nwd);
    unsigned int r   = (ird - FAEST_EM_256F_LAMBDA - j * 32 * FAEST_EM_256F_Nwd) % 32 / 8;
    unsigned int c = ((ird - FAEST_EM_256F_LAMBDA - j * 32 * FAEST_EM_256F_Nwd - 8 * r) / 32 + r) %
                     FAEST_EM_256F_Nwd;

    if (FAEST_EM_256F_Nwd == 8 && r >= 2) {
      c = (c + 1) % FAEST_EM_256F_Nwd;
    }

    bf256_t bf_z_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      memcpy(bf_z_tilde + i, get_vole_aes_256(vbb, chunk_idx + i), sizeof(bf256_t));
    }

    bf256_t bf_y_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf256_add(bf256_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                bf_z_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
    bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);

    y_out[4 * FAEST_EM_256F_Nwd * j + 4 * c + r] = bf256_byte_combine(bf_y_tilde);
  }
}

static void em_enc_backward_256_vbb_linear_access(vbb_t* vbb, const bf256_t* bf_x, vbb_t* vbb_out,
                                                  uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                                  bf256_t* y_out) {
  // Step: 1
  const bf256_t bf_delta = delta ? bf256_load(delta) : bf256_zero();
  const bf256_t factor =
      bf256_mul_bit(bf256_add(bf256_mul_bit(bf_delta, Mkey), bf256_from_bit(1 ^ Mkey)), 1 ^ Mtag);

  for (unsigned int chunk_idx = 0; chunk_idx < FAEST_EM_256F_LAMBDA; chunk_idx += 8) {
    unsigned int j   = FAEST_EM_256F_R - 1;
    unsigned int ird = chunk_idx + 32 * FAEST_EM_256F_Nwd * (j + 1);
    unsigned int r   = chunk_idx % 32 / 8;
    unsigned int c   = ((chunk_idx - 8 * r) / 32 + r) % FAEST_EM_256F_Nwd;

    if (FAEST_EM_256F_Nwd == 8 && r >= 2) {
      c = (c + 1) % FAEST_EM_256F_Nwd;
    }

    bf256_t bf_z_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_z_tilde[i] = *get_vole_aes_256(vbb_out, chunk_idx + i);
      if (bf_x) {
        bf_z_tilde[i] = bf256_add(bf_z_tilde[i], bf_x[ird + i]);
      }
    }

    bf256_t bf_y_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf256_add(bf256_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                bf_z_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
    bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);

    y_out[4 * FAEST_EM_256F_Nwd * j + 4 * c + r] = bf256_byte_combine(bf_y_tilde);
  }

  for (unsigned int chunk_idx = FAEST_EM_256F_LAMBDA;
       chunk_idx < FAEST_EM_256F_R * FAEST_EM_256F_Nwd * 4 * 8; chunk_idx += 8) {
    unsigned int ird = chunk_idx;
    unsigned int j   = (ird - FAEST_EM_256F_LAMBDA) / (32 * FAEST_EM_256F_Nwd);
    unsigned int r   = (ird - FAEST_EM_256F_LAMBDA - j * 32 * FAEST_EM_256F_Nwd) % 32 / 8;
    unsigned int c = ((ird - FAEST_EM_256F_LAMBDA - j * 32 * FAEST_EM_256F_Nwd - 8 * r) / 32 + r) %
                     FAEST_EM_256F_Nwd;

    if (FAEST_EM_256F_Nwd == 8 && r >= 2) {
      c = (c + 1) % FAEST_EM_256F_Nwd;
    }

    bf256_t bf_z_tilde[8];
    for (unsigned int i = 0; i < 8; i++) {
      memcpy(bf_z_tilde + i, get_vole_aes_256(vbb, chunk_idx + i), sizeof(bf256_t));
    }

    bf256_t bf_y_tilde[8];
    for (unsigned int i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf256_add(bf256_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                bf_z_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf256_add(bf_y_tilde[0], factor);
    bf_y_tilde[2] = bf256_add(bf_y_tilde[2], factor);

    y_out[4 * FAEST_EM_256F_Nwd * j + 4 * c + r] = bf256_byte_combine(bf_y_tilde);
  }
}

static void em_enc_constraints_Mkey_0_256(const uint8_t* out, const uint8_t* x, const uint8_t* w,
                                          vbb_t* vbb, zk_hash_256_ctx* a0_ctx,
                                          zk_hash_256_ctx* a1_ctx) {
  // Step 6
  uint8_t w_out[FAEST_EM_256F_LAMBDA / 8];
  xor_u8_array(out, w, w_out, sizeof(w_out));

  bf256_t bf_s[FAEST_EM_256F_Senc];
  bf256_t bf_vs[FAEST_EM_256F_Senc];
  bf256_t bf_s_dash[FAEST_EM_256F_Senc];
  bf256_t bf_vs_dash[FAEST_EM_256F_Senc];
  em_enc_forward_256_1(w, x, bf_s);
  em_enc_forward_256_vbb(vbb, NULL, bf_vs);
  em_enc_backward_256_1(w, x, w_out, bf_s_dash);
  em_enc_backward_256_vbb_linear_access(vbb, NULL, vbb, 1, 0, NULL, bf_vs_dash);

  for (unsigned int j = 0; j < FAEST_EM_256F_Senc; j++) {
    const bf256_t tmp = bf256_mul(bf_vs[j], bf_vs_dash[j]);
    zk_hash_256_update(a0_ctx, tmp);
    zk_hash_256_update(a1_ctx,
                       bf256_add(bf256_add(bf256_mul(bf256_add(bf_s[j], bf_vs[j]),
                                                     bf256_add(bf_s_dash[j], bf_vs_dash[j])),
                                           tmp),
                                 bf256_one()));
  }
}

static void em_enc_constraints_Mkey_1_256(const uint8_t* out, const uint8_t* x, vbb_t* vbb,
                                          const uint8_t* delta, zk_hash_256_ctx* b0_ctx) {
  // Step: 18, 19
  // TODO: compute these on demand in em_enc_backward_256
  const bf256_t bf_delta = bf256_load(delta);
  bf256_t* bf_x = faest_aligned_alloc(BF256_ALIGN, sizeof(bf256_t) * 256 * (FAEST_EM_256F_R + 1));
  for (unsigned int i = 0; i < 256 * (FAEST_EM_256F_R + 1); i++) {
    bf_x[i] = bf256_mul_bit(bf_delta, ptr_get_bit(x, i));
  }

  // Step 21
  bf256_t* bf_q_out = faest_aligned_alloc(BF256_ALIGN, sizeof(bf256_t) * FAEST_EM_256F_LAMBDA);
  for (unsigned int i = 0; i < FAEST_EM_256F_LAMBDA; i++) {
    bf_q_out[i] =
        bf256_add(bf256_mul_bit(bf_delta, ptr_get_bit(out, i)), *get_vole_aes_256(vbb, i));
  }

  bf256_t bf_qs[FAEST_EM_256F_Senc];
  bf256_t bf_qs_dash[FAEST_EM_256F_Senc];
  em_enc_forward_256_vbb(vbb, bf_x, bf_qs);
  em_enc_backward_256_linear_access_verify(vbb, bf_x, bf_q_out, 0, 1, delta, bf_qs_dash);
  faest_aligned_free(bf_q_out);
  faest_aligned_free(bf_x);

  // Step: 13..14
  bf256_t minus_part = bf256_mul(bf_delta, bf_delta);
  for (unsigned int j = 0; j < FAEST_EM_256F_Senc; j++) {
    zk_hash_256_update(b0_ctx, bf256_add(bf256_mul(bf_qs[j], bf_qs_dash[j]), minus_part));
  }
}

static void em_prove_256(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
                         const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde) {
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

  zk_hash_256_ctx a0_ctx;
  zk_hash_256_ctx a1_ctx;

  zk_hash_256_init(&a0_ctx, chall);
  zk_hash_256_init(&a1_ctx, chall);
  em_enc_constraints_Mkey_0_256(out, x, w, vbb, &a0_ctx, &a1_ctx);

  zk_hash_256_finalize(a_tilde, &a1_ctx, bf256_load(get_vole_u(vbb) + FAEST_EM_256F_Lenc / 8));
  zk_hash_256_finalize(b_tilde, &a0_ctx, bf256_sum_poly_vbb(vbb, FAEST_EM_256F_Lenc, 0));
}

static uint8_t* em_verify_256(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                              const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out) {
  const uint8_t* delta = chall_3;

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
  em_enc_constraints_Mkey_1_256(out, x, vbb, delta, &b0_ctx);

  uint8_t* q_tilde = malloc(FAEST_EM_256F_LAMBDA / 8);
  zk_hash_256_finalize(q_tilde, &b0_ctx, bf256_sum_poly_vbb(vbb, FAEST_EM_256F_Lenc, 0));

  bf256_t bf_qtilde = bf256_load(q_tilde);
  bf256_store(q_tilde, bf256_add(bf_qtilde, bf256_mul(bf256_load(a_tilde), bf256_load(delta))));

  return q_tilde;
}

// dispatchers

void aes_prove(const uint8_t* w, vbb_t* vbb, const uint8_t* in, const uint8_t* out,
               const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params) {
  switch (params->faest_param.lambda) {
  case 256:
    if (params->faest_param.Lke) {
      aes_prove_256(w, vbb, in, out, chall, a_tilde, b_tilde, params);
    } else {
      em_prove_256(w, vbb, in, out, chall, a_tilde, b_tilde);
    }
    break;
  case 192:
    if (params->faest_param.Lke) {
      aes_prove_192(w, vbb, in, out, chall, a_tilde, b_tilde, params);
    } else {
      em_prove_192(w, vbb, in, out, chall, a_tilde, b_tilde);
    }
    break;
  default:
    if (params->faest_param.Lke) {
      aes_prove_128(w, vbb, in, out, chall, a_tilde, b_tilde, params);
    } else {
      em_prove_128(w, vbb, in, out, chall, a_tilde, b_tilde);
    }
  }
}

uint8_t* aes_verify(vbb_t* vbb, const uint8_t* chall_2, const uint8_t* chall_3,
                    const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out,
                    const faest_paramset_t* params) {
  switch (params->faest_param.lambda) {
  case 256:
    if (params->faest_param.Lke) {
      return aes_verify_256(vbb, chall_2, chall_3, a_tilde, in, out);
    } else {
      return em_verify_256(vbb, chall_2, chall_3, a_tilde, in, out);
    }
  case 192:
    if (params->faest_param.Lke) {
      return aes_verify_192(vbb, chall_2, chall_3, a_tilde, in, out);
    } else {
      return em_verify_192(vbb, chall_2, chall_3, a_tilde, in, out);
    }
  default:
    if (params->faest_param.Lke) {
      return aes_verify_128(vbb, chall_2, chall_3, a_tilde, in, out);
    } else {
      return em_verify_128(vbb, chall_2, chall_3, a_tilde, in, out);
    }
  }
}
