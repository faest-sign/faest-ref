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

#include <string.h>
#include <stdlib.h>

static const bf8_t Rcon[30] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a,
    0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91,
};

ATTR_CONST static inline bf8_t get_bit(bf8_t in, uint8_t index) {
  return (in >> index) & 0x01;
}

ATTR_CONST static inline bf8_t set_bit(bf8_t in, uint8_t index) {
  return (in << index);
}

static inline uint8_t ptr_get_bit(const uint8_t* in, unsigned int index) {
  return (in[index / 8] >> (index % 8)) & 1;
}

static inline void ptr_set_bit(uint8_t* dst, uint8_t in, unsigned int index) {
  dst[index / 8] |= in << (index % 8);
}

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

// lambda == 128 implementation

static void aes_key_schedule_forward_128(uint32_t m, const uint8_t* x, const bf128_t* v,
                                         uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                         uint8_t* out, bf128_t* bf_out,
                                         const faest_paramset_t* params) {
  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return;
  }

  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int R           = params->cipher_param.numRounds;
  const unsigned int Nwd         = params->faest_param.Nwd;
  const unsigned int lambdaBytes = lambda / 8;

  if (m == 1) {
    const unsigned int out_len = (R + 1) * 128 / 8;
    // Step 3
    memcpy(out, x, lambdaBytes);
    memset(out + lambdaBytes, 0, out_len - lambdaBytes);

    // Step: 4
    uint32_t i_wd = lambda;
    // Step: 5..10
    for (uint32_t j = Nwd; j < 4 * (R + 1); j++) {
      if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {
        memcpy(out + 32 * j / 8, x + i_wd / 8, 4);
        i_wd += 32;
      } else {
        for (uint32_t i = 0; i < 32; i += 8) {
          // bit spliced
          out[(32 * j + i) / 8] |= out[(32 * (j - Nwd) + i) / 8] ^ out[(32 * (j - 1) + i) / 8];
        }
      }
    }
    return;
  }

  for (uint32_t i = 0; i < lambdaBytes; i++) {
    for (uint32_t j = 0; j < 8; j++) {
      bf_out[(i * 8) + j] = v[i * 8 + j];
    }
  }

  // Step: 4
  uint32_t i_wd = lambda;
  // Step: 5..10
  for (uint32_t j = Nwd; j < 4 * (R + 1); j++) {
    if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {
      for (uint32_t i = (j * 32); i <= ((j * 32) + 31); i++, i_wd++) {
        bf_out[i] = v[i_wd];
      }
    } else {
      for (uint32_t i = 0; i < 32; i++) {
        bf_out[(32 * j) + i] = bf128_add(bf_out[32 * (j - Nwd) + i], bf_out[32 * (j - 1) + i]);
      }
    }
  }
}

static void aes_key_schedule_backward_128(uint32_t m, const uint8_t* x, const bf128_t* v,
                                          const uint8_t* xk, const bf128_t* Vk, uint8_t Mtag,
                                          uint8_t Mkey, const uint8_t* delta, uint8_t* out,
                                          bf128_t* bf_out, const faest_paramset_t* params) {
  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return;
  }

  const unsigned int lambda = params->faest_param.lambda;
  const unsigned int Ske    = params->faest_param.Ske;

  // STep: 2
  if (m == 1) {
    assert(delta == NULL);

    uint32_t iwd   = 0;
    uint32_t c     = 0;
    bool rmvRcon   = true;
    uint32_t ircon = 0;

    for (uint32_t j = 0; j < Ske; j++) {
      // Step 7 (bit sliced)
      uint8_t x_tilde = x[j] ^ xk[(iwd + 8 * c) / 8];

      // Step 8
      if (Mtag == 0 && rmvRcon == true && c == 0) {
        uint8_t rcon = Rcon[ircon];
        ircon        = ircon + 1;
        // Steps 12 and 13, bitsliced; delta is always 0
        x_tilde ^= rcon;
      }

      // Step: 15..19
      uint8_t y_tilde = 0;
      for (uint32_t i = 0; i < 8; ++i) {
        y_tilde ^= set_bit(get_bit(x_tilde, (i + 7) % 8) ^ get_bit(x_tilde, (i + 5) % 8) ^
                               get_bit(x_tilde, (i + 2) % 8),
                           i);
      }
      y_tilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 0);
      y_tilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 2);

      // Step: 19
      out[j] = y_tilde;
      // Step: 20
      c = c + 1;

      if (c == 4) {
        c = 0;
        if (lambda == 192) {
          iwd += 192;
        } else {
          iwd += 128;
          if (lambda == 256) {
            rmvRcon = !rmvRcon;
          }
        }
      }
    }
    return;
  }

  bf128_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf128_zero();
  } else {
    bf_delta = bf128_load(delta);
  }

  uint32_t iwd   = 0;
  uint32_t c     = 0;
  bool rmvRcon   = true;
  uint32_t ircon = 0;

  bf128_t bf_mkey             = bf128_from_bit(Mkey);
  bf128_t bf_minus_mkey       = bf128_from_bit(1 ^ Mkey);
  bf128_t bf_minus_mtag       = bf128_from_bit(1 ^ Mtag);
  bf128_t bf_mkey_times_delta = bf128_mul(bf_mkey, bf_delta);
  bf_mkey_times_delta         = bf128_add(bf_mkey_times_delta, bf_minus_mkey);

  for (uint32_t j = 0; j < Ske; j++) {
    // Step 7
    bf128_t bf_x_tilde[8];
    for (uint32_t i = 0; i < 8; i++) {
      bf_x_tilde[i] = bf128_add(v[8 * j + i], Vk[iwd + 8 * c + i]);
    }

    if (Mtag == 0 && rmvRcon == true && c == 0) {
      // Step 9
      uint8_t r = Rcon[ircon];
      ircon     = ircon + 1;

      bf128_t bf_r[8];
      for (uint32_t i = 0; i < 8; i++) {
        // Step 12
        bf_r[i] = bf128_mul(bf128_from_bit(get_bit(r, i)), bf_mkey_times_delta);
        // Step 13
        bf_x_tilde[i] = bf128_add(bf_x_tilde[i], bf_r[i]);
      }
    }

    bf128_t bf_y_tilde[8];
    for (uint32_t i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf128_add(bf128_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                bf_x_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf128_add(bf_y_tilde[0], bf128_mul(bf_minus_mtag, bf_mkey_times_delta));
    bf_y_tilde[2] = bf128_add(bf_y_tilde[2], bf128_mul(bf_minus_mtag, bf_mkey_times_delta));

    // TODO: get rid of this copy
    for (uint32_t i = 0; i < 8; i++) {
      bf_out[(8 * j) + i] = bf_y_tilde[i];
    }
    c = c + 1;

    if (c == 4) {
      c = 0;
      if (lambda == 192) {
        iwd += 192;
      } else {
        iwd += 128;
        if (lambda == 256) {
          rmvRcon = !rmvRcon;
        }
      }
    }
  }
}

static void aes_key_schedule_constraints_128(const uint8_t* w, const bf128_t* v, const uint8_t Mkey,
                                             const bf128_t* q, const uint8_t* delta, bf128_t* A0,
                                             bf128_t* A1, uint8_t* k, bf128_t* vk, bf128_t* B,
                                             bf128_t* qk, const faest_paramset_t* params) {
  const unsigned int lambda     = params->faest_param.lambda;
  const unsigned int Nwd        = params->faest_param.Nwd;
  const unsigned int Ske        = params->faest_param.Ske;
  const unsigned int lambdaByte = lambda / 8;

  if (Mkey == 0) {
    assert(delta == NULL);

    // Step: 2
    aes_key_schedule_forward_128(1, w, NULL, 0, 0, NULL, k, NULL, params);

    // Step: 3
    aes_key_schedule_forward_128(lambda, NULL, v, 1, 0, NULL, NULL, vk, params);

    // Step: 4
    uint8_t* w_dash = malloc(Ske);
    aes_key_schedule_backward_128(1, w + lambdaByte, NULL, k, NULL, 0, 0, NULL, w_dash, NULL,
                                  params);

    // Step: 5
    bf128_t* v_w_dash = malloc(Ske * 8 * sizeof(bf128_t));
    aes_key_schedule_backward_128(lambda, NULL, v + lambda, NULL, vk, 1, 0, NULL, NULL, v_w_dash,
                                  params);

    // Step: 6..8
    uint32_t iwd = 32 * (Nwd - 1);
    for (uint32_t j = 0; j < Ske / 4; j++) {
      bf128_t bf_k_hat[4];
      bf128_t bf_v_k_hat[4];
      bf128_t bf_w_dash_hat[4];
      bf128_t bf_v_w_dash_hat[4];
      for (uint32_t r = 0; r <= 3; r++) {
        // Step: 10..11
        bf_k_hat[(r + 3) % 4]   = bf128_byte_combine_bits(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat[(r + 3) % 4] = bf128_byte_combine(vk + (iwd + 8 * r));
        bf_w_dash_hat[r]        = bf128_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
        bf_v_w_dash_hat[r]      = bf128_byte_combine(v_w_dash + (32 * j + 8 * r));
      }
      // Step: 13..17
      for (uint32_t r = 0; r <= 3; r++) {
        // check that the constraint actually holds
        bf128_t should_be_zero = bf128_add(bf128_mul(bf_k_hat[r], bf_w_dash_hat[r]), bf128_one());
        assert(should_be_zero.values[0] == 0 && should_be_zero.values[1] == 0);
        A0[4 * j + r] = bf128_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
        A1[4 * j + r] =
            bf128_add(bf128_add(bf128_mul(bf128_add(bf_k_hat[r], bf_v_k_hat[r]),
                                          bf128_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                bf128_from_bf8(bf8_one())),
                      A0[4 * j + r]);
      }
      if (lambda == 192) {
        iwd = iwd + 192;
      } else {
        iwd = iwd + 128;
      }
    }
    free(v_w_dash);
    free(w_dash);
    return;
  }

  assert(delta);

  // Step: 19..20
  aes_key_schedule_forward_128(lambda, NULL, q, 0, 1, delta, NULL, qk, params);
  bf128_t* q_w_dash = malloc(Ske * 8 * sizeof(bf128_t));
  aes_key_schedule_backward_128(lambda, NULL, &q[lambda], NULL, qk, 0, 1, delta, NULL, q_w_dash,
                                params);

  const bf128_t bf_delta = bf128_load(delta);

  // Step 23..24
  uint32_t iwd = 32 * (Nwd - 1);
  for (uint32_t j = 0; j < Ske / 4; j++) {
    bf128_t bf_q_hat_k[4];
    bf128_t bf_q_hat_w_dash[4];
    for (uint32_t r = 0; r <= 3; r++) {
      // Step: 25..26
      bf_q_hat_k[(r + 3) % 4] = bf128_byte_combine(qk + ((iwd + 8 * r)));
      bf_q_hat_w_dash[r]      = bf128_byte_combine(q_w_dash + ((32 * j + 8 * r)));
    }
    // STep: 27
    for (uint32_t r = 0; r <= 3; r++) {
      bf128_t bf_tmp = bf128_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      B[4 * j + r]   = bf128_add(bf_tmp, bf128_mul(bf_delta, bf_delta));
    }
    if (lambda == 192) {
      iwd = iwd + 192;
    } else {
      iwd = iwd + 128;
    }
  }
  free(q_w_dash);
}

static void aes_enc_forward_128(uint32_t m, const uint8_t* x, const bf128_t* bf_x,
                                const uint8_t* xk, const bf128_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf128_t* bf_y,
                                const faest_paramset_t* params) {
  const unsigned int R = params->cipher_param.numRounds;

  bf128_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf128_zero();
  } else {
    bf_delta = bf128_load(delta);
  }

  if (m == 1) {
    // Step: 2
    for (uint32_t i = 0; i < 16; i++) {
      uint8_t xin = 0;
      // STep: 3
      for (uint32_t j = 0; j < 8; j++) {
        // STep: 4
        // TODO: check bit order
        xin |= set_bit(get_bit(in[i], j) & (1 ^ Mtag) & (1 ^ Mkey), j);
      }
      // STep: 5
      bf_y[i] = bf128_add(bf128_byte_combine_bits(xin), bf128_byte_combine_bits(xk[i]));
    }
    uint32_t ix, ik, iy;
    for (uint32_t j = 1; j < R; j++) {
      for (uint32_t c = 0; c <= 3; c++) {
        ix = 128 * (j - 1) + 32 * c;
        ik = 128 * j + 32 * c;
        iy = 16 * j + 4 * c;
        bf128_t bf_x_hat[4];
        bf128_t bf_xk_hat[4];
        for (uint32_t r = 0; r <= 3; r++) {
          // Step: 12..13
          bf_x_hat[r]  = bf128_byte_combine_bits(x[(ix + 8 * r) / 8]);
          bf_xk_hat[r] = bf128_byte_combine_bits(xk[(ik + 8 * r) / 8]);
        }
        bf128_t bf_one   = bf128_one();
        bf128_t bf_two   = bf128_byte_combine_bits(2);
        bf128_t bf_three = bf128_byte_combine_bits(3);
        // Step : 14
        bf_y[iy + 0] = bf128_add(bf_xk_hat[0], bf128_mul(bf_x_hat[0], bf_two));
        bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[1], bf_three));
        bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[2], bf_one));
        bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[3], bf_one));

        // STep: 15
        bf_y[iy + 1] = bf128_add(bf_xk_hat[1], bf128_mul(bf_x_hat[0], bf_one));
        bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[1], bf_two));
        bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[2], bf_three));
        bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[3], bf_one));

        // STep: 16
        bf_y[iy + 2] = bf128_add(bf_xk_hat[2], bf128_mul(bf_x_hat[0], bf_one));
        bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[1], bf_one));
        bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[2], bf_two));
        bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[3], bf_three));

        // Step: 17
        bf_y[iy + 3] = bf128_add(bf_xk_hat[3], bf128_mul(bf_x_hat[0], bf_three));
        bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[1], bf_one));
        bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[2], bf_one));
        bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[3], bf_two));
      }
    }
    return;
  }

  bf128_t bf_minus_mtag = bf128_from_bit(1 ^ Mtag);
  bf128_t bf_minus_mkey = bf128_from_bit(1 ^ Mkey);
  bf128_t bf_mkey       = bf128_from_bit(Mkey);

  // Step: 2..4
  for (uint32_t i = 0; i < 16; i++) {
    bf128_t bf_xin[8];
    for (uint32_t j = 0; j < 8; j++) {
      bf_xin[j] = bf128_mul(bf128_mul(bf128_from_bit(get_bit(in[i], j)), bf_minus_mtag),
                            bf128_add(bf128_mul(bf_mkey, bf_delta), bf_minus_mkey));
    }
    // Step: 5
    bf_y[i] = bf128_add(bf128_byte_combine(bf_xin), bf128_byte_combine(bf_xk + (8 * i)));
  }
  uint32_t ix, ik, iy;
  for (uint32_t j = 1; j < R; j++) {
    for (uint32_t c = 0; c <= 3; c++) {
      ix = 128 * (j - 1) + 32 * c;
      ik = 128 * j + 32 * c;
      iy = 16 * j + 4 * c;
      bf128_t bf_x_hat[4];
      bf128_t bf_xk_hat[4];
      for (uint32_t r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf128_byte_combine(bf_x + (ix + 8 * r));
        bf_xk_hat[r] = bf128_byte_combine(bf_xk + (ik + 8 * r));
      }
      bf128_t bf_one   = bf128_one();
      bf128_t bf_two   = bf128_byte_combine_bits(2);
      bf128_t bf_three = bf128_byte_combine_bits(3);
      // Step : 14
      bf_y[iy + 0] = bf128_add(bf_xk_hat[0], bf128_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[2], bf_one));
      bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[3], bf_one));

      // STep: 15
      bf_y[iy + 1] = bf128_add(bf_xk_hat[1], bf128_mul(bf_x_hat[0], bf_one));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[3], bf_one));

      // STep: 16
      bf_y[iy + 2] = bf128_add(bf_xk_hat[2], bf128_mul(bf_x_hat[0], bf_one));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[1], bf_one));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf128_add(bf_xk_hat[3], bf128_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[1], bf_one));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[2], bf_one));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[3], bf_two));
    }
  }
}

static void aes_enc_backward_128(uint32_t m, const uint8_t* x, const bf128_t* bf_x,
                                 const uint8_t* xk, const bf128_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf128_t* y_out, const faest_paramset_t* params) {
  const unsigned int R = params->cipher_param.numRounds;

  if (m == 1) {
    uint8_t xtilde;
    // Step:2..4
    for (uint32_t j = 0; j < R; j++) {
      for (uint32_t c = 0; c <= 3; c++) {
        for (uint32_t r = 0; r <= 3; r++) {
          // Step: 5..6
          unsigned int ird = (128 * j) + (32 * ((c - r) % 4)) + (8 * r);
          if (j < (R - 1)) {
            // Step: 7
            xtilde = x[ird / 8];
          } else {
            // Step: 9
            uint8_t xout = 0;
            for (uint32_t i = 0; i < 8; i++) {
              // Step: 10..11
              // delta is always \bot if called with m == 1
              // TODO bit splice
              xout |=
                  set_bit(get_bit(out[(ird - 128 * (R - 1)) / 8], i) & (1 ^ Mtag) & (1 ^ Mkey), i);
            }
            xtilde = xout ^ xk[(128 + ird) / 8];
          }
          // Step: 12..17
          uint8_t ytilde = 0;
          for (uint32_t i = 0; i < 8; ++i) {
            ytilde ^= set_bit(get_bit(xtilde, (i + 7) % 8) ^ get_bit(xtilde, (i + 5) % 8) ^
                                  get_bit(xtilde, (i + 2) % 8),
                              i);
          }
          ytilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 0);
          ytilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 2);

          // Step: 18
          y_out[16 * j + 4 * c + r] = bf128_byte_combine_bits(ytilde);
        }
      }
    }
    return;
  }

  // Step: 1
  bf128_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf128_zero();
  } else {
    bf_delta = bf128_load(delta);
  }

  const bf128_t factor =
      bf128_mul(bf128_from_bit(1 ^ Mtag),
                bf128_add(bf128_mul(bf128_from_bit(Mkey), bf_delta), bf128_from_bit(1 ^ Mkey)));

  // STep: 2..4
  for (uint32_t j = 0; j < R; j++) {
    for (uint32_t c = 0; c <= 3; c++) {
      for (uint32_t r = 0; r <= 3; r++) {
        bf128_t bf_x_tilde[8];
        // Step: 5
        unsigned int ird = (128 * j) + (32 * ((c - r) % 4)) + (8 * r);
        // Step: 6
        if (j < (R - 1)) {
          // Step: 7
          for (uint32_t i = 0; i < 8; i++) {
            bf_x_tilde[i] = bf_x[ird + i];
          }
        } else {
          bf128_t bf_xout[8];
          // Step: 10
          for (uint32_t i = 0; i < 8; ++i) {
            // Step: 11
            // TODO: check 0 extension
            bf_xout[i] = bf128_from_bit(get_bit(out[(ird - 128 * (R - 1)) / 8], i));
            bf_xout[i] = bf128_mul(bf_xout[i], factor);
            // Step: 12
            bf_x_tilde[i] = bf128_add(bf_xout[i], bf_xk[128 + ird + i]);
          }
        }
        // STep: 13..17
        bf128_t bf_y_tilde[8];
        for (uint32_t i = 0; i < 8; ++i) {
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
                                    const uint8_t* delta, bf128_t* A0, bf128_t* A1, bf128_t* B,
                                    const faest_paramset_t* params) {
  const unsigned int lambda = params->faest_param.lambda;
  const unsigned int Senc   = params->faest_param.Senc;

  if (Mkey == 0) {
    bf128_t* s       = malloc(sizeof(bf128_t) * Senc);
    bf128_t* vs      = malloc(sizeof(bf128_t) * Senc);
    bf128_t* s_dash  = malloc(sizeof(bf128_t) * Senc);
    bf128_t* vs_dash = malloc(sizeof(bf128_t) * Senc);
    aes_enc_forward_128(1, w, NULL, k, NULL, in, 0, 0, NULL, s, params);
    aes_enc_forward_128(lambda, NULL, v, NULL, vk, in, 1, 0, NULL, vs, params);
    aes_enc_backward_128(1, w, NULL, k, NULL, 0, 0, NULL, out, s_dash, params);
    aes_enc_backward_128(lambda, NULL, v, NULL, vk, 1, 0, NULL, out, vs_dash, params);

    for (uint32_t j = 0; j < Senc; j++) {
      // check that the constraint actually holds
      bf128_t should_be_zero = bf128_add(bf128_mul(s[j], s_dash[j]), bf128_one());
      assert(should_be_zero.values[0] == 0 && should_be_zero.values[1] == 0);
      A0[j] = bf128_mul(vs[j], vs_dash[j]);
      A1[j] = bf128_add(
          bf128_add(bf128_mul(bf128_add(s[j], vs[j]), bf128_add(s_dash[j], vs_dash[j])), A0[j]),
          bf128_one());
    }

    free(vs_dash);
    free(s_dash);
    free(vs);
    free(s);
  } else {
    // Step: 11..12
    bf128_t* qs      = malloc(sizeof(bf128_t) * Senc);
    bf128_t* qs_dash = malloc(sizeof(bf128_t) * Senc);
    aes_enc_forward_128(lambda, NULL, q, NULL, qk, in, 0, 1, delta, qs, params);
    aes_enc_backward_128(lambda, NULL, q, NULL, qk, 0, 1, delta, out, qs_dash, params);

    // Step: 13..14
    bf128_t minus_part = bf128_mul(bf128_load(delta), bf128_load(delta));
    for (uint32_t j = 0; j < Senc; j++) {
      B[j] = bf128_add(bf128_mul(qs[j], qs_dash[j]), minus_part);
    }
    free(qs);
    free(qs_dash);
  }
}

static void aes_prove_128(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
                          const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde,
                          uint8_t* b_tilde, const faest_paramset_t* params) {
  const unsigned int beta = params->faest_param.beta;
  const unsigned int l    = params->faest_param.l;
  const unsigned int Lke  = params->faest_param.Lke;
  const unsigned int Lenc = params->faest_param.Lenc;
  const unsigned int R    = params->cipher_param.numRounds;
  const unsigned int Ske  = params->faest_param.Ske;
  const unsigned int Senc = params->faest_param.Senc;

  // Step: 1..2
  bf128_t* bf_v = column_to_row_major_and_shrink_V_128(V, l);

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7
  const unsigned int length_a = Ske + (beta * Senc) + 1;
  bf128_t* A0                 = malloc(sizeof(bf128_t) * length_a);
  bf128_t* A1                 = malloc(sizeof(bf128_t) * length_a);
  uint8_t* k                  = malloc((R + 1) * 128 / 8);
  bf128_t* vk                 = malloc(sizeof(bf128_t) * ((R + 1) * 128));
  bf128_t* qk                 = malloc(sizeof(bf128_t) * ((R + 1) * 128));
  if (Lke > 0) {
    aes_key_schedule_constraints_128(w, bf_v, 0, NULL, NULL, A0, A1, k, vk, NULL, qk, params);
  }

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  aes_enc_constraints_128(in, out, w + Lke / 8, bf_v + Lke, k, vk, 0, NULL, NULL, NULL, A0 + Ske,
                          A1 + Ske, NULL, params);
  // Step: 12
  if (beta == 2) {
    // Step: 13
    // Step: 14, 15
    aes_enc_constraints_128(in + 16, out + 16, w + (Lke + Lenc) / 8, bf_v + Lke + Lenc, k, vk, 0,
                            NULL, NULL, NULL, A0 + (Ske + Senc), A1 + (Ske + Senc), NULL, params);
  }
  free(qk);
  free(vk);
  free(k);

  // Step: 16..18
  A1[length_a - 1] = bf128_load(u + l / 8);
  A0[length_a - 1] = bf128_sum_poly(bf_v + l);
  free(bf_v);

  zk_hash_128(a_tilde, chall, A1, length_a - 1);
  zk_hash_128(b_tilde, chall, A0, length_a - 1);

  free(A0);
  free(A1);
}

static uint8_t* aes_verify_128(uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                               const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                               const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int tau         = params->faest_param.tau;
  const unsigned int t0          = params->faest_param.t0;
  const unsigned int k0          = params->faest_param.k0;
  const unsigned int t1          = params->faest_param.t1;
  const unsigned int k1          = params->faest_param.k1;
  const unsigned int beta        = params->faest_param.beta;
  const unsigned int l           = params->faest_param.l;
  const unsigned int Lke         = params->faest_param.Lke;
  const unsigned int Lenc        = params->faest_param.Lenc;
  const unsigned int R           = params->cipher_param.numRounds;
  const unsigned int Ske         = params->faest_param.Ske;
  const unsigned int Senc        = params->faest_param.Senc;
  const unsigned int lambdaBytes = lambda / 8;

  // Step: 1
  const uint8_t* delta = chall_3;
  // Step: 2,3
  // do nothing

  // Step: 4..10
  for (uint32_t i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t* fancy_d   = malloc(depth);
    ChalDec(chall_3, i, k0, t0, k1, t1, fancy_d);
    for (uint32_t j = 0; j < depth; j++, ++col) {
      if (fancy_d[j] == 1) {
        xorUint8Arr(d, Q[col], Q[col], (l + 7) / 8);
      }
    }
    free(fancy_d);
  }

  // Step: 11..12
  bf128_t* bf_q = column_to_row_major_and_shrink_V_128(Q, l);

  // Step: 13
  const unsigned int length_b = Ske + (beta * Senc) + 1;
  uint8_t* k                  = malloc((R + 1) * 128);
  bf128_t* vk                 = malloc(sizeof(bf128_t) * ((R + 1) * 128));
  bf128_t* qk                 = malloc(sizeof(bf128_t) * ((R + 1) * 128));
  bf128_t* B_0                = malloc(sizeof(bf128_t) * length_b);
  if (Lke > 0) {
    aes_key_schedule_constraints_128(NULL, NULL, 1, bf_q, delta, NULL, NULL, k, vk, B_0, qk,
                                     params);
  }

  // Step: 14
  bf128_t* B_1 = B_0 + Ske;
  aes_enc_constraints_128(in, out, NULL, NULL, NULL, NULL, 1, bf_q + Lke, qk, delta, NULL, NULL,
                          B_1, params);

  if (beta == 2) {
    // Step: 18
    bf128_t* B_2 = B_0 + (Ske + Senc);
    aes_enc_constraints_128(in + 16, out + 16, NULL, NULL, NULL, NULL, 1, bf_q + (Lke + Lenc), qk,
                            delta, NULL, NULL, B_2, params);
  }
  free(qk);
  free(vk);
  free(k);

  // Step: 20
  B_0[length_b - 1] = bf128_sum_poly(bf_q + l);
  free(bf_q);

  // Step 21
  uint8_t* q_tilde = malloc(lambdaBytes);
  zk_hash_128(q_tilde, chall_2, B_0, length_b - 1);
  free(B_0);

  bf128_t bf_qtilde = bf128_load(q_tilde);
  bf128_store(q_tilde, bf128_add(bf_qtilde, bf128_mul(bf128_load(a_tilde), bf128_load(delta))));

  return q_tilde;
}

// lambda == 192 implementation

static void aes_key_schedule_forward_192(uint32_t m, const uint8_t* x, const bf192_t* v,
                                         uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                         uint8_t* out, bf192_t* bf_out,
                                         const faest_paramset_t* params) {
  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return;
  }

  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int R           = params->cipher_param.numRounds;
  const unsigned int Nwd         = params->faest_param.Nwd;
  const unsigned int lambdaBytes = lambda / 8;

  if (m == 1) {
    const unsigned int out_len = (R + 1) * 128 / 8;
    // Step 3
    memcpy(out, x, lambdaBytes);
    memset(out + lambdaBytes, 0, out_len - lambdaBytes);

    // Step: 4
    uint32_t i_wd = lambda;
    // Step: 5..10
    for (uint32_t j = Nwd; j < 4 * (R + 1); j++) {
      if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {
        memcpy(out + 32 * j / 8, x + i_wd / 8, 4);
        i_wd += 32;
      } else {
        for (uint32_t i = 0; i < 32; i += 8) {
          // bit spliced
          out[(32 * j + i) / 8] |= out[(32 * (j - Nwd) + i) / 8] ^ out[(32 * (j - 1) + i) / 8];
        }
      }
    }
    return;
  }

  for (uint32_t i = 0; i < lambdaBytes; i++) {
    for (uint32_t j = 0; j < 8; j++) {
      bf_out[(i * 8) + j] = v[i * 8 + j];
    }
  }

  // Step: 4
  uint32_t i_wd = lambda;
  // Step: 5..10
  for (uint32_t j = Nwd; j < 4 * (R + 1); j++) {
    if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {
      for (uint32_t i = (j * 32); i <= ((j * 32) + 31); i++, i_wd++) {
        bf_out[i] = v[i_wd];
      }
    } else {
      for (uint32_t i = 0; i < 32; i++) {
        bf_out[(32 * j) + i] = bf192_add(bf_out[32 * (j - Nwd) + i], bf_out[32 * (j - 1) + i]);
      }
    }
  }
}

static void aes_key_schedule_backward_192(uint32_t m, const uint8_t* x, const bf192_t* v,
                                          const uint8_t* xk, const bf192_t* Vk, uint8_t Mtag,
                                          uint8_t Mkey, const uint8_t* delta, uint8_t* out,
                                          bf192_t* bf_out, const faest_paramset_t* params) {
  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return;
  }

  const unsigned int lambda = params->faest_param.lambda;
  const unsigned int Ske    = params->faest_param.Ske;

  // STep: 2
  if (m == 1) {
    assert(delta == NULL);

    uint32_t iwd   = 0;
    uint32_t c     = 0;
    bool rmvRcon   = true;
    uint32_t ircon = 0;

    for (uint32_t j = 0; j < Ske; j++) {
      // Step 7 (bit sliced)
      uint8_t x_tilde = x[j] ^ xk[(iwd + 8 * c) / 8];

      // Step 8
      if (Mtag == 0 && rmvRcon == true && c == 0) {
        uint8_t rcon = Rcon[ircon];
        ircon        = ircon + 1;
        // Steps 12 and 13, bitsliced; delta is always 0
        x_tilde ^= rcon;
      }

      // Step: 15..19
      uint8_t y_tilde = 0;
      for (uint32_t i = 0; i < 8; ++i) {
        y_tilde ^= set_bit(get_bit(x_tilde, (i + 7) % 8) ^ get_bit(x_tilde, (i + 5) % 8) ^
                               get_bit(x_tilde, (i + 2) % 8),
                           i);
      }
      y_tilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 0);
      y_tilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 2);

      // Step: 19
      out[j] = y_tilde;
      // Step: 20
      c = c + 1;

      if (c == 4) {
        c = 0;
        if (lambda == 192) {
          iwd += 192;
        } else {
          iwd += 128;
          if (lambda == 256) {
            rmvRcon = !rmvRcon;
          }
        }
      }
    }
    return;
  }

  bf192_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf192_zero();
  } else {
    bf_delta = bf192_load(delta);
  }

  uint32_t iwd   = 0;
  uint32_t c     = 0;
  bool rmvRcon   = true;
  uint32_t ircon = 0;

  bf192_t bf_mkey             = bf192_from_bit(Mkey);
  bf192_t bf_minus_mkey       = bf192_from_bit(1 ^ Mkey);
  bf192_t bf_minus_mtag       = bf192_from_bit(1 ^ Mtag);
  bf192_t bf_mkey_times_delta = bf192_mul(bf_mkey, bf_delta);
  bf_mkey_times_delta         = bf192_add(bf_mkey_times_delta, bf_minus_mkey);

  for (uint32_t j = 0; j < Ske; j++) {
    // Step 7
    bf192_t bf_x_tilde[8];
    for (uint32_t i = 0; i < 8; i++) {
      bf_x_tilde[i] = bf192_add(v[8 * j + i], Vk[iwd + 8 * c + i]);
    }

    if (Mtag == 0 && rmvRcon == true && c == 0) {
      // Step 9
      uint8_t r = Rcon[ircon];
      ircon     = ircon + 1;

      bf192_t bf_r[8];
      for (uint32_t i = 0; i < 8; i++) {
        // Step 12
        bf_r[i] = bf192_mul(bf192_from_bit(get_bit(r, i)), bf_mkey_times_delta);
        // Step 13
        bf_x_tilde[i] = bf192_add(bf_x_tilde[i], bf_r[i]);
      }
    }

    bf192_t bf_y_tilde[8];
    for (uint32_t i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf192_add(bf192_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                bf_x_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf192_add(bf_y_tilde[0], bf192_mul(bf_minus_mtag, bf_mkey_times_delta));
    bf_y_tilde[2] = bf192_add(bf_y_tilde[2], bf192_mul(bf_minus_mtag, bf_mkey_times_delta));

    // TODO: get rid of this copy
    for (uint32_t i = 0; i < 8; i++) {
      bf_out[(8 * j) + i] = bf_y_tilde[i];
    }
    c = c + 1;

    if (c == 4) {
      c = 0;
      if (lambda == 192) {
        iwd += 192;
      } else {
        iwd += 128;
        if (lambda == 256) {
          rmvRcon = !rmvRcon;
        }
      }
    }
  }
}

static void aes_key_schedule_constraints_192(const uint8_t* w, const bf192_t* v, const uint8_t Mkey,
                                             const bf192_t* q, const uint8_t* delta, bf192_t* A0,
                                             bf192_t* A1, uint8_t* k, bf192_t* vk, bf192_t* B,
                                             bf192_t* qk, const faest_paramset_t* params) {
  const unsigned int lambda     = params->faest_param.lambda;
  const unsigned int Nwd        = params->faest_param.Nwd;
  const unsigned int Ske        = params->faest_param.Ske;
  const unsigned int lambdaByte = lambda / 8;

  if (Mkey == 0) {
    assert(delta == NULL);

    // Step: 2
    aes_key_schedule_forward_192(1, w, NULL, 0, 0, NULL, k, NULL, params);

    // Step: 3
    aes_key_schedule_forward_192(lambda, NULL, v, 1, 0, NULL, NULL, vk, params);

    // Step: 4
    uint8_t* w_dash = malloc(Ske);
    aes_key_schedule_backward_192(1, w + lambdaByte, NULL, k, NULL, 0, 0, NULL, w_dash, NULL,
                                  params);

    // Step: 5
    bf192_t* v_w_dash = malloc(Ske * 8 * sizeof(bf192_t));
    aes_key_schedule_backward_192(lambda, NULL, v + lambda, NULL, vk, 1, 0, NULL, NULL, v_w_dash,
                                  params);

    // Step: 6..8
    uint32_t iwd = 32 * (Nwd - 1);
    for (uint32_t j = 0; j < Ske / 4; j++) {
      bf192_t bf_k_hat[4];
      bf192_t bf_v_k_hat[4];
      bf192_t bf_w_dash_hat[4];
      bf192_t bf_v_w_dash_hat[4];
      for (uint32_t r = 0; r <= 3; r++) {
        // Step: 10..11
        bf_k_hat[(r + 3) % 4]   = bf192_byte_combine_bits(k[(iwd + 8 * r) / 8]);
        bf_v_k_hat[(r + 3) % 4] = bf192_byte_combine(vk + (iwd + 8 * r));
        bf_w_dash_hat[r]        = bf192_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
        bf_v_w_dash_hat[r]      = bf192_byte_combine(v_w_dash + (32 * j + 8 * r));
      }
      // Step: 13..17
      for (uint32_t r = 0; r <= 3; r++) {
        // check that the constraint actually holds
        bf192_t should_be_zero = bf192_add(bf192_mul(bf_k_hat[r], bf_w_dash_hat[r]), bf192_one());
        assert(should_be_zero.values[0] == 0 && should_be_zero.values[1] == 0 &&
               should_be_zero.values[2] == 0);
        A0[4 * j + r] = bf192_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
        A1[4 * j + r] =
            bf192_add(bf192_add(bf192_mul(bf192_add(bf_k_hat[r], bf_v_k_hat[r]),
                                          bf192_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                bf192_from_bf8(bf8_one())),
                      A0[4 * j + r]);
      }
      if (lambda == 192) {
        iwd = iwd + 192;
      } else {
        iwd = iwd + 128;
      }
    }
    free(v_w_dash);
    free(w_dash);
    return;
  }

  assert(delta);

  // Step: 19..20
  aes_key_schedule_forward_192(lambda, NULL, q, 0, 1, delta, NULL, qk, params);
  bf192_t* q_w_dash = malloc(Ske * 8 * sizeof(bf192_t));
  aes_key_schedule_backward_192(lambda, NULL, &q[lambda], NULL, qk, 0, 1, delta, NULL, q_w_dash,
                                params);

  const bf192_t bf_delta = bf192_load(delta);

  // Step 23..24
  uint32_t iwd = 32 * (Nwd - 1);
  for (uint32_t j = 0; j < Ske / 4; j++) {
    bf192_t bf_q_hat_k[4];
    bf192_t bf_q_hat_w_dash[4];
    for (uint32_t r = 0; r <= 3; r++) {
      // Step: 25..26
      bf_q_hat_k[(r + 3) % 4] = bf192_byte_combine(qk + ((iwd + 8 * r)));
      bf_q_hat_w_dash[r]      = bf192_byte_combine(q_w_dash + ((32 * j + 8 * r)));
    }
    // STep: 27
    for (uint32_t r = 0; r <= 3; r++) {
      bf192_t bf_tmp = bf192_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      B[4 * j + r]   = bf192_add(bf_tmp, bf192_mul(bf_delta, bf_delta));
    }
    if (lambda == 192) {
      iwd = iwd + 192;
    } else {
      iwd = iwd + 128;
    }
  }
  free(q_w_dash);
}

static void aes_enc_forward_192(uint32_t m, const uint8_t* x, const bf192_t* bf_x,
                                const uint8_t* xk, const bf192_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf192_t* bf_y,
                                const faest_paramset_t* params) {
  const unsigned int R = params->cipher_param.numRounds;

  bf192_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf192_zero();
  } else {
    bf_delta = bf192_load(delta);
  }

  if (m == 1) {
    // Step: 2
    for (uint32_t i = 0; i < 16; i++) {
      uint8_t xin = 0;
      // STep: 3
      for (uint32_t j = 0; j < 8; j++) {
        // STep: 4
        // TODO: check bit order
        xin |= set_bit(get_bit(in[i], j) & (1 ^ Mtag) & (1 ^ Mkey), j);
      }
      // STep: 5
      bf_y[i] = bf192_add(bf192_byte_combine_bits(xin), bf192_byte_combine_bits(xk[i]));
    }
    uint32_t ix, ik, iy;
    for (uint32_t j = 1; j < R; j++) {
      for (uint32_t c = 0; c <= 3; c++) {
        ix = 128 * (j - 1) + 32 * c;
        ik = 128 * j + 32 * c;
        iy = 16 * j + 4 * c;
        bf192_t bf_x_hat[4];
        bf192_t bf_xk_hat[4];
        for (uint32_t r = 0; r <= 3; r++) {
          // Step: 12..13
          bf_x_hat[r]  = bf192_byte_combine_bits(x[(ix + 8 * r) / 8]);
          bf_xk_hat[r] = bf192_byte_combine_bits(xk[(ik + 8 * r) / 8]);
        }
        bf192_t bf_one   = bf192_one();
        bf192_t bf_two   = bf192_byte_combine_bits(2);
        bf192_t bf_three = bf192_byte_combine_bits(3);
        // Step : 14
        bf_y[iy + 0] = bf192_add(bf_xk_hat[0], bf192_mul(bf_x_hat[0], bf_two));
        bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf192_mul(bf_x_hat[1], bf_three));
        bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf192_mul(bf_x_hat[2], bf_one));
        bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf192_mul(bf_x_hat[3], bf_one));

        // STep: 15
        bf_y[iy + 1] = bf192_add(bf_xk_hat[1], bf192_mul(bf_x_hat[0], bf_one));
        bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[1], bf_two));
        bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[2], bf_three));
        bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[3], bf_one));

        // STep: 16
        bf_y[iy + 2] = bf192_add(bf_xk_hat[2], bf192_mul(bf_x_hat[0], bf_one));
        bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[1], bf_one));
        bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[2], bf_two));
        bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[3], bf_three));

        // Step: 17
        bf_y[iy + 3] = bf192_add(bf_xk_hat[3], bf192_mul(bf_x_hat[0], bf_three));
        bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_x_hat[1], bf_one));
        bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_x_hat[2], bf_one));
        bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_x_hat[3], bf_two));
      }
    }
    return;
  }

  bf192_t bf_minus_mtag = bf192_from_bit(1 ^ Mtag);
  bf192_t bf_minus_mkey = bf192_from_bit(1 ^ Mkey);
  bf192_t bf_mkey       = bf192_from_bit(Mkey);

  // Step: 2..4
  for (uint32_t i = 0; i < 16; i++) {
    bf192_t bf_xin[8];
    for (uint32_t j = 0; j < 8; j++) {
      bf_xin[j] = bf192_mul(bf192_mul(bf192_from_bit(get_bit(in[i], j)), bf_minus_mtag),
                            bf192_add(bf192_mul(bf_mkey, bf_delta), bf_minus_mkey));
    }
    // Step: 5
    bf_y[i] = bf192_add(bf192_byte_combine(bf_xin), bf192_byte_combine(bf_xk + (8 * i)));
  }
  uint32_t ix, ik, iy;
  for (uint32_t j = 1; j < R; j++) {
    for (uint32_t c = 0; c <= 3; c++) {
      ix = 128 * (j - 1) + 32 * c;
      ik = 128 * j + 32 * c;
      iy = 16 * j + 4 * c;
      bf192_t bf_x_hat[4];
      bf192_t bf_xk_hat[4];
      for (uint32_t r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf192_byte_combine(bf_x + (ix + 8 * r));
        bf_xk_hat[r] = bf192_byte_combine(bf_xk + (ik + 8 * r));
      }
      bf192_t bf_one   = bf192_one();
      bf192_t bf_two   = bf192_byte_combine_bits(2);
      bf192_t bf_three = bf192_byte_combine_bits(3);
      // Step : 14
      bf_y[iy + 0] = bf192_add(bf_xk_hat[0], bf192_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf192_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf192_mul(bf_x_hat[2], bf_one));
      bf_y[iy + 0] = bf192_add(bf_y[iy + 0], bf192_mul(bf_x_hat[3], bf_one));

      // STep: 15
      bf_y[iy + 1] = bf192_add(bf_xk_hat[1], bf192_mul(bf_x_hat[0], bf_one));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf192_add(bf_y[iy + 1], bf192_mul(bf_x_hat[3], bf_one));

      // STep: 16
      bf_y[iy + 2] = bf192_add(bf_xk_hat[2], bf192_mul(bf_x_hat[0], bf_one));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[1], bf_one));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf192_add(bf_y[iy + 2], bf192_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf192_add(bf_xk_hat[3], bf192_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_x_hat[1], bf_one));
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_x_hat[2], bf_one));
      bf_y[iy + 3] = bf192_add(bf_y[iy + 3], bf192_mul(bf_x_hat[3], bf_two));
    }
  }
}

static void aes_enc_backward_192(uint32_t m, const uint8_t* x, const bf192_t* bf_x,
                                 const uint8_t* xk, const bf192_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf192_t* y_out, const faest_paramset_t* params) {
  const unsigned int R = params->cipher_param.numRounds;

  if (m == 1) {
    uint8_t xtilde;
    // Step:2..4
    for (uint32_t j = 0; j < R; j++) {
      for (uint32_t c = 0; c <= 3; c++) {
        for (uint32_t r = 0; r <= 3; r++) {
          // Step: 5..6
          unsigned int ird = (128 * j) + (32 * ((c - r) % 4)) + (8 * r);
          if (j < (R - 1)) {
            // Step: 7
            xtilde = x[ird / 8];
          } else {
            // Step: 9
            uint8_t xout = 0;
            for (uint32_t i = 0; i < 8; i++) {
              // Step: 10..11
              // delta is always \bot if called with m == 1
              // TODO bit splice
              xout |=
                  set_bit(get_bit(out[(ird - 128 * (R - 1)) / 8], i) & (1 ^ Mtag) & (1 ^ Mkey), i);
            }
            xtilde = xout ^ xk[(128 + ird) / 8];
          }
          // Step: 12..17
          uint8_t ytilde = 0;
          for (uint32_t i = 0; i < 8; ++i) {
            ytilde ^= set_bit(get_bit(xtilde, (i + 7) % 8) ^ get_bit(xtilde, (i + 5) % 8) ^
                                  get_bit(xtilde, (i + 2) % 8),
                              i);
          }
          ytilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 0);
          ytilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 2);

          // Step: 18
          y_out[16 * j + 4 * c + r] = bf192_byte_combine_bits(ytilde);
        }
      }
    }
    return;
  }

  // Step: 1
  bf192_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf192_zero();
  } else {
    bf_delta = bf192_load(delta);
  }

  const bf192_t factor =
      bf192_mul(bf192_from_bit(1 ^ Mtag),
                bf192_add(bf192_mul(bf192_from_bit(Mkey), bf_delta), bf192_from_bit(1 ^ Mkey)));

  // STep: 2..4
  for (uint32_t j = 0; j < R; j++) {
    for (uint32_t c = 0; c <= 3; c++) {
      for (uint32_t r = 0; r <= 3; r++) {
        bf192_t bf_x_tilde[8];
        // Step: 5
        unsigned int ird = (128 * j) + (32 * ((c - r) % 4)) + (8 * r);
        // Step: 6
        if (j < (R - 1)) {
          // Step: 7
          for (uint32_t i = 0; i < 8; i++) {
            bf_x_tilde[i] = bf_x[ird + i];
          }
        } else {
          bf192_t bf_xout[8];
          // Step: 10
          for (uint32_t i = 0; i < 8; ++i) {
            // Step: 11
            // TODO: check 0 extension
            bf_xout[i] = bf192_from_bit(get_bit(out[(ird - 128 * (R - 1)) / 8], i));
            bf_xout[i] = bf192_mul(bf_xout[i], factor);
            // Step: 12
            bf_x_tilde[i] = bf192_add(bf_xout[i], bf_xk[128 + ird + i]);
          }
        }
        // STep: 13..17
        bf192_t bf_y_tilde[8];
        for (uint32_t i = 0; i < 8; ++i) {
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
                                    const uint8_t* delta, bf192_t* A0, bf192_t* A1, bf192_t* B,
                                    const faest_paramset_t* params) {
  const unsigned int lambda = params->faest_param.lambda;
  const unsigned int Senc   = params->faest_param.Senc;

  if (Mkey == 0) {
    bf192_t* s       = malloc(sizeof(bf192_t) * Senc);
    bf192_t* vs      = malloc(sizeof(bf192_t) * Senc);
    bf192_t* s_dash  = malloc(sizeof(bf192_t) * Senc);
    bf192_t* vs_dash = malloc(sizeof(bf192_t) * Senc);
    aes_enc_forward_192(1, w, NULL, k, NULL, in, 0, 0, NULL, s, params);
    aes_enc_forward_192(lambda, NULL, v, NULL, vk, in, 1, 0, NULL, vs, params);
    aes_enc_backward_192(1, w, NULL, k, NULL, 0, 0, NULL, out, s_dash, params);
    aes_enc_backward_192(lambda, NULL, v, NULL, vk, 1, 0, NULL, out, vs_dash, params);

    for (uint32_t j = 0; j < Senc; j++) {
      // check that the constraint actually holds
      bf192_t should_be_zero = bf192_add(bf192_mul(s[j], s_dash[j]), bf192_one());
      assert(should_be_zero.values[0] == 0 && should_be_zero.values[1] == 0 &&
             should_be_zero.values[2] == 0);
      A0[j] = bf192_mul(vs[j], vs_dash[j]);
      A1[j] = bf192_add(
          bf192_add(bf192_mul(bf192_add(s[j], vs[j]), bf192_add(s_dash[j], vs_dash[j])), A0[j]),
          bf192_one());
    }

    free(vs_dash);
    free(s_dash);
    free(vs);
    free(s);
  } else {
    // Step: 11..12
    bf192_t* qs      = malloc(sizeof(bf192_t) * Senc);
    bf192_t* qs_dash = malloc(sizeof(bf192_t) * Senc);
    aes_enc_forward_192(lambda, NULL, q, NULL, qk, in, 0, 1, delta, qs, params);
    aes_enc_backward_192(lambda, NULL, q, NULL, qk, 0, 1, delta, out, qs_dash, params);

    // Step: 13..14
    bf192_t minus_part = bf192_mul(bf192_load(delta), bf192_load(delta));
    for (uint32_t j = 0; j < Senc; j++) {
      B[j] = bf192_add(bf192_mul(qs[j], qs_dash[j]), minus_part);
    }
    free(qs);
    free(qs_dash);
  }
}

static void aes_prove_192(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
                          const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde,
                          uint8_t* b_tilde, const faest_paramset_t* params) {
  const unsigned int beta = params->faest_param.beta;
  const unsigned int l    = params->faest_param.l;
  const unsigned int Lke  = params->faest_param.Lke;
  const unsigned int Lenc = params->faest_param.Lenc;
  const unsigned int R    = params->cipher_param.numRounds;
  const unsigned int Ske  = params->faest_param.Ske;
  const unsigned int Senc = params->faest_param.Senc;

  // Step: 1..2
  bf192_t* bf_v = column_to_row_major_and_shrink_V_192(V, l);

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7
  const unsigned int length_a = Ske + (beta * Senc) + 1;
  bf192_t* A0                 = malloc(sizeof(bf192_t) * length_a);
  bf192_t* A1                 = malloc(sizeof(bf192_t) * length_a);
  uint8_t* k                  = malloc((R + 1) * 128 / 8);
  bf192_t* vk                 = malloc(sizeof(bf192_t) * ((R + 1) * 128));
  bf192_t* qk                 = malloc(sizeof(bf192_t) * ((R + 1) * 128));
  if (Lke > 0) {
    aes_key_schedule_constraints_192(w, bf_v, 0, NULL, NULL, A0, A1, k, vk, NULL, qk, params);
  }

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  aes_enc_constraints_192(in, out, w + Lke / 8, bf_v + Lke, k, vk, 0, NULL, NULL, NULL, A0 + Ske,
                          A1 + Ske, NULL, params);
  // Step: 12
  if (beta == 2) {
    // Step: 13
    // Step: 14, 15
    aes_enc_constraints_192(in + 16, out + 16, w + (Lke + Lenc) / 8, bf_v + Lke + Lenc, k, vk, 0,
                            NULL, NULL, NULL, A0 + (Ske + Senc), A1 + (Ske + Senc), NULL, params);
  }
  free(qk);
  free(vk);
  free(k);

  // Step: 16..18
  A1[length_a - 1] = bf192_load(u + l / 8);
  A0[length_a - 1] = bf192_sum_poly(bf_v + l);
  free(bf_v);

  zk_hash_192(a_tilde, chall, A1, length_a - 1);
  zk_hash_192(b_tilde, chall, A0, length_a - 1);

  free(A0);
  free(A1);
}

static uint8_t* aes_verify_192(uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                               const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                               const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int tau         = params->faest_param.tau;
  const unsigned int t0          = params->faest_param.t0;
  const unsigned int k0          = params->faest_param.k0;
  const unsigned int t1          = params->faest_param.t1;
  const unsigned int k1          = params->faest_param.k1;
  const unsigned int beta        = params->faest_param.beta;
  const unsigned int l           = params->faest_param.l;
  const unsigned int Lke         = params->faest_param.Lke;
  const unsigned int Lenc        = params->faest_param.Lenc;
  const unsigned int R           = params->cipher_param.numRounds;
  const unsigned int Ske         = params->faest_param.Ske;
  const unsigned int Senc        = params->faest_param.Senc;
  const unsigned int lambdaBytes = lambda / 8;

  // Step: 1
  const uint8_t* delta = chall_3;
  // Step: 2,3
  // do nothing

  // Step: 4..10
  for (uint32_t i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t* fancy_d   = malloc(depth);
    ChalDec(chall_3, i, k0, t0, k1, t1, fancy_d);
    for (uint32_t j = 0; j < depth; j++, ++col) {
      if (fancy_d[j] == 1) {
        xorUint8Arr(d, Q[col], Q[col], (l + 7) / 8);
      }
    }
    free(fancy_d);
  }

  // Step: 11..12
  bf192_t* bf_q = column_to_row_major_and_shrink_V_192(Q, l);

  // Step: 13
  const unsigned int length_b = Ske + (beta * Senc) + 1;
  uint8_t* k                  = malloc((R + 1) * 128);
  bf192_t* vk                 = malloc(sizeof(bf192_t) * ((R + 1) * 128));
  bf192_t* qk                 = malloc(sizeof(bf192_t) * ((R + 1) * 128));
  bf192_t* B_0                = malloc(sizeof(bf192_t) * length_b);
  if (Lke > 0) {
    aes_key_schedule_constraints_192(NULL, NULL, 1, bf_q, delta, NULL, NULL, k, vk, B_0, qk,
                                     params);
  }

  // Step: 14
  bf192_t* B_1 = B_0 + Ske;
  aes_enc_constraints_192(in, out, NULL, NULL, NULL, NULL, 1, bf_q + Lke, qk, delta, NULL, NULL,
                          B_1, params);

  if (beta == 2) {
    // Step: 18
    bf192_t* B_2 = B_0 + (Ske + Senc);
    aes_enc_constraints_192(in + 16, out + 16, NULL, NULL, NULL, NULL, 1, bf_q + (Lke + Lenc), qk,
                            delta, NULL, NULL, B_2, params);
  }
  free(qk);
  free(vk);
  free(k);

  // Step: 20
  B_0[length_b - 1] = bf192_sum_poly(bf_q + l);
  free(bf_q);

  // Step 21
  uint8_t* q_tilde = malloc(lambdaBytes);
  zk_hash_192(q_tilde, chall_2, B_0, length_b - 1);
  free(B_0);

  bf192_t bf_qtilde = bf192_load(q_tilde);
  bf192_store(q_tilde, bf192_add(bf_qtilde, bf192_mul(bf192_load(a_tilde), bf192_load(delta))));

  return q_tilde;
}

// lambda == 256 implementation

static void aes_key_schedule_forward_256(uint32_t m, const uint8_t* x, const bf256_t* v,
                                         uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                                         uint8_t* out, bf256_t* bf_out,
                                         const faest_paramset_t* params) {
  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return;
  }

  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int R           = params->cipher_param.numRounds;
  const unsigned int Nwd         = params->faest_param.Nwd;
  const unsigned int lambdaBytes = lambda / 8;

  if (m == 1) {
    const unsigned int out_len = (R + 1) * 128 / 8;
    // Step 3
    memcpy(out, x, lambdaBytes);
    memset(out + lambdaBytes, 0, out_len - lambdaBytes);

    // Step: 4
    uint32_t i_wd = lambda;
    // Step: 5..10
    for (uint32_t j = Nwd; j < 4 * (R + 1); j++) {
      if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {
        memcpy(out + 32 * j / 8, x + i_wd / 8, 4);
        i_wd += 32;
      } else {
        for (uint32_t i = 0; i < 32; i += 8) {
          // bit spliced
          out[(32 * j + i) / 8] |= out[(32 * (j - Nwd) + i) / 8] ^ out[(32 * (j - 1) + i) / 8];
        }
      }
    }
    return;
  }

  // TODO: replace with one memcpy
  for (uint32_t i = 0; i < lambdaBytes; i++) {
    for (uint32_t j = 0; j < 8; j++) {
      bf_out[(i * 8) + j] = v[i * 8 + j];
    }
  }

  // Step: 4
  uint32_t i_wd = lambda;
  // Step: 5..10
  for (uint32_t j = Nwd; j < 4 * (R + 1); j++) {
    if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {
      for (uint32_t i = (j * 32); i <= ((j * 32) + 31); i++, i_wd++) {
        bf_out[i] = v[i_wd];
      }
    } else {
      for (uint32_t i = 0; i < 32; i++) {
        bf_out[(32 * j) + i] = bf256_add(bf_out[32 * (j - Nwd) + i], bf_out[32 * (j - 1) + i]);
      }
    }
  }
}

static void aes_key_schedule_backward_256(uint32_t m, const uint8_t* x, const bf256_t* v,
                                          const uint8_t* xk, const bf256_t* Vk, uint8_t Mtag,
                                          uint8_t Mkey, const uint8_t* delta, uint8_t* out,
                                          bf256_t* bf_out, const faest_paramset_t* params) {
  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return;
  }

  const unsigned int lambda = params->faest_param.lambda;
  const unsigned int Ske    = params->faest_param.Ske;

  // STep: 2
  if (m == 1) {
    assert(delta == NULL);

    uint32_t iwd   = 0;
    uint32_t c     = 0;
    bool rmvRcon   = true;
    uint32_t ircon = 0;

    for (uint32_t j = 0; j < Ske; j++) {
      // Step 7 (bit sliced)
      uint8_t x_tilde = x[j] ^ xk[(iwd + 8 * c) / 8];

      // Step 8
      if (Mtag == 0 && rmvRcon == true && c == 0) {
        uint8_t rcon = Rcon[ircon];
        ircon        = ircon + 1;
        // Steps 12 and 13, bitsliced; delta is always 0, Mkey is always 0
        x_tilde ^= rcon;
      }

      // Step: 15..19
      uint8_t y_tilde = 0;
      for (uint32_t i = 0; i < 8; ++i) {
        y_tilde ^= set_bit(get_bit(x_tilde, (i + 7) % 8) ^ get_bit(x_tilde, (i + 5) % 8) ^
                               get_bit(x_tilde, (i + 2) % 8),
                           i);
      }
      // delta is always 0
      y_tilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 0);
      y_tilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 2);

      // Step: 19
      out[j] = y_tilde;
      // Step: 20
      c = c + 1;

      if (c == 4) {
        c = 0;
        if (lambda == 192) {
          iwd += 192;
        } else {
          iwd += 128;
          if (lambda == 256) {
            rmvRcon = !rmvRcon;
          }
        }
      }
    }
    return;
  }

  bf256_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf256_zero();
  } else {
    bf_delta = bf256_load(delta);
  }

  uint32_t iwd   = 0;
  uint32_t c     = 0;
  bool rmvRcon   = true;
  uint32_t ircon = 0;

  bf256_t bf_mkey             = bf256_from_bit(Mkey);
  bf256_t bf_minus_mkey       = bf256_from_bit(1 ^ Mkey);
  bf256_t bf_minus_mtag       = bf256_from_bit(1 ^ Mtag);
  bf256_t bf_mkey_times_delta = bf256_mul(bf_mkey, bf_delta);
  bf_mkey_times_delta         = bf256_add(bf_mkey_times_delta, bf_minus_mkey);

  for (uint32_t j = 0; j < Ske; j++) {
    // Step 7
    bf256_t bf_x_tilde[8];
    for (uint32_t i = 0; i < 8; i++) {
      bf_x_tilde[i] = bf256_add(v[8 * j + i], Vk[iwd + 8 * c + i]);
    }

    if (Mtag == 0 && rmvRcon == true && c == 0) {
      // Step 9
      uint8_t r = Rcon[ircon];
      ircon     = ircon + 1;

      bf256_t bf_r[8];
      for (uint32_t i = 0; i < 8; i++) {
        // Step 12
        bf_r[i] = bf256_mul(bf256_from_bit(get_bit(r, i)), bf_mkey_times_delta);
        // Step 13
        bf_x_tilde[i] = bf256_add(bf_x_tilde[i], bf_r[i]);
      }
    }

    bf256_t bf_y_tilde[8];
    for (uint32_t i = 0; i < 8; ++i) {
      bf_y_tilde[i] = bf256_add(bf256_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                                bf_x_tilde[(i + 2) % 8]);
    }
    bf_y_tilde[0] = bf256_add(bf_y_tilde[0], bf256_mul(bf_minus_mtag, bf_mkey_times_delta));
    bf_y_tilde[2] = bf256_add(bf_y_tilde[2], bf256_mul(bf_minus_mtag, bf_mkey_times_delta));

    // TODO: get rid of this copy
    for (uint32_t i = 0; i < 8; i++) {
      bf_out[(8 * j) + i] = bf_y_tilde[i];
    }
    c = c + 1;

    if (c == 4) {
      c = 0;
      if (lambda == 192) {
        iwd += 192;
      } else {
        iwd += 128;
        if (lambda == 256) {
          rmvRcon = !rmvRcon;
        }
      }
    }
  }
}

static void aes_key_schedule_constraints_256(const uint8_t* w, const bf256_t* v, const uint8_t Mkey,
                                             const bf256_t* q, const uint8_t* delta, bf256_t* A0,
                                             bf256_t* A1, uint8_t* k, bf256_t* vk, bf256_t* B,
                                             bf256_t* qk, const faest_paramset_t* params) {
  const unsigned int lambda     = params->faest_param.lambda;
  const unsigned int Nwd        = params->faest_param.Nwd;
  const unsigned int Ske        = params->faest_param.Ske;
  const unsigned int lambdaByte = lambda / 8;

  bool rotate_word = true;

  if (Mkey == 0) {
    assert(delta == NULL);

    // Step: 2
    aes_key_schedule_forward_256(1, w, NULL, 0, 0, NULL, k, NULL, params);

    // Step: 3
    aes_key_schedule_forward_256(lambda, NULL, v, 1, 0, NULL, NULL, vk, params);

    // Step: 4
    uint8_t* w_dash = malloc(Ske);
    aes_key_schedule_backward_256(1, w + lambdaByte, NULL, k, NULL, 0, 0, NULL, w_dash, NULL,
                                  params);

    // Step: 5
    bf256_t* v_w_dash = malloc(Ske * 8 * sizeof(bf256_t));
    aes_key_schedule_backward_256(lambda, NULL, v + lambda, NULL, vk, 1, 0, NULL, NULL, v_w_dash,
                                  params);

    // Step: 6..8
    uint32_t iwd = 32 * (Nwd - 1);
    for (uint32_t j = 0; j < Ske / 4; j++) {
      bf256_t bf_k_hat[4];
      bf256_t bf_v_k_hat[4];
      bf256_t bf_w_dash_hat[4];
      bf256_t bf_v_w_dash_hat[4];
      for (uint32_t r = 0; r <= 3; r++) {
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
      for (uint32_t r = 0; r <= 3; r++) {
        // check that the constraint actually holds
        bf256_t should_be_zero = bf256_add(bf256_mul(bf_k_hat[r], bf_w_dash_hat[r]), bf256_one());
        assert(should_be_zero.values[0] == 0 && should_be_zero.values[1] == 0 &&
               should_be_zero.values[2] == 0 && should_be_zero.values[3] == 0);
        A0[4 * j + r] = bf256_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
        A1[4 * j + r] =
            bf256_add(bf256_add(bf256_mul(bf256_add(bf_k_hat[r], bf_v_k_hat[r]),
                                          bf256_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r])),
                                bf256_from_bf8(bf8_one())),
                      A0[4 * j + r]);
      }
      if (lambda == 192) {
        iwd = iwd + 192;
      } else {
        iwd = iwd + 128;
        if (lambda == 256) {
          rotate_word = !rotate_word;
        }
      }
    }
    free(v_w_dash);
    free(w_dash);
    return;
  }

  assert(delta);

  // Step: 19..20
  aes_key_schedule_forward_256(lambda, NULL, q, 0, 1, delta, NULL, qk, params);
  bf256_t* q_w_dash = malloc(Ske * 8 * sizeof(bf256_t));
  aes_key_schedule_backward_256(lambda, NULL, &q[lambda], NULL, qk, 0, 1, delta, NULL, q_w_dash,
                                params);

  const bf256_t bf_delta = bf256_load(delta);

  // Step 23..24
  uint32_t iwd = 32 * (Nwd - 1);
  for (uint32_t j = 0; j < Ske / 4; j++) {
    bf256_t bf_q_hat_k[4];
    bf256_t bf_q_hat_w_dash[4];
    for (uint32_t r = 0; r <= 3; r++) {
      // Step: 25..26
      if (rotate_word) {
        bf_q_hat_k[(r + 3) % 4] = bf256_byte_combine(qk + ((iwd + 8 * r)));
        bf_q_hat_w_dash[r]      = bf256_byte_combine(q_w_dash + ((32 * j + 8 * r)));
      } else {
        bf_q_hat_k[r]      = bf256_byte_combine(qk + ((iwd + 8 * r)));
        bf_q_hat_w_dash[r] = bf256_byte_combine(q_w_dash + ((32 * j + 8 * r)));
      }
    }
    // STep: 27
    for (uint32_t r = 0; r <= 3; r++) {
      bf256_t bf_tmp = bf256_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      B[4 * j + r]   = bf256_add(bf_tmp, bf256_mul(bf_delta, bf_delta));
    }
    if (lambda == 192) {
      iwd = iwd + 192;
    } else {
      iwd = iwd + 128;
      if (lambda == 256) {
        rotate_word = !rotate_word;
      }
    }
  }
  free(q_w_dash);
}

static void aes_enc_forward_256(uint32_t m, const uint8_t* x, const bf256_t* bf_x,
                                const uint8_t* xk, const bf256_t* bf_xk, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf256_t* bf_y,
                                const faest_paramset_t* params) {
  const unsigned int R = params->cipher_param.numRounds;

  bf256_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf256_zero();
  } else {
    bf_delta = bf256_load(delta);
  }

  if (m == 1) {
    // Step: 2
    for (uint32_t i = 0; i < 16; i++) {
      uint8_t xin = 0;
      // STep: 3
      for (uint32_t j = 0; j < 8; j++) {
        // STep: 4
        // TODO: check bit order
        xin |= set_bit(get_bit(in[i], j) & (1 ^ Mtag) & (1 ^ Mkey), j);
      }
      // STep: 5
      bf_y[i] = bf256_add(bf256_byte_combine_bits(xin), bf256_byte_combine_bits(xk[i]));
    }
    uint32_t ix, ik, iy;
    for (uint32_t j = 1; j < R; j++) {
      for (uint32_t c = 0; c <= 3; c++) {
        ix = 128 * (j - 1) + 32 * c;
        ik = 128 * j + 32 * c;
        iy = 16 * j + 4 * c;
        bf256_t bf_x_hat[4];
        bf256_t bf_xk_hat[4];
        for (uint32_t r = 0; r <= 3; r++) {
          // Step: 12..13
          bf_x_hat[r]  = bf256_byte_combine_bits(x[(ix + 8 * r) / 8]);
          bf_xk_hat[r] = bf256_byte_combine_bits(xk[(ik + 8 * r) / 8]);
        }
        bf256_t bf_one   = bf256_one();
        bf256_t bf_two   = bf256_byte_combine_bits(2);
        bf256_t bf_three = bf256_byte_combine_bits(3);
        // Step : 14
        bf_y[iy + 0] = bf256_add(bf_xk_hat[0], bf256_mul(bf_x_hat[0], bf_two));
        bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf256_mul(bf_x_hat[1], bf_three));
        bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf256_mul(bf_x_hat[2], bf_one));
        bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf256_mul(bf_x_hat[3], bf_one));

        // STep: 15
        bf_y[iy + 1] = bf256_add(bf_xk_hat[1], bf256_mul(bf_x_hat[0], bf_one));
        bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[1], bf_two));
        bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[2], bf_three));
        bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[3], bf_one));

        // STep: 16
        bf_y[iy + 2] = bf256_add(bf_xk_hat[2], bf256_mul(bf_x_hat[0], bf_one));
        bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[1], bf_one));
        bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[2], bf_two));
        bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[3], bf_three));

        // Step: 17
        bf_y[iy + 3] = bf256_add(bf_xk_hat[3], bf256_mul(bf_x_hat[0], bf_three));
        bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_x_hat[1], bf_one));
        bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_x_hat[2], bf_one));
        bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_x_hat[3], bf_two));
      }
    }
    return;
  }

  bf256_t bf_minus_mtag = bf256_from_bit(1 ^ Mtag);
  bf256_t bf_minus_mkey = bf256_from_bit(1 ^ Mkey);
  bf256_t bf_mkey       = bf256_from_bit(Mkey);

  // Step: 2..4
  for (uint32_t i = 0; i < 16; i++) {
    bf256_t bf_xin[8];
    for (uint32_t j = 0; j < 8; j++) {
      bf_xin[j] = bf256_mul(bf256_mul(bf256_from_bit(get_bit(in[i], j)), bf_minus_mtag),
                            bf256_add(bf256_mul(bf_mkey, bf_delta), bf_minus_mkey));
    }
    // Step: 5
    bf_y[i] = bf256_add(bf256_byte_combine(bf_xin), bf256_byte_combine(bf_xk + (8 * i)));
  }
  uint32_t ix, ik, iy;
  for (uint32_t j = 1; j < R; j++) {
    for (uint32_t c = 0; c <= 3; c++) {
      ix = 128 * (j - 1) + 32 * c;
      ik = 128 * j + 32 * c;
      iy = 16 * j + 4 * c;
      bf256_t bf_x_hat[4];
      bf256_t bf_xk_hat[4];
      for (uint32_t r = 0; r <= 3; r++) {
        // Step: 12..13
        bf_x_hat[r]  = bf256_byte_combine(bf_x + (ix + 8 * r));
        bf_xk_hat[r] = bf256_byte_combine(bf_xk + (ik + 8 * r));
      }
      bf256_t bf_one   = bf256_one();
      bf256_t bf_two   = bf256_byte_combine_bits(2);
      bf256_t bf_three = bf256_byte_combine_bits(3);
      // Step : 14
      bf_y[iy + 0] = bf256_add(bf_xk_hat[0], bf256_mul(bf_x_hat[0], bf_two));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf256_mul(bf_x_hat[1], bf_three));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf256_mul(bf_x_hat[2], bf_one));
      bf_y[iy + 0] = bf256_add(bf_y[iy + 0], bf256_mul(bf_x_hat[3], bf_one));

      // STep: 15
      bf_y[iy + 1] = bf256_add(bf_xk_hat[1], bf256_mul(bf_x_hat[0], bf_one));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[1], bf_two));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[2], bf_three));
      bf_y[iy + 1] = bf256_add(bf_y[iy + 1], bf256_mul(bf_x_hat[3], bf_one));

      // STep: 16
      bf_y[iy + 2] = bf256_add(bf_xk_hat[2], bf256_mul(bf_x_hat[0], bf_one));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[1], bf_one));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[2], bf_two));
      bf_y[iy + 2] = bf256_add(bf_y[iy + 2], bf256_mul(bf_x_hat[3], bf_three));

      // Step: 17
      bf_y[iy + 3] = bf256_add(bf_xk_hat[3], bf256_mul(bf_x_hat[0], bf_three));
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_x_hat[1], bf_one));
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_x_hat[2], bf_one));
      bf_y[iy + 3] = bf256_add(bf_y[iy + 3], bf256_mul(bf_x_hat[3], bf_two));
    }
  }
}

static void aes_enc_backward_256(uint32_t m, const uint8_t* x, const bf256_t* bf_x,
                                 const uint8_t* xk, const bf256_t* bf_xk, uint8_t Mtag,
                                 uint8_t Mkey, const uint8_t* delta, const uint8_t* out,
                                 bf256_t* y_out, const faest_paramset_t* params) {
  const unsigned int R = params->cipher_param.numRounds;

  if (m == 1) {
    uint8_t xtilde;
    // Step:2..4
    for (uint32_t j = 0; j < R; j++) {
      for (uint32_t c = 0; c <= 3; c++) {
        for (uint32_t r = 0; r <= 3; r++) {
          // Step: 5..6
          unsigned int ird = (128 * j) + (32 * ((c - r) % 4)) + (8 * r);
          if (j < (R - 1)) {
            // Step: 7
            xtilde = x[ird / 8];
          } else {
            // Step: 9
            uint8_t xout = 0;
            for (uint32_t i = 0; i < 8; i++) {
              // Step: 10..11
              // delta is always \bot if called with m == 1
              // TODO bit splice
              xout |=
                  set_bit(get_bit(out[(ird - 128 * (R - 1)) / 8], i) & (1 ^ Mtag) & (1 ^ Mkey), i);
            }
            xtilde = xout ^ xk[(128 + ird) / 8];
          }
          // Step: 12..17
          uint8_t ytilde = 0;
          for (uint32_t i = 0; i < 8; ++i) {
            ytilde ^= set_bit(get_bit(xtilde, (i + 7) % 8) ^ get_bit(xtilde, (i + 5) % 8) ^
                                  get_bit(xtilde, (i + 2) % 8),
                              i);
          }
          ytilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 0);
          ytilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 2);

          // Step: 18
          y_out[16 * j + 4 * c + r] = bf256_byte_combine_bits(ytilde);
        }
      }
    }
    return;
  }

  // Step: 1
  bf256_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf256_zero();
  } else {
    bf_delta = bf256_load(delta);
  }

  const bf256_t factor =
      bf256_mul(bf256_from_bit(1 ^ Mtag),
                bf256_add(bf256_mul(bf256_from_bit(Mkey), bf_delta), bf256_from_bit(1 ^ Mkey)));

  // STep: 2..4
  for (uint32_t j = 0; j < R; j++) {
    for (uint32_t c = 0; c <= 3; c++) {
      for (uint32_t r = 0; r <= 3; r++) {
        bf256_t bf_x_tilde[8];
        // Step: 5
        unsigned int ird = (128 * j) + (32 * ((c - r) % 4)) + (8 * r);
        // Step: 6
        if (j < (R - 1)) {
          // Step: 7
          for (uint32_t i = 0; i < 8; i++) {
            bf_x_tilde[i] = bf_x[ird + i];
          }
        } else {
          bf256_t bf_xout[8];
          // Step: 10
          for (uint32_t i = 0; i < 8; ++i) {
            // Step: 11
            // TODO: check 0 extension
            bf_xout[i] = bf256_from_bit(get_bit(out[(ird - 128 * (R - 1)) / 8], i));
            bf_xout[i] = bf256_mul(bf_xout[i], factor);
            // Step: 12
            bf_x_tilde[i] = bf256_add(bf_xout[i], bf_xk[128 + ird + i]);
          }
        }
        // STep: 13..17
        bf256_t bf_y_tilde[8];
        for (uint32_t i = 0; i < 8; ++i) {
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
                                    const uint8_t* delta, bf256_t* A0, bf256_t* A1, bf256_t* B,
                                    const faest_paramset_t* params) {
  const unsigned int lambda = params->faest_param.lambda;
  const unsigned int Senc   = params->faest_param.Senc;

  if (Mkey == 0) {
    bf256_t* s       = malloc(sizeof(bf256_t) * Senc);
    bf256_t* vs      = malloc(sizeof(bf256_t) * Senc);
    bf256_t* s_dash  = malloc(sizeof(bf256_t) * Senc);
    bf256_t* vs_dash = malloc(sizeof(bf256_t) * Senc);
    aes_enc_forward_256(1, w, NULL, k, NULL, in, 0, 0, NULL, s, params);
    aes_enc_forward_256(lambda, NULL, v, NULL, vk, in, 1, 0, NULL, vs, params);
    aes_enc_backward_256(1, w, NULL, k, NULL, 0, 0, NULL, out, s_dash, params);
    aes_enc_backward_256(lambda, NULL, v, NULL, vk, 1, 0, NULL, out, vs_dash, params);

    for (uint32_t j = 0; j < Senc; j++) {
      // check that the constraint actually holds
      bf256_t should_be_zero = bf256_add(bf256_mul(s[j], s_dash[j]), bf256_one());
      assert(should_be_zero.values[0] == 0 && should_be_zero.values[1] == 0 &&
             should_be_zero.values[2] == 0 && should_be_zero.values[3] == 0);
      A0[j] = bf256_mul(vs[j], vs_dash[j]);
      A1[j] = bf256_add(
          bf256_add(bf256_mul(bf256_add(s[j], vs[j]), bf256_add(s_dash[j], vs_dash[j])), A0[j]),
          bf256_one());
    }

    free(vs_dash);
    free(s_dash);
    free(vs);
    free(s);
  } else {
    // Step: 11..12
    bf256_t* qs      = malloc(sizeof(bf256_t) * Senc);
    bf256_t* qs_dash = malloc(sizeof(bf256_t) * Senc);
    aes_enc_forward_256(lambda, NULL, q, NULL, qk, in, 0, 1, delta, qs, params);
    aes_enc_backward_256(lambda, NULL, q, NULL, qk, 0, 1, delta, out, qs_dash, params);

    // Step: 13..14
    bf256_t minus_part = bf256_mul(bf256_load(delta), bf256_load(delta));
    for (uint32_t j = 0; j < Senc; j++) {
      B[j] = bf256_add(bf256_mul(qs[j], qs_dash[j]), minus_part);
    }
    free(qs);
    free(qs_dash);
  }
}

static void aes_prove_256(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
                          const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde,
                          uint8_t* b_tilde, const faest_paramset_t* params) {
  const unsigned int beta = params->faest_param.beta;
  const unsigned int l    = params->faest_param.l;
  const unsigned int Lke  = params->faest_param.Lke;
  const unsigned int Lenc = params->faest_param.Lenc;
  const unsigned int R    = params->cipher_param.numRounds;
  const unsigned int Ske  = params->faest_param.Ske;
  const unsigned int Senc = params->faest_param.Senc;

  // Step: 1..2
  bf256_t* bf_v = column_to_row_major_and_shrink_V_256(V, l);

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7
  const unsigned int length_a = Ske + (beta * Senc) + 1;
  bf256_t* A0                 = malloc(sizeof(bf256_t) * length_a);
  bf256_t* A1                 = malloc(sizeof(bf256_t) * length_a);
  uint8_t* k                  = malloc((R + 1) * 128 / 8);
  bf256_t* vk                 = malloc(sizeof(bf256_t) * ((R + 1) * 128));
  bf256_t* qk                 = malloc(sizeof(bf256_t) * ((R + 1) * 128));
  if (Lke > 0) {
    aes_key_schedule_constraints_256(w, bf_v, 0, NULL, NULL, A0, A1, k, vk, NULL, qk, params);
  }

  // Step: Skipping 8 in implementation
  // Step: 9

  // Step: 10,11
  aes_enc_constraints_256(in, out, w + Lke / 8, bf_v + Lke, k, vk, 0, NULL, NULL, NULL, A0 + Ske,
                          A1 + Ske, NULL, params);
  // Step: 12
  if (beta == 2) {
    // Step: 13
    // Step: 14, 15
    aes_enc_constraints_256(in + 16, out + 16, w + (Lke + Lenc) / 8, bf_v + Lke + Lenc, k, vk, 0,
                            NULL, NULL, NULL, A0 + (Ske + Senc), A1 + (Ske + Senc), NULL, params);
  }
  free(qk);
  free(vk);
  free(k);

  // Step: 16..18
  A1[length_a - 1] = bf256_load(u + l / 8);
  A0[length_a - 1] = bf256_sum_poly(bf_v + l);
  free(bf_v);

  zk_hash_256(a_tilde, chall, A1, length_a - 1);
  zk_hash_256(b_tilde, chall, A0, length_a - 1);

  free(A0);
  free(A1);
}

static uint8_t* aes_verify_256(uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                               const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* in,
                               const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int tau         = params->faest_param.tau;
  const unsigned int t0          = params->faest_param.t0;
  const unsigned int k0          = params->faest_param.k0;
  const unsigned int t1          = params->faest_param.t1;
  const unsigned int k1          = params->faest_param.k1;
  const unsigned int beta        = params->faest_param.beta;
  const unsigned int l           = params->faest_param.l;
  const unsigned int Lke         = params->faest_param.Lke;
  const unsigned int Lenc        = params->faest_param.Lenc;
  const unsigned int R           = params->cipher_param.numRounds;
  const unsigned int Ske         = params->faest_param.Ske;
  const unsigned int Senc        = params->faest_param.Senc;
  const unsigned int lambdaBytes = lambda / 8;

  // Step: 1
  const uint8_t* delta = chall_3;
  // Step: 2,3
  // do nothing

  // Step: 4..10
  for (uint32_t i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t* fancy_d   = malloc(depth);
    ChalDec(chall_3, i, k0, t0, k1, t1, fancy_d);
    for (uint32_t j = 0; j < depth; j++, ++col) {
      if (fancy_d[j] == 1) {
        xorUint8Arr(d, Q[col], Q[col], (l + 7) / 8);
      }
    }
    free(fancy_d);
  }

  // Step: 11..12
  bf256_t* bf_q = column_to_row_major_and_shrink_V_256(Q, l);

  // Step: 13
  const unsigned int length_b = Ske + (beta * Senc) + 1;
  uint8_t* k                  = malloc((R + 1) * 128);
  bf256_t* vk                 = malloc(sizeof(bf256_t) * ((R + 1) * 128));
  bf256_t* qk                 = malloc(sizeof(bf256_t) * ((R + 1) * 128));
  bf256_t* B_0                = malloc(sizeof(bf256_t) * length_b);
  if (Lke > 0) {
    aes_key_schedule_constraints_256(NULL, NULL, 1, bf_q, delta, NULL, NULL, k, vk, B_0, qk,
                                     params);
  }

  // Step: 14
  bf256_t* B_1 = B_0 + Ske;
  aes_enc_constraints_256(in, out, NULL, NULL, NULL, NULL, 1, bf_q + Lke, qk, delta, NULL, NULL,
                          B_1, params);

  if (beta == 2) {
    // Step: 18
    bf256_t* B_2 = B_0 + (Ske + Senc);
    aes_enc_constraints_256(in + 16, out + 16, NULL, NULL, NULL, NULL, 1, bf_q + (Lke + Lenc), qk,
                            delta, NULL, NULL, B_2, params);
  }
  free(qk);
  free(vk);
  free(k);

  // Step: 20
  B_0[length_b - 1] = bf256_sum_poly(bf_q + l);
  free(bf_q);

  // Step 21
  uint8_t* q_tilde = malloc(lambdaBytes);
  zk_hash_256(q_tilde, chall_2, B_0, length_b - 1);
  free(B_0);

  bf256_t bf_qtilde = bf256_load(q_tilde);
  bf256_store(q_tilde, bf256_add(bf_qtilde, bf256_mul(bf256_load(a_tilde), bf256_load(delta))));

  return q_tilde;
}

// EM-128

static void em_enc_forward_128(uint32_t m, const uint8_t* z, const uint8_t* x, const bf128_t* bf_x,
                               const uint8_t* xk, const bf128_t* bf_xk, const uint8_t* in,
                               uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, bf128_t* bf_y,
                               const faest_paramset_t* params) {

  const unsigned int R = params->cipher_param.numRounds;

  bf128_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf128_zero();
  } else {
    bf_delta = bf128_load(delta);
  }

  if (m == 1) {
    // Step: 2
    for (uint32_t j = 0; j < 16; j++) {
      bf_y[j] = bf128_add(bf128_byte_combine_bits(z[j]), bf128_byte_combine_bits(x[j]));
    }

    uint32_t i, iy;

    for (uint32_t j = 1; j < R; j++) {

      for (uint32_t c = 0; c <= 3; c++) {

        i  = 128 * j + 32 * c;
        iy = 16 * j + 4 * c;

        bf128_t bf_x_hat[4];
        bf128_t bf_z_hat[4];

        for (uint32_t r = 0; r <= 3; r++) {

          // Step: 12..13
          bf_z_hat[r] = bf128_byte_combine_bits(z[(i + 8 * r) / 8]);
          bf_x_hat[r] = bf128_byte_combine_bits(x[(i + 8 * r) / 8]);
        }

        bf128_t bf_one   = bf128_one();
        bf128_t bf_two   = bf128_byte_combine_bits(2);
        bf128_t bf_three = bf128_byte_combine_bits(3);

        bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_two), bf128_mul(bf_z_hat[1], bf_three));
        bf_y[iy + 1] = bf128_add(bf_y[iy + 0], bf128_mul(bf_z_hat[2], bf_one));
        bf_y[iy + 2] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[3], bf_one));
        bf_y[iy + 3] = bf128_add(bf_y[iy + 0], bf_x_hat[0]);

        bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_one), bf128_mul(bf_z_hat[1], bf_two));
        bf_y[iy + 1] = bf128_add(bf_y[iy + 0], bf128_mul(bf_z_hat[2], bf_three));
        bf_y[iy + 2] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[3], bf_one));
        bf_y[iy + 3] = bf128_add(bf_y[iy + 0], bf_x_hat[1]);

        bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_one), bf128_mul(bf_z_hat[1], bf_one));
        bf_y[iy + 1] = bf128_add(bf_y[iy + 0], bf128_mul(bf_z_hat[2], bf_two));
        bf_y[iy + 2] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[3], bf_three));
        bf_y[iy + 3] = bf128_add(bf_y[iy + 0], bf_x_hat[2]);

        bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_three), bf128_mul(bf_z_hat[1], bf_one));
        bf_y[iy + 1] = bf128_add(bf_y[iy + 0], bf128_mul(bf_z_hat[2], bf_one));
        bf_y[iy + 2] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[3], bf_two));
        bf_y[iy + 3] = bf128_add(bf_y[iy + 0], bf_x_hat[3]);
      }
    }
    return;
  }

  bf128_t bf_minus_mtag = bf128_from_bit(1 ^ Mtag);
  bf128_t bf_minus_mkey = bf128_from_bit(1 ^ Mkey);
  bf128_t bf_mkey       = bf128_from_bit(Mkey);

  // Step: 2
  for (uint32_t j = 0; j < 16; j++) {
    bf_y[j] = bf128_add(bf128_byte_combine_bits(z[j]), bf128_byte_combine_bits(x[j]));
  }

  uint32_t i, iy;

  for (uint32_t j = 1; j < R; j++) {

    for (uint32_t c = 0; c <= 3; c++) {

      i  = 128 * j + 32 * c;
      iy = 16 * j + 4 * c;

      bf128_t bf_x_hat[4];
      bf128_t bf_z_hat[4];

      for (uint32_t r = 0; r <= 3; r++) {

        // Step: 12..13
        bf_z_hat[r] = bf128_byte_combine(z[(i + 8 * r) / 8]);
        bf_x_hat[r] = bf128_byte_combine(x[(i + 8 * r) / 8]);
      }

      bf128_t bf_one   = bf128_one();
      bf128_t bf_two   = bf128_byte_combine_bits(2);
      bf128_t bf_three = bf128_byte_combine_bits(3);

      bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_two), bf128_mul(bf_z_hat[1], bf_three));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 0], bf128_mul(bf_z_hat[2], bf_one));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[3], bf_one));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 0], bf_x_hat[0]);

      bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_one), bf128_mul(bf_z_hat[1], bf_two));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 0], bf128_mul(bf_z_hat[2], bf_three));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[3], bf_one));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 0], bf_x_hat[1]);

      bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_one), bf128_mul(bf_z_hat[1], bf_one));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 0], bf128_mul(bf_z_hat[2], bf_two));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[3], bf_three));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 0], bf_x_hat[2]);

      bf_y[iy + 0] = bf128_add(bf128_mul(bf_z_hat[0], bf_three), bf128_mul(bf_z_hat[1], bf_one));
      bf_y[iy + 1] = bf128_add(bf_y[iy + 0], bf128_mul(bf_z_hat[2], bf_one));
      bf_y[iy + 2] = bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[3], bf_two));
      bf_y[iy + 3] = bf128_add(bf_y[iy + 0], bf_x_hat[3]);
    }
  }
}

static void em_enc_backward_128(uint32_t m, const uint8_t* z, const bf128_t* bf_z, const uint8_t* x,
                                const bf128_t* bf_x, uint8_t Mtag, uint8_t Mkey,
                                const uint8_t* delta, const uint8_t* z_out, bf128_t* y_out,
                                const faest_paramset_t* params) {
  const unsigned int R = params->cipher_param.numRounds;

  if (m == 1) {
    uint8_t z_tilde;

    for (uint32_t j = 0; j < R; j++) {

      for (uint32_t c = 0; c <= 3; c++) {

        for (uint32_t r = 0; r <= 3; r++) {

          unsigned int ird = (128 * j) + (32 * ((c - r) % 4)) + (8 * r);

          if (j < (R - 1)) {

            z_tilde = z[ird / 8];

          } else {

            uint8_t z_tilde_out = 0;

            for (uint32_t i = 0; i < 8; i++) {

              // delta is always \bot if called with m == 1
              // TODO bit splice
              z_tilde_out |= set_bit(get_bit(z_out[(ird - 128 * (R - 1)) / 8], i), i);
            }
            z_tilde = z_tilde_out ^ x[(128 + ird) / 8];
          }

          uint8_t y_tilde = 0;

          for (uint32_t i = 0; i < 8; ++i) {
            y_tilde ^= set_bit(get_bit(z_tilde, (i + 7) % 8) ^ get_bit(z_tilde, (i + 5) % 8) ^
                                   get_bit(z_tilde, (i + 2) % 8),
                               i);
          }
          y_tilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 0);
          y_tilde ^= set_bit(((1 ^ Mtag) & (1 ^ Mkey)), 2);

          // Step: 18
          y_out[16 * j + 4 * c + r] = bf128_byte_combine_bits(y_tilde);
        }
      }
    }
    return;
  }

  // Step: 1
  bf128_t bf_delta;
  if (delta == NULL) {
    bf_delta = bf128_zero();
  } else {
    bf_delta = bf128_load(delta);
  }

  const bf128_t factor =
      bf128_mul(bf128_from_bit(1 ^ Mtag),
                bf128_add(bf128_mul(bf128_from_bit(Mkey), bf_delta), bf128_from_bit(1 ^ Mkey)));

  for (uint32_t j = 0; j < R; j++) {

    for (uint32_t c = 0; c <= 3; c++) {

      for (uint32_t r = 0; r <= 3; r++) {

        bf128_t bf_z_tilde[8];

        unsigned int ird = (128 * j) + (32 * ((c - r) % 4)) + (8 * r);

        if (j < (R - 1)) {

          for (uint32_t i = 0; i < 8; i++) {
            bf_z_tilde[i] = bf_z[ird + i];
          }

        } else {
          bf128_t bf_z_tilde_out[8];

          for (uint32_t i = 0; i < 8; ++i) {

            bf_z_tilde_out[i] = bf128_from_bit(get_bit(z_out[(ird - 128 * (R - 1)) / 8], i));
            bf_z_tilde_out[i] = bf128_mul(bf_z_tilde_out[i], factor);
            // Step: 12
            bf_z_tilde[i] = bf128_add(bf_z_tilde_out[i], bf_x[128 + ird + i]);
          }
        }

        bf128_t bf_y_tilde[8];

        for (uint32_t i = 0; i < 8; ++i) {
          bf_y_tilde[i] = bf128_add(bf128_add(bf_z_tilde[(i + 7) % 8], bf_z_tilde[(i + 5) % 8]),
                                    bf_z_tilde[(i + 2) % 8]);
        }
        bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
        bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

        // Step: 18
        y_out[16 * j + 4 * c + r] = bf128_byte_combine(bf_y_tilde);
      }
    }
  }
}

// dispatchers

void aes_prove(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
               const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params) {
  switch (params->faest_param.lambda) {
  case 256:
    aes_prove_256(w, u, V, in, out, chall, a_tilde, b_tilde, params);
    break;
  case 192:
    aes_prove_192(w, u, V, in, out, chall, a_tilde, b_tilde, params);
    break;
  default:
    aes_prove_128(w, u, V, in, out, chall, a_tilde, b_tilde, params);
  }
}

uint8_t* aes_verify(uint8_t* d, uint8_t** Q, const uint8_t* chall_2, const uint8_t* chall_3,
                    const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out,
                    const faest_paramset_t* params) {
  switch (params->faest_param.lambda) {
  case 256:
    return aes_verify_256(d, Q, chall_2, chall_3, a_tilde, in, out, params);
  case 192:
    return aes_verify_192(d, Q, chall_2, chall_3, a_tilde, in, out, params);
  default:
    return aes_verify_128(d, Q, chall_2, chall_3, a_tilde, in, out, params);
  }
}
