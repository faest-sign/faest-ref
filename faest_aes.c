#include "faest_aes.h"
// TODO: Do not pass lambdaBytes everywhere, compute it in the function....
// TODO: change q to Q where applicable

static uint8_t Rcon[10] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x26};

// TODO: generalize bf128_t to bf(lambda)_t
int aes_key_schedule_forward(uint32_t lambda, uint32_t R, uint32_t Nwd, uint32_t Lke, uint32_t m,
                             const uint8_t* x, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                             uint8_t* y_out) {

  uint32_t lambdaBytes = lambda / 8;
  uint32_t x_len       = (Lke / 8) * ((m + 7) / 8);
  uint32_t delta_len   = lambdaBytes;
  uint32_t xk_len      = (Lke / 8) * ((m + 7) / 8);
  uint32_t y_out_len   = (R + 1) * 128 * ((m + 7) / 8);

  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return 0;
  }

  // STep: 2..3
  if (m == 1) {
    bf8_t* bf_y = malloc(y_out_len);
    for (uint32_t i = 0; i < lambdaBytes; i++) {
      for (uint32_t j = 0; j < 8; j++) {
        bf_y[(i * 8) + j] = bf8_from_bit(getBit(x + i, j));
      }
    }
    // Step: 4
    uint32_t i_wd = lambda;

    // Step: 5..10
    for (uint32_t j = Nwd; j < 4 * (R + 1); j++) {

      if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {

        for (uint32_t i = (j * 32); i <= ((j * 32) + 31); i++) {
          // TODO: hopefully this magic works fine here
          bf_y[i] = bf8_from_bit(getBit(x + (i / 8), i % 8));
        }
        i_wd += 32;
      } else {
        for (uint32_t i = 0; i < 32; i++) {
          bf_y[(32 * j) + i] = bf8_add(bf_y[32 * (j - Nwd) + i], bf_y[32 * (j - 1) + i]);
        }
      }
    }

    y_out = malloc(y_out_len);
    for (uint32_t i = 0; i < y_out_len; i++) {
      bf8_store(y_out + i, bf_y[i]);
    }
  } else {
    bf128_t* bf_y = malloc(y_out_len);
    for (uint32_t i = 0; i < lambdaBytes; i++) {
      for (uint32_t j = 0; j < 8; j++) {
        bf_y[(i * 8) + j] = bf128_from_bit(getBit(x + i, j));
      }
    }
    // Step: 4
    uint32_t i_wd = lambda;

    // Step: 5..10
    for (uint32_t j = Nwd; j < 4 * (R + 1); j++) {

      if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {

        for (uint32_t i = (j * 32); i <= ((j * 32) + 31); i++) {
          // TODO: hopefully this magic works fine here
          bf_y[i] = bf128_from_bit(getBit(x + (i / 8), i % 8));
        }
        i_wd += 32;
      } else {
        for (uint32_t i = 0; i < 32; i++) {
          bf_y[(32 * j) + i] = bf128_add(bf_y[32 * (j - Nwd) + i], bf_y[32 * (j - 1) + i]);
        }
      }
    }
    y_out        = malloc(y_out_len);
    uint32_t idx = 0;
    for (uint32_t i = 0; i < y_out_len; i += lambdaBytes) {
      bf128_store(y_out + i, bf_y[idx]);
      idx += 1;
    }
  }
}

int aes_key_schedule_backward(uint32_t lambda, uint32_t R, uint32_t Nwd, uint32_t Ske, uint8_t Lke,
                              uint32_t m, const uint8_t* x, const uint8_t* xk, uint8_t Mtag,
                              uint8_t Mkey, const uint8_t* delta, uint8_t* y_out) {

  uint32_t lambdaBytes = lambda / 8;
  uint32_t x_len       = (Lke / 8) * ((m + 7) / 8);
  uint32_t delta_len   = lambdaBytes;
  uint32_t xk_len      = (Lke / 8) * ((m + 7) / 8);
  uint32_t y_out_len   = 8 * Ske * ((m + 7) / 8);

  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return 0;
  }

  bf128_t bf_delta = bf128_load(delta);

  // STep: 2
  if (m == 1) {

    // Mkey should be 0, else we hav ethe multiplication with delta problem
    if (Mkey != 0) {
      return -1;
    }

    bf8_t* bf_y    = malloc(8 * Ske);
    uint32_t iwd   = 0;
    uint32_t c     = 0;
    bool rmvRcon   = true;
    uint32_t ircon = 0;

    bf8_t bf_mkey       = bf8_from_bit(Mkey);
    bf8_t bf_minus_mkey = bf8_from_bit(1 - Mkey);
    bf8_t bf_minus_mtag = bf8_from_bit(1 - Mtag);

    bf8_t* bf_x_tilde = malloc(8);

    for (uint32_t j = 0; j < Ske; j++) {

      for (uint32_t i = 0; i < 8; i++) {
        bf_x_tilde[i] = bf8_add(bf8_from_bit(getBit(x + j, i)),
                                bf8_from_bit(getBit(xk + ((iwd + 8 * c) / 8), i)));
      }

      if (Mtag == 0 && rmvRcon == true && c == 0) {
        // TODO: different r's !!!!
        uint8_t r   = Rcon[ircon];
        ircon       = ircon + 1;
        bf8_t* bf_r = malloc(8);
        for (uint32_t i = 0; i < 8; i++) {
          bf_r[i]       = bf8_from_bit(getBit(r, i));
          bf_r[i]       = bf8_mul(bf_r[i], bf_minus_mkey);
          bf_x_tilde[i] = bf8_add(bf_x_tilde[i], bf_r[i]);
        }
      }

      bf8_t* bf_y_tilde = malloc(8);

      bf_y_tilde[7] = bf8_add(bf8_add(bf8_add(bf_x_tilde[5], bf_x_tilde[2]), bf_x_tilde[0]),
                              bf8_mul(bf_minus_mtag, bf_minus_mkey));
      bf_y_tilde[6] = bf8_add(bf8_add(bf_x_tilde[7], bf_x_tilde[4]), bf_x_tilde[1]);
      bf_y_tilde[5] = bf8_add(bf8_add(bf8_add(bf_x_tilde[6], bf_x_tilde[3]), bf_x_tilde[0]),
                              bf8_mul(bf_minus_mtag, bf_minus_mkey));
      bf_y_tilde[4] = bf8_add(bf8_add(bf_x_tilde[7], bf_x_tilde[5]), bf_x_tilde[2]);
      bf_y_tilde[3] = bf8_add(bf8_add(bf_x_tilde[6], bf_x_tilde[4]), bf_x_tilde[1]);
      bf_y_tilde[2] = bf8_add(bf8_add(bf_x_tilde[5], bf_x_tilde[3]), bf_x_tilde[0]);
      bf_y_tilde[1] = bf8_add(bf8_add(bf_x_tilde[7], bf_x_tilde[4]), bf_x_tilde[2]);
      bf_y_tilde[0] = bf8_add(bf8_add(bf_x_tilde[6], bf_x_tilde[3]), bf_x_tilde[1]);

      for (uint32_t i = 0; i < 8; i++) {
        bf_y[(8 * j) + i] = bf_y_tilde[i];
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
    y_out = malloc(y_out_len);
    for (uint32_t i = 0; i < y_out_len; i++) {
      bf8_store(y_out + i, bf_y[i]);
    }
  } else {
    bf128_t* bf_y  = malloc(sizeof(bf128_t) * 8 * Ske);
    uint32_t iwd   = 0;
    uint32_t c     = 0;
    bool rmvRcon   = true;
    uint32_t ircon = 0;

    bf128_t bf_mkey       = bf128_from_bit(Mkey);
    bf128_t bf_minus_mkey = bf128_from_bit(1 - Mkey);
    bf128_t bf_minus_mtag = bf128_from_bit(1 - Mtag);

    bf128_t* bf_x_tilde = malloc(sizeof(bf128_t) * 8);

    for (uint32_t j = 0; j < Ske; j++) {

      for (uint32_t i = 0; i < 8; i++) {
        bf_x_tilde[i] = bf128_add(bf128_from_bit(getBit(x + j, i)),
                                  bf128_from_bit(getBit(xk + ((iwd + 8 * c) / 8), i)));
      }

      if (Mtag == 0 && rmvRcon == true && c == 0) {
        uint8_t r     = Rcon[ircon];
        ircon         = ircon + 1;
        bf128_t* bf_r = malloc(sizeof(bf128_t) * 8);

        for (uint32_t i = 0; i < 8; i++) {
          // TODO: hopefully the endianness is fine here
          bf_r[i]                   = bf128_from_bit(getBit(r, i));
          bf128_t bf_zeros_r_concat = bf128_zero();
          bf_zeros_r_concat         = bf128_add(bf_r[i], bf_zeros_r_concat);
          bf_r[i]                   = bf128_mul(bf_zeros_r_concat,
                                                bf128_add(bf128_mul(bf_mkey, bf128_load(delta)), bf_minus_mkey));
          bf_x_tilde[i]             = bf128_add(bf_x_tilde[i], bf_r[i]);
        }
      }

      bf128_t* bf_y_tilde = malloc(sizeof(bf128_t) * 8);

      bf_y_tilde[7] = bf128_add(bf128_add(bf128_add(bf_x_tilde[5], bf_x_tilde[2]), bf_x_tilde[0]),
                                bf128_mul(bf_minus_mtag, bf_minus_mkey));
      bf_y_tilde[6] = bf128_add(bf128_add(bf_x_tilde[7], bf_x_tilde[4]), bf_x_tilde[1]);
      bf_y_tilde[5] = bf128_add(bf128_add(bf128_add(bf_x_tilde[6], bf_x_tilde[3]), bf_x_tilde[0]),
                                bf128_mul(bf_minus_mtag, bf_minus_mkey));
      bf_y_tilde[4] = bf128_add(bf128_add(bf_x_tilde[7], bf_x_tilde[5]), bf_x_tilde[2]);
      bf_y_tilde[3] = bf128_add(bf128_add(bf_x_tilde[6], bf_x_tilde[4]), bf_x_tilde[1]);
      bf_y_tilde[2] = bf128_add(bf128_add(bf_x_tilde[5], bf_x_tilde[3]), bf_x_tilde[0]);
      bf_y_tilde[1] = bf128_add(bf128_add(bf_x_tilde[7], bf_x_tilde[4]), bf_x_tilde[2]);
      bf_y_tilde[0] = bf128_add(bf128_add(bf_x_tilde[6], bf_x_tilde[3]), bf_x_tilde[1]);

      for (uint32_t i = 0; i < 8; i++) {
        bf_y[(8 * j) + i] = bf_y_tilde[i];
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
    y_out        = malloc(y_out_len);
    uint32_t idx = 0;
    for (uint32_t i = 0; i < y_out_len; i += lambdaBytes) {
      bf128_store(y_out + i, bf_y[idx]);
      idx += 1;
    }
  }
  return 1;
}

void aes_key_schedule_constraints(uint32_t lambda, uint32_t R, uint32_t Nwd, uint32_t Ske,
                                  uint32_t Lke, const uint8_t* w, const uint8_t* v,
                                  const uint8_t Mkey, const uint8_t* q, const uint8_t* delta,
                                  uint8_t* A0, uint8_t* A1, uint8_t* k, uint8_t* vk, uint8_t* B,
                                  uint8_t* qk) {
  uint32_t lambdaByte = lambda / 8;
  uint32_t w_len      = Lke / 8;
  uint32_t v_len      = lambdaByte * (Lke / 8);
  uint32_t q_len      = lambdaByte * (Lke / 8);
  uint32_t delta_len  = lambdaByte;

  if (Mkey == 0) {
    // If m == 1, Mkey should be 0, because we then cancel delta in line 10
    if (Mkey != 0) {
      return -1;
    }

    // STep: 2
    uint8_t *k, vk, w_dash, v_w_dash;
    aes_key_schedule_forward(lambda, R, Nwd, Lke, 1, w, 0, 0, NULL, k);

    // Step: 3
    aes_key_schedule_forward(lambda, R, Nwd, Lke, lambda, v, 1, 0, NULL, vk);

    // Step: 4
    uint8_t* w_lambda = malloc(w_len - lambdaByte);
    memcpy(w_lambda, w + lambdaByte, w_len - lambdaByte);
    aes_key_schedule_backward(lambda, R, Nwd, Ske, Lke, 1, w_lambda, k, 0, 0, NULL, w_dash);

    // Step: 5
    uint8_t* v_lambda = malloc(sizeof(v) - lambdaByte);
    memcpy(v_lambda, v + lambdaByte, sizeof(v) - lambdaByte);
    aes_key_schedule_backward(lambda, R, Nwd, Ske, Lke, lambda, v_lambda, vk, 1, 0, NULL, v_w_dash);

    // TODO: correctly initilaized ??
    A0 = malloc(lambdaByte * (Ske / 8));
    A1 = malloc(lambdaByte * (Ske / 8));

    // Step: 6..8
    uint32_t iwd = 32 * (Nwd - 1);
    for (uint32_t j = 0; j < Ske / 4; j++) {

      bf128_t* bf_k_hat        = malloc(sizeof(bf128_t) * 4);
      bf128_t* bf_v_k_hat      = malloc(sizeof(bf128_t) * 4);
      bf128_t* bf_w_dash_hat   = malloc(sizeof(bf128_t) * 4);
      bf128_t* bf_v_w_dash_hat = malloc(sizeof(bf128_t) * 4);
      for (uint32_t r = 0; r <= 3; r) {
        // TODO: ByteCombine
        bf_k_hat[(r + 3) % 4];
        bf_k_hat[(r + 3) % 4];
      }
      // Step: 13..15
      for (uint32_t r = 0; r <= 3; r) {
        // TODO: check
        // TODO: unsure what is happeninig here with A0, and A1
        bf128_mul(bf_v_k_hat[r], bf_v_w_dash_hat[r]);
        // TODO: xor with 1_{F_{2^8}} ??
        bf128_mul(bf128_add(bf_k_hat[r], bf_v_k_hat[r]),
                  bf128_add(bf_w_dash_hat[r], bf_v_w_dash_hat[r]));
      }
      if (lambda == 192) {
        iwd = iwd + 192;
      } else {
        iwd = iwd + 128;
      }
    }
  } else {
    // Step: 19..20
    aes_key_schedule_forward(lambda, R, Nwd, Lke, lambda, q, 0, 1, delta, qk);
    uint8_t* q_w_dash;
    uint8_t* q_lambda = malloc(q_len - lambdaByte);
    memcpy(q_lambda, q + lambdaByte, q_len - lambdaByte);
    aes_key_schedule_backward(lambda, R, Nwd, Ske, Lke, lambda, q_lambda, qk, 0, 1, delta,
                              q_w_dash);

    B = malloc(lambdaByte * (Ske / 8));
    // Step 20..22
    uint32_t iwd = 32 * (Nwd - 1);
    for (uint32_t j = 0; j < Ske / 4; j++) {

      bf128_t* bf_q_hat_k      = malloc(sizeof(bf128_t) * 4);
      bf128_t* bf_q_hat_w_dash = malloc(sizeof(bf128_t) * 4);

      for (uint32_t r = 0; r <= 3; r++) {
        // TODO: ByteCombine
        bf_q_hat_k[(r + 3) % 4];
        bf_q_hat_w_dash[r];
      }
      // STep: 27
      for (uint32_t r = 0; r <= 3; r++) {
        // TODO: mul with 1_{F_{2^8}} ??
        bf128_mul(bf_q_hat_k[r], bf_q_hat_w_dash[r]);
      }
      if (lambda == 192) {
        iwd = iwd + 192;
      } else {
        iwd = iwd + 128;
      }
    }
  }

  return 1;
}

int aes_enc_forward(uint32_t lambda, uint32_t R, uint32_t m, uint32_t Lenc, const uint8_t* x,
                    uint8_t* xk, uint8_t* in, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                    uint8_t* y_out) {

  uint32_t lambdaBytes = lambda / 8;
  uint32_t x_len       = lambdaBytes * (Lenc / 8);
  uint32_t delta_len   = lambdaBytes;
  uint32_t in_len      = 16;
  uint32_t y_out_len   = lambdaBytes * R * 16;

  uint32_t bf_y_len = sizeof(bf128_t) * R * 16;
  bf128_t* bf_y     = malloc(bf_y_len);

  uint32_t ird;

  if (m == 1) {
    // If m == 1, Mkey should be 0, because we then cancel delta in line 10
    if (Mkey != 0) {
      return -1;
    }

    bf8_t bf_minus_mtag = bf8_from_bit(1 - Mtag);
    bf8_t bf_minus_mkey = bf8_from_bit(1 - Mkey);

    for (uint32_t i = 0; i < 16; i++) {
      bf8_t* bf_xin = malloc(8);
      for (uint32_t j = 0; j < 8; j++) {
        bf_xin[j] =
            bf8_mul(bf8_mul(bf8_from_bit(getBit(in + (8 * i), j)), bf_minus_mtag), bf_minus_mkey);
      }
      // TODO: ByteCombine
      bf_y[i];
    }
    uint32_t ix, ik, iy;
    for (uint32_t j = 0; j < R; j++) {
      for (uint32_t c = 0; c < 3; c++) {
        ix                 = 128 * (j - 1) + 32 * c;
        ik                 = 128 * j + 32 * c;
        iy                 = 16 * j + 4 * c;
        bf128_t* bf_x_hat  = malloc(sizeof(bf128_t) * 3);
        bf128_t* bf_xk_hat = malloc(sizeof(bf128_t) * 3);
        for (uint32_t r = 0; r < 3; r++) {
          // TODO: Bytecombine
          bf_x_hat[r];
          bf_xk_hat[r];
        }
        bf8_t one   = bf8_one();
        bf8_t two   = bf8_add(one, one);
        bf8_t three = bf8_add(two, one);
        // TODO: multiplication by {...}_{F_{2^8}}
        bf_y[iy + 0] = bf_x_hat[0];
        bf_y[iy + 1] = bf_x_hat[0];
        bf_y[iy + 2] = bf_x_hat[0];
        bf_y[iy + 3] = bf_x_hat[0];
      }
    }
    y_out        = malloc(bf_y_len);
    uint32_t idx = 0;
    for (uint32_t i = 0; i < bf_y_len; i += lambdaBytes) {
      bf128_store(y_out + i, bf_y[idx]);
      idx += 1;
    }
  } else {
    bf128_t bf_minus_mtag = bf128_from_bit(1 - Mtag);
    bf128_t bf_minus_mkey = bf128_from_bit(1 - Mkey);
    bf128_t bf_mkey       = bf128_from_bit(Mkey);

    for (uint32_t i = 0; i < 16; i++) {
      bf128_t* bf_xin = malloc(sizeof(bf128_t) * 8);
      for (uint32_t j = 0; j < 8; j++) {
        bf_xin[j] = bf128_mul(bf128_mul(bf128_from_bit(getBit(in + (8 * i), j)), bf_minus_mtag),
                              bf128_add(bf128_mul(bf_mkey, bf128_load(delta)), bf_minus_mkey));
      }
      // TODO: ByteCombine
      bf_y[i];
    }
    uint32_t ix, ik, iy;
    for (uint32_t j = 0; j < R; j++) {
      for (uint32_t c = 0; c < 3; c++) {
        ix                 = 128 * (j - 1) + 32 * c;
        ik                 = 128 * j + 32 * c;
        iy                 = 16 * j + 4 * c;
        bf128_t* bf_x_hat  = malloc(sizeof(bf128_t) * 3);
        bf128_t* bf_xk_hat = malloc(sizeof(bf128_t) * 3);
        for (uint32_t r = 0; r < 3; r++) {
          // TODO: Bytecombine
          bf_x_hat[r];
          bf_xk_hat[r];
        }
        bf8_t one   = bf8_one();
        bf8_t two   = bf8_add(one, one);
        bf8_t three = bf8_add(two, one);
        // TODO: multiplication by {...}_{F_{2^8}}
        bf_y[iy + 0] = bf_x_hat[0];
        bf_y[iy + 1] = bf_x_hat[0];
        bf_y[iy + 2] = bf_x_hat[0];
        bf_y[iy + 3] = bf_x_hat[0];
      }
    }
    y_out        = malloc(bf_y_len);
    uint32_t idx = 0;
    for (uint32_t i = 0; i < bf_y_len; i += lambdaBytes) {
      bf128_store(y_out + i, bf_y[idx]);
      idx += 1;
    }
  }
  return 1;
}

int aes_enc_backward(uint32_t lambda, uint32_t R, uint32_t m, uint32_t Lenc, const uint8_t* x,
                     uint8_t* xk, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, uint8_t* out,
                     uint8_t* y_out) {

  uint32_t lambdaBytes = lambda / 8;
  uint32_t x_len       = lambdaBytes * (Lenc / 8);
  uint32_t delta_len   = lambdaBytes;
  uint32_t out_len     = 16;
  uint32_t y_out_len   = lambdaBytes * R * 16;

  uint32_t bf_y_len = sizeof(bf128_t) * R * 16;
  bf128_t* bf_y     = malloc(bf_y_len);

  uint32_t ird;

  if (m == 1) {

    // If m == 1, Mkey should be 0, because we then cancel delta in line 10
    if (Mkey != 0) {
      return -1;
    }

    bf8_t* bf_x_tilde   = malloc(8);
    bf8_t* bf_xout      = malloc(8);
    bf8_t bf_minus_mtag = bf8_from_bit(1 - Mtag);
    bf8_t bf_minus_mkey = bf8_from_bit(1 - Mkey);

    for (uint32_t j = 0; j < R; j++) {
      for (uint32_t c = 0; c < 3; c++) {
        for (uint32_t r = 0; r < 3; r++) {
          ird = (128 * j) + (32 * (c - r % 4)) + (8 * r);

          if (j < (R - 1)) {

            for (uint32_t i = 0; i < 8; i++) {
              // converting from bit index to byte idx
              bf_x_tilde[i] = bf8_from_bit(getBit(x + (ird / 8), i));
            }

          } else {
            for (uint32_t i = 0; i < 8; i++) {
              // converting bit idx to byte idx
              bf_xout[i] =
                  bf8_mul(bf8_mul(bf8_from_bit(getBit(out[(ird - 1152) / 8], i)), bf_minus_mtag),
                          bf_minus_mkey);
              bf_x_tilde[i] = bf8_add(bf_xout[i], bf8_from_bit(getBit(xk + (128 + ird) / 8, i)));
            }
          }
          bf8_t* bf_y_tilde = malloc(8);
          // TODO: multiplication with different field size ???
          bf_y_tilde[7] = bf8_add(bf8_add(bf8_add(bf_x_tilde[5], bf_x_tilde[2]), bf_x_tilde[0]),
                                  bf8_mul(bf_minus_mtag, bf_minus_mkey));
          bf_y_tilde[6] = bf8_add(bf8_add(bf_x_tilde[7], bf_x_tilde[4]), bf_x_tilde[1]);
          bf_y_tilde[5] = bf8_add(bf8_add(bf8_add(bf_x_tilde[6], bf_x_tilde[3]), bf_x_tilde[0]),
                                  bf8_mul(bf_minus_mtag, bf_minus_mkey));
          bf_y_tilde[4] = bf8_add(bf8_add(bf_x_tilde[7], bf_x_tilde[5]), bf_x_tilde[2]);
          bf_y_tilde[3] = bf8_add(bf8_add(bf_x_tilde[6], bf_x_tilde[4]), bf_x_tilde[1]);
          bf_y_tilde[2] = bf8_add(bf8_add(bf_x_tilde[5], bf_x_tilde[3]), bf_x_tilde[0]);
          bf_y_tilde[1] = bf8_add(bf8_add(bf_x_tilde[7], bf_x_tilde[4]), bf_x_tilde[2]);
          bf_y_tilde[0] = bf8_add(bf8_add(bf_x_tilde[6], bf_x_tilde[3]), bf_x_tilde[1]);

          // TODO: ByteCombine is missing !!!!
          bf_y[16 * j + 4 * c + r];
        }
      }
    }
    y_out        = malloc(y_out_len);
    uint32_t idx = 0;
    for (uint32_t i = 0; i < y_out_len; i += lambdaBytes) {
      bf128_store(y_out + i, bf_y[idx]);
      idx++;
    }
  } else {

    bf128_t* bf_x_tilde   = malloc(sizeof(bf128_t) * 8);
    bf128_t* bf_xout      = malloc(sizeof(bf128_t) * 8);
    bf128_t bf_minus_mtag = bf128_from_bit(1 - Mtag);
    bf128_t bf_minus_mkey = bf128_from_bit(1 - Mkey);
    bf128_t bf_mkey       = bf128_from_bit(Mkey);

    for (uint32_t j = 0; j < R; j++) {
      for (uint32_t c = 0; c < 3; c++) {
        for (uint32_t r = 0; r < 3; r++) {
          ird = (128 * j) + (32 * (c - r % 4)) + (8 * r);

          if (j < (R - 1)) {

            for (uint32_t i = 0; i < 8; i++) {
              // TODO: check here if getting the bit makes sense and the converting to bf128 !!
              bf_x_tilde[i] = bf128_from_bit(getBit(x + (ird / 8), i));
            }

          } else {
            for (uint32_t i = 0; i < 8; i++) {
              // TODO: hopefully the endianness is fine
              uint8_t* out_zeros_concat = malloc(lambdaBytes);
              out_zeros_concat[0]       = (getBit(out[(ird - 1152) / 8], i) << 7);
              memset(out_zeros_concat + 1, 0, lambdaBytes - 1);

              bf_xout[i] =
                  bf128_mul(bf128_mul(bf128_load(out_zeros_concat), bf_minus_mtag),
                            bf128_add(bf128_mul(bf_mkey, bf128_load(delta)), bf_minus_mkey));

              // TODO: check here if getting the bit makes sense and the converting to bf128 !!
              bf_x_tilde[i] =
                  bf128_add(bf_xout[i], bf128_from_bit(getBit(xk + ((128 + ird) / 8), i)));
            }
          }
          bf128_t* bf_y_tilde = malloc(8);
          // TODO: multiplication with different field size ???
          bf_y_tilde[7] =
              bf128_add(bf128_add(bf128_add(bf_x_tilde[5], bf_x_tilde[2]), bf_x_tilde[0]),
                        bf128_mul(bf_minus_mtag, bf_minus_mkey));
          bf_y_tilde[6] = bf128_add(bf128_add(bf_x_tilde[7], bf_x_tilde[4]), bf_x_tilde[1]);
          bf_y_tilde[5] =
              bf128_add(bf128_add(bf128_add(bf_x_tilde[6], bf_x_tilde[3]), bf_x_tilde[0]),
                        bf128_mul(bf_minus_mtag, bf_minus_mkey));
          bf_y_tilde[4] = bf128_add(bf128_add(bf_x_tilde[7], bf_x_tilde[5]), bf_x_tilde[2]);
          bf_y_tilde[3] = bf128_add(bf128_add(bf_x_tilde[6], bf_x_tilde[4]), bf_x_tilde[1]);
          bf_y_tilde[2] = bf128_add(bf128_add(bf_x_tilde[5], bf_x_tilde[3]), bf_x_tilde[0]);
          bf_y_tilde[1] = bf128_add(bf128_add(bf_x_tilde[7], bf_x_tilde[4]), bf_x_tilde[2]);
          bf_y_tilde[0] = bf128_add(bf128_add(bf_x_tilde[6], bf_x_tilde[3]), bf_x_tilde[1]);

          // TODO: ByteCombine is missing !!!!
          bf_y[16 * j + 4 * c + r];
        }
      }
    }
    y_out        = malloc(y_out_len);
    uint32_t idx = 0;
    for (uint32_t i = 0; i < y_out_len; i += lambdaBytes) {
      bf128_store(y_out + i, bf_y[idx]);
      idx++;
    }
  }
  return 1;
}

int aes_enc_constraints(uint32_t lambda, uint32_t R, uint32_t Lenc, uint32_t Senc,
                        const uint8_t* in, const uint8_t* out, const uint8_t* w, const uint8_t* v,
                        const uint8_t* k, const uint8_t* vk, uint8_t Mkey, const uint8_t* q,
                        const uint8_t* qk, const uint8_t* delta, uint8_t* A0, uint8_t* A1,
                        uint8_t* B) {

  // TODO: unused vars, remove them in the end !!
  uint32_t lambdaBytes = lambda / 8;
  uint32_t in_len      = 16;
  uint32_t out_len     = 16;
  uint32_t w_len       = Lenc;
  uint32_t v_len       = Lenc;
  uint32_t k_len       = (R + 1) * 128;
  uint32_t vk_len      = ((R + 1) * 128) * lambdaBytes;
  uint32_t q_len       = Lenc * lambdaBytes;
  uint32_t qk_len      = ((R + 1) * 128) * lambdaBytes;
  uint32_t delta_len   = lambdaBytes;

  if (Mkey == 0) {
    uint8_t *s, vs, s_dash, vs_dash;
    aes_enc_forward(lambda, R, 1, Lenc, w, k, in, 0, 0, NULL, s);
    aes_enc_forward(lambda, R, lambda, Lenc, v, vk, in, 1, 0, NULL, vs);
    aes_enc_backward(lambda, R, 1, Lenc, w, k, 0, 0, delta, out, s_dash);
    aes_enc_backward(lambda, R, lambda, Lenc, v, vk, 1, 0, delta, out, vs_dash);

    bf128_t* bf_s       = malloc(sizeof(bf128_t) * Senc);
    bf128_t* bf_vs      = malloc(sizeof(bf128_t) * Senc);
    bf128_t* bf_s_dash  = malloc(sizeof(bf128_t) * Senc);
    bf128_t* bf_vs_dash = malloc(sizeof(bf128_t) * Senc);
    uint32_t idx        = 0;
    for (uint8_t i = 0; i < Senc; i++) {
      bf_s[i]       = bf128_load(s + idx);
      bf_vs[i]      = bf128_load(vs + idx);
      bf_s_dash[i]  = bf128_load(s_dash + idx);
      bf_vs_dash[i] = bf128_load(vs_dash + idx);
      idx += lambdaBytes;
    }
    bf128_t* bf_a0 = malloc(sizeof(bf128_t) * Senc);
    bf128_t* bf_a1 = malloc(sizeof(bf128_t) * Senc);
    A0             = malloc(lambdaBytes * Senc);
    A1             = malloc(lambdaBytes * Senc);
    for (uint32_t j = 0; j < Senc; j++) {
      bf_a0[j] = bf128_mul(bf_vs[j], bf_vs_dash[j]);
      bf128_store(A0 + (lambdaBytes * j), bf_a0[j]);
      // TODO : Subtraction by 1_{F_{2^8}}
      bf_a1[j] =
          bf128_add(bf128_mul(bf128_add(bf_s[j], bf_vs[j]), bf128_add(bf_s_dash[j], bf_vs_dash[j])),
                    bf_a0[j]);
      bf128_store(A1 + (lambdaBytes * j), bf_a1[j]);
    }
  } else {
    uint8_t *qs, qs_dash;
    aes_enc_forward(lambda, R, lambda, Lenc, q, qk, in, 0, 1, delta, qs);
    aes_enc_backward(lambda, R, lambda, Lenc, q, qk, 0, 1, delta, out, qs_dash);

    bf128_t* bf_qs      = malloc(sizeof(bf128_t) * Senc);
    bf128_t* bf_qs_dash = malloc(sizeof(bf128_t) * Senc);
    uint32_t idx        = 0;
    for (uint32_t i = 0; i < Senc; i++) {
      bf_qs[i]      = bf128_load(qs + idx);
      bf_qs_dash[i] = bf128_load(qs_dash + idx);
      idx += lambdaBytes;
    }
    bf128_t* bf_b = malloc(sizeof(bf128_t) * Senc);
    B             = malloc(lambdaBytes * Senc);
    for (uint32_t j = 0; j < Senc; j++) {
      // TODO: multiplication by 1_{F_{2^8}}
      bf_b[j] = bf128_mul(bf_qs[j], bf_qs_dash[j]);
      bf128_store(B + (lambdaBytes * j), bf_b[j]);
    }
  }
}

void aes_prove(uint8_t* w, uint8_t* u, uint8_t** V, uint8_t* in, uint8_t* out, uint8_t* chal,
               uint32_t lambda, uint32_t R, uint32_t tau, uint32_t l, uint32_t beta, uint32_t Lke,
               uint32_t Lenc, uint32_t C, uint32_t Nwd, uint32_t Ske, uint32_t Senc,
               uint8_t* a_tilde, uint8_t* b_tilde) {

  uint32_t lambdaBytes = lambda / 8;
  uint32_t w_len       = l;
  uint32_t u_len       = lambda;
  uint32_t in_len      = beta * 128;
  uint32_t out_len     = beta * 128;
  uint32_t chal        = (3 * lambda) + 64;

  uint8_t* in_0  = malloc(16);
  uint8_t* out_0 = malloc(16);
  memcpy(in_0, in, 16);
  memcpy(out_0, out, 16);

  uint8_t* in_1  = malloc(16);
  uint8_t* out_1 = malloc(16);
  if (beta == 2) {
    memcpy(in_1, in + 16, 16);
    memcpy(out_1, out + 16, 16);
  }

  uint8_t* w_tilde = malloc(Lke / 8);
  uint8_t* v_tilde = malloc(Lke / 8);
  memcpy(w_tilde, w, Lke / 8);
  // TODO: fix V copy
  memcpy(v_tilde, V, Lke / 8);

  uint8_t *A0, *A1, *k, *vk, *B, *qk;
  aes_key_schedule_constraints(lambda, R, Nwd, Ske, Lke, w_tilde, v_tilde, 0, NULL, NULL, A0, A1, k,
                               vk, B, qk);

  uint32_t lByte_enc;
  memcpy(w_tilde, w + Lke, Lenc);
  // TODO: fix V copy
  memcpy(v_tilde, V + Lke, Lenc);

  uint8_t *A0_1, *A1_1;

  aes_enc_constraints(lambda, R, Lenc, Senc, in_0, out_0, w_tilde, v_tilde, k, vk, 0, NULL, NULL,
                      NULL, A0_1, A1_1, B);

  A0 = realloc(A0, lambdaBytes * (Ske / 8) + (lambdaBytes * Senc));
  A1 = realloc(A1, lambdaBytes * (Ske / 8) + (lambdaBytes * Senc));
  memcpy(A0 + (lambdaBytes * (Ske / 8)), A0_1, lambdaBytes * Senc);
  memcpy(A1 + (lambdaBytes * (Ske / 8)), A1_1, lambdaBytes * Senc);

  if (beta == 2) {
    uint8_t *A0_2, *A1_2;
    memcpy(w_tilde, w + Lke + Lenc, l - (Lke + Lenc));
    memcpy(v_tilde, V + Lke + Lenc, l - (Lke + Lenc));
    aes_enc_constraints(lambda, R, Lenc, Senc, in_1, out_1, w_tilde, v_tilde, k, vk, 0, NULL, NULL,
                        NULL, A0_2, A1_2, B);

    A0 = realloc(A0, (lambdaBytes * (Ske / 8)) + (2 * (lambdaBytes * Senc)));
    A1 = realloc(A1, (lambdaBytes * (Ske / 8)) + (2 * (lambdaBytes * Senc)));
    memcpy(A0 + (lambdaBytes * (Ske / 8)) + (lambdaBytes * Senc), A0_2, lambdaBytes * Senc);
    memcpy(A1 + (lambdaBytes * (Ske / 8)) + (lambdaBytes * Senc), A1_2, lambdaBytes * Senc);
  }

  bf128_t* bf_u = malloc(sizeof(bf128_t) * lambda);
  uint32_t idx  = 0;
  for (uint32_t i = l; i < l + lambda; i++) {
    bf_u[i] = bf128_load(u + idx);
    idx += lambdaBytes;
  }
  // TODO: unsure what is happeninig in line 18
  // TODO: check if the sizes have been assigned correctly
  uint8_t* us = malloc(lambdaBytes * lambdaBytes);
  uint8_t* vs = malloc(lambdaBytes * lambdaBytes);

  uint8_t* a1_us_concat = malloc(((lambdaBytes * (Ske / 8)) + (2 * (lambdaBytes * Senc))) +
                                 sizeof(lambdaBytes * lambdaBytes));
  uint8_t* a0_vs_concat = malloc(((lambdaBytes * (Ske / 8)) + (2 * (lambdaBytes * Senc))) +
                                 sizeof(lambdaBytes * lambdaBytes));

  // TODO: what was ell, r and t ?
  zk_hash_128(a_tilde, r, chal, t, a1_us_concat, l);
  zk_hash_128(b_tilde, r, chal, t, a0_vs_concat, l);
}

bool aes_verify(uint8_t* d, uint8_t** Q, uint8_t* chal_2, uint8_t* chal_3, uint8_t* a_tilde,
                uint8_t* b_tilde, uint8_t* in, uint8_t* out, uint32_t lambda, uint32_t tau,
                uint32_t l, uint32_t beta, uint32_t R, uint32_t Nwd, uint32_t Ske, uint32_t Lke,
                uint32_t Lenc, uint32_t Senc, uint32_t C, uint32_t k0, uint32_t k1, uint32_t t0,
                uint32_t t1) {

  uint32_t lambdaBytes = lambda / 8;
  uint32_t d_len       = l;
  uint32_t Q_len       = tau;
  uint32_t Qi_len      = (l + lambda) + lambda;
  uint32_t chal_2_len  = 3 * lambda + 64;
  uint32_t chal_3_len  = lambda;
  uint32_t in_len      = beta * 128;
  uint32_t out_len     = beta * 128;

  bf128_t bf_delta = bf128_load(chal_3);
  uint8_t* delta   = malloc(16);
  memcpy(delta, chal_3, 16);

  uint8_t* in_0  = malloc(16);
  uint8_t* out_0 = malloc(16);
  memcpy(in_0, in, 16);
  memcpy(out_0, out, 16);

  uint8_t* in_1  = malloc(16);
  uint8_t* out_1 = malloc(16);
  if (beta == 2) {
    memcpy(in_1, in + 16, 16);
    memcpy(out_1, out + 16, 16);
  }

  uint32_t kb;
  for (uint32_t i = 0; i < tau; i++) {
    if (i < t0) {
      kb = k0;
    } else {
      kb = k1;
    }
    uint8_t* d = malloc(kb);
    ChalDec(chal_3, i, k0, t0, k1, t1, d);

    for (uint32_t j = 0; j < kb; j++) {
      // TODO: have a look what is happeninig here !!
    }
  }

  bf128_t* bf_q  = malloc(sizeof(bf128_t) * (l + tau));
  uint8_t* q_lke = malloc(lambdaBytes * (Lke));
  for (uint32_t i = 0; i < l + tau; i++) {
    bf_q[i] = bf128_load(Q[i]);
  }
  memcpy(q_lke, bf_q, sizeof(bf128_t) * Lke);

  bf128_t* bf_b = malloc(sizeof(bf128_t) * C);

  uint8_t *A0, *A1, *k, *vk, *B_0, *qk;
  aes_key_schedule_constraints(lambda, R, Nwd, Ske, Lke, NULL, NULL, 1, q_lke, delta, A0, A1, k, vk,
                               B_0, qk);

  uint8_t *A0_1, *A1_1, *B_1;
  uint8_t* q_lke_lenc = malloc(lambdaBytes * (Lenc - Lke));
  memcpy(q_lke_lenc, bf_q + Lke, sizeof(bf128_t) * (Lenc - Lke));
  aes_enc_constraints(lambda, R, Lenc, Senc, in_0, out_0, NULL, NULL, NULL, NULL, 1, q_lke_lenc, qk,
                      delta, A0_1, A1_1, B_1);

  uint8_t* B;
  if (beta == 1) {
    B = malloc(16 * R * 2);
    memcpy(B, B_0, 16 * R);
    memcpy(B + (R * 16), B_1, 16 * R);
  } else {
    uint8_t *A0_2, *A1_2, *B_2;
    uint8_t* q_lke_enc_l = malloc(lambdaBytes * (l - (Lke + Lenc)));
    memcpy(q_lke_enc_l, bf_b + (Lke + Lenc), sizeof(bf128_t) * (l - (Lke + Lenc)));
    aes_enc_constraints(lambda, R, Lenc, Senc, in_1, out_1, NULL, NULL, NULL, NULL, 1, q_lke_enc_l,
                        qk, delta, A0_2, A1_2, B_2);
    memcpy(B, B_0, (R * 16));
    memcpy(B + (R * 16), B_1, (R * 16));
    memcpy(B + (R * 16 * 2), B_2, (R * 16));
  }

  // TODO: line 18 and do till end !!
  // TODO: what is ell ??
  uint8_t* qs;
  uint8_t* q_tilde;
  uint8_t* b_qs_conat;
  uint8_t *r, t;
  uint32_t ell;
  switch (lambda) {
  case 256:
    zk_hash_256(q_tilde, chal_2, b_qs_conat, ell);
    break;
  case 192:
    zk_hash_192(q_tilde, chal_2, b_qs_conat, ell);
    break;
  default:
    zk_hash_128(q_tilde, chal_2, b_qs_conat, ell);
    break;
  }

  // TODO: Do field operation here !!
  uint8_t* q_check = malloc(sizeof(a_tilde));
  for (uint32_t i = 0; i < sizeof(q_check); i++) {
    q_check[i] = a_tilde[i] + (b_tilde[i] * delta[i]);
  }

  if (memcmp(q_tilde, q_check, sizeof(q_check)) == 0) {
    return true;
  } else {
    return false;
  }
}