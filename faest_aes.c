#include "faest_aes.h"
// TODO: Do not pass lambdaBytes everywhere, compute it in the function....
// TODO: change q to Q where applicable

// Bwd -> block_words,,, for EM mode
void aes_extend_witness(uint32_t lambda, uint32_t R, uint32_t Nwd, uint32_t Bwd, uint32_t l,
                        uint32_t Ske, uint32_t beta, const uint8_t* key, const uint8_t** in,
                        uint8_t* w_out) {

  // Step 1..5
  uint32_t w_byte_idx = 0;
  w_out               = malloc(l + 7 / 8);
  aes_round_key_t* round_keys;
  expand_key(round_keys, key, Nwd, 4, R);
  memcpy(w_out, round_keys, Nwd);
  w_byte_idx += Nwd;
  uint32_t ik = Nwd;
  // Step: 6..9
  for (uint32_t j = 0; j < Ske / 4; j++) {
    memcpy(w_out + w_byte_idx, round_keys[ik], 4);
    w_byte_idx += 4;
    if (lambda == 192) {
      ik += 6;
    } else {
      ik += 4;
    }
  }
  // Step: 10..19
  for (uint32_t b = 0; b < beta; b++) {
    uint8_t* state = malloc(16);
    memcpy(state, in[b], 16);
    add_round_key(0, state, round_keys, Bwd);
    for (uint32_t j = 1; j < R; j++) {
      sub_bytes(state, Bwd);
      shift_row(state, Bwd);
      memcpy(w_out + w_byte_idx, state, sizeof(state));
      w_byte_idx += sizeof(state);
      mix_column(state, Bwd);
      add_round_key(j, state, round_keys, Bwd);
    }
  }
}

// TODO: generalize bf128_t to bf(lambda)_t
int aes_key_schedule_forward(uint32_t lambda, uint32_t R, uint32_t Nwd, uint32_t m,
                             const uint8_t* x, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                             uint8_t* y_out) {

  uint32_t lambdaByte = lambda / 8;

  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return 0;
  }

  // STep: 2..3
  if (m == 1) {
    bf8_t* bf_y = malloc(sizeof(bf8_t) * ((R + 1) * 128));
    for (uint32_t i = 0; i < lambda; i++) {
      bf_y[i] = bf8_load(x + i);
    }
    // Step: 4
    uint32_t i_wd = lambdaByte;

    // Step: 5..10
    for (uint32_t j = Nwd; j < 4 * (R + 1); j++) {
      if (j % Nwd == 0 || (Nwd > 6 && j % Nwd == 4)) {
        for (uint32_t i = (j * 32); i < ((j * 32) + 31); i++) {
          bf_y[i] = bf8_load(x + i);
        }
        i_wd += 32;
      } else {
        for (uint32_t i = 0; i < 32; i++) {
          bf_y[32 * j + i] = bf8_add(bf_y[32 * (j - Nwd) + i], bf_y[32 * (j - 1) + i]);
        }
      }
    }
    y_out = malloc(sizeof(bf_y));
    for (uint32_t i = 0; i < sizeof(bf_y); i++) {
      bf8_store(y_out + i, bf_y[i]);
    }
    return 1;

  } else {
    bf128_t* bf_y = malloc(sizeof(bf128_t) * ((R + 1) * 128));
    for (uint32_t i = 0; i < lambda; i++) {
      bf_y[i] = bf128_load(x + (i * lambdaByte));
    }
    // Step: 4
    uint32_t i_wd = lambdaByte;

    // Step: 5..10
    for (uint32_t j = Nwd; j < 4 * (R + 1); j++) {
      if (j % Nwd == 0 || (Nwd > 6 && j % Nwd == 4)) {
        for (uint32_t i = (j * 32); i < ((j * 32) + 31); i++) {
          bf_y[i] = bf128_load(x + (lambda * i));
        }
        i_wd += 32;
      } else {
        for (uint32_t i = 0; i < 32; i++) {
          bf_y[32 * j + i] = bf128_add(bf_y[32 * (j - Nwd) + i], bf_y[32 * (j - 1) + i]);
        }
      }
    }
    y_out = malloc(sizeof(bf_y));
    for (uint32_t i = 0; i < sizeof(bf_y) / sizeof(bf_y[0]); i++) {
      bf128_store(y_out + (i * sizeof(bf_y[0])), bf_y[i]);
    }
    return 1;
  }
}

// TODO: generalize bf128_t to bf(lambda)_t
int aes_key_schedule_backward(uint32_t lambda, uint32_t R, uint32_t Nwd, uint32_t Ske, uint8_t Lke,
                              uint32_t m, const uint8_t* x, const uint8_t* xk, uint8_t Mtag,
                              uint8_t Mkey, const uint8_t* delta, uint8_t* y_out) {

  uint32_t lambdaByte = lambda / 8;

  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return 0;
  }

  bf128_t bf_delta = bf128_load(delta);

  // STep: 2
  if (m == 1) {
    bf8_t* bf_y    = malloc(8 * Ske);
    uint32_t iwd   = 0;
    uint32_t c     = 0;
    bool rmvRcon   = true;
    uint32_t ircon = 0;

    bf8_t* bf_x = malloc(Lke);
    for (uint32_t i = 0; i < Lke; i++) {
      bf_x[i] = bf8_load(x + i);
    }
    bf8_t* bf_xk = malloc(Lke);
    for (uint32_t i = 0; i < Lke; i++) {
      bf_xk[i] = bf8_load(xk + i);
    }
    bf8_t* bf_x_tilde = malloc(8);

    for (uint32_t j = 0; j < Ske; j++) {
      bf_x_tilde[0] = bf8_add(bf_x[8 * j], bf_xk[iwd + 8 * c]);
      bf_x_tilde[1] = bf8_add(bf_x[8 * j + 1], bf_xk[iwd + (8 * c) + 1]);
      bf_x_tilde[2] = bf8_add(bf_x[8 * j + 2], bf_xk[iwd + (8 * c) + 2]);
      bf_x_tilde[3] = bf8_add(bf_x[8 * j + 3], bf_xk[iwd + (8 * c) + 3]);
      bf_x_tilde[4] = bf8_add(bf_x[8 * j + 4], bf_xk[iwd + (8 * c) + 4]);
      bf_x_tilde[5] = bf8_add(bf_x[8 * j + 5], bf_xk[iwd + (8 * c) + 5]);
      bf_x_tilde[6] = bf8_add(bf_x[8 * j + 6], bf_xk[iwd + (8 * c) + 6]);
      bf_x_tilde[7] = bf8_add(bf_x[8 * j + 7], bf_xk[iwd + (8 * c) + 7]);

      if (Mtag == 0 && rmvRcon == true && c == 0) {
        // TODO: step:9, not sure what is happening here
        // rcon = ToBits(...) ??
        bf8_t* rcon = malloc(8);
        ircon       = ircon + 1;
        bf8_t* bf_r = malloc(8);
        for (uint32_t i = 0; i < 8; i++) {
          // TODO:: complet this after clarification !!
          bf8_t bf_Mkey;
          bf8_t bf_t;
          if (Mkey == 1) {
            bf_Mkey = bf8_one();
            bf_t    = bf8_zero();
          } else {
            bf_Mkey = bf8_zero();
            bf_t    = bf8_one();
          }
          // TODO: Delta is in F_(2^lambda), but bf_Mkey is F_2 ??
          bf_r[i] = bf8_mul(rcon[i], bf8_add(bf8_mul(bf_Mkey, bf_delta), bf_t));
        }
        for (uint32_t i = 0; i < 8; i++) {
          bf_x_tilde[i] = bf8_add(bf_x_tilde[i], bf_r[i]);
        }
      }
      bf8_t* bf_y_tilde = malloc(8);
      // TODO: Delta is in F_(2^lambda), but bf_Mkey is F_2 ??
      bf_y_tilde[0] = bf8_add(bf8_add(bf_x_tilde[2], bf_x_tilde[5]), bf_x_tilde[7]);
      bf_y_tilde[1] = bf8_add(bf8_add(bf_x_tilde[0], bf_x_tilde[3]), bf_x_tilde[6]);
      // TODO: Delta is in F_(2^lambda), but bf_Mkey is F_2 ??
      bf_y_tilde[2] = bf8_add(bf8_add(bf_x_tilde[1], bf_x_tilde[4]), bf_x_tilde[7]);
      bf_y_tilde[3] = bf8_add(bf8_add(bf_x_tilde[0], bf_x_tilde[2]), bf_x_tilde[5]);
      bf_y_tilde[4] = bf8_add(bf8_add(bf_x_tilde[1], bf_x_tilde[3]), bf_x_tilde[6]);
      bf_y_tilde[5] = bf8_add(bf8_add(bf_x_tilde[2], bf_x_tilde[4]), bf_x_tilde[7]);
      bf_y_tilde[6] = bf8_add(bf8_add(bf_x_tilde[0], bf_x_tilde[3]), bf_x_tilde[5]);
      bf_y_tilde[7] = bf8_add(bf8_add(bf_x_tilde[1], bf_x_tilde[4]), bf_x_tilde[6]);

      bf_y[8 * j]     = bf_y_tilde[0];
      bf_y[8 * j + 1] = bf_y_tilde[1];
      bf_y[8 * j + 2] = bf_y_tilde[2];
      bf_y[8 * j + 3] = bf_y_tilde[3];
      bf_y[8 * j + 4] = bf_y_tilde[4];
      bf_y[8 * j + 5] = bf_y_tilde[5];
      bf_y[8 * j + 6] = bf_y_tilde[6];
      bf_y[8 * j + 7] = bf_y_tilde[7];

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
    y_out = malloc(sizeof(bf_y));
    for (uint32_t i = 0; i < sizeof(bf_y); i++) {
      bf8_store(y_out, bf_y[i]);
    }
    return 1;
  } else {
    bf128_t* bf_y  = malloc(8 * Ske * sizeof(bf128_t));
    uint32_t iwd   = 0;
    uint32_t c     = 0;
    bool rmvRcon   = true;
    uint32_t ircon = 0;

    bf128_t* bf_x = malloc(Lke * sizeof(bf128_t));
    for (uint32_t i = 0; i < Lke; i++) {
      bf_x[i] = bf128_load(x + (i * lambdaByte));
    }
    bf128_t* bf_xk = malloc(Lke * sizeof(bf128_t));
    for (uint32_t i = 0; i < Lke; i++) {
      bf_xk[i] = bf128_load(xk + (i * lambdaByte));
    }
    bf128_t* bf_x_tilde = malloc(8 * sizeof(bf128_t));

    for (uint32_t j = 0; j < Ske; j++) {
      bf_x_tilde[0] = bf128_add(bf_x[8 * j], bf_xk[iwd + 8 * c]);
      bf_x_tilde[1] = bf128_add(bf_x[8 * j + 1], bf_xk[iwd + (8 * c) + 1]);
      bf_x_tilde[2] = bf128_add(bf_x[8 * j + 2], bf_xk[iwd + (8 * c) + 2]);
      bf_x_tilde[3] = bf128_add(bf_x[8 * j + 3], bf_xk[iwd + (8 * c) + 3]);
      bf_x_tilde[4] = bf128_add(bf_x[8 * j + 4], bf_xk[iwd + (8 * c) + 4]);
      bf_x_tilde[5] = bf128_add(bf_x[8 * j + 5], bf_xk[iwd + (8 * c) + 5]);
      bf_x_tilde[6] = bf128_add(bf_x[8 * j + 6], bf_xk[iwd + (8 * c) + 6]);
      bf_x_tilde[7] = bf128_add(bf_x[8 * j + 7], bf_xk[iwd + (8 * c) + 7]);

      if (Mtag == 0 && rmvRcon == true && c == 0) {
        // TODO: step:9, not sure what is happening here
        // rcon = ToBits(...) ??
        bf128_t* rcon = malloc(8 * sizeof(bf128_t));
        ircon         = ircon + 1;
        bf128_t* bf_r = malloc(8 * sizeof(bf128_t));
        for (uint32_t i = 0; i < 8; i++) {
          // TODO:: complet this after clarification !!
          bf128_t bf_Mkey;
          bf128_t bf_t;
          if (Mkey == 1) {
            bf_Mkey = bf128_one();
            bf_t    = bf128_zero();
          } else {
            bf_Mkey = bf128_zero();
            bf_t    = bf128_one();
          }
          // TODO: Delta is in F_(2^lambda), but bf_Mkey is F_2 ??
          bf_r[i] = bf128_mul(rcon[i], bf128_add(bf128_mul(bf_Mkey, bf_delta), bf_t));
        }
        for (uint32_t i = 0; i < 8; i++) {
          bf_x_tilde[i] = bf128_add(bf_x_tilde[i], bf_r[i]);
        }
      }
      bf128_t* bf_y_tilde = malloc(8 * sizeof(bf128_t));
      // TODO: Delta is in F_(2^lambda), but bf_Mkey is F_2 ??
      bf_y_tilde[0] = bf128_add(bf128_add(bf_x_tilde[2], bf_x_tilde[5]), bf_x_tilde[7]);
      bf_y_tilde[1] = bf128_add(bf128_add(bf_x_tilde[0], bf_x_tilde[3]), bf_x_tilde[6]);
      // TODO: Delta is in F_(2^lambda), but bf_Mkey is F_2 ??
      bf_y_tilde[2] = bf128_add(bf128_add(bf_x_tilde[1], bf_x_tilde[4]), bf_x_tilde[7]);
      bf_y_tilde[3] = bf128_add(bf128_add(bf_x_tilde[0], bf_x_tilde[2]), bf_x_tilde[5]);
      bf_y_tilde[4] = bf128_add(bf128_add(bf_x_tilde[1], bf_x_tilde[3]), bf_x_tilde[6]);
      bf_y_tilde[5] = bf128_add(bf128_add(bf_x_tilde[2], bf_x_tilde[4]), bf_x_tilde[7]);
      bf_y_tilde[6] = bf128_add(bf128_add(bf_x_tilde[0], bf_x_tilde[3]), bf_x_tilde[5]);
      bf_y_tilde[7] = bf128_add(bf128_add(bf_x_tilde[1], bf_x_tilde[4]), bf_x_tilde[6]);

      bf_y[8 * j]     = bf_y_tilde[0];
      bf_y[8 * j + 1] = bf_y_tilde[1];
      bf_y[8 * j + 2] = bf_y_tilde[2];
      bf_y[8 * j + 3] = bf_y_tilde[3];
      bf_y[8 * j + 4] = bf_y_tilde[4];
      bf_y[8 * j + 5] = bf_y_tilde[5];
      bf_y[8 * j + 6] = bf_y_tilde[6];
      bf_y[8 * j + 7] = bf_y_tilde[7];

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
    y_out = malloc(sizeof(bf_y));
    for (uint32_t i = 0; i < sizeof(bf_y) / sizeof(bf_y[0]); i++) {
      bf128_store(y_out + (i * sizeof(bf_y[0])), bf_y[i]);
    }
    return 1;
  }
}

// TODO
void aes_key_schedule_constraints(uint32_t lambda, uint32_t R, uint32_t Nwd, uint32_t Ske,
                                  uint32_t Lke, const uint8_t* w, const uint8_t* v,
                                  const uint8_t Mkey, const uint8_t* q, const uint8_t* delta,
                                  uint8_t** A, uint8_t* k, uint8_t* vk, uint8_t* B, uint8_t* qk) {
  uint32_t lambdaByte = lambda / 8;
  if (Mkey == 0) {
    // STep: 2
    uint8_t *k, vk, w_dash, v_w_dash;
    aes_key_schedule_forward(lambda, R, Nwd, 1, w, 0, 0, NULL, k);

    // Step: 3
    aes_key_schedule_forward(lambda, R, Nwd, lambda, v, 1, 0, NULL, vk);

    // Step: 4
    uint8_t* w_lambda = malloc(sizeof(w) - lambdaByte);
    memcpy(w_lambda, w + lambdaByte, sizeof(w) - lambdaByte);
    aes_key_schedule_backward(lambda, R, Nwd, Ske, Lke, 1, w_lambda, k, 0, 0, NULL, w_dash);

    // Step: 5
    uint8_t* v_lambda = malloc(sizeof(v) - lambdaByte);
    memcpy(v_lambda, v + lambdaByte, sizeof(v) - lambdaByte);
    aes_key_schedule_backward(lambda, R, Nwd, Ske, Lke, lambda, v_lambda, vk, 1, 0, NULL, v_w_dash);

    // Step: 6..8
    uint32_t iwd = 32 * (Nwd - 1);
    for (uint32_t j = 0; j < Ske / 4; j++) {

      for (uint32_t r = 0; r <= 3; r) {
        // TODO: line 9..12
      }
      // Step: 13..15
      for (uint32_t r = 0; r <= 3; r) {
        // TODO: line 14..15
      }
      if (lambda == 192) {
        iwd = iwd + 192;
      } else {
        iwd = iwd + 128;
      }
    }
  } else {
    // Step: 18..19
    aes_key_schedule_forward(lambda, R, Nwd, lambda, q, 0, 1, delta, qk);
    uint8_t* q_w_dash;
    uint8_t* q_lambda = malloc(sizeof(q) - lambdaByte);
    memcpy(q_lambda, q + lambdaByte, lambdaByte);
    aes_key_schedule_backward(lambda, R, Nwd, Ske, Lke, lambda, q_lambda, qk, 0, 1, delta,
                              q_w_dash);
    // Step 20..22
    uint32_t iwd = 32 * (Nwd - 1);
    for (uint32_t j = 0; j < Ske / 4; j++) {
      for (uint32_t r = 0; r <= 3; r++) {
        // TODO: steps 25..26
      }
      // STep: 27
      for (uint32_t r = 0; r <= 3; r++) {
        // TODO: step 28
      }
      if (lambda == 192) {
        iwd = iwd + 192;
      } else {
        iwd = iwd + 128;
      }
    }
  }
}

int aes_enc_forward(uint32_t lambda, uint32_t R, uint32_t m, uint32_t Lenc, const uint8_t* x,
                    uint8_t* xk, uint8_t* in, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                    uint8_t* y_out) {

  uint32_t lambdaBytes = lambda / 8;
  uint32_t x_len       = lambdaBytes * Lenc;
  uint32_t delta_len   = lambdaBytes;
  uint32_t in_len      = 16;
  uint32_t y_out_len   = lambdaBytes * R * 16;

  uint32_t bf_y_len = sizeof(bf128_t) * R * 16;
  bf128_t* bf_y     = malloc(bf_y_len);

  uint32_t ird;

  if (m == 1) {
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
}

int aes_enc_backward(uint32_t lambda, uint32_t R, uint32_t m, uint32_t Lenc, const uint8_t* x,
                     uint8_t* xk, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, uint8_t* out,
                     uint8_t* y_out) {

  uint32_t lambdaBytes = lambda / 8;
  uint32_t x_len       = lambdaBytes * Lenc;
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
              uint8_t* out_zeros_concat = malloc(m / 8);
              out_zeros_concat[0]       = (getBit(out[(ird - 1152) / 8], i) << 7);
              memset(out_zeros_concat + 1, 0, (m / 8) - 1);

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
}

int aes_enc_constraints() {}

// w (extended witness bits), u (vole masks), v (vole tags), in (faest public-key input block
// bits), out (faest public-key output block bits), chal (quicksilver challenge), a_tilde and
// b_tilde (quicksilver response)
void aes_prove(uint8_t* w, uint8_t* u, uint8_t** v, uint8_t* in, uint8_t* out, uint8_t* chal,
               uint32_t lambda, uint32_t tau, uint32_t l, uint8_t* a_tilde, uint8_t* b_tilde) {

  uint32_t lByte      = l / 8;
  uint32_t lambdaByte = lambda / 8;

  // Step: 1
  // TODO: Is to field, equivalent to the function bf_load ??
  // TODO: Are we really passing it bit by bit to ToField function ?
  uint8_t* w_field = malloc(lByte);
  //   bf8_t* bf_w = malloc(lByte);
  //   for (uint32_t i = 0; lByte < l; i++) {
  //     bf_w[i] = bf8_load(w[i]);
  //   }

  // Step: 2
  // TODO: Is to field, equivalent to the function bf_load ??
  // TODO: Are we really passing it bit by bit to ToField function ?
  // TODO: Unsure what to do here
  uint8_t v_field = malloc(lByte + lambdaByte);
  //   for (uint32_t i = 0; i < lByte + lambdaBytes; i++) {
  //   }

  // Step: 3
  uint8_t* in_0  = malloc(16);
  uint8_t* out_0 = malloc(16);
  memcpy(in_0, in, 16);
  memcpy(out_0, out, 16);

  // Step: 4
  // TODO: What is B ?
  uint32_t B     = 2;
  uint8_t* in_1  = malloc(16);
  uint8_t* out_1 = malloc(16);
  if (B == 2) {
    memcpy(in_1, in + 16, 16);
    memcpy(out_1, out + 16, 16);
  }

  // Step: 5
  // TODO: Sampling some element from the field ??
  uint8_t *a0, a1;

  // Step: 6
  // TODO: What is l_ke ??
  uint32_t lByte_ke      = 16;
  uint8_t* w_field_tilde = malloc(lByte_ke);
  memcpy(w_field_tilde, w_field, lByte_ke);
  uint8_t* v_field_tilde = malloc(lByte_ke);
  memcpy(v_field_tilde, v_field, lByte_ke);

  // Step: 7
  uint8_t MKey = 0;
  uint8_t* q;
  uint8_t* delta;
  uint8_t* a0_tilde;
  uint8_t* a1_tilde;
  uint8_t* k;
  uint8_t* vk;
  aes_key_schedule_constrains(lambda, w_field_tilde, v_field_tilde, MKey, q, delta, a0_tilde,
                              a1_tilde, k, vk);

  // Step: 8
  // TODO: resize it before memcpy
  memcpy(a0, a0_tilde, sizeof(a0_tilde));
  memcpy(a1, a1_tilde, sizeof(a1_tilde));

  // Step: 9
  // TODO
  uint32_t lByte_enc;
  memcpy(w_field_tilde, w_field + lByte_ke, lByte_enc);
  memcpy(v_field_tilde, v_field + lByte_ke, lByte_enc);

  // Step: 10
  // TODO
  uint8_t* qk;
  aes_cipher_constrains(lambda, in_0, out_0, w_field_tilde, v_field_tilde, k, vk, MKey, q, qk,
                        delta, a0_tilde, a1_tilde);

  // Step: 11
  // TODO
  memcpy(a0, a0_tilde, sizeof(a0_tilde));
  memcpy(a1, a1_tilde, sizeof(a1_tilde));

  // STep: 12..15
  // TODO
  if (B == 2) {
    memcpy(w_field_tilde, w_field + (lByte_ke + lByte_enc), (lByte) - (lByte_ke + lByte_enc));
    memcpy(v_field_tilde, v_field + (lByte_ke + lByte_enc), (lByte) - (lByte_ke + lByte_enc));
    aes_cipher_constrains(lambda, in_1, out_1, w_field_tilde, v_field_tilde, k, vk, MKey, q, qk,
                          delta, a0_tilde, a1_tilde);
  }

  // Step: 16..17
  for (uint32_t i = lByte; i < lByte + lambdaByte; i++) {
    // TODO: unsure what is happening here...
  }

  // STep: 18
  // TODO: unsure
  uint8_t* us;
  uint8_t* vs;

  // Step: 19..20
  // TODO
  uint8_t* a1_us_concat = malloc(sizeof(a1) + sizeof(us));
  uint8_t* a0_vs_concat = malloc(sizeof(a0) + sizeof(vs));
  switch (lambda) {
  case 256:
    zk_hash_256(a_tilde, r, chal, t, a1_us_concat, ell);
    zk_hash_256(b_tilde, r, chal, t, a0_vs_concat, ell);
    break;
  case 192:
    zk_hash_192(a_tilde, r, chal, t, a1_us_concat, ell);
    zk_hash_192(b_tilde, r, chal, t, a0_vs_concat, ell);
    break;
  default:
    zk_hash_128(a_tilde, r, chal, t, a1_us_concat, ell);
    zk_hash_128(b_tilde, r, chal, t, a0_vs_concat, ell);
    break;
  }
}

bool aes_verify(uint8_t* d, uint8_t* Q, uint8_t* chal_2, uint8_t* chal_3, uint8_t* a_tilde,
                uint8_t* b_tilde, uint8_t* in, uint8_t* out, uint32_t lambda, uint32_t tau,
                uint32_t l, uint32_t k0, uint32_t k1) {

  uint32_t lByte = l / 8;

  // Step: 1
  // TODO: confirm if this is a check or we actually assign it here
  uint32_t t0 = lambda % tau;
  uint32_t t1 = (lambda - (k0 * t0)) / k1;

  // STep: 2
  // TODO: field operation
  uint8_t* delta;

  // Step: 3
  uint8_t* in0 = malloc(16);
  memcpy(in0, in, 16);
  uint8_t* out0 = malloc(16);
  memcpy(out0, out, 16);

  // STep: 4
  // TODO: What is B ?
  uint32_t B = 2;
  uint8_t *in1, out1;
  if (B == 2) {
    in1  = malloc(16);
    out1 = malloc(16);
  }

  // Step: 5..9
  // TODO: from where comes the small b ?
  // TODO:
  uint8_t* b;
  uint32_t kb = 0;
  uint8_t** q;
  for (uint32_t i = 0; i < tau; i++) {
    if (i < t0) {
      b  = 0;
      kb = k0;
    } else {
      b  = 1;
      kb = k1;
    }
    uint8_t* chalout;
    ChalDec(chal_3, i, k0, t0, k1, t1, chalout);
    for (uint32_t j = 0; j < kb; j++) {
      // TODO: check it, very likely it is wrong, need to do bf multiplication I guess
      q[i][j] = q[i][j] ^ (chalout[j] * d[j]);
    }
  }
  // Step: 10..12
  for (uint32_t i = 0; i < lByte + tau; i++) {
    // TODO: ToField problem....
  }
  // TODO: sampling b

  // Step: 13
  // TODO: pass and check for null instead of 0 ?
  // TODO:
  uint32_t Mkey = 1;
  uint32_t lByte_ke;
  uint8_t* q_;
  uint8_t* b1;
  uint8_t* qk;
  memcpy(q_, q, lByte_ke);
  aes_key_schedule_constrains(lambda, NULL, NULL, Mkey, q_, delta, b1, qk);

  // Step: 14
  // TODO
  uint8_t* q_1;
  uint32_t lByte_enc;
  uint8_t* b2;
  memcpy(q_1, q + lByte_ke, lByte_enc);
  aes_cipher_constrains(lambda, in0, out0, NULL, NULL, NULL, NULL, Mkey, q_1, qk, delta, b2);

  // Step: 15..17
  // TODO
  if (B == 1) {
    memcpy(b, b1, sizeof(b1));
    memcpy(b + sizeof(b1), b2, sizeof(b2));
  } else {
    uint8_t* b3;
    uint8_t* q_2;
    memcpy(q_2, q + (lByte_ke + lByte_enc), lByte - (lByte_ke + lByte_enc));
    aes_cipher_constrains(lambda, in1, out1, NULL, NULL, NULL, NULL, Mkey, q_2, qk, delta, b3);
    memcpy(b, b1, sizeof(b1));
    memcpy(b + sizeof(b1), b2, sizeof(b2));
    memcpy(b + sizeof(b1) + sizeof(b2), b3, sizeof(b3));
  }

  // Step: 18
  // TODO : unsure what to do
  uint8_t* qs;

  // Step: 19
  // TODO
  uint8_t* q_tilde;
  uint8_t* b_qs_conat;
  uint8_t *r, t;
  uint32_t ell;
  switch (lambda) {
  case 256:
    zk_hash_256(q_tilde, r, chal_2, t, b_qs_conat, ell);
    break;
  case 192:
    zk_hash_192(q_tilde, r, chal_2, t, b_qs_conat, ell);
    break;
  default:
    zk_hash_128(q_tilde, r, chal_2, t, b_qs_conat, ell);
    break;
  }
  uint8_t* q_check = malloc(sizeof(a_tilde));
  // TODO: its bf operation i guess
  for (uint32_t i = 0; i < sizeof(q_check); i++) {
    q_check[i] = a_tilde[i] + (b_tilde[i] * delta[i]);
  }

  if (memcmp(q_tilde, q_check, sizeof(q_check)) == 0) {
    return true;
  } else {
    return false;
  }
}