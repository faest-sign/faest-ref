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


// Arranged in first occurance in EncSctrnts basis
// ADD ROUND KEY
static void aes_128_add_round_key_prover(bf128_t* out, bf128_t* out_tag, const bf128_t* in, const bf128_t* in_tag, const bf128_t* k, const bf128_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out[i] = bf128_add(in[i], k[i]);
    out_tag[i] = bf128_add(in_tag[i], k_tag[i]);
  }
}
static void aes_128_add_round_key_verifier(bf128_t* out_tag, const bf128_t* in_tag, const bf128_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out_tag[i] = bf128_add(in_tag[i], k_tag[i]);
  }
}

static void aes_192_add_round_key_prover(bf192_t* out, bf192_t* out_tag, const bf192_t* in, const bf192_t* in_tag, const bf192_t* k, const bf192_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out[i] = bf192_add(in[i], k[i]);
    out_tag[i] = bf192_add(in_tag[i], k_tag[i]);
  }
}
static void aes_192_add_round_key_verifier(bf192_t* out_tag, const bf192_t* in_tag, const bf192_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out_tag[i] = bf192_add(in_tag[i], k_tag[i]);
  }
}

static void aes_256_add_round_key_prover(bf256_t* out, bf256_t* out_tag, const bf256_t* in, const bf256_t* in_tag, const bf256_t* k, const bf256_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out[i] = bf256_add(in[i], k[i]);
    out_tag[i] = bf256_add(in_tag[i], k_tag[i]);
  }
}
static void aes_256_add_round_key_verifier(bf256_t* out_tag, const bf256_t* in_tag, const bf256_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out_tag[i] = bf256_add(in_tag[i], k_tag[i]);
  }
}

// F256/F2.CONJUGATES
static void aes_128_f256_f2_conjugates_1(bf128_t* y, const uint8_t* state) {
  for (unsigned int i = 0; i != 4; ++i) {
    uint8_t x0 = state[i];
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf128_byte_combine_bits(x0);
      x0           = bits_square(x0);
    }
    y[i * 8 + 7] = bf128_byte_combine_bits(x0);
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
    y[i * 8 + 7] = bf128_byte_combine(x);
  }
}

static void aes_192_f256_f2_conjugates_1(bf192_t* y, const uint8_t* state) {
  for (unsigned int i = 0; i != 4; ++i) {
    uint8_t x0 = state[i];
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf192_byte_combine_bits(x0);
      x0           = bits_square(x0);
    }
    y[i * 8 + 7] = bf192_byte_combine_bits(x0);
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
    y[i * 8 + 7] = bf192_byte_combine(x);
  }
}

static void aes_256_f256_f2_conjugates_1(bf256_t* y, const uint8_t* state) {
  for (unsigned int i = 0; i != 4; ++i) {
    uint8_t x0 = state[i];
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf256_byte_combine_bits(x0);
      x0           = bits_square(x0);
    }
    y[i * 8 + 7] = bf256_byte_combine_bits(x0);
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
    y[i * 8 + 7] = bf256_byte_combine(x);
  }
}

static void aes_em_192_f256_f2_conjugates_1(bf192_t* y, const uint8_t* state) {
  for (unsigned int i = 0; i != 6; ++i) {
    uint8_t x0 = state[i];
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf192_byte_combine_bits(x0);
      x0           = bits_square(x0);
    }
    y[i * 8 + 7] = bf192_byte_combine_bits(x0);
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
    y[i * 8 + 7] = bf192_byte_combine(x);
  }
}
static void aes_em_256_f256_f2_conjugates_1(bf256_t* y, const uint8_t* state) {
  for (unsigned int i = 0; i != 8; ++i) {
    uint8_t x0 = state[i];
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf256_byte_combine_bits(x0);
      x0           = bits_square(x0);
    }
    y[i * 8 + 7] = bf256_byte_combine_bits(x0);
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

// INV NORM TO CONJUGATES
static void aes_128_inv_norm_to_conjugates_prover(bf128_t* y, uint8_t x) {
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

// INV NORM CONSTRAINTS
void aes_128_inv_norm_constraints_prover(bf128_t* z0, bf128_t* z1, const bf128_t* state_bits, const bf128_t* state_bits_tag, const uint8_t* y, const bf128_t* y_tag) {
  
    z0[0] = bf128_add(
              bf128_mul(bf128_mul(
                                bf128_byte_combine_bits(y), 
                                state_bits[1]), 
                        state_bits[4]),
              state_bits[0]);

    z1[0] = bf128_add(
              bf128_mul(bf128_mul(
                                bf128_byte_combine_bits(y_tag), 
                                state_bits_tag[1]), 
                        state_bits_tag[4]),
              state_bits_tag[0]);

}
void aes_128_inv_norm_constraints_verifier(bf128_t* z1, const bf128_t* state_bits_key, const bf128_t* y_key) {
  z1[0] = bf128_add(
              bf128_mul(bf128_mul(
                                bf128_byte_combine_bits(y_key), 
                                state_bits_key[1]), 
                        state_bits_key[4]),
              state_bits_key[0]);
}

void aes_192_inv_norm_constraints_prover(bf192_t* z0, bf192_t* z1, const bf192_t* state_bits, const bf192_t* state_bits_tag, const uint8_t* y, const bf192_t* y_tag) {
  
    z0[0] = bf192_add(
              bf192_mul(bf192_mul(
                                bf192_byte_combine_bits(y), 
                                state_bits[1]), 
                        state_bits[4]),
              state_bits[0]);

    z1[0] = bf192_add(
              bf192_mul(bf192_mul(
                                bf192_byte_combine_bits(y_tag), 
                                state_bits_tag[1]), 
                        state_bits_tag[4]),
              state_bits_tag[0]);

}
void aes_192_inv_norm_constraints_verifier(bf192_t* z1, const bf192_t* state_bits_key, const bf192_t* y_key) {
  z1[0] = bf192_add(
              bf192_mul(bf192_mul(
                                bf192_byte_combine_bits(y_key), 
                                state_bits_key[1]), 
                        state_bits_key[4]),
              state_bits_key[0]);
}

void aes_256_inv_norm_constraints_prover(bf256_t* z0, bf256_t* z1, const bf256_t* state_bits, const bf256_t* state_bits_tag, const uint8_t* y, const bf256_t* y_tag) {
  
    z0[0] = bf256_add(
              bf256_mul(bf256_mul(
                                bf256_byte_combine_bits(y), 
                                state_bits[1]), 
                        state_bits[4]),
              state_bits[0]);

    z1[0] = bf256_add(
              bf256_mul(bf256_mul(
                                bf256_byte_combine_bits(y_tag), 
                                state_bits_tag[1]), 
                        state_bits_tag[4]),
              state_bits_tag[0]);

}
void aes_256_inv_norm_constraints_verifier(bf256_t* z1, const bf256_t* state_bits_key, const bf256_t* y_key) {
  z1[0] = bf256_add(
              bf256_mul(bf256_mul(
                                bf256_byte_combine_bits(y_key), 
                                state_bits_key[1]), 
                        state_bits_key[4]),
              state_bits_key[0]);
}

// STATE TO BYTES
void aes_128_state_to_bytes_prover(bf128_t* out, bf128_t* out_tag, const uint8_t* state, const bf128_t* state_tag) {
  for (unsigned int i = 0; i < 16; i++) {
    out[i] = bf128_byte_combine_bits(state[i]);
    out_tag[i] = bf128_byte_combine(state_tag + i*8);
  }
}
void aes_128_state_to_bytes_verifier(bf128_t* out_key, const bf128_t* state_key) {
  for (unsigned int i = 0; i < 16; i++) {
    out_key[i] = bf128_byte_combine(state_key + i*8);
  }
}

void aes_192_state_to_bytes_prover(bf192_t* out, bf192_t* out_tag, const uint8_t* state, const bf192_t* state_tag) {
  for (unsigned int i = 0; i < 16; i++) {
    out[i] = bf192_byte_combine_bits(state[i]);
    out_tag[i] = bf192_byte_combine(state_tag + i*8);
  }
}
void aes_192_state_to_bytes_verifier(bf192_t* out_key, const bf192_t* state_key) {
  for (unsigned int i = 0; i < 16; i++) {
    out_key[i] = bf192_byte_combine(state_key + i*8);
  }
}

void aes_256_state_to_bytes_prover(bf256_t* out, bf256_t* out_tag, const uint8_t* state, const bf256_t* state_tag) {
  for (unsigned int i = 0; i < 16; i++) {
    out[i] = bf256_byte_combine_bits(state[i]);
    out_tag[i] = bf256_byte_combine(state_tag + i*8);
  }
}
void aes_256_state_to_bytes_verifier(bf256_t* out_key, const bf256_t* state_key) {
  for (unsigned int i = 0; i < 16; i++) {
    out_key[i] = bf256_byte_combine(state_key + i*8);
  }
}

// SBOX AFFINE
static void aes_128_sbox_affine_prover(bf128_t* out, bf128_t* out_tag, const bf128_t* in, const bf128_t* in_tag, bool dosq, const faest_paramset_t* params) {

  unsigned int Nst_bytes = params->faest_param.lambda/8;

  bf128_t C[9];
  uint8_t x[9] = {0x05,0x09,0xf9,0x25,0xf4,0x01,0xb5,0x8f,0x63};

  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf128_byte_combine_bits(x[i]);
      C[i] = bf128_mul(C[i], C[i]);
    }
  } else {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf128_byte_combine_bits(x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 9; Cidx++) {
      out[i] += bf128_mul(C[Cidx], bf128_mul(in[i], in[i]));
      out_tag[i] += bf128_mul(C[Cidx], bf128_mul(in_tag[i], in_tag[i]));
    }
  }
}
static void aes_128_sbox_affine_verify(bf128_t* out_tag, const bf128_t* in_tag, bool dosq, const faest_paramset_t* params) {

  unsigned int Nst_bytes = params->faest_param.lambda/8;

  bf128_t C[9];
  uint8_t x[9] = {0x05,0x09,0xf9,0x25,0xf4,0x01,0xb5,0x8f,0x63};

  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf128_byte_combine_bits(x[i]);
      C[i] = bf128_mul(C[i], C[i]);
    }
  } else {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf128_byte_combine_bits(x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 9; Cidx++) {
      out_tag[i] += bf128_mul(C[Cidx], bf128_mul(in_tag[i], in_tag[i]));
    }
  }
}

static void aes_192_sbox_affine_prover(bf192_t* out, bf192_t* out_tag, const bf192_t* in, const bf192_t* in_tag, bool dosq, const faest_paramset_t* params) {

  unsigned int Nst_bytes = params->faest_param.lambda/8;

  bf192_t C[9];
  uint8_t x[9] = {0x05,0x09,0xf9,0x25,0xf4,0x01,0xb5,0x8f,0x63};

  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf192_byte_combine_bits(x[i]);
      C[i] = bf192_mul(C[i], C[i]);
    }
  } else {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf192_byte_combine_bits(x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 9; Cidx++) {
      out[i] += bf192_mul(C[Cidx], bf192_mul(in[i], in[i]));
      out_tag[i] += bf192_mul(C[Cidx], bf192_mul(in_tag[i], in_tag[i]));
    }
  }
}
static void aes_192_sbox_affine_verify(bf192_t* out_tag, const bf192_t* in_tag, bool dosq, const faest_paramset_t* params) {

  unsigned int Nst_bytes = params->faest_param.lambda/8;

  bf192_t C[9];
  uint8_t x[9] = {0x05,0x09,0xf9,0x25,0xf4,0x01,0xb5,0x8f,0x63};

  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf192_byte_combine_bits(x[i]);
      C[i] = bf192_mul(C[i], C[i]);
    }
  } else {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf192_byte_combine_bits(x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 9; Cidx++) {
      out_tag[i] += bf192_mul(C[Cidx], bf192_mul(in_tag[i], in_tag[i]));
    }
  }
}

static void aes_256_sbox_affine_prover(bf256_t* out, bf256_t* out_tag, const bf256_t* in, const bf256_t* in_tag, bool dosq, const faest_paramset_t* params) {

  unsigned int Nst_bytes = params->faest_param.lambda/8;

  bf256_t C[9];
  uint8_t x[9] = {0x05,0x09,0xf9,0x25,0xf4,0x01,0xb5,0x8f,0x63};

  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf256_byte_combine_bits(x[i]);
      C[i] = bf256_mul(C[i], C[i]);
    }
  } else {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf256_byte_combine_bits(x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 9; Cidx++) {
      out[i] += bf256_mul(C[Cidx], bf256_mul(in[i], in[i]));
      out_tag[i] += bf256_mul(C[Cidx], bf256_mul(in_tag[i], in_tag[i]));
    }
  }
}
static void aes_256_sbox_affine_verify(bf256_t* out_tag, const bf256_t* in_tag, bool dosq, const faest_paramset_t* params) {

  unsigned int Nst_bytes = params->faest_param.lambda/8;

  bf256_t C[9];
  uint8_t x[9] = {0x05,0x09,0xf9,0x25,0xf4,0x01,0xb5,0x8f,0x63};

  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf256_byte_combine_bits(x[i]);
      C[i] = bf256_mul(C[i], C[i]);
    }
  } else {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf256_byte_combine_bits(x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 9; Cidx++) {
      out_tag[i] += bf256_mul(C[Cidx], bf256_mul(in_tag[i], in_tag[i]));
    }
  }
}

// SHIFT ROWS
static void aes_128_shiftrows_prover(uint8_t* out, bf128_t* out_tag, const uint8_t* in, const bf128_t* in_tag, const faest_paramset_t* params) {
  unsigned int Nst = 4;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if (r <= 1) {
        out[4*c + r] = in[4*((c + r) % 4) + r];
        out_tag[4*c + r] = in_tag[(4*((c + r) % 4) + r)];
      } 
      else {
        out[4*c + r] = in[4*((c + r) % 4) + r];
        out_tag[4*c + r] = in_tag[(4*((c + r + 1) % 4) + r)];
      }
    }
  }
}
static void aes_128_shiftrows_verifier(bf128_t* out_tag, const bf128_t* in_tag, const faest_paramset_t* params) {
  unsigned int Nst = 4;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if (r <= 1) {
        out_tag[4*c + r] = in_tag[(4*((c + r) % 4) + r)];
      } 
      else {
        out_tag[4*c + r] = in_tag[(4*((c + r + 1) % 4) + r)];
      }
    }
  }
}

static void aes_192_shiftrows_prover(uint8_t* out, bf192_t* out_tag, const uint8_t* in, const bf192_t* in_tag, const faest_paramset_t* params) {
  unsigned int Nst = 4;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if (r <= 1) {
        out[4*c + r] = in[4*((c + r) % 4) + r];
        out_tag[4*c + r] = in_tag[(4*((c + r) % 4) + r)];
      } 
      else {
        out[4*c + r] = in[4*((c + r) % 4) + r];
        out_tag[4*c + r] = in_tag[(4*((c + r + 1) % 4) + r)];
      }
    }
  }
}
static void aes_192_shiftrows_verifier(bf192_t* out_tag, const bf192_t* in_tag, const faest_paramset_t* params) {
  unsigned int Nst = 4;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if (r <= 1) {
        out_tag[4*c + r] = in_tag[(4*((c + r) % 4) + r)];
      } 
      else {
        out_tag[4*c + r] = in_tag[(4*((c + r + 1) % 4) + r)];
      }
    }
  }
}

static void aes_256_shiftrows_prover(uint8_t* out, bf256_t* out_tag, const uint8_t* in, const bf256_t* in_tag, const faest_paramset_t* params) {
  unsigned int Nst = 4;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if (r <= 1) {
        out[4*c + r] = in[4*((c + r) % 4) + r];
        out_tag[4*c + r] = in_tag[(4*((c + r) % 4) + r)];
      } 
      else {
        out[4*c + r] = in[4*((c + r) % 4) + r];
        out_tag[4*c + r] = in_tag[(4*((c + r + 1) % 4) + r)];
      }
    }
  }
}
static void aes_256_shiftrows_verifier(bf256_t* out_tag, const bf256_t* in_tag, const faest_paramset_t* params) {
  unsigned int Nst = 4;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if (r <= 1) {
        out_tag[4*c + r] = in_tag[(4*((c + r) % 4) + r)];
      } 
      else {
        out_tag[4*c + r] = in_tag[(4*((c + r + 1) % 4) + r)];
      }
    }
  }
}

// MIX COLOUMNS
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
static void aes_128_mix_coloumns_verifier(bf128_t* y_key, const uint8_t* in_key, bool dosq, const faest_paramset_t* params) {
  
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
      bf128_t tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i0]), v2);
    bf128_t tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i1]), v3);
    bf128_t tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i2]), v1);
    bf128_t tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i3]), v1);
    y_key[i0] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    // ::8
    tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i0]), v1);
    tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i1]), v2);
    tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i2]), v3);
    tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i3]), v1);
    y_key[i1] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    // ::9
    tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i0]), v1);
    tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i1]), v1);
    tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i2]), v2);
    tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i3]), v3);
    y_key[i2] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    // ::10
    tmp1_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i0]), v1);
    tmp2_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i1]), v1);
    tmp3_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i2]), v2);
    tmp4_tag = bf128_mul(bf128_byte_combine_bits_sq(in_key[i3]), v3);
    y_key[i3] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));
  
  }

}

static void aes_192_mix_coloumns_prover(bf192_t* y, bf192_t* y_tag, const uint8_t* in, const uint8_t* in_tag, bool dosq, const faest_paramset_t* params) {
  
  unsigned int Nst = 3;
  
  //  ::2-4
  bf192_t v1 = bf192_byte_combine_bits(1);
  bf192_t v2 = bf192_byte_combine_bits(2);
  bf192_t v3 = bf192_byte_combine_bits(3);
  if (dosq) {
    v1 = bf192_mul(v1, v1);
    v2 = bf192_mul(v2, v2);
    v3 = bf192_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {

    unsigned int i0 = 4*c;
    unsigned int i1 = 4*c + 1;
    unsigned int i2 = 4*c + 2;
    unsigned int i3 = 4*c + 3;

    // ::7
    bf192_t tmp1_val = bf192_mul(bf192_byte_combine_bits_sq(in[i0]), v2);
    bf192_t tmp2_val = bf192_mul(bf192_byte_combine_bits_sq(in[i1]), v3);
    bf192_t tmp3_val = bf192_mul(bf192_byte_combine_bits_sq(in[i2]), v1);
    bf192_t tmp4_val = bf192_mul(bf192_byte_combine_bits_sq(in[i3]), v1);
    y[i0] = bf192_add(bf192_add(tmp1_val, tmp2_val), bf192_add(tmp3_val, tmp4_val));

    bf192_t tmp1_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i0]), v2);
    bf192_t tmp2_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i1]), v3);
    bf192_t tmp3_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i2]), v1);
    bf192_t tmp4_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i3]), v1);
    y_tag[i0] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));

    // ::8
    tmp1_val = bf192_mul(bf192_byte_combine_bits_sq(in[i0]), v1);
    tmp2_val = bf192_mul(bf192_byte_combine_bits_sq(in[i1]), v2);
    tmp3_val = bf192_mul(bf192_byte_combine_bits_sq(in[i2]), v3);
    tmp4_val = bf192_mul(bf192_byte_combine_bits_sq(in[i3]), v1);
    y[i1] = bf192_add(bf192_add(tmp1_val, tmp2_val), bf192_add(tmp3_val, tmp4_val));

    tmp1_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i1]), v2);
    tmp3_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i2]), v3);
    tmp4_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i3]), v1);
    y_tag[i1] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));

    // ::9
    tmp1_val = bf192_mul(bf192_byte_combine_bits_sq(in[i0]), v1);
    tmp2_val = bf192_mul(bf192_byte_combine_bits_sq(in[i1]), v1);
    tmp3_val = bf192_mul(bf192_byte_combine_bits_sq(in[i2]), v2);
    tmp4_val = bf192_mul(bf192_byte_combine_bits_sq(in[i3]), v3);
    y[i2] = bf192_add(bf192_add(tmp1_val, tmp2_val), bf192_add(tmp3_val, tmp4_val));

    tmp1_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i1]), v1);
    tmp3_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i2]), v2);
    tmp4_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i3]), v3);
    y_tag[i2] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));

    // ::10
    tmp1_val = bf192_mul(bf192_byte_combine_bits_sq(in[i0]), v3);
    tmp2_val = bf192_mul(bf192_byte_combine_bits_sq(in[i1]), v1);
    tmp3_val = bf192_mul(bf192_byte_combine_bits_sq(in[i2]), v1);
    tmp4_val = bf192_mul(bf192_byte_combine_bits_sq(in[i3]), v2);
    y[i3] = bf192_add(bf192_add(tmp1_val, tmp2_val), bf192_add(tmp3_val, tmp4_val));

    tmp1_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i1]), v1);
    tmp3_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i2]), v2);
    tmp4_tag = bf192_mul(bf192_byte_combine_bits_sq(in_tag[i3]), v3);
    y_tag[i3] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));
  
  }

}
static void aes_192_mix_coloumns_verifier(bf192_t* y_key, const uint8_t* in_key, bool dosq, const faest_paramset_t* params) {
  
  unsigned int Nst = 3;
  
  //  ::2-4
  bf192_t v1 = bf192_byte_combine_bits(1);
  bf192_t v2 = bf192_byte_combine_bits(2);
  bf192_t v3 = bf192_byte_combine_bits(3);
  if (dosq) {
    v1 = bf192_mul(v1, v1);
    v2 = bf192_mul(v2, v2);
    v3 = bf192_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {

    unsigned int i0 = 4*c;
    unsigned int i1 = 4*c + 1;
    unsigned int i2 = 4*c + 2;
    unsigned int i3 = 4*c + 3;

    // ::7
      bf192_t tmp1_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i0]), v2);
    bf192_t tmp2_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i1]), v3);
    bf192_t tmp3_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i2]), v1);
    bf192_t tmp4_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i3]), v1);
    y_key[i0] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));

    // ::8
    tmp1_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i0]), v1);
    tmp2_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i1]), v2);
    tmp3_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i2]), v3);
    tmp4_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i3]), v1);
    y_key[i1] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));

    // ::9
    tmp1_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i0]), v1);
    tmp2_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i1]), v1);
    tmp3_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i2]), v2);
    tmp4_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i3]), v3);
    y_key[i2] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));

    // ::10
    tmp1_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i0]), v1);
    tmp2_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i1]), v1);
    tmp3_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i2]), v2);
    tmp4_tag = bf192_mul(bf192_byte_combine_bits_sq(in_key[i3]), v3);
    y_key[i3] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));
  
  }

}

static void aes_256_mix_coloumns_prover(bf256_t* y, bf256_t* y_tag, const uint8_t* in, const uint8_t* in_tag, bool dosq, const faest_paramset_t* params) {
  
  unsigned int Nst = 3;
  
  //  ::2-4
  bf256_t v1 = bf256_byte_combine_bits(1);
  bf256_t v2 = bf256_byte_combine_bits(2);
  bf256_t v3 = bf256_byte_combine_bits(3);
  if (dosq) {
    v1 = bf256_mul(v1, v1);
    v2 = bf256_mul(v2, v2);
    v3 = bf256_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {

    unsigned int i0 = 4*c;
    unsigned int i1 = 4*c + 1;
    unsigned int i2 = 4*c + 2;
    unsigned int i3 = 4*c + 3;

    // ::7
    bf256_t tmp1_val = bf256_mul(bf256_byte_combine_bits_sq(in[i0]), v2);
    bf256_t tmp2_val = bf256_mul(bf256_byte_combine_bits_sq(in[i1]), v3);
    bf256_t tmp3_val = bf256_mul(bf256_byte_combine_bits_sq(in[i2]), v1);
    bf256_t tmp4_val = bf256_mul(bf256_byte_combine_bits_sq(in[i3]), v1);
    y[i0] = bf256_add(bf256_add(tmp1_val, tmp2_val), bf256_add(tmp3_val, tmp4_val));

    bf256_t tmp1_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i0]), v2);
    bf256_t tmp2_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i1]), v3);
    bf256_t tmp3_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i2]), v1);
    bf256_t tmp4_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i3]), v1);
    y_tag[i0] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));

    // ::8
    tmp1_val = bf256_mul(bf256_byte_combine_bits_sq(in[i0]), v1);
    tmp2_val = bf256_mul(bf256_byte_combine_bits_sq(in[i1]), v2);
    tmp3_val = bf256_mul(bf256_byte_combine_bits_sq(in[i2]), v3);
    tmp4_val = bf256_mul(bf256_byte_combine_bits_sq(in[i3]), v1);
    y[i1] = bf256_add(bf256_add(tmp1_val, tmp2_val), bf256_add(tmp3_val, tmp4_val));

    tmp1_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i1]), v2);
    tmp3_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i2]), v3);
    tmp4_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i3]), v1);
    y_tag[i1] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));

    // ::9
    tmp1_val = bf256_mul(bf256_byte_combine_bits_sq(in[i0]), v1);
    tmp2_val = bf256_mul(bf256_byte_combine_bits_sq(in[i1]), v1);
    tmp3_val = bf256_mul(bf256_byte_combine_bits_sq(in[i2]), v2);
    tmp4_val = bf256_mul(bf256_byte_combine_bits_sq(in[i3]), v3);
    y[i2] = bf256_add(bf256_add(tmp1_val, tmp2_val), bf256_add(tmp3_val, tmp4_val));

    tmp1_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i1]), v1);
    tmp3_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i2]), v2);
    tmp4_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i3]), v3);
    y_tag[i2] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));

    // ::10
    tmp1_val = bf256_mul(bf256_byte_combine_bits_sq(in[i0]), v3);
    tmp2_val = bf256_mul(bf256_byte_combine_bits_sq(in[i1]), v1);
    tmp3_val = bf256_mul(bf256_byte_combine_bits_sq(in[i2]), v1);
    tmp4_val = bf256_mul(bf256_byte_combine_bits_sq(in[i3]), v2);
    y[i3] = bf256_add(bf256_add(tmp1_val, tmp2_val), bf256_add(tmp3_val, tmp4_val));

    tmp1_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i0]), v1);
    tmp2_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i1]), v1);
    tmp3_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i2]), v2);
    tmp4_tag = bf256_mul(bf256_byte_combine_bits_sq(in_tag[i3]), v3);
    y_tag[i3] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));
  
  }

}
static void aes_256_mix_coloumns_verifier(bf256_t* y_key, const uint8_t* in_key, bool dosq, const faest_paramset_t* params) {
  
  unsigned int Nst = 3;
  
  //  ::2-4
  bf256_t v1 = bf256_byte_combine_bits(1);
  bf256_t v2 = bf256_byte_combine_bits(2);
  bf256_t v3 = bf256_byte_combine_bits(3);
  if (dosq) {
    v1 = bf256_mul(v1, v1);
    v2 = bf256_mul(v2, v2);
    v3 = bf256_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {

    unsigned int i0 = 4*c;
    unsigned int i1 = 4*c + 1;
    unsigned int i2 = 4*c + 2;
    unsigned int i3 = 4*c + 3;

    // ::7
      bf256_t tmp1_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i0]), v2);
    bf256_t tmp2_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i1]), v3);
    bf256_t tmp3_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i2]), v1);
    bf256_t tmp4_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i3]), v1);
    y_key[i0] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));

    // ::8
    tmp1_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i0]), v1);
    tmp2_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i1]), v2);
    tmp3_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i2]), v3);
    tmp4_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i3]), v1);
    y_key[i1] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));

    // ::9
    tmp1_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i0]), v1);
    tmp2_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i1]), v1);
    tmp3_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i2]), v2);
    tmp4_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i3]), v3);
    y_key[i2] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));

    // ::10
    tmp1_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i0]), v1);
    tmp2_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i1]), v1);
    tmp3_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i2]), v2);
    tmp4_tag = bf256_mul(bf256_byte_combine_bits_sq(in_key[i3]), v3);
    y_key[i3] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));
  
  }

}

// ADD ROUND KEY BYTES
static void aes_128_add_round_key_bytes_prover(bf128_t* out, bf128_t* out_tag, const bf128_t* in, const bf128_t* in_tag, const bf128_t* k, const bf128_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out[i] = bf128_add(in[i], k[i]);
    out_tag[i] = bf128_add(in_tag[i], k_tag[i]);
  }
}
static void aes_128_add_round_key_bytes_verifier(bf128_t* out_tag, const bf128_t* in_tag, const bf128_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out_tag[i] = bf128_add(in_tag[i], k_tag[i]);
  }
}

static void aes_192_add_round_key_bytes_prover(bf192_t* out, bf192_t* out_tag, const bf192_t* in, const bf192_t* in_tag, const bf192_t* k, const bf192_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out[i] = bf192_add(in[i], k[i]);
    out_tag[i] = bf192_add(in_tag[i], k_tag[i]);
  }
}
static void aes_192_add_round_key_bytes_verifier(bf192_t* out_tag, const bf192_t* in_tag, const bf192_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out_tag[i] = bf192_add(in_tag[i], k_tag[i]);
  }
}

static void aes_256_add_round_key_bytes_prover(bf256_t* out, bf256_t* out_tag, const bf256_t* in, const bf256_t* in_tag, const bf256_t* k, const bf256_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out[i] = bf256_add(in[i], k[i]);
    out_tag[i] = bf256_add(in_tag[i], k_tag[i]);
  }
}
static void aes_256_add_round_key_bytes_verifier(bf256_t* out_tag, const bf256_t* in_tag, const bf256_t* k_tag, const faest_paramset_t* params) {

  unsigned int Nst = 4;
  unsigned int Nstbits = Nst*32;

  for (unsigned int i = 0; i < Nstbits/8; i++) {
    out_tag[i] = bf256_add(in_tag[i], k_tag[i]);
  }
}

// INVERSE SHIFT ROWS
static void aes_128_inverse_shiftrows_prover(uint8_t* out, bf128_t* out_tag, const uint8_t* in, const bf128_t* in_tag, const faest_paramset_t* params) {
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

      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out[8*(4*c + r) + byte_idx] = in[8*i + byte_idx];
        out_tag[8*(4*c + r) + byte_idx] = in_tag[8*i + byte_idx];
      }
    }
  }
}
static void aes_128_inverse_shiftrows_verifier(bf128_t* out_tag, const bf128_t* in_tag, const faest_paramset_t* params) {
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
      
      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out_tag[8*(4*c + r) + byte_idx] = in_tag[8*i + byte_idx];
      }
    }
  }
}

static void aes_192_inverse_shiftrows_prover(uint8_t* out, bf192_t* out_tag, const uint8_t* in, const bf192_t* in_tag, const faest_paramset_t* params) {
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

      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out[8*(4*c + r) + byte_idx] = in[8*i + byte_idx];
        out_tag[8*(4*c + r) + byte_idx] = in_tag[8*i + byte_idx];
      }
    }
  }
}
static void aes_192_inverse_shiftrows_verifier(bf192_t* out_tag, const bf192_t* in_tag, const faest_paramset_t* params) {
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
      
      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out_tag[8*(4*c + r) + byte_idx] = in_tag[8*i + byte_idx];
      }
    }
  }
}

static void aes_256_inverse_shiftrows_prover(uint8_t* out, bf256_t* out_tag, const uint8_t* in, const bf256_t* in_tag, const faest_paramset_t* params) {
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

      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out[8*(4*c + r) + byte_idx] = in[8*i + byte_idx];
        out_tag[8*(4*c + r) + byte_idx] = in_tag[8*i + byte_idx];
      }
    }
  }
}
static void aes_256_inverse_shiftrows_verifier(bf256_t* out_tag, const bf256_t* in_tag, const faest_paramset_t* params) {
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
      
      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out_tag[8*(4*c + r) + byte_idx] = in_tag[8*i + byte_idx];
      }
    }
  }
}

// INVERSE AFFINE
// DONE: Looks good
static void aes_128_inverse_affine_prover(uint8_t* y, bf128_t* y_tag, const uint8_t* x, const bf128_t* x_tag, const faest_paramset_t* params) {
  
  unsigned int state_size_bytes = 16;

  for (unsigned int i = 0; i < state_size_bytes; i++) {
    uint8_t x_bits[8];
    bf128_t x_bits_tag[8];
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_bits[bit_i] = (x[i] >> bit_i) & 1;
      x_bits_tag[bit_i] = x_tag[i*8 + bit_i];
    }

    uint8_t y_bits[8];
    bf128_t y_bits_tag[8];
    unsigned int c = 0;
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      if (bit_i == 0 || bit_i == 2) {
        c = 1;
      } else {
        c = 0;
      }
      bf128_t c_tag;
      constant_to_vole_128_prover(&c_tag, c);
      y_bits[bit_i] = x_bits[(bit_i - 1 + 8)%8] ^ x_bits[(bit_i - 3 + 8)%8] ^ x_bits[(bit_i - 6 + 8)%8] ^ c;
      y_bits_tag[bit_i] = bf128_add(
                                    bf128_add(x_bits_tag[(bit_i - 1 + 8)%8], x_bits_tag[(bit_i - 3 + 8)%8]),
                                    bf128_add(x_bits_tag[(bit_i - 6 + 8)%8], c_tag));
    }

    for (unsigned int bit_i = 0; bit_i < state_size_bytes; bit_i++) {
      y[i] &= (y_bits[bit_i] << bit_i);
      y_tag[i * 8 + bit_i] = y_bits_tag[bit_i];
    }
  }
}
// DONE: Looks good
static void aes_128_inverse_affine_verifier(bf128_t* y_key, const bf128_t* x_key, bf128_t delta, const faest_paramset_t* params) {
  
  unsigned int state_size_bytes = 16;

  for (unsigned int i = 0; i < state_size_bytes; i++) {
    bf128_t x_bits_key[8];
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_bits_key[bit_i] = x_key[i*8 + bit_i];
    }

    bf128_t y_bits_key[8];
    unsigned int c = 0;
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      if (bit_i == 0 || bit_i == 2) {
        c = 1;
      } else {
        c = 0;
      }
      bf128_t c_tag;
      constant_to_vole_128_verifier(&c_tag, c, delta);
      y_bits_key[bit_i] = bf128_add(
                                    bf128_add(x_bits_key[(bit_i - 1 + 8)%8], x_bits_key[(bit_i - 3 + 8)%8]),
                                    bf128_add(x_bits_key[(bit_i - 6 + 8)%8], c_tag));
    }

    for (unsigned int bit_i = 0; bit_i < state_size_bytes; bit_i++) {
      y_key[i * 8 + bit_i] = y_bits_key[bit_i];
    }
  }
}

static void aes_192_inverse_affine_prover(uint8_t* y, bf192_t* y_tag, const uint8_t* x, const bf192_t* x_tag, const faest_paramset_t* params) {
  
  unsigned int state_size_bytes = 16;

  for (unsigned int i = 0; i < state_size_bytes; i++) {
    uint8_t x_bits[8];
    bf192_t x_bits_tag[8];
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_bits[bit_i] = (x[i] >> bit_i) & 1;
      x_bits_tag[bit_i] = x_tag[i*8 + bit_i];
    }

    uint8_t y_bits[8];
    bf192_t y_bits_tag[8];
    unsigned int c = 0;
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      if (bit_i == 0 || bit_i == 2) {
        c = 1;
      } else {
        c = 0;
      }
      bf192_t c_tag;  // we don't really need to call constant to cole but keeping it for similarity with spec
      constant_to_vole_192_prover(&c_tag, c);
      y_bits[bit_i] = x_bits[(bit_i - 1 + 8)%8] ^ x_bits[(bit_i - 3 + 8)%8] ^ x_bits[(bit_i - 6 + 8)%8] ^ c;
      y_bits_tag[bit_i] = bf192_add(
                                    bf192_add(x_bits_tag[(bit_i - 1 + 8)%8], x_bits_tag[(bit_i - 3 + 8)%8]),
                                    x_bits_tag[(bit_i - 6 + 8)%8] ^ c_tag);
    }

    for (unsigned int bit_i = 0; bit_i < state_size_bytes; bit_i++) {
      y[i] &= (y_bits[bit_i] << bit_i);
      y_tag[i * 8 + bit_i] = y_bits_tag[bit_i];
    }
  }
}
static void aes_192_inverse_affine_verifier(bf192_t* y_key, const bf192_t* x_key, const faest_paramset_t* params) {
  
  unsigned int state_size_bytes = 16;

  for (unsigned int i = 0; i < state_size_bytes; i++) {
    bf192_t x_bits_key[8];
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_bits_key[bit_i] = x_key[i*8 + bit_i];
    }

    bf192_t y_bits_key[8];
    unsigned int c = 0;
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      if (bit_i == 0 || bit_i == 2) {
        c = 1;
      } else {
        c = 0;
      }
      bf192_t c_tag;  // we don't really need to call constant to cole but keeping it for similarity with spec
      constant_to_vole_192_prover(&c_tag, c);
      y_bits_key[bit_i] = bf192_add(
                                    bf192_add(x_bits_key[(bit_i - 1 + 8)%8], x_bits_key[(bit_i - 3 + 8)%8]),
                                    x_bits_key[(bit_i - 6 + 8)%8] ^ c_tag);
    }

    for (unsigned int bit_i = 0; bit_i < state_size_bytes; bit_i++) {
      y_key[i * 8 + bit_i] = y_bits_key[bit_i];
    }
  }
}

static void aes_256_inverse_affine_prover(uint8_t* y, bf256_t* y_tag, const uint8_t* x, const bf256_t* x_tag, const faest_paramset_t* params) {
  
  unsigned int state_size_bytes = 16;

  for (unsigned int i = 0; i < state_size_bytes; i++) {
    uint8_t x_bits[8];
    bf256_t x_bits_tag[8];
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_bits[bit_i] = (x[i] >> bit_i) & 1;
      x_bits_tag[bit_i] = x_tag[i*8 + bit_i];
    }

    uint8_t y_bits[8];
    bf256_t y_bits_tag[8];
    unsigned int c = 0;
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      if (bit_i == 0 || bit_i == 2) {
        c = 1;
      } else {
        c = 0;
      }
      bf256_t c_tag;  // we don't really need to call constant to cole but keeping it for similarity with spec
      constant_to_vole_256_prover(&c_tag, c);
      y_bits[bit_i] = x_bits[(bit_i - 1 + 8)%8] ^ x_bits[(bit_i - 3 + 8)%8] ^ x_bits[(bit_i - 6 + 8)%8] ^ c;
      y_bits_tag[bit_i] = bf256_add(
                                    bf256_add(x_bits_tag[(bit_i - 1 + 8)%8], x_bits_tag[(bit_i - 3 + 8)%8]),
                                    x_bits_tag[(bit_i - 6 + 8)%8] ^ c_tag);
    }

    for (unsigned int bit_i = 0; bit_i < state_size_bytes; bit_i++) {
      y[i] &= (y_bits[bit_i] << bit_i);
      y_tag[i * 8 + bit_i] = y_bits_tag[bit_i];
    }
  }
}
static void aes_256_inverse_affine_verifier(bf256_t* y_key, const bf256_t* x_key, const faest_paramset_t* params) {
  
  unsigned int state_size_bytes = 16;

  for (unsigned int i = 0; i < state_size_bytes; i++) {
    bf256_t x_bits_key[8];
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_bits_key[bit_i] = x_key[i*8 + bit_i];
    }

    bf256_t y_bits_key[8];
    unsigned int c = 0;
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      if (bit_i == 0 || bit_i == 2) {
        c = 1;
      } else {
        c = 0;
      }
      bf256_t c_tag;  // we don't really need to call constant to cole but keeping it for similarity with spec
      constant_to_vole_256_prover(&c_tag, c);
      y_bits_key[bit_i] = bf256_add(
                                    bf256_add(x_bits_key[(bit_i - 1 + 8)%8], x_bits_key[(bit_i - 3 + 8)%8]),
                                    x_bits_key[(bit_i - 6 + 8)%8] ^ c_tag);
    }

    for (unsigned int bit_i = 0; bit_i < state_size_bytes; bit_i++) {
      y_key[i * 8 + bit_i] = y_bits_key[bit_i];
    }
  }
}

// BITWISE MIX COLOUMNS
static void aes_128_bitwise_mix_coloumn_prover(bf128_t* out, bf128_t* out_tag, uint8_t* s, bf128_t* s_tag) {

  unsigned int Nst = 4;

  for (unsigned int c = 0; c < Nst; c++) {

    uint8_t a_bits[4*8];
    bf128_t a_bits_tag[4*8];

    uint8_t b_bits[4*8];
    bf128_t b_bits_tag[4*8];

    // ::1
    for(unsigned int r = 0; r < 4; r++) {
      // :2
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits[r*8 + bit_i] = (s[(32*c+8*r)/8] >> bit_i) & 1;
        a_bits_tag[r*8 + bit_i] = s_tag[32*c+8*r+bit_i];
      }
      // :5
      b_bits[r*8 + 0] = a_bits[r*8 + 7];
      b_bits[r*8 + 1] = a_bits[r*8 + 0] + a_bits[r*8 + 7];
      b_bits[r*8 + 2] = a_bits[r*8 + 1];
      b_bits[r*8 + 3] = a_bits[r*8 + 2] + a_bits[r*8 + 7];
      b_bits[r*8 + 4] = a_bits[r*8 + 3] + a_bits[r*8 + 7];
      b_bits[r*8 + 5] = a_bits[r*8 + 4];
      b_bits[r*8 + 6] = a_bits[r*8 + 5];
      b_bits[r*8 + 7] = a_bits[r*8 + 6];

      b_bits_tag[r*8 + 0] = a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 1] = a_bits_tag[r*8 + 0] + a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 2] = a_bits_tag[r*8 + 1];
      b_bits_tag[r*8 + 3] = a_bits_tag[r*8 + 2] + a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 4] = a_bits_tag[r*8 + 3] + a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 5] = a_bits_tag[r*8 + 4];
      b_bits_tag[r*8 + 6] = a_bits_tag[r*8 + 5];
      b_bits_tag[r*8 + 7] = a_bits_tag[r*8 + 6];

    }

    uint8_t a[4];
    uint8_t b[4];

    bf128_t a_bf[4];
    bf128_t a_tag_bf[4];
    bf128_t b_bf[4];
    bf128_t b_tag_bf[4];
    for (unsigned int round = 0; round < 4; round++) {
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        a[round] &= (a_bits[round*8 + bit_i] << bit_i);
        b[round] &= (b_bits[round*8 + bit_i] << bit_i);
      }
      a_bf[round] = bf128_byte_combine_bits(a[round]);
      b_bf[round] = bf128_byte_combine_bits(b[round]);

      a_tag_bf[round] = bf128_byte_combine(a_bits_tag + round*8);
      b_tag_bf[round] = bf128_byte_combine(b_bits_tag + round*8);
    }

    // ::6-9
    out[c*4] = bf128_add(
                        bf128_add(
                                bf128_add(b_bf[0], a_bf[3]), bf128_add(a_bf[2], b_bf[1])
                                ), a_bf[1]);
    out[c*4 + 1] = bf128_add(
                        bf128_add(
                                bf128_add(b_bf[1], a_bf[0]), bf128_add(a_bf[3], b_bf[2])
                                ), a_bf[2]);
    out[c*4 + 2] = bf128_add(
                        bf128_add(
                                bf128_add(b_bf[2], a_bf[1]), bf128_add(a_bf[0], b_bf[3])
                                ), a_bf[3]);
    out[c*4 + 3] = bf128_add(
                        bf128_add(
                                bf128_add(b_bf[3], a_bf[2]), bf128_add(a_bf[1], b_bf[0])
                                ), a_bf[0]);


    out_tag[c*4] = bf128_add(
                        bf128_add(
                                bf128_add(b_tag_bf[0], a_tag_bf[3]), bf128_add(a_tag_bf[2], b_tag_bf[1])
                                ), a_tag_bf[1]);
    out_tag[c*4 + 1] = bf128_add(
                        bf128_add(
                                bf128_add(b_tag_bf[1], a_tag_bf[0]), bf128_add(a_tag_bf[3], b_tag_bf[2])
                                ), a_tag_bf[2]);
    out_tag[c*4 + 2] = bf128_add(
                        bf128_add(
                                bf128_add(b_tag_bf[2], a_tag_bf[1]), bf128_add(a_tag_bf[0], b_tag_bf[3])
                                ), a_tag_bf[3]);
    out_tag[c*4 + 3] = bf128_add(
                        bf128_add(
                                bf128_add(b_tag_bf[3], a_tag_bf[2]), bf128_add(a_tag_bf[1], b_tag_bf[0])
                                ), a_tag_bf[0]);

  }
}
static void aes_128_bitwise_mix_coloumn_verifier(bf128_t* out_key, bf128_t* s_key) {

  unsigned int Nst = 4;

  for (unsigned int c = 0; c < Nst; c++) {

    bf128_t a_bits_key[4*8];
    bf128_t b_bits_key[4*8];

    // ::1
    for(unsigned int r = 0; r < 4; r++) {
      // :2
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits_key[r*8 + bit_i] = s_key[32*c+8*r+bit_i];
      }
      // :5
      b_bits_key[r*8 + 0] = a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 1] = a_bits_key[r*8 + 0] + a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 2] = a_bits_key[r*8 + 1];
      b_bits_key[r*8 + 3] = a_bits_key[r*8 + 2] + a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 4] = a_bits_key[r*8 + 3] + a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 5] = a_bits_key[r*8 + 4];
      b_bits_key[r*8 + 6] = a_bits_key[r*8 + 5];
      b_bits_key[r*8 + 7] = a_bits_key[r*8 + 6];

    }

    uint8_t a[4];
    uint8_t b[4];

    bf128_t a_bf[4];
    bf128_t a_key_bf[4];
    bf128_t b_bf[4];
    bf128_t b_key_bf[4];
    for (unsigned int round = 0; round < 4; round++) {
      a_bf[round] = bf128_byte_combine_bits(a[round]);
      b_bf[round] = bf128_byte_combine_bits(b[round]);

      a_key_bf[round] = bf128_byte_combine(a_bits_key + round*8);
      b_key_bf[round] = bf128_byte_combine(b_bits_key + round*8);
    }

    // ::6-9
    out_key[c*4] = bf128_add(
                        bf128_add(
                                bf128_add(b_key_bf[0], a_key_bf[3]), bf128_add(a_key_bf[2], b_key_bf[1])
                                ), a_key_bf[1]);
    out_key[c*4 + 1] = bf128_add(
                        bf128_add(
                                bf128_add(b_key_bf[1], a_key_bf[0]), bf128_add(a_key_bf[3], b_key_bf[2])
                                ), a_key_bf[2]);
    out_key[c*4 + 2] = bf128_add(
                        bf128_add(
                                bf128_add(b_key_bf[2], a_key_bf[1]), bf128_add(a_key_bf[0], b_key_bf[3])
                                ), a_key_bf[3]);
    out_key[c*4 + 3] = bf128_add(
                        bf128_add(
                                bf128_add(b_key_bf[3], a_key_bf[2]), bf128_add(a_key_bf[1], b_key_bf[0])
                                ), a_key_bf[0]);

  }
}

static void aes_192_bitwise_mix_coloumn_prover(bf192_t* out, bf192_t* out_tag, uint8_t* s, bf192_t* s_tag) {

  unsigned int Nst = 4;

  for (unsigned int c = 0; c < Nst; c++) {

    uint8_t a_bits[4*8];
    bf192_t a_bits_tag[4*8];

    uint8_t b_bits[4*8];
    bf192_t b_bits_tag[4*8];

    // ::1
    for(unsigned int r = 0; r < 4; r++) {
      // :2
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits[r*8 + bit_i] = (s[(32*c+8*r)/8] >> bit_i) & 1;
        a_bits_tag[r*8 + bit_i] = s_tag[32*c+8*r+bit_i];
      }
      // :5
      b_bits[r*8 + 0] = a_bits[r*8 + 7];
      b_bits[r*8 + 1] = a_bits[r*8 + 0] + a_bits[r*8 + 7];
      b_bits[r*8 + 2] = a_bits[r*8 + 1];
      b_bits[r*8 + 3] = a_bits[r*8 + 2] + a_bits[r*8 + 7];
      b_bits[r*8 + 4] = a_bits[r*8 + 3] + a_bits[r*8 + 7];
      b_bits[r*8 + 5] = a_bits[r*8 + 4];
      b_bits[r*8 + 6] = a_bits[r*8 + 5];
      b_bits[r*8 + 7] = a_bits[r*8 + 6];

      b_bits_tag[r*8 + 0] = a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 1] = a_bits_tag[r*8 + 0] + a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 2] = a_bits_tag[r*8 + 1];
      b_bits_tag[r*8 + 3] = a_bits_tag[r*8 + 2] + a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 4] = a_bits_tag[r*8 + 3] + a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 5] = a_bits_tag[r*8 + 4];
      b_bits_tag[r*8 + 6] = a_bits_tag[r*8 + 5];
      b_bits_tag[r*8 + 7] = a_bits_tag[r*8 + 6];

    }

    uint8_t a[4];
    uint8_t b[4];

    bf192_t a_bf[4];
    bf192_t a_tag_bf[4];
    bf192_t b_bf[4];
    bf192_t b_tag_bf[4];
    for (unsigned int round = 0; round < 4; round++) {
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        a[round] &= (a_bits[round*8 + bit_i] << bit_i);
        b[round] &= (b_bits[round*8 + bit_i] << bit_i);
      }
      a_bf[round] = bf192_byte_combine_bits(a[round]);
      b_bf[round] = bf192_byte_combine_bits(b[round]);

      a_tag_bf[round] = bf192_byte_combine(a_bits_tag + round*8);
      b_tag_bf[round] = bf192_byte_combine(b_bits_tag + round*8);
    }

    // ::6-9
    out[c*4] = bf192_add(
                        bf192_add(
                                bf192_add(b_bf[0], a_bf[3]), bf192_add(a_bf[2], b_bf[1])
                                ), a_bf[1]);
    out[c*4 + 1] = bf192_add(
                        bf192_add(
                                bf192_add(b_bf[1], a_bf[0]), bf192_add(a_bf[3], b_bf[2])
                                ), a_bf[2]);
    out[c*4 + 2] = bf192_add(
                        bf192_add(
                                bf192_add(b_bf[2], a_bf[1]), bf192_add(a_bf[0], b_bf[3])
                                ), a_bf[3]);
    out[c*4 + 3] = bf192_add(
                        bf192_add(
                                bf192_add(b_bf[3], a_bf[2]), bf192_add(a_bf[1], b_bf[0])
                                ), a_bf[0]);


    out_tag[c*4] = bf192_add(
                        bf192_add(
                                bf192_add(b_tag_bf[0], a_tag_bf[3]), bf192_add(a_tag_bf[2], b_tag_bf[1])
                                ), a_tag_bf[1]);
    out_tag[c*4 + 1] = bf192_add(
                        bf192_add(
                                bf192_add(b_tag_bf[1], a_tag_bf[0]), bf192_add(a_tag_bf[3], b_tag_bf[2])
                                ), a_tag_bf[2]);
    out_tag[c*4 + 2] = bf192_add(
                        bf192_add(
                                bf192_add(b_tag_bf[2], a_tag_bf[1]), bf192_add(a_tag_bf[0], b_tag_bf[3])
                                ), a_tag_bf[3]);
    out_tag[c*4 + 3] = bf192_add(
                        bf192_add(
                                bf192_add(b_tag_bf[3], a_tag_bf[2]), bf192_add(a_tag_bf[1], b_tag_bf[0])
                                ), a_tag_bf[0]);

  }
}
static void aes_192_bitwise_mix_coloumn_verifier(bf192_t* out_key, bf192_t* s_key) {

  unsigned int Nst = 4;

  for (unsigned int c = 0; c < Nst; c++) {

    bf192_t a_bits_key[4*8];
    bf192_t b_bits_key[4*8];

    // ::1
    for(unsigned int r = 0; r < 4; r++) {
      // :2
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits_key[r*8 + bit_i] = s_key[32*c+8*r+bit_i];
      }
      // :5
      b_bits_key[r*8 + 0] = a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 1] = a_bits_key[r*8 + 0] + a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 2] = a_bits_key[r*8 + 1];
      b_bits_key[r*8 + 3] = a_bits_key[r*8 + 2] + a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 4] = a_bits_key[r*8 + 3] + a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 5] = a_bits_key[r*8 + 4];
      b_bits_key[r*8 + 6] = a_bits_key[r*8 + 5];
      b_bits_key[r*8 + 7] = a_bits_key[r*8 + 6];

    }

    uint8_t a[4];
    uint8_t b[4];

    bf192_t a_bf[4];
    bf192_t a_key_bf[4];
    bf192_t b_bf[4];
    bf192_t b_key_bf[4];
    for (unsigned int round = 0; round < 4; round++) {
      a_bf[round] = bf192_byte_combine_bits(a[round]);
      b_bf[round] = bf192_byte_combine_bits(b[round]);

      a_key_bf[round] = bf192_byte_combine(a_bits_key + round*8);
      b_key_bf[round] = bf192_byte_combine(b_bits_key + round*8);
    }

    // ::6-9
    out_key[c*4] = bf192_add(
                        bf192_add(
                                bf192_add(b_key_bf[0], a_key_bf[3]), bf192_add(a_key_bf[2], b_key_bf[1])
                                ), a_key_bf[1]);
    out_key[c*4 + 1] = bf192_add(
                        bf192_add(
                                bf192_add(b_key_bf[1], a_key_bf[0]), bf192_add(a_key_bf[3], b_key_bf[2])
                                ), a_key_bf[2]);
    out_key[c*4 + 2] = bf192_add(
                        bf192_add(
                                bf192_add(b_key_bf[2], a_key_bf[1]), bf192_add(a_key_bf[0], b_key_bf[3])
                                ), a_key_bf[3]);
    out_key[c*4 + 3] = bf192_add(
                        bf192_add(
                                bf192_add(b_key_bf[3], a_key_bf[2]), bf192_add(a_key_bf[1], b_key_bf[0])
                                ), a_key_bf[0]);

  }
}

static void aes_256_bitwise_mix_coloumn_prover(bf256_t* out, bf256_t* out_tag, uint8_t* s, bf256_t* s_tag) {

  unsigned int Nst = 4;

  for (unsigned int c = 0; c < Nst; c++) {

    uint8_t a_bits[4*8];
    bf256_t a_bits_tag[4*8];

    uint8_t b_bits[4*8];
    bf256_t b_bits_tag[4*8];

    // ::1
    for(unsigned int r = 0; r < 4; r++) {
      // :2
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits[r*8 + bit_i] = (s[(32*c+8*r)/8] >> bit_i) & 1;
        a_bits_tag[r*8 + bit_i] = s_tag[32*c+8*r+bit_i];
      }
      // :5
      b_bits[r*8 + 0] = a_bits[r*8 + 7];
      b_bits[r*8 + 1] = a_bits[r*8 + 0] + a_bits[r*8 + 7];
      b_bits[r*8 + 2] = a_bits[r*8 + 1];
      b_bits[r*8 + 3] = a_bits[r*8 + 2] + a_bits[r*8 + 7];
      b_bits[r*8 + 4] = a_bits[r*8 + 3] + a_bits[r*8 + 7];
      b_bits[r*8 + 5] = a_bits[r*8 + 4];
      b_bits[r*8 + 6] = a_bits[r*8 + 5];
      b_bits[r*8 + 7] = a_bits[r*8 + 6];

      b_bits_tag[r*8 + 0] = a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 1] = a_bits_tag[r*8 + 0] + a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 2] = a_bits_tag[r*8 + 1];
      b_bits_tag[r*8 + 3] = a_bits_tag[r*8 + 2] + a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 4] = a_bits_tag[r*8 + 3] + a_bits_tag[r*8 + 7];
      b_bits_tag[r*8 + 5] = a_bits_tag[r*8 + 4];
      b_bits_tag[r*8 + 6] = a_bits_tag[r*8 + 5];
      b_bits_tag[r*8 + 7] = a_bits_tag[r*8 + 6];

    }

    uint8_t a[4];
    uint8_t b[4];

    bf256_t a_bf[4];
    bf256_t a_tag_bf[4];
    bf256_t b_bf[4];
    bf256_t b_tag_bf[4];
    for (unsigned int round = 0; round < 4; round++) {
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        a[round] &= (a_bits[round*8 + bit_i] << bit_i);
        b[round] &= (b_bits[round*8 + bit_i] << bit_i);
      }
      a_bf[round] = bf256_byte_combine_bits(a[round]);
      b_bf[round] = bf256_byte_combine_bits(b[round]);

      a_tag_bf[round] = bf256_byte_combine(a_bits_tag + round*8);
      b_tag_bf[round] = bf256_byte_combine(b_bits_tag + round*8);
    }

    // ::6-9
    out[c*4] = bf256_add(
                        bf256_add(
                                bf256_add(b_bf[0], a_bf[3]), bf256_add(a_bf[2], b_bf[1])
                                ), a_bf[1]);
    out[c*4 + 1] = bf256_add(
                        bf256_add(
                                bf256_add(b_bf[1], a_bf[0]), bf256_add(a_bf[3], b_bf[2])
                                ), a_bf[2]);
    out[c*4 + 2] = bf256_add(
                        bf256_add(
                                bf256_add(b_bf[2], a_bf[1]), bf256_add(a_bf[0], b_bf[3])
                                ), a_bf[3]);
    out[c*4 + 3] = bf256_add(
                        bf256_add(
                                bf256_add(b_bf[3], a_bf[2]), bf256_add(a_bf[1], b_bf[0])
                                ), a_bf[0]);


    out_tag[c*4] = bf256_add(
                        bf256_add(
                                bf256_add(b_tag_bf[0], a_tag_bf[3]), bf256_add(a_tag_bf[2], b_tag_bf[1])
                                ), a_tag_bf[1]);
    out_tag[c*4 + 1] = bf256_add(
                        bf256_add(
                                bf256_add(b_tag_bf[1], a_tag_bf[0]), bf256_add(a_tag_bf[3], b_tag_bf[2])
                                ), a_tag_bf[2]);
    out_tag[c*4 + 2] = bf256_add(
                        bf256_add(
                                bf256_add(b_tag_bf[2], a_tag_bf[1]), bf256_add(a_tag_bf[0], b_tag_bf[3])
                                ), a_tag_bf[3]);
    out_tag[c*4 + 3] = bf256_add(
                        bf256_add(
                                bf256_add(b_tag_bf[3], a_tag_bf[2]), bf256_add(a_tag_bf[1], b_tag_bf[0])
                                ), a_tag_bf[0]);

  }
}
static void aes_256_bitwise_mix_coloumn_verifier(bf256_t* out_key, bf256_t* s_key) {

  unsigned int Nst = 4;

  for (unsigned int c = 0; c < Nst; c++) {

    bf256_t a_bits_key[4*8];
    bf256_t b_bits_key[4*8];

    // ::1
    for(unsigned int r = 0; r < 4; r++) {
      // :2
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits_key[r*8 + bit_i] = s_key[32*c+8*r+bit_i];
      }
      // :5
      b_bits_key[r*8 + 0] = a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 1] = a_bits_key[r*8 + 0] + a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 2] = a_bits_key[r*8 + 1];
      b_bits_key[r*8 + 3] = a_bits_key[r*8 + 2] + a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 4] = a_bits_key[r*8 + 3] + a_bits_key[r*8 + 7];
      b_bits_key[r*8 + 5] = a_bits_key[r*8 + 4];
      b_bits_key[r*8 + 6] = a_bits_key[r*8 + 5];
      b_bits_key[r*8 + 7] = a_bits_key[r*8 + 6];

    }

    uint8_t a[4];
    uint8_t b[4];

    bf256_t a_bf[4];
    bf256_t a_key_bf[4];
    bf256_t b_bf[4];
    bf256_t b_key_bf[4];
    for (unsigned int round = 0; round < 4; round++) {
      a_bf[round] = bf256_byte_combine_bits(a[round]);
      b_bf[round] = bf256_byte_combine_bits(b[round]);

      a_key_bf[round] = bf256_byte_combine(a_bits_key + round*8);
      b_key_bf[round] = bf256_byte_combine(b_bits_key + round*8);
    }

    // ::6-9
    out_key[c*4] = bf256_add(
                        bf256_add(
                                bf256_add(b_key_bf[0], a_key_bf[3]), bf256_add(a_key_bf[2], b_key_bf[1])
                                ), a_key_bf[1]);
    out_key[c*4 + 1] = bf256_add(
                        bf256_add(
                                bf256_add(b_key_bf[1], a_key_bf[0]), bf256_add(a_key_bf[3], b_key_bf[2])
                                ), a_key_bf[2]);
    out_key[c*4 + 2] = bf256_add(
                        bf256_add(
                                bf256_add(b_key_bf[2], a_key_bf[1]), bf256_add(a_key_bf[0], b_key_bf[3])
                                ), a_key_bf[3]);
    out_key[c*4 + 3] = bf256_add(
                        bf256_add(
                                bf256_add(b_key_bf[3], a_key_bf[2]), bf256_add(a_key_bf[1], b_key_bf[0])
                                ), a_key_bf[0]);

  }
}

// CONSTANT TO VOLE
// DONE: Looks good
static void constant_to_vole_128_prover(bf128_t* tag, const uint8_t* val) {
  // the val stay the same as the val is a pub const!
  for (unsigned int i = 0; i < 128; i++) {
    tag[i] = bf128_zero();  // for constant values the tag is zero
  }   
}
static void constant_to_vole_128_verifier(bf128_t* key, const uint8_t* val, bf128_t delta) {
  for (unsigned int i = 0; i < 128; i++) {
    key[i] = bf128_mul(bf128_from_bit(get_bit(val[i/8], i%8)), delta);  // multiply delta with each bit
  }  
}
// DONE: Looks good
static void constant_to_vole_192_prover(bf192_t* tag, const uint8_t* val) {
  // the val stay the same as the val is a pub const!
  for (unsigned int i = 0; i < 192; i++) {
    tag[i] = bf192_zero();  // for constant values the tag is zero
  }   
}
static void constant_to_vole_192_verifier(bf192_t* key, const uint8_t* val, bf192_t delta) {
  for (unsigned int i = 0; i < 192; i++) {
    key[i] = bf192_mul(bf192_from_bit(get_bit(val[i/8], i%8)), delta);  // multiply delta with each bit
  }  
}
// DONE: Looks good
static void constant_to_vole_256_prover(bf256_t* tag, const uint8_t* val) {
  // the val stay the same as the val is a pub const!
  for (unsigned int i = 0; i < 256; i++) {
    tag[i] = bf256_zero();  // for constant values the tag is zero
  }   
}
static void constant_to_vole_256_verifier(bf256_t* key, const uint8_t* val, bf256_t delta) {
  for (unsigned int i = 0; i < 256; i++) {
    key[i] = bf256_mul(bf256_from_bit(get_bit(val[i/8], i%8)), delta);  // multiply delta with each bit
  }  
}

// DEG 2 TO 3
static void aes_128_deg2to3_prover(bf128_t* deg1, bf128_t* deg2, bf128_t tag, bf128_t val) {
  deg1[0] = tag;
  deg2[0] = val;
}
static void aes_128_deg2to3_verifier(bf128_t* deg1, bf128_t key, bf128_t delta) {
  deg1[0] = bf128_mul(key, delta);
}

static void aes_192_deg2to3_prover(bf192_t* deg1, bf192_t* deg2, bf192_t tag, bf192_t val) {
  deg1[0] = tag;
  deg2[0] = val;
}
static void aes_192_deg2to3_verifier(bf192_t* deg1, bf192_t key, bf192_t delta) {
  deg1[0] = bf192_mul(key, delta);
}

static void aes_256_deg2to3_prover(bf256_t* deg1, bf256_t* deg2, bf256_t tag, bf256_t val) {
  deg1[0] = tag;
  deg2[0] = val;
}
static void aes_256_deg2to3_verifier(bf256_t* deg1, bf256_t key, bf256_t delta) {
  deg1[0] = bf256_mul(key, delta);
}
// EncSctrnts internal functions end!!

// COLOUM TO ROW MAJOR
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
static bf192_t* column_to_row_major_and_shrink_V_192(const uint8_t** v, unsigned int ell) {
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
static bf256_t* column_to_row_major_and_shrink_V_256(const uint8_t** v, unsigned int ell) {
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

// KEY EXP BKWD
// DONE: Looks good
static void aes_128_keyexp_backward_prover(uint8_t* y, bf128_t* y_tag, const uint8_t* x, const bf128_t* x_tag, uint8_t* key, bf128_t* key_tag,
                                const faest_paramset_t* params) {

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
    x_tilde = x[j] ^ key[iwd/8 + j%4];  // for the witness
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_tilde_tag[bit_i] = bf128_add(x_tag[j*8 + bit_i], key_tag[iwd + (j%4)*8 + bit_i]); // for the tags of each witness bit
    }

    // ::8-10
    if (rmvRcon == true && j % 4 == 0) {
      // adding round constant to the witness
      x_tilde = x_tilde ^ Rcon[ircon];
      // adding round constant to the tags
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // for prover, no multiplication with delta
        bf128_t rcon_tag;
        constant_to_vole_128_prover(&rcon_tag, (Rcon[ircon] >> bit_i) & 1); // TODO: in spec there should be call to ConstantToVOLE() call
        x_tilde_tag[bit_i] = bf128_add(x_tilde_tag[bit_i], rcon_tag);

      }
      ++ircon;
    }
    // ::11
    aes_128_inverse_affine_prover(y, y_tag, x_tilde, x_tilde_tag, params);

    // ::12-16 lines only relavant for aes-128
    if (j%4 == 4) {
      iwd += 128;
    }
  }
}
// DONE: Looks good
static void aes_128_keyexp_backward_verifier(bf128_t* y_key, const bf128_t* x_key, bf128_t* key_key, const faest_paramset_t* params, bf128_t delta) {

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
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_tilde_tag[bit_i] = bf128_add(x_key[j*8 + bit_i], key_key[iwd + (j%4)*8 + bit_i]); // for the tags of each witness bit
    }

    // ::8-10
    if (rmvRcon == true && j % 4 == 0) {
      // adding round constant to the tags
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // for prover, no multiplication with delta
        bf128_t rcon_key;
        constant_to_vole_128_verifier(&rcon_key, (Rcon[ircon] >> bit_i) & 1, delta); // TODO: in spec there should be call to ConstantToVOLE() call
        x_tilde_tag[bit_i] = bf128_add(x_tilde_tag[bit_i], rcon_key);

      }
      ++ircon;
    }
    // ::11
    aes_128_inverse_affine_verifier(y_key, x_key, delta, params);

    // ::12-16 lines only relavant for aes-128
    if (j%4 == 4) {
      iwd += 128;
    }
  }
}

// KEY EXP FWD
// DONE: Looks good
static void aes_128_keyexp_forward_prover(uint8_t* y, bf128_t* y_tag, const uint8_t* w, const bf128_t* w_tag, const faest_paramset_t* params) {

  unsigned int lambda = params->faest_param.lambda;
  unsigned int Nk = lambda/8;
  unsigned int R = params->faest_param.R;

  // ::1-2
  for (unsigned int i = 0; i < Nk; i++) {
    y[i] = w[i];  // unit8 contains 8 bits
  }
  for (unsigned int i = 0; i < lambda; i++) {
    y_tag[i] = w_tag[i]; 
  }
  // ::3
  unsigned int i_wd = Nk;
  // ::4-10
  for (unsigned int j = Nk; j < 4*(R + 1); j++) {
    // ::5
    if (j % Nk == 0) {
      // ::6
      for (unsigned int word_idx = i_wd; word_idx < 4; word_idx++) {
        y[4*j + word_idx] = w[i_wd + word_idx];   // storing byte by byte
      }
      for (unsigned int word_bit_idx = i_wd*8; word_bit_idx < 4*8; word_bit_idx++) {
        y_tag[4*8*j + word_bit_idx] = w_tag[i_wd*8 + word_bit_idx]; // storing bit tags, "bit by bit"
      }
      // ::7
      i_wd += 4;    // 32 bits -> 4 words
    // ::8
    } else {
      // ::9-10
      for (unsigned int word_idx = 0; word_idx < 4; word_idx++) {
        y[4*j + word_idx] = y[4*(j - Nk) + word_idx] ^ y[4*(j - 1) + word_idx]; // adding bitwise
      }
      for (unsigned int word_bit_idx = 0; word_bit_idx < 32; word_bit_idx++) {
        y_tag[4*8*j + word_bit_idx] = bf128_add(y_tag[4*8*(j - Nk) + word_bit_idx], y_tag[4*8*(j -1) + word_bit_idx]);  // adding the tags in F2lambda "bitwwise"
      }
    }
  }
}
// DONE: Looks good
static void aes_128_keyexp_forward_verifier(bf128_t* y_tag, const bf128_t* w_tag, const faest_paramset_t* params, bf128_t delta) {

  unsigned int lambda = params->faest_param.lambda;
  unsigned int Nk = lambda/8;
  unsigned int R = params->faest_param.R;

  // ::1-2
  for (unsigned int i = 0; i < lambda; i++) {
    y_tag[i] = w_tag[i]; 
  }
  // ::3
  unsigned int i_wd = Nk;
  // ::4-10
  for (unsigned int j = Nk; j < 4*(R + 1); j++) {
    // ::5
    if (j % Nk == 0) {
      // ::6
      for (unsigned int word_bit_idx = i_wd*8; word_bit_idx < 4*8; word_bit_idx++) {
        y_tag[4*8*j + word_bit_idx] = w_tag[i_wd*8 + word_bit_idx]; // storing bit tags, "bit by bit"
      }
      // ::7
      i_wd += 4;    // 32 bits -> 4 words
    // ::8
    } else {
      // ::9-10
      for (unsigned int word_bit_idx = 0; word_bit_idx < 32; word_bit_idx++) {
        y_tag[4*8*j + word_bit_idx] = bf128_add(y_tag[4*8*(j - Nk) + word_bit_idx], y_tag[4*8*(j -1) + word_bit_idx]);  // adding the tags in F2lambda "bitwwise"
      }
    }
  }
}

// KEY EXP CSTRNTS
static void aes_128_expkey_constraints_prover(bf128_t* z_deg0, bf128_t* z_deg1, uint8_t* k, bf128_t* k_tag, const uint8_t* w, const bf128_t* w_tag, 
                                    const faest_paramset_t* params) {

  unsigned int Ske = params->faest_param.Ske;
  unsigned int R = params->faest_param.R;
  unsigned int lambda = params->faest_param.lambda;
  unsigned int Nk = lambda/32;

  // ::1
  aes_128_keyexp_forward_prover(k, k_tag, w, w_tag, params);
  // ::2
  uint8_t w_flat[Ske];
  bf128_t w_flat_tag[8*Ske];
  aes_128_keyexp_backward_prover(w_flat, w_flat_tag, w, w_tag, k, k_tag, params);

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
static void aes_128_expkey_constraints_verifier(bf128_t* z0, bf128_t* k_key, const bf128_t* w_key, bf128_t delta,
                                    const faest_paramset_t* params) {


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

// ENC CSTRNTS
static aes_128_enc_constraints_prover(bf128_t* z0, bf128_t* z1, bf128_t* z2, const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* w, 
                                      const bf128_t* w_tag, const faest_paramset_t* params) {

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

    }

}
static aes_128_enc_constraints_verifier(bf128_t* z0, bf128_t* z1, bf128_t* z2, const uint8_t* owf_in, const uint8_t* owf_out, 
                                        bf128_t* w_key, bf128_t delta, const faest_paramset_t* params) {

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

    }

}

// DONE: Looks good
// OWF CONSTRAINTS
static void aes_128_constraints_prover(bf128_t* z_deg0, bf128_t* z_deg1, bf128_t* z_deg2, const uint8_t* w, const bf128_t* w_tag, const uint8_t* owf_in, 
                                        const uint8_t* owf_out, const faest_paramset_t* params) {

  unsigned int lambda = params->faest_param.lambda;
  unsigned int R = params->faest_param.R;
  unsigned int Ske = params->faest_param.Ske;
  unsigned int Lke = lambda + 8*Ske;
  unsigned int Lenc = params->faest_param.Lenc;
  unsigned int Senc = params->faest_param.Senc;
  // ::1-3 owf_in, owf_out, z and z_tag

  // ::4-5
  aes_128_deg2to3_prover(z_deg1, z_deg2, bf128_from_bit((w[0]&1) & ((w[0]>>1)&1)), bf128_mul(w_tag[0], w_tag[1]));

  // jump to ::13 for AES
  bf128_t owf_in_tag[lambda];
  constant_to_vole_128_prover(owf_in_tag, owf_in);
  // ::14
  bf128_t owf_out_tag[lambda];
  constant_to_vole_128_prover(owf_out_tag, owf_out);
  // ::15 skiped as B = 1
  // ::16
  // TODO: check if the sizes are correct!!
  bf128_t z_tilde_deg0_tag[FAEST_128F_Ske / 4 + 2*4];
  bf128_t z_tilde_deg1_val[FAEST_128F_Ske / 4 + 2*4];
  uint8_t k[(R+1)*lambda/8];
  bf128_t k_tag[(R+1)*lambda];
  aes_128_expkey_constraints_prover(z_tilde_deg0_tag, z_tilde_deg1_val, k, k_tag, w, w_tag, params);

  // ::17
  for (unsigned int i = 0; i < FAEST_128F_Ske/4+2*4; i++) {
    aes_128_deg2to3_prover(z_deg1 + 1, z_deg2 + 1, z_tilde_deg0_tag[i], z_tilde_deg1_val[i]);
  }

  // ::18 b = 0
  // ::19
  uint8_t w_tilde[Lenc/8];
  bf128_t w_tilde_tag[Lenc];
  for (unsigned int i = 0; i < Lenc/8; i++) {
    w_tilde[i] = w[Lke/8 + i];  // copying 8 bits at a time
  }
  for (unsigned int i = 0; i < Lenc; i++) {
    w_tilde_tag[i] = w_tag[Lke + i];  // copying 1 bit's tag at a time
  }
  // ::20 not needed for AES-128
  // ::21
  bf128_t z_tilde_deg0[Senc];
  bf128_t z_tilde_deg1[Senc];
  bf128_t z_tilde_deg2[Senc];
  aes_128_enc_constraints_prover(z_tilde_deg0, z_tilde_deg1, z_tilde_deg2, owf_in, owf_out, w, w_tag, params);

  // :22
  for (unsigned int i = 0; i < Senc; i++) {
    z_deg0[1+(FAEST_128F_Ske/4+2*4) + i] = z_tilde_deg0[i];
    z_deg1[1+(FAEST_128F_Ske/4+2*4) + i] = z_tilde_deg1[i];
    z_deg2[1+(FAEST_128F_Ske/4+2*4) + i] = z_tilde_deg2[i];
  }

}
// DONE: Looks good
static void aes_128_constraints_verifier(bf128_t* z_deg0, bf128_t* z_deg1, bf128_t* z_deg2, const bf128_t* w_key, const uint8_t* owf_in, 
                                        const uint8_t* owf_out, bf128_t delta, const faest_paramset_t* params) {

  unsigned int lambda = params->faest_param.lambda;
  unsigned int R = params->faest_param.R;
  unsigned int Ske = params->faest_param.Ske;
  unsigned int Lke = lambda + 8*Ske;
  unsigned int Lenc = params->faest_param.Lenc;
  unsigned int Senc = params->faest_param.Senc;
  // ::1-3 owf_in, owf_out, z and z_tag

  // ::4-5
  aes_128_deg2to3_verifier(z_deg1, bf128_mul(w_key[0], w_key[1]), delta);
  
  // jump to ::13 for AES
  bf128_t owf_in_key[128];
  constant_to_vole_128_verifier(owf_in_key, owf_in, delta);

  // ::14
  bf128_t owf_out_key[128];
  constant_to_vole_128_verifier(owf_out_key, owf_out, delta);
  // ::15 skiped as B = 1
  // ::16
  bf128_t z_tilde_deg0_key[FAEST_128F_Ske / 4 + 2*4];
  bf128_t z_tilde_deg1_val[FAEST_128F_Ske / 4 + 2*4];
  bf128_t k_key[(R+1)*lambda];
  aes_128_expkey_constraints_verifier(z_tilde_deg0_key, k_key, w_key, delta, params);
  
  // ::17
  for (unsigned int i = 0; i < FAEST_128F_Ske / 4 + 2*4; i++) {
    aes_128_deg2to3_verifier(z_deg1 + 1+i, z_tilde_deg0_key[i], delta);
  }
  // ::18 b = 0
  // ::19
  bf128_t w_tilde_key[Lenc];
  for (unsigned int i = 0; i < Lenc; i++) {
    w_tilde_key[i] = w_key[Lke + i];  // copying 1 bit's key tag at a time
  }
  // ::20 not needed for aes128
  // ::21
  bf128_t z_tilde_deg0[Senc];
  bf128_t z_tilde_deg1[Senc];
  bf128_t z_tilde_deg2[Senc];
  aes_128_enc_constraints_verifier(z_tilde_deg0, z_tilde_deg1, z_tilde_deg2, owf_in, owf_out, w_key, delta, params);

  for (unsigned int i = 0; i < Senc; i++) {
    z_deg0[1+(FAEST_128F_Ske/4+2*4) + i] = z_tilde_deg0[i];
    z_deg1[1+(FAEST_128F_Ske/4+2*4) + i] = z_tilde_deg1[i];
    z_deg2[1+(FAEST_128F_Ske/4+2*4) + i] = z_tilde_deg2[i];
  }

}

// DONE: Looks good
// OWF PROVE VERIFY
static void aes_128_prover(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde, const uint8_t* w, const uint8_t* u, 
                          const uint8_t** V, const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* chall_2, const faest_paramset_t* params) {

  unsigned int lambda = params->faest_param.lambda;
  unsigned int ske = params->faest_param.Ske;
  unsigned int senc = params->faest_param.Senc;
  // TODO: CHANGE THIS FOR OTHER SETTING WHEN COPY PASTING!!!!!
  unsigned int beta = 1;
  unsigned int c = 2*ske + (3/2)*senc + 1;  // TODO: how is c affected if we have 4 bits in 1 unit8_t

  // ::1-5
  // V becomes the w_tag
  bf128_t* w_tag = column_to_row_major_and_shrink_V_128(V, FAEST_128F_ELL); // This is the tag for w

  // ::6-7 embed VOLE masks
  bf128_t bf_u_star_0 = bf128_load(u);
  bf128_t bf_u_star_1 = bf128_load(u + lambda/8);
  // ::8-9
  bf128_t bf_v_star_0 = bf128_sum_poly(w_tag);
  bf128_t bf_v_star_1 = bf128_sum_poly(w_tag + lambda);

  // ::10-12
  bf128_t z0_tag[c]; // this contains the bf tag
  bf128_t z1_val[c]; // this contains the bf val
  bf128_t z2_gamma[c]; // this contains the bf gamma
  aes_128_constraints_prover(z0_tag, z1_val, z2_gamma, w, w_tag, owf_in, owf_out, params);
 
  // Step: 13-18
  zk_hash_128_ctx a0_ctx;
  zk_hash_128_ctx a1_ctx;
  zk_hash_128_ctx a2_ctx;
  zk_hash_128_init(&a0_ctx, chall_2);
  zk_hash_128_init(&a1_ctx, chall_2);
  zk_hash_128_init(&a2_ctx, chall_2);

  for (unsigned int i = 0; i < c*3; i++) {
    zk_hash_128_update(&a0_ctx, z0_tag[i]);
    zk_hash_128_update(&a1_ctx, z1_val[i]);
    zk_hash_128_update(&a2_ctx, z2_gamma[i]);
  }

  zk_hash_128_finalize(a0_tilde, &a0_ctx, bf_u_star_0);
  zk_hash_128_finalize(a1_tilde, &a1_ctx, bf128_add(bf_v_star_0, bf_u_star_1));
  zk_hash_128_finalize(a2_tilde, &a2_ctx, bf_v_star_1);

}

// DONE: Looks good
static uint8_t* aes_128_verifier(const uint8_t* d, const uint8_t** Q, const uint8_t* owf_in, const uint8_t* owf_out,
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
  unsigned int ell = params->faest_param.L; // TODO: I hope L is ELL, that is l_ke + l_enc from Fig. 1.5

  // ::1
  bf128_t bf_delta = bf128_load(chall_3);
  bf128_t bf_delta_sq = bf128_mul(bf_delta, bf_delta);

  // ::2-6
  bf128_t* bf_Q = column_to_row_major_and_shrink_V_128(Q, ell);

  // ::7-9
  bf128_t bf_q_star_0 = bf128_sum_poly(bf_Q + ell);
  bf128_t bf_q_star_1 = bf128_sum_poly(bf_Q + ell + lambda);

  // ::10
  bf128_t bf_q_star = bf128_add(bf_q_star_0, bf128_mul(bf_delta, bf_q_star_1));

  // ::11-12
  bf128_t bf_z0_tag[c];
  bf128_t bf_z1_val[c];
  bf128_t bf_z2_gamma[c];
  bf128_t w_key[ell];
  for (unsigned int i = 0; i < ell; i++) {
    w_key[i] = bf128_add(
                          bf_Q[i], 
                          bf128_mul(bf128_from_bit(get_bit(d[i/8], i%8)),  // TODO: of course, here we have the annoying 4 bit witness problem urghhhh!!!!
                                    bf_delta));
  }
  aes_128_enc_constraints_verifier(bf_z0_tag, bf_z1_val, bf_z2_gamma, owf_in, owf_out, w_key, bf_delta, params);

  // ::13-14
  zk_hash_128_ctx b_ctx;
  zk_hash_128_init(&b_ctx, chall_2);
  for (unsigned int i = 0; i < c; i++) {
    zk_hash_128_update(&b_ctx, bf_z2_gamma[i]);
  }
  uint8_t q_tilde[lambda/8*c];
  zk_hash_128_finalize(q_tilde, &b_ctx, bf_q_star);

  // ::16
  bf128_t tmp1 = bf128_mul(bf128_load(a1_tilde), bf_delta);
  bf128_t tmp2 = bf128_mul(bf128_load(a2_tilde), bf_delta_sq);
  bf128_t tmp3 = bf128_add(tmp1, tmp2);
  bf128_t ret = bf128_add(bf128_load(q_tilde), tmp3);

  uint8_t* a0_tilde = malloc(lambda/8*c);
  bf128_store(a0_tilde, ret);
  return a0_tilde;

}


// AES dispatchers
void aes_prove(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde, const uint8_t* w, const uint8_t* u, const uint8_t** V, 
                const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* chall_2, const faest_paramset_t* params) {
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
      aes_128_prover(a0_tilde, a1_tilde, a2_tilde, w, u, V, owf_in, owf_out, chall_2, params);
    } else {
      // em_prove_128(w, u, V, owf_in, owf_out, chall_2, a0_tilde, a12_tilde);
    }
  }
}

uint8_t* aes_verify(const uint8_t* d, const uint8_t** Q, const uint8_t* chall_2, const uint8_t* chall_3,
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
      return aes_128_verifier(d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
    } else {
      // return em_verify_128(d, Q, chall_2, chall_3, a_tilde, in, out, params);
    }
  }
}

