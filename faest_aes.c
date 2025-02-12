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

#include <stdio.h>

static_assert(FAEST_128_LAMBDA == FAEST_128S_LAMBDA, "Invalid parameters");
static_assert(FAEST_128F_ELL == FAEST_128S_ELL, "Invalid parameters");
static_assert(FAEST_128F_LAMBDA == FAEST_128S_LAMBDA, "Invalid parameters");
static_assert(FAEST_128F_Lke == FAEST_128S_Lke, "Invalid parameters");
static_assert(FAEST_128F_Nst == FAEST_128S_Nst, "Invalid parameters");
static_assert(FAEST_128F_R == FAEST_128S_R, "Invalid parameters");
static_assert(FAEST_128F_Senc == FAEST_128S_Senc, "Invalid parameters");
static_assert(FAEST_128F_Ske == FAEST_128S_Ske, "Invalid parameters");
static_assert(FAEST_128F_C == FAEST_128S_C, "Invalid parameters");

static_assert(FAEST_192_LAMBDA == FAEST_192S_LAMBDA, "Invalid parameters");
static_assert(FAEST_192F_ELL == FAEST_192S_ELL, "Invalid parameters");
static_assert(FAEST_192F_LAMBDA == FAEST_192S_LAMBDA, "Invalid parameters");
static_assert(FAEST_192F_Lke == FAEST_192S_Lke, "Invalid parameters");
static_assert(FAEST_192F_Nst == FAEST_192S_Nst, "Invalid parameters");
static_assert(FAEST_192F_R == FAEST_192S_R, "Invalid parameters");
static_assert(FAEST_192F_Senc == FAEST_192S_Senc, "Invalid parameters");
static_assert(FAEST_192F_Ske == FAEST_192S_Ske, "Invalid parameters");
static_assert(FAEST_192F_C == FAEST_192S_C, "Invalid parameters");

static_assert(FAEST_256_LAMBDA == FAEST_256S_LAMBDA, "Invalid parameters");
static_assert(FAEST_256F_ELL == FAEST_256S_ELL, "Invalid parameters");
static_assert(FAEST_256F_LAMBDA == FAEST_256S_LAMBDA, "Invalid parameters");
static_assert(FAEST_256F_Lke == FAEST_256S_Lke, "Invalid parameters");
static_assert(FAEST_256F_Nst == FAEST_256S_Nst, "Invalid parameters");
static_assert(FAEST_256F_R == FAEST_256S_R, "Invalid parameters");
static_assert(FAEST_256F_Senc == FAEST_256S_Senc, "Invalid parameters");
static_assert(FAEST_256F_Ske == FAEST_256S_Ske, "Invalid parameters");
static_assert(FAEST_256F_C == FAEST_256S_C, "Invalid parameters");

static_assert(FAEST_128_LAMBDA == FAEST_EM_128S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_128F_LAMBDA == FAEST_EM_128S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_128F_Lenc == FAEST_EM_128S_Lenc, "Invalid parameters");
static_assert(FAEST_EM_128F_Nst == FAEST_EM_128S_Nst, "Invalid parameters");
static_assert(FAEST_EM_128F_R == FAEST_EM_128S_R, "Invalid parameters");
static_assert(FAEST_EM_128F_Senc == FAEST_EM_128S_Senc, "Invalid parameters");
static_assert(FAEST_EM_128F_C == FAEST_EM_128S_C, "Invalid parameters");
// for scan-build
static_assert(FAEST_EM_128F_LAMBDA * (FAEST_EM_128F_R + 1) / 8 ==
                  sizeof(aes_word_t) * FAEST_EM_128F_Nst * (FAEST_EM_128F_R + 1),
              "Invalid parameters");

static_assert(FAEST_192_LAMBDA == FAEST_EM_192S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_192F_LAMBDA == FAEST_EM_192S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_192F_Lenc == FAEST_EM_192S_Lenc, "Invalid parameters");
static_assert(FAEST_EM_192F_Nst == FAEST_EM_192S_Nst, "Invalid parameters");
static_assert(FAEST_EM_192F_R == FAEST_EM_192S_R, "Invalid parameters");
static_assert(FAEST_EM_192F_Senc == FAEST_EM_192S_Senc, "Invalid parameters");
static_assert(FAEST_EM_192F_C == FAEST_EM_192S_C, "Invalid parameters");
// for scan-build
static_assert(FAEST_EM_192F_LAMBDA * (FAEST_EM_192F_R + 1) / 8 ==
                  sizeof(aes_word_t) * FAEST_EM_192F_Nst * (FAEST_EM_192F_R + 1),
              "Invalid parameters");

static_assert(FAEST_256_LAMBDA == FAEST_EM_256S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_256F_LAMBDA == FAEST_EM_256S_LAMBDA, "Invalid parameters");
static_assert(FAEST_EM_256F_Lenc == FAEST_EM_256S_Lenc, "Invalid parameters");
static_assert(FAEST_EM_256F_Nst == FAEST_EM_256S_Nst, "Invalid parameters");
static_assert(FAEST_EM_256F_R == FAEST_EM_256S_R, "Invalid parameters");
static_assert(FAEST_EM_256F_Senc == FAEST_EM_256S_Senc, "Invalid parameters");
static_assert(FAEST_EM_256F_C == FAEST_EM_256S_C, "Invalid parameters");
// for scan-build
static_assert(FAEST_EM_256F_LAMBDA * (FAEST_EM_256F_R + 1) / 8 ==
                  sizeof(aes_word_t) * FAEST_EM_256F_Nst * (FAEST_EM_256F_R + 1),
              "Invalid parameters");

/* pad sizes to multiples of some value; assumes that a is a power of 2 */
#define PAD_TO(s, a) (((s) + (a) - 1) & ~((a) - 1))

#define BF128_ALLOC(s) faest_aligned_alloc(BF128_ALIGN, PAD_TO((s) * sizeof(bf128_t), BF128_ALIGN))
#define BF192_ALLOC(s) faest_aligned_alloc(BF192_ALIGN, PAD_TO((s) * sizeof(bf192_t), BF192_ALIGN))
#define BF256_ALLOC(s) faest_aligned_alloc(BF256_ALIGN, PAD_TO((s) * sizeof(bf256_t), BF256_ALIGN))

static const bf8_t Rcon[30] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a,
    0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91,
};

// ADD ROUND KEY
/*
Called in EncCstrnts, takes in owf_in (bits) and their tags (0 for prover, owf_in (bit) * delta for
verifier)
*/
static void aes_128_add_round_key_prover(uint8_t* out, bf128_t* out_tag, const uint8_t* in,
                                         const bf128_t* in_tag, const uint8_t* k,
                                         const bf128_t* k_tag, const faest_paramset_t* params) {
  const unsigned int Nst     = params->Nst;
  const unsigned int Nstbits = Nst * 32;
  for (unsigned int i = 0; i < Nstbits; i++) {
    out[i]     = in[i] ^ k[i];
    out_tag[i] = bf128_add(in_tag[i], k_tag[i]);
  }
}

static void aes_192_add_round_key_prover(uint8_t* out, bf192_t* out_tag, const uint8_t* in,
                                         const bf192_t* in_tag, const uint8_t* k,
                                         const bf192_t* k_tag, const faest_paramset_t* params) {
  const unsigned int Nst     = params->Nst;
  const unsigned int Nstbits = Nst * 32;
  for (unsigned int i = 0; i < Nstbits; i++) {
    out[i]     = in[i] ^ k[i];
    out_tag[i] = bf192_add(in_tag[i], k_tag[i]);
  }
}

static void aes_256_add_round_key_prover(uint8_t* out, bf256_t* out_tag, const uint8_t* in,
                                         const bf256_t* in_tag, const uint8_t* k,
                                         const bf256_t* k_tag, const faest_paramset_t* params) {
  const unsigned int Nst     = params->Nst;
  const unsigned int Nstbits = Nst * 32;
  for (unsigned int i = 0; i < Nstbits; i++) {
    out[i]     = in[i] ^ k[i];
    out_tag[i] = bf256_add(in_tag[i], k_tag[i]);
  }
}

static void aes_128_add_round_key_verifier(bf128_t* out_key, const bf128_t* in_key,
                                           const bf128_t* k_key, const faest_paramset_t* params) {
  const unsigned int Nst     = params->Nst;
  const unsigned int Nstbits = Nst * 32;

  for (unsigned int i = 0; i < Nstbits; i++) {
    out_key[i] = bf128_add(in_key[i], k_key[i]);
  }
}

static void aes_192_add_round_key_verifier(bf192_t* out_key, const bf192_t* in_key,
                                           const bf192_t* k_key, const faest_paramset_t* params) {
  const unsigned int Nst     = params->Nst;
  const unsigned int Nstbits = Nst * 32;

  for (unsigned int i = 0; i < Nstbits; i++) {
    out_key[i] = bf192_add(in_key[i], k_key[i]);
  }
}

static void aes_256_add_round_key_verifier(bf256_t* out_key, const bf256_t* in_key,
                                           const bf256_t* k_key, const faest_paramset_t* params) {
  const unsigned int Nst     = params->Nst;
  const unsigned int Nstbits = Nst * 32;

  for (unsigned int i = 0; i < Nstbits; i++) {
    out_key[i] = bf256_add(in_key[i], k_key[i]);
  }
}

// F256/F2.CONJUGATES
static void aes_128_f256_f2_conjugates_1(bf128_t* y, const uint8_t* state,
                                         const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i != Nst_bytes; ++i) {
    uint8_t x0[8];
    // array of bits containing byte i
    for (unsigned int j = 0; j != 8; ++j) {
      x0[j] = state[i * 8 + j];
    }

    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf128_byte_combine_bits_array(x0);
      bits_sq_array(x0);
    }
    y[i * 8 + 7] = bf128_byte_combine_bits_array(x0);
  }
}

static void aes_192_f256_f2_conjugates_1(bf192_t* y, const uint8_t* state,
                                         const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i != Nst_bytes; ++i) {
    uint8_t x0[8];
    // array of bits containing byte i
    for (unsigned int j = 0; j != 8; ++j) {
      x0[j] = state[i * 8 + j];
    }

    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf192_byte_combine_bits_array(x0);
      bits_sq_array(x0);
    }
    y[i * 8 + 7] = bf192_byte_combine_bits_array(x0);
  }
}

static void aes_256_f256_f2_conjugates_1(bf256_t* y, const uint8_t* state,
                                         const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i != Nst_bytes; ++i) {
    uint8_t x0[8];
    // array of bits containing byte i
    for (unsigned int j = 0; j != 8; ++j) {
      x0[j] = state[i * 8 + j];
    }

    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf256_byte_combine_bits_array(x0);
      bits_sq_array(x0);
    }
    y[i * 8 + 7] = bf256_byte_combine_bits_array(x0);
  }
}

static void aes_128_f256_f2_conjugates_128(bf128_t* y, const bf128_t* state,
                                           const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i != Nst_bytes; ++i) {
    bf128_t x[8];
    memcpy(x, state + i * 8, sizeof(x));
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf128_byte_combine(x);
      bf128_t tmp[8];
      memcpy(tmp, x, sizeof(x));
      bf128_sq_bit(x, tmp);
    }
    y[i * 8 + 7] = bf128_byte_combine(x);
  }
}

static void aes_192_f256_f2_conjugates_192(bf192_t* y, const bf192_t* state,
                                           const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i != Nst_bytes; ++i) {
    bf192_t x[8];
    memcpy(x, state + (i * 8), sizeof(x));
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf192_byte_combine(x);
      bf192_t tmp[8];
      memcpy(tmp, x, sizeof(x));
      bf192_sq_bit(x, tmp);
    }
    y[i * 8 + 7] = bf192_byte_combine(x);
  }
}

static void aes_256_f256_f2_conjugates_256(bf256_t* y, const bf256_t* state,
                                           const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i != Nst_bytes; ++i) {
    bf256_t x[8];
    memcpy(x, state + (i * 8), sizeof(x));
    for (unsigned int j = 0; j != 7; ++j) {
      y[i * 8 + j] = bf256_byte_combine(x);
      bf256_t tmp[8];
      memcpy(tmp, x, sizeof(x));
      bf256_sq_bit(x, tmp);
    }
    y[i * 8 + 7] = bf256_byte_combine(x);
  }
}

// INV NORM TO CONJUGATES

static const bf128_t bf128_beta_squares[5] = {
    BF128C(UINT64_C(0xaca8c2a7190f676b), UINT64_C(0xdb4a932e2cae3d8b)),
    BF128C(UINT64_C(0x34d2f7fba603e341), UINT64_C(0x500317bd159d73bb)),
    BF128C(UINT64_C(0xcbc26e38bdbd6c62), UINT64_C(0xf210539fd8dd2772)),
    BF128C(UINT64_C(0x53b85b6402b1e849), UINT64_C(0x7959d70ce1ee6942)),
    BF128C(UINT64_C(0xaca8c2a7190f676b), UINT64_C(0xdb4a932e2cae3d8b)),
};
static const bf128_t bf128_beta_cubes[4] = {
    BF128C(UINT64_C(0x53b85b6402b1e849), UINT64_C(0x7959d70ce1ee6942)),
    BF128C(UINT64_C(0xaca8c2a7190f676b), UINT64_C(0xdb4a932e2cae3d8b)),
    BF128C(UINT64_C(0x34d2f7fba603e341), UINT64_C(0x500317bd159d73bb)),
    BF128C(UINT64_C(0xcbc26e38bdbd6c62), UINT64_C(0xf210539fd8dd2772)),
};

static const bf192_t bf192_beta_squares[5] = {
    BF192C(UINT64_C(0x4a2fe80a53fe16e7), UINT64_C(0x89bfb5055f8ff266),
           UINT64_C(0x6c4aecc3fdd0f812)),
    BF192C(UINT64_C(0xf061255c52e359d4), UINT64_C(0x125a337e5e808579),
           UINT64_C(0x54ed13a2d09da6a4)),
    BF192C(UINT64_C(0x941080db88d09584), UINT64_C(0x62217a15d25ec520),
           UINT64_C(0xee0fc537c5a9cb74)),
    BF192C(UINT64_C(0x2e5e4d8d89cddab6), UINT64_C(0xf9c4fc6ed351b23f),
           UINT64_C(0xd6a83a56e8e495c2)),
    BF192C(UINT64_C(0x4a2fe80a53fe16e7), UINT64_C(0x89bfb5055f8ff266),
           UINT64_C(0x6c4aecc3fdd0f812)),
};
static const bf192_t bf192_beta_cubes[4] = {
    BF192C(UINT64_C(0x2e5e4d8d89cddab6), UINT64_C(0xf9c4fc6ed351b23f),
           UINT64_C(0xd6a83a56e8e495c2)),
    BF192C(UINT64_C(0x4a2fe80a53fe16e7), UINT64_C(0x89bfb5055f8ff266),
           UINT64_C(0x6c4aecc3fdd0f812)),
    BF192C(UINT64_C(0xf061255c52e359d4), UINT64_C(0x125a337e5e808579),
           UINT64_C(0x54ed13a2d09da6a4)),
    BF192C(UINT64_C(0x941080db88d09584), UINT64_C(0x62217a15d25ec520),
           UINT64_C(0xee0fc537c5a9cb74)),
};

static const bf256_t bf256_beta_squares[5] = {
    BF256C(UINT64_C(0x4074aaa06e5faade), UINT64_C(0xd9ca3e577cc14c4a), UINT64_C(0xda065d89ba07bf41),
           UINT64_C(0x96ff78f99e8ca1cb)),
    BF256C(UINT64_C(0xfc8046a8eaa315e8), UINT64_C(0x3cd447937f790879), UINT64_C(0x2adc5a08c13f1cfb),
           UINT64_C(0x67dbaeb2ec4abc7c)),
    BF256C(UINT64_C(0x4074aaa16e5eabeb), UINT64_C(0xd9ca3e577cc14c4b), UINT64_C(0xda065d89ba07bf40),
           UINT64_C(0x96ff78f99e8ca1cb)),
    BF256C(UINT64_C(0xfc8046a9eaa214dc), UINT64_C(0x3cd447937f790878), UINT64_C(0x2adc5a08c13f1cfa),
           UINT64_C(0x67dbaeb2ec4abc7c)),
    BF256C(UINT64_C(0x4074aaa06e5faade), UINT64_C(0xd9ca3e577cc14c4a), UINT64_C(0xda065d89ba07bf41),
           UINT64_C(0x96ff78f99e8ca1cb)),
};
static const bf256_t bf256_beta_cubes[4] = {
    BF256C(UINT64_C(0xfc8046a9eaa214dc), UINT64_C(0x3cd447937f790878), UINT64_C(0x2adc5a08c13f1cfa),
           UINT64_C(0x67dbaeb2ec4abc7c)),
    BF256C(UINT64_C(0x4074aaa06e5faade), UINT64_C(0xd9ca3e577cc14c4a), UINT64_C(0xda065d89ba07bf41),
           UINT64_C(0x96ff78f99e8ca1cb)),
    BF256C(UINT64_C(0xfc8046a8eaa315e8), UINT64_C(0x3cd447937f790879), UINT64_C(0x2adc5a08c13f1cfb),
           UINT64_C(0x67dbaeb2ec4abc7c)),
    BF256C(UINT64_C(0x4074aaa16e5eabeb), UINT64_C(0xd9ca3e577cc14c4b), UINT64_C(0xda065d89ba07bf40),
           UINT64_C(0x96ff78f99e8ca1cb)),
};

static void aes_128_inv_norm_to_conjugates_prover(bf128_t* y_val, bf128_t* y_tag,
                                                  const uint8_t* x_val, const bf128_t* x_tag) {
  // :1-2
  for (unsigned int i = 0; i != 4; ++i) {
    y_val[i] = bf128_add(bf128_add(bf128_mul_bit(bf128_one(), x_val[0]),
                                   bf128_mul_bit(bf128_beta_squares[i], x_val[1])),
                         bf128_add(bf128_mul_bit(bf128_beta_squares[i + 1], x_val[2]),
                                   bf128_mul_bit(bf128_beta_cubes[i], x_val[3])));
    y_tag[i] = bf128_add(
        bf128_add(bf128_mul(bf128_one(), x_tag[0]), bf128_mul(bf128_beta_squares[i], x_tag[1])),
        bf128_add(bf128_mul(bf128_beta_squares[i + 1], x_tag[2]),
                  bf128_mul(bf128_beta_cubes[i], x_tag[3])));
  }
}

static void aes_192_inv_norm_to_conjugates_prover(bf192_t* y_val, bf192_t* y_tag,
                                                  const uint8_t* x_val, const bf192_t* x_tag) {
  // :1-2
  for (unsigned int i = 0; i != 4; ++i) {
    y_val[i] = bf192_add(bf192_add(bf192_mul_bit(bf192_one(), x_val[0]),
                                   bf192_mul_bit(bf192_beta_squares[i], x_val[1])),
                         bf192_add(bf192_mul_bit(bf192_beta_squares[i + 1], x_val[2]),
                                   bf192_mul_bit(bf192_beta_cubes[i], x_val[3])));
    y_tag[i] = bf192_add(
        bf192_add(bf192_mul(bf192_one(), x_tag[0]), bf192_mul(bf192_beta_squares[i], x_tag[1])),
        bf192_add(bf192_mul(bf192_beta_squares[i + 1], x_tag[2]),
                  bf192_mul(bf192_beta_cubes[i], x_tag[3])));
  }
}

static void aes_256_inv_norm_to_conjugates_prover(bf256_t* y_val, bf256_t* y_tag,
                                                  const uint8_t* x_val, const bf256_t* x_tag) {
  // :1-2
  for (unsigned int i = 0; i != 4; ++i) {
    y_val[i] = bf256_add(bf256_add(bf256_mul_bit(bf256_one(), x_val[0]),
                                   bf256_mul_bit(bf256_beta_squares[i], x_val[1])),
                         bf256_add(bf256_mul_bit(bf256_beta_squares[i + 1], x_val[2]),
                                   bf256_mul_bit(bf256_beta_cubes[i], x_val[3])));
    y_tag[i] = bf256_add(
        bf256_add(bf256_mul(bf256_one(), x_tag[0]), bf256_mul(bf256_beta_squares[i], x_tag[1])),
        bf256_add(bf256_mul(bf256_beta_squares[i + 1], x_tag[2]),
                  bf256_mul(bf256_beta_cubes[i], x_tag[3])));
  }
}

static void aes_128_inv_norm_to_conjugates_verifier(bf128_t* y_eval, const bf128_t* x_eval) {
  // :1-2
  for (unsigned int i = 0; i != 4; ++i) {
    y_eval[i] = bf128_add(
        bf128_add(bf128_mul(bf128_one(), x_eval[0]), bf128_mul(bf128_beta_squares[i], x_eval[1])),
        bf128_add(bf128_mul(bf128_beta_squares[i + 1], x_eval[2]),
                  bf128_mul(bf128_beta_cubes[i], x_eval[3])));
  }
}

static void aes_192_inv_norm_to_conjugates_verifier(bf192_t* y_eval, const bf192_t* x_eval) {
  // :1-2
  for (unsigned int i = 0; i != 4; ++i) {
    y_eval[i] = bf192_add(
        bf192_add(bf192_mul(bf192_one(), x_eval[0]), bf192_mul(bf192_beta_squares[i], x_eval[1])),
        bf192_add(bf192_mul(bf192_beta_squares[i + 1], x_eval[2]),
                  bf192_mul(bf192_beta_cubes[i], x_eval[3])));
  }
}

static void aes_256_inv_norm_to_conjugates_verifier(bf256_t* y_eval, const bf256_t* x_eval) {
  // :1-2
  for (unsigned int i = 0; i != 4; ++i) {
    y_eval[i] = bf256_add(
        bf256_add(bf256_mul(bf256_one(), x_eval[0]), bf256_mul(bf256_beta_squares[i], x_eval[1])),
        bf256_add(bf256_mul(bf256_beta_squares[i + 1], x_eval[2]),
                  bf256_mul(bf256_beta_cubes[i], x_eval[3])));
  }
}

// // INV NORM CONSTRAINTS
static void aes_128_inv_norm_constraints_prover(zk_hash_128_3_ctx* hasher,
                                                const bf128_t* conjugates,
                                                const bf128_t* conjugates_tag, const bf128_t* y,
                                                const bf128_t* y_tag) {
  zk_hash_128_3_update(
      hasher, bf128_mul(bf128_mul(*y_tag, conjugates_tag[1]), conjugates_tag[4]),
      bf128_add(bf128_add(bf128_mul(bf128_mul(*y, conjugates_tag[1]), conjugates_tag[4]),
                          bf128_mul(bf128_mul(*y_tag, conjugates_tag[1]), conjugates[4])),
                bf128_mul(bf128_mul(*y_tag, conjugates[1]), conjugates_tag[4])),
      bf128_add(bf128_add(bf128_add(bf128_mul(bf128_mul(*y, conjugates[1]), conjugates_tag[4]),
                                    bf128_mul(bf128_mul(*y, conjugates_tag[1]), conjugates[4])),
                          bf128_mul(bf128_mul(*y_tag, conjugates[1]), conjugates[4])),
                conjugates_tag[0]));
}

static void aes_192_inv_norm_constraints_prover(bf192_t* z_deg0, bf192_t* z_deg1, bf192_t* z_deg2,
                                                const bf192_t* conjugates,
                                                const bf192_t* conjugates_tag, const bf192_t* y,
                                                const bf192_t* y_tag) {
  z_deg0[0] = bf192_mul(bf192_mul(*y_tag, conjugates_tag[1]), conjugates_tag[4]);

  z_deg1[0] = bf192_add(bf192_add(bf192_mul(bf192_mul(*y, conjugates_tag[1]), conjugates_tag[4]),
                                  bf192_mul(bf192_mul(*y_tag, conjugates_tag[1]), conjugates[4])),
                        bf192_mul(bf192_mul(*y_tag, conjugates[1]), conjugates_tag[4]));

  z_deg2[0] =
      bf192_add(bf192_add(bf192_add(bf192_mul(bf192_mul(*y, conjugates[1]), conjugates_tag[4]),
                                    bf192_mul(bf192_mul(*y, conjugates_tag[1]), conjugates[4])),
                          bf192_mul(bf192_mul(*y_tag, conjugates[1]), conjugates[4])),
                conjugates_tag[0]);
}

static void aes_256_inv_norm_constraints_prover(bf256_t* z_deg0, bf256_t* z_deg1, bf256_t* z_deg2,
                                                const bf256_t* conjugates,
                                                const bf256_t* conjugates_tag, const bf256_t* y,
                                                const bf256_t* y_tag) {
  z_deg0[0] = bf256_mul(bf256_mul(*y_tag, conjugates_tag[1]), conjugates_tag[4]);

  z_deg1[0] = bf256_add(bf256_add(bf256_mul(bf256_mul(*y, conjugates_tag[1]), conjugates_tag[4]),
                                  bf256_mul(bf256_mul(*y_tag, conjugates_tag[1]), conjugates[4])),
                        bf256_mul(bf256_mul(*y_tag, conjugates[1]), conjugates_tag[4]));

  z_deg2[0] =
      bf256_add(bf256_add(bf256_add(bf256_mul(bf256_mul(*y, conjugates[1]), conjugates_tag[4]),
                                    bf256_mul(bf256_mul(*y, conjugates_tag[1]), conjugates[4])),
                          bf256_mul(bf256_mul(*y_tag, conjugates[1]), conjugates[4])),
                conjugates_tag[0]);
}

static void aes_128_inv_norm_constraints_verifier(zk_hash_128_ctx* hasher,
                                                  const bf128_t* conjugates_eval,
                                                  const bf128_t* y_eval, const bf128_t delta) {
  zk_hash_128_update(
      hasher, bf128_add(bf128_mul(bf128_mul(*y_eval, conjugates_eval[1]), conjugates_eval[4]),
                        bf128_mul(conjugates_eval[0], bf128_mul(delta, delta))));
}

static void aes_192_inv_norm_constraints_verifier(bf192_t* z_eval, const bf192_t* conjugates_eval,
                                                  const bf192_t* y_eval, const bf192_t delta) {
  z_eval[0] = bf192_add(bf192_mul(bf192_mul(*y_eval, conjugates_eval[1]), conjugates_eval[4]),
                        bf192_mul(conjugates_eval[0], bf192_mul(delta, delta)));
}

static void aes_256_inv_norm_constraints_verifier(bf256_t* z_eval, const bf256_t* conjugates_eval,
                                                  const bf256_t* y_eval, const bf256_t delta) {
  z_eval[0] = bf256_add(bf256_mul(bf256_mul(*y_eval, conjugates_eval[1]), conjugates_eval[4]),
                        bf256_mul(conjugates_eval[0], bf256_mul(delta, delta)));
}

// STATE TO BYTES
static void aes_128_state_to_bytes_prover(bf128_t* out, bf128_t* out_tag, const uint8_t* k,
                                          const bf128_t* k_tag, const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    out[i]     = bf128_byte_combine_bits_array(k + i * 8);
    out_tag[i] = bf128_byte_combine(k_tag + i * 8);
  }
}

static void aes_192_state_to_bytes_prover(bf192_t* out, bf192_t* out_tag, const uint8_t* k,
                                          const bf192_t* k_tag, const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    out[i]     = bf192_byte_combine_bits_array(k + i * 8);
    out_tag[i] = bf192_byte_combine(k_tag + i * 8);
  }
}

static void aes_256_state_to_bytes_prover(bf256_t* out, bf256_t* out_tag, const uint8_t* k,
                                          const bf256_t* k_tag, const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    out[i]     = bf256_byte_combine_bits_array(k + i * 8);
    out_tag[i] = bf256_byte_combine(k_tag + i * 8);
  }
}

static void aes_128_state_to_bytes_verifier(bf128_t* out_key, const bf128_t* k_key,
                                            const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    out_key[i] = bf128_byte_combine(k_key + i * 8);
  }
}

static void aes_192_state_to_bytes_verifier(bf192_t* out_key, const bf192_t* k_key,
                                            const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    out_key[i] = bf192_byte_combine(k_key + i * 8);
  }
}

static void aes_256_state_to_bytes_verifier(bf256_t* out_key, const bf256_t* k_key,
                                            const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    out_key[i] = bf256_byte_combine(k_key + i * 8);
  }
}

// SBOX AFFINE
static const uint8_t sbox_x[9]    = {0x05, 0x09, 0xf9, 0x25, 0xf4, 0x01, 0xb5, 0x8f, 0x63};
static const uint8_t sbox_x_sq[9] = {0x11, 0x41, 0x07, 0x7d, 0x56, 0x01, 0xfc, 0xcf, 0xc2};

static void aes_128_sbox_affine_prover(bf128_t* out_deg0, bf128_t* out_deg1, bf128_t* out_deg2,
                                       const bf128_t* in_deg0, const bf128_t* in_deg1,
                                       const bf128_t* in_deg2, bool dosq,
                                       const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  bf128_t C[9];
  uint8_t t;
  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf128_byte_combine_bits(sbox_x_sq[i]);
    }
    t = 1;
  } else {
    t = 0;
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf128_byte_combine_bits(sbox_x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 8; Cidx++) {
      out_deg2[i] = bf128_add(out_deg2[i], bf128_mul(C[Cidx], in_deg2[i * 8 + (Cidx + t) % 8]));
      out_deg1[i] = bf128_add(out_deg1[i], bf128_mul(C[Cidx], in_deg1[i * 8 + (Cidx + t) % 8]));
      out_deg0[i] = bf128_add(out_deg0[i], bf128_mul(C[Cidx], in_deg0[i * 8 + (Cidx + t) % 8]));
    }
    // add the constant C[8] to the highest coefficient
    out_deg2[i] = bf128_add(out_deg2[i], C[8]);
  }
}

static void aes_192_sbox_affine_prover(bf192_t* out_deg0, bf192_t* out_deg1, bf192_t* out_deg2,
                                       const bf192_t* in_deg0, const bf192_t* in_deg1,
                                       const bf192_t* in_deg2, bool dosq,
                                       const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  bf192_t C[9];
  uint8_t t;
  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf192_byte_combine_bits(sbox_x_sq[i]);
    }
    t = 1;
  } else {
    t = 0;
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf192_byte_combine_bits(sbox_x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 8; Cidx++) {
      out_deg2[i] = bf192_add(out_deg2[i], bf192_mul(C[Cidx], in_deg2[i * 8 + (Cidx + t) % 8]));
      out_deg1[i] = bf192_add(out_deg1[i], bf192_mul(C[Cidx], in_deg1[i * 8 + (Cidx + t) % 8]));
      out_deg0[i] = bf192_add(out_deg0[i], bf192_mul(C[Cidx], in_deg0[i * 8 + (Cidx + t) % 8]));
    }
    // add the constant C[8] to the highest coefficient
    out_deg2[i] = bf192_add(out_deg2[i], C[8]);
  }
}

static void aes_256_sbox_affine_prover(bf256_t* out_deg0, bf256_t* out_deg1, bf256_t* out_deg2,
                                       const bf256_t* in_deg0, const bf256_t* in_deg1,
                                       const bf256_t* in_deg2, bool dosq,
                                       const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  bf256_t C[9];
  uint8_t t;
  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf256_byte_combine_bits(sbox_x_sq[i]);
    }
    t = 1;
  } else {
    t = 0;
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf256_byte_combine_bits(sbox_x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 8; Cidx++) {
      out_deg2[i] = bf256_add(out_deg2[i], bf256_mul(C[Cidx], in_deg2[i * 8 + (Cidx + t) % 8]));
      out_deg1[i] = bf256_add(out_deg1[i], bf256_mul(C[Cidx], in_deg1[i * 8 + (Cidx + t) % 8]));
      out_deg0[i] = bf256_add(out_deg0[i], bf256_mul(C[Cidx], in_deg0[i * 8 + (Cidx + t) % 8]));
    }
    // add the constant C[8] to the highest coefficient
    out_deg2[i] = bf256_add(out_deg2[i], C[8]);
  }
}

static void aes_128_sbox_affine_verifier(bf128_t* out_deg1, const bf128_t* in_deg1, bf128_t delta,
                                         bool dosq, const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  bf128_t C[9];
  uint8_t t;
  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf128_byte_combine_bits(sbox_x_sq[i]);
    }
    t = 1;
  } else {
    t = 0;
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf128_byte_combine_bits(sbox_x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 8; Cidx++) {
      out_deg1[i] = bf128_add(out_deg1[i], bf128_mul(C[Cidx], in_deg1[i * 8 + (Cidx + t) % 8]));
    }
    // add the constant C[8] by multiplying with delta^2
    out_deg1[i] = bf128_add(out_deg1[i], bf128_mul(C[8], bf128_mul(delta, delta)));
  }
}

static void aes_192_sbox_affine_verifier(bf192_t* out_deg1, const bf192_t* in_deg1, bf192_t delta,
                                         bool dosq, const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  bf192_t C[9];
  uint8_t t;
  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf192_byte_combine_bits(sbox_x_sq[i]);
    }
    t = 1;
  } else {
    t = 0;
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf192_byte_combine_bits(sbox_x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 8; Cidx++) {
      out_deg1[i] = bf192_add(out_deg1[i], bf192_mul(C[Cidx], in_deg1[i * 8 + (Cidx + t) % 8]));
    }
    // add the constant C[8] by multiplying with delta^2
    out_deg1[i] = bf192_add(out_deg1[i], bf192_mul(C[8], bf192_mul(delta, delta)));
  }
}

static void aes_256_sbox_affine_verifier(bf256_t* out_deg1, const bf256_t* in_deg1, bf256_t delta,
                                         bool dosq, const faest_paramset_t* params) {
  const unsigned int Nst_bytes = params->Nst * 4;

  bf256_t C[9];
  uint8_t t;
  // ::5-6
  if (dosq) {
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf256_byte_combine_bits(sbox_x_sq[i]);
    }
    t = 1;
  } else {
    t = 0;
    for (unsigned i = 0; i < 9; i++) {
      C[i] = bf256_byte_combine_bits(sbox_x[i]);
    }
  }

  for (unsigned int i = 0; i < Nst_bytes; i++) {
    for (unsigned int Cidx = 0; Cidx < 8; Cidx++) {
      out_deg1[i] = bf256_add(out_deg1[i], bf256_mul(C[Cidx], in_deg1[i * 8 + (Cidx + t) % 8]));
    }
    // add the constant C[8] by multiplying with delta^2
    out_deg1[i] = bf256_add(out_deg1[i], bf256_mul(C[8], bf256_mul(delta, delta)));
  }
}

// SHIFT ROWS
static void aes_128_shiftrows_prover(bf128_t* out_deg0, bf128_t* out_deg1, bf128_t* out_deg2,
                                     const bf128_t* in_deg0, const bf128_t* in_deg1,
                                     const bf128_t* in_deg2, const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if ((Nst != 8) || (r <= 1)) {
        out_deg2[4 * c + r] = in_deg2[4 * ((c + r) % Nst) + r];
        out_deg1[4 * c + r] = in_deg1[4 * ((c + r) % Nst) + r];
        out_deg0[4 * c + r] = in_deg0[4 * ((c + r) % Nst) + r];
      } else {
        out_deg2[4 * c + r] = in_deg2[4 * ((c + r + 1) % Nst) + r];
        out_deg1[4 * c + r] = in_deg1[4 * ((c + r + 1) % Nst) + r];
        out_deg0[4 * c + r] = in_deg0[4 * ((c + r + 1) % Nst) + r];
      }
    }
  }
}

static void aes_192_shiftrows_prover(bf192_t* out_deg0, bf192_t* out_deg1, bf192_t* out_deg2,
                                     const bf192_t* in_deg0, const bf192_t* in_deg1,
                                     const bf192_t* in_deg2, const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if ((Nst != 8) || (r <= 1)) {
        out_deg2[4 * c + r] = in_deg2[4 * ((c + r) % Nst) + r];
        out_deg1[4 * c + r] = in_deg1[4 * ((c + r) % Nst) + r];
        out_deg0[4 * c + r] = in_deg0[4 * ((c + r) % Nst) + r];
      } else {
        out_deg2[4 * c + r] = in_deg2[4 * ((c + r + 1) % Nst) + r];
        out_deg1[4 * c + r] = in_deg1[4 * ((c + r + 1) % Nst) + r];
        out_deg0[4 * c + r] = in_deg0[4 * ((c + r + 1) % Nst) + r];
      }
    }
  }
}

static void aes_256_shiftrows_prover(bf256_t* out_deg0, bf256_t* out_deg1, bf256_t* out_deg2,
                                     const bf256_t* in_deg0, const bf256_t* in_deg1,
                                     const bf256_t* in_deg2, const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if ((Nst != 8) || (r <= 1)) {
        out_deg2[4 * c + r] = in_deg2[4 * ((c + r) % Nst) + r];
        out_deg1[4 * c + r] = in_deg1[4 * ((c + r) % Nst) + r];
        out_deg0[4 * c + r] = in_deg0[4 * ((c + r) % Nst) + r];
      } else {
        out_deg2[4 * c + r] = in_deg2[4 * ((c + r + 1) % Nst) + r];
        out_deg1[4 * c + r] = in_deg1[4 * ((c + r + 1) % Nst) + r];
        out_deg0[4 * c + r] = in_deg0[4 * ((c + r + 1) % Nst) + r];
      }
    }
  }
}

static void aes_128_shiftrows_verifier(bf128_t* out_deg1, const bf128_t* in_deg1,
                                       const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if ((Nst != 8) || (r <= 1)) {
        out_deg1[4 * c + r] = in_deg1[(4 * ((c + r) % Nst) + r)];
      } else {
        out_deg1[4 * c + r] = in_deg1[(4 * ((c + r + 1) % Nst) + r)];
      }
    }
  }
}

static void aes_192_shiftrows_verifier(bf192_t* out_deg1, const bf192_t* in_deg1,
                                       const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if ((Nst != 8) || (r <= 1)) {
        out_deg1[4 * c + r] = in_deg1[(4 * ((c + r) % Nst) + r)];
      } else {
        out_deg1[4 * c + r] = in_deg1[(4 * ((c + r + 1) % Nst) + r)];
      }
    }
  }
}

static void aes_256_shiftrows_verifier(bf256_t* out_deg1, const bf256_t* in_deg1,
                                       const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      if ((Nst != 8) || (r <= 1)) {
        out_deg1[4 * c + r] = in_deg1[(4 * ((c + r) % Nst) + r)];
      } else {
        out_deg1[4 * c + r] = in_deg1[(4 * ((c + r + 1) % Nst) + r)];
      }
    }
  }
}

// MIX COLOUMNS
static void aes_128_mix_columns_prover(bf128_t* y_deg0, bf128_t* y_deg1, bf128_t* y_deg2,
                                       const bf128_t* in_deg0, const bf128_t* in_deg1,
                                       const bf128_t* in_deg2, bool dosq,
                                       const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  //  ::2-4
  bf128_t v2 = bf128_byte_combine_bits(2);
  bf128_t v3 = bf128_byte_combine_bits(3);
  if (dosq) {
    v2 = bf128_mul(v2, v2);
    v3 = bf128_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {
    unsigned int i0 = 4 * c;
    unsigned int i1 = 4 * c + 1;
    unsigned int i2 = 4 * c + 2;
    unsigned int i3 = 4 * c + 3;
    bf128_t tmp1, tmp2, tmp3, tmp4;

    // ::7
    tmp1       = bf128_mul(in_deg2[i0], v2);
    tmp2       = bf128_mul(in_deg2[i1], v3);
    tmp3       = in_deg2[i2];
    tmp4       = in_deg2[i3];
    y_deg2[i0] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    tmp1       = bf128_mul(in_deg1[i0], v2);
    tmp2       = bf128_mul(in_deg1[i1], v3);
    tmp3       = in_deg1[i2];
    tmp4       = in_deg1[i3];
    y_deg1[i0] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    tmp1       = bf128_mul(in_deg0[i0], v2);
    tmp2       = bf128_mul(in_deg0[i1], v3);
    tmp3       = in_deg0[i2];
    tmp4       = in_deg0[i3];
    y_deg0[i0] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    // ::8
    tmp1       = in_deg2[i0];
    tmp2       = bf128_mul(in_deg2[i1], v2);
    tmp3       = bf128_mul(in_deg2[i2], v3);
    tmp4       = in_deg2[i3];
    y_deg2[i1] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    tmp1       = in_deg1[i0];
    tmp2       = bf128_mul(in_deg1[i1], v2);
    tmp3       = bf128_mul(in_deg1[i2], v3);
    tmp4       = in_deg1[i3];
    y_deg1[i1] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    tmp1       = in_deg0[i0];
    tmp2       = bf128_mul(in_deg0[i1], v2);
    tmp3       = bf128_mul(in_deg0[i2], v3);
    tmp4       = in_deg0[i3];
    y_deg0[i1] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    // ::9
    tmp1       = in_deg2[i0];
    tmp2       = in_deg2[i1];
    tmp3       = bf128_mul(in_deg2[i2], v2);
    tmp4       = bf128_mul(in_deg2[i3], v3);
    y_deg2[i2] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    tmp1       = in_deg1[i0];
    tmp2       = in_deg1[i1];
    tmp3       = bf128_mul(in_deg1[i2], v2);
    tmp4       = bf128_mul(in_deg1[i3], v3);
    y_deg1[i2] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    tmp1       = in_deg0[i0];
    tmp2       = in_deg0[i1];
    tmp3       = bf128_mul(in_deg0[i2], v2);
    tmp4       = bf128_mul(in_deg0[i3], v3);
    y_deg0[i2] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    // ::10
    tmp1       = bf128_mul(in_deg2[i0], v3);
    tmp2       = in_deg2[i1];
    tmp3       = in_deg2[i2];
    tmp4       = bf128_mul(in_deg2[i3], v2);
    y_deg2[i3] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    tmp1       = bf128_mul(in_deg1[i0], v3);
    tmp2       = in_deg1[i1];
    tmp3       = in_deg1[i2];
    tmp4       = bf128_mul(in_deg1[i3], v2);
    y_deg1[i3] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));

    tmp1       = bf128_mul(in_deg0[i0], v3);
    tmp2       = in_deg0[i1];
    tmp3       = in_deg0[i2];
    tmp4       = bf128_mul(in_deg0[i3], v2);
    y_deg0[i3] = bf128_add(bf128_add(tmp1, tmp2), bf128_add(tmp3, tmp4));
  }
}

static void aes_192_mix_columns_prover(bf192_t* y_deg0, bf192_t* y_deg1, bf192_t* y_deg2,
                                       const bf192_t* in_deg0, const bf192_t* in_deg1,
                                       const bf192_t* in_deg2, bool dosq,
                                       const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  //  ::2-4
  bf192_t v2 = bf192_byte_combine_bits(2);
  bf192_t v3 = bf192_byte_combine_bits(3);
  if (dosq) {
    v2 = bf192_mul(v2, v2);
    v3 = bf192_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {
    unsigned int i0 = 4 * c;
    unsigned int i1 = 4 * c + 1;
    unsigned int i2 = 4 * c + 2;
    unsigned int i3 = 4 * c + 3;
    bf192_t tmp1, tmp2, tmp3, tmp4;

    // ::7
    tmp1       = bf192_mul(in_deg2[i0], v2);
    tmp2       = bf192_mul(in_deg2[i1], v3);
    tmp3       = in_deg2[i2];
    tmp4       = in_deg2[i3];
    y_deg2[i0] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    tmp1       = bf192_mul(in_deg1[i0], v2);
    tmp2       = bf192_mul(in_deg1[i1], v3);
    tmp3       = in_deg1[i2];
    tmp4       = in_deg1[i3];
    y_deg1[i0] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    tmp1       = bf192_mul(in_deg0[i0], v2);
    tmp2       = bf192_mul(in_deg0[i1], v3);
    tmp3       = in_deg0[i2];
    tmp4       = in_deg0[i3];
    y_deg0[i0] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    // ::8
    tmp1       = in_deg2[i0];
    tmp2       = bf192_mul(in_deg2[i1], v2);
    tmp3       = bf192_mul(in_deg2[i2], v3);
    tmp4       = in_deg2[i3];
    y_deg2[i1] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    tmp1       = in_deg1[i0];
    tmp2       = bf192_mul(in_deg1[i1], v2);
    tmp3       = bf192_mul(in_deg1[i2], v3);
    tmp4       = in_deg1[i3];
    y_deg1[i1] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    tmp1       = in_deg0[i0];
    tmp2       = bf192_mul(in_deg0[i1], v2);
    tmp3       = bf192_mul(in_deg0[i2], v3);
    tmp4       = in_deg0[i3];
    y_deg0[i1] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    // ::9
    tmp1       = in_deg2[i0];
    tmp2       = in_deg2[i1];
    tmp3       = bf192_mul(in_deg2[i2], v2);
    tmp4       = bf192_mul(in_deg2[i3], v3);
    y_deg2[i2] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    tmp1       = in_deg1[i0];
    tmp2       = in_deg1[i1];
    tmp3       = bf192_mul(in_deg1[i2], v2);
    tmp4       = bf192_mul(in_deg1[i3], v3);
    y_deg1[i2] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    tmp1       = in_deg0[i0];
    tmp2       = in_deg0[i1];
    tmp3       = bf192_mul(in_deg0[i2], v2);
    tmp4       = bf192_mul(in_deg0[i3], v3);
    y_deg0[i2] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    // ::10
    tmp1       = bf192_mul(in_deg2[i0], v3);
    tmp2       = in_deg2[i1];
    tmp3       = in_deg2[i2];
    tmp4       = bf192_mul(in_deg2[i3], v2);
    y_deg2[i3] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    tmp1       = bf192_mul(in_deg1[i0], v3);
    tmp2       = in_deg1[i1];
    tmp3       = in_deg1[i2];
    tmp4       = bf192_mul(in_deg1[i3], v2);
    y_deg1[i3] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));

    tmp1       = bf192_mul(in_deg0[i0], v3);
    tmp2       = in_deg0[i1];
    tmp3       = in_deg0[i2];
    tmp4       = bf192_mul(in_deg0[i3], v2);
    y_deg0[i3] = bf192_add(bf192_add(tmp1, tmp2), bf192_add(tmp3, tmp4));
  }
}

static void aes_256_mix_columns_prover(bf256_t* y_deg0, bf256_t* y_deg1, bf256_t* y_deg2,
                                       const bf256_t* in_deg0, const bf256_t* in_deg1,
                                       const bf256_t* in_deg2, bool dosq,
                                       const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  //  ::2-4
  bf256_t v2 = bf256_byte_combine_bits(2);
  bf256_t v3 = bf256_byte_combine_bits(3);
  if (dosq) {
    v2 = bf256_mul(v2, v2);
    v3 = bf256_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {
    unsigned int i0 = 4 * c;
    unsigned int i1 = 4 * c + 1;
    unsigned int i2 = 4 * c + 2;
    unsigned int i3 = 4 * c + 3;
    bf256_t tmp1, tmp2, tmp3, tmp4;

    // ::7
    tmp1       = bf256_mul(in_deg2[i0], v2);
    tmp2       = bf256_mul(in_deg2[i1], v3);
    tmp3       = in_deg2[i2];
    tmp4       = in_deg2[i3];
    y_deg2[i0] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    tmp1       = bf256_mul(in_deg1[i0], v2);
    tmp2       = bf256_mul(in_deg1[i1], v3);
    tmp3       = in_deg1[i2];
    tmp4       = in_deg1[i3];
    y_deg1[i0] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    tmp1       = bf256_mul(in_deg0[i0], v2);
    tmp2       = bf256_mul(in_deg0[i1], v3);
    tmp3       = in_deg0[i2];
    tmp4       = in_deg0[i3];
    y_deg0[i0] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    // ::8
    tmp1       = in_deg2[i0];
    tmp2       = bf256_mul(in_deg2[i1], v2);
    tmp3       = bf256_mul(in_deg2[i2], v3);
    tmp4       = in_deg2[i3];
    y_deg2[i1] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    tmp1       = in_deg1[i0];
    tmp2       = bf256_mul(in_deg1[i1], v2);
    tmp3       = bf256_mul(in_deg1[i2], v3);
    tmp4       = in_deg1[i3];
    y_deg1[i1] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    tmp1       = in_deg0[i0];
    tmp2       = bf256_mul(in_deg0[i1], v2);
    tmp3       = bf256_mul(in_deg0[i2], v3);
    tmp4       = in_deg0[i3];
    y_deg0[i1] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    // ::9
    tmp1       = in_deg2[i0];
    tmp2       = in_deg2[i1];
    tmp3       = bf256_mul(in_deg2[i2], v2);
    tmp4       = bf256_mul(in_deg2[i3], v3);
    y_deg2[i2] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    tmp1       = in_deg1[i0];
    tmp2       = in_deg1[i1];
    tmp3       = bf256_mul(in_deg1[i2], v2);
    tmp4       = bf256_mul(in_deg1[i3], v3);
    y_deg1[i2] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    tmp1       = in_deg0[i0];
    tmp2       = in_deg0[i1];
    tmp3       = bf256_mul(in_deg0[i2], v2);
    tmp4       = bf256_mul(in_deg0[i3], v3);
    y_deg0[i2] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    // ::10
    tmp1       = bf256_mul(in_deg2[i0], v3);
    tmp2       = in_deg2[i1];
    tmp3       = in_deg2[i2];
    tmp4       = bf256_mul(in_deg2[i3], v2);
    y_deg2[i3] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    tmp1       = bf256_mul(in_deg1[i0], v3);
    tmp2       = in_deg1[i1];
    tmp3       = in_deg1[i2];
    tmp4       = bf256_mul(in_deg1[i3], v2);
    y_deg1[i3] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));

    tmp1       = bf256_mul(in_deg0[i0], v3);
    tmp2       = in_deg0[i1];
    tmp3       = in_deg0[i2];
    tmp4       = bf256_mul(in_deg0[i3], v2);
    y_deg0[i3] = bf256_add(bf256_add(tmp1, tmp2), bf256_add(tmp3, tmp4));
  }
}

static void aes_128_mix_columns_verifier(bf128_t* y_deg1, const bf128_t* in_deg1, bool dosq,
                                         const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  //  ::2-4
  bf128_t v2 = bf128_byte_combine_bits(2);
  bf128_t v3 = bf128_byte_combine_bits(3);
  if (dosq) {
    v2 = bf128_mul(v2, v2);
    v3 = bf128_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {
    unsigned int i0 = 4 * c;
    unsigned int i1 = 4 * c + 1;
    unsigned int i2 = 4 * c + 2;
    unsigned int i3 = 4 * c + 3;

    bf128_t tmp1_tag = bf128_mul(in_deg1[i0], v2);
    bf128_t tmp2_tag = bf128_mul(in_deg1[i1], v3);
    bf128_t tmp3_tag = in_deg1[i2];
    bf128_t tmp4_tag = in_deg1[i3];
    y_deg1[i0]       = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    tmp1_tag   = in_deg1[i0];
    tmp2_tag   = bf128_mul(in_deg1[i1], v2);
    tmp3_tag   = bf128_mul(in_deg1[i2], v3);
    tmp4_tag   = in_deg1[i3];
    y_deg1[i1] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    tmp1_tag   = in_deg1[i0];
    tmp2_tag   = in_deg1[i1];
    tmp3_tag   = bf128_mul(in_deg1[i2], v2);
    tmp4_tag   = bf128_mul(in_deg1[i3], v3);
    y_deg1[i2] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));

    tmp1_tag   = bf128_mul(in_deg1[i0], v3);
    tmp2_tag   = in_deg1[i1];
    tmp3_tag   = in_deg1[i2];
    tmp4_tag   = bf128_mul(in_deg1[i3], v2);
    y_deg1[i3] = bf128_add(bf128_add(tmp1_tag, tmp2_tag), bf128_add(tmp3_tag, tmp4_tag));
  }
}

static void aes_192_mix_columns_verifier(bf192_t* y_deg1, const bf192_t* in_deg1, bool dosq,
                                         const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  //  ::2-4
  bf192_t v2 = bf192_byte_combine_bits(2);
  bf192_t v3 = bf192_byte_combine_bits(3);
  if (dosq) {
    v2 = bf192_mul(v2, v2);
    v3 = bf192_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {
    unsigned int i0 = 4 * c;
    unsigned int i1 = 4 * c + 1;
    unsigned int i2 = 4 * c + 2;
    unsigned int i3 = 4 * c + 3;

    bf192_t tmp1_tag = bf192_mul(in_deg1[i0], v2);
    bf192_t tmp2_tag = bf192_mul(in_deg1[i1], v3);
    bf192_t tmp3_tag = in_deg1[i2];
    bf192_t tmp4_tag = in_deg1[i3];
    y_deg1[i0]       = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));

    tmp1_tag   = in_deg1[i0];
    tmp2_tag   = bf192_mul(in_deg1[i1], v2);
    tmp3_tag   = bf192_mul(in_deg1[i2], v3);
    tmp4_tag   = in_deg1[i3];
    y_deg1[i1] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));

    tmp1_tag   = in_deg1[i0];
    tmp2_tag   = in_deg1[i1];
    tmp3_tag   = bf192_mul(in_deg1[i2], v2);
    tmp4_tag   = bf192_mul(in_deg1[i3], v3);
    y_deg1[i2] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));

    tmp1_tag   = bf192_mul(in_deg1[i0], v3);
    tmp2_tag   = in_deg1[i1];
    tmp3_tag   = in_deg1[i2];
    tmp4_tag   = bf192_mul(in_deg1[i3], v2);
    y_deg1[i3] = bf192_add(bf192_add(tmp1_tag, tmp2_tag), bf192_add(tmp3_tag, tmp4_tag));
  }
}

static void aes_256_mix_columns_verifier(bf256_t* y_deg1, const bf256_t* in_deg1, bool dosq,
                                         const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  //  ::2-4
  bf256_t v2 = bf256_byte_combine_bits(2);
  bf256_t v3 = bf256_byte_combine_bits(3);
  if (dosq) {
    v2 = bf256_mul(v2, v2);
    v3 = bf256_mul(v3, v3);
  }

  for (unsigned int c = 0; c < Nst; c++) {
    unsigned int i0 = 4 * c;
    unsigned int i1 = 4 * c + 1;
    unsigned int i2 = 4 * c + 2;
    unsigned int i3 = 4 * c + 3;

    bf256_t tmp1_tag = bf256_mul(in_deg1[i0], v2);
    bf256_t tmp2_tag = bf256_mul(in_deg1[i1], v3);
    bf256_t tmp3_tag = in_deg1[i2];
    bf256_t tmp4_tag = in_deg1[i3];
    y_deg1[i0]       = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));

    tmp1_tag   = in_deg1[i0];
    tmp2_tag   = bf256_mul(in_deg1[i1], v2);
    tmp3_tag   = bf256_mul(in_deg1[i2], v3);
    tmp4_tag   = in_deg1[i3];
    y_deg1[i1] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));

    tmp1_tag   = in_deg1[i0];
    tmp2_tag   = in_deg1[i1];
    tmp3_tag   = bf256_mul(in_deg1[i2], v2);
    tmp4_tag   = bf256_mul(in_deg1[i3], v3);
    y_deg1[i2] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));

    tmp1_tag   = bf256_mul(in_deg1[i0], v3);
    tmp2_tag   = in_deg1[i1];
    tmp3_tag   = in_deg1[i2];
    tmp4_tag   = bf256_mul(in_deg1[i3], v2);
    y_deg1[i3] = bf256_add(bf256_add(tmp1_tag, tmp2_tag), bf256_add(tmp3_tag, tmp4_tag));
  }
}

// ADD ROUND KEY BYTES
// on degree-2 state and degree-2 key
static void aes_128_add_round_key_bytes_prover_degree_2(
    bf128_t* y_deg0, bf128_t* y_deg1, bf128_t* y_deg2, const bf128_t* in_deg0,
    const bf128_t* in_deg1, const bf128_t* in_deg2, const bf128_t* k_deg0, const bf128_t* k_deg1,
    const bf128_t* k_deg2, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    y_deg2[i] = bf128_add(in_deg2[i], k_deg2[i]);
    y_deg1[i] = bf128_add(in_deg1[i], k_deg1[i]);
    y_deg0[i] = bf128_add(in_deg0[i], k_deg0[i]);
  }
}

static void aes_192_add_round_key_bytes_prover_degree_2(
    bf192_t* y_deg0, bf192_t* y_deg1, bf192_t* y_deg2, const bf192_t* in_deg0,
    const bf192_t* in_deg1, const bf192_t* in_deg2, const bf192_t* k_deg0, const bf192_t* k_deg1,
    const bf192_t* k_deg2, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    y_deg2[i] = bf192_add(in_deg2[i], k_deg2[i]);
    y_deg1[i] = bf192_add(in_deg1[i], k_deg1[i]);
    y_deg0[i] = bf192_add(in_deg0[i], k_deg0[i]);
  }
}

static void aes_256_add_round_key_bytes_prover_degree_2(
    bf256_t* y_deg0, bf256_t* y_deg1, bf256_t* y_deg2, const bf256_t* in_deg0,
    const bf256_t* in_deg1, const bf256_t* in_deg2, const bf256_t* k_deg0, const bf256_t* k_deg1,
    const bf256_t* k_deg2, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    y_deg2[i] = bf256_add(in_deg2[i], k_deg2[i]);
    y_deg1[i] = bf256_add(in_deg1[i], k_deg1[i]);
    y_deg0[i] = bf256_add(in_deg0[i], k_deg0[i]);
  }
}

static void
aes_128_add_round_key_bytes_prover_degree_1(bf128_t* y_deg0, bf128_t* y_deg1, bf128_t* y_deg2,
                                            const bf128_t* in_deg0, const bf128_t* in_deg1,
                                            const bf128_t* in_deg2, const bf128_t* k_deg1,
                                            const bf128_t* k_deg2, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    y_deg2[i] = bf128_add(in_deg2[i], k_deg2[i]);
    y_deg1[i] = bf128_add(in_deg1[i], k_deg1[i]);
    y_deg0[i] = in_deg0[i];
  }
}

static void
aes_192_add_round_key_bytes_prover_degree_1(bf192_t* y_deg0, bf192_t* y_deg1, bf192_t* y_deg2,
                                            const bf192_t* in_deg0, const bf192_t* in_deg1,
                                            const bf192_t* in_deg2, const bf192_t* k_deg1,
                                            const bf192_t* k_deg2, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    y_deg2[i] = bf192_add(in_deg2[i], k_deg2[i]);
    y_deg1[i] = bf192_add(in_deg1[i], k_deg1[i]);
    y_deg0[i] = in_deg0[i];
  }
}

static void
aes_256_add_round_key_bytes_prover_degree_1(bf256_t* y_deg0, bf256_t* y_deg1, bf256_t* y_deg2,
                                            const bf256_t* in_deg0, const bf256_t* in_deg1,
                                            const bf256_t* in_deg2, const bf256_t* k_deg1,
                                            const bf256_t* k_deg2, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    y_deg2[i] = bf256_add(in_deg2[i], k_deg2[i]);
    y_deg1[i] = bf256_add(in_deg1[i], k_deg1[i]);
    y_deg0[i] = in_deg0[i];
  }
}

// Use shift_tag if key is degree-1 instead of degree-2
static void aes_128_add_round_key_bytes_verifier(bf128_t* y_deg1, const bf128_t* in_tag,
                                                 const bf128_t* k_tag, bf128_t delta,
                                                 bool shift_tag, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    if (shift_tag) {
      // Multiply tag by delta to align degrees
      y_deg1[i] = bf128_add(in_tag[i], bf128_mul(k_tag[i], delta));
    } else {
      y_deg1[i] = bf128_add(in_tag[i], k_tag[i]);
    }
  }
}

static void aes_192_add_round_key_bytes_verifier(bf192_t* y_deg1, const bf192_t* in_tag,
                                                 const bf192_t* k_tag, bf192_t delta,
                                                 bool shift_tag, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    if (shift_tag) {
      // Multiply tag by delta to align degrees
      y_deg1[i] = bf192_add(in_tag[i], bf192_mul(k_tag[i], delta));
    } else {
      y_deg1[i] = bf192_add(in_tag[i], k_tag[i]);
    }
  }
}

static void aes_256_add_round_key_bytes_verifier(bf256_t* y_deg1, const bf256_t* in_tag,
                                                 const bf256_t* k_tag, bf256_t delta,
                                                 bool shift_tag, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    if (shift_tag) {
      // Multiply tag by delta to align degrees
      y_deg1[i] = bf256_add(in_tag[i], bf256_mul(k_tag[i], delta));
    } else {
      y_deg1[i] = bf256_add(in_tag[i], k_tag[i]);
    }
  }
}

// INVERSE SHIFT ROWS
static void aes_128_inverse_shiftrows_prover(uint8_t* out, bf128_t* out_tag, const uint8_t* in,
                                             const bf128_t* in_tag,
                                             const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      unsigned int i;
      if ((Nst != 8) || (r <= 1)) {
        i = 4 * ((c + Nst - r) % Nst) + r;
      } else {
        i = 4 * ((c + Nst - r - 1) % Nst) + r;
      }
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        out[8 * (4 * c + r) + bit_i]     = in[8 * i + bit_i]; // moving bitwise
        out_tag[8 * (4 * c + r) + bit_i] = in_tag[8 * i + bit_i];
      }
    }
  }
}

static void aes_192_inverse_shiftrows_prover(uint8_t* out, bf192_t* out_tag, const uint8_t* in,
                                             const bf192_t* in_tag,
                                             const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      unsigned int i;
      if ((Nst != 8) || (r <= 1)) {
        i = 4 * ((c + Nst - r) % Nst) + r;
      } else {
        i = 4 * ((c + Nst - r - 1) % Nst) + r;
      }

      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out[8 * (4 * c + r) + byte_idx]     = in[8 * i + byte_idx];
        out_tag[8 * (4 * c + r) + byte_idx] = in_tag[8 * i + byte_idx];
      }
    }
  }
}

static void aes_256_inverse_shiftrows_prover(uint8_t* out, bf256_t* out_tag, const uint8_t* in,
                                             const bf256_t* in_tag,
                                             const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      unsigned int i;
      if ((Nst != 8) || (r <= 1)) {
        i = 4 * ((c + Nst - r) % Nst) + r;
      } else {
        i = 4 * ((c + Nst - r - 1) % Nst) + r;
      }

      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out[8 * (4 * c + r) + byte_idx]     = in[8 * i + byte_idx];
        out_tag[8 * (4 * c + r) + byte_idx] = in_tag[8 * i + byte_idx];
      }
    }
  }
}

static void aes_128_inverse_shiftrows_verifier(bf128_t* out_tag, const bf128_t* in_tag,
                                               const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      unsigned int i;
      if ((Nst != 8) || (r <= 1)) {
        i = 4 * ((c + Nst - r) % Nst) + r;
      } else {
        i = 4 * ((c + Nst - r - 1) % Nst) + r;
      }

      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out_tag[8 * (4 * c + r) + byte_idx] = in_tag[8 * i + byte_idx];
      }
    }
  }
}

static void aes_192_inverse_shiftrows_verifier(bf192_t* out_tag, const bf192_t* in_tag,
                                               const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      unsigned int i;
      if ((Nst != 8) || (r <= 1)) {
        i = 4 * ((c + Nst - r) % Nst) + r;
      } else {
        i = 4 * ((c + Nst - r - 1) % Nst) + r;
      }

      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out_tag[8 * (4 * c + r) + byte_idx] = in_tag[8 * i + byte_idx];
      }
    }
  }
}

static void aes_256_inverse_shiftrows_verifier(bf256_t* out_tag, const bf256_t* in_tag,
                                               const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int r = 0; r < 4; r++) {
    for (unsigned int c = 0; c < Nst; c++) {
      unsigned int i;
      if ((Nst != 8) || (r <= 1)) {
        i = 4 * ((c + Nst - r) % Nst) + r;
      } else {
        i = 4 * ((c + Nst - r - 1) % Nst) + r;
      }

      for (unsigned int byte_idx = 0; byte_idx < 8; byte_idx++) {
        out_tag[8 * (4 * c + r) + byte_idx] = in_tag[8 * i + byte_idx];
      }
    }
  }
}

// BITWISE MIX COLUMNS
static void aes_128_bitwise_mix_column_prover(uint8_t* out, bf128_t* out_tag, uint8_t* s,
                                              bf128_t* s_tag, const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int c = 0; c < Nst; c++) {
    uint8_t a_bits[4 * 8];
    bf128_t a_bits_tag[4 * 8];

    uint8_t b_bits[4 * 8];
    bf128_t b_bits_tag[4 * 8];

    // ::1
    for (unsigned int r = 0; r < 4; r++) {
      // :2
      // assign a^(r), vector of 8 bits
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits[r * 8 + bit_i]     = s[32 * c + 8 * r + bit_i];
        a_bits_tag[r * 8 + bit_i] = s_tag[32 * c + 8 * r + bit_i];
      }
      // :5
      // b^(r), vector of 8 bits
      b_bits[r * 8 + 0] = a_bits[r * 8 + 7];
      b_bits[r * 8 + 1] = a_bits[r * 8 + 0] ^ a_bits[r * 8 + 7];
      b_bits[r * 8 + 2] = a_bits[r * 8 + 1];
      b_bits[r * 8 + 3] = a_bits[r * 8 + 2] ^ a_bits[r * 8 + 7];
      b_bits[r * 8 + 4] = a_bits[r * 8 + 3] ^ a_bits[r * 8 + 7];
      b_bits[r * 8 + 5] = a_bits[r * 8 + 4];
      b_bits[r * 8 + 6] = a_bits[r * 8 + 5];
      b_bits[r * 8 + 7] = a_bits[r * 8 + 6];

      b_bits_tag[r * 8 + 0] = a_bits_tag[r * 8 + 7];
      b_bits_tag[r * 8 + 1] = bf128_add(a_bits_tag[r * 8 + 0], a_bits_tag[r * 8 + 7]);
      b_bits_tag[r * 8 + 2] = a_bits_tag[r * 8 + 1];
      b_bits_tag[r * 8 + 3] = bf128_add(a_bits_tag[r * 8 + 2], a_bits_tag[r * 8 + 7]);
      b_bits_tag[r * 8 + 4] = bf128_add(a_bits_tag[r * 8 + 3], a_bits_tag[r * 8 + 7]);
      b_bits_tag[r * 8 + 5] = a_bits_tag[r * 8 + 4];
      b_bits_tag[r * 8 + 6] = a_bits_tag[r * 8 + 5];
      b_bits_tag[r * 8 + 7] = a_bits_tag[r * 8 + 6];
    }
    for (uint16_t i_bit = 0; i_bit < 8; ++i_bit) {
      // out[c*4] = b_bits[0] ^ a_bits[3] ^ a_bits[2] ^ b_bits[1] ^ a_bits[1];
      out[8 * (c * 4) + i_bit] = b_bits[0 * 8 + i_bit] ^ a_bits[3 * 8 + i_bit] ^
                                 a_bits[2 * 8 + i_bit] ^ b_bits[1 * 8 + i_bit] ^
                                 a_bits[1 * 8 + i_bit];
      // out[c*4 + 1] = b_bits[1] ^ a_bits[0] ^ a_bits[3] ^ b_bits[2] ^ a_bits[2];
      out[8 * (c * 4 + 1) + i_bit] = b_bits[1 * 8 + i_bit] ^ a_bits[0 * 8 + i_bit] ^
                                     a_bits[3 * 8 + i_bit] ^ b_bits[2 * 8 + i_bit] ^
                                     a_bits[2 * 8 + i_bit];

      out[8 * (c * 4 + 2) + i_bit] = b_bits[2 * 8 + i_bit] ^ a_bits[1 * 8 + i_bit] ^
                                     a_bits[0 * 8 + i_bit] ^ b_bits[3 * 8 + i_bit] ^
                                     a_bits[3 * 8 + i_bit];

      out[8 * (c * 4 + 3) + i_bit] = b_bits[3 * 8 + i_bit] ^ a_bits[2 * 8 + i_bit] ^
                                     a_bits[1 * 8 + i_bit] ^ b_bits[0 * 8 + i_bit] ^
                                     a_bits[0 * 8 + i_bit];

      out_tag[8 * (c * 4) + i_bit] =
          bf128_add(bf128_add(bf128_add(b_bits_tag[0 * 8 + i_bit], a_bits_tag[3 * 8 + i_bit]),
                              bf128_add(a_bits_tag[2 * 8 + i_bit], b_bits_tag[1 * 8 + i_bit])),
                    a_bits_tag[1 * 8 + i_bit]);
      out_tag[8 * (c * 4 + 1) + i_bit] =
          bf128_add(bf128_add(bf128_add(b_bits_tag[1 * 8 + i_bit], a_bits_tag[0 * 8 + i_bit]),
                              bf128_add(a_bits_tag[3 * 8 + i_bit], b_bits_tag[2 * 8 + i_bit])),
                    a_bits_tag[2 * 8 + i_bit]);
      out_tag[8 * (c * 4 + 2) + i_bit] =
          bf128_add(bf128_add(bf128_add(b_bits_tag[2 * 8 + i_bit], a_bits_tag[1 * 8 + i_bit]),
                              bf128_add(a_bits_tag[0 * 8 + i_bit], b_bits_tag[3 * 8 + i_bit])),
                    a_bits_tag[3 * 8 + i_bit]);
      out_tag[8 * (c * 4 + 3) + i_bit] =
          bf128_add(bf128_add(bf128_add(b_bits_tag[3 * 8 + i_bit], a_bits_tag[2 * 8 + i_bit]),
                              bf128_add(a_bits_tag[1 * 8 + i_bit], b_bits_tag[0 * 8 + i_bit])),
                    a_bits_tag[0 * 8 + i_bit]);
    }
  }
}

static void aes_192_bitwise_mix_column_prover(uint8_t* out, bf192_t* out_tag, uint8_t* s,
                                              bf192_t* s_tag, const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int c = 0; c < Nst; c++) {
    uint8_t a_bits[4 * 8];
    bf192_t a_bits_tag[4 * 8];

    uint8_t b_bits[4 * 8];
    bf192_t b_bits_tag[4 * 8];

    // ::1
    for (unsigned int r = 0; r < 4; r++) {
      // :2
      // assign a^(r), vector of 8 bits
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits[r * 8 + bit_i]     = s[32 * c + 8 * r + bit_i];
        a_bits_tag[r * 8 + bit_i] = s_tag[32 * c + 8 * r + bit_i];
      }
      // :5
      // b^(r), vector of 8 bits
      b_bits[r * 8 + 0] = a_bits[r * 8 + 7];
      b_bits[r * 8 + 1] = a_bits[r * 8 + 0] ^ a_bits[r * 8 + 7];
      b_bits[r * 8 + 2] = a_bits[r * 8 + 1];
      b_bits[r * 8 + 3] = a_bits[r * 8 + 2] ^ a_bits[r * 8 + 7];
      b_bits[r * 8 + 4] = a_bits[r * 8 + 3] ^ a_bits[r * 8 + 7];
      b_bits[r * 8 + 5] = a_bits[r * 8 + 4];
      b_bits[r * 8 + 6] = a_bits[r * 8 + 5];
      b_bits[r * 8 + 7] = a_bits[r * 8 + 6];

      b_bits_tag[r * 8 + 0] = a_bits_tag[r * 8 + 7];
      b_bits_tag[r * 8 + 1] = bf192_add(a_bits_tag[r * 8 + 0], a_bits_tag[r * 8 + 7]);
      b_bits_tag[r * 8 + 2] = a_bits_tag[r * 8 + 1];
      b_bits_tag[r * 8 + 3] = bf192_add(a_bits_tag[r * 8 + 2], a_bits_tag[r * 8 + 7]);
      b_bits_tag[r * 8 + 4] = bf192_add(a_bits_tag[r * 8 + 3], a_bits_tag[r * 8 + 7]);
      b_bits_tag[r * 8 + 5] = a_bits_tag[r * 8 + 4];
      b_bits_tag[r * 8 + 6] = a_bits_tag[r * 8 + 5];
      b_bits_tag[r * 8 + 7] = a_bits_tag[r * 8 + 6];
    }
    for (uint16_t i_bit = 0; i_bit < 8; ++i_bit) {
      // out[c*4] = b_bits[0] ^ a_bits[3] ^ a_bits[2] ^ b_bits[1] ^ a_bits[1];
      out[8 * (c * 4) + i_bit] = b_bits[0 * 8 + i_bit] ^ a_bits[3 * 8 + i_bit] ^
                                 a_bits[2 * 8 + i_bit] ^ b_bits[1 * 8 + i_bit] ^
                                 a_bits[1 * 8 + i_bit];
      // out[c*4 + 1] = b_bits[1] ^ a_bits[0] ^ a_bits[3] ^ b_bits[2] ^ a_bits[2];
      out[8 * (c * 4 + 1) + i_bit] = b_bits[1 * 8 + i_bit] ^ a_bits[0 * 8 + i_bit] ^
                                     a_bits[3 * 8 + i_bit] ^ b_bits[2 * 8 + i_bit] ^
                                     a_bits[2 * 8 + i_bit];

      out[8 * (c * 4 + 2) + i_bit] = b_bits[2 * 8 + i_bit] ^ a_bits[1 * 8 + i_bit] ^
                                     a_bits[0 * 8 + i_bit] ^ b_bits[3 * 8 + i_bit] ^
                                     a_bits[3 * 8 + i_bit];

      out[8 * (c * 4 + 3) + i_bit] = b_bits[3 * 8 + i_bit] ^ a_bits[2 * 8 + i_bit] ^
                                     a_bits[1 * 8 + i_bit] ^ b_bits[0 * 8 + i_bit] ^
                                     a_bits[0 * 8 + i_bit];

      out_tag[8 * (c * 4) + i_bit] =
          bf192_add(bf192_add(bf192_add(b_bits_tag[0 * 8 + i_bit], a_bits_tag[3 * 8 + i_bit]),
                              bf192_add(a_bits_tag[2 * 8 + i_bit], b_bits_tag[1 * 8 + i_bit])),
                    a_bits_tag[1 * 8 + i_bit]);
      out_tag[8 * (c * 4 + 1) + i_bit] =
          bf192_add(bf192_add(bf192_add(b_bits_tag[1 * 8 + i_bit], a_bits_tag[0 * 8 + i_bit]),
                              bf192_add(a_bits_tag[3 * 8 + i_bit], b_bits_tag[2 * 8 + i_bit])),
                    a_bits_tag[2 * 8 + i_bit]);
      out_tag[8 * (c * 4 + 2) + i_bit] =
          bf192_add(bf192_add(bf192_add(b_bits_tag[2 * 8 + i_bit], a_bits_tag[1 * 8 + i_bit]),
                              bf192_add(a_bits_tag[0 * 8 + i_bit], b_bits_tag[3 * 8 + i_bit])),
                    a_bits_tag[3 * 8 + i_bit]);
      out_tag[8 * (c * 4 + 3) + i_bit] =
          bf192_add(bf192_add(bf192_add(b_bits_tag[3 * 8 + i_bit], a_bits_tag[2 * 8 + i_bit]),
                              bf192_add(a_bits_tag[1 * 8 + i_bit], b_bits_tag[0 * 8 + i_bit])),
                    a_bits_tag[0 * 8 + i_bit]);
    }
  }
}

static void aes_256_bitwise_mix_column_prover(uint8_t* out, bf256_t* out_tag, uint8_t* s,
                                              bf256_t* s_tag, const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int c = 0; c < Nst; c++) {
    uint8_t a_bits[4 * 8];
    bf256_t a_bits_tag[4 * 8];

    uint8_t b_bits[4 * 8];
    bf256_t b_bits_tag[4 * 8];

    // ::1
    for (unsigned int r = 0; r < 4; r++) {
      // :2
      // assign a^(r), vector of 8 bits
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits[r * 8 + bit_i]     = s[32 * c + 8 * r + bit_i];
        a_bits_tag[r * 8 + bit_i] = s_tag[32 * c + 8 * r + bit_i];
      }
      // :5
      // b^(r), vector of 8 bits
      b_bits[r * 8 + 0] = a_bits[r * 8 + 7];
      b_bits[r * 8 + 1] = a_bits[r * 8 + 0] ^ a_bits[r * 8 + 7];
      b_bits[r * 8 + 2] = a_bits[r * 8 + 1];
      b_bits[r * 8 + 3] = a_bits[r * 8 + 2] ^ a_bits[r * 8 + 7];
      b_bits[r * 8 + 4] = a_bits[r * 8 + 3] ^ a_bits[r * 8 + 7];
      b_bits[r * 8 + 5] = a_bits[r * 8 + 4];
      b_bits[r * 8 + 6] = a_bits[r * 8 + 5];
      b_bits[r * 8 + 7] = a_bits[r * 8 + 6];

      b_bits_tag[r * 8 + 0] = a_bits_tag[r * 8 + 7];
      b_bits_tag[r * 8 + 1] = bf256_add(a_bits_tag[r * 8 + 0], a_bits_tag[r * 8 + 7]);
      b_bits_tag[r * 8 + 2] = a_bits_tag[r * 8 + 1];
      b_bits_tag[r * 8 + 3] = bf256_add(a_bits_tag[r * 8 + 2], a_bits_tag[r * 8 + 7]);
      b_bits_tag[r * 8 + 4] = bf256_add(a_bits_tag[r * 8 + 3], a_bits_tag[r * 8 + 7]);
      b_bits_tag[r * 8 + 5] = a_bits_tag[r * 8 + 4];
      b_bits_tag[r * 8 + 6] = a_bits_tag[r * 8 + 5];
      b_bits_tag[r * 8 + 7] = a_bits_tag[r * 8 + 6];
    }
    for (uint16_t i_bit = 0; i_bit < 8; ++i_bit) {
      // out[c*4] = b_bits[0] ^ a_bits[3] ^ a_bits[2] ^ b_bits[1] ^ a_bits[1];
      out[8 * (c * 4) + i_bit] = b_bits[0 * 8 + i_bit] ^ a_bits[3 * 8 + i_bit] ^
                                 a_bits[2 * 8 + i_bit] ^ b_bits[1 * 8 + i_bit] ^
                                 a_bits[1 * 8 + i_bit];
      // out[c*4 + 1] = b_bits[1] ^ a_bits[0] ^ a_bits[3] ^ b_bits[2] ^ a_bits[2];
      out[8 * (c * 4 + 1) + i_bit] = b_bits[1 * 8 + i_bit] ^ a_bits[0 * 8 + i_bit] ^
                                     a_bits[3 * 8 + i_bit] ^ b_bits[2 * 8 + i_bit] ^
                                     a_bits[2 * 8 + i_bit];

      out[8 * (c * 4 + 2) + i_bit] = b_bits[2 * 8 + i_bit] ^ a_bits[1 * 8 + i_bit] ^
                                     a_bits[0 * 8 + i_bit] ^ b_bits[3 * 8 + i_bit] ^
                                     a_bits[3 * 8 + i_bit];

      out[8 * (c * 4 + 3) + i_bit] = b_bits[3 * 8 + i_bit] ^ a_bits[2 * 8 + i_bit] ^
                                     a_bits[1 * 8 + i_bit] ^ b_bits[0 * 8 + i_bit] ^
                                     a_bits[0 * 8 + i_bit];

      out_tag[8 * (c * 4) + i_bit] =
          bf256_add(bf256_add(bf256_add(b_bits_tag[0 * 8 + i_bit], a_bits_tag[3 * 8 + i_bit]),
                              bf256_add(a_bits_tag[2 * 8 + i_bit], b_bits_tag[1 * 8 + i_bit])),
                    a_bits_tag[1 * 8 + i_bit]);
      out_tag[8 * (c * 4 + 1) + i_bit] =
          bf256_add(bf256_add(bf256_add(b_bits_tag[1 * 8 + i_bit], a_bits_tag[0 * 8 + i_bit]),
                              bf256_add(a_bits_tag[3 * 8 + i_bit], b_bits_tag[2 * 8 + i_bit])),
                    a_bits_tag[2 * 8 + i_bit]);
      out_tag[8 * (c * 4 + 2) + i_bit] =
          bf256_add(bf256_add(bf256_add(b_bits_tag[2 * 8 + i_bit], a_bits_tag[1 * 8 + i_bit]),
                              bf256_add(a_bits_tag[0 * 8 + i_bit], b_bits_tag[3 * 8 + i_bit])),
                    a_bits_tag[3 * 8 + i_bit]);
      out_tag[8 * (c * 4 + 3) + i_bit] =
          bf256_add(bf256_add(bf256_add(b_bits_tag[3 * 8 + i_bit], a_bits_tag[2 * 8 + i_bit]),
                              bf256_add(a_bits_tag[1 * 8 + i_bit], b_bits_tag[0 * 8 + i_bit])),
                    a_bits_tag[0 * 8 + i_bit]);
    }
  }
}

static void aes_128_bitwise_mix_column_verifier(bf128_t* out_key, bf128_t* s_keys_tag,
                                                const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int c = 0; c < Nst; c++) {
    bf128_t a_bits_key[4 * 8];
    bf128_t b_bits_key[4 * 8];

    // ::1
    for (unsigned int r = 0; r < 4; r++) {
      // :2
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits_key[r * 8 + bit_i] = s_keys_tag[32 * c + 8 * r + bit_i];
      }
      // :5
      b_bits_key[r * 8 + 0] = a_bits_key[r * 8 + 7];
      b_bits_key[r * 8 + 1] = bf128_add(a_bits_key[r * 8 + 0], a_bits_key[r * 8 + 7]);
      b_bits_key[r * 8 + 2] = a_bits_key[r * 8 + 1];
      b_bits_key[r * 8 + 3] = bf128_add(a_bits_key[r * 8 + 2], a_bits_key[r * 8 + 7]);
      b_bits_key[r * 8 + 4] = bf128_add(a_bits_key[r * 8 + 3], a_bits_key[r * 8 + 7]);
      b_bits_key[r * 8 + 5] = a_bits_key[r * 8 + 4];
      b_bits_key[r * 8 + 6] = a_bits_key[r * 8 + 5];
      b_bits_key[r * 8 + 7] = a_bits_key[r * 8 + 6];
    }

    // ::6-9
    for (uint16_t i_bit = 0; i_bit < 8; ++i_bit) {
      out_key[8 * (c * 4) + i_bit] =
          bf128_add(bf128_add(bf128_add(b_bits_key[0 * 8 + i_bit], a_bits_key[3 * 8 + i_bit]),
                              bf128_add(a_bits_key[2 * 8 + i_bit], b_bits_key[1 * 8 + i_bit])),
                    a_bits_key[1 * 8 + i_bit]);
      out_key[8 * (c * 4 + 1) + i_bit] =
          bf128_add(bf128_add(bf128_add(b_bits_key[1 * 8 + i_bit], a_bits_key[0 * 8 + i_bit]),
                              bf128_add(a_bits_key[3 * 8 + i_bit], b_bits_key[2 * 8 + i_bit])),
                    a_bits_key[2 * 8 + i_bit]);
      out_key[8 * (c * 4 + 2) + i_bit] =
          bf128_add(bf128_add(bf128_add(b_bits_key[2 * 8 + i_bit], a_bits_key[1 * 8 + i_bit]),
                              bf128_add(a_bits_key[0 * 8 + i_bit], b_bits_key[3 * 8 + i_bit])),
                    a_bits_key[3 * 8 + i_bit]);
      out_key[8 * (c * 4 + 3) + i_bit] =
          bf128_add(bf128_add(bf128_add(b_bits_key[3 * 8 + i_bit], a_bits_key[2 * 8 + i_bit]),
                              bf128_add(a_bits_key[1 * 8 + i_bit], b_bits_key[0 * 8 + i_bit])),
                    a_bits_key[0 * 8 + i_bit]);
    }
  }
}

static void aes_192_bitwise_mix_column_verifier(bf192_t* out_key, bf192_t* s_keys_tag,
                                                const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int c = 0; c < Nst; c++) {
    bf192_t a_bits_key[4 * 8];
    bf192_t b_bits_key[4 * 8];

    // ::1
    for (unsigned int r = 0; r < 4; r++) {
      // :2
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits_key[r * 8 + bit_i] = s_keys_tag[32 * c + 8 * r + bit_i];
      }
      // :5
      b_bits_key[r * 8 + 0] = a_bits_key[r * 8 + 7];
      b_bits_key[r * 8 + 1] = bf192_add(a_bits_key[r * 8 + 0], a_bits_key[r * 8 + 7]);
      b_bits_key[r * 8 + 2] = a_bits_key[r * 8 + 1];
      b_bits_key[r * 8 + 3] = bf192_add(a_bits_key[r * 8 + 2], a_bits_key[r * 8 + 7]);
      b_bits_key[r * 8 + 4] = bf192_add(a_bits_key[r * 8 + 3], a_bits_key[r * 8 + 7]);
      b_bits_key[r * 8 + 5] = a_bits_key[r * 8 + 4];
      b_bits_key[r * 8 + 6] = a_bits_key[r * 8 + 5];
      b_bits_key[r * 8 + 7] = a_bits_key[r * 8 + 6];
    }

    // ::6-9
    for (uint16_t i_bit = 0; i_bit < 8; ++i_bit) {
      out_key[8 * (c * 4) + i_bit] =
          bf192_add(bf192_add(bf192_add(b_bits_key[0 * 8 + i_bit], a_bits_key[3 * 8 + i_bit]),
                              bf192_add(a_bits_key[2 * 8 + i_bit], b_bits_key[1 * 8 + i_bit])),
                    a_bits_key[1 * 8 + i_bit]);
      out_key[8 * (c * 4 + 1) + i_bit] =
          bf192_add(bf192_add(bf192_add(b_bits_key[1 * 8 + i_bit], a_bits_key[0 * 8 + i_bit]),
                              bf192_add(a_bits_key[3 * 8 + i_bit], b_bits_key[2 * 8 + i_bit])),
                    a_bits_key[2 * 8 + i_bit]);
      out_key[8 * (c * 4 + 2) + i_bit] =
          bf192_add(bf192_add(bf192_add(b_bits_key[2 * 8 + i_bit], a_bits_key[1 * 8 + i_bit]),
                              bf192_add(a_bits_key[0 * 8 + i_bit], b_bits_key[3 * 8 + i_bit])),
                    a_bits_key[3 * 8 + i_bit]);
      out_key[8 * (c * 4 + 3) + i_bit] =
          bf192_add(bf192_add(bf192_add(b_bits_key[3 * 8 + i_bit], a_bits_key[2 * 8 + i_bit]),
                              bf192_add(a_bits_key[1 * 8 + i_bit], b_bits_key[0 * 8 + i_bit])),
                    a_bits_key[0 * 8 + i_bit]);
    }
  }
}

static void aes_256_bitwise_mix_column_verifier(bf256_t* out_key, bf256_t* s_keys_tag,
                                                const faest_paramset_t* params) {
  const unsigned int Nst = params->Nst;

  for (unsigned int c = 0; c < Nst; c++) {
    bf256_t a_bits_key[4 * 8];
    bf256_t b_bits_key[4 * 8];

    // ::1
    for (unsigned int r = 0; r < 4; r++) {
      // :2
      for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
        // :3
        a_bits_key[r * 8 + bit_i] = s_keys_tag[32 * c + 8 * r + bit_i];
      }
      // :5
      b_bits_key[r * 8 + 0] = a_bits_key[r * 8 + 7];
      b_bits_key[r * 8 + 1] = bf256_add(a_bits_key[r * 8 + 0], a_bits_key[r * 8 + 7]);
      b_bits_key[r * 8 + 2] = a_bits_key[r * 8 + 1];
      b_bits_key[r * 8 + 3] = bf256_add(a_bits_key[r * 8 + 2], a_bits_key[r * 8 + 7]);
      b_bits_key[r * 8 + 4] = bf256_add(a_bits_key[r * 8 + 3], a_bits_key[r * 8 + 7]);
      b_bits_key[r * 8 + 5] = a_bits_key[r * 8 + 4];
      b_bits_key[r * 8 + 6] = a_bits_key[r * 8 + 5];
      b_bits_key[r * 8 + 7] = a_bits_key[r * 8 + 6];
    }

    // ::6-9
    for (uint16_t i_bit = 0; i_bit < 8; ++i_bit) {
      out_key[8 * (c * 4) + i_bit] =
          bf256_add(bf256_add(bf256_add(b_bits_key[0 * 8 + i_bit], a_bits_key[3 * 8 + i_bit]),
                              bf256_add(a_bits_key[2 * 8 + i_bit], b_bits_key[1 * 8 + i_bit])),
                    a_bits_key[1 * 8 + i_bit]);
      out_key[8 * (c * 4 + 1) + i_bit] =
          bf256_add(bf256_add(bf256_add(b_bits_key[1 * 8 + i_bit], a_bits_key[0 * 8 + i_bit]),
                              bf256_add(a_bits_key[3 * 8 + i_bit], b_bits_key[2 * 8 + i_bit])),
                    a_bits_key[2 * 8 + i_bit]);
      out_key[8 * (c * 4 + 2) + i_bit] =
          bf256_add(bf256_add(bf256_add(b_bits_key[2 * 8 + i_bit], a_bits_key[1 * 8 + i_bit]),
                              bf256_add(a_bits_key[0 * 8 + i_bit], b_bits_key[3 * 8 + i_bit])),
                    a_bits_key[3 * 8 + i_bit]);
      out_key[8 * (c * 4 + 3) + i_bit] =
          bf256_add(bf256_add(bf256_add(b_bits_key[3 * 8 + i_bit], a_bits_key[2 * 8 + i_bit]),
                              bf256_add(a_bits_key[1 * 8 + i_bit], b_bits_key[0 * 8 + i_bit])),
                    a_bits_key[0 * 8 + i_bit]);
    }
  }
}

// CONSTANT TO VOLE
static void constant_to_vole_128_prover(bf128_t* tag, unsigned int n) {
  // the val stay the same as the val is a pub const!
  for (unsigned int i = 0; i < n; i++) {
    tag[i] = bf128_zero(); // for constant values the tag is zero
  }
}

static void constant_to_vole_192_prover(bf192_t* tag, unsigned int n) {
  // the val stay the same as the val is a pub const!
  for (unsigned int i = 0; i < n; i++) {
    tag[i] = bf192_zero(); // for constant values the tag is zero
  }
}

static void constant_to_vole_256_prover(bf256_t* tag, unsigned int n) {
  // the val stay the same as the val is a pub const!
  for (unsigned int i = 0; i < n; i++) {
    tag[i] = bf256_zero(); // for constant values the tag is zero
  }
}

static void constant_to_vole_128_verifier(bf128_t* key, const uint8_t* val, bf128_t delta,
                                          unsigned int n) {
  for (unsigned int i = 0; i < n; i++) {
    key[i] = bf128_mul_bit(delta, get_bit(val[i / 8], i % 8));
  }
}

static void constant_to_vole_192_verifier(bf192_t* key, const uint8_t* val, bf192_t delta,
                                          unsigned int n) {
  for (unsigned int i = 0; i < n; i++) {
    key[i] = bf192_mul_bit(delta, get_bit(val[i / 8], i % 8));
  }
}

static void constant_to_vole_256_verifier(bf256_t* key, const uint8_t* val, bf256_t delta,
                                          unsigned int n) {
  for (unsigned int i = 0; i < n; i++) {
    key[i] = bf256_mul_bit(delta, get_bit(val[i / 8], i % 8));
  }
}

// // INVERSE AFFINE
static void aes_128_inverse_affine_byte_prover(uint8_t* y_bits, bf128_t* y_bits_tag,
                                               const uint8_t* x_bits, const bf128_t* x_bits_tag) {
  uint8_t c = 0;

  for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
    if (bit_i == 0 || bit_i == 2) {
      c = 1;
    } else {
      c = 0;
    }
    y_bits[bit_i] =
        x_bits[(bit_i - 1 + 8) % 8] ^ x_bits[(bit_i - 3 + 8) % 8] ^ x_bits[(bit_i - 6 + 8) % 8] ^ c;
    y_bits_tag[bit_i] =
        bf128_add(bf128_add(x_bits_tag[(bit_i - 1 + 8) % 8], x_bits_tag[(bit_i - 3 + 8) % 8]),
                  x_bits_tag[(bit_i - 6 + 8) % 8]);
  }
}

static void aes_192_inverse_affine_byte_prover(uint8_t* y_bits, bf192_t* y_bits_tag,
                                               const uint8_t* x_bits, const bf192_t* x_bits_tag) {
  uint8_t c = 0;

  for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
    if (bit_i == 0 || bit_i == 2) {
      c = 1;
    } else {
      c = 0;
    }
    y_bits[bit_i] =
        x_bits[(bit_i - 1 + 8) % 8] ^ x_bits[(bit_i - 3 + 8) % 8] ^ x_bits[(bit_i - 6 + 8) % 8] ^ c;
    y_bits_tag[bit_i] =
        bf192_add(bf192_add(x_bits_tag[(bit_i - 1 + 8) % 8], x_bits_tag[(bit_i - 3 + 8) % 8]),
                  x_bits_tag[(bit_i - 6 + 8) % 8]);
  }
}

static void aes_256_inverse_affine_byte_prover(uint8_t* y_bits, bf256_t* y_bits_tag,
                                               const uint8_t* x_bits, const bf256_t* x_bits_tag) {
  uint8_t c = 0;

  for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
    if (bit_i == 0 || bit_i == 2) {
      c = 1;
    } else {
      c = 0;
    }
    y_bits[bit_i] =
        x_bits[(bit_i - 1 + 8) % 8] ^ x_bits[(bit_i - 3 + 8) % 8] ^ x_bits[(bit_i - 6 + 8) % 8] ^ c;
    y_bits_tag[bit_i] =
        bf256_add(bf256_add(x_bits_tag[(bit_i - 1 + 8) % 8], x_bits_tag[(bit_i - 3 + 8) % 8]),
                  x_bits_tag[(bit_i - 6 + 8) % 8]);
  }
}

static void aes_128_inverse_affine_prover(uint8_t* y, bf128_t* y_tag, const uint8_t* x,
                                          const bf128_t* x_tag, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    aes_128_inverse_affine_byte_prover(y + i * 8, y_tag + i * 8, x + i * 8, x_tag + i * 8);
  }
}

static void aes_192_inverse_affine_prover(uint8_t* y, bf192_t* y_tag, const uint8_t* x,
                                          const bf192_t* x_tag, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    aes_192_inverse_affine_byte_prover(y + i * 8, y_tag + i * 8, x + i * 8, x_tag + i * 8);
  }
}

static void aes_256_inverse_affine_prover(uint8_t* y, bf256_t* y_tag, const uint8_t* x,
                                          const bf256_t* x_tag, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    aes_256_inverse_affine_byte_prover(y + i * 8, y_tag + i * 8, x + i * 8, x_tag + i * 8);
  }
}

static void aes_128_inverse_affine_byte_verifier(bf128_t* y_bits_key, const bf128_t* x_bits_key,
                                                 bf128_t delta) {
  uint8_t c = 0;
  for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
    if (bit_i == 0 || bit_i == 2) {
      c = 1;
    } else {
      c = 0;
    }
    y_bits_key[bit_i] =
        bf128_add(bf128_add(x_bits_key[(bit_i - 1 + 8) % 8], x_bits_key[(bit_i - 3 + 8) % 8]),
                  bf128_add(x_bits_key[(bit_i - 6 + 8) % 8], bf128_mul_bit(delta, c)));
  }
}

static void aes_192_inverse_affine_byte_verifier(bf192_t* y_bits_key, const bf192_t* x_bits_key,
                                                 bf192_t delta) {
  uint8_t c = 0;
  for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
    if (bit_i == 0 || bit_i == 2) {
      c = 1;
    } else {
      c = 0;
    }
    y_bits_key[bit_i] =
        bf192_add(bf192_add(x_bits_key[(bit_i - 1 + 8) % 8], x_bits_key[(bit_i - 3 + 8) % 8]),
                  bf192_add(x_bits_key[(bit_i - 6 + 8) % 8], bf192_mul_bit(delta, c)));
  }
}

static void aes_256_inverse_affine_byte_verifier(bf256_t* y_bits_key, const bf256_t* x_bits_key,
                                                 bf256_t delta) {
  uint8_t c = 0;
  for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
    if (bit_i == 0 || bit_i == 2) {
      c = 1;
    } else {
      c = 0;
    }
    y_bits_key[bit_i] =
        bf256_add(bf256_add(x_bits_key[(bit_i - 1 + 8) % 8], x_bits_key[(bit_i - 3 + 8) % 8]),
                  bf256_add(x_bits_key[(bit_i - 6 + 8) % 8], bf256_mul_bit(delta, c)));
  }
}

static void aes_128_inverse_affine_verifier(bf128_t* y_key, const bf128_t* x_key, bf128_t delta,
                                            const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    aes_128_inverse_affine_byte_verifier(y_key + i * 8, x_key + i * 8, delta);
  }
}

static void aes_192_inverse_affine_verifier(bf192_t* y_key, const bf192_t* x_key, bf192_t delta,
                                            const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    aes_192_inverse_affine_byte_verifier(y_key + i * 8, x_key + i * 8, delta);
  }
}

static void aes_256_inverse_affine_verifier(bf256_t* y_key, const bf256_t* x_key, bf256_t delta,
                                            const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbytes = Nst * 4;

  for (unsigned int i = 0; i < Nstbytes; i++) {
    aes_256_inverse_affine_byte_verifier(y_key + i * 8, x_key + i * 8, delta);
  }
}

// EncSctrnts internal functions end!!

// COLOUM TO ROW MAJOR
static bf128_t* column_to_row_major_and_shrink_V_128(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + 2
  // \lambda rows and \lambda columns storing in row-major order
  bf128_t* new_v = BF128_ALLOC(ell + FAEST_128F_LAMBDA * 2);
  for (unsigned int row = 0; row != ell + FAEST_128F_LAMBDA * 2; ++row) {
    uint8_t new_row[BF128_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_128F_LAMBDA; ++column) {
      ptr_set_bit(new_row, column, ptr_get_bit(v[column], row));
    }
    new_v[row] = bf128_load(new_row);
  }
  return new_v;
}

static bf192_t* column_to_row_major_and_shrink_V_192(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf192_t* new_v = BF192_ALLOC(ell + FAEST_192F_LAMBDA * 2);
  for (unsigned int row = 0; row != ell + FAEST_192F_LAMBDA * 2; ++row) {
    uint8_t new_row[BF192_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_192F_LAMBDA; ++column) {
      ptr_set_bit(new_row, column, ptr_get_bit(v[column], row));
    }
    new_v[row] = bf192_load(new_row);
  }

  return new_v;
}

static bf256_t* column_to_row_major_and_shrink_V_256(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf256_t* new_v = BF256_ALLOC(ell + FAEST_256F_LAMBDA * 2);
  for (unsigned int row = 0; row != ell + FAEST_256F_LAMBDA * 2; ++row) {
    uint8_t new_row[BF256_NUM_BYTES] = {0};
    for (unsigned int column = 0; column != FAEST_256F_LAMBDA; ++column) {
      ptr_set_bit(new_row, column, ptr_get_bit(v[column], row));
    }
    new_v[row] = bf256_load(new_row);
  }

  return new_v;
}

// // KEY EXP FWD / BKWD
static void aes_128_keyexp_backward_prover(uint8_t* y, bf128_t* y_tag, const uint8_t* x,
                                           const bf128_t* x_tag, const uint8_t* key,
                                           const bf128_t* key_tag, const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_128_LAMBDA;

  // ::2
  uint8_t x_tilde[8];
  bf128_t x_tilde_tag[8];
  // ::3
  unsigned int iwd = 0;
  // ::4
  bool rmvRcon = true;
  // ::5-6
  for (unsigned int j = 0; j < Ske; j++) {
    // ::7-10
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {

      x_tilde[bit_i] = x[j * 8 + bit_i] ^ key[iwd + (j % 4) * 8 + bit_i]; // for the witness
      x_tilde_tag[bit_i] =
          bf128_add(x_tag[j * 8 + bit_i],
                    key_tag[iwd + (j % 4) * 8 + bit_i]); // for the tags of each witness bit

      if (rmvRcon == true && j % 4 == 0) {
        // adding round constant to the witness
        x_tilde[bit_i] = x_tilde[bit_i] ^ get_bit(Rcon[j / 4], bit_i);
      }
    }

    // ::11
    aes_128_inverse_affine_byte_prover(y + 8 * j, y_tag + 8 * j, x_tilde,
                                       x_tilde_tag); // working in bit per uint8

    // ::12-16 lines only relavant for aes-128
    if (j % 4 == 3) {
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

static void aes_192_keyexp_backward_prover(uint8_t* y, bf192_t* y_tag, const uint8_t* x,
                                           const bf192_t* x_tag, const uint8_t* key,
                                           const bf192_t* key_tag, const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_192_LAMBDA;

  // ::2
  uint8_t x_tilde[8];
  bf192_t x_tilde_tag[8];
  // ::3
  unsigned int iwd = 0;
  // ::4
  bool rmvRcon = true;
  // ::5-6
  for (unsigned int j = 0; j < Ske; j++) {
    // ::7-10
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {

      x_tilde[bit_i] = x[j * 8 + bit_i] ^ key[iwd + (j % 4) * 8 + bit_i]; // for the witness
      x_tilde_tag[bit_i] =
          bf192_add(x_tag[j * 8 + bit_i],
                    key_tag[iwd + (j % 4) * 8 + bit_i]); // for the tags of each witness bit

      if (rmvRcon == true && j % 4 == 0) {
        // adding round constant to the witness
        x_tilde[bit_i] = x_tilde[bit_i] ^ get_bit(Rcon[j / 4], bit_i);
      }
    }

    // ::11
    aes_192_inverse_affine_byte_prover(y + 8 * j, y_tag + 8 * j, x_tilde,
                                       x_tilde_tag); // working in bit per uint8

    // ::12-16
    if (j % 4 == 3) {
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

static void aes_256_keyexp_backward_prover(uint8_t* y, bf256_t* y_tag, const uint8_t* x,
                                           const bf256_t* x_tag, const uint8_t* key,
                                           const bf256_t* key_tag, const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_256_LAMBDA;

  // ::2
  uint8_t x_tilde[8];
  bf256_t x_tilde_tag[8];
  // ::3
  unsigned int iwd = 0;
  // ::4
  bool rmvRcon = true;
  // ::5-6
  for (unsigned int j = 0; j < Ske; j++) {
    // ::7-10
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {

      x_tilde[bit_i] = x[j * 8 + bit_i] ^ key[iwd + (j % 4) * 8 + bit_i]; // for the witness
      x_tilde_tag[bit_i] =
          bf256_add(x_tag[j * 8 + bit_i],
                    key_tag[iwd + (j % 4) * 8 + bit_i]); // for the tags of each witness bit

      if (rmvRcon == true && j % 4 == 0) {
        // adding round constant to the witness
        x_tilde[bit_i] = x_tilde[bit_i] ^ get_bit(Rcon[j / 8], bit_i);
      }
    }

    // ::11
    aes_256_inverse_affine_byte_prover(y + 8 * j, y_tag + 8 * j, x_tilde,
                                       x_tilde_tag); // working in bit per uint8

    // ::12-16
    if (j % 4 == 3) {
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

static void aes_128_keyexp_backward_verifier(bf128_t* y_key, const bf128_t* x_key,
                                             const bf128_t* key_key, bf128_t delta,
                                             const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_128_LAMBDA;

  // ::2
  bf128_t x_tilde_key[8];
  // ::3
  unsigned int iwd = 0;
  // ::4
  bool rmvRcon = true;
  // ::5-6
  for (unsigned int j = 0; j < Ske; j++) {
    // ::7
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_tilde_key[bit_i] =
          bf128_add(x_key[j * 8 + bit_i],
                    key_key[iwd + (j % 4) * 8 + bit_i]); // for the tags of each witness bit
      // ::8-10
      if (rmvRcon == true && j % 4 == 0) {
        bf128_t rcon_key;
        const uint8_t c = (Rcon[j / 4] >> bit_i) & 1;
        constant_to_vole_128_verifier(&rcon_key, &c, delta, 1);
        x_tilde_key[bit_i] = bf128_add(x_tilde_key[bit_i], rcon_key);
      }
    }
    // ::11
    aes_128_inverse_affine_byte_verifier(y_key + 8 * j, x_tilde_key, delta);

    // ::12-16 lines only relavant for aes-128
    if (j % 4 == 3) {
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

static void aes_192_keyexp_backward_verifier(bf192_t* y_key, const bf192_t* x_key, bf192_t* key_key,
                                             bf192_t delta, const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_192_LAMBDA;

  // ::2
  bf192_t x_tilde_key[8];
  // ::3
  unsigned int iwd = 0;
  // ::4
  bool rmvRcon = true;
  // ::5-6
  for (unsigned int j = 0; j < Ske; j++) {
    // ::7
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_tilde_key[bit_i] =
          bf192_add(x_key[j * 8 + bit_i],
                    key_key[iwd + (j % 4) * 8 + bit_i]); // for the tags of each witness bit
      // ::8-10
      if (rmvRcon == true && j % 4 == 0) {
        bf192_t rcon_key;
        const uint8_t c = (Rcon[j / 4] >> bit_i) & 1;
        constant_to_vole_192_verifier(&rcon_key, &c, delta, 1);
        x_tilde_key[bit_i] = bf192_add(x_tilde_key[bit_i], rcon_key);
      }
    }
    // ::11
    aes_192_inverse_affine_byte_verifier(y_key + 8 * j, x_tilde_key, delta);

    // ::12-16
    if (j % 4 == 3) {
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

static void aes_256_keyexp_backward_verifier(bf256_t* y_key, const bf256_t* x_key, bf256_t* key_key,
                                             bf256_t delta, const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_256_LAMBDA;

  // ::2
  bf256_t x_tilde_key[8];
  // ::3
  unsigned int iwd = 0;
  // ::4
  bool rmvRcon = true;
  // ::5-6
  for (unsigned int j = 0; j < Ske; j++) {
    // ::7
    for (unsigned int bit_i = 0; bit_i < 8; bit_i++) {
      x_tilde_key[bit_i] =
          bf256_add(x_key[j * 8 + bit_i],
                    key_key[iwd + (j % 4) * 8 + bit_i]); // for the tags of each witness bit
      // ::8-10
      if (rmvRcon == true && j % 4 == 0) {
        bf256_t rcon_key;
        const uint8_t c = (Rcon[j / 8] >> bit_i) & 1;
        constant_to_vole_256_verifier(&rcon_key, &c, delta, 1);
        x_tilde_key[bit_i] = bf256_add(x_tilde_key[bit_i], rcon_key);
      }
    }
    // ::11
    aes_256_inverse_affine_byte_verifier(y_key + 8 * j, x_tilde_key, delta);

    // ::12-16
    if (j % 4 == 3) {
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

static void aes_128_keyexp_forward_prover(uint8_t* y, bf128_t* y_tag, const uint8_t* w,
                                          const bf128_t* w_tag, const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_128_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  const unsigned int R      = params->R;

  // ::1-2
  for (unsigned int i = 0; i < lambda; i++) {
    y[i]     = w[i];
    y_tag[i] = w_tag[i];
  }
  // ::3
  unsigned int i_wd = lambda;
  // ::4-10
  for (unsigned int j = Nk; j < 4 * (R + 1); j++) {
    // ::5
    if ((j % Nk == 0) || ((Nk > 6) && (j % Nk == 4))) {
      // ::6
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y[32 * j + word_idx]     = w[i_wd + word_idx];     // storing bit by bit
        y_tag[32 * j + word_idx] = w_tag[i_wd + word_idx]; // storing tags
      }
      // ::7
      i_wd += 32;
      // ::8
    } else {
      // ::9-10
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y[32 * j + word_idx] = y[32 * (j - Nk) + word_idx] ^ y[32 * (j - 1) + word_idx];
        y_tag[32 * j + word_idx] =
            bf128_add(y_tag[32 * (j - Nk) + word_idx], y_tag[32 * (j - 1) + word_idx]);
      }
    }
  }
}

static void aes_192_keyexp_forward_prover(uint8_t* y, bf192_t* y_tag, const uint8_t* w,
                                          const bf192_t* w_tag, const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_192_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  const unsigned int R      = params->R;

  // ::1-2
  for (unsigned int i = 0; i < lambda; i++) {
    y[i]     = w[i];
    y_tag[i] = w_tag[i];
  }
  // ::3
  unsigned int i_wd = lambda;
  // ::4-10
  for (unsigned int j = Nk; j < 4 * (R + 1); j++) {
    // ::5
    if ((j % Nk == 0) || ((Nk > 6) && (j % Nk == 4))) {
      // ::6
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y[32 * j + word_idx]     = w[i_wd + word_idx];     // storing bit by bit
        y_tag[32 * j + word_idx] = w_tag[i_wd + word_idx]; // storing tags
      }
      // ::7
      i_wd += 32;
      // ::8
    } else {
      // ::9-10
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y[32 * j + word_idx] = y[32 * (j - Nk) + word_idx] ^ y[32 * (j - 1) + word_idx];
        y_tag[32 * j + word_idx] =
            bf192_add(y_tag[32 * (j - Nk) + word_idx], y_tag[32 * (j - 1) + word_idx]);
      }
    }
  }
}

static void aes_256_keyexp_forward_prover(uint8_t* y, bf256_t* y_tag, const uint8_t* w,
                                          const bf256_t* w_tag, const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_256_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  const unsigned int R      = params->R;

  // ::1-2
  for (unsigned int i = 0; i < lambda; i++) {
    y[i]     = w[i];
    y_tag[i] = w_tag[i];
  }
  // ::3
  unsigned int i_wd = lambda;
  // ::4-10
  for (unsigned int j = Nk; j < 4 * (R + 1); j++) {
    // ::5
    if ((j % Nk == 0) || ((Nk > 6) && (j % Nk == 4))) {
      // ::6
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y[32 * j + word_idx]     = w[i_wd + word_idx];     // storing bit by bit
        y_tag[32 * j + word_idx] = w_tag[i_wd + word_idx]; // storing tags
      }
      // ::7
      i_wd += 32;
      // ::8
    } else {
      // ::9-10
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y[32 * j + word_idx] = y[32 * (j - Nk) + word_idx] ^ y[32 * (j - 1) + word_idx];
        y_tag[32 * j + word_idx] =
            bf256_add(y_tag[32 * (j - Nk) + word_idx], y_tag[32 * (j - 1) + word_idx]);
      }
    }
  }
}

static void aes_128_keyexp_forward_verifier(bf128_t* y_key, const bf128_t* w_key,
                                            const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_128_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  const unsigned int R      = params->R;

  // ::1-2
  for (unsigned int i = 0; i < lambda; i++) {
    y_key[i] = w_key[i];
  }
  // ::3
  unsigned int i_wd = lambda;
  // ::4-10
  for (unsigned int j = Nk; j < 4 * (R + 1); j++) {
    // ::5
    if ((j % Nk == 0) || ((Nk > 6) && (j % Nk == 4))) {
      // ::6
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y_key[32 * j + word_idx] = w_key[i_wd + word_idx];
      }
      // ::7
      i_wd += 32; // 32 bits -> 4 words
      // ::8
    } else {
      // ::9-10
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y_key[32 * j + word_idx] =
            bf128_add(y_key[32 * (j - Nk) + word_idx], y_key[32 * (j - 1) + word_idx]);
      }
    }
  }
}

static void aes_192_keyexp_forward_verifier(bf192_t* y_key, const bf192_t* w_key,
                                            const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_192_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  const unsigned int R      = params->R;

  // ::1-2
  for (unsigned int i = 0; i < lambda; i++) {
    y_key[i] = w_key[i];
  }
  // ::3
  unsigned int i_wd = lambda;
  // ::4-10
  for (unsigned int j = Nk; j < 4 * (R + 1); j++) {
    // ::5
    if ((j % Nk == 0) || ((Nk > 6) && (j % Nk == 4))) {
      // ::6
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y_key[32 * j + word_idx] = w_key[i_wd + word_idx];
      }
      // ::7
      i_wd += 32; // 32 bits -> 4 words
      // ::8
    } else {
      // ::9-10
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y_key[32 * j + word_idx] =
            bf192_add(y_key[32 * (j - Nk) + word_idx], y_key[32 * (j - 1) + word_idx]);
      }
    }
  }
}

static void aes_256_keyexp_forward_verifier(bf256_t* y_key, const bf256_t* w_key,
                                            const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_256_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  const unsigned int R      = params->R;

  // ::1-2
  for (unsigned int i = 0; i < lambda; i++) {
    y_key[i] = w_key[i];
  }
  // ::3
  unsigned int i_wd = lambda;
  // ::4-10
  for (unsigned int j = Nk; j < 4 * (R + 1); j++) {
    // ::5
    if ((j % Nk == 0) || ((Nk > 6) && (j % Nk == 4))) {
      // ::6
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y_key[32 * j + word_idx] = w_key[i_wd + word_idx];
      }
      // ::7
      i_wd += 32; // 32 bits -> 4 words
      // ::8
    } else {
      // ::9-10
      for (unsigned int word_idx = 0; word_idx < 32; word_idx++) {
        y_key[32 * j + word_idx] =
            bf256_add(y_key[32 * (j - Nk) + word_idx], y_key[32 * (j - 1) + word_idx]);
      }
    }
  }
}

// // KEY EXP CSTRNTS
static void aes_128_expkey_constraints_prover(zk_hash_128_3_ctx* hasher, uint8_t* k, bf128_t* k_tag,
                                              const uint8_t* w, const bf128_t* w_tag,
                                              const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_128_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  unsigned int r_prime;

  bool do_rot_word = true;

  // ::1
  aes_128_keyexp_forward_prover(k, k_tag, w, w_tag, params);
  // ::2
  uint8_t* w_flat     = malloc(8 * Ske * sizeof(uint8_t));
  bf128_t* w_flat_tag = BF128_ALLOC(8 * Ske);
  aes_128_keyexp_backward_prover(w_flat, w_flat_tag, w + lambda, w_tag + lambda, k, k_tag, params);

  // ::3-5
  unsigned int iwd = 32 * (Nk - 1);
  // ::6 Used only on AES-256
  // ::7
  for (unsigned int j = 0; j < Ske / 4; j++) {
    // ::8
    bf128_t k_hat[4];    // expnaded key witness
    bf128_t w_hat[4];    // inverse output
    bf128_t k_hat_sq[4]; // expanded key witness sq
    bf128_t w_hat_sq[4]; // inverse output sq

    bf128_t k_hat_tag[4];    // expanded key witness tag
    bf128_t w_hat_tag[4];    // inverse output tag
    bf128_t k_hat_tag_sq[4]; // expanded key tag sq
    bf128_t w_hat_tag_sq[4]; // inverser output tag sq

    // ::9
    for (unsigned int r = 0; r < 4; r++) {
      // ::10
      r_prime = r;
      // ::11
      if (do_rot_word) {
        r_prime = (r + 3) % 4;
      }
      // ::12-15
      k_hat[r_prime] = bf128_byte_combine_bits_array(&k[(iwd + 8 * r)]); // lifted key witness
      k_hat_sq[r_prime] =
          bf128_byte_combine_bits_array_sq(&k[(iwd + 8 * r)]); // lifted key witness sq

      w_hat[r]    = bf128_byte_combine_bits_array(&w_flat[(32 * j + 8 * r)]);    // lifted output
      w_hat_sq[r] = bf128_byte_combine_bits_array_sq(&w_flat[(32 * j + 8 * r)]); // lifted output sq

      // done by both prover and verifier
      k_hat_tag[r_prime]    = bf128_byte_combine(k_tag + (iwd + 8 * r));    // lifted key tag
      k_hat_tag_sq[r_prime] = bf128_byte_combine_sq(k_tag + (iwd + 8 * r)); // lifted key tag sq

      w_hat_tag[r] = bf128_byte_combine(w_flat_tag + ((32 * j + 8 * r))); // lifted output tag
      w_hat_tag_sq[r] =
          bf128_byte_combine_sq(w_flat_tag + (32 * j + 8 * r)); // lifted output tag sq
    }
    // ::16 used only for AES-256
    if (lambda == 256) {
      do_rot_word = !do_rot_word;
    }
    // ::17
    for (unsigned int r = 0; r < 4; r++) {
      // ::18-19
      zk_hash_128_3_raise_and_update(hasher, bf128_mul(k_hat_tag_sq[r], w_hat_tag[r]),
                                     bf128_add(bf128_add(bf128_mul(k_hat_sq[r], w_hat_tag[r]),
                                                         bf128_mul(k_hat_tag_sq[r], w_hat[r])),
                                               k_hat_tag[r]));
      zk_hash_128_3_raise_and_update(hasher, bf128_mul(k_hat_tag[r], w_hat_tag_sq[r]),
                                     bf128_add(bf128_add(bf128_mul(k_hat[r], w_hat_tag_sq[r]),
                                                         bf128_mul(k_hat_tag[r], w_hat_sq[r])),
                                               w_hat_tag[r]));
    }
    if (lambda == 192) {
      iwd += 192;
    } else {
      iwd += 128;
    }
  }
  free(w_flat);
  faest_aligned_free(w_flat_tag);
}

static void aes_192_expkey_constraints_prover(bf192_t* z_deg0, bf192_t* z_deg1, uint8_t* k,
                                              bf192_t* k_tag, const uint8_t* w,
                                              const bf192_t* w_tag,
                                              const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_192_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  unsigned int r_prime;

  bool do_rot_word = true;

  // ::1
  aes_192_keyexp_forward_prover(k, k_tag, w, w_tag, params);
  // ::2
  uint8_t* w_flat     = malloc(8 * Ske * sizeof(uint8_t));
  bf192_t* w_flat_tag = BF192_ALLOC(8 * Ske);
  aes_192_keyexp_backward_prover(w_flat, w_flat_tag, w + lambda, w_tag + lambda, k, k_tag, params);

  // ::3-5
  unsigned int iwd = 32 * (Nk - 1);
  // ::6 Used only on AES-256
  // ::7
  for (unsigned int j = 0; j < Ske / 4; j++) {
    // ::8
    bf192_t k_hat[4];    // expnaded key witness
    bf192_t w_hat[4];    // inverse output
    bf192_t k_hat_sq[4]; // expanded key witness sq
    bf192_t w_hat_sq[4]; // inverse output sq

    bf192_t k_hat_tag[4];    // expanded key witness tag
    bf192_t w_hat_tag[4];    // inverse output tag
    bf192_t k_hat_tag_sq[4]; // expanded key tag sq
    bf192_t w_hat_tag_sq[4]; // inverser output tag sq

    // ::9
    for (unsigned int r = 0; r < 4; r++) {
      // ::10
      r_prime = r;
      // ::11
      if (do_rot_word) {
        r_prime = (r + 3) % 4;
      }
      // ::12-15
      k_hat[r_prime] = bf192_byte_combine_bits_array(&k[(iwd + 8 * r)]); // lifted key witness
      k_hat_sq[r_prime] =
          bf192_byte_combine_bits_array_sq(&k[(iwd + 8 * r)]); // lifted key witness sq

      w_hat[r]    = bf192_byte_combine_bits_array(&w_flat[(32 * j + 8 * r)]);    // lifted output
      w_hat_sq[r] = bf192_byte_combine_bits_array_sq(&w_flat[(32 * j + 8 * r)]); // lifted output sq

      // done by both prover and verifier
      k_hat_tag[r_prime]    = bf192_byte_combine(k_tag + (iwd + 8 * r));    // lifted key tag
      k_hat_tag_sq[r_prime] = bf192_byte_combine_sq(k_tag + (iwd + 8 * r)); // lifted key tag sq

      w_hat_tag[r] = bf192_byte_combine(w_flat_tag + ((32 * j + 8 * r))); // lifted output tag
      w_hat_tag_sq[r] =
          bf192_byte_combine_sq(w_flat_tag + (32 * j + 8 * r)); // lifted output tag sq
    }
    // ::16 used only for AES-256
    if (lambda == 256) {
      do_rot_word = !do_rot_word;
    }
    // ::17
    for (unsigned int r = 0; r < 4; r++) {
      // ::18-19
      z_deg1[8 * j + 2 * r] = bf192_add(
          bf192_add(bf192_mul(k_hat_sq[r], w_hat_tag[r]), bf192_mul(k_hat_tag_sq[r], w_hat[r])),
          k_hat_tag[r]);

      z_deg1[8 * j + 2 * r + 1] = bf192_add(
          bf192_add(bf192_mul(k_hat[r], w_hat_tag_sq[r]), bf192_mul(k_hat_tag[r], w_hat_sq[r])),
          w_hat_tag[r]);

      z_deg0[8 * j + 2 * r]     = bf192_mul(k_hat_tag_sq[r], w_hat_tag[r]);
      z_deg0[8 * j + 2 * r + 1] = bf192_mul(k_hat_tag[r], w_hat_tag_sq[r]);

      // z_deg1[8*j + 2*r + 1] = bf192_add(bf192_mul(k_hat[r], w_hat_sq[r]), w_hat[r]);
      // z_deg0[8*j + 2*r] = bf192_add(bf192_mul(k_hat_tag_sq[r], w_hat_tag[r]), k_hat_tag[r]);
      // z_deg0[8*j + 2*r + 1] = bf192_add(bf192_mul(k_hat_tag[r], w_hat_tag_sq[r]), k_hat_tag[r]);
    }
    if (lambda == 192) {
      iwd += 192;
    } else {
      iwd += 128;
    }
  }
  free(w_flat);
  faest_aligned_free(w_flat_tag);
}

static void aes_256_expkey_constraints_prover(bf256_t* z_deg0, bf256_t* z_deg1, uint8_t* k,
                                              bf256_t* k_tag, const uint8_t* w,
                                              const bf256_t* w_tag,
                                              const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_256_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  unsigned int r_prime;

  bool do_rot_word = true;

  // ::1
  aes_256_keyexp_forward_prover(k, k_tag, w, w_tag, params);
  // ::2
  uint8_t* w_flat     = malloc(8 * Ske * sizeof(uint8_t));
  bf256_t* w_flat_tag = BF256_ALLOC(8 * Ske);
  aes_256_keyexp_backward_prover(w_flat, w_flat_tag, w + lambda, w_tag + lambda, k, k_tag, params);

  // ::3-5
  unsigned int iwd = 32 * (Nk - 1);
  // ::6 Used only on AES-256
  // ::7
  for (unsigned int j = 0; j < Ske / 4; j++) {
    // ::8
    bf256_t k_hat[4];    // expnaded key witness
    bf256_t w_hat[4];    // inverse output
    bf256_t k_hat_sq[4]; // expanded key witness sq
    bf256_t w_hat_sq[4]; // inverse output sq

    bf256_t k_hat_tag[4];    // expanded key witness tag
    bf256_t w_hat_tag[4];    // inverse output tag
    bf256_t k_hat_tag_sq[4]; // expanded key tag sq
    bf256_t w_hat_tag_sq[4]; // inverser output tag sq

    // ::9
    for (unsigned int r = 0; r < 4; r++) {
      // ::10
      r_prime = r;
      // ::11
      if (do_rot_word) {
        r_prime = (r + 3) % 4;
      }
      // ::12-15
      k_hat[r_prime] = bf256_byte_combine_bits_array(&k[(iwd + 8 * r)]); // lifted key witness
      k_hat_sq[r_prime] =
          bf256_byte_combine_bits_array_sq(&k[(iwd + 8 * r)]); // lifted key witness sq

      w_hat[r]    = bf256_byte_combine_bits_array(&w_flat[(32 * j + 8 * r)]);    // lifted output
      w_hat_sq[r] = bf256_byte_combine_bits_array_sq(&w_flat[(32 * j + 8 * r)]); // lifted output sq

      // done by both prover and verifier
      k_hat_tag[r_prime]    = bf256_byte_combine(k_tag + (iwd + 8 * r));    // lifted key tag
      k_hat_tag_sq[r_prime] = bf256_byte_combine_sq(k_tag + (iwd + 8 * r)); // lifted key tag sq

      w_hat_tag[r] = bf256_byte_combine(w_flat_tag + ((32 * j + 8 * r))); // lifted output tag
      w_hat_tag_sq[r] =
          bf256_byte_combine_sq(w_flat_tag + (32 * j + 8 * r)); // lifted output tag sq
    }
    // ::16 used only for AES-256
    if (lambda == 256) {
      do_rot_word = !do_rot_word;
    }
    // ::17
    for (unsigned int r = 0; r < 4; r++) {
      // ::18-19
      z_deg1[8 * j + 2 * r] = bf256_add(
          bf256_add(bf256_mul(k_hat_sq[r], w_hat_tag[r]), bf256_mul(k_hat_tag_sq[r], w_hat[r])),
          k_hat_tag[r]);

      z_deg1[8 * j + 2 * r + 1] = bf256_add(
          bf256_add(bf256_mul(k_hat[r], w_hat_tag_sq[r]), bf256_mul(k_hat_tag[r], w_hat_sq[r])),
          w_hat_tag[r]);

      z_deg0[8 * j + 2 * r]     = bf256_mul(k_hat_tag_sq[r], w_hat_tag[r]);
      z_deg0[8 * j + 2 * r + 1] = bf256_mul(k_hat_tag[r], w_hat_tag_sq[r]);

      // z_deg1[8*j + 2*r + 1] = bf256_add(bf256_mul(k_hat[r], w_hat_sq[r]), w_hat[r]);
      // z_deg0[8*j + 2*r] = bf256_add(bf256_mul(k_hat_tag_sq[r], w_hat_tag[r]), k_hat_tag[r]);
      // z_deg0[8*j + 2*r + 1] = bf256_add(bf256_mul(k_hat_tag[r], w_hat_tag_sq[r]), k_hat_tag[r]);
    }
    if (lambda == 192) {
      iwd += 192;
    } else {
      iwd += 128;
    }
  }
  free(w_flat);
  faest_aligned_free(w_flat_tag);
}

static void aes_128_expkey_constraints_verifier(zk_hash_128_ctx* hasher, bf128_t* k_key,
                                                const bf128_t* w_key, bf128_t delta,
                                                const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_128_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  unsigned int r_prime;

  bool do_rot_word = true;

  // ::1
  aes_128_keyexp_forward_verifier(k_key, w_key, params);
  // ::2
  bf128_t* w_flat_key = BF128_ALLOC(8 * Ske);
  aes_128_keyexp_backward_verifier(w_flat_key, w_key + lambda, k_key, delta, params);

  // ::3-5
  unsigned int iwd = 32 * (Nk - 1); // as 1 unit8 has 8 bits
  // ::6 Used only on AES-256
  // ::7
  for (unsigned int j = 0; j < Ske / 4; j++) {
    // ::8
    bf128_t k_hat_key[4];    // expanded key witness tag
    bf128_t w_hat_key[4];    // inverse output tag
    bf128_t k_hat_key_sq[4]; // expanded key tag sq
    bf128_t w_hat_key_sq[4]; // inverser output tag sq

    // ::9
    for (unsigned int r = 0; r < 4; r++) {
      // ::10
      r_prime = r;
      // ::11
      if (do_rot_word) {
        r_prime = (r + 3) % 4;
      }
      // ::12-15
      k_hat_key[r_prime]    = bf128_byte_combine(k_key + (iwd + 8 * r));    // lifted key tag
      k_hat_key_sq[r_prime] = bf128_byte_combine_sq(k_key + (iwd + 8 * r)); // lifted key tag sq

      w_hat_key[r] = bf128_byte_combine(w_flat_key + ((32 * j + 8 * r))); // lifted output tag
      w_hat_key_sq[r] =
          bf128_byte_combine_sq(w_flat_key + (32 * j + 8 * r)); // lifted output tag sq
    }
    // ::16 used only for AES-256
    if (lambda == 256) {
      do_rot_word = !do_rot_word;
    }
    // ::17-20
    for (unsigned int r = 0; r < 4; r++) {
      // also raise degree
      zk_hash_128_update(hasher,
                         bf128_mul(delta, bf128_add(bf128_mul(k_hat_key_sq[r], w_hat_key[r]),
                                                    bf128_mul(delta, k_hat_key[r]))));
      zk_hash_128_update(hasher,
                         bf128_mul(delta, bf128_add(bf128_mul(k_hat_key[r], w_hat_key_sq[r]),
                                                    bf128_mul(delta, w_hat_key[r]))));
    }
    if (lambda == 192) {
      iwd += 192;
    } else {
      iwd += 128;
    }
  }
  faest_aligned_free(w_flat_key);
}

static void aes_192_expkey_constraints_verifier(bf192_t* z_deg1, bf192_t* k_key,
                                                const bf192_t* w_key, bf192_t delta,
                                                const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_192_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  unsigned int r_prime;

  bool do_rot_word = true;

  // ::1
  aes_192_keyexp_forward_verifier(k_key, w_key, params);
  // ::2
  bf192_t* w_flat_key = BF192_ALLOC(8 * Ske);
  aes_192_keyexp_backward_verifier(w_flat_key, w_key + lambda, k_key, delta, params);

  // ::3-5
  unsigned int iwd = 32 * (Nk - 1); // as 1 unit8 has 8 bits
  // ::6 Used only on AES-256
  // ::7
  for (unsigned int j = 0; j < Ske / 4; j++) {
    // ::8
    bf192_t k_hat_key[4];    // expanded key witness tag
    bf192_t w_hat_key[4];    // inverse output tag
    bf192_t k_hat_key_sq[4]; // expanded key tag sq
    bf192_t w_hat_key_sq[4]; // inverser output tag sq

    // ::9
    for (unsigned int r = 0; r < 4; r++) {
      // ::10
      r_prime = r;
      // ::11
      if (do_rot_word) {
        r_prime = (r + 3) % 4;
      }
      // ::12-15
      k_hat_key[r_prime]    = bf192_byte_combine(k_key + (iwd + 8 * r));    // lifted key tag
      k_hat_key_sq[r_prime] = bf192_byte_combine_sq(k_key + (iwd + 8 * r)); // lifted key tag sq

      w_hat_key[r] = bf192_byte_combine(w_flat_key + ((32 * j + 8 * r))); // lifted output tag
      w_hat_key_sq[r] =
          bf192_byte_combine_sq(w_flat_key + (32 * j + 8 * r)); // lifted output tag sq
    }
    // ::16 used only for AES-256
    if (lambda == 256) {
      do_rot_word = !do_rot_word;
    }
    // ::17-20
    for (unsigned int r = 0; r < 4; r++) {
      z_deg1[8 * j + 2 * r] =
          bf192_add(bf192_mul(k_hat_key_sq[r], w_hat_key[r]), bf192_mul(delta, k_hat_key[r]));
      z_deg1[8 * j + 2 * r + 1] =
          bf192_add(bf192_mul(k_hat_key[r], w_hat_key_sq[r]), bf192_mul(delta, w_hat_key[r]));
    }
    if (lambda == 192) {
      iwd += 192;
    } else {
      iwd += 128;
    }
  }
  faest_aligned_free(w_flat_key);
}

static void aes_256_expkey_constraints_verifier(bf256_t* z_deg1, bf256_t* k_key,
                                                const bf256_t* w_key, bf256_t delta,
                                                const faest_paramset_t* params) {
  const unsigned int Ske    = params->Ske;
  const unsigned int lambda = FAEST_256_LAMBDA;
  const unsigned int Nk     = lambda / 32;
  unsigned int r_prime;

  bool do_rot_word = true;

  // ::1
  aes_256_keyexp_forward_verifier(k_key, w_key, params);
  // ::2
  bf256_t* w_flat_key = BF256_ALLOC(8 * Ske);
  aes_256_keyexp_backward_verifier(w_flat_key, w_key + lambda, k_key, delta, params);

  // ::3-5
  unsigned int iwd = 32 * (Nk - 1); // as 1 unit8 has 8 bits
  // ::6 Used only on AES-256
  // ::7
  for (unsigned int j = 0; j < Ske / 4; j++) {
    // ::8
    bf256_t k_hat_key[4];    // expanded key witness tag
    bf256_t w_hat_key[4];    // inverse output tag
    bf256_t k_hat_key_sq[4]; // expanded key tag sq
    bf256_t w_hat_key_sq[4]; // inverser output tag sq

    // ::9
    for (unsigned int r = 0; r < 4; r++) {
      // ::10
      r_prime = r;
      // ::11
      if (do_rot_word) {
        r_prime = (r + 3) % 4;
      }
      // ::12-15
      k_hat_key[r_prime]    = bf256_byte_combine(k_key + (iwd + 8 * r));    // lifted key tag
      k_hat_key_sq[r_prime] = bf256_byte_combine_sq(k_key + (iwd + 8 * r)); // lifted key tag sq

      w_hat_key[r] = bf256_byte_combine(w_flat_key + ((32 * j + 8 * r))); // lifted output tag
      w_hat_key_sq[r] =
          bf256_byte_combine_sq(w_flat_key + (32 * j + 8 * r)); // lifted output tag sq
    }
    // ::16 used only for AES-256
    if (lambda == 256) {
      do_rot_word = !do_rot_word;
    }
    // ::17-20
    for (unsigned int r = 0; r < 4; r++) {
      z_deg1[8 * j + 2 * r] =
          bf256_add(bf256_mul(k_hat_key_sq[r], w_hat_key[r]), bf256_mul(delta, k_hat_key[r]));
      z_deg1[8 * j + 2 * r + 1] =
          bf256_add(bf256_mul(k_hat_key[r], w_hat_key_sq[r]), bf256_mul(delta, w_hat_key[r]));
    }
    if (lambda == 192) {
      iwd += 192;
    } else {
      iwd += 128;
    }
  }
  faest_aligned_free(w_flat_key);
}

// // ENC CSTRNTS
static void aes_128_enc_constraints_prover(zk_hash_128_3_ctx* hasher, const uint8_t* owf_in,
                                           const bf128_t* owf_in_tag, const uint8_t* owf_out,
                                           const bf128_t* owf_out_tag, const uint8_t* w,
                                           const bf128_t* w_tag, const uint8_t* k,
                                           const bf128_t* k_tag, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbits  = 32 * Nst;
  const unsigned int R        = params->R;
  const unsigned int Nstbytes = Nstbits / 8;

  /// ::1 AddRoundKey
  uint8_t* state_bits     = malloc(Nstbits * sizeof(uint8_t));
  bf128_t* state_bits_tag = BF128_ALLOC(Nstbits);

  aes_128_add_round_key_prover(state_bits, state_bits_tag, owf_in, owf_in_tag, k, k_tag, params);

  // for conjugates of state and s-box outputs
  bf128_t* state_conj     = BF128_ALLOC(8 * Nstbytes);
  bf128_t* state_conj_tag = BF128_ALLOC(8 * Nstbytes);
  bf128_t* st_dash_deg2   = BF128_ALLOC(8 * Nstbytes);
  bf128_t* st_dash_deg1   = BF128_ALLOC(8 * Nstbytes);
  bf128_t* st_dash_deg0   = BF128_ALLOC(8 * Nstbytes);

  bf128_t* state_bytewise_deg2    = BF128_ALLOC(Nstbytes);
  bf128_t* state_bytewise_deg1    = BF128_ALLOC(Nstbytes);
  bf128_t* state_bytewise_deg0    = BF128_ALLOC(Nstbytes);
  bf128_t* state_sq_bytewise_deg2 = BF128_ALLOC(Nstbytes);
  bf128_t* state_sq_bytewise_deg1 = BF128_ALLOC(Nstbytes);
  bf128_t* state_sq_bytewise_deg0 = BF128_ALLOC(Nstbytes);

  // ::2
  for (unsigned int r = 0; r < R / 2; r++) {
    // ::3-4
    aes_128_f256_f2_conjugates_1(state_conj, state_bits, params);
    aes_128_f256_f2_conjugates_128(state_conj_tag, state_bits_tag, params);

    // ::5-6 : start of norms in witness
    const uint8_t* norms_ptr     = w + 3 * Nstbits * r / 2;
    const bf128_t* norm_tags_ptr = w_tag + 3 * Nstbits * r / 2;

    // ::7
    for (unsigned int i = 0; i < Nstbytes; i++) {
      // ::8-9
      bf128_t y[4];
      bf128_t y_tag[4];
      aes_128_inv_norm_to_conjugates_prover(y, y_tag, norms_ptr + 4 * i, norm_tags_ptr + 4 * i);

      // ::10-11
      aes_128_inv_norm_constraints_prover(hasher, state_conj + 8 * i, state_conj_tag + 8 * i, y,
                                          y_tag);

      // ::12
      for (unsigned int j = 0; j < 8; j++) {
        // ::13-14
        uint32_t conj_index     = i * 8 + ((j + 4) % 8);
        uint32_t y_index        = j % 4;
        st_dash_deg2[i * 8 + j] = bf128_mul(state_conj[conj_index], y[y_index]);
        st_dash_deg1[i * 8 + j] = bf128_add(bf128_mul(state_conj[conj_index], y_tag[y_index]),
                                            bf128_mul(state_conj_tag[conj_index], y[y_index]));
        st_dash_deg0[i * 8 + j] = bf128_mul(state_conj_tag[conj_index], y_tag[y_index]);
      }
    }

    // ::15-16
    bf128_t k_0_deg0[32];
    bf128_t k_0_deg1[32];
    aes_128_state_to_bytes_prover(k_0_deg1, k_0_deg0, // k is in bits
                                  k + (2 * r + 1) * Nstbits, k_tag + (2 * r + 1) * Nstbits, params);
    // ::17
    bf128_t k_1_deg0[32];
    bf128_t k_1_deg1[32];
    bf128_t k_1_deg2[32];
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      k_1_deg0[byte_i] = bf128_mul(k_0_deg0[byte_i], k_0_deg0[byte_i]);
      k_1_deg1[byte_i] = bf128_zero();
      k_1_deg2[byte_i] = bf128_mul(k_0_deg1[byte_i], k_0_deg1[byte_i]);
    }

    // ::18
    bf128_t st_b_deg0[2][32];
    bf128_t st_b_deg1[2][32];
    bf128_t st_b_deg2[2][32];
    bf128_t st_b_deg0_tmp[2][32];
    bf128_t st_b_deg1_tmp[2][32];
    bf128_t st_b_deg2_tmp[2][32];
    memset(st_b_deg0, 0x00, sizeof(st_b_deg0));
    memset(st_b_deg1, 0x00, sizeof(st_b_deg1));
    memset(st_b_deg2, 0x00, sizeof(st_b_deg2));
    memset(st_b_deg0_tmp, 0x00, sizeof(st_b_deg0_tmp));
    memset(st_b_deg1_tmp, 0x00, sizeof(st_b_deg1_tmp));
    memset(st_b_deg2_tmp, 0x00, sizeof(st_b_deg2_tmp));

    for (unsigned int b = 0; b < 2; b++) {
      // ::19
      aes_128_sbox_affine_prover(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b], st_dash_deg0,
                                 st_dash_deg1, st_dash_deg2, b, params);
      // ::20
      aes_128_shiftrows_prover(st_b_deg0_tmp[b], st_b_deg1_tmp[b], st_b_deg2_tmp[b], st_b_deg0[b],
                               st_b_deg1[b], st_b_deg2[b], params);
      memcpy(st_b_deg0[b], st_b_deg0_tmp[b], sizeof(bf128_t) * Nstbytes);
      memcpy(st_b_deg1[b], st_b_deg1_tmp[b], sizeof(bf128_t) * Nstbytes);
      memcpy(st_b_deg2[b], st_b_deg2_tmp[b], sizeof(bf128_t) * Nstbytes);
      // ::21
      aes_128_mix_columns_prover(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b], st_b_deg0_tmp[b],
                                 st_b_deg1_tmp[b], st_b_deg2_tmp[b], b, params);
      // ::22
      if (b == 0) {
        aes_128_add_round_key_bytes_prover_degree_1(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    k_0_deg0, k_0_deg1, params);
      } else {
        aes_128_add_round_key_bytes_prover_degree_2(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    k_1_deg0, k_1_deg1, k_1_deg2, params);
      }
    }
    // ::23-24
    uint8_t* s_tilde     = malloc(Nstbits * sizeof(uint8_t));
    bf128_t* s_tilde_tag = BF128_ALLOC(Nstbits);
    if (r == R / 2 - 1) {
      // ::25
      aes_128_add_round_key_prover(s_tilde, s_tilde_tag, owf_out, owf_out_tag, k + R * Nstbits,
                                   k_tag + R * Nstbits, params);
    } else {
      // ::27-28
      unsigned int idx = 0;
      for (unsigned int i = ((Nstbits / 2) + (Nstbits / 2) * 3 * r);
           i < ((Nstbits / 2 * 3 * r) + (Nstbits / 2 * 3)); i++) {
        s_tilde[idx]     = w[i];
        s_tilde_tag[idx] = w_tag[i];
        idx++;
      }
    }
    // ::29
    uint8_t* s_dash_dash     = malloc(Nstbits * sizeof(uint8_t));
    bf128_t* s_dash_dash_tag = BF128_ALLOC(Nstbits);
    aes_128_inverse_shiftrows_prover(s_dash_dash, s_dash_dash_tag, s_tilde, s_tilde_tag, params);
    // ::30
    uint8_t* s     = malloc(Nstbits * sizeof(uint8_t));
    bf128_t* s_tag = BF128_ALLOC(Nstbits);
    aes_128_inverse_affine_prover(s, s_tag, s_dash_dash, s_dash_dash_tag, params);

    // ::31
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      bf128_t s_deg1;
      bf128_t s_deg0;
      bf128_t s_sq_deg1;
      bf128_t s_sq_deg0;
      // ::32
      s_deg1 = bf128_byte_combine_bits_array(s + 8 * byte_i); // taking 8 uint8 "bits"
      s_deg0 = bf128_byte_combine(s_tag + 8 * byte_i);
      // ::33
      s_sq_deg1 = bf128_byte_combine_bits_array_sq(s + 8 * byte_i); // taking 8 uint8 "bits"
      s_sq_deg0 = bf128_byte_combine_sq(s_tag + 8 * byte_i);

      // ::36
      // compute <s^sq>^1 * <st_{0,i}>^2 - <s>^1
      //    deg0: s_sq[0] * st[0]
      //    deg1: s_sq[0] * st[1] + s_sq[1] * st[0]
      //    deg2: s_sq[0] * st[2] + s_sq[1] * st[1] + s[0]
      //
      zk_hash_128_3_update(hasher, bf128_mul(s_sq_deg0, st_b_deg0[0][byte_i]),
                           bf128_add(bf128_mul(s_sq_deg0, st_b_deg1[0][byte_i]),
                                     bf128_mul(s_sq_deg1, st_b_deg0[0][byte_i])),
                           bf128_add(bf128_add(bf128_mul(s_sq_deg0, st_b_deg2[0][byte_i]),
                                               bf128_mul(s_sq_deg1, st_b_deg1[0][byte_i])),
                                     s_deg0));

      // ::37
      // compute <s>^1 * <st_{1,i}>^2 - <st_{0,i}>^2
      //    deg0: s[0] * st_{1,i}[0]
      //    deg1: s[0] * st_{1,i}[1] + s[1] * st_{1,i}[0] + st_{0,i}[0]
      //    deg2: s[0] * st_{1,i}[2] + s[1] * st_{1,i}[1] + st_{0,i}[1]
      //
      zk_hash_128_3_update(hasher, bf128_mul(s_deg0, st_b_deg0[1][byte_i]),
                           bf128_add(bf128_add(bf128_mul(s_deg0, st_b_deg1[1][byte_i]),
                                               bf128_mul(s_deg1, st_b_deg0[1][byte_i])),
                                     st_b_deg0[0][byte_i]),
                           bf128_add(bf128_add(bf128_mul(s_deg0, st_b_deg2[1][byte_i]),
                                               bf128_mul(s_deg1, st_b_deg1[1][byte_i])),
                                     st_b_deg1[0][byte_i]));
    }
    if (r != (R / 2) - 1) {
      uint8_t* tmp_state     = malloc(Nstbits * sizeof(uint8_t));
      bf128_t* tmp_state_tag = BF128_ALLOC(Nstbits);
      aes_128_bitwise_mix_column_prover(tmp_state, tmp_state_tag, s_tilde, s_tilde_tag, params);
      aes_128_add_round_key_prover(state_bits, state_bits_tag, tmp_state, tmp_state_tag,
                                   k + (2 * r + 2) * Nstbits, k_tag + (2 * r + 2) * Nstbits,
                                   params);
      faest_aligned_free(tmp_state_tag);
      free(tmp_state);
    }

    faest_aligned_free(s_tilde_tag);
    faest_aligned_free(s_dash_dash_tag);
    faest_aligned_free(s_tag);
    free(s_tilde);
    free(s_dash_dash);
    free(s);
  }

  faest_aligned_free(state_sq_bytewise_deg0);
  faest_aligned_free(state_sq_bytewise_deg1);
  faest_aligned_free(state_sq_bytewise_deg2);
  faest_aligned_free(state_bytewise_deg0);
  faest_aligned_free(state_bytewise_deg1);
  faest_aligned_free(state_bytewise_deg2);
  faest_aligned_free(st_dash_deg0);
  faest_aligned_free(st_dash_deg1);
  faest_aligned_free(st_dash_deg2);
  faest_aligned_free(state_conj_tag);
  faest_aligned_free(state_conj);
  faest_aligned_free(state_bits_tag);
  free(state_bits);
}

static void aes_192_enc_constraints_prover(bf192_t* z_deg0, bf192_t* z_deg1, bf192_t* z_deg2,
                                           const uint8_t* owf_in, const bf192_t* owf_in_tag,
                                           const uint8_t* owf_out, const bf192_t* owf_out_tag,
                                           const uint8_t* w, const bf192_t* w_tag, const uint8_t* k,
                                           const bf192_t* k_tag, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbits  = 32 * Nst;
  const unsigned int R        = params->R;
  const unsigned int Nstbytes = Nstbits / 8;

  /// ::1 AddRoundKey
  uint8_t* state_bits     = malloc(Nstbits * sizeof(uint8_t));
  bf192_t* state_bits_tag = BF192_ALLOC(Nstbits);

  aes_192_add_round_key_prover(state_bits, state_bits_tag, owf_in, owf_in_tag, k, k_tag, params);

  // for conjugates of state and s-box outputs
  bf192_t* state_conj     = BF192_ALLOC(8 * Nstbytes);
  bf192_t* state_conj_tag = BF192_ALLOC(8 * Nstbytes);
  bf192_t* st_dash_deg2   = BF192_ALLOC(8 * Nstbytes);
  bf192_t* st_dash_deg1   = BF192_ALLOC(8 * Nstbytes);
  bf192_t* st_dash_deg0   = BF192_ALLOC(8 * Nstbytes);

  bf192_t* state_bytewise_deg2    = BF192_ALLOC(Nstbytes);
  bf192_t* state_bytewise_deg1    = BF192_ALLOC(Nstbytes);
  bf192_t* state_bytewise_deg0    = BF192_ALLOC(Nstbytes);
  bf192_t* state_sq_bytewise_deg2 = BF192_ALLOC(Nstbytes);
  bf192_t* state_sq_bytewise_deg1 = BF192_ALLOC(Nstbytes);
  bf192_t* state_sq_bytewise_deg0 = BF192_ALLOC(Nstbytes);

  // ::2
  for (unsigned int r = 0; r < R / 2; r++) {
    // ::3-4
    aes_192_f256_f2_conjugates_1(state_conj, state_bits, params);
    aes_192_f256_f2_conjugates_192(state_conj_tag, state_bits_tag, params);

    // ::5-6 : start of norms in witness
    const uint8_t* norms_ptr     = w + 3 * Nstbits * r / 2;
    const bf192_t* norm_tags_ptr = w_tag + 3 * Nstbits * r / 2;
    // ::7
    for (unsigned int i = 0; i < Nstbytes; i++) {
      // ::8-9
      bf192_t y[4];
      bf192_t y_tag[4];
      aes_192_inv_norm_to_conjugates_prover(y, y_tag, norms_ptr + 4 * i, norm_tags_ptr + 4 * i);

      // ::10-11
      aes_192_inv_norm_constraints_prover(
          z_deg0 + 3 * r * Nstbytes + i, z_deg1 + 3 * r * Nstbytes + i,
          z_deg2 + 3 * r * Nstbytes + i, state_conj + 8 * i, state_conj_tag + 8 * i, y, y_tag);

      // ::12
      for (unsigned int j = 0; j < 8; j++) {
        // ::13-14
        uint32_t conj_index     = i * 8 + ((j + 4) % 8);
        uint32_t y_index        = j % 4;
        st_dash_deg2[i * 8 + j] = bf192_mul(state_conj[conj_index], y[y_index]);
        st_dash_deg1[i * 8 + j] = bf192_add(bf192_mul(state_conj[conj_index], y_tag[y_index]),
                                            bf192_mul(state_conj_tag[conj_index], y[y_index]));
        st_dash_deg0[i * 8 + j] = bf192_mul(state_conj_tag[conj_index], y_tag[y_index]);
      }
    }

    // ::15-16
    bf192_t k_0_deg0[32];
    bf192_t k_0_deg1[32];
    aes_192_state_to_bytes_prover(k_0_deg1, k_0_deg0, // k is in bits
                                  k + (2 * r + 1) * Nstbits, k_tag + (2 * r + 1) * Nstbits, params);
    // ::17
    bf192_t k_1_deg0[32];
    bf192_t k_1_deg1[32];
    bf192_t k_1_deg2[32];
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      k_1_deg0[byte_i] = bf192_mul(k_0_deg0[byte_i], k_0_deg0[byte_i]);
      k_1_deg1[byte_i] = bf192_zero();
      k_1_deg2[byte_i] = bf192_mul(k_0_deg1[byte_i], k_0_deg1[byte_i]);
    }

    // ::18
    bf192_t st_b_deg0[2][32];
    bf192_t st_b_deg1[2][32];
    bf192_t st_b_deg2[2][32];
    bf192_t st_b_deg0_tmp[2][32];
    bf192_t st_b_deg1_tmp[2][32];
    bf192_t st_b_deg2_tmp[2][32];
    memset(st_b_deg0, 0x00, sizeof(st_b_deg0));
    memset(st_b_deg1, 0x00, sizeof(st_b_deg1));
    memset(st_b_deg2, 0x00, sizeof(st_b_deg2));
    memset(st_b_deg0_tmp, 0x00, sizeof(st_b_deg0_tmp));
    memset(st_b_deg1_tmp, 0x00, sizeof(st_b_deg1_tmp));
    memset(st_b_deg2_tmp, 0x00, sizeof(st_b_deg2_tmp));

    for (unsigned int b = 0; b < 2; b++) {
      // ::19
      aes_192_sbox_affine_prover(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b], st_dash_deg0,
                                 st_dash_deg1, st_dash_deg2, b, params);
      // ::20
      aes_192_shiftrows_prover(st_b_deg0_tmp[b], st_b_deg1_tmp[b], st_b_deg2_tmp[b], st_b_deg0[b],
                               st_b_deg1[b], st_b_deg2[b], params);
      memcpy(st_b_deg0[b], st_b_deg0_tmp[b], sizeof(bf192_t) * Nstbytes);
      memcpy(st_b_deg1[b], st_b_deg1_tmp[b], sizeof(bf192_t) * Nstbytes);
      memcpy(st_b_deg2[b], st_b_deg2_tmp[b], sizeof(bf192_t) * Nstbytes);
      // ::21
      aes_192_mix_columns_prover(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b], st_b_deg0_tmp[b],
                                 st_b_deg1_tmp[b], st_b_deg2_tmp[b], b, params);
      // ::22
      if (b == 0) {
        aes_192_add_round_key_bytes_prover_degree_1(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    k_0_deg0, k_0_deg1, params);
      } else {
        aes_192_add_round_key_bytes_prover_degree_2(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    k_1_deg0, k_1_deg1, k_1_deg2, params);
      }
    }
    // ::23-24
    uint8_t* s_tilde     = malloc(Nstbits * sizeof(uint8_t));
    bf192_t* s_tilde_tag = BF192_ALLOC(Nstbits);
    if (r == R / 2 - 1) {
      // ::25
      aes_192_add_round_key_prover(s_tilde, s_tilde_tag, owf_out, owf_out_tag, k + R * Nstbits,
                                   k_tag + R * Nstbits, params);
    } else {
      // ::27-28
      unsigned int idx = 0;
      for (unsigned int i = ((Nstbits / 2) + (Nstbits / 2) * 3 * r);
           i < ((Nstbits / 2 * 3 * r) + (Nstbits / 2 * 3)); i++) {
        s_tilde[idx]     = w[i];
        s_tilde_tag[idx] = w_tag[i];
        idx++;
      }
    }
    // ::29
    uint8_t* s_dash_dash     = malloc(Nstbits * sizeof(uint8_t));
    bf192_t* s_dash_dash_tag = BF192_ALLOC(Nstbits);
    aes_192_inverse_shiftrows_prover(s_dash_dash, s_dash_dash_tag, s_tilde, s_tilde_tag, params);
    // ::30
    uint8_t* s     = malloc(Nstbits * sizeof(uint8_t));
    bf192_t* s_tag = BF192_ALLOC(Nstbits);
    aes_192_inverse_affine_prover(s, s_tag, s_dash_dash, s_dash_dash_tag, params);

    // ::31
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      bf192_t s_deg1;
      bf192_t s_deg0;
      bf192_t s_sq_deg1;
      bf192_t s_sq_deg0;
      // ::32
      s_deg1 = bf192_byte_combine_bits_array(s + 8 * byte_i); // taking 8 uint8 "bits"
      s_deg0 = bf192_byte_combine(s_tag + 8 * byte_i);
      // ::33
      s_sq_deg1 = bf192_byte_combine_bits_array_sq(s + 8 * byte_i); // taking 8 uint8 "bits"
      s_sq_deg0 = bf192_byte_combine_sq(s_tag + 8 * byte_i);

      // ::36
      // compute <s^sq>^1 * <st_{0,i}>^2 - <s>^1
      //    deg0: s_sq[0] * st[0]
      //    deg1: s_sq[0] * st[1] + s_sq[1] * st[0]
      //    deg2: s_sq[0] * st[2] + s_sq[1] * st[1] + s[0]
      //

      z_deg0[(3 * r + 1) * Nstbytes + 2 * byte_i] = bf192_mul(s_sq_deg0, st_b_deg0[0][byte_i]);

      z_deg1[(3 * r + 1) * Nstbytes + 2 * byte_i] = bf192_add(
          bf192_mul(s_sq_deg0, st_b_deg1[0][byte_i]), bf192_mul(s_sq_deg1, st_b_deg0[0][byte_i]));
      z_deg2[(3 * r + 1) * Nstbytes + 2 * byte_i] =
          bf192_add(bf192_add(bf192_mul(s_sq_deg0, st_b_deg2[0][byte_i]),
                              bf192_mul(s_sq_deg1, st_b_deg1[0][byte_i])),
                    s_deg0);
      // ::37
      // compute <s>^1 * <st_{1,i}>^2 - <st_{0,i}>^2
      //    deg0: s[0] * st_{1,i}[0]
      //    deg1: s[0] * st_{1,i}[1] + s[1] * st_{1,i}[0] + st_{0,i}[0]
      //    deg2: s[0] * st_{1,i}[2] + s[1] * st_{1,i}[1] + st_{0,i}[1]
      //
      z_deg0[(3 * r + 1) * Nstbytes + 2 * byte_i + 1] = bf192_mul(s_deg0, st_b_deg0[1][byte_i]);
      z_deg1[(3 * r + 1) * Nstbytes + 2 * byte_i + 1] =
          bf192_add(bf192_add(bf192_mul(s_deg0, st_b_deg1[1][byte_i]),
                              bf192_mul(s_deg1, st_b_deg0[1][byte_i])),
                    st_b_deg0[0][byte_i]);
      z_deg2[(3 * r + 1) * Nstbytes + 2 * byte_i + 1] =
          bf192_add(bf192_add(bf192_mul(s_deg0, st_b_deg2[1][byte_i]),
                              bf192_mul(s_deg1, st_b_deg1[1][byte_i])),
                    st_b_deg1[0][byte_i]);
    }
    if (r != (R / 2) - 1) {
      uint8_t* tmp_state     = malloc(Nstbits * sizeof(uint8_t));
      bf192_t* tmp_state_tag = BF192_ALLOC(Nstbits);
      aes_192_bitwise_mix_column_prover(tmp_state, tmp_state_tag, s_tilde, s_tilde_tag, params);
      aes_192_add_round_key_prover(state_bits, state_bits_tag, tmp_state, tmp_state_tag,
                                   k + (2 * r + 2) * Nstbits, k_tag + (2 * r + 2) * Nstbits,
                                   params);
      faest_aligned_free(tmp_state_tag);
      free(tmp_state);
    }

    faest_aligned_free(s_tilde_tag);
    faest_aligned_free(s_dash_dash_tag);
    faest_aligned_free(s_tag);
    free(s_tilde);
    free(s_dash_dash);
    free(s);
  }

  faest_aligned_free(state_sq_bytewise_deg0);
  faest_aligned_free(state_sq_bytewise_deg1);
  faest_aligned_free(state_sq_bytewise_deg2);
  faest_aligned_free(state_bytewise_deg0);
  faest_aligned_free(state_bytewise_deg1);
  faest_aligned_free(state_bytewise_deg2);
  faest_aligned_free(st_dash_deg0);
  faest_aligned_free(st_dash_deg1);
  faest_aligned_free(st_dash_deg2);
  faest_aligned_free(state_conj_tag);
  faest_aligned_free(state_conj);
  faest_aligned_free(state_bits_tag);
  free(state_bits);
}

static void aes_256_enc_constraints_prover(bf256_t* z_deg0, bf256_t* z_deg1, bf256_t* z_deg2,
                                           const uint8_t* owf_in, const bf256_t* owf_in_tag,
                                           const uint8_t* owf_out, const bf256_t* owf_out_tag,
                                           const uint8_t* w, const bf256_t* w_tag, const uint8_t* k,
                                           const bf256_t* k_tag, const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbits  = 32 * Nst;
  const unsigned int R        = params->R;
  const unsigned int Nstbytes = Nstbits / 8;

  /// ::1 AddRoundKey
  uint8_t* state_bits     = malloc(Nstbits * sizeof(uint8_t));
  bf256_t* state_bits_tag = BF256_ALLOC(Nstbits);

  aes_256_add_round_key_prover(state_bits, state_bits_tag, owf_in, owf_in_tag, k, k_tag, params);

  // for conjugates of state and s-box outputs
  bf256_t* state_conj     = BF256_ALLOC(8 * Nstbytes);
  bf256_t* state_conj_tag = BF256_ALLOC(8 * Nstbytes);
  bf256_t* st_dash_deg2   = BF256_ALLOC(8 * Nstbytes);
  bf256_t* st_dash_deg1   = BF256_ALLOC(8 * Nstbytes);
  bf256_t* st_dash_deg0   = BF256_ALLOC(8 * Nstbytes);

  bf256_t* state_bytewise_deg2    = BF256_ALLOC(Nstbytes);
  bf256_t* state_bytewise_deg1    = BF256_ALLOC(Nstbytes);
  bf256_t* state_bytewise_deg0    = BF256_ALLOC(Nstbytes);
  bf256_t* state_sq_bytewise_deg2 = BF256_ALLOC(Nstbytes);
  bf256_t* state_sq_bytewise_deg1 = BF256_ALLOC(Nstbytes);
  bf256_t* state_sq_bytewise_deg0 = BF256_ALLOC(Nstbytes);

  // ::2
  for (unsigned int r = 0; r < R / 2; r++) {
    // ::3-4
    aes_256_f256_f2_conjugates_1(state_conj, state_bits, params);
    aes_256_f256_f2_conjugates_256(state_conj_tag, state_bits_tag, params);

    // ::5-6 : start of norms in witness
    const uint8_t* norms_ptr     = w + 3 * Nstbits * r / 2;
    const bf256_t* norm_tags_ptr = w_tag + 3 * Nstbits * r / 2;

    // ::7
    for (unsigned int i = 0; i < Nstbytes; i++) {
      // ::8-9
      bf256_t y[4];
      bf256_t y_tag[4];
      aes_256_inv_norm_to_conjugates_prover(y, y_tag, norms_ptr + 4 * i, norm_tags_ptr + 4 * i);

      // ::10-11
      aes_256_inv_norm_constraints_prover(
          z_deg0 + 3 * r * Nstbytes + i, z_deg1 + 3 * r * Nstbytes + i,
          z_deg2 + 3 * r * Nstbytes + i, state_conj + 8 * i, state_conj_tag + 8 * i, y, y_tag);

      // ::12
      for (unsigned int j = 0; j < 8; j++) {
        // ::13-14
        uint32_t conj_index     = i * 8 + ((j + 4) % 8);
        uint32_t y_index        = j % 4;
        st_dash_deg2[i * 8 + j] = bf256_mul(state_conj[conj_index], y[y_index]);
        st_dash_deg1[i * 8 + j] = bf256_add(bf256_mul(state_conj[conj_index], y_tag[y_index]),
                                            bf256_mul(state_conj_tag[conj_index], y[y_index]));
        st_dash_deg0[i * 8 + j] = bf256_mul(state_conj_tag[conj_index], y_tag[y_index]);
      }
    }

    // ::15-16
    bf256_t k_0_deg0[32];
    bf256_t k_0_deg1[32];
    aes_256_state_to_bytes_prover(k_0_deg1, k_0_deg0, k + (2 * r + 1) * Nstbits,
                                  k_tag + (2 * r + 1) * Nstbits, params);
    // ::17
    bf256_t k_1_deg0[32];
    bf256_t k_1_deg1[32];
    bf256_t k_1_deg2[32];
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      k_1_deg0[byte_i] = bf256_mul(k_0_deg0[byte_i], k_0_deg0[byte_i]);
      k_1_deg1[byte_i] = bf256_zero();
      k_1_deg2[byte_i] = bf256_mul(k_0_deg1[byte_i], k_0_deg1[byte_i]);
    }

    // ::18
    bf256_t st_b_deg0[2][32];
    bf256_t st_b_deg1[2][32];
    bf256_t st_b_deg2[2][32];
    bf256_t st_b_deg0_tmp[2][32];
    bf256_t st_b_deg1_tmp[2][32];
    bf256_t st_b_deg2_tmp[2][32];
    memset(st_b_deg0, 0x00, sizeof(st_b_deg0));
    memset(st_b_deg1, 0x00, sizeof(st_b_deg1));
    memset(st_b_deg2, 0x00, sizeof(st_b_deg2));
    memset(st_b_deg0_tmp, 0x00, sizeof(st_b_deg0_tmp));
    memset(st_b_deg1_tmp, 0x00, sizeof(st_b_deg1_tmp));
    memset(st_b_deg2_tmp, 0x00, sizeof(st_b_deg2_tmp));

    for (unsigned int b = 0; b < 2; b++) {
      // ::19
      aes_256_sbox_affine_prover(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b], st_dash_deg0,
                                 st_dash_deg1, st_dash_deg2, b, params);
      // ::20
      aes_256_shiftrows_prover(st_b_deg0_tmp[b], st_b_deg1_tmp[b], st_b_deg2_tmp[b], st_b_deg0[b],
                               st_b_deg1[b], st_b_deg2[b], params);
      memcpy(st_b_deg0[b], st_b_deg0_tmp[b], sizeof(bf256_t) * Nstbytes);
      memcpy(st_b_deg1[b], st_b_deg1_tmp[b], sizeof(bf256_t) * Nstbytes);
      memcpy(st_b_deg2[b], st_b_deg2_tmp[b], sizeof(bf256_t) * Nstbytes);
      // ::21
      aes_256_mix_columns_prover(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b], st_b_deg0_tmp[b],
                                 st_b_deg1_tmp[b], st_b_deg2_tmp[b], b, params);
      // ::22
      if (b == 0) {
        aes_256_add_round_key_bytes_prover_degree_1(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    k_0_deg0, k_0_deg1, params);
      } else {
        aes_256_add_round_key_bytes_prover_degree_2(st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    st_b_deg0[b], st_b_deg1[b], st_b_deg2[b],
                                                    k_1_deg0, k_1_deg1, k_1_deg2, params);
      }
    }
    // ::23-24
    uint8_t* s_tilde     = malloc(Nstbits * sizeof(uint8_t));
    bf256_t* s_tilde_tag = BF256_ALLOC(Nstbits);
    if (r == R / 2 - 1) {
      // ::25
      aes_256_add_round_key_prover(s_tilde, s_tilde_tag, owf_out, owf_out_tag, k + R * Nstbits,
                                   k_tag + R * Nstbits, params);
    } else {
      // ::27-28
      unsigned int idx = 0;
      for (unsigned int i = ((Nstbits / 2) + (Nstbits / 2) * 3 * r);
           i < ((Nstbits / 2 * 3 * r) + (Nstbits / 2 * 3)); i++) {
        s_tilde[idx]     = w[i];
        s_tilde_tag[idx] = w_tag[i];
        idx++;
      }
    }
    // ::29
    uint8_t* s_dash_dash     = malloc(Nstbits * sizeof(uint8_t));
    bf256_t* s_dash_dash_tag = BF256_ALLOC(Nstbits);
    aes_256_inverse_shiftrows_prover(s_dash_dash, s_dash_dash_tag, s_tilde, s_tilde_tag, params);
    // ::30
    uint8_t* s     = malloc(Nstbits * sizeof(uint8_t));
    bf256_t* s_tag = BF256_ALLOC(Nstbits);
    aes_256_inverse_affine_prover(s, s_tag, s_dash_dash, s_dash_dash_tag, params);

    // ::31
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      bf256_t s_deg1;
      bf256_t s_deg0;
      bf256_t s_sq_deg1;
      bf256_t s_sq_deg0;
      // ::32
      s_deg1 = bf256_byte_combine_bits_array(s + 8 * byte_i); // taking 8 uint8 "bits"
      s_deg0 = bf256_byte_combine(s_tag + 8 * byte_i);
      // ::33
      s_sq_deg1 = bf256_byte_combine_bits_array_sq(s + 8 * byte_i); // taking 8 uint8 "bits"
      s_sq_deg0 = bf256_byte_combine_sq(s_tag + 8 * byte_i);

      // ::36
      // compute <s^sq>^1 * <st_{0,i}>^2 - <s>^1
      //    deg0: s_sq[0] * st[0]
      //    deg1: s_sq[0] * st[1] + s_sq[1] * st[0]
      //    deg2: s_sq[0] * st[2] + s_sq[1] * st[1] + s[0]
      //

      z_deg0[(3 * r + 1) * Nstbytes + 2 * byte_i] = bf256_mul(s_sq_deg0, st_b_deg0[0][byte_i]);

      z_deg1[(3 * r + 1) * Nstbytes + 2 * byte_i] = bf256_add(
          bf256_mul(s_sq_deg0, st_b_deg1[0][byte_i]), bf256_mul(s_sq_deg1, st_b_deg0[0][byte_i]));
      z_deg2[(3 * r + 1) * Nstbytes + 2 * byte_i] =
          bf256_add(bf256_add(bf256_mul(s_sq_deg0, st_b_deg2[0][byte_i]),
                              bf256_mul(s_sq_deg1, st_b_deg1[0][byte_i])),
                    s_deg0);
      // ::37
      // compute <s>^1 * <st_{1,i}>^2 - <st_{0,i}>^2
      //    deg0: s[0] * st_{1,i}[0]
      //    deg1: s[0] * st_{1,i}[1] + s[1] * st_{1,i}[0] + st_{0,i}[0]
      //    deg2: s[0] * st_{1,i}[2] + s[1] * st_{1,i}[1] + st_{0,i}[1]
      //
      z_deg0[(3 * r + 1) * Nstbytes + 2 * byte_i + 1] = bf256_mul(s_deg0, st_b_deg0[1][byte_i]);
      z_deg1[(3 * r + 1) * Nstbytes + 2 * byte_i + 1] =
          bf256_add(bf256_add(bf256_mul(s_deg0, st_b_deg1[1][byte_i]),
                              bf256_mul(s_deg1, st_b_deg0[1][byte_i])),
                    st_b_deg0[0][byte_i]);
      z_deg2[(3 * r + 1) * Nstbytes + 2 * byte_i + 1] =
          bf256_add(bf256_add(bf256_mul(s_deg0, st_b_deg2[1][byte_i]),
                              bf256_mul(s_deg1, st_b_deg1[1][byte_i])),
                    st_b_deg1[0][byte_i]);
    }
    if (r != (R / 2) - 1) {
      uint8_t* tmp_state     = malloc(Nstbits * sizeof(uint8_t));
      bf256_t* tmp_state_tag = BF256_ALLOC(Nstbits);
      aes_256_bitwise_mix_column_prover(tmp_state, tmp_state_tag, s_tilde, s_tilde_tag, params);
      aes_256_add_round_key_prover(state_bits, state_bits_tag, tmp_state, tmp_state_tag,
                                   k + (2 * r + 2) * Nstbits, k_tag + (2 * r + 2) * Nstbits,
                                   params);
      faest_aligned_free(tmp_state_tag);
      free(tmp_state);
    }

    faest_aligned_free(s_tilde_tag);
    faest_aligned_free(s_dash_dash_tag);
    faest_aligned_free(s_tag);
    free(s_tilde);
    free(s_dash_dash);
    free(s);
  }

  faest_aligned_free(state_sq_bytewise_deg0);
  faest_aligned_free(state_sq_bytewise_deg1);
  faest_aligned_free(state_sq_bytewise_deg2);
  faest_aligned_free(state_bytewise_deg0);
  faest_aligned_free(state_bytewise_deg1);
  faest_aligned_free(state_bytewise_deg2);
  faest_aligned_free(st_dash_deg0);
  faest_aligned_free(st_dash_deg1);
  faest_aligned_free(st_dash_deg2);
  faest_aligned_free(state_conj_tag);
  faest_aligned_free(state_conj);
  faest_aligned_free(state_bits_tag);
  free(state_bits);
}

static void aes_128_enc_constraints_verifier(zk_hash_128_ctx* hasher, const bf128_t* owf_in_key,
                                             const bf128_t* owf_out_key, const bf128_t* w_key,
                                             const bf128_t* rkeys_key, const bf128_t delta,
                                             const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbits  = 32 * Nst;
  const unsigned int R        = params->R;
  const unsigned int Nstbytes = Nstbits / 8;

  bf128_t* state_bits_key = BF128_ALLOC(Nstbits);

  /// ::1 AddRoundKey
  aes_128_add_round_key_verifier(state_bits_key, owf_in_key, rkeys_key, params);

  // for conjugates of state and s-box outputs
  bf128_t* state_conj_key = BF128_ALLOC(8 * Nstbytes);
  bf128_t* st_dash_key    = BF128_ALLOC(8 * Nstbytes);

  // ::2
  for (unsigned int r = 0; r < R / 2; r++) {
    // ::3-4
    aes_128_f256_f2_conjugates_128(state_conj_key, state_bits_key, params);

    // ::5-6 : start of norms in witness
    const bf128_t* norm_keys_ptr = w_key + 3 * Nstbits * r / 2;

    // ::7
    for (unsigned int i = 0; i < Nstbytes; i++) {
      // ::8-9
      bf128_t y_key[4];
      aes_128_inv_norm_to_conjugates_verifier(y_key, norm_keys_ptr + 4 * i);

      // ::10-11
      aes_128_inv_norm_constraints_verifier(hasher, state_conj_key + 8 * i, y_key, delta);

      // ::12
      for (unsigned int j = 0; j < 8; j++) {
        // ::13-14
        uint32_t conj_index    = i * 8 + ((j + 4) % 8);
        uint32_t y_index       = j % 4;
        st_dash_key[i * 8 + j] = bf128_mul(state_conj_key[conj_index], y_key[y_index]);
      }
    }

    // ::15-16
    bf128_t k_0_key[32];
    bf128_t k_1_key[32];
    aes_128_state_to_bytes_verifier(k_0_key, rkeys_key + (2 * r + 1) * Nstbits, params);
    // ::17
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      k_1_key[byte_i] = bf128_mul(k_0_key[byte_i], k_0_key[byte_i]);
    }

    // ::18
    bf128_t st_b_key[2][32];
    bf128_t st_b_tmp_key[2][32];
    memset(st_b_key, 0x00, sizeof(st_b_key));
    memset(st_b_tmp_key, 0x00, sizeof(st_b_tmp_key));

    for (unsigned int b = 0; b < 2; b++) {
      // ::19
      aes_128_sbox_affine_verifier(st_b_key[b], st_dash_key, delta, b, params);
      // ::20
      aes_128_shiftrows_verifier(st_b_tmp_key[b], st_b_key[b], params);
      memcpy(st_b_key[b], st_b_tmp_key[b], sizeof(bf128_t) * Nstbytes);
      // ::21
      aes_128_mix_columns_verifier(st_b_tmp_key[b], st_b_key[b], b, params);
      memcpy(st_b_key[b], st_b_tmp_key[b], sizeof(bf128_t) * Nstbytes);
      // ::22
      if (b == 0) {
        aes_128_add_round_key_bytes_verifier(st_b_key[b], st_b_key[b], k_0_key, delta, true,
                                             params);
      } else {
        aes_128_add_round_key_bytes_verifier(st_b_key[b], st_b_key[b], k_1_key, delta, false,
                                             params);
      }
    }
    // ::23-24
    bf128_t* s_tilde_key = BF128_ALLOC(Nstbits);
    if (r == R / 2 - 1) {
      // ::25
      aes_128_add_round_key_verifier(s_tilde_key, owf_out_key, rkeys_key + R * Nstbits, params);
    } else {
      // ::27-28
      unsigned int idx = 0;
      for (unsigned int i = ((Nstbits / 2) + (Nstbits / 2) * 3 * r);
           i < ((Nstbits / 2 * 3 * r) + (Nstbits / 2 * 3)); i++) {
        s_tilde_key[idx] = w_key[i];
        idx++;
      }
    }
    // ::29
    bf128_t* s_dash_dash_key = BF128_ALLOC(Nstbits);
    aes_128_inverse_shiftrows_verifier(s_dash_dash_key, s_tilde_key, params);
    // ::30
    bf128_t* s_state_key = BF128_ALLOC(Nstbits);
    aes_128_inverse_affine_verifier(s_state_key, s_dash_dash_key, delta, params);

    // ::31
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      bf128_t s_key;
      bf128_t s_sq_key;
      // ::32
      s_key = bf128_byte_combine(s_state_key + 8 * byte_i);
      // ::33
      s_sq_key = bf128_byte_combine_sq(s_state_key + 8 * byte_i);

      // ::36
      // compute <s^sq>^1 * <st_{0,i}>^2 - <s>^1
      //    s_sq * st_{0,i} + delta^2 * s
      //
      zk_hash_128_update(hasher, bf128_add(bf128_mul(s_sq_key, st_b_key[0][byte_i]),
                                           bf128_mul(delta, bf128_mul(delta, s_key))));
      // ::37
      // compute <s>^1 * <st_{1,i}>^2 - <st_{0,i}>^2
      //    s * st_{1,i} + delta * st_{0,i}
      //
      zk_hash_128_update(hasher, bf128_add(bf128_mul(s_key, st_b_key[1][byte_i]),
                                           bf128_mul(delta, st_b_key[0][byte_i])));
    }
    if (r != (R / 2) - 1) {
      bf128_t* tmp_state_key = BF128_ALLOC(Nstbits);
      aes_128_bitwise_mix_column_verifier(tmp_state_key, s_tilde_key, params);
      aes_128_add_round_key_verifier(state_bits_key, tmp_state_key,
                                     rkeys_key + (2 * r + 2) * Nstbits, params);
      faest_aligned_free(tmp_state_key);
    }
    faest_aligned_free(s_tilde_key);
    faest_aligned_free(s_dash_dash_key);
    faest_aligned_free(s_state_key);
  }
  faest_aligned_free(st_dash_key);
  faest_aligned_free(state_conj_key);
  faest_aligned_free(state_bits_key);
}

static void aes_192_enc_constraints_verifier(bf192_t* z_key, const bf192_t* owf_in_key,
                                             const bf192_t* owf_out_key, const bf192_t* w_key,
                                             const bf192_t* rkeys_key, const bf192_t delta,
                                             const faest_paramset_t* params) {
  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbits  = 32 * Nst;
  const unsigned int R        = params->R;
  const unsigned int Nstbytes = Nstbits / 8;

  bf192_t* state_bits_key = BF192_ALLOC(Nstbits);

  /// ::1 AddRoundKey
  aes_192_add_round_key_verifier(state_bits_key, owf_in_key, rkeys_key, params);

  // for conjugates of state and s-box outputs
  bf192_t* state_conj_key = BF192_ALLOC(8 * Nstbytes);
  bf192_t* st_dash_key    = BF192_ALLOC(8 * Nstbytes);

  // ::2
  for (unsigned int r = 0; r < R / 2; r++) {
    // ::3-4
    aes_192_f256_f2_conjugates_192(state_conj_key, state_bits_key, params);

    // ::5-6 : start of norms in witness
    const bf192_t* norm_keys_ptr = w_key + 3 * Nstbits * r / 2;

    // ::7
    for (unsigned int i = 0; i < Nstbytes; i++) {
      // ::8-9
      bf192_t y_key[4];
      aes_192_inv_norm_to_conjugates_verifier(y_key, norm_keys_ptr + 4 * i);

      // ::10-11
      aes_192_inv_norm_constraints_verifier(z_key + 3 * r * Nstbytes + i, state_conj_key + 8 * i,
                                            y_key, delta);

      // ::12
      for (unsigned int j = 0; j < 8; j++) {
        // ::13-14
        uint32_t conj_index    = i * 8 + ((j + 4) % 8);
        uint32_t y_index       = j % 4;
        st_dash_key[i * 8 + j] = bf192_mul(state_conj_key[conj_index], y_key[y_index]);
      }
    }

    // ::15-16
    bf192_t k_0_key[32];
    bf192_t k_1_key[32];
    aes_192_state_to_bytes_verifier(k_0_key, rkeys_key + (2 * r + 1) * Nstbits, params);
    // ::17
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      k_1_key[byte_i] = bf192_mul(k_0_key[byte_i], k_0_key[byte_i]);
    }

    // ::18
    bf192_t st_b_key[2][32];
    bf192_t st_b_tmp_key[2][32];
    memset(st_b_key, 0x00, sizeof(st_b_key));
    memset(st_b_tmp_key, 0x00, sizeof(st_b_tmp_key));

    for (unsigned int b = 0; b < 2; b++) {
      // ::19
      aes_192_sbox_affine_verifier(st_b_key[b], st_dash_key, delta, b, params);
      // ::20
      aes_192_shiftrows_verifier(st_b_tmp_key[b], st_b_key[b], params);
      memcpy(st_b_key[b], st_b_tmp_key[b], sizeof(bf192_t) * Nstbytes);
      // ::21
      aes_192_mix_columns_verifier(st_b_tmp_key[b], st_b_key[b], b, params);
      memcpy(st_b_key[b], st_b_tmp_key[b], sizeof(bf192_t) * Nstbytes);
      // ::22
      if (b == 0) {
        aes_192_add_round_key_bytes_verifier(st_b_key[b], st_b_key[b], k_0_key, delta, true,
                                             params);
      } else {
        aes_192_add_round_key_bytes_verifier(st_b_key[b], st_b_key[b], k_1_key, delta, false,
                                             params);
      }
    }
    // ::23-24
    bf192_t* s_tilde_key = BF192_ALLOC(Nstbits);
    if (r == R / 2 - 1) {
      // ::25
      aes_192_add_round_key_verifier(s_tilde_key, owf_out_key, rkeys_key + R * Nstbits, params);
    } else {
      // ::27-28
      unsigned int idx = 0;
      for (unsigned int i = ((Nstbits / 2) + (Nstbits / 2) * 3 * r);
           i < ((Nstbits / 2 * 3 * r) + (Nstbits / 2 * 3)); i++) {
        s_tilde_key[idx] = w_key[i];
        idx++;
      }
    }
    // ::29
    bf192_t* s_dash_dash_key = BF192_ALLOC(Nstbits);
    aes_192_inverse_shiftrows_verifier(s_dash_dash_key, s_tilde_key, params);
    // ::30
    bf192_t* s_state_key = BF192_ALLOC(Nstbits);
    aes_192_inverse_affine_verifier(s_state_key, s_dash_dash_key, delta, params);

    // ::31
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      bf192_t s_key;
      bf192_t s_sq_key;
      // ::32
      s_key = bf192_byte_combine(s_state_key + 8 * byte_i);
      // ::33
      s_sq_key = bf192_byte_combine_sq(s_state_key + 8 * byte_i);

      // ::36
      // compute <s^sq>^1 * <st_{0,i}>^2 - <s>^1
      //    s_sq * st_{0,i} + delta^2 * s
      //
      z_key[(3 * r + 1) * Nstbytes + 2 * byte_i] = bf192_add(
          bf192_mul(s_sq_key, st_b_key[0][byte_i]), bf192_mul(delta, bf192_mul(delta, s_key)));
      // ::37
      // compute <s>^1 * <st_{1,i}>^2 - <st_{0,i}>^2
      //    s * st_{1,i} + delta * st_{0,i}
      //
      z_key[(3 * r + 1) * Nstbytes + 2 * byte_i + 1] =
          bf192_add(bf192_mul(s_key, st_b_key[1][byte_i]), bf192_mul(delta, st_b_key[0][byte_i]));
    }
    if (r != (R / 2) - 1) {
      bf192_t* tmp_state_key = BF192_ALLOC(Nstbits);
      aes_192_bitwise_mix_column_verifier(tmp_state_key, s_tilde_key, params);
      aes_192_add_round_key_verifier(state_bits_key, tmp_state_key,
                                     rkeys_key + (2 * r + 2) * Nstbits, params);
      faest_aligned_free(tmp_state_key);
    }
    faest_aligned_free(s_tilde_key);
    faest_aligned_free(s_dash_dash_key);
    faest_aligned_free(s_state_key);
  }
  faest_aligned_free(st_dash_key);
  faest_aligned_free(state_conj_key);
  faest_aligned_free(state_bits_key);
}

static void aes_256_enc_constraints_verifier(bf256_t* z_key, const bf256_t* owf_in_key,
                                             const bf256_t* owf_out_key, const bf256_t* w_key,
                                             const bf256_t* rkeys_key, const bf256_t delta,
                                             const faest_paramset_t* params) {

  const unsigned int Nst      = params->Nst;
  const unsigned int Nstbits  = 32 * Nst;
  const unsigned int R        = params->R;
  const unsigned int Nstbytes = Nstbits / 8;

  bf256_t* state_bits_key = BF256_ALLOC(Nstbits);

  /// ::1 AddRoundKey
  aes_256_add_round_key_verifier(state_bits_key, owf_in_key, rkeys_key, params);

  // for conjugates of state and s-box outputs
  bf256_t* state_conj_key = BF256_ALLOC(8 * Nstbytes);
  bf256_t* st_dash_key    = BF256_ALLOC(8 * Nstbytes);

  // ::2
  for (unsigned int r = 0; r < R / 2; r++) {
    // ::3-4
    aes_256_f256_f2_conjugates_256(state_conj_key, state_bits_key, params);

    // ::5-6 : start of norms in witness
    const bf256_t* norm_keys_ptr = w_key + 3 * Nstbits * r / 2;

    // ::7
    for (unsigned int i = 0; i < Nstbytes; i++) {
      // ::8-9
      bf256_t y_key[4];
      aes_256_inv_norm_to_conjugates_verifier(y_key, norm_keys_ptr + 4 * i);

      // ::10-11
      aes_256_inv_norm_constraints_verifier(z_key + 3 * r * Nstbytes + i, state_conj_key + 8 * i,
                                            y_key, delta);

      // ::12
      for (unsigned int j = 0; j < 8; j++) {
        // ::13-14
        uint32_t conj_index    = i * 8 + ((j + 4) % 8);
        uint32_t y_index       = j % 4;
        st_dash_key[i * 8 + j] = bf256_mul(state_conj_key[conj_index], y_key[y_index]);
      }
    }

    // ::15-16
    bf256_t k_0_key[32];
    bf256_t k_1_key[32];
    aes_256_state_to_bytes_verifier(k_0_key, rkeys_key + (2 * r + 1) * Nstbits, params);
    // ::17
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      k_1_key[byte_i] = bf256_mul(k_0_key[byte_i], k_0_key[byte_i]);
    }

    // ::18
    bf256_t st_b_key[2][32];
    bf256_t st_b_tmp_key[2][32];
    memset(st_b_key, 0x00, sizeof(st_b_key));
    memset(st_b_tmp_key, 0x00, sizeof(st_b_tmp_key));

    for (unsigned int b = 0; b < 2; b++) {
      // ::19
      aes_256_sbox_affine_verifier(st_b_key[b], st_dash_key, delta, b, params);
      // ::20
      aes_256_shiftrows_verifier(st_b_tmp_key[b], st_b_key[b], params);
      memcpy(st_b_key[b], st_b_tmp_key[b], sizeof(bf256_t) * Nstbytes);
      // ::21
      aes_256_mix_columns_verifier(st_b_tmp_key[b], st_b_key[b], b, params);
      memcpy(st_b_key[b], st_b_tmp_key[b], sizeof(bf256_t) * Nstbytes);
      // ::22
      if (b == 0) {
        aes_256_add_round_key_bytes_verifier(st_b_key[b], st_b_key[b], k_0_key, delta, true,
                                             params);
      } else {
        aes_256_add_round_key_bytes_verifier(st_b_key[b], st_b_key[b], k_1_key, delta, false,
                                             params);
      }
    }
    // ::23-24
    bf256_t* s_tilde_key = BF256_ALLOC(Nstbits);
    if (r == R / 2 - 1) {
      // ::25
      aes_256_add_round_key_verifier(s_tilde_key, owf_out_key, rkeys_key + R * Nstbits, params);
    } else {
      // ::27-28
      unsigned int idx = 0;
      for (unsigned int i = ((Nstbits / 2) + (Nstbits / 2) * 3 * r);
           i < ((Nstbits / 2 * 3 * r) + (Nstbits / 2 * 3)); i++) {
        s_tilde_key[idx] = w_key[i];
        idx++;
      }
    }
    // ::29
    bf256_t* s_dash_dash_key = BF256_ALLOC(Nstbits);
    aes_256_inverse_shiftrows_verifier(s_dash_dash_key, s_tilde_key, params);
    // ::30
    bf256_t* s_state_key = BF256_ALLOC(Nstbits);
    aes_256_inverse_affine_verifier(s_state_key, s_dash_dash_key, delta, params);

    // ::31
    for (unsigned int byte_i = 0; byte_i < Nstbytes; byte_i++) {
      bf256_t s_key;
      bf256_t s_sq_key;
      // ::32
      s_key = bf256_byte_combine(s_state_key + 8 * byte_i);
      // ::33
      s_sq_key = bf256_byte_combine_sq(s_state_key + 8 * byte_i);

      // ::36
      // compute <s^sq>^1 * <st_{0,i}>^2 - <s>^1
      //    s_sq * st_{0,i} + delta^2 * s
      //
      z_key[(3 * r + 1) * Nstbytes + 2 * byte_i] = bf256_add(
          bf256_mul(s_sq_key, st_b_key[0][byte_i]), bf256_mul(delta, bf256_mul(delta, s_key)));
      // ::37
      // compute <s>^1 * <st_{1,i}>^2 - <st_{0,i}>^2
      //    s * st_{1,i} + delta * st_{0,i}
      //
      z_key[(3 * r + 1) * Nstbytes + 2 * byte_i + 1] =
          bf256_add(bf256_mul(s_key, st_b_key[1][byte_i]), bf256_mul(delta, st_b_key[0][byte_i]));
    }
    if (r != (R / 2) - 1) {
      bf256_t* tmp_state_key = BF256_ALLOC(Nstbits);
      aes_256_bitwise_mix_column_verifier(tmp_state_key, s_tilde_key, params);
      aes_256_add_round_key_verifier(state_bits_key, tmp_state_key,
                                     rkeys_key + (2 * r + 2) * Nstbits, params);
      faest_aligned_free(tmp_state_key);
    }
    faest_aligned_free(s_tilde_key);
    faest_aligned_free(s_dash_dash_key);
    faest_aligned_free(s_state_key);
  }
  faest_aligned_free(st_dash_key);
  faest_aligned_free(state_conj_key);
  faest_aligned_free(state_bits_key);
}

// OWF CONSTRAINTS
static void aes_128_constraints_prover(zk_hash_128_3_ctx* hasher, const uint8_t* w,
                                       const bf128_t* w_tag, const uint8_t* owf_in,
                                       const uint8_t* owf_out, const faest_paramset_t* params) {
  const unsigned int lambda    = FAEST_128_LAMBDA;
  const unsigned int R         = params->R;
  const unsigned int Lke       = params->Lke;
  const unsigned int Lenc      = params->Lenc;
  const unsigned int Nk        = lambda / 32;
  const unsigned int Nst       = params->Nst;
  const unsigned int blocksize = 32 * params->Nst;
  const unsigned int beta      = (lambda + blocksize - 1) / blocksize;
  // ::1-3 owf_in, owf_out, z and z_tag

  // ::4-5
  zk_hash_128_3_raise_and_update(
      hasher, bf128_mul(w_tag[0], w_tag[1]),
      bf128_add(bf128_mul_bit(w_tag[0], w[1]), bf128_mul_bit(w_tag[1], w[0])));

  // ::7-8
  uint8_t* in        = malloc(blocksize);
  bf128_t* in_tag    = BF128_ALLOC(blocksize);
  uint8_t* out       = malloc(blocksize);
  bf128_t* out_tag   = BF128_ALLOC(beta * blocksize);
  uint8_t* rkeys     = malloc((R + 1) * blocksize * sizeof(uint8_t));
  bf128_t* rkeys_tag = BF128_ALLOC((R + 1) * blocksize);

  if (faest_is_em(params)) {
    aes_round_keys_t round_keys;
    expand_key(&round_keys, owf_in, Nk, Nk, R);

    unsigned int idx = 0;
    for (unsigned int r = 0; r < R + 1; r++) {
      for (unsigned int n = 0; n < Nst; n++) {
        for (unsigned int i = 0; i < 4; i++) {
          for (unsigned int j = 0; j < 8; j++) {
            rkeys[8 * idx + j]     = get_bit(round_keys.round_keys[r][n][i], j);
            rkeys_tag[8 * idx + j] = bf128_zero();
          }
          idx++;
        }
      }
    }
    // ::10
    for (unsigned int i = 0; i < blocksize; i++) {
      in[i]     = w[i];
      in_tag[i] = w_tag[i];
    }
    // ::11
    for (unsigned int i = 0; i < blocksize; i++) {
      out[i]     = w[i] ^ ((owf_out[i / 8] >> (i % 8)) & 1);
      out_tag[i] = w_tag[i];
    }
  } else {
    // jump to ::13 for AES
    for (unsigned int i = 0; i < blocksize; i++) {
      in[i] = (owf_in[i / 8] >> (i % 8)) & 1;
    }
    constant_to_vole_128_prover(in_tag, blocksize);

    // ::14
    for (unsigned int i = 0; i < beta * blocksize; i++) {
      out[i] = (owf_out[i / 8] >> (i % 8)) & 1;
    }
    constant_to_vole_128_prover(out_tag, blocksize);

    // ::15 skiped as B = 1
    // ::16
    aes_128_expkey_constraints_prover(hasher, rkeys, rkeys_tag, w, w_tag, params);
  }

  uint8_t* w_tilde     = malloc(Lenc * sizeof(uint8_t));
  bf128_t* w_tilde_tag = BF128_ALLOC(Lenc);

  // ::18-20
  for (unsigned int b = 0; b < beta; b++) {
    for (unsigned int i = 0; i < Lenc; i++) {
      w_tilde[i]     = w[Lke + b * Lenc + i];
      w_tilde_tag[i] = w_tag[Lke + b * Lenc + i];
    }

    if (b == 1) {
      in[0]     = in[0] ^ 0x01;
      in_tag[0] = bf128_add(in_tag[0], bf128_one());
      out_tag += blocksize;
    }

    aes_128_enc_constraints_prover(hasher, in, in_tag, out, out_tag, w_tilde, w_tilde_tag, rkeys,
                                   rkeys_tag, params);
  }

  free(in);
  faest_aligned_free(in_tag);
  free(out);
  faest_aligned_free(out_tag);
  free(rkeys);
  faest_aligned_free(rkeys_tag);

  free(w_tilde);
  faest_aligned_free(w_tilde_tag);
}

static void aes_192_constraints_prover(bf192_t* z_deg0, bf192_t* z_deg1, bf192_t* z_deg2,
                                       const uint8_t* w, const bf192_t* w_tag,
                                       const uint8_t* owf_in, const uint8_t* owf_out,
                                       const faest_paramset_t* params) {
  const unsigned int lambda              = FAEST_192_LAMBDA;
  const unsigned int R                   = params->R;
  const unsigned int Ske                 = params->Ske;
  const unsigned int Lke                 = params->Lke;
  const unsigned int Lenc                = params->Lenc;
  const unsigned int Senc                = params->Senc;
  const unsigned int Nk                  = lambda / 32;
  const unsigned int Nst                 = params->Nst;
  const unsigned int num_enc_constraints = 3 * Senc / 2;
  const unsigned int blocksize           = 32 * params->Nst;
  const unsigned int beta                = (lambda + blocksize - 1) / blocksize;
  // ::1-3 owf_in, owf_out, z and z_tag

  // ::4-5
  z_deg0[0] = bf192_zero();
  z_deg1[0] = bf192_mul(w_tag[0], w_tag[1]);
  z_deg2[0] = bf192_add(bf192_mul_bit(w_tag[0], w[1]), bf192_mul_bit(w_tag[1], w[0]));
  // ::7-8
  uint8_t* in        = malloc(blocksize);
  bf192_t* in_tag    = BF192_ALLOC(blocksize);
  uint8_t* out       = malloc(beta * blocksize);
  bf192_t* out_tag   = BF192_ALLOC(beta * blocksize);
  uint8_t* rkeys     = malloc((R + 1) * blocksize * sizeof(uint8_t));
  bf192_t* rkeys_tag = BF192_ALLOC((R + 1) * blocksize);

  if (faest_is_em(params)) {
    aes_round_keys_t round_keys;
    expand_key(&round_keys, owf_in, Nk, Nk, R);

    unsigned int idx = 0;
    for (unsigned int r = 0; r < R + 1; r++) {
      for (unsigned int n = 0; n < Nst; n++) {
        for (unsigned int i = 0; i < 4; i++) {
          for (unsigned int j = 0; j < 8; j++) {
            rkeys[8 * idx + j]     = get_bit(round_keys.round_keys[r][n][i], j);
            rkeys_tag[8 * idx + j] = bf192_zero();
          }
          idx++;
        }
      }
    }
    // ::10
    for (unsigned int i = 0; i < blocksize; i++) {
      in[i]     = w[i];
      in_tag[i] = w_tag[i];
    }
    // ::11
    for (unsigned int i = 0; i < blocksize; i++) {
      out[i]     = w[i] ^ ((owf_out[i / 8] >> (i % 8)) & 1);
      out_tag[i] = w_tag[i];
    }
  } else {
    // jump to ::13 for AES
    for (unsigned int i = 0; i < blocksize; i++) {
      in[i] = (owf_in[i / 8] >> (i % 8)) & 1;
    }
    constant_to_vole_192_prover(in_tag, blocksize);

    // ::14
    for (unsigned int i = 0; i < beta * blocksize; i++) {
      out[i] = (owf_out[i / 8] >> (i % 8)) & 1;
    }
    constant_to_vole_192_prover(out_tag, beta * blocksize);

    // ::15 skiped as B = 1
    // ::16
    bf192_t* z_tilde_deg0_tag = BF192_ALLOC(2 * Ske);
    bf192_t* z_tilde_deg1_val = BF192_ALLOC(2 * Ske);
    aes_192_expkey_constraints_prover(z_tilde_deg0_tag, z_tilde_deg1_val, rkeys, rkeys_tag, w,
                                      w_tag, params); // w is bit per uint8

    // ::17 raise degree
    for (unsigned int i = 0; i < 2 * Ske; i++) {
      z_deg0[1 + i] = bf192_zero();
      z_deg1[1 + i] = z_tilde_deg0_tag[i];
      z_deg2[1 + i] = z_tilde_deg1_val[i];
      // aes_192_deg2to3_prover(z_deg1 + 1, z_deg2 + 1, z_tilde_deg0_tag[i], z_tilde_deg1_val[i]);
    }

    faest_aligned_free(z_tilde_deg0_tag);
    faest_aligned_free(z_tilde_deg1_val);
  }
  uint8_t* w_tilde     = malloc(Lenc * sizeof(uint8_t));
  bf192_t* w_tilde_tag = BF192_ALLOC(Lenc);

  bf192_t* z_tilde_deg0 = BF192_ALLOC(num_enc_constraints);
  bf192_t* z_tilde_deg1 = BF192_ALLOC(num_enc_constraints);
  bf192_t* z_tilde_deg2 = BF192_ALLOC(num_enc_constraints);

  // ::18-20
  for (unsigned int b = 0; b < beta; b++) {
    for (unsigned int i = 0; i < Lenc; i++) {
      w_tilde[i]     = w[Lke + b * Lenc + i];
      w_tilde_tag[i] = w_tag[Lke + b * Lenc + i];
    }

    memset(z_tilde_deg0, 0, num_enc_constraints * sizeof(bf192_t));
    memset(z_tilde_deg1, 0, num_enc_constraints * sizeof(bf192_t));
    memset(z_tilde_deg2, 0, num_enc_constraints * sizeof(bf192_t));

    if (b == 1) {
      in[0] = in[0] ^ 0x01;
    }

    aes_192_enc_constraints_prover(z_tilde_deg0, z_tilde_deg1, z_tilde_deg2, in, in_tag,
                                   out + b * blocksize, out_tag + b * blocksize, w_tilde,
                                   w_tilde_tag, rkeys, rkeys_tag, params);

    // :22
    for (unsigned int i = 0; i < num_enc_constraints; i++) {
      z_deg0[1 + (2 * Ske) + b * num_enc_constraints + i] = z_tilde_deg0[i];
      z_deg1[1 + (2 * Ske) + b * num_enc_constraints + i] = z_tilde_deg1[i];
      z_deg2[1 + (2 * Ske) + b * num_enc_constraints + i] = z_tilde_deg2[i];
    }
  }
  free(in);
  faest_aligned_free(in_tag);
  free(out);
  faest_aligned_free(out_tag);
  free(rkeys);
  faest_aligned_free(rkeys_tag);

  free(w_tilde);
  faest_aligned_free(w_tilde_tag);

  faest_aligned_free(z_tilde_deg0);
  faest_aligned_free(z_tilde_deg1);
  faest_aligned_free(z_tilde_deg2);
}

static void aes_256_constraints_prover(bf256_t* z_deg0, bf256_t* z_deg1, bf256_t* z_deg2,
                                       const uint8_t* w, const bf256_t* w_tag,
                                       const uint8_t* owf_in, const uint8_t* owf_out,
                                       const faest_paramset_t* params) {
  const unsigned int lambda              = FAEST_256_LAMBDA;
  const unsigned int R                   = params->R;
  const unsigned int Ske                 = params->Ske;
  const unsigned int Lke                 = params->Lke;
  const unsigned int Lenc                = params->Lenc;
  const unsigned int Senc                = params->Senc;
  const unsigned int Nk                  = lambda / 32;
  const unsigned int Nst                 = params->Nst;
  const unsigned int num_enc_constraints = 3 * Senc / 2;
  const unsigned int blocksize           = 32 * params->Nst;
  const unsigned int beta                = (lambda + blocksize - 1) / blocksize;
  // ::1-3 owf_in, owf_out, z and z_tag

  // ::4-5
  z_deg0[0] = bf256_zero();
  z_deg1[0] = bf256_mul(w_tag[0], w_tag[1]);
  z_deg2[0] = bf256_add(bf256_mul_bit(w_tag[0], w[1]), bf256_mul_bit(w_tag[1], w[0]));
  // ::7-8
  uint8_t* in        = malloc(blocksize);
  bf256_t* in_tag    = BF256_ALLOC(blocksize);
  uint8_t* out       = malloc(beta * blocksize);
  bf256_t* out_tag   = BF256_ALLOC(beta * blocksize);
  uint8_t* rkeys     = malloc((R + 1) * blocksize * sizeof(uint8_t));
  bf256_t* rkeys_tag = BF256_ALLOC((R + 1) * blocksize);

  if (faest_is_em(params)) {
    aes_round_keys_t round_keys;
    expand_key(&round_keys, owf_in, Nk, Nk, R);

    unsigned int idx = 0;
    for (unsigned int r = 0; r < R + 1; r++) {
      for (unsigned int n = 0; n < Nst; n++) {
        for (unsigned int i = 0; i < 4; i++) {
          for (unsigned int j = 0; j < 8; j++) {
            rkeys[8 * idx + j]     = get_bit(round_keys.round_keys[r][n][i], j);
            rkeys_tag[8 * idx + j] = bf256_zero();
          }
          idx++;
        }
      }
    }
    // ::10
    for (unsigned int i = 0; i < blocksize; i++) {
      in[i]     = w[i];
      in_tag[i] = w_tag[i];
    }
    // ::11
    for (unsigned int i = 0; i < blocksize; i++) {
      out[i]     = w[i] ^ ((owf_out[i / 8] >> (i % 8)) & 1);
      out_tag[i] = w_tag[i];
    }
  } else {
    // jump to ::13 for AES
    for (unsigned int i = 0; i < blocksize; i++) {
      in[i] = (owf_in[i / 8] >> (i % 8)) & 1;
    }
    constant_to_vole_256_prover(in_tag, blocksize);

    // ::14
    for (unsigned int i = 0; i < beta * blocksize; i++) {
      out[i] = (owf_out[i / 8] >> (i % 8)) & 1;
    }
    constant_to_vole_256_prover(out_tag, beta * blocksize);

    // ::15 skiped as B = 1
    // ::16
    bf256_t* z_tilde_deg0_tag = BF256_ALLOC(2 * Ske);
    bf256_t* z_tilde_deg1_val = BF256_ALLOC(2 * Ske);
    aes_256_expkey_constraints_prover(z_tilde_deg0_tag, z_tilde_deg1_val, rkeys, rkeys_tag, w,
                                      w_tag, params); // w is bit per uint8

    // ::17 raise degree
    for (unsigned int i = 0; i < 2 * Ske; i++) {
      z_deg0[1 + i] = bf256_zero();
      z_deg1[1 + i] = z_tilde_deg0_tag[i];
      z_deg2[1 + i] = z_tilde_deg1_val[i];
      // aes_256_deg2to3_prover(z_deg1 + 1, z_deg2 + 1, z_tilde_deg0_tag[i], z_tilde_deg1_val[i]);
    }

    faest_aligned_free(z_tilde_deg0_tag);
    faest_aligned_free(z_tilde_deg1_val);
  }
  uint8_t* w_tilde     = malloc(Lenc * sizeof(uint8_t));
  bf256_t* w_tilde_tag = BF256_ALLOC(Lenc);

  bf256_t* z_tilde_deg0 = BF256_ALLOC(num_enc_constraints);
  bf256_t* z_tilde_deg1 = BF256_ALLOC(num_enc_constraints);
  bf256_t* z_tilde_deg2 = BF256_ALLOC(num_enc_constraints);

  // ::18-20
  for (unsigned int b = 0; b < beta; b++) {
    for (unsigned int i = 0; i < Lenc; i++) {
      w_tilde[i]     = w[Lke + b * Lenc + i];
      w_tilde_tag[i] = w_tag[Lke + b * Lenc + i];
    }

    memset(z_tilde_deg0, 0, num_enc_constraints * sizeof(bf256_t));
    memset(z_tilde_deg1, 0, num_enc_constraints * sizeof(bf256_t));
    memset(z_tilde_deg2, 0, num_enc_constraints * sizeof(bf256_t));

    if (b == 1) {
      in[0] = in[0] ^ 0x01;
    }

    aes_256_enc_constraints_prover(z_tilde_deg0, z_tilde_deg1, z_tilde_deg2, in, in_tag,
                                   out + b * blocksize, out_tag + b * blocksize, w_tilde,
                                   w_tilde_tag, rkeys, rkeys_tag, params);

    // :22
    for (unsigned int i = 0; i < num_enc_constraints; i++) {
      z_deg0[1 + (2 * Ske) + b * num_enc_constraints + i] = z_tilde_deg0[i];
      z_deg1[1 + (2 * Ske) + b * num_enc_constraints + i] = z_tilde_deg1[i];
      z_deg2[1 + (2 * Ske) + b * num_enc_constraints + i] = z_tilde_deg2[i];
    }
  }
  free(in);
  faest_aligned_free(in_tag);
  free(out);
  faest_aligned_free(out_tag);
  free(rkeys);
  faest_aligned_free(rkeys_tag);

  free(w_tilde);
  faest_aligned_free(w_tilde_tag);

  faest_aligned_free(z_tilde_deg0);
  faest_aligned_free(z_tilde_deg1);
  faest_aligned_free(z_tilde_deg2);
}

// OWF CONSTRAINTS VERIFIER
static void aes_128_constraints_verifier(zk_hash_128_ctx* hasher, const bf128_t* w_key,
                                         const uint8_t* owf_in, const uint8_t* owf_out,
                                         bf128_t delta, const faest_paramset_t* params) {
  const unsigned int lambda    = FAEST_128_LAMBDA;
  const unsigned int R         = params->R;
  const unsigned int Lke       = params->Lke;
  const unsigned int Lenc      = params->Lenc;
  const unsigned int Nk        = lambda / 32;
  const unsigned int Nst       = params->Nst;
  const unsigned int blocksize = 32 * Nst;
  const unsigned int beta      = (lambda + blocksize - 1) / blocksize;
  // ::1-3 owf_in, owf_out, z and z_tag

  // ::4-5
  zk_hash_128_update(hasher, bf128_mul(delta, bf128_mul(w_key[0], w_key[1])));

  // ::7-8
  bf128_t* rkeys_key = BF128_ALLOC((R + 1) * blocksize);
  bf128_t* in_key    = BF128_ALLOC(blocksize);
  bf128_t* out_key   = BF128_ALLOC(beta * blocksize);

  if (faest_is_em(params)) {
    aes_round_keys_t round_keys;
    expand_key(&round_keys, owf_in, Nk, Nk, R);

    unsigned int idx = 0;
    for (unsigned int r = 0; r < R + 1; r++) {
      for (unsigned int n = 0; n < Nst; n++) {
        for (unsigned int i = 0; i < 4; i++) {
          uint8_t rk_byte = round_keys.round_keys[r][n][i];
          for (unsigned int j = 0; j < 8; j++) {
            rkeys_key[8 * idx + j] = bf128_mul_bit(delta, get_bit(rk_byte, j));
          }
          idx++;
        }
      }
    }
    // ::10-11
    for (unsigned int i = 0; i < blocksize; i++) {
      in_key[i]  = w_key[i];
      out_key[i] = bf128_add(w_key[i], bf128_mul_bit(delta, ptr_get_bit(owf_out, i)));
    }
  } else {
    // jump to ::13 for AES
    constant_to_vole_128_verifier(in_key, owf_in, delta, blocksize);

    // ::14-15
    // if beta=2, load both public key blocks
    constant_to_vole_128_verifier(out_key, owf_out, delta, beta * blocksize);

    // ::16
    aes_128_expkey_constraints_verifier(hasher, rkeys_key, w_key, delta, params);
  }
  // ::18-20
  bf128_t* w_tilde_key = BF128_ALLOC(Lenc);

  for (unsigned int b = 0; b < beta; b++) {
    for (unsigned int i = 0; i < Lenc; i++) {
      w_tilde_key[i] = w_key[Lke + b * Lenc + i];
    }
    // ::21
    if (b == 1) {
      in_key[0] = bf128_add(in_key[0], delta); // adding one
      out_key += blocksize;
    }
    aes_128_enc_constraints_verifier(hasher, in_key, out_key, w_tilde_key, rkeys_key, delta,
                                     params);
  }

  faest_aligned_free(rkeys_key);
  faest_aligned_free(in_key);
  faest_aligned_free(out_key);
  faest_aligned_free(w_tilde_key);
}

static void aes_192_constraints_verifier(bf192_t* z_key, const bf192_t* w_key,
                                         const uint8_t* owf_in, const uint8_t* owf_out,
                                         bf192_t delta, const faest_paramset_t* params) {
  const unsigned int lambda              = FAEST_192_LAMBDA;
  const unsigned int R                   = params->R;
  const unsigned int Lke                 = params->Lke;
  const unsigned int Lenc                = params->Lenc;
  const unsigned int Senc                = params->Senc;
  const unsigned int Ske                 = params->Ske;
  const unsigned int Nk                  = lambda / 32;
  const unsigned int Nst                 = params->Nst;
  const unsigned int num_enc_constraints = 3 * Senc / 2;
  const unsigned int num_ks_constraints  = 2 * Ske;
  const unsigned int blocksize           = 32 * Nst;
  const unsigned int beta                = (lambda + blocksize - 1) / blocksize;
  // ::1-3 owf_in, owf_out, z and z_tag

  // ::4-5
  z_key[0] = bf192_mul(delta, bf192_mul(w_key[0], w_key[1]));

  // ::7-8
  bf192_t* rkeys_key = BF192_ALLOC((R + 1) * blocksize);
  bf192_t* in_key    = BF192_ALLOC(blocksize);
  bf192_t* out_key   = BF192_ALLOC(beta * blocksize);

  if (faest_is_em(params)) {
    aes_round_keys_t round_keys;
    expand_key(&round_keys, owf_in, Nk, Nk, R);

    unsigned int idx = 0;
    for (unsigned int r = 0; r < R + 1; r++) {
      for (unsigned int n = 0; n < Nst; n++) {
        for (unsigned int i = 0; i < 4; i++) {
          uint8_t rk_byte = round_keys.round_keys[r][n][i];
          for (unsigned int j = 0; j < 8; j++) {
            rkeys_key[8 * idx + j] = bf192_mul_bit(delta, get_bit(rk_byte, j));
          }
          idx++;
        }
      }
    }
    // ::10-11
    for (unsigned int i = 0; i < blocksize; i++) {
      in_key[i]  = w_key[i];
      out_key[i] = bf192_add(w_key[i], bf192_mul_bit(delta, ptr_get_bit(owf_out, i)));
    }
  } else {
    // jump to ::13 for AES
    constant_to_vole_192_verifier(in_key, owf_in, delta, blocksize);

    // ::14-15
    // if beta=2, load both public key blocks
    constant_to_vole_192_verifier(out_key, owf_out, delta, beta * blocksize);

    // ::16
    aes_192_expkey_constraints_verifier(z_key + 1, rkeys_key, w_key, delta, params);

    // ::17 raise degree
    for (unsigned int i = 0; i < num_ks_constraints; i++) {
      z_key[1 + i] = bf192_mul(delta, z_key[1 + i]);
    }
  }
  // ::18-20
  bf192_t* w_tilde_key = BF192_ALLOC(Lenc);

  for (unsigned int b = 0; b < beta; b++) {
    for (unsigned int i = 0; i < Lenc; i++) {
      w_tilde_key[i] = w_key[Lke + b * Lenc + i];
    }
    // ::21
    if (b == 1) {
      in_key[0] = bf192_add(in_key[0], delta); // adding one
    }
    aes_192_enc_constraints_verifier(z_key + 1 + num_ks_constraints + b * num_enc_constraints,
                                     in_key, out_key + b * blocksize, w_tilde_key, rkeys_key, delta,
                                     params);
  }

  faest_aligned_free(rkeys_key);
  faest_aligned_free(in_key);
  faest_aligned_free(out_key);
  faest_aligned_free(w_tilde_key);
}

static void aes_256_constraints_verifier(bf256_t* z_key, const bf256_t* w_key,
                                         const uint8_t* owf_in, const uint8_t* owf_out,
                                         bf256_t delta, const faest_paramset_t* params) {
  const unsigned int lambda              = FAEST_256_LAMBDA;
  const unsigned int R                   = params->R;
  const unsigned int Lke                 = params->Lke;
  const unsigned int Lenc                = params->Lenc;
  const unsigned int Senc                = params->Senc;
  const unsigned int Ske                 = params->Ske;
  const unsigned int Nk                  = lambda / 32;
  const unsigned int Nst                 = params->Nst;
  const unsigned int num_enc_constraints = 3 * Senc / 2;
  const unsigned int num_ks_constraints  = 2 * Ske;
  const unsigned int blocksize           = 32 * Nst;
  const unsigned int beta                = (lambda + blocksize - 1) / blocksize;
  // ::1-3 owf_in, owf_out, z and z_tag

  // ::4-5
  z_key[0] = bf256_mul(delta, bf256_mul(w_key[0], w_key[1]));

  // ::7-8
  bf256_t* rkeys_key = BF256_ALLOC((R + 1) * blocksize);
  bf256_t* in_key    = BF256_ALLOC(blocksize);
  bf256_t* out_key   = BF256_ALLOC(beta * blocksize);

  if (faest_is_em(params)) {
    aes_round_keys_t round_keys;
    expand_key(&round_keys, owf_in, Nk, Nk, R);

    unsigned int idx = 0;
    for (unsigned int r = 0; r < R + 1; r++) {
      for (unsigned int n = 0; n < Nst; n++) {
        for (unsigned int i = 0; i < 4; i++) {
          uint8_t rk_byte = round_keys.round_keys[r][n][i];
          for (unsigned int j = 0; j < 8; j++) {
            rkeys_key[8 * idx + j] = bf256_mul_bit(delta, get_bit(rk_byte, j));
          }
          idx++;
        }
      }
    }
    // ::10-11
    for (unsigned int i = 0; i < blocksize; i++) {
      in_key[i]  = w_key[i];
      out_key[i] = bf256_add(w_key[i], bf256_mul_bit(delta, ptr_get_bit(owf_out, i)));
    }
  } else {
    // jump to ::13 for AES
    constant_to_vole_256_verifier(in_key, owf_in, delta, blocksize);

    // ::14-15
    // if beta=2, load both public key blocks
    constant_to_vole_256_verifier(out_key, owf_out, delta, beta * blocksize);

    // ::16
    aes_256_expkey_constraints_verifier(z_key + 1, rkeys_key, w_key, delta, params);

    // ::17 raise degree
    for (unsigned int i = 0; i < num_ks_constraints; i++) {
      z_key[1 + i] = bf256_mul(delta, z_key[1 + i]);
    }
  }
  // ::18-20
  bf256_t* w_tilde_key = BF256_ALLOC(Lenc);

  for (unsigned int b = 0; b < beta; b++) {
    for (unsigned int i = 0; i < Lenc; i++) {
      w_tilde_key[i] = w_key[Lke + b * Lenc + i];
    }
    // ::21
    if (b == 1) {
      in_key[0] = bf256_add(in_key[0], delta); // adding one
    }
    aes_256_enc_constraints_verifier(z_key + 1 + num_ks_constraints + b * num_enc_constraints,
                                     in_key, out_key + b * blocksize, w_tilde_key, rkeys_key, delta,
                                     params);
  }

  faest_aligned_free(rkeys_key);
  faest_aligned_free(in_key);
  faest_aligned_free(out_key);
  faest_aligned_free(w_tilde_key);
}

// OWF PROVER
static void aes_128_prover(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde,
                           const uint8_t* w_bits, const uint8_t* u, uint8_t** V,
                           const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* chall_2,
                           const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_128_LAMBDA;
  const unsigned int c      = params->C;
  const unsigned int ell    = params->l;

  // ::1-5
  // V becomes the w_tag: ell + 2*lambda field elements
  bf128_t* w_tag = column_to_row_major_and_shrink_V_128(V, ell); // This is the tag for w

  // ::6-7 embed VOLE masks
  bf128_t bf_u_star_0 = bf128_sum_poly_bits(u);
  bf128_t bf_u_star_1 = bf128_sum_poly_bits(u + lambda / 8);

  // ::8-9
  bf128_t bf_v_star_0 = bf128_sum_poly(w_tag + ell);
  bf128_t bf_v_star_1 = bf128_sum_poly(w_tag + ell + lambda);

  // Step: 13-18
  zk_hash_128_3_ctx hasher;
  zk_hash_128_3_init(&hasher, chall_2);

  aes_128_constraints_prover(&hasher, w_bits, w_tag, owf_in, owf_out, params);

  zk_hash_128_3_finalize(a0_tilde, a1_tilde, a2_tilde, &hasher, bf_v_star_0,
                         bf128_add(bf_u_star_0, bf_v_star_1), bf_u_star_1);

  faest_aligned_free(w_tag);
}

static void aes_192_prover(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde,
                           const uint8_t* w_bits, const uint8_t* u, uint8_t** V,
                           const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* chall_2,
                           const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_192_LAMBDA;
  const unsigned int c      = params->C;
  const unsigned int ell    = params->l;

  // ::1-5
  // V becomes the w_tag
  bf192_t* w_tag = column_to_row_major_and_shrink_V_192(V, ell); // This is the tag for w

  // ::6-7 embed VOLE masks
  bf192_t bf_u_star_0 = bf192_sum_poly_bits(u);
  bf192_t bf_u_star_1 = bf192_sum_poly_bits(u + lambda / 8);

  // ::8-9
  bf192_t bf_v_star_0 = bf192_sum_poly(w_tag + ell);
  bf192_t bf_v_star_1 = bf192_sum_poly(w_tag + ell + lambda);

  // ::10-12
  bf192_t* z0_tag   = BF192_ALLOC(c); // this contains the bf tag
  bf192_t* z1_val   = BF192_ALLOC(c); // this contains the bf val
  bf192_t* z2_gamma = BF192_ALLOC(c); // this contains the bf gamma
  memset(z0_tag, 0, c * sizeof(bf192_t));
  memset(z1_val, 0, c * sizeof(bf192_t));
  memset(z2_gamma, 0, c * sizeof(bf192_t));

  aes_192_constraints_prover(z0_tag, z1_val, z2_gamma, w_bits, w_tag, owf_in, owf_out, params);

  // Step: 13-18
  zk_hash_192_ctx a0_ctx;
  zk_hash_192_ctx a1_ctx;
  zk_hash_192_ctx a2_ctx;
  zk_hash_192_init(&a0_ctx, chall_2);
  zk_hash_192_init(&a1_ctx, chall_2);
  zk_hash_192_init(&a2_ctx, chall_2);

  for (unsigned int i = 0; i < c; i++) {
    zk_hash_192_update(&a0_ctx, z0_tag[i]);
    zk_hash_192_update(&a1_ctx, z1_val[i]);
    zk_hash_192_update(&a2_ctx, z2_gamma[i]);
  }

  faest_aligned_free(z0_tag);
  faest_aligned_free(z1_val);
  faest_aligned_free(z2_gamma);

  zk_hash_192_finalize(a0_tilde, &a0_ctx, bf_v_star_0);
  zk_hash_192_finalize(a1_tilde, &a1_ctx, bf192_add(bf_u_star_0, bf_v_star_1));
  zk_hash_192_finalize(a2_tilde, &a2_ctx, bf_u_star_1);

  faest_aligned_free(w_tag);
}

static void aes_256_prover(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde,
                           const uint8_t* w_bits, const uint8_t* u, uint8_t** V,
                           const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* chall_2,
                           const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_256_LAMBDA;
  const unsigned int c      = params->C;
  const unsigned int ell    = params->l;

  // ::1-5
  // V becomes the w_tag
  bf256_t* w_tag = column_to_row_major_and_shrink_V_256(V, ell); // This is the tag for w

  // ::6-7 embed VOLE masks
  bf256_t bf_u_star_0 = bf256_sum_poly_bits(u);
  bf256_t bf_u_star_1 = bf256_sum_poly_bits(u + lambda / 8);

  // ::8-9
  bf256_t bf_v_star_0 = bf256_sum_poly(w_tag + ell);
  bf256_t bf_v_star_1 = bf256_sum_poly(w_tag + ell + lambda);

  // ::10-12
  bf256_t* z0_tag   = BF256_ALLOC(c); // this contains the bf tag
  bf256_t* z1_val   = BF256_ALLOC(c); // this contains the bf val
  bf256_t* z2_gamma = BF256_ALLOC(c); // this contains the bf gamma
  memset(z0_tag, 0, c * sizeof(bf256_t));
  memset(z1_val, 0, c * sizeof(bf256_t));
  memset(z2_gamma, 0, c * sizeof(bf256_t));

  aes_256_constraints_prover(z0_tag, z1_val, z2_gamma, w_bits, w_tag, owf_in, owf_out, params);

  // Step: 13-18
  zk_hash_256_ctx a0_ctx;
  zk_hash_256_ctx a1_ctx;
  zk_hash_256_ctx a2_ctx;
  zk_hash_256_init(&a0_ctx, chall_2);
  zk_hash_256_init(&a1_ctx, chall_2);
  zk_hash_256_init(&a2_ctx, chall_2);

  for (unsigned int i = 0; i < c; i++) {
    zk_hash_256_update(&a0_ctx, z0_tag[i]);
    zk_hash_256_update(&a1_ctx, z1_val[i]);
    zk_hash_256_update(&a2_ctx, z2_gamma[i]);
  }

  faest_aligned_free(z0_tag);
  faest_aligned_free(z1_val);
  faest_aligned_free(z2_gamma);

  zk_hash_256_finalize(a0_tilde, &a0_ctx, bf_v_star_0);
  zk_hash_256_finalize(a1_tilde, &a1_ctx, bf256_add(bf_u_star_0, bf_v_star_1));
  zk_hash_256_finalize(a2_tilde, &a2_ctx, bf_u_star_1);

  faest_aligned_free(w_tag);
}

// OWF VERIFIER
static void aes_128_verifier(uint8_t* a0_tilde, const uint8_t* d, uint8_t** Q,
                             const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* chall_2,
                             const uint8_t* chall_3, const uint8_t* a1_tilde,
                             const uint8_t* a2_tilde, const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_128_LAMBDA;
  const unsigned int ell    = params->l;

  // ::1
  bf128_t bf_delta    = bf128_load(chall_3);
  bf128_t bf_delta_sq = bf128_mul(bf_delta, bf_delta);

  // ::2-6
  bf128_t* q_key = column_to_row_major_and_shrink_V_128(Q, ell);

  // ::7-9
  bf128_t q_star_0 = bf128_sum_poly(q_key + ell);
  bf128_t q_star_1 = bf128_sum_poly(q_key + ell + lambda);

  // ::10
  bf128_t q_star = bf128_add(q_star_0, bf128_mul(bf_delta, q_star_1));

  // ::13-14
  zk_hash_128_ctx b_ctx;
  zk_hash_128_init(&b_ctx, chall_2);

  bf128_t* w_key = BF128_ALLOC(ell);
  for (unsigned int i = 0; i < ell; i++) {
    w_key[i] = bf128_add(q_key[i], bf128_mul_bit(bf_delta, ptr_get_bit(d, i)));
  }

  // ::11-12
  aes_128_constraints_verifier(&b_ctx, w_key, owf_in, owf_out, bf_delta, params);

  // ::13-14
  uint8_t q_tilde[FAEST_128_LAMBDA / 8];
  zk_hash_128_finalize(q_tilde, &b_ctx, q_star);

  faest_aligned_free(w_key);

  // ::16
  bf128_t tmp1 = bf128_mul(bf128_load(a1_tilde), bf_delta);
  bf128_t tmp2 = bf128_mul(bf128_load(a2_tilde), bf_delta_sq);
  bf128_t tmp3 = bf128_add(tmp1, tmp2);
  bf128_t ret  = bf128_add(bf128_load(q_tilde), tmp3);

  faest_aligned_free(q_key);

  bf128_store(a0_tilde, ret);
}

static void aes_192_verifier(uint8_t* a0_tilde, const uint8_t* d, uint8_t** Q,
                             const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* chall_2,
                             const uint8_t* chall_3, const uint8_t* a1_tilde,
                             const uint8_t* a2_tilde, const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_192_LAMBDA;
  const unsigned int c      = params->C;
  const unsigned int ell    = params->l;

  // ::1
  bf192_t bf_delta    = bf192_load(chall_3);
  bf192_t bf_delta_sq = bf192_mul(bf_delta, bf_delta);

  // ::2-6
  bf192_t* q_key = column_to_row_major_and_shrink_V_192(Q, ell);

  // ::7-9
  bf192_t q_star_0 = bf192_sum_poly(q_key + ell);
  bf192_t q_star_1 = bf192_sum_poly(q_key + ell + lambda);

  // ::10
  bf192_t q_star = bf192_add(q_star_0, bf192_mul(bf_delta, q_star_1));

  // ::11-12
  bf192_t* z2_key = BF192_ALLOC(c);
  bf192_t* w_key  = BF192_ALLOC(ell);
  for (unsigned int i = 0; i < ell; i++) {
    w_key[i] = bf192_add(q_key[i], bf192_mul_bit(bf_delta, ptr_get_bit(d, i)));
  }
  memset(z2_key, 0, c * sizeof(bf192_t));
  aes_192_constraints_verifier(z2_key, w_key, owf_in, owf_out, bf_delta, params);

  // ::13-14
  zk_hash_192_ctx b_ctx;
  zk_hash_192_init(&b_ctx, chall_2);
  for (unsigned int i = 0; i < c; i++) {
    zk_hash_192_update(&b_ctx, z2_key[i]);
  }

  uint8_t q_tilde[FAEST_192_LAMBDA / 8];
  zk_hash_192_finalize(q_tilde, &b_ctx, q_star);

  faest_aligned_free(z2_key);
  faest_aligned_free(w_key);

  // ::16
  bf192_t tmp1 = bf192_mul(bf192_load(a1_tilde), bf_delta);
  bf192_t tmp2 = bf192_mul(bf192_load(a2_tilde), bf_delta_sq);
  bf192_t tmp3 = bf192_add(tmp1, tmp2);
  bf192_t ret  = bf192_add(bf192_load(q_tilde), tmp3);

  faest_aligned_free(q_key);

  bf192_store(a0_tilde, ret);
}

static void aes_256_verifier(uint8_t* a0_tilde, const uint8_t* d, uint8_t** Q,
                             const uint8_t* owf_in, const uint8_t* owf_out, const uint8_t* chall_2,
                             const uint8_t* chall_3, const uint8_t* a1_tilde,
                             const uint8_t* a2_tilde, const faest_paramset_t* params) {
  const unsigned int lambda = FAEST_256_LAMBDA;
  const unsigned int c      = params->C;
  const unsigned int ell    = params->l;

  // ::1
  bf256_t bf_delta    = bf256_load(chall_3);
  bf256_t bf_delta_sq = bf256_mul(bf_delta, bf_delta);

  // ::2-6
  bf256_t* q_key = column_to_row_major_and_shrink_V_256(Q, ell);

  // ::7-9
  bf256_t q_star_0 = bf256_sum_poly(q_key + ell);
  bf256_t q_star_1 = bf256_sum_poly(q_key + ell + lambda);

  // ::10
  bf256_t q_star = bf256_add(q_star_0, bf256_mul(bf_delta, q_star_1));

  // ::11-12
  bf256_t* z2_key = BF256_ALLOC(c);
  bf256_t* w_key  = BF256_ALLOC(ell);
  for (unsigned int i = 0; i < ell; i++) {
    w_key[i] = bf256_add(q_key[i], bf256_mul_bit(bf_delta, ptr_get_bit(d, i)));
  }
  memset(z2_key, 0, c * sizeof(bf256_t));
  aes_256_constraints_verifier(z2_key, w_key, owf_in, owf_out, bf_delta, params);

  // ::13-14
  zk_hash_256_ctx b_ctx;
  zk_hash_256_init(&b_ctx, chall_2);
  for (unsigned int i = 0; i < c; i++) {
    zk_hash_256_update(&b_ctx, z2_key[i]);
  }

  uint8_t q_tilde[FAEST_256_LAMBDA];
  zk_hash_256_finalize(q_tilde, &b_ctx, q_star);

  faest_aligned_free(z2_key);
  faest_aligned_free(w_key);

  // ::16
  bf256_t tmp1 = bf256_mul(bf256_load(a1_tilde), bf_delta);
  bf256_t tmp2 = bf256_mul(bf256_load(a2_tilde), bf_delta_sq);
  bf256_t tmp3 = bf256_add(tmp1, tmp2);
  bf256_t ret  = bf256_add(bf256_load(q_tilde), tmp3);

  faest_aligned_free(q_key);

  bf256_store(a0_tilde, ret);
}

// AES(-EM) OWF dispatchers
void aes_prove(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde, const uint8_t* w_bits,
               const uint8_t* u, uint8_t** V, const uint8_t* owf_in, const uint8_t* owf_out,
               const uint8_t* chall_2, const faest_paramset_t* params) {
  switch (params->lambda) {
  case 256:
    aes_256_prover(a0_tilde, a1_tilde, a2_tilde, w_bits, u, V, owf_in, owf_out, chall_2, params);
    break;
  case 192:
    aes_192_prover(a0_tilde, a1_tilde, a2_tilde, w_bits, u, V, owf_in, owf_out, chall_2, params);
    break;
  default:
    aes_128_prover(a0_tilde, a1_tilde, a2_tilde, w_bits, u, V, owf_in, owf_out, chall_2, params);
  }
}

void aes_verify(uint8_t* a0_tilde, const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                const uint8_t* chall_3, const uint8_t* a1_tilde, const uint8_t* a2_tilde,
                const uint8_t* owf_in, const uint8_t* owf_out, const faest_paramset_t* params) {
  switch (params->lambda) {
  case 256:
    aes_256_verifier(a0_tilde, d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
    break;
  case 192:
    aes_192_verifier(a0_tilde, d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
    break;
  default:
    aes_128_verifier(a0_tilde, d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
  }
}
