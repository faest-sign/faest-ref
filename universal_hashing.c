/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "instances.h"
#include "universal_hashing.h"
#include "utils.h"

#include <assert.h>
#include <string.h>

static bf64_t compute_h1(const uint8_t* t, const uint8_t* x, unsigned int lambda,
                         unsigned int ell) {
  const bf64_t b_t = bf64_load(t);

  unsigned int lambdaBytes         = lambda / 8;
  const unsigned int length_lambda = (ell + lambda + lambda - 1) / lambda;

  uint8_t tmp[MAX_LAMBDA_BYTES] = {0};
  memcpy(tmp, x + (length_lambda - 1) * lambdaBytes,
         (ell + lambda) % lambda == 0 ? lambdaBytes : ((ell + lambda) % lambda) / 8);

  bf64_t h1        = bf64_zero();
  bf64_t running_t = bf64_one();
  unsigned int i   = 0;
  for (; i < lambdaBytes; i += 8, running_t = bf64_mul(running_t, b_t)) {
    h1 = bf64_add(h1, bf64_mul(running_t, bf64_load(tmp + (lambdaBytes - i - 8))));
  }
  for (; i < length_lambda * lambdaBytes; i += 8, running_t = bf64_mul(running_t, b_t)) {
    h1 = bf64_add(h1, bf64_mul(running_t, bf64_load(x + (length_lambda * lambdaBytes - i - 8))));
  }

  return h1;
}

void vole_hash_128(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell) {
  const uint8_t* r0 = sd;
  const uint8_t* r1 = sd + 1 * sizeof(bf128_t);
  const uint8_t* r2 = sd + 2 * sizeof(bf128_t);
  const uint8_t* r3 = sd + 3 * sizeof(bf128_t);
  const uint8_t* s  = sd + 4 * sizeof(bf128_t);
  const uint8_t* t  = sd + 5 * sizeof(bf128_t);
  const uint8_t* x1 = x + (ell + sizeof(bf128_t) * 8) / 8;

  const unsigned int length_lambda = (ell + 2 * sizeof(bf128_t) * 8 - 1) / (sizeof(bf128_t) * 8);

  uint8_t tmp[sizeof(bf128_t)] = {0};
  memcpy(tmp, x + (length_lambda - 1) * sizeof(bf128_t),
         (ell + sizeof(bf128_t) * 8) % (sizeof(bf128_t) * 8) == 0
             ? sizeof(bf128_t)
             : ((ell + sizeof(bf128_t) * 8) % (sizeof(bf128_t) * 8)) / 8);
  bf128_t h0 = bf128_load(tmp);

  const bf128_t b_s = bf128_load(s);
  bf128_t running_s = b_s;
  for (unsigned int i = 1; i != length_lambda; ++i, running_s = bf128_mul(running_s, b_s)) {
    h0 = bf128_add(h0,
                   bf128_mul(running_s, bf128_load(x + (length_lambda - 1 - i) * sizeof(bf128_t))));
  }

  bf128_t h1p = bf128_from_bf64(compute_h1(t, x, sizeof(bf128_t) * 8, ell));
  bf128_t h2  = bf128_add(bf128_mul(bf128_load(r0), h0), bf128_mul(bf128_load(r1), h1p));
  bf128_t h3  = bf128_add(bf128_mul(bf128_load(r2), h0), bf128_mul(bf128_load(r3), h1p));

  bf128_store(h, h2);
  bf128_store(tmp, h3);
  memcpy(h + sizeof(bf128_t), tmp, UNIVERSAL_HASH_B);
  xor_u8_array(h, x1, h, sizeof(bf128_t) + UNIVERSAL_HASH_B);
}

void vole_hash_192(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell) {
  const uint8_t* r0 = sd;
  const uint8_t* r1 = sd + 1 * sizeof(bf192_t);
  const uint8_t* r2 = sd + 2 * sizeof(bf192_t);
  const uint8_t* r3 = sd + 3 * sizeof(bf192_t);
  const uint8_t* s  = sd + 4 * sizeof(bf192_t);
  const uint8_t* t  = sd + 5 * sizeof(bf192_t);
  const uint8_t* x1 = x + (ell + sizeof(bf192_t) * 8) / 8;

  const unsigned int length_lambda = (ell + 2 * sizeof(bf192_t) * 8 - 1) / (sizeof(bf192_t) * 8);

  uint8_t tmp[sizeof(bf192_t)] = {0};
  memcpy(tmp, x + (length_lambda - 1) * sizeof(bf192_t),
         (ell + sizeof(bf192_t) * 8) % (sizeof(bf192_t) * 8) == 0
             ? sizeof(bf192_t)
             : ((ell + sizeof(bf192_t) * 8) % (sizeof(bf192_t) * 8)) / 8);
  bf192_t h0 = bf192_load(tmp);

  const bf192_t b_s = bf192_load(s);
  bf192_t running_s = b_s;
  for (unsigned int i = 1; i != length_lambda; ++i, running_s = bf192_mul(running_s, b_s)) {
    h0 = bf192_add(h0,
                   bf192_mul(running_s, bf192_load(x + (length_lambda - 1 - i) * sizeof(bf192_t))));
  }

  bf192_t h1p = bf192_from_bf64(compute_h1(t, x, sizeof(bf192_t) * 8, ell));
  bf192_t h2  = bf192_add(bf192_mul(bf192_load(r0), h0), bf192_mul(bf192_load(r1), h1p));
  bf192_t h3  = bf192_add(bf192_mul(bf192_load(r2), h0), bf192_mul(bf192_load(r3), h1p));

  bf192_store(h, h2);
  bf192_store(tmp, h3);
  memcpy(h + sizeof(bf192_t), tmp, UNIVERSAL_HASH_B);
  xor_u8_array(h, x1, h, sizeof(bf192_t) + UNIVERSAL_HASH_B);
}

void vole_hash_256(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell) {
  const uint8_t* r0 = sd;
  const uint8_t* r1 = sd + 1 * sizeof(bf256_t);
  const uint8_t* r2 = sd + 2 * sizeof(bf256_t);
  const uint8_t* r3 = sd + 3 * sizeof(bf256_t);
  const uint8_t* s  = sd + 4 * sizeof(bf256_t);
  const uint8_t* t  = sd + 5 * sizeof(bf256_t);
  const uint8_t* x1 = x + (ell + sizeof(bf256_t) * 8) / 8;

  const unsigned int length_lambda = (ell + 2 * sizeof(bf256_t) * 8 - 1) / (sizeof(bf256_t) * 8);

  uint8_t tmp[sizeof(bf256_t)] = {0};
  memcpy(tmp, x + (length_lambda - 1) * sizeof(bf256_t),
         (ell + sizeof(bf256_t) * 8) % (sizeof(bf256_t) * 8) == 0
             ? sizeof(bf256_t)
             : ((ell + sizeof(bf256_t) * 8) % (sizeof(bf256_t) * 8)) / 8);
  bf256_t h0 = bf256_load(tmp);

  const bf256_t b_s = bf256_load(s);
  bf256_t running_s = b_s;
  for (unsigned int i = 1; i != length_lambda; ++i, running_s = bf256_mul(running_s, b_s)) {
    h0 = bf256_add(h0,
                   bf256_mul(running_s, bf256_load(x + (length_lambda - 1 - i) * sizeof(bf256_t))));
  }

  bf256_t h1p = bf256_from_bf64(compute_h1(t, x, sizeof(bf256_t) * 8, ell));
  bf256_t h2  = bf256_add(bf256_mul(bf256_load(r0), h0), bf256_mul(bf256_load(r1), h1p));
  bf256_t h3  = bf256_add(bf256_mul(bf256_load(r2), h0), bf256_mul(bf256_load(r3), h1p));

  bf256_store(h, h2);
  bf256_store(tmp, h3);
  memcpy(h + sizeof(bf256_t), tmp, UNIVERSAL_HASH_B);
  xor_u8_array(h, x1, h, sizeof(bf256_t) + UNIVERSAL_HASH_B);
}

void vole_hash(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell, uint32_t lambda) {
  switch (lambda) {
  case 256:
    vole_hash_256(h, sd, x, ell);
    break;
  case 192:
    vole_hash_192(h, sd, x, ell);
    break;
  default:
    vole_hash_128(h, sd, x, ell);
    break;
  }
}

void zk_hash_128(uint8_t* h, const uint8_t* sd, const bf128_t* x, unsigned int ell) {
  const uint8_t* r0 = sd;
  const uint8_t* r1 = sd + 128 / 8;
  const uint8_t* s  = sd + 2 * 128 / 8;
  const uint8_t* t  = sd + 3 * 128 / 8;
  const bf128_t* x1 = x + ell;

  bf128_t b_s       = bf128_load(s);
  bf128_t b_t       = bf128_from_bf64(bf64_load(t));
  bf128_t running_s = bf128_one();
  bf128_t running_t = bf128_one();
  bf128_t h0        = bf128_zero();
  bf128_t h1        = bf128_zero();
  for (unsigned int i = 0; i != ell;
       ++i, running_s = bf128_mul(running_s, b_s), running_t = bf128_mul(running_t, b_t)) {
    h0 = bf128_add(h0, bf128_mul(running_s, x[ell - 1 - i]));
    h1 = bf128_add(h1, bf128_mul(running_t, x[ell - 1 - i]));
  }

  h0 = bf128_add(bf128_add(bf128_mul(bf128_load(r0), h0), bf128_mul(bf128_load(r1), h1)), *x1);
  bf128_store(h, h0);
}

void zk_hash_192(uint8_t* h, const uint8_t* sd, const bf192_t* x, unsigned int ell) {
  const uint8_t* r0 = sd;
  const uint8_t* r1 = sd + 192 / 8;
  const uint8_t* s  = sd + 2 * 192 / 8;
  const uint8_t* t  = sd + 3 * 192 / 8;
  const bf192_t* x1 = x + ell;

  bf192_t b_s       = bf192_load(s);
  bf192_t b_t       = bf192_from_bf64(bf64_load(t));
  bf192_t running_s = bf192_one();
  bf192_t running_t = bf192_one();
  bf192_t h0        = bf192_zero();
  bf192_t h1        = bf192_zero();
  for (unsigned int i = 0; i != ell;
       ++i, running_s = bf192_mul(running_s, b_s), running_t = bf192_mul(running_t, b_t)) {
    h0 = bf192_add(h0, bf192_mul(running_s, x[ell - 1 - i]));
    h1 = bf192_add(h1, bf192_mul(running_t, x[ell - 1 - i]));
  }

  h0 = bf192_add(bf192_add(bf192_mul(bf192_load(r0), h0), bf192_mul(bf192_load(r1), h1)), *x1);
  bf192_store(h, h0);
}

void zk_hash_256(uint8_t* h, const uint8_t* sd, const bf256_t* x, unsigned int ell) {
  const uint8_t* r0 = sd;
  const uint8_t* r1 = sd + 256 / 8;
  const uint8_t* s  = sd + 2 * 256 / 8;
  const uint8_t* t  = sd + 3 * 256 / 8;
  const bf256_t* x1 = x + ell;

  bf256_t b_s       = bf256_load(s);
  bf256_t b_t       = bf256_from_bf64(bf64_load(t));
  bf256_t running_s = bf256_one();
  bf256_t running_t = bf256_one();
  bf256_t h0        = bf256_zero();
  bf256_t h1        = bf256_zero();
  for (unsigned int i = 0; i != ell;
       ++i, running_s = bf256_mul(running_s, b_s), running_t = bf256_mul(running_t, b_t)) {
    h0 = bf256_add(h0, bf256_mul(running_s, x[ell - 1 - i]));
    h1 = bf256_add(h1, bf256_mul(running_t, x[ell - 1 - i]));
  }

  h0 = bf256_add(bf256_add(bf256_mul(bf256_load(r0), h0), bf256_mul(bf256_load(r1), h1)), *x1);
  bf256_store(h, h0);
}
