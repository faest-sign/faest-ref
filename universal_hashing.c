/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "universal_hashing.h"

#include <assert.h>
#include <string.h>

static bf64_t compute_h1(const uint8_t* t, const uint8_t* x, size_t ell) {
  bf64_t b_t = bf64_load(t);
  bf64_t h1  = bf64_zero();

  const size_t length_64 = ell / 64;
  assert(length_64 * 64 == ell);

  bf64_t running_t = bf64_one();
  for (size_t i = 0; i != length_64; ++i, running_t = bf64_mul(running_t, b_t)) {
    h1 = bf64_add(h1, bf64_mul(running_t, bf64_load(x + (length_64 - 1 - i) * 64 / 8)));
  }

  return h1;
}

void vole_hash_128(uint8_t* h, const uint8_t* r0, const uint8_t* r1, const uint8_t* s,
                   const uint8_t* t, const uint8_t* x, size_t ell) {
  const size_t length_lambda = ell / 128;
  assert(length_lambda * 128 == ell);

  bf128_t b_s       = bf128_load(s);
  bf128_t running_s = bf128_one();
  bf128_t h0        = bf128_zero();
  for (size_t i = 0; i != length_lambda; ++i, running_s = bf128_mul(running_s, b_s)) {
    h0 = bf128_add(h0, bf128_mul(running_s, bf128_load(x + (length_lambda - 1 - i) * 128 / 8)));
  }

  // TODO: check with spec
  uint8_t tmp[128 / 8] = {0};
  bf64_store(tmp, compute_h1(t, x, ell));
  bf128_t h1p = bf128_load(tmp);

  h0  = bf128_add(h0, bf128_mul(bf128_load(r0), h1p));
  h1p = bf128_mul(bf128_load(r1), h1p);

  bf128_store(h, h0);
  bf128_store(tmp, h1p);
  memcpy(h + 128 / 8, tmp, UNIVERSAL_HASH_B);
}

void vole_hash_192(uint8_t* h, const uint8_t* r0, const uint8_t* r1, const uint8_t* s,
                   const uint8_t* t, const uint8_t* x, size_t ell) {
  const size_t length_lambda = ell / 192;
  assert(length_lambda * 192 == ell);

  bf192_t b_s       = bf192_load(s);
  bf192_t running_s = bf192_one();
  bf192_t h0        = bf192_zero();
  for (size_t i = 0; i != length_lambda; ++i, running_s = bf192_mul(running_s, b_s)) {
    h0 = bf192_add(h0, bf192_mul(running_s, bf192_load(x + (length_lambda - 1 - i) * 192 / 8)));
  }

  // TODO: check with spec
  uint8_t tmp[192 / 8] = {0};
  bf64_store(tmp, compute_h1(t, x, ell));
  bf192_t h1p = bf192_load(tmp);

  h0  = bf192_add(h0, bf192_mul(bf192_load(r0), h1p));
  h1p = bf192_mul(bf192_load(r1), h1p);

  bf192_store(h, h0);
  bf192_store(tmp, h1p);
  memcpy(h + 192 / 8, tmp, UNIVERSAL_HASH_B);
}

void vole_hash_256(uint8_t* h, const uint8_t* r0, const uint8_t* r1, const uint8_t* s,
                   const uint8_t* t, const uint8_t* x, size_t ell) {
  const size_t length_lambda = ell / 256;
  assert(length_lambda * 256 == ell);

  bf256_t b_s       = bf256_load(s);
  bf256_t running_s = bf256_one();
  bf256_t h0        = bf256_zero();
  for (size_t i = 0; i != length_lambda; ++i, running_s = bf256_mul(running_s, b_s)) {
    h0 = bf256_add(h0, bf256_mul(running_s, bf256_load(x + (length_lambda - 1 - i) * 256 / 8)));
  }

  // TODO: check with spec
  uint8_t tmp[256 / 8] = {0};
  bf64_store(tmp, compute_h1(t, x, ell));
  bf256_t h1p = bf256_load(tmp);

  h0  = bf256_add(h0, bf256_mul(bf256_load(r0), h1p));
  h1p = bf256_mul(bf256_load(r1), h1p);

  bf256_store(h, h0);
  bf256_store(tmp, h1p);
  memcpy(h + 256 / 8, tmp, UNIVERSAL_HASH_B);
}

void zk_hash_128(uint8_t* h, const uint8_t* r, const uint8_t* s, const uint8_t* t, const bf128_t* x,
                 size_t ell) {
  // TODO: check with spec
  uint8_t tmp[128 / 8] = {0};
  memcpy(tmp + (128 - 64) / 8, t, 64 / 8);

  bf128_t b_s       = bf128_load(s);
  bf128_t b_t       = bf128_load(tmp);
  bf128_t running_s = bf128_one();
  bf128_t running_t = bf128_one();
  bf128_t h0        = bf128_zero();
  bf128_t h1        = bf128_zero();
  for (size_t i = 0; i != ell;
       ++i, running_s = bf128_mul(running_s, b_s), running_s = bf128_mul(running_t, b_t)) {
    h0 = bf128_add(h0, bf128_mul(running_s, x[ell - 1 - i]));
    h1 = bf128_add(h1, bf128_mul(running_t, x[ell - 1 - i]));
  }

  h0 = bf128_add(h0, bf128_mul(bf128_load(r), h1));
  bf128_store(h, h0);
}

void zk_hash_192(uint8_t* h, const uint8_t* r, const uint8_t* s, const uint8_t* t, const bf192_t* x,
                 size_t ell) {
  // TODO: check with spec
  uint8_t tmp[192 / 8] = {0};
  memcpy(tmp + (192 - 64) / 8, t, 64 / 8);

  bf192_t b_s       = bf192_load(s);
  bf192_t b_t       = bf192_load(tmp);
  bf192_t running_s = bf192_one();
  bf192_t running_t = bf192_one();
  bf192_t h0        = bf192_zero();
  bf192_t h1        = bf192_zero();
  for (size_t i = 0; i != ell;
       ++i, running_s = bf192_mul(running_s, b_s), running_s = bf192_mul(running_t, b_t)) {
    h0 = bf192_add(h0, bf192_mul(running_s, x[ell - 1 - i]));
    h1 = bf192_add(h1, bf192_mul(running_t, x[ell - 1 - i]));
  }

  h0 = bf192_add(h0, bf192_mul(bf192_load(r), h1));
  bf192_store(h, h0);
}

void zk_hash_256(uint8_t* h, const uint8_t* r, const uint8_t* s, const uint8_t* t, const bf256_t* x,
                 size_t ell) {
  // TODO: check with spec
  uint8_t tmp[256 / 8] = {0};
  memcpy(tmp + (256 - 64) / 8, t, 64 / 8);

  bf256_t b_s       = bf256_load(s);
  bf256_t b_t       = bf256_load(tmp);
  bf256_t running_s = bf256_one();
  bf256_t running_t = bf256_one();
  bf256_t h0        = bf256_zero();
  bf256_t h1        = bf256_zero();
  for (size_t i = 0; i != ell;
       ++i, running_s = bf256_mul(running_s, b_s), running_s = bf256_mul(running_t, b_t)) {
    h0 = bf256_add(h0, bf256_mul(running_s, x[ell - 1 - i]));
    h1 = bf256_add(h1, bf256_mul(running_t, x[ell - 1 - i]));
  }

  h0 = bf256_add(h0, bf256_mul(bf256_load(r), h1));
  bf256_store(h, h0);
}
