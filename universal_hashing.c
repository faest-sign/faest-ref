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

  unsigned int lambda_bytes        = lambda / 8;
  const unsigned int length_lambda = (ell + 3 * lambda - 1) / lambda;

  uint8_t tmp[MAX_LAMBDA_BYTES] = {0};
  memcpy(tmp, x + (length_lambda - 1) * lambda_bytes,
         (ell + lambda) % lambda == 0 ? lambda_bytes : ((ell + lambda) % lambda) / 8);

  bf64_t h1        = bf64_zero();
  bf64_t running_t = bf64_one();
  unsigned int i   = 0;
  for (; i < lambda_bytes; i += 8, running_t = bf64_mul(running_t, b_t)) {
    h1 = bf64_add(h1, bf64_mul(running_t, bf64_load(tmp + (lambda_bytes - i - 8))));
  }
  for (; i < length_lambda * lambda_bytes; i += 8, running_t = bf64_mul(running_t, b_t)) {
    h1 = bf64_add(h1, bf64_mul(running_t, bf64_load(x + (length_lambda * lambda_bytes - i - 8))));
  }

  return h1;
}

void vole_hash_128(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell) {
  const uint8_t* r0 = sd;
  const uint8_t* r1 = sd + 1 * BF128_NUM_BYTES;
  const uint8_t* r2 = sd + 2 * BF128_NUM_BYTES;
  const uint8_t* r3 = sd + 3 * BF128_NUM_BYTES;
  const uint8_t* s  = sd + 4 * BF128_NUM_BYTES;
  const uint8_t* t  = sd + 5 * BF128_NUM_BYTES;
  const uint8_t* x1 = x + (ell + 2 * BF128_NUM_BYTES * 8) / 8;

  const unsigned int length_lambda = (ell + 3 * BF128_NUM_BYTES * 8 - 1) / (BF128_NUM_BYTES * 8);

  uint8_t tmp[BF128_NUM_BYTES] = {0};
  memcpy(tmp, x + (length_lambda - 1) * BF128_NUM_BYTES,
         (ell + BF128_NUM_BYTES * 8) % (BF128_NUM_BYTES * 8) == 0
             ? BF128_NUM_BYTES
             : ((ell + BF128_NUM_BYTES * 8) % (BF128_NUM_BYTES * 8)) / 8);
  bf128_t h0 = bf128_load(tmp);

  const bf128_t b_s = bf128_load(s);
  bf128_t running_s = b_s;
  for (unsigned int i = 1; i != length_lambda; ++i, running_s = bf128_mul(running_s, b_s)) {
    h0 = bf128_add(h0,
                   bf128_mul(running_s, bf128_load(x + (length_lambda - 1 - i) * BF128_NUM_BYTES)));
  }

  bf64_t h1  = compute_h1(t, x, BF128_NUM_BYTES * 8, ell);
  bf128_t h2 = bf128_add(bf128_mul(bf128_load(r0), h0), bf128_mul_64(bf128_load(r1), h1));
  bf128_t h3 = bf128_add(bf128_mul(bf128_load(r2), h0), bf128_mul_64(bf128_load(r3), h1));

  bf128_store(h, h2);
  bf128_store(tmp, h3);
  memcpy(h + BF128_NUM_BYTES, tmp, UNIVERSAL_HASH_B);
  xor_u8_array(h, x1, h, BF128_NUM_BYTES + UNIVERSAL_HASH_B);
}

void vole_hash_192(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell) {
  const uint8_t* r0 = sd;
  const uint8_t* r1 = sd + 1 * BF192_NUM_BYTES;
  const uint8_t* r2 = sd + 2 * BF192_NUM_BYTES;
  const uint8_t* r3 = sd + 3 * BF192_NUM_BYTES;
  const uint8_t* s  = sd + 4 * BF192_NUM_BYTES;
  const uint8_t* t  = sd + 5 * BF192_NUM_BYTES;
  const uint8_t* x1 = x + (ell + 2 * BF192_NUM_BYTES * 8) / 8;

  const unsigned int length_lambda = (ell + 3 * BF192_NUM_BYTES * 8 - 1) / (BF192_NUM_BYTES * 8);

  uint8_t tmp[BF192_NUM_BYTES] = {0};
  memcpy(tmp, x + (length_lambda - 1) * BF192_NUM_BYTES,
         (ell + BF192_NUM_BYTES * 8) % (BF192_NUM_BYTES * 8) == 0
             ? BF192_NUM_BYTES
             : ((ell + BF192_NUM_BYTES * 8) % (BF192_NUM_BYTES * 8)) / 8);
  bf192_t h0 = bf192_load(tmp);

  const bf192_t b_s = bf192_load(s);
  bf192_t running_s = b_s;
  for (unsigned int i = 1; i != length_lambda; ++i, running_s = bf192_mul(running_s, b_s)) {
    h0 = bf192_add(h0,
                   bf192_mul(running_s, bf192_load(x + (length_lambda - 1 - i) * BF192_NUM_BYTES)));
  }

  bf64_t h1  = compute_h1(t, x, BF192_NUM_BYTES * 8, ell);
  bf192_t h2 = bf192_add(bf192_mul(bf192_load(r0), h0), bf192_mul_64(bf192_load(r1), h1));
  bf192_t h3 = bf192_add(bf192_mul(bf192_load(r2), h0), bf192_mul_64(bf192_load(r3), h1));

  bf192_store(h, h2);
  bf192_store(tmp, h3);
  memcpy(h + BF192_NUM_BYTES, tmp, UNIVERSAL_HASH_B);
  xor_u8_array(h, x1, h, BF192_NUM_BYTES + UNIVERSAL_HASH_B);
}

void vole_hash_256(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int ell) {
  const uint8_t* r0 = sd;
  const uint8_t* r1 = sd + 1 * BF256_NUM_BYTES;
  const uint8_t* r2 = sd + 2 * BF256_NUM_BYTES;
  const uint8_t* r3 = sd + 3 * BF256_NUM_BYTES;
  const uint8_t* s  = sd + 4 * BF256_NUM_BYTES;
  const uint8_t* t  = sd + 5 * BF256_NUM_BYTES;
  const uint8_t* x1 = x + (ell + 2 * BF256_NUM_BYTES * 8) / 8;

  const unsigned int length_lambda = (ell + 3 * BF256_NUM_BYTES * 8 - 1) / (BF256_NUM_BYTES * 8);

  uint8_t tmp[BF256_NUM_BYTES] = {0};
  memcpy(tmp, x + (length_lambda - 1) * BF256_NUM_BYTES,
         (ell + BF256_NUM_BYTES * 8) % (BF256_NUM_BYTES * 8) == 0
             ? BF256_NUM_BYTES
             : ((ell + BF256_NUM_BYTES * 8) % (BF256_NUM_BYTES * 8)) / 8);
  bf256_t h0 = bf256_load(tmp);

  const bf256_t b_s = bf256_load(s);
  bf256_t running_s = b_s;
  for (unsigned int i = 1; i != length_lambda; ++i, running_s = bf256_mul(running_s, b_s)) {
    h0 = bf256_add(h0,
                   bf256_mul(running_s, bf256_load(x + (length_lambda - 1 - i) * BF256_NUM_BYTES)));
  }

  bf64_t h1  = compute_h1(t, x, BF256_NUM_BYTES * 8, ell);
  bf256_t h2 = bf256_add(bf256_mul(bf256_load(r0), h0), bf256_mul_64(bf256_load(r1), h1));
  bf256_t h3 = bf256_add(bf256_mul(bf256_load(r2), h0), bf256_mul_64(bf256_load(r3), h1));

  bf256_store(h, h2);
  bf256_store(tmp, h3);
  memcpy(h + BF256_NUM_BYTES, tmp, UNIVERSAL_HASH_B);
  xor_u8_array(h, x1, h, BF256_NUM_BYTES + UNIVERSAL_HASH_B);
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

void zk_hash_128_init(zk_hash_128_ctx* ctx, const uint8_t* sd) {
  const uint8_t* s = sd + 2 * BF128_NUM_BYTES;
  const uint8_t* t = sd + 3 * BF128_NUM_BYTES;

  ctx->h0 = bf128_zero();
  ctx->h1 = bf128_zero();
  ctx->s  = bf128_load(s);
  ctx->t  = bf64_load(t);
  ctx->sd = sd;
}

void zk_hash_128_update(zk_hash_128_ctx* ctx, bf128_t v) {
  ctx->h0 = bf128_add(bf128_mul(ctx->h0, ctx->s), v);
  ctx->h1 = bf128_add(bf128_mul_64(ctx->h1, ctx->t), v);
}

void zk_hash_128_finalize(uint8_t* h, zk_hash_128_ctx* ctx, bf128_t x1) {
  const uint8_t* r0 = ctx->sd;
  const uint8_t* r1 = ctx->sd + BF128_NUM_BYTES;

  bf128_store(h, bf128_add(bf128_add(bf128_mul(bf128_load(r0), ctx->h0),
                                     bf128_mul(bf128_load(r1), ctx->h1)),
                           x1));
}

void zk_hash_128_3_init(zk_hash_128_3_ctx* ctx, const uint8_t* sd) {
  const uint8_t* s = sd + 2 * BF128_NUM_BYTES;
  const uint8_t* t = sd + 3 * BF128_NUM_BYTES;

  ctx->h0[0] = bf128_zero();
  ctx->h0[1] = bf128_zero();
  ctx->h0[2] = bf128_zero();
  ctx->h1[0] = bf128_zero();
  ctx->h1[1] = bf128_zero();
  ctx->h1[2] = bf128_zero();
  ctx->s     = bf128_load(s);
  ctx->t     = bf64_load(t);
  ctx->sd    = sd;
}

void zk_hash_128_3_update(zk_hash_128_3_ctx* ctx, bf128_t v_0, bf128_t v_1, bf128_t v_2) {
  ctx->h0[0] = bf128_add(bf128_mul(ctx->h0[0], ctx->s), v_0);
  ctx->h1[0] = bf128_add(bf128_mul_64(ctx->h1[0], ctx->t), v_0);
  ctx->h0[1] = bf128_add(bf128_mul(ctx->h0[1], ctx->s), v_1);
  ctx->h1[1] = bf128_add(bf128_mul_64(ctx->h1[1], ctx->t), v_1);
  ctx->h0[2] = bf128_add(bf128_mul(ctx->h0[2], ctx->s), v_2);
  ctx->h1[2] = bf128_add(bf128_mul_64(ctx->h1[2], ctx->t), v_2);
}

void zk_hash_128_3_raise_and_update(zk_hash_128_3_ctx* ctx, bf128_t v_1, bf128_t v_2) {
  zk_hash_128_3_update(ctx, bf128_zero(), v_1, v_2);
}

void zk_hash_128_3_finalize(uint8_t* h_0, uint8_t* h_1, uint8_t* h_2, zk_hash_128_3_ctx* ctx,
                            bf128_t x1_0, bf128_t x1_1, bf128_t x1_2) {
  const bf128_t r0 = bf128_load(ctx->sd);
  const bf128_t r1 = bf128_load(ctx->sd + BF128_NUM_BYTES);

  bf128_store(h_0,
              bf128_add(bf128_add(bf128_mul(r0, ctx->h0[0]), bf128_mul(r1, ctx->h1[0])), x1_0));
  bf128_store(h_1,
              bf128_add(bf128_add(bf128_mul(r0, ctx->h0[1]), bf128_mul(r1, ctx->h1[1])), x1_1));
  bf128_store(h_2,
              bf128_add(bf128_add(bf128_mul(r0, ctx->h0[2]), bf128_mul(r1, ctx->h1[2])), x1_2));
}

void zk_hash_192_init(zk_hash_192_ctx* ctx, const uint8_t* sd) {
  const uint8_t* s = sd + 2 * BF192_NUM_BYTES;
  const uint8_t* t = sd + 3 * BF192_NUM_BYTES;

  ctx->h0 = bf192_zero();
  ctx->h1 = bf192_zero();
  ctx->s  = bf192_load(s);
  ctx->t  = bf64_load(t);
  ctx->sd = sd;
}

void zk_hash_192_update(zk_hash_192_ctx* ctx, bf192_t v) {
  ctx->h0 = bf192_add(bf192_mul(ctx->h0, ctx->s), v);
  ctx->h1 = bf192_add(bf192_mul_64(ctx->h1, ctx->t), v);
}

void zk_hash_192_finalize(uint8_t* h, zk_hash_192_ctx* ctx, bf192_t x1) {
  const uint8_t* r0 = ctx->sd;
  const uint8_t* r1 = ctx->sd + BF192_NUM_BYTES;

  bf192_store(h, bf192_add(bf192_add(bf192_mul(bf192_load(r0), ctx->h0),
                                     bf192_mul(bf192_load(r1), ctx->h1)),
                           x1));
}

void zk_hash_192_3_init(zk_hash_192_3_ctx* ctx, const uint8_t* sd) {
  const uint8_t* s = sd + 2 * BF192_NUM_BYTES;
  const uint8_t* t = sd + 3 * BF192_NUM_BYTES;

  ctx->h0[0] = bf192_zero();
  ctx->h0[1] = bf192_zero();
  ctx->h0[2] = bf192_zero();
  ctx->h1[0] = bf192_zero();
  ctx->h1[1] = bf192_zero();
  ctx->h1[2] = bf192_zero();
  ctx->s     = bf192_load(s);
  ctx->t     = bf64_load(t);
  ctx->sd    = sd;
}

void zk_hash_192_3_update(zk_hash_192_3_ctx* ctx, bf192_t v_0, bf192_t v_1, bf192_t v_2) {
  ctx->h0[0] = bf192_add(bf192_mul(ctx->h0[0], ctx->s), v_0);
  ctx->h1[0] = bf192_add(bf192_mul_64(ctx->h1[0], ctx->t), v_0);
  ctx->h0[1] = bf192_add(bf192_mul(ctx->h0[1], ctx->s), v_1);
  ctx->h1[1] = bf192_add(bf192_mul_64(ctx->h1[1], ctx->t), v_1);
  ctx->h0[2] = bf192_add(bf192_mul(ctx->h0[2], ctx->s), v_2);
  ctx->h1[2] = bf192_add(bf192_mul_64(ctx->h1[2], ctx->t), v_2);
}

void zk_hash_192_3_raise_and_update(zk_hash_192_3_ctx* ctx, bf192_t v_1, bf192_t v_2) {
  zk_hash_192_3_update(ctx, bf192_zero(), v_1, v_2);
}

void zk_hash_192_3_finalize(uint8_t* h_0, uint8_t* h_1, uint8_t* h_2, zk_hash_192_3_ctx* ctx,
                            bf192_t x1_0, bf192_t x1_1, bf192_t x1_2) {
  const bf192_t r0 = bf192_load(ctx->sd);
  const bf192_t r1 = bf192_load(ctx->sd + BF192_NUM_BYTES);

  bf192_store(h_0,
              bf192_add(bf192_add(bf192_mul(r0, ctx->h0[0]), bf192_mul(r1, ctx->h1[0])), x1_0));
  bf192_store(h_1,
              bf192_add(bf192_add(bf192_mul(r0, ctx->h0[1]), bf192_mul(r1, ctx->h1[1])), x1_1));
  bf192_store(h_2,
              bf192_add(bf192_add(bf192_mul(r0, ctx->h0[2]), bf192_mul(r1, ctx->h1[2])), x1_2));
}

void zk_hash_256_init(zk_hash_256_ctx* ctx, const uint8_t* sd) {
  const uint8_t* s = sd + 2 * BF256_NUM_BYTES;
  const uint8_t* t = sd + 3 * BF256_NUM_BYTES;

  ctx->h0 = bf256_zero();
  ctx->h1 = bf256_zero();
  ctx->s  = bf256_load(s);
  ctx->t  = bf64_load(t);
  ctx->sd = sd;
}

void zk_hash_256_update(zk_hash_256_ctx* ctx, bf256_t v) {
  ctx->h0 = bf256_add(bf256_mul(ctx->h0, ctx->s), v);
  ctx->h1 = bf256_add(bf256_mul_64(ctx->h1, ctx->t), v);
}

void zk_hash_256_finalize(uint8_t* h, zk_hash_256_ctx* ctx, bf256_t x1) {
  const uint8_t* r0 = ctx->sd;
  const uint8_t* r1 = ctx->sd + BF256_NUM_BYTES;

  bf256_store(h, bf256_add(bf256_add(bf256_mul(bf256_load(r0), ctx->h0),
                                     bf256_mul(bf256_load(r1), ctx->h1)),
                           x1));
}

void zk_hash_256_3_init(zk_hash_256_3_ctx* ctx, const uint8_t* sd) {
  const uint8_t* s = sd + 2 * BF256_NUM_BYTES;
  const uint8_t* t = sd + 3 * BF256_NUM_BYTES;

  ctx->h0[0] = bf256_zero();
  ctx->h0[1] = bf256_zero();
  ctx->h0[2] = bf256_zero();
  ctx->h1[0] = bf256_zero();
  ctx->h1[1] = bf256_zero();
  ctx->h1[2] = bf256_zero();
  ctx->s     = bf256_load(s);
  ctx->t     = bf64_load(t);
  ctx->sd    = sd;
}

void zk_hash_256_3_update(zk_hash_256_3_ctx* ctx, bf256_t v_0, bf256_t v_1, bf256_t v_2) {
  ctx->h0[0] = bf256_add(bf256_mul(ctx->h0[0], ctx->s), v_0);
  ctx->h1[0] = bf256_add(bf256_mul_64(ctx->h1[0], ctx->t), v_0);
  ctx->h0[1] = bf256_add(bf256_mul(ctx->h0[1], ctx->s), v_1);
  ctx->h1[1] = bf256_add(bf256_mul_64(ctx->h1[1], ctx->t), v_1);
  ctx->h0[2] = bf256_add(bf256_mul(ctx->h0[2], ctx->s), v_2);
  ctx->h1[2] = bf256_add(bf256_mul_64(ctx->h1[2], ctx->t), v_2);
}

void zk_hash_256_3_raise_and_update(zk_hash_256_3_ctx* ctx, bf256_t v_1, bf256_t v_2) {
  zk_hash_256_3_update(ctx, bf256_zero(), v_1, v_2);
}

void zk_hash_256_3_finalize(uint8_t* h_0, uint8_t* h_1, uint8_t* h_2, zk_hash_256_3_ctx* ctx,
                            bf256_t x1_0, bf256_t x1_1, bf256_t x1_2) {
  const bf256_t r0 = bf256_load(ctx->sd);
  const bf256_t r1 = bf256_load(ctx->sd + BF256_NUM_BYTES);

  bf256_store(h_0,
              bf256_add(bf256_add(bf256_mul(r0, ctx->h0[0]), bf256_mul(r1, ctx->h1[0])), x1_0));
  bf256_store(h_1,
              bf256_add(bf256_add(bf256_mul(r0, ctx->h0[1]), bf256_mul(r1, ctx->h1[1])), x1_1));
  bf256_store(h_2,
              bf256_add(bf256_add(bf256_mul(r0, ctx->h0[2]), bf256_mul(r1, ctx->h1[2])), x1_2));
}

#if defined(FAEST_TESTS)
void zk_hash_128(uint8_t* h, const uint8_t* sd, const bf128_t* x, unsigned int ell) {
  zk_hash_128_ctx ctx;
  zk_hash_128_init(&ctx, sd);
  for (unsigned int i = 0; i != ell; ++i) {
    zk_hash_128_update(&ctx, x[i]);
  }
  zk_hash_128_finalize(h, &ctx, x[ell]);
}

void zk_hash_192(uint8_t* h, const uint8_t* sd, const bf192_t* x, unsigned int ell) {
  zk_hash_192_ctx ctx;
  zk_hash_192_init(&ctx, sd);
  for (unsigned int i = 0; i != ell; ++i) {
    zk_hash_192_update(&ctx, x[i]);
  }
  zk_hash_192_finalize(h, &ctx, x[ell]);
}

void zk_hash_256(uint8_t* h, const uint8_t* sd, const bf256_t* x, unsigned int ell) {
  zk_hash_256_ctx ctx;
  zk_hash_256_init(&ctx, sd);
  for (unsigned int i = 0; i != ell; ++i) {
    zk_hash_256_update(&ctx, x[i]);
  }
  zk_hash_256_finalize(h, &ctx, x[ell]);
}
#endif

void leaf_hash_128(uint8_t* h, const uint8_t* uhash, const uint8_t* x) {
  const uint8_t* x0 = x;
  const uint8_t* x1 = x + BF128_NUM_BYTES;

  bf384_store(h, bf384_add(bf384_mul_128(bf384_load(uhash), bf128_load(x0)), bf384_load(x1)));
}

void leaf_hash_192(uint8_t* h, const uint8_t* uhash, const uint8_t* x) {
  const uint8_t* x0 = x;
  const uint8_t* x1 = x + BF192_NUM_BYTES;

  bf576_store(h, bf576_add(bf576_mul_192(bf576_load(uhash), bf192_load(x0)), bf576_load(x1)));
}

void leaf_hash_256(uint8_t* h, const uint8_t* uhash, const uint8_t* x) {
  const uint8_t* x0 = x;
  const uint8_t* x1 = x + BF256_NUM_BYTES;

  bf768_store(h, bf768_add(bf768_mul_256(bf768_load(uhash), bf256_load(x0)), bf768_load(x1)));
}

void leaf_hash(uint8_t* h, const uint8_t* sd, const uint8_t* x, unsigned int lambda) {
  switch (lambda) {
  case 256:
    leaf_hash_256(h, sd, x);
    break;
  case 192:
    leaf_hash_192(h, sd, x);
    break;
  default:
    leaf_hash_128(h, sd, x);
    break;
  }
}
