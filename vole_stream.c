#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "vole_stream.h"
#include "aes.h"
#include "utils.h"
#include "random_oracle.h"
#include "vole.h"
#include "vc.h"

#include <stdbool.h>
#include <string.h>
#include <stdio.h>

static void ConstructVoleCMO(const uint8_t* iv, stream_vec_com_t* sVecCom, unsigned int lambda,
                             unsigned int outLenBytes, uint8_t* u, uint8_t* v, uint8_t* h,
                             unsigned int begin, unsigned int end) {
  unsigned int depth               = sVecCom->depth;
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;
  unsigned int len                 = end - begin;

#define V_CMO(idx) (v + ((idx)-begin) * outLenBytes)

  uint8_t* sd  = malloc(lambda_bytes);
  uint8_t* com = malloc(lambda_bytes * 2);
  H1_context_t* h1_ctx;
  if (h != NULL) {
    h1_ctx = malloc(sizeof(H1_context_t));
    H1_init(h1_ctx, lambda);
  }

  uint8_t* r = malloc(outLenBytes);

  // Clear initial memory
  if (u != NULL) {
    memset(u, 0, outLenBytes);
  }
  if (v != NULL) {
    memset(v, 0, len * outLenBytes);
  }

  for (unsigned int i = 0; i < num_instances; i++) {
    get_sd_com(sVecCom, iv, lambda, i, sd, com);
    if (h != NULL) {
      H1_update(h1_ctx, com, lambda_bytes * 2);
    }

    // Seed expansion
    prg(sd, iv, r, lambda, outLenBytes);
    if (u != NULL) {
      xor_u8_array(u, r, u, outLenBytes);
    }
    if (v != NULL) {
      for (unsigned int j = begin; j < end; j++) {
        // Apply r if the j'th bit is set
        if ((i >> j) & 1) {
          xor_u8_array(V_CMO(j), r, V_CMO(j), outLenBytes);
        }
      }
    }
  }

  if (h != NULL) {
    H1_final(h1_ctx, h, lambda_bytes * 2);
  }

  if (h != NULL) {
    free(h1_ctx);
  }
  free(sd);
  free(com);
  free(r);
}

// NOTE - Assumes v is cleared (initially)!!
static void ConstructVoleRMO(const uint8_t* iv, unsigned int start, unsigned int len,
                             stream_vec_com_t* sVecCom, unsigned int lambda,
                             unsigned int outLenBytes, uint8_t* v, unsigned int offset) {
  unsigned int depth               = sVecCom->depth;
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;

  uint8_t* sd              = malloc(lambda_bytes);
  uint8_t* com             = malloc(lambda_bytes * 2);
  uint8_t* r               = malloc(outLenBytes);
  unsigned int bit_offset  = (offset % 8);
  unsigned int byte_offset = (offset / 8);

  for (unsigned int i = 0; i < num_instances; i++) {
    get_sd_com(sVecCom, iv, lambda, i, sd, com);
    prg(sd, iv, r, lambda, outLenBytes);

    for (unsigned int row_idx = 0; row_idx < len; row_idx++) {
      unsigned int byte_idx = (row_idx + start) / 8;
      unsigned int bit_idx  = (row_idx + start) % 8;
      uint8_t bit           = (r[byte_idx] >> (bit_idx)) & 1;
      if (bit == 0) {
        continue;
      }
      unsigned int base_idx = row_idx * lambda_bytes + byte_offset;
      unsigned int amount   = (bit_offset + depth + 7) / 8;
      // Avoid carry by breaking into two steps
      v[base_idx] ^= i << bit_offset;
      for (unsigned int j = 1; j < amount; j++) {
        v[base_idx + j] ^= i >> (j * 8 - bit_offset);
      }
    }
  }
  free(sd);
  free(com);
  free(r);
}

void partial_vole_commit_cmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                             const faest_paramset_t* params, stream_vec_com_t* sVecCom, uint8_t* v,
                             unsigned int start, unsigned int len, uint8_t* u, uint8_t* hcom,
                             uint8_t* c) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;
  unsigned int max_depth    = MAX(k0, k1);

  uint8_t* expanded_keys = malloc(tau * lambda_bytes);
  prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);
  uint8_t* path = malloc(lambda_bytes * max_depth);

  H1_context_t h1_ctx;
  uint8_t* h = NULL;
  if (hcom != NULL) {
    H1_init(&h1_ctx, lambda);
    h = malloc(lambda_bytes * 2);
  }

  unsigned int end        = start + len;
  unsigned int tree_start = 0;
  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;

    if (tree_start + depth > start) {
      stream_vector_commitment(expanded_keys + i * lambda_bytes, lambda, &sVecCom[i], depth);

      unsigned int lbegin = 0;
      if (start > tree_start) {
        lbegin = start - tree_start;
      }

      unsigned int lend  = MIN(depth, end - tree_start);
      unsigned int v_idx = 0;
      if (tree_start > start) {
        v_idx = tree_start - start;
      }

      // This is for optimization when len = ell_hat
      uint8_t* u_ptr = NULL;
      if (u != NULL && i == 0) {
        u_ptr = u;
      } else if (c != NULL && i != 0) {
        u_ptr = c + (i - 1) * ellhat_bytes;
      }

      sVecCom[i].path = path;
      ConstructVoleCMO(iv, &sVecCom[i], lambda, ellhat_bytes, u_ptr, v + ellhat_bytes * v_idx, h,
                       lbegin, lend);
      sVecCom[i].path = NULL;

      if (c != NULL && u != NULL && i >= 1) {
        xor_u8_array(u, u_ptr, u_ptr, ellhat_bytes);
      }

      if (hcom != NULL) {
        H1_update(&h1_ctx, h, lambda_bytes * 2);
      }
    }

    tree_start += depth;
    if (tree_start >= end) {
      break;
    }
  }
  if (hcom != NULL) {
    H1_final(&h1_ctx, hcom, lambda_bytes * 2);
    free(h);
  }

  free(expanded_keys);
  free(path);
}

void vole_commit_u_hcom_c(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                          const faest_paramset_t* params, uint8_t* hcom, stream_vec_com_t* sVecCom,
                          uint8_t* c, uint8_t* u) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;
  unsigned int max_depth    = MAX(k0, k1);

  uint8_t* expanded_keys = malloc(tau * lambda_bytes);
  prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);

  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  uint8_t* h    = malloc(lambda_bytes * 2);
  uint8_t* path = malloc(lambda_bytes * max_depth);

  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;
    stream_vector_commitment(expanded_keys + i * lambda_bytes, lambda, &sVecCom[i], depth);

    uint8_t* u_ptr = NULL;
    if (i == 0) {
      u_ptr = u;
    } else {
      u_ptr = c + (i - 1) * ellhat_bytes;
    }

    sVecCom[i].path = path;
    ConstructVoleCMO(iv, &sVecCom[i], lambda, ellhat_bytes, u_ptr, NULL, h, 0, depth);
    sVecCom[i].path = NULL;

    if (i != 0) {
      xor_u8_array(u, u_ptr, u_ptr, ellhat_bytes);
    }
    H1_update(&h1_ctx, h, lambda_bytes * 2);
  }

  H1_final(&h1_ctx, hcom, lambda_bytes * 2);
  free(expanded_keys);
  free(h);
  free(path);
}

void partial_vole_commit_rmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int start,
                             unsigned int len, const faest_paramset_t* params,
                             stream_vec_com_t* sVecCom, uint8_t* v) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ell_hat =
      params->faest_param.l + params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ell_hat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;
  unsigned int max_depth    = MAX(k0, k1);

  uint8_t* expanded_keys = malloc(tau * lambda_bytes);
  prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);

  uint8_t* path = malloc(lambda_bytes * max_depth);

  memset(v, 0, ((size_t)len) * (size_t)lambda_bytes);

  unsigned int col_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;

    stream_vector_commitment(expanded_keys + i * lambda_bytes, lambda, &sVecCom[i], depth);
    sVecCom[i].path = path;
    ConstructVoleRMO(iv, start, len, &sVecCom[i], lambda, ellhat_bytes, v, col_idx);
    sVecCom[i].path = NULL;

    col_idx += depth;
  }

  free(expanded_keys);
  free(path);
}

// reconstruction
static void ReconstructVoleCMO(const uint8_t* iv, stream_vec_com_rec_t* sVecComRec,
                               unsigned int lambda, unsigned int outLenBytes, uint8_t* q,
                               uint8_t* h, unsigned int begin, unsigned int end) {
  unsigned int depth               = sVecComRec->depth;
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;
  unsigned int len                 = end - begin;

#define Q_CMO(idx) (q + ((idx)-begin) * outLenBytes)

  H1_context_t h1_ctx;
  if (h != NULL) {
    H1_init(&h1_ctx, lambda);
  }
  uint8_t* sd  = malloc(lambda_bytes);
  uint8_t* com = malloc(lambda_bytes * 2);
  uint8_t* r;

  unsigned int offset = NumRec(depth, sVecComRec->b);

  if (q != NULL) {
    r = malloc(outLenBytes);
    memset(q, 0, len * outLenBytes);
  }

  for (unsigned int i = 0; i < num_instances; i++) {
    unsigned int offset_index = i ^ offset;
    if (offset_index == 0) {
      if (h != NULL) {
        H1_update(&h1_ctx, sVecComRec->com_j, lambda_bytes * 2);
      }
      continue;
    }

    get_sd_com_rec(sVecComRec, iv, lambda, i, sd, com);
    if (h != NULL) {
      H1_update(&h1_ctx, com, lambda_bytes * 2);
    }

    if (q == NULL) {
      continue;
    }

    prg(sd, iv, r, lambda, outLenBytes);
    for (unsigned int j = begin; j < end; j++) {
      // Apply r if the j'th bit is set
      if ((offset_index >> j) & 1) {
        xor_u8_array(Q_CMO(j), r, Q_CMO(j), outLenBytes);
      }
    }
  }
  if (q != NULL) {
    free(r);
  }
  free(sd);
  free(com);
  if (h != NULL) {
    H1_final(&h1_ctx, h, lambda_bytes * 2);
  }
}

void partial_vole_reconstruct_cmo(const uint8_t* iv, const uint8_t* chall,
                                  const uint8_t* const* pdec, const uint8_t* const* com_j,
                                  uint8_t* hcom, uint8_t* q, unsigned int ellhat,
                                  const faest_paramset_t* params, unsigned int start,
                                  unsigned int len) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int tau1         = params->faest_param.t1;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;

  H1_context_t h1_ctx;
  if (hcom != NULL) {
    H1_init(&h1_ctx, lambda);
  }

  stream_vec_com_rec_t sVecComRec;
  unsigned int max_depth = MAX(k0, k1);
  sVecComRec.b           = alloca(max_depth * sizeof(uint8_t));
  sVecComRec.nodes       = calloc(max_depth, lambda_bytes);
  sVecComRec.com_j       = alloca(lambda_bytes * 2);
  sVecComRec.path        = alloca(lambda_bytes * (max_depth - 1));

  uint8_t* h = NULL;
  if (hcom != NULL) {
    h = alloca(lambda_bytes * 2);
  }

  unsigned int end        = start + len;
  unsigned int tree_start = 0;

  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;
    if (tree_start + depth > start) {
      uint8_t chalout[MAX_DEPTH];
      ChalDec(chall, i, k0, tau0, k1, tau1, chalout);
      stream_vector_reconstruction(pdec[i], com_j[i], chalout, lambda, depth, &sVecComRec);

      unsigned int lbegin = 0;
      if (start > tree_start) {
        lbegin = start - tree_start;
      }

      unsigned int lend  = MIN(depth, end - tree_start);
      unsigned int q_idx = 0;
      if (tree_start > start) {
        q_idx = tree_start - start;
      }

      ReconstructVoleCMO(iv, &sVecComRec, lambda, ellhat_bytes, q + q_idx * ellhat_bytes, h, lbegin,
                         lend);
      if (hcom != NULL) {
        H1_update(&h1_ctx, h, lambda_bytes * 2);
      }
    }

    tree_start += depth;
    if (tree_start >= end) {
      break;
    }
  }

  free(sVecComRec.nodes);
  if (hcom != NULL) {
    H1_final(&h1_ctx, hcom, lambda_bytes * 2);
  }
}

void vole_reconstruct_hcom(const uint8_t* iv, const uint8_t* chall, const uint8_t* const* pdec,
                           const uint8_t* const* com_j, uint8_t* hcom, unsigned int ellhat,
                           const faest_paramset_t* params) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int tau1         = params->faest_param.t1;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;

  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);

  stream_vec_com_rec_t sVecComRec;
  unsigned int max_depth = MAX(k0, k1);
  sVecComRec.b           = alloca(max_depth * sizeof(uint8_t));
  sVecComRec.nodes       = calloc(max_depth, lambda_bytes);
  sVecComRec.com_j       = alloca(lambda_bytes * 2);
  sVecComRec.path        = alloca(lambda_bytes * (max_depth - 1));

  uint8_t* h              = alloca(lambda_bytes * 2);
  unsigned int tree_start = 0;

  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;

    uint8_t chalout[MAX_DEPTH];
    ChalDec(chall, i, k0, tau0, k1, tau1, chalout);
    stream_vector_reconstruction(pdec[i], com_j[i], chalout, lambda, depth, &sVecComRec);

    ReconstructVoleCMO(iv, &sVecComRec, lambda, ellhat_bytes, NULL, h, 0, depth);

    H1_update(&h1_ctx, h, lambda_bytes * 2);

    tree_start += depth;
  }

  free(sVecComRec.nodes);

  H1_final(&h1_ctx, hcom, lambda_bytes * 2);
}
