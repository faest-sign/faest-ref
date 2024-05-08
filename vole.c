/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "vole.h"
#include "aes.h"
#include "utils.h"
#include "random_oracle.h"

#include <stdbool.h>
#include <string.h>

int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout) {
  if (i >= t0 + t1) {
    return 0;
  }

  unsigned int lo;
  unsigned int hi;
  if (i < t0) {
    lo = i * k0;
    hi = ((i + 1) * k0);
  } else {
    unsigned int t = i - t0;
    lo             = (t0 * k0) + (t * k1);
    hi             = (t0 * k0) + ((t + 1) * k1);
  }

  assert(hi - lo == k0 || hi - lo == k1);
  for (unsigned int j = lo; j < hi; ++j) {
    // set_bit(chalout, i - lo, get_bit(chal, i));
    chalout[j - lo] = ptr_get_bit(chal, j);
  }
  return 1;
}

// NOTE - Assumes v is cleared (initially)!!
static void ConstructVoleRMO(const uint8_t* iv, unsigned int start, unsigned int len,
                             vec_com_t* vec_com, unsigned int lambda, unsigned int out_len_bytes,
                             uint8_t* v, unsigned int offset) {
  unsigned int depth               = vec_com->depth;
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;

  uint8_t* sd              = malloc(lambda_bytes);
  uint8_t* com             = malloc(lambda_bytes * 2);
  uint8_t* r               = malloc(out_len_bytes);
  unsigned int bit_offset  = (offset % 8);
  unsigned int byte_offset = (offset / 8);

  for (unsigned int i = 0; i < num_instances; i++) {
    extract_sd_com(vec_com, iv, lambda, i, sd, com);
    prg(sd, iv, r, lambda, out_len_bytes);

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
                             unsigned int start, unsigned int end,
                             sign_vole_mode_ctx_t vole_mode, const faest_paramset_t* params) {
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

  H1_context_t hcom_ctx;
  H1_context_t com_ctx;
  uint8_t* h = NULL;
  if (vole_mode.mode != EXCLUDE_U_HCOM_C) {
    h = malloc(lambda_bytes * 2);
    H1_init(&hcom_ctx, lambda);
  }

  // STEP 1: To commit to [start,end] we first compute which trees we need to consider
  unsigned int depth_tau_0     = tau0 * k0;

  unsigned int k0_trees_begin  = (start < depth_tau_0) ? start / k0 : tau0;
  unsigned int k1_trees_begin  = (start < depth_tau_0) ? 0 : (start - depth_tau_0) / k1;
  unsigned int k0_trees_end    = (end < depth_tau_0) ? (end + (k0-1)) / k0 : tau0; // ceiled
  unsigned int k1_trees_end    = (end < depth_tau_0) ? 0 : (end - depth_tau_0 + (k1-1)) / k1;  // ceiled

  unsigned int tree_start = k0_trees_begin+k1_trees_begin;
  unsigned int tree_end   = k0_trees_end+k1_trees_end;

  // Compute the cummulative sum of the tree depths until the requested start
  unsigned int v_progress = k0 * k0_trees_begin + k1 * k1_trees_begin;
  
  for (unsigned int t = tree_start; t < tree_end; t++) {
    bool is_first_tree      = (t == 0);
    unsigned int tree_depth = t < tau0 ? k0 : k1;

    // v_cache_offset is used to compute the index we should write v to relative to our cache
    unsigned int v_cache_offset  = (v_progress > start) ? v_progress - start : 0; // (i.e. MAX(v_progress-start, 0))
    // [v_start, v_end] is the v's that t provides (capped by requested start/end)
    unsigned int v_start         = MAX(v_progress, start); 
    unsigned int v_end           = MIN(end, v_progress+tree_depth);

    // (Setup for STEP 2)
    const unsigned int num_seeds = 1 << tree_depth;
    uint8_t sd[MAX_LAMBDA_BYTES];
    uint8_t com[2*MAX_LAMBDA_BYTES];
    uint8_t* r = malloc(ellhat_bytes);

    uint8_t* u_ptr = NULL;
    if (vole_mode.mode != EXCLUDE_U_HCOM_C) {
      H1_init(&com_ctx, lambda);
      u_ptr = is_first_tree ? vole_mode.u : vole_mode.c + (t - 1) * ellhat_bytes;
      memset(u_ptr, 0, ellhat_bytes);
    }
    if (vole_mode.mode != EXCLUDE_V) {
      unsigned int v_count = v_end-v_start;
      memset(vole_mode.v+(v_cache_offset*ellhat_bytes), 0, v_count*ellhat_bytes);
    }

    vec_com_t vec_com;
    vector_commitment(expanded_keys + t * lambda_bytes, lambda, tree_depth, path, &vec_com);

    // STEP 2: For this tree, extract all seeds and commitments and compute according to the VOLE-mode
    for (unsigned int i = 0; i < num_seeds; i++) {
      extract_sd_com(&vec_com, iv, lambda, i, sd, com);
      prg(sd, iv, r, lambda, ellhat_bytes); // Seed expansion

      if (vole_mode.mode != EXCLUDE_U_HCOM_C) {
        H1_update(&com_ctx, com, lambda_bytes * 2);
        xor_u8_array(u_ptr, r, u_ptr, ellhat_bytes);
      }
      if (vole_mode.mode != EXCLUDE_V) {
        for (unsigned int j = v_start; j < v_end; j++) {
          // Instead of writing v_j at V[j], use the v_cache_offset
          uint8_t* write_idx = (vole_mode.v+(j-v_start+v_cache_offset) * ellhat_bytes);
          unsigned int t_v = j-v_progress; // That is; t provides depth num of v's where t_v reflects the current v \in [0, depth]
          // Apply r if the i/2^t_v is odd
          if ((i >> t_v) & 1) {
            xor_u8_array(write_idx, r, write_idx, ellhat_bytes);
          }
        }
      }
    }

    free(r);

    if (vole_mode.mode != EXCLUDE_U_HCOM_C) {
      if (!is_first_tree) {
        xor_u8_array(vole_mode.u, u_ptr, u_ptr, ellhat_bytes); // Correction values
      }
      H1_final(&com_ctx, h, lambda_bytes * 2);
      H1_update(&hcom_ctx, h, lambda_bytes * 2);
    }

    v_progress += tree_depth;
  }

  if (vole_mode.mode != EXCLUDE_U_HCOM_C) {
    H1_final(&hcom_ctx, vole_mode.hcom, lambda_bytes * 2);
    free(h);
  }

  free(expanded_keys);
  free(path);
}

void partial_vole_commit_rmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int start,
                             unsigned int len, const faest_paramset_t* params, uint8_t* v) {
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
  vec_com_t vec_com;
  memset(v, 0, ((size_t)len) * (size_t)lambda_bytes);

  unsigned int col_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;

    vector_commitment(expanded_keys + i * lambda_bytes, lambda, depth, path, &vec_com);
    ConstructVoleRMO(iv, start, len, &vec_com, lambda, ellhat_bytes, v, col_idx);

    col_idx += depth;
  }
  free(expanded_keys);
  free(path);
}

// reconstruction
static void ReconstructVoleCMO(const uint8_t* iv, vec_com_rec_t* vec_com_rec, unsigned int lambda,
                               unsigned int out_len_bytes, uint8_t* q, uint8_t* h, unsigned int begin,
                               unsigned int end) {
  unsigned int depth               = vec_com_rec->depth;
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;
  unsigned int len                 = end - begin;

#define Q_CMO(idx) (q + ((idx)-begin) * out_len_bytes)

  H1_context_t h1_ctx;
  if (h != NULL) {
    H1_init(&h1_ctx, lambda);
  }
  uint8_t* sd  = malloc(lambda_bytes);
  uint8_t* com = malloc(lambda_bytes * 2);
  uint8_t* r;

  unsigned int offset = NumRec(depth, vec_com_rec->b);

  if (q != NULL) {
    r = malloc(out_len_bytes);
    memset(q, 0, len * out_len_bytes);
  }

  for (unsigned int i = 0; i < num_instances; i++) {
    unsigned int offset_index = i ^ offset;
    if (offset_index == 0) {
      if (h != NULL) {
        H1_update(&h1_ctx, vec_com_rec->com_j, lambda_bytes * 2);
      }
      continue;
    }

    extract_sd_com_rec(vec_com_rec, iv, lambda, i, sd, com);
    if (h != NULL) {
      H1_update(&h1_ctx, com, lambda_bytes * 2);
    }

    if (q == NULL) {
      continue;
    }

    prg(sd, iv, r, lambda, out_len_bytes);
    for (unsigned int j = begin; j < end; j++) {
      // Apply r if the j'th bit is set
      if ((offset_index >> j) & 1) {
        xor_u8_array(Q_CMO(j), r, Q_CMO(j), out_len_bytes);
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

void partial_vole_reconstruct_cmo(const uint8_t* iv, const uint8_t* chall, const uint8_t* const* pdec, 
                                  const uint8_t* const* com_j, unsigned int ellhat,
                                  unsigned int start, unsigned int len,
                                  verify_vole_mode_ctx_t vole_mode, 
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
  if (vole_mode.mode != EXCLUDE_HCOM) {
    H1_init(&h1_ctx, lambda);
  }

  unsigned int max_depth = MAX(k0, k1);
  vec_com_rec_t vec_com_rec;
  vec_com_rec.b       = malloc(max_depth * sizeof(uint8_t));
  vec_com_rec.nodes   = calloc(max_depth, lambda_bytes);
  vec_com_rec.com_j   = malloc(lambda_bytes * 2);
  uint8_t* tree_nodes = malloc(lambda_bytes * (max_depth - 1));

  uint8_t* h = NULL;
  if (vole_mode.mode != EXCLUDE_HCOM) {
    h = malloc(lambda_bytes * 2);
  }

  unsigned int end        = start + len;
  unsigned int tree_start = 0;

  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;
    if (tree_start + depth > start) {
      unsigned int lbegin = 0;
      if (start > tree_start) {
        lbegin = start - tree_start;
      }

      unsigned int lend  = MIN(depth, end - tree_start);
      unsigned int q_idx = 0;
      if (tree_start > start) {
        q_idx = tree_start - start;
      }

      uint8_t* q_ptr = NULL;
      if (vole_mode.mode != EXCLUDE_Q) {
        q_ptr = vole_mode.q + q_idx * ellhat_bytes;
      }

      uint8_t chalout[MAX_DEPTH];
      ChalDec(chall, i, k0, tau0, k1, tau1, chalout);
      vector_reconstruction(pdec[i], com_j[i], chalout, lambda, depth, tree_nodes, &vec_com_rec);
      ReconstructVoleCMO(iv, &vec_com_rec, lambda, ellhat_bytes, q_ptr, h,
                         lbegin, lend);
      if (vole_mode.mode != EXCLUDE_HCOM) {
        H1_update(&h1_ctx, h, lambda_bytes * 2);
      }
    }

    tree_start += depth;
    if (tree_start >= end) {
      break;
    }
  }

  free(vec_com_rec.b);
  free(vec_com_rec.nodes);
  free(vec_com_rec.com_j);
  free(tree_nodes);
  if (vole_mode.mode != EXCLUDE_HCOM) {
    free(h);
    H1_final(&h1_ctx, vole_mode.hcom, lambda_bytes * 2);
  }
}

static void ReconstructVoleRMO(const uint8_t* iv, vec_com_rec_t* vec_com_rec, unsigned int lambda,
                               unsigned int out_len_bytes, uint8_t* q, unsigned int start,
                               unsigned int len, unsigned int col_idx) {
  unsigned int depth               = vec_com_rec->depth;
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;

  uint8_t* sd              = malloc(lambda_bytes);
  uint8_t* com             = malloc(lambda_bytes * 2);
  uint8_t* r               = malloc(out_len_bytes);
  unsigned int bit_offset  = (col_idx % 8);
  unsigned int byte_offset = (col_idx / 8);

  unsigned int offset = NumRec(depth, vec_com_rec->b);

  for (unsigned int i = 0; i < num_instances; i++) {
    unsigned int offset_index = i ^ offset;
    if (offset_index == 0) {
      continue;
    }
    extract_sd_com_rec(vec_com_rec, iv, lambda, i, sd, com);
    prg(sd, iv, r, lambda, out_len_bytes);

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
      q[base_idx] ^= offset_index << bit_offset;
      for (unsigned int j = 1; j < amount; j++) {
        q[base_idx + j] ^= offset_index >> (j * 8 - bit_offset);
      }
    }
  }
  free(sd);
  free(com);
  free(r);
}

void partial_vole_reconstruct_rmo(const uint8_t* iv, const uint8_t* chall,
                                  const uint8_t* const* pdec, const uint8_t* const* com_j,
                                  uint8_t* q, unsigned int ellhat, const faest_paramset_t* params,
                                  unsigned int start, unsigned int len) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int tau1         = params->faest_param.t1;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;

  unsigned int max_depth = MAX(k0, k1);
  vec_com_rec_t vec_com_rec;
  vec_com_rec.b       = malloc(max_depth * sizeof(uint8_t));
  vec_com_rec.nodes   = calloc(max_depth, lambda_bytes);
  vec_com_rec.com_j   = malloc(lambda_bytes * 2);
  uint8_t* tree_nodes = malloc(lambda_bytes * (max_depth - 1));

  memset(q, 0, len * lambda_bytes);

  unsigned int col_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;
    uint8_t chalout[MAX_DEPTH];
    ChalDec(chall, i, k0, tau0, k1, tau1, chalout);
    vector_reconstruction(pdec[i], com_j[i], chalout, lambda, depth, tree_nodes, &vec_com_rec);
    ReconstructVoleRMO(iv, &vec_com_rec, lambda, ellhat_bytes, q, start, len, col_idx);
    col_idx += depth;
  }

  free(vec_com_rec.b);
  free(vec_com_rec.nodes);
  free(vec_com_rec.com_j);
  free(tree_nodes);
}
