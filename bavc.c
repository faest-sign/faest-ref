/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bavc.h"
#include "random_oracle.h"
#include "compat.h"
#include "aes.h"
#include "instances.h"
#include "universal_hashing.h"

#include <string.h>

typedef struct tree_t {
  uint8_t* nodes;   /* The data for each node */
  size_t numNodes;  /* The total number of nodes in the tree */
  size_t numLeaves; /* The total number of leaves in the tree */
} tree_t;

#define NODE(nodes, node, lambda_bytes) (&nodes[(node) * (lambda_bytes)])

static void expand_seeds(uint8_t* nodes, const uint8_t* iv, const faest_paramset_t* params) {
  const unsigned int lambda_bytes = params->faest_param.lambda / 8;

  for (unsigned int alpha = 0; alpha < params->faest_param.L - 1; ++alpha) {
    // the nodes are located other in memory consecutively
    prg(NODE(nodes, alpha, lambda_bytes), iv, alpha, NODE(nodes, 2 * alpha + 1, lambda_bytes),
        params->faest_param.lambda, lambda_bytes * 2);
  }
}

static uint8_t* generate_seeds(const uint8_t* root_seed, const uint8_t* iv,
                               const faest_paramset_t* params) {
  unsigned int lambda_bytes = params->faest_param.lambda / 8;
  uint8_t* nodes            = calloc(2 * params->faest_param.L - 1, lambda_bytes);

  memcpy(NODE(nodes, 0, lambda_bytes), root_seed, lambda_bytes);
  expand_seeds(nodes, iv, params);

  return nodes;
}

/* Gets the bit string of a node according to its position in the binary tree */
/* idx -> 2 -> {0,1},, Little Endian */
int BitDec(unsigned int leafIndex, unsigned int depth, uint8_t* out) {
  if (leafIndex >= (1u << depth)) {
    return -1;
  }
  for (unsigned int j = 0; j < depth; j++, leafIndex /= 2) {
    out[j] = leafIndex % 2;
  }
  return 1;
}

unsigned int NumRec(unsigned int depth, const uint8_t* bi) {
  unsigned int out = 0;
  for (unsigned int i = 0; i < depth; i++) {
    out += ((unsigned int)bi[i]) << i;
  }
  return out;
}

// FAEST.LeafCommit
static void faest_leaf_commit(uint8_t* sd, uint8_t* com, const uint8_t* key, const uint8_t* iv,
                              uint32_t tweak, const uint8_t* uhash, unsigned int lambda) {
  const unsigned int lambda_bytes = lambda / 8;

  uint8_t buffer[MAX_LAMBDA_BYTES * 4];
  prg(key, iv, tweak, buffer, lambda, lambda_bytes * 4);
  leaf_hash(com, uhash, buffer, lambda);
  memcpy(sd, buffer, lambda_bytes);
}

// FAEST-EM.LeafCommit
static void faest_em_leaf_commit(uint8_t* sd, uint8_t* com, const uint8_t* key, const uint8_t* iv,
                                 uint32_t tweak, unsigned int lambda) {
  const unsigned int lambda_bytes = lambda / 8;

  memcpy(sd, key, lambda_bytes);
  prg(key, iv, tweak, com, lambda, lambda_bytes * 2);
}

// BAVC.PosInTree
static inline unsigned int pos_in_tree(unsigned int i, unsigned int j,
                                       const faest_paramset_t* params) {
  const unsigned int tmp = 1 << (params->faest_param.k - 1);
  if (j < tmp) {
    return params->faest_param.L - 1 + params->faest_param.tau * j + i;
  }
  // mod 2^(k-1) is the same as & 2^(k-1)-1
  const unsigned int mask = tmp - 1;
  return params->faest_param.L - 1 + params->faest_param.tau * tmp +
         params->faest_param.tau1 * (j & mask) + i;
}

// BAVC.Commit for FAEST
static void bavc_commit_faest(const uint8_t* rootKey, const uint8_t* iv,
                              const faest_paramset_t* params, vec_com_t* vecCom) {
  const unsigned int lambda       = params->faest_param.lambda;
  const unsigned int L            = params->faest_param.L;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int com_size     = lambda_bytes * 3; // size of com_ij

  H0_context_t uhash_ctx;
  H0_init(&uhash_ctx, lambda);
  H0_update(&uhash_ctx, iv, IV_SIZE);
  H0_final_for_squeeze(&uhash_ctx);

  H1_context_t h1_com_ctx;
  H1_init(&h1_com_ctx, lambda);

  // Generating the tree (k)
  uint8_t* nodes = generate_seeds(rootKey, iv, params);

  // Initialzing stuff
  vecCom->h   = malloc(lambda_bytes * 2);
  vecCom->com = malloc(L * com_size);
  vecCom->sd  = malloc(L * lambda_bytes);

  // Step: 1..3
  vecCom->k = NODE(nodes, 0, lambda_bytes);

  // Step: 4..5
  // compute commitments for remaining instances
  for (unsigned int i = 0, offset = 0; i < params->faest_param.tau; ++i) {
    uint8_t uhash[MAX_LAMBDA_BYTES * 3];
    H0_squeeze(&uhash_ctx, uhash, 3 * lambda_bytes);

    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda);

    const unsigned int N_i =
        bavc_max_node_index(i, params->faest_param.tau1, params->faest_param.k);
    for (unsigned int j = 0; j < N_i; ++j, ++offset) {
      const unsigned int alpha = pos_in_tree(i, j, params);
      faest_leaf_commit(vecCom->sd + offset * lambda_bytes, vecCom->com + offset * com_size,
                        NODE(nodes, alpha, lambda_bytes), iv, i + L - 1, uhash, lambda);
      H1_update(&h1_ctx, vecCom->com + offset * com_size, com_size);
    }

    uint8_t hi[MAX_LAMBDA_BYTES * 2];
    // Step 11
    H1_final(&h1_ctx, hi, lambda_bytes * 2);
    // Step 12
    H1_update(&h1_com_ctx, hi, lambda_bytes * 2);
  }
  H0_clear(&uhash_ctx);

  // Step 12
  H1_final(&h1_com_ctx, vecCom->h, lambda_bytes * 2);
}

// BAVC.Commit for FAEST-EM
static void bavc_commit_faest_em(const uint8_t* rootKey, const uint8_t* iv,
                                 const faest_paramset_t* params, vec_com_t* vecCom) {
  const unsigned int lambda       = params->faest_param.lambda;
  const unsigned int L            = params->faest_param.L;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int com_size     = lambda_bytes * 2; // size of com_ij

  H1_context_t h1_com_ctx;
  H1_init(&h1_com_ctx, lambda);

  // Generating the tree
  uint8_t* nodes = generate_seeds(rootKey, iv, params);

  // Initialzing stuff
  vecCom->h   = malloc(lambda_bytes * 2);
  vecCom->com = malloc(L * com_size);
  vecCom->sd  = malloc(L * lambda_bytes);

  // Step: 1..3
  vecCom->k = NODE(nodes, 0, lambda_bytes);

  // Step: 4..5
  // compute commitments for remaining instances
  for (unsigned int i = 0, offset = 0; i < params->faest_param.tau; ++i) {
    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda);

    const unsigned int N_i =
        bavc_max_node_index(i, params->faest_param.tau1, params->faest_param.k);
    for (unsigned int j = 0; j < N_i; ++j, ++offset) {
      const unsigned int alpha = pos_in_tree(i, j, params);
      faest_em_leaf_commit(vecCom->sd + offset * lambda_bytes, vecCom->com + offset * com_size,
                           NODE(nodes, alpha, lambda_bytes), iv, i + L - 1, lambda);
      H1_update(&h1_ctx, vecCom->com + offset * com_size, com_size);
    }

    uint8_t hi[MAX_LAMBDA_BYTES * 2];
    // Step 11
    H1_final(&h1_ctx, hi, lambda_bytes * 2);
    // Step 12
    H1_update(&h1_com_ctx, hi, lambda_bytes * 2);
  }

  // Step 12
  H1_final(&h1_com_ctx, vecCom->h, lambda_bytes * 2);
}

void bavc_commit(const uint8_t* rootKey, const uint8_t* iv, const faest_paramset_t* params,
                 vec_com_t* vecCom) {
  if (faest_is_em(params)) {
    bavc_commit_faest_em(rootKey, iv, params, vecCom);
  } else {
    bavc_commit_faest(rootKey, iv, params, vecCom);
  }
}

bool bavc_open(const vec_com_t* vc, const uint16_t* i_delta, uint8_t* decom_i,
               const faest_paramset_t* params) {
  const unsigned int lambda       = params->faest_param.lambda;
  const unsigned int L            = params->faest_param.L;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int k            = params->faest_param.k;
  const unsigned int tau          = params->faest_param.tau;
  const unsigned int tau_1        = params->faest_param.tau1;
  const unsigned int com_size     = faest_is_em(params) ? (2 * lambda_bytes) : (3 * lambda_bytes);

  uint8_t* decom_i_end = decom_i + com_size * tau + params->faest_param.T_open * lambda_bytes;

  // Step 5
  uint8_t* s = calloc((2 * L - 1 + 7) / 8, 1);
  // Step 6
  unsigned int nh = 0;

  // Step 7..15
  for (unsigned int i = 0; i < tau; ++i) {
    unsigned int alpha = pos_in_tree(i, i_delta[i], params);
    ptr_set_bit(s, 1, alpha);
    ++nh;

    while (alpha > 0 && ptr_get_bit(s, (alpha - 1) / 2) == 0) {
      alpha = (alpha - 1) / 2;
      ptr_set_bit(s, 1, alpha);
      ++nh;
    }
  }

  // Step 16..17
  if (nh - 2 * tau + 1 > params->faest_param.T_open) {
    free(s);
    return false;
  }

  // Step 3
  const uint8_t* com = vc->com;
  for (unsigned int i = 0; i < tau; ++i) {
    memcpy(decom_i, com + i_delta[i] * com_size, com_size);
    com += bavc_max_node_index(i, tau_1, k) * com_size;
    decom_i += com_size;
  }

  // Step 19..25
  for (int i = L - 2; i >= 0; --i) {
    ptr_set_bit(s, ptr_get_bit(s, 2 * i + 1) | ptr_get_bit(s, 2 * i + 2), i);
    if ((ptr_get_bit(s, 2 * i + 1) ^ ptr_get_bit(s, 2 * i + 2)) == 1) {
      const unsigned int alpha = 2 * i + 1 + ptr_get_bit(s, 2 * i + 1);
      memcpy(decom_i, NODE(vc->k, alpha, lambda_bytes), lambda_bytes);
      decom_i += lambda_bytes;
    }
  }

  memset(decom_i, 0, decom_i_end - decom_i);

  free(s);
  return true;
}

static bool reconstruct_keys(uint8_t* s, uint8_t* keys, const uint8_t* decom_i,
                             const uint16_t* i_delta, const uint8_t* iv,
                             const faest_paramset_t* params) {
  const unsigned int lambda       = params->faest_param.lambda;
  const unsigned int L            = params->faest_param.L;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int tau          = params->faest_param.tau;

  const uint8_t* nodes = decom_i + (faest_is_em(params) ? 2 : 3) * tau * lambda_bytes;
  const uint8_t* end   = nodes + params->faest_param.T_open * lambda_bytes;

  // Step 7..10
  for (unsigned int i = 0; i < tau; ++i) {
    unsigned int alpha = pos_in_tree(i, i_delta[i], params);
    ptr_set_bit(s, 1, alpha);
  }

  // Step 12.12
  for (int i = L - 2; i >= 0; --i) {
    ptr_set_bit(s, ptr_get_bit(s, 2 * i + 1) | ptr_get_bit(s, 2 * i + 2), i);
    if ((ptr_get_bit(s, 2 * i + 1) ^ ptr_get_bit(s, 2 * i + 2)) == 1) {
      if (nodes == end) {
        return false;
      }

      const unsigned int alpha = 2 * i + 1 + ptr_get_bit(s, 2 * i + 1);
      memcpy(keys + alpha * lambda_bytes, nodes, lambda_bytes);
      nodes += lambda_bytes;
    }
  }

  for (; nodes != end; ++nodes) {
    if (*nodes) {
      return false;
    }
  }

  for (unsigned int i = 0; i != L - 1; ++i) {
    if (!ptr_get_bit(s, i)) {
      prg(keys + i * lambda_bytes, iv, i, keys + (2 * i + 1) * lambda_bytes, lambda,
          2 * lambda_bytes);
    }
  }

  return true;
}

static bool bavc_reconstruct_faest(const uint8_t* decom_i, const uint16_t* i_delta,
                                   const uint8_t* iv, const faest_paramset_t* params,
                                   vec_com_rec_t* vecComRec) {
  // Initializing
  const unsigned int lambda       = params->faest_param.lambda;
  const unsigned int L            = params->faest_param.L;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int k            = params->faest_param.k;
  const unsigned int tau          = params->faest_param.tau;
  const unsigned int tau_1        = params->faest_param.tau1;
  const unsigned int com_size     = lambda_bytes * 3; // size of com_ij

  // Step 6
  uint8_t* s    = calloc((2 * L - 1 + 7) >> 3, 1);
  uint8_t* keys = calloc(2 * params->faest_param.L - 1, lambda_bytes);

  if (!reconstruct_keys(s, keys, decom_i, i_delta, iv, params)) {
    free(keys);
    free(s);
    return false;
  }

  H0_context_t uhash_ctx;
  H0_init(&uhash_ctx, lambda);
  H0_update(&uhash_ctx, iv, IV_SIZE);
  H0_final_for_squeeze(&uhash_ctx);

  H1_context_t h1_com_ctx;
  H1_init(&h1_com_ctx, lambda);

  for (unsigned int i = 0, offset = 0; i != tau; ++i) {
    uint8_t uhash[MAX_LAMBDA_BYTES * 3];
    H0_squeeze(&uhash_ctx, uhash, com_size);

    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda);

    const unsigned int N_i = bavc_max_node_index(i, tau_1, k);
    for (unsigned int j = 0; j != N_i; ++j) {
      const unsigned int alpha = pos_in_tree(i, j, params);
      if (ptr_get_bit(s, alpha)) {
        H1_update(&h1_ctx, decom_i + i * com_size, com_size);
      } else {
        uint8_t com[3 * MAX_LAMBDA_BYTES];
        faest_leaf_commit(vecComRec->s + offset * lambda_bytes, com, keys + alpha * lambda_bytes,
                          iv, i + L - 1, uhash, lambda);
        ++offset;
        H1_update(&h1_ctx, com, com_size);
      }
    }

    uint8_t hi[MAX_LAMBDA_BYTES * 2];
    H1_final(&h1_ctx, hi, lambda_bytes * 2);
    H1_update(&h1_com_ctx, hi, lambda_bytes * 2);
  }
  H0_clear(&uhash_ctx);

  H1_final(&h1_com_ctx, vecComRec->h, lambda_bytes * 2);

  free(keys);
  free(s);
  return true;
}

static bool bavc_reconstruct_faest_em(const uint8_t* decom_i, const uint16_t* i_delta,
                                      const uint8_t* iv, const faest_paramset_t* params,
                                      vec_com_rec_t* vecComRec) {
  // Initializing
  const unsigned int lambda       = params->faest_param.lambda;
  const unsigned int L            = params->faest_param.L;
  const unsigned int lambda_bytes = lambda / 8;
  const unsigned int k            = params->faest_param.k;
  const unsigned int tau          = params->faest_param.tau;
  const unsigned int tau_1        = params->faest_param.tau1;
  const unsigned int com_size     = lambda_bytes * 2; // size of com_ij

  // Step 6
  uint8_t* s    = calloc((2 * L - 1 + 7) / 8, 1);
  uint8_t* keys = calloc(2 * params->faest_param.L - 1, lambda_bytes);

  // Step 7..10
  if (!reconstruct_keys(s, keys, decom_i, i_delta, iv, params)) {
    free(keys);
    free(s);
    return false;
  }

  H1_context_t h1_com_ctx;
  H1_init(&h1_com_ctx, lambda);

  for (unsigned int i = 0, offset = 0; i != tau; ++i) {
    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda);

    const unsigned int N_i = bavc_max_node_index(i, tau_1, k);
    for (unsigned int j = 0; j != N_i; ++j) {
      const unsigned int alpha = pos_in_tree(i, j, params);
      if (ptr_get_bit(s, alpha)) {
        H1_update(&h1_ctx, decom_i + i * com_size, com_size);
      } else {
        uint8_t com[2 * MAX_LAMBDA_BYTES];
        faest_em_leaf_commit(vecComRec->s + offset * lambda_bytes, com, keys + alpha * lambda_bytes,
                             iv, i + L - 1, lambda);
        ++offset;
        H1_update(&h1_ctx, com, com_size);
      }
    }

    uint8_t hi[MAX_LAMBDA_BYTES * 2];
    H1_final(&h1_ctx, hi, lambda_bytes * 2);
    H1_update(&h1_com_ctx, hi, lambda_bytes * 2);
  }

  H1_final(&h1_com_ctx, vecComRec->h, lambda_bytes * 2);

  free(keys);
  free(s);
  return true;
}

bool bavc_reconstruct(const uint8_t* decom_i, const uint16_t* i_delta, const uint8_t* iv,
                      const faest_paramset_t* params, vec_com_rec_t* vecComRec) {
  return faest_is_em(params) ? bavc_reconstruct_faest_em(decom_i, i_delta, iv, params, vecComRec)
                             : bavc_reconstruct_faest(decom_i, i_delta, iv, params, vecComRec);
}

void vec_com_clear(vec_com_t* com) {
  free(com->sd);
  free(com->com);
  free(com->k);
  free(com->h);
}

void vec_com_rec_clear(vec_com_rec_t* rec) {
  free(rec->s);
  free(rec->h);
}
