/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vc.h"
#include "random_oracle.h"
#include "compat.h"
#include "aes.h"
#include "instances.h"

#include <assert.h>
#include <string.h>

typedef struct tree_t {
  uint8_t* nodes;   /* The data for each node */
  size_t numNodes;  /* The total number of nodes in the tree */
  size_t numLeaves; /* The total number of leaves in the tree */
} tree_t;

#define NODE(tree, node, lambda_bytes) (&(tree).nodes[(node) * (lambda_bytes)])

static ATTR_CONST int is_left_child(size_t node) {
  assert(node != 0);
  return (node % 2 == 1);
}

static tree_t create_tree(const faest_paramset_t* params, unsigned int depth) {
  tree_t tree;
  uint32_t lambdaBytes = params->faest_param.lambda / 8;

  tree.numNodes  = getBinaryTreeNodeCount(depth);
  tree.numLeaves = 1 << depth;
  tree.nodes     = calloc(tree.numNodes, lambdaBytes);

  return tree;
}

static ATTR_CONST size_t get_parent(size_t node) {
  assert(node != 0);

  if (is_left_child(node)) {
    return (node - 1) / 2;
  }
  return (node - 2) / 2;
}

static void expand_seeds(tree_t* tree, const uint8_t* iv, const faest_paramset_t* params) {
  const unsigned int lambda_bytes = params->faest_param.lambda / 8;

  /* Walk the tree, expanding seeds where possible. Compute children of
   * non-leaf nodes. */
  size_t lastNonLeaf = get_parent(tree->numNodes - 1);
  // for scan build
  assert(2 * lastNonLeaf + 2 < tree->numNodes);

  for (size_t i = 0; i <= lastNonLeaf; i++) {
    // the nodes are located other in memory consecutively
    prg(NODE(*tree, i, lambda_bytes), iv, NODE(*tree, 2 * i + 1, lambda_bytes),
        params->faest_param.lambda, lambda_bytes * 2);
  }
}

static tree_t generate_seeds(const uint8_t* rootSeed, const uint8_t* iv,
                             const faest_paramset_t* params, unsigned int depth) {
  uint32_t lambdaBytes = params->faest_param.lambda / 8;
  tree_t tree          = create_tree(params, depth);

  memcpy(NODE(tree, 0, lambdaBytes), rootSeed, lambdaBytes);
  expand_seeds(&tree, iv, params);

  return tree;
}

/* Gets how many nodes will be there in the tree in total including root node */
uint64_t getBinaryTreeNodeCount(unsigned int depth) {
  return (1 << (depth + 1)) - 1;
}

/* Calculates the flat array index of the binary tree position */
uint64_t getNodeIndex(uint64_t depth, uint64_t levelIndex) {
  if (depth == 0) {
    return 0;
  }
  return (((2 << (depth - 1)) - 2) + (levelIndex + 1));
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

void vector_commitment(const uint8_t* rootKey, const uint8_t* iv, const faest_paramset_t* params,
                       uint32_t lambda, vec_com_t* vecCom, uint32_t depth) {
  const unsigned int lambdaBytes      = lambda / 8;
  const unsigned int numVoleInstances = 1 << depth;

  // Generating the tree
  tree_t tree = generate_seeds(rootKey, iv, params, depth);

  // Initialzing stuff
  vecCom->h   = malloc(lambdaBytes * 2);
  vecCom->com = malloc(numVoleInstances * lambdaBytes * 2);
  vecCom->sd  = malloc(numVoleInstances * lambdaBytes);

  // Step: 1..3
  vecCom->k = NODE(tree, 0, lambdaBytes);

  // Step: 4..5
  const unsigned int base_index = tree.numNodes - tree.numLeaves;
  unsigned int i                = 0;
  // compute commitments for 4 instances in parallel
  for (; i < numVoleInstances / 4 * 4; i += 4) {
    H0_context_x4_t h0_ctx;
    H0_x4_init(&h0_ctx, lambda);
    H0_x4_update(&h0_ctx, NODE(tree, base_index + i, lambdaBytes),
                 NODE(tree, base_index + i + 1, lambdaBytes),
                 NODE(tree, base_index + i + 2, lambdaBytes),
                 NODE(tree, base_index + i + 3, lambdaBytes), lambdaBytes);
    H0_x4_update(&h0_ctx, iv, iv, iv, iv, IV_SIZE);
    H0_x4_final(&h0_ctx, vecCom->sd + i * lambdaBytes, vecCom->sd + (i + 1) * lambdaBytes,
                vecCom->sd + (i + 2) * lambdaBytes, vecCom->sd + (i + 3) * lambdaBytes, lambdaBytes,
                vecCom->com + i * lambdaBytes * 2, vecCom->com + (i + 1) * lambdaBytes * 2,
                vecCom->com + (i + 2) * lambdaBytes * 2, vecCom->com + (i + 3) * lambdaBytes * 2,
                (lambdaBytes * 2));
  }
  // compute commitments for remaining instances
  for (; i < numVoleInstances; i++) {
    H0_context_t h0_ctx;
    H0_init(&h0_ctx, lambda);
    H0_update(&h0_ctx, NODE(tree, base_index + i, lambdaBytes), lambdaBytes);
    H0_update(&h0_ctx, iv, 16);
    H0_final(&h0_ctx, vecCom->sd + (i * lambdaBytes), lambdaBytes,
             vecCom->com + (i * (lambdaBytes * 2)), (lambdaBytes * 2));
  }

  tree.nodes = NULL;

  // Step: 6
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  for (uint32_t j = 0; j < numVoleInstances; j++) {
    H1_update(&h1_ctx, vecCom->com + (j * (lambdaBytes * 2)), (lambdaBytes * 2));
  }
  H1_final(&h1_ctx, vecCom->h, lambdaBytes * 2);
}

void vector_open(const uint8_t* k, const uint8_t* com, const uint8_t* b, uint8_t* cop,
                 uint8_t* com_j, uint32_t depth, uint32_t lambdaBytes) {
  // Step: 1
  uint64_t leafIndex = NumRec(depth, b);

  // Step: 3..6
  uint32_t a = 0;
  for (uint32_t i = 0; i < depth; i++) {
    memcpy(cop + (lambdaBytes * i),
           k + (lambdaBytes * getNodeIndex(i + 1, (2 * a) + !b[depth - 1 - i])), lambdaBytes);
    a = (2 * a) + b[depth - 1 - i];
  }

  // Step: 7
  memcpy(com_j, com + (leafIndex * lambdaBytes * 2), lambdaBytes * 2);
}

void vector_reconstruction(const uint8_t* iv, const uint8_t* cop, const uint8_t* com_j,
                           const uint8_t* b, uint32_t lambda, uint32_t depth,
                           vec_com_rec_t* vecComRec) {
  // Initializing
  const unsigned int lambdaBytes      = lambda / 8;
  const unsigned int numVoleInstances = 1 << depth;
  const uint64_t leafIndex            = NumRec(depth, b);

  // Step: 3..9
  uint32_t a = 0;
  for (uint32_t i = 1; i <= depth; i++) {
    memcpy(vecComRec->k + (lambdaBytes * getNodeIndex(i, 2 * a + !b[depth - i])),
           cop + (lambdaBytes * (i - 1)), lambdaBytes);
    memset(vecComRec->k + (lambdaBytes * getNodeIndex(i, 2 * a + b[depth - i])), 0, lambdaBytes);

    const uint32_t current_depth = (1 << (i - 1));
    for (uint32_t j = 0; j < current_depth; j++) {
      if (j == a) {
        continue;
      }

      uint8_t out[2 * MAX_LAMBDA_BYTES];
      prg(vecComRec->k + (lambdaBytes * getNodeIndex(i - 1, j)), iv, out, lambda, lambdaBytes * 2);
      memcpy(vecComRec->k + (lambdaBytes * getNodeIndex(i, 2 * j)), out, lambdaBytes);
      memcpy(vecComRec->k + (lambdaBytes * getNodeIndex(i, (2 * j) + 1)), out + lambdaBytes,
             lambdaBytes);
    }

    a = a * 2 + b[depth - i];
  }

  // Step: 10..11
  unsigned int j = 0;
  // reconstruct commitments for 4 instances in parallel
  for (; j < leafIndex / 4 * 4; j += 4) {
    H0_context_x4_t h0_ctx;
    H0_x4_init(&h0_ctx, lambda);
    H0_x4_update(&h0_ctx, vecComRec->k + (getNodeIndex(depth, j) * lambdaBytes),
                 vecComRec->k + (getNodeIndex(depth, j + 1) * lambdaBytes),
                 vecComRec->k + (getNodeIndex(depth, j + 2) * lambdaBytes),
                 vecComRec->k + (getNodeIndex(depth, j + 3) * lambdaBytes), lambdaBytes);
    H0_x4_update(&h0_ctx, iv, iv, iv, iv, IV_SIZE);
    H0_x4_final(&h0_ctx, vecComRec->s + j * lambdaBytes, vecComRec->s + (j + 1) * lambdaBytes,
                vecComRec->s + (j + 2) * lambdaBytes, vecComRec->s + (j + 3) * lambdaBytes,
                lambdaBytes, vecComRec->com + j * lambdaBytes * 2,
                vecComRec->com + (j + 1) * lambdaBytes * 2,
                vecComRec->com + (j + 2) * lambdaBytes * 2,
                vecComRec->com + (j + 3) * lambdaBytes * 2, lambdaBytes * 2);
  }
  // reconstruct commitments up until the leafIndex
  for (; j < leafIndex; ++j) {
    H0_context_t h0_ctx;
    H0_init(&h0_ctx, lambda);
    H0_update(&h0_ctx, vecComRec->k + getNodeIndex(depth, j) * lambdaBytes, lambdaBytes);
    H0_update(&h0_ctx, iv, IV_SIZE);
    H0_final(&h0_ctx, vecComRec->s + j * lambdaBytes, lambdaBytes,
             vecComRec->com + j * lambdaBytes * 2, lambdaBytes * 2);
  }
  // skip leafIndex
  ++j;
  // reconstruct until index is divisible by 4 again
  for (; j < numVoleInstances && j % 4; ++j) {
    H0_context_t h0_ctx;
    H0_init(&h0_ctx, lambda);
    H0_update(&h0_ctx, vecComRec->k + getNodeIndex(depth, j) * lambdaBytes, lambdaBytes);
    H0_update(&h0_ctx, iv, IV_SIZE);
    H0_final(&h0_ctx, vecComRec->s + j * lambdaBytes, lambdaBytes,
             vecComRec->com + j * lambdaBytes * 2, lambdaBytes * 2);
  }
  // reconstruct 4 instances in parallel
  for (; j < numVoleInstances / 4 * 4; j += 4) {
    H0_context_x4_t h0_ctx;
    H0_x4_init(&h0_ctx, lambda);
    H0_x4_update(&h0_ctx, vecComRec->k + (getNodeIndex(depth, j) * lambdaBytes),
                 vecComRec->k + (getNodeIndex(depth, j + 1) * lambdaBytes),
                 vecComRec->k + (getNodeIndex(depth, j + 2) * lambdaBytes),
                 vecComRec->k + (getNodeIndex(depth, j + 3) * lambdaBytes), lambdaBytes);
    H0_x4_update(&h0_ctx, iv, iv, iv, iv, IV_SIZE);
    H0_x4_final(&h0_ctx, vecComRec->s + j * lambdaBytes, vecComRec->s + (j + 1) * lambdaBytes,
                vecComRec->s + (j + 2) * lambdaBytes, vecComRec->s + (j + 3) * lambdaBytes,
                lambdaBytes, vecComRec->com + j * lambdaBytes * 2,
                vecComRec->com + (j + 1) * lambdaBytes * 2,
                vecComRec->com + (j + 2) * lambdaBytes * 2,
                vecComRec->com + (j + 3) * lambdaBytes * 2, lambdaBytes * 2);
  }
  // reconstruct remaining instances
  for (; j < numVoleInstances; ++j) {
    H0_context_t h0_ctx;
    H0_init(&h0_ctx, lambda);
    H0_update(&h0_ctx, vecComRec->k + getNodeIndex(depth, j) * lambdaBytes, lambdaBytes);
    H0_update(&h0_ctx, iv, IV_SIZE);
    H0_final(&h0_ctx, vecComRec->s + j * lambdaBytes, lambdaBytes,
             vecComRec->com + j * lambdaBytes * 2, lambdaBytes * 2);
  }

  // Step: 12..13
  memcpy(vecComRec->com + (lambdaBytes * 2 * leafIndex), com_j, lambdaBytes * 2);
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  H1_update(&h1_ctx, vecComRec->com, lambdaBytes * 2 * numVoleInstances);
  H1_final(&h1_ctx, vecComRec->h, lambdaBytes * 2);
}

#if defined(FAEST_TESTS)
int vector_verify(const uint8_t* iv, const uint8_t* pdec, const uint8_t* com_j, const uint8_t* b,
                  uint32_t lambda, uint32_t depth, vec_com_rec_t* rec, const uint8_t* vecComH) {
  const unsigned int lambdaBytes      = lambda / 8;
  const unsigned int numVoleInstances = 1 << depth;

  vec_com_rec_t vecComRec;
  vecComRec.h   = malloc(lambdaBytes * 2);
  vecComRec.k   = calloc(getBinaryTreeNodeCount(depth), lambdaBytes);
  vecComRec.com = malloc(numVoleInstances * lambdaBytes * 2);
  vecComRec.s   = malloc(numVoleInstances * lambdaBytes);

  // Step: 2
  vector_reconstruction(iv, pdec, com_j, b, lambda, depth, &vecComRec);

  // Step: 3
  int ret = memcmp(vecComH, vecComRec.h, lambdaBytes * 2);
  if (!rec || ret) {
    vec_com_rec_clear(&vecComRec);
  }

  if (ret == 0) {
    if (rec) {
      *rec = vecComRec;
    }
    return 1;
  } else {
    return 0;
  }
}
#endif

void vec_com_clear(vec_com_t* com) {
  free(com->sd);
  free(com->com);
  free(com->k);
  free(com->h);
}

void vec_com_rec_clear(vec_com_rec_t* rec) {
  free(rec->s);
  free(rec->com);
  free(rec->k);
  free(rec->h);
}
