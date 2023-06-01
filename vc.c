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

#include <assert.h>
#include <string.h>

#define MAX_SEED_SIZE (256 / 8)

typedef struct tree_t {
  size_t depth;      /* The depth of the tree */
  uint8_t** nodes;   /* The data for each node */
  uint8_t* haveNode; /* If we have the data (seed or hash) for node i, haveSeed[i] is 1 */
  uint8_t* exists;   /* Since the tree is not always complete, nodes marked 0 don't exist */
  size_t numNodes;   /* The total number of nodes in the tree */
  size_t numLeaves;  /* The total number of leaves in the tree */
} tree_t;

static int exists(tree_t* tree, size_t i) {
  if (i >= tree->numNodes) {
    return 0;
  }
  if (tree->exists[i]) {
    return 1;
  }
  return 0;
}

static uint8_t** getLeaves(tree_t* tree) {
  return &tree->nodes[tree->numNodes - tree->numLeaves];
}

static int isLeftChild(size_t node) {
  assert(node != 0);
  return (node % 2 == 1);
}

static tree_t createTree(const faest_paramset_t* params, uint32_t num_nodes) {
  assert(num_nodes > 0);
  tree_t tree;
  uint32_t lambdaBytes = params->faest_param.lambda / 8;

  tree.depth = ceil_log2(num_nodes) + 1;
  tree.numNodes =
      ((1 << (tree.depth)) - 1) -
      ((1 << (tree.depth - 1)) - num_nodes); /* Num nodes in complete - number of missing leaves */
  tree.numLeaves = num_nodes;
  tree.nodes     = malloc(tree.numNodes * sizeof(uint8_t*));

  uint8_t* slab = calloc(tree.numNodes, lambdaBytes);

  for (size_t i = 0; i < tree.numNodes; i++) {
    tree.nodes[i] = slab;
    slab += lambdaBytes;
  }

  tree.haveNode = calloc(tree.numNodes, 1);

  /* Depending on the number of leaves, the tree may not be complete */
  tree.exists = calloc(tree.numNodes, 1);
  memset(tree.exists + tree.numNodes - tree.numLeaves, 1, tree.numLeaves); /* Set leaves */
  for (int i = tree.numNodes - tree.numLeaves; i > 0; i--) {
    if (exists(&tree, 2 * i + 1) || exists(&tree, 2 * i + 2)) {
      tree.exists[i] = 1;
    }
  }
  tree.exists[0] = 1;

  return tree;
}

static void freeTree(tree_t* tree) {
  if (tree != NULL) {
    free(tree->nodes[0]);
    free(tree->nodes);
    free(tree->haveNode);
    free(tree->exists);
  }
}

static size_t getParent(size_t node) {
  assert(node != 0);

  if (isLeftChild(node)) {
    return (node - 1) / 2;
  }
  return (node - 2) / 2;
}

static void expandSeeds(tree_t* tree, const uint8_t* iv, const faest_paramset_t* params) {
  uint32_t lambdaBytes = params->faest_param.lambda / 8;

  uint8_t out[2 * MAX_SEED_SIZE];
  assert(lambdaBytes <= sizeof(out));

  /* Walk the tree, expanding seeds where possible. Compute children of
   * non-leaf nodes. */
  size_t lastNonLeaf = getParent(tree->numNodes - 1);

  for (size_t i = 0; i <= lastNonLeaf; i++) {
    if (!tree->haveNode[i]) {
      continue;
    }

    prg(tree->nodes[i], iv, out, params->faest_param.lambda, params->faest_param.lambda / 4);

    if (!tree->haveNode[2 * i + 1]) {
      memcpy(tree->nodes[2 * i + 1], out, lambdaBytes);
      tree->haveNode[2 * i + 1] = 1;
    }

    /* The last non-leaf node will only have a left child when there are an odd number of leaves */
    if (exists(tree, 2 * i + 2) && !tree->haveNode[2 * i + 2]) {
      memcpy(tree->nodes[2 * i + 2], out + lambdaBytes, lambdaBytes);
      tree->haveNode[2 * i + 2] = 1;
    }
  }
}

static tree_t generateSeeds(const uint8_t* rootSeed, const uint8_t* iv,
                            const faest_paramset_t* params, uint32_t numVoleInstances) {

  uint32_t lambdaBytes = params->faest_param.lambda / 8;
  tree_t tree          = createTree(params, numVoleInstances);

  memcpy(tree.nodes[0], rootSeed, lambdaBytes);
  tree.haveNode[0] = 1;
  expandSeeds(&tree, iv, params);

#if 0
  printTree("tree", tree);
  printTreeInfo("tree_info", tree);
#endif

  return tree;
}

/* Gets how many nodes will be there in the tree in total including root node */
uint64_t getBinaryTreeNodeCount(uint32_t numVoleInstances) {
  uint32_t depth = ceil_log2(numVoleInstances) + 1;
  return ((1 << depth) - 1) - ((1 << (depth - 1)) - numVoleInstances);

  // uint64_t out = 0;
  // for (uint64_t i = depth; i > 0; i--) {
  //   out += (1 << i);
  // }
  // out += 1;
  // return out;
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
int BitDec(uint32_t leafIndex, uint32_t depth, uint8_t* out) {
  uint32_t i = leafIndex;
  if (leafIndex >= (uint32_t)(1 << depth)) {
    return -1;
  }
  for (uint32_t j = 0; j < depth; j++) {
    out[j] = i % 2;
    i      = (i - out[j]) / 2;
  }
  return 1;
}

uint64_t NumRec(uint32_t depth, const uint8_t* bi) {
  uint64_t out = 0;
  for (uint32_t i = 0; i < depth; i++) {
    out = out + ((uint64_t)bi[i] << i);
  }
  return out;
}

void vector_commitment(const uint8_t* rootKey, const uint8_t* iv, const faest_paramset_t* params,
                       uint32_t lambda, uint32_t lambdaBytes, vec_com_t* vecCom,
                       uint32_t numVoleInstances) {

  // Generating the tree
  tree_t tree = generateSeeds(rootKey, iv, params, numVoleInstances);

  // Initialzing stuff
  vecCom->h   = malloc(lambdaBytes * 2);
  vecCom->k   = malloc(tree.numNodes * lambdaBytes);
  vecCom->com = malloc(numVoleInstances * lambdaBytes * 2);
  vecCom->sd  = malloc(numVoleInstances * lambdaBytes);

  // Step: 1..3
  for (uint32_t i = 0; i < tree.numNodes; i++) {
    memcpy(vecCom->k + (i * lambdaBytes), tree.nodes[i], lambdaBytes);
  }

  // Step: 4..5
  uint8_t** leaves = getLeaves(&tree);
  for (uint32_t i = 0; i < numVoleInstances; i++) {
    H0_context_t h0_ctx;
    H0_init(&h0_ctx, lambda);
    H0_update(&h0_ctx, leaves[i], lambdaBytes);
    H0_update(&h0_ctx, iv, 16);
    H0_final(&h0_ctx, vecCom->sd + (i * lambdaBytes), lambdaBytes,
             vecCom->com + (i * (lambdaBytes * 2)), (lambdaBytes * 2));
  }
  freeTree(&tree);

  // Step: 6
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  for (uint32_t i = 0; i < numVoleInstances; i++) {
    H1_update(&h1_ctx, vecCom->com + (i * (lambdaBytes * 2)), (lambdaBytes * 2));
  }
  H1_final(&h1_ctx, vecCom->h, lambdaBytes * 2);
  /* printHex("hashed_leaves", vecCom->com, 1 << 10); */
}

void vector_open(const uint8_t* k, const uint8_t* com, const uint8_t* b, uint8_t* cop,
                 uint8_t* com_j, uint32_t numVoleInstances, uint32_t lambdaBytes) {
  uint32_t depth = ceil_log2(numVoleInstances);

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
                           const uint8_t* b, uint32_t lambda, uint32_t lambdaBytes,
                           uint32_t numVoleInstances, uint32_t depth, vec_com_rec_t* vecComRec) {
  // Initializing
  uint64_t leafIndex = NumRec(depth, b);

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

      uint8_t out[2 * MAX_SEED_SIZE];
      prg(vecComRec->k + (lambdaBytes * getNodeIndex(i - 1, j)), iv, out, lambda, lambdaBytes * 2);
      memcpy(vecComRec->k + (lambdaBytes * getNodeIndex(i, 2 * j)), out, lambdaBytes);
      memcpy(vecComRec->k + (lambdaBytes * getNodeIndex(i, (2 * j) + 1)), out + lambdaBytes,
             lambdaBytes);
    }

    a = a * 2 + b[depth - i];
  }

  // Step: 10..11
  for (uint32_t j = 0; j < numVoleInstances; j++) {
    /* Reconstruct the coms and the m from the ks while keeping k_j* secret */
    if (j != leafIndex) {
      H0_context_t h0_ctx;
      H0_init(&h0_ctx, lambda);
      H0_update(&h0_ctx, vecComRec->k + (getNodeIndex(depth, j) * lambdaBytes), lambdaBytes);
      H0_update(&h0_ctx, iv, 16);
      H0_final(&h0_ctx, vecComRec->s + (lambdaBytes * j), lambdaBytes,
               vecComRec->com + (lambdaBytes * 2 * j), lambdaBytes * 2);
    }
  }
  // Step: 12..13
  memcpy(vecComRec->com + (lambdaBytes * 2 * leafIndex), com_j, lambdaBytes * 2);
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  H1_update(&h1_ctx, vecComRec->com, lambdaBytes * 2 * numVoleInstances);
  H1_final(&h1_ctx, vecComRec->h, lambdaBytes * 2);
}

int vector_verify(const uint8_t* iv, const uint8_t* pdec, const uint8_t* com_j, const uint8_t* b,
                  uint32_t lambda, uint32_t lambdaBytes, uint32_t numVoleInstances, uint32_t depth,
                  vec_com_rec_t* rec, const uint8_t* vecComH) {
  vec_com_rec_t vecComRec;
  vecComRec.h   = malloc(lambdaBytes * 2);
  vecComRec.k   = calloc(getBinaryTreeNodeCount(numVoleInstances), lambdaBytes);
  vecComRec.com = malloc(numVoleInstances * lambdaBytes * 2);
  vecComRec.s   = malloc(numVoleInstances * lambdaBytes);

  // Step: 2
  vector_reconstruction(iv, pdec, com_j, b, lambda, lambdaBytes, numVoleInstances, depth,
                        &vecComRec);

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
