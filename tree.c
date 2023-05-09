#include <assert.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "hash_shake.h"
#include "tree.h"
#include "aes.h"
#include "utils.h"

static int contains(size_t* list, size_t len, size_t value) {
  for (size_t i = 0; i < len; i++) {
    if (list[i] == value) {
      return 1;
    }
  }
  return 0;
}

int exists(tree_t* tree, size_t i) {
  if (i >= tree->numNodes) {
    return 0;
  }
  if (tree->exists[i]) {
    return 1;
  }
  return 0;
}

tree_t* createTree(faest_paramset_t* params) {
  tree_t* tree = malloc(sizeof(tree_t));

  tree->depth    = ceil_log2(params->faest_param.t) + 1;
  tree->numNodes = ((1 << (tree->depth)) - 1) -
                   ((1 << (tree->depth - 1)) -
                    params->faest_param.t); /* Num nodes in complete - number of missing leaves */
  tree->numLeaves = params->faest_param.t;
  tree->dataSize  = params->faest_param.seedSizeBytes;
  tree->nodes     = malloc(tree->numNodes * sizeof(uint8_t*));

  uint8_t* slab = calloc(tree->numNodes, params->faest_param.seedSizeBytes);

  for (size_t i = 0; i < tree->numNodes; i++) {
    tree->nodes[i] = slab;
    slab += params->faest_param.seedSizeBytes;
  }

  tree->haveNode = calloc(tree->numNodes, 1);

  /* Depending on the number of leaves, the tree may not be complete */
  tree->exists = calloc(tree->numNodes, 1);
  memset(tree->exists + tree->numNodes - tree->numLeaves, 1, tree->numLeaves); /* Set leaves */
  for (int i = tree->numNodes - tree->numLeaves; i > 0; i--) {
    if (exists(tree, 2 * i + 1) || exists(tree, 2 * i + 2)) {
      tree->exists[i] = 1;
    }
  }
  tree->exists[0] = 1;

  return tree;
}

void freeTree(tree_t* tree) {
  if (tree != NULL) {
    free(tree->nodes[0]);
    free(tree->nodes);
    free(tree->haveNode);
    free(tree->exists);
    free(tree);
  }
}

int isLeftChild(size_t node) {
  assert(node != 0);
  return (node % 2 == 1);
}

int hasRightChild(tree_t* tree, size_t node) {
  return (2 * node + 2 < tree->numNodes && exists(tree, node));
}

int hasLeftChild(tree_t* tree, size_t node) {
  return (2 * node + 1 < tree->numNodes);
}

size_t getParent(size_t node) {
  assert(node != 0);

  if (isLeftChild(node)) {
    return (node - 1) / 2;
  }
  return (node - 2) / 2;
}

uint8_t** getLeaves(tree_t* tree) {
  return &tree->nodes[tree->numNodes - tree->numLeaves];
}

uint8_t* getLeaf(tree_t* tree, size_t leafIndex) {
  assert(leafIndex < tree->numLeaves);
  size_t firstLeaf = tree->numNodes - tree->numLeaves;
  return tree->nodes[firstLeaf + leafIndex];
}

void hashSeed(uint8_t* digest, const uint8_t* inputSeed, uint8_t* salt, size_t repIndex,
              size_t nodeIndex, faest_paramset_t* params) {
  hash_context ctx;

  hash_init(&ctx, params->cipher_param.stateSizeBits);
  hash_update(&ctx, inputSeed, params->faest_param.seedSizeBytes);
  hash_update(&ctx, salt, params->faest_param.saltSizeBytes);
  hash_update_uint16_le(&ctx, (uint16_t)repIndex);
  hash_update_uint16_le(&ctx, (uint16_t)nodeIndex);
  hash_final(&ctx);
  hash_squeeze(&ctx, digest, params->faest_param.seedSizeBytes);

  // HashInit(&ctx, params, hashPrefix);
  // HashUpdate(&ctx, inputSeed, params->seedSizeBytes);
  // HashUpdate(&ctx, salt, params->saltSizeBytes);
  // HashUpdateIntLE(&ctx, (uint16_t)repIndex);
  // HashUpdateIntLE(&ctx, (uint16_t)nodeIndex);
  // HashFinal(&ctx);
  // HashSqueeze(&ctx, digest, 2 * params->seedSizeBytes);
}

void expandSeeds(tree_t* tree, faest_paramset_t* params) {
  uint8_t out[2 * MAX_SEED_SIZE_BYTES];

  /* Walk the tree, expanding seeds where possible. Compute children of
   * non-leaf nodes. */
  size_t lastNonLeaf = getParent(tree->numNodes - 1);

  for (size_t i = 0; i <= lastNonLeaf; i++) {
    if (!tree->haveNode[i]) {
      continue;
    }

    // Always statring with zeor IV ?
    uint8_t iv[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

    // Here we use the AES ctr PRG to get the nodes, starting from root and
    // assign it to the tree
    aes_prg(tree->nodes[i], iv, out, params->faest_param.seclvl);
    if (!tree->haveNode[2 * i + 1]) {
      /* left child = H_left(seed_i || salt || t || i) */
      memcpy(tree->nodes[2 * i + 1], out, params->faest_param.seedSizeBytes);
      tree->haveNode[2 * i + 1] = 1;
    }

    aes_prg(tree->nodes[i], iv, out, params->faest_param.seclvl);
    /* The last non-leaf node will only have a left child when there are an odd number of leaves */
    if (exists(tree, 2 * i + 2) && !tree->haveNode[2 * i + 2]) {
      /* right child = H_right(seed_i || salt || t || i)  */
      memcpy(tree->nodes[2 * i + 2], out + params->faest_param.seedSizeBytes,
             params->faest_param.seedSizeBytes);
      tree->haveNode[2 * i + 2] = 1;
    }
  }
}

tree_t* generateSeeds(uint8_t* rootSeed, faest_paramset_t* params) {
  tree_t* tree = createTree(params);

  memcpy(tree->nodes[0], rootSeed, params->faest_param.seedSizeBytes);
  tree->haveNode[0] = 1;
  expandSeeds(tree, params);

  return tree;
}

int isLeafNode(tree_t* tree, size_t node) {
  return (2 * node + 1 >= tree->numNodes);
}

int hasSibling(tree_t* tree, size_t node) {
  if (!exists(tree, node)) {
    return 0;
  }

  if (isLeftChild(node) && !exists(tree, node + 1)) {
    return 0;
  }

  return 1;
}

size_t getSibling(tree_t* tree, size_t node) {
  assert(node < tree->numNodes);
  assert(node != 0);
  assert(hasSibling(tree, node));

  if (isLeftChild(node)) {
    if (node + 1 < tree->numNodes) {
      return node + 1;
    } else {
      assert(!"getSibling: request for node with not sibling");
      return 0;
    }
  } else {
    return node - 1;
  }
}

void printSeeds(uint8_t* seedsBuf, size_t seedLen, size_t numSeeds) {
  for (size_t i = 0; i < numSeeds; i++) {
    printf("seed %lu", i);
    printHex("", seedsBuf, seedLen);
    seedsBuf += seedLen;
  }
}

void printLeaves(tree_t* tree) {
  size_t firstLeaf = tree->numNodes - tree->numLeaves;

  printSeeds(tree->nodes[firstLeaf], tree->dataSize, tree->numLeaves);
}

/* Returns the number of bytes written to output */
static size_t* getRevealedNodes(tree_t* tree, uint16_t* hideList, size_t hideListSize,
                                size_t* outputSize) {
  /* Compute paths up from hideList to root, store as sets of nodes */
  size_t pathLen = tree->depth - 1;

  /* pathSets[i][0...hideListSize] stores the nodes in the path at depth i
   * for each of the leaf nodes in hideListSize */
  size_t** pathSets = malloc(pathLen * sizeof(size_t*));
  size_t* slab      = malloc(hideListSize * pathLen * sizeof(size_t));

  for (size_t i = 0; i < pathLen; i++) {
    pathSets[i] = slab;
    slab += hideListSize;
  }

  /* Compute the paths back to the root */
  for (size_t i = 0; i < hideListSize; i++) {
    size_t pos = 0;
    size_t node =
        hideList[i] +
        (tree->numNodes - tree->numLeaves); /* input lists leaf indexes, translate to nodes */
    pathSets[pos][i] = node;
    pos++;
    while ((node = getParent(node)) != 0) {
      pathSets[pos][i] = node;
      pos++;
    }
  }

  /* Determine seeds to reveal */
  size_t* revealed   = malloc(tree->numLeaves * sizeof(size_t));
  size_t revealedPos = 0;
  for (size_t d = 0; d < pathLen; d++) {
    for (size_t i = 0; i < hideListSize; i++) {
      if (!hasSibling(tree, pathSets[d][i])) {
        continue;
      }
      size_t sibling = getSibling(tree, pathSets[d][i]);
      if (!contains(pathSets[d], hideListSize, sibling)) {
        // Determine the seed to reveal
        while (!hasRightChild(tree, sibling) && !isLeafNode(tree, sibling)) {
          sibling = 2 * sibling + 1; // sibling = leftChild(sibling)
        }

        // Only reveal if we haven't already
        if (!contains(revealed, revealedPos, sibling)) {
          revealed[revealedPos] = sibling;
          revealedPos++;
        }
      }
    }
  }

  free(pathSets[0]);
  free(pathSets);

  *outputSize = revealedPos;
  return revealed;
}

size_t revealSeedsSize(uint16_t* hideList, size_t hideListSize, faest_paramset_t* params) {
  tree_t* tree            = createTree(params);
  size_t numNodesRevealed = 0;
  size_t* revealed        = getRevealedNodes(tree, hideList, hideListSize, &numNodesRevealed);

  freeTree(tree);
  free(revealed);
  return numNodesRevealed * params->faest_param.seedSizeBytes;
}

size_t revealSeeds(tree_t* tree, uint16_t* hideList, size_t hideListSize, uint8_t* output,
                   size_t outputSize, faest_paramset_t* params) {
  uint8_t* outputBase = output;
  size_t revealedSize = 0;

  if (outputSize > INT_MAX) {
    return -1;
  }
  int outLen = (int)outputSize;

  size_t* revealed = getRevealedNodes(tree, hideList, hideListSize, &revealedSize);
  for (size_t i = 0; i < revealedSize; i++) {
    outLen -= params->faest_param.seedSizeBytes;
    if (outLen < 0) {
      assert(!"Insufficient sized buffer provided to revealSeeds");
      free(revealed);
      return 0;
    }
    memcpy(output, tree->nodes[revealed[i]], params->faest_param.seedSizeBytes);
    output += params->faest_param.seedSizeBytes;
  }

  free(revealed);
  return output - outputBase;
}

int reconstructSeeds(tree_t* tree, uint16_t* hideList, size_t hideListSize, uint8_t* input,
                     size_t inputLen, faest_paramset_t* params) {
  int ret = 0;

  if (inputLen > INT_MAX) {
    return -1;
  }
  int inLen = (int)inputLen;

  size_t revealedSize = 0;
  size_t* revealed    = getRevealedNodes(tree, hideList, hideListSize, &revealedSize);
  for (size_t i = 0; i < revealedSize; i++) {
    inLen -= params->faest_param.seedSizeBytes;
    if (inLen < 0) {
      ret = -1;
      goto Exit;
    }
    memcpy(tree->nodes[revealed[i]], input, params->faest_param.seedSizeBytes);
    tree->haveNode[revealed[i]] = 1;
    input += params->faest_param.seedSizeBytes;
  }

  expandSeeds(tree, params);

Exit:
  free(revealed);
  return ret;
}

static void computeParentHash(tree_t* tree, size_t child, uint8_t* salt, faest_paramset_t* params) {
  if (!exists(tree, child)) {
    return;
  }

  size_t parent = getParent(child);

  if (tree->haveNode[parent]) {
    return;
  }

  /* Compute the hash for parent, if we have everything */
  if (!tree->haveNode[2 * parent + 1]) {
    return;
  }

  if (exists(tree, 2 * parent + 2) && !tree->haveNode[2 * parent + 2]) {
    return;
  }

  /* Compute parent data = H(left child data || [right child data] || salt || parent idx) */
  hash_context ctx;
  hash_init(&ctx, params->cipher_param.stateSizeBits);
  hash_update(&ctx, tree->nodes[2 * parent + 1], params->faest_param.h1digestSizeBytes);
  if (hasRightChild(tree, parent)) {
    /* One node may not have a right child when there's an odd number of leaves */
    hash_update(&ctx, tree->nodes[2 * parent + 2], params->faest_param.h1digestSizeBytes);
  }
  hash_update(&ctx, salt, params->faest_param.saltSizeBytes);
  hash_update_uint16_le(&ctx, (uint16_t)parent);
  hash_final(&ctx);
  hash_squeeze(&ctx, tree->nodes[parent], params->faest_param.h1digestSizeBytes);
  tree->haveNode[parent] = 1;

  // HashInit(&ctx, params, HASH_PREFIX_3);
  // HashUpdate(&ctx, tree->nodes[2 * parent + 1], params->digestSizeBytes);
  // if (hasRightChild(tree, parent)) {
  //     /* One node may not have a right child when there's an odd number of leaves */
  //     HashUpdate(&ctx, tree->nodes[2 * parent + 2], params->digestSizeBytes);
  // }

  // HashUpdate(&ctx, salt, params->saltSizeBytes);
  // HashUpdateIntLE(&ctx, (uint16_t)parent);
  // HashFinal(&ctx);
  // HashSqueeze(&ctx, tree->nodes[parent], params->digestSizeBytes);
  // tree->haveNode[parent] = 1;
}

/* Create a Merkle tree by hashing up all nodes.
 * leafData must have length tree->numNodes, but some may be NULL. */
void buildMerkleTree(tree_t* tree, uint8_t** leafData, uint8_t* salt, faest_paramset_t* params) {
  size_t firstLeaf = tree->numNodes - tree->numLeaves;

  /* Copy data to the leaves. The actual data being committed to has already been
   * hashed, according to the spec. */
  for (size_t i = 0; i < tree->numLeaves; i++) {
    if (leafData[i] != NULL) {
      memcpy(tree->nodes[firstLeaf + i], leafData[i], tree->dataSize);
      tree->haveNode[firstLeaf + i] = 1;
    }
  }
  /* Starting at the leaves, work up the tree, computing the hashes for intermediate nodes */
  for (int i = (int)tree->numNodes; i > 0; i--) {
    computeParentHash(tree, i, salt, params);
  }
}

/* Note that we never output the root node */
static size_t* getRevealedMerkleNodes(tree_t* tree, uint16_t* missingLeaves,
                                      size_t missingLeavesSize, size_t* outputSize) {
  size_t firstLeaf      = tree->numNodes - tree->numLeaves;
  uint8_t* missingNodes = calloc(tree->numNodes, 1);

  /* Mark leaves that are missing */
  for (size_t i = 0; i < missingLeavesSize; i++) {
    missingNodes[firstLeaf + missingLeaves[i]] = 1;
  }

  /* For the nonleaf nodes, if both leaves are missing, mark it as missing too */
  int lastNonLeaf = getParent(tree->numNodes - 1);
  for (int i = lastNonLeaf; i > 0; i--) {
    if (!exists(tree, i)) {
      continue;
    }
    if (exists(tree, 2 * i + 2)) {
      if (missingNodes[2 * i + 1] && missingNodes[2 * i + 2]) {
        missingNodes[i] = 1;
      }
    } else {
      if (missingNodes[2 * i + 1]) {
        missingNodes[i] = 1;
      }
    }
  }

  /* For each missing leaf node, add the highest missing node on the path
   * back to the root to the set to be revealed */
  size_t* revealed = malloc(tree->numLeaves * sizeof(size_t));
  size_t pos       = 0;
  for (size_t i = 0; i < missingLeavesSize; i++) {
    size_t node = missingLeaves[i] + firstLeaf; /* input is leaf indexes, translate to nodes */
    do {
      if (!missingNodes[getParent(node)]) {
        if (!contains(revealed, pos, node)) {
          revealed[pos] = node;
          pos++;
        }
        break;
      }
    } while ((node = getParent(node)) != 0);
  }

  free(missingNodes);
  *outputSize = pos;
  return revealed;
}

size_t openMerkleTreeSize(uint16_t* missingLeaves, size_t missingLeavesSize,
                          faest_paramset_t* params) {

  tree_t* tree        = createTree(params);
  size_t revealedSize = 0;
  size_t* revealed = getRevealedMerkleNodes(tree, missingLeaves, missingLeavesSize, &revealedSize);

  freeTree(tree);
  free(revealed);

  return revealedSize * params->faest_param.h1digestSizeBytes;
}

/* Serialze the missing nodes that the verifier will require to check commitments for non-missing
 * leaves */
uint8_t* openMerkleTree(tree_t* tree, uint16_t* missingLeaves, size_t missingLeavesSize,
                        size_t* outputSizeBytes) {
  size_t revealedSize = 0;
  size_t* revealed = getRevealedMerkleNodes(tree, missingLeaves, missingLeavesSize, &revealedSize);

  /* Serialize output */
  *outputSizeBytes    = revealedSize * tree->dataSize;
  uint8_t* output     = malloc(*outputSizeBytes);
  uint8_t* outputBase = output;

  for (size_t i = 0; i < revealedSize; i++) {
    memcpy(output, tree->nodes[revealed[i]], tree->dataSize);
    output += tree->dataSize;
  }

  free(revealed);

  return outputBase;
}

/* addMerkleNodes: deserialize and add the data for nodes provided by the committer */
int addMerkleNodes(tree_t* tree, uint16_t* missingLeaves, size_t missingLeavesSize, uint8_t* input,
                   size_t inputSize) {
  int ret = 0;

  assert(missingLeavesSize < tree->numLeaves);

  if (inputSize > INT_MAX) {
    return -1;
  }
  int intLen = (int)inputSize;

  size_t revealedSize = 0;
  size_t* revealed = getRevealedMerkleNodes(tree, missingLeaves, missingLeavesSize, &revealedSize);
  assert(!contains(revealed, revealedSize, 0));

  /* Deserialize input */
  for (size_t i = 0; i < revealedSize; i++) {
    intLen -= tree->dataSize;
    if (intLen < 0) {
      ret = -1;
      goto Exit;
    }
    memcpy(tree->nodes[revealed[i]], input, tree->dataSize);
    input += tree->dataSize;
    tree->haveNode[revealed[i]] = 1;
  }

  if (intLen != 0) {
    ret = -1;
    goto Exit;
  }

Exit:

  free(revealed);

  return ret;
}

/* verifyMerkleTree: verify for each leaf that is set */
int verifyMerkleTree(tree_t* tree, /* uint16_t* missingLeaves, size_t missingLeavesSize, */
                     uint8_t** leafData, uint8_t* salt, faest_paramset_t* params) {
  size_t firstLeaf = tree->numNodes - tree->numLeaves;

  /* Copy the leaf data, where we have it. The actual data being committed to has already been
   * hashed, according to the spec. */
  for (size_t i = 0; i < tree->numLeaves; i++) {
    if (leafData[i] != NULL) {
      if (tree->haveNode[firstLeaf + i] == 1) {
        return -1; /* A leaf was assigned from the prover for a node we've recomputed */
      }

      if (leafData[i] != NULL) {
        memcpy(tree->nodes[firstLeaf + i], leafData[i], tree->dataSize);
        tree->haveNode[firstLeaf + i] = 1;
      }
    }
  }

  /* At this point the tree has some of the leaves, and some intermediate nodes
   * Work up the tree, computing all nodes we don't have that are missing. */
  for (int i = (int)tree->numNodes; i > 0; i--) {
    computeParentHash(tree, i, salt, params);
  }

  /* Fail if the root was not computed. */
  if (!tree->haveNode[0]) {
    return -1;
  }

  return 0;
}

/* Gets how many nodes will be there in the tree in total including root node */
uint64_t getBinaryTreeNodeCount(uint64_t depth) {
  uint64_t out = 0;
  for (uint64_t i = depth; i < 1; i--) {
    out += (1 << i);
  }
  out += 1;
  return out;
}

/* Gets the number of nodes+leaves (excluding the root node) for a tree with a given depth */
uint64_t getBinaryTreeNodeCountFromIndex(uint64_t depth) {
  uint64_t count = 0;
  for (uint32_t i = 1; i < depth; i++) {
    count += (1 << i);
  }
  return count;
}

/* Calculates the flat array index of the binary tree position */
uint64_t getNodeIndex(uint64_t depth, uint64_t pos) {

  // always between [0...(2^depth)-1]
  return (uint64_t)(((2 << depth) - 1) + pos);
}