#include "instances.h"

typedef struct tree_t {
  size_t depth;      /* The depth of the tree */
  uint8_t** nodes;   /* The data for each node */
  size_t dataSize;   /* The size data at each node, in bytes */
  uint8_t* haveNode; /* If we have the data (seed or hash) for node i, haveSeed[i] is 1 */
  uint8_t* exists;   /* Since the tree is not always complete, nodes marked 0 don't exist */
  size_t numNodes;   /* The total number of nodes in the tree */
  size_t numLeaves;  /* The total number of leaves in the tree */
} tree_t;

/* The largest seed size is 256 bits, for the Picnic3-L5-FS parameter set. */
#define MAX_SEED_SIZE_BYTES (32)

tree_t* createTree(const faest_paramset_t* params, uint32_t numVoleInstances);
void freeTree(tree_t* tree);
uint8_t** getLeaves(tree_t* tree);
/* Get one leaf, leafIndex must be in [0, tree->numLeaves -1] */
uint8_t* getLeaf(tree_t* tree, size_t leafIndex);
void printLeaves(tree_t* tree);

/* Functions for trees used to derive seeds.
 *    Signer's usage:   generateSeeds -> revealSeeds -> freeTree
 *    Verifier's usage: createTree -> reconstructSeeds -> freeTree
 */

/* Returns the number of bytes written to output.  A safe number of bytes for
 * callers to allocate is numLeaves*params->seedSizeBytes, or call revealSeedsSize. */
tree_t* generateSeeds(const uint8_t* rootSeed, const faest_paramset_t* params,
                      uint32_t numVoleInstances);
size_t revealSeeds(tree_t* tree, uint16_t* hideList, size_t hideListSize, uint8_t* output,
                   size_t outputLen, faest_paramset_t* params);
size_t revealSeedsSize(uint16_t* hideList, size_t hideListSize, faest_paramset_t* params);
int reconstructSeeds(tree_t* tree, uint16_t* hideList, size_t hideListSize, uint8_t* input,
                     size_t inputLen, faest_paramset_t* params);

/* Functions for Merkle hash trees used for commitments.
 *
 * Signer call sequence:
 *     1. createTree
 *     2. buildMerkleTree  with all commitments as leaf nodes
 *     3. openMerkleTree   with missingLeaves - list of commitments the verifier won't recompute
 *     4. freeTree
 *  Verifier call sequence
 *      1. createTree
 *      2. addMerkleNodes       with the output of the signer
 *      3. verifyMerkleTree     Checks that all leaf nodes present are correct commitments
 *      4. freeTree
 */
void buildMerkleTree(tree_t* tree, uint8_t** leafData, uint8_t* salt, faest_paramset_t* params);
uint8_t* openMerkleTree(tree_t* tree, uint16_t* missingLeaves, size_t missingLeavesSize,
                        size_t* outputSizeBytes);
size_t openMerkleTreeSize(uint16_t* missingLeaves, size_t missingLeavesSize,
                          faest_paramset_t* params);
int addMerkleNodes(tree_t* tree, uint16_t* missingLeaves, size_t missingLeavesSize, uint8_t* input,
                   size_t inputSize);
int verifyMerkleTree(tree_t* tree, uint8_t** leafData, uint8_t* salt, faest_paramset_t* params);
uint64_t getBinaryTreeNodeCount(uint32_t numVoleInstances);
uint64_t getNodeIndex(uint64_t depth, uint64_t levelIndex);
void printTree(const char* label, tree_t* tree);
void printTreeInfo(const char* label, tree_t* tree);
int contains(size_t* list, size_t len, size_t value);