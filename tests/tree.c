#include <assert.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>

#include "../hash_shake.h"
#include "../tree.h"
#include "../instances.h"
#include "../utils.h"
#include "../compat.h"

int runSeedTest(uint16_t* hideList, size_t hideListSize, uint8_t* rootKey,
                faest_paramset_t* params) {
  int freeHideList   = 0;
  int ret            = 1;
  size_t numVoleInst = params->faest_param.t; // number of leaves

  if (numVoleInst < hideListSize - 1) {
    printf("%s invalid input (numLeaves = %lu, hideListSize = %lu)\n", __func__, numVoleInst,
           hideListSize);
    return 0;
  }

  if (hideList == NULL) {
    hideList     = malloc(hideListSize * sizeof(uint16_t));
    freeHideList = 1;
    uint16_t val;
    for (size_t i = 0; i < hideListSize; i++) {
      do {
        val = ((uint16_t)rand()) % numVoleInst;
      } while (contains(hideList, i, val));
      hideList[i] = val;
    }
  }

  // printf("%s: Generating seeds\n", __func__);
  tree_t* tree  = generateSeeds(rootKey, params);
  tree_t* tree2 = createTree(params);

#if 0
  printTree("tree", tree);
#endif

  size_t initialOutputSize = (tree->numLeaves) * params->faest_param.seedSizeBytes;
  uint8_t* output          = malloc(initialOutputSize);

  size_t expectedOutputLen = revealSeedsSize(hideList, hideListSize, params);
  if (hideListSize > 0 && expectedOutputLen == 0) {
    printf("Failed to get exepctedOutputLen\n");
    ret = 0;
    goto Exit;
  }
  if (expectedOutputLen % params->faest_param.seedSizeBytes != 0) {
    printf("ExepctedOutputLen is not a multiple of the seed length\n");
    ret = 0;
    goto Exit;
  }

  // printf("%s: Revealing seeds\n", __func__);
  size_t outputLen = revealSeeds(tree, hideList, hideListSize, output, initialOutputSize, params);
  if (outputLen == 0) {
    printf("Failed to revealSeeds, output buffer too small\n");
    ret = 0;
    goto Exit;
  }

  if (outputLen != expectedOutputLen) {
    printf("Expected output lengthd doesn't match output length\n");
    ret = 0;
    goto Exit;
  }

  // printf("%s: numLeaves = %lu, revealed %lu\n", __func__, tree->numLeaves,
  // outputLen/tree->dataSize);

  if (params->faest_param.numOpenRounds *
          ceil_log2(params->faest_param.t / params->faest_param.numOpenRounds) <
      outputLen / tree->dataSize) {
    printf("%s: Output length is larger than expected\n", __func__);
    ret = 0;
    goto Exit;
  }

  // printf("%s: Reconstructing seeds\n", __func__);
  int res = reconstructSeeds(tree2, hideList, hideListSize, output, outputLen, params);
  if (res != 0) {
    printf("%s: Reconstructing seeds FAILED\n", __func__);
    ret = 0;
    goto Exit;
  }

  // printf("seeds in reconstructed tree:\n");
  // printSeeds(tree2->nodes[0], params->seedSizeBytes, 15 );

  // Check that we have the correct seeds, and that they match
  size_t firstLeaf = tree->numNodes - tree->numLeaves;
  for (size_t i = firstLeaf; i < tree->numNodes; i++) {
    if (contains(hideList, hideListSize, i - firstLeaf)) {
      if (tree2->haveNode[i]) {
        printf("%s FAIL: reconstructed tree contains a seed that should have been hidden, node %lu "
               "(leaf node %lu)\n",
               __func__, i, i - firstLeaf);
        printHex("tree->nodes[i] ", tree->nodes[i], params->faest_param.seedSizeBytes);
        printHex("tree2->nodes[i]", tree2->nodes[i], params->faest_param.seedSizeBytes);
        ret = 0;
        goto Exit;
      }
    } else {

      if (!tree2->haveNode[i]) {
        printf("%s FAIL: expected to have seed for node %lu, but don't\n", __func__, i);
        ret = 0;
        goto Exit;
      }
      if (!tree->haveNode[i]) {
        printf("%s FAIL: initial tree is missing node %lu -- not contructed properly?\n", __func__,
               i);
        // printTreeInfo("tree", tree);
        ret = 0;
        goto Exit;
      }

      if (memcmp(tree->nodes[i], tree2->nodes[i], params->faest_param.seedSizeBytes) != 0) {
        printf("%s FAIL: reconstructed tree has an incorrect seed node %lu\n", __func__, i);
        ret = 0;
        goto Exit;
      }
    }
  }

Exit:
  if (freeHideList) {
    free(hideList);
  }
  free(output);
  freeTree(tree);
  freeTree(tree2);

  return ret;
}

int main(void) {

  size_t tests  = 0;
  size_t passed = 0;

  printf("Running seed tree tests\n");

  /* We can use the uint8 32 root, the 128 and 192 aes turncate the root key appropriately in the
  implementation */
  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  for (faest_paramid_t p = FAEST_128S; p <= FAEST_256F; p++) {
    faest_paramset_t params = faest_get_paramset(p);

    passed += runSeedTest(NULL, params.faest_param.numOpenRounds, rootKey, &params);
    tests++;

    printf("Done, %lu of %lu tests passed\n", passed, tests);
  }
  return 0;
}
