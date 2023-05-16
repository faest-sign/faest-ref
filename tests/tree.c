#include <assert.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>

#include "../hash_shake.h"
#include "../tree.h"
#include "../instances.h"
#include "../utils.h"
#include "../compat.h"

int runSeedTest(uint8_t* rootKey, faest_paramset_t* params, uint32_t numVoleInstances) {

  tree_t* tree = generateSeeds(rootKey, params, numVoleInstances);

#if 0
  printTreeInfo("tree", tree);
  printTree("tree", tree);
  printLeaves(tree);
#endif

  // TODO: make tests
  return 1;
}

int main(void) {

  size_t tests  = 6;
  size_t passed = 0;

  printf("Running seed tree tests\n");

  /* We can use the uint8 32 root, the 128 and 192 aes turncate the root key appropriately in the
  implementation */
  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  for (faest_paramid_t p = FAEST_128S; p <= FAEST_256F; p++) {
    faest_paramset_t params = faest_get_paramset(p);

    passed += runSeedTest(rootKey, &params, params.faest_param.k0);
  }
  printf("Done, %lu of %lu tests passed\n", passed, tests);
  return 0;
}
