#include "vc.h"

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
    out = out + ((uint64_t)bi[i] * (1 << i));
  }
  return out;
}

void vector_commitment(const uint8_t* rootKey, const faest_paramset_t* params, uint32_t lambda,
                       uint32_t lambdaBytes, vec_com_t* vecCom, uint32_t numVoleInstances) {

  // Generating the tree
  tree_t* tree = generateSeeds(rootKey, params, numVoleInstances);

  // Initialzing stuff
  vecCom->h   = malloc(lambdaBytes * 2);
  vecCom->k   = malloc(tree->numNodes * lambdaBytes);
  vecCom->com = malloc(numVoleInstances * lambdaBytes * 2);
  vecCom->sd  = malloc(numVoleInstances * lambdaBytes);

  // Step: 1..3
  for (uint32_t i = 0; i < tree->numNodes; i++) {
    memcpy(vecCom->k + (i * lambdaBytes), tree->nodes[i], lambdaBytes);
  }

  // Step: 4..5
  uint8_t** leaves = getLeaves(tree);
  for (uint32_t i = 0; i < numVoleInstances; i++) {
    H0_context_t h0_ctx;
    H0_init(&h0_ctx, lambda);
    H0_update(&h0_ctx, leaves[i], lambdaBytes);
    H0_final(&h0_ctx, vecCom->sd + (i * lambdaBytes), lambdaBytes,
             vecCom->com + (i * (lambdaBytes * 2)), (lambdaBytes * 2));
  }

  // Step: 6
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  for (uint32_t i = 0; i < numVoleInstances; i++) {
    H1_update(&h1_ctx, vecCom->com + (i * (lambdaBytes * 2)), (lambdaBytes * 2));
  }
  H1_final(&h1_ctx, vecCom->h, lambdaBytes * 2);
}

void vector_open(const uint8_t* k, const uint8_t* com, const uint8_t* b, uint8_t* pdec,
                 uint8_t* com_j, uint32_t numVoleInstances, uint32_t lambdaBytes) {

  uint32_t depth = ceil_log2(numVoleInstances);

  // Step: 1
  uint64_t leafIndex = NumRec(depth, b);

  // Step: 3..6
  uint32_t a = 0;
  for (uint32_t d = 0; d < depth; d++) {
    memcpy(pdec + (lambdaBytes * d),
           k + (lambdaBytes * getNodeIndex(d + 1, (2 * a) + !b[(depth - 1) - d])), lambdaBytes);
    a = (2 * a) + b[(depth - 1) - d];
  }

  // Step: 7
  memcpy(com_j, com + (leafIndex * lambdaBytes * 2), lambdaBytes * 2);
}

void vector_reconstruction(const uint8_t* pdec, const uint8_t* com_j, const uint8_t* b,
                           uint32_t lambda, uint32_t lambdaBytes, uint32_t numVoleInstances,
                           vec_com_rec_t* vecComRec) {

  // Initializing
  uint8_t depth      = ceil_log2(numVoleInstances);
  uint64_t leafIndex = NumRec(depth, b);
  vecComRec->h       = malloc(lambdaBytes * 2);
  vecComRec->com     = malloc(numVoleInstances * lambdaBytes * 2);
  vecComRec->m       = malloc(numVoleInstances * lambdaBytes);
  vecComRec->k       = malloc(getBinaryTreeNodeCount(numVoleInstances) * lambdaBytes);

  // Step: 3..9
  uint8_t* zeros = malloc(lambdaBytes);
  memset(zeros, 0, lambdaBytes);
  memset(vecComRec->k, 0, lambdaBytes);
  uint32_t a = 0;
  for (uint32_t d = 1; d <= depth; d++) {
    memcpy(vecComRec->k + (lambdaBytes * getNodeIndex(d, 2 * a + !b[depth - d])),
           pdec + (lambdaBytes * (d - 1)), lambdaBytes);
    memset(vecComRec->k + (lambdaBytes * getNodeIndex(d, 2 * a + b[depth - d])), 0, lambdaBytes);
    a = a * 2 + b[depth - d];

    for (uint32_t j = 0; j < (1 << d); j++) {
      if (memcmp(vecComRec->k + (lambdaBytes * getNodeIndex(d, j)), zeros, lambdaBytes) == 0) {
        continue;
      }
      uint8_t iv[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
      uint8_t* out   = malloc(lambdaBytes * 2);
      prg(vecComRec->k + (lambdaBytes * getNodeIndex(d, j)), iv, out, lambda, lambdaBytes * 2);
      memcpy(vecComRec->k + (lambdaBytes * getNodeIndex(d + 1, 2 * j)), out, lambdaBytes);
      memcpy(vecComRec->k + (lambdaBytes * getNodeIndex(d + 1, (2 * j) + 1)), out + lambdaBytes,
             lambdaBytes);
    }
  }

  // Step: 10..11
  for (uint32_t j = 0; j < numVoleInstances; j++) {
    /* Reconstruct the coms and the m from the ks while keeping k_j* secret */
    if (j != leafIndex) {
      H0_context_t h0_ctx;
      H0_init(&h0_ctx, lambda);
      H0_update(&h0_ctx, vecComRec->k + (getNodeIndex(depth, j) * lambdaBytes), lambdaBytes);
      H0_final(&h0_ctx, vecComRec->m + (lambdaBytes * j), lambdaBytes,
               vecComRec->com + (lambdaBytes * 2 * j), lambdaBytes * 2);
    }
  }
  // Step: 12
  memcpy(vecComRec->com + (lambdaBytes * 2 * leafIndex), com_j, lambdaBytes * 2);
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  H1_update(&h1_ctx, vecComRec->com, lambdaBytes * 2 * numVoleInstances);
  H1_final(&h1_ctx, vecComRec->h, lambdaBytes * 2);
}

int vector_verify(const uint8_t* pdec, const uint8_t* com_j, const uint8_t* b, uint32_t lambda,
                  uint32_t lambdaBytes, uint32_t numVoleInstances, vec_com_rec_t* vecComRec,
                  uint8_t* vecComH) {

  vector_reconstruction(pdec, com_j, b, lambda, lambdaBytes, numVoleInstances, vecComRec);
  if (memcmp(vecComH, vecComRec->h, lambdaBytes * 2) == 0) {
    return 1;
  } else {
    return 0;
  }
}