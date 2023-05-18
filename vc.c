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
                       uint32_t lambdaBytes, vec_com_t* vecCom, tree_t* tree,
                       uint32_t numVoleInstances) {

  /* Generating the tree with a rootkey */
  *tree = *generateSeeds(rootKey, params, numVoleInstances);

  vecCom->h   = malloc(lambdaBytes * 2);
  vecCom->k   = malloc(tree->numNodes * lambdaBytes);
  vecCom->com = malloc(numVoleInstances * lambdaBytes * 2);
  vecCom->sd  = malloc(numVoleInstances * lambdaBytes);

  /* Saving the tree nodes in K */
  for (uint32_t i = 0; i < tree->numNodes; i++) {
    memcpy(vecCom->k + (i * lambdaBytes), tree->nodes[i], lambdaBytes);
  }

  uint8_t** leaves = getLeaves(tree);

  for (uint32_t i = 0; i < numVoleInstances; i++) {
    H0_context_t h0_ctx;
    switch (params->faest_param.lambda) {
    case 128:
      H0_init(&h0_ctx, 128);
      break;
    default:
      H0_init(&h0_ctx, 256);
      break;
    }
    /* Generating the sd messages and com commitments from each leaf */
    H0_update(&h0_ctx, leaves[i], lambdaBytes);
    H0_final(&h0_ctx, vecCom->sd + (i * lambdaBytes), lambdaBytes,
             vecCom->com + (i * (lambdaBytes * 2)), (lambdaBytes * 2));
  }

  H1_context_t h1_ctx;
  switch (params->faest_param.lambda) {
  case 128:
    H1_init(&h1_ctx, 128);
    break;
  default:
    H1_init(&h1_ctx, 256);
    break;
  }
  for (uint32_t i = 0; i < numVoleInstances; i++) {
    H1_update(&h1_ctx, vecCom->com + (i * (lambdaBytes * 2)), (lambdaBytes * 2));
  }
  /* Generating final commitment from all the com commitments */
  H1_final(&h1_ctx, vecCom->h, params->faest_param.lambda / 4);
}

void vector_open(const uint8_t* k, const uint8_t* com, const uint8_t* b, uint8_t* pdec,
                 uint8_t* com_j, uint32_t numVoleInstances, uint32_t lambda, uint32_t lambdaBytes,
                 vec_com_rec_t* vecComRec, const vec_com_t* vecCom) {

  uint32_t depth     = ceil_log2(numVoleInstances);
  uint64_t leafIndex = NumRec(depth, b);

  memcpy(com_j, com + (leafIndex * lambdaBytes * 2), lambdaBytes * 2);
  uint32_t a = 0;
  for (uint32_t d = 0; d < depth; d++) {
    memcpy(pdec + (lambdaBytes * d),
           k + (lambdaBytes * getNodeIndex(d + 1, (2 * a) + b[(depth - 1) - d])), lambdaBytes);
    a = (2 * a) + b[(depth - 1) - d];
  }

  // TODO
  /* Unsure but perhaps we open the nodes (k) of the tree exlcuding the pdec nodes here and not in
   * vec reconstruct */
  vecComRec->k = malloc(getBinaryTreeNodeCount(numVoleInstances) * lambdaBytes);
  memset(vecComRec->k, 0, lambdaBytes); // setting root node to 0
  memcpy(vecComRec->k + lambdaBytes, vecCom->k + lambdaBytes,
         (getBinaryTreeNodeCount(numVoleInstances) - 1) * lambdaBytes);
  a              = 0;
  uint8_t* zeros = malloc(lambdaBytes);
  memset(zeros, 0, lambdaBytes);
  for (uint32_t d = 0; d < depth; d++) {
    memcpy(vecComRec->k + (lambdaBytes * (getNodeIndex((d + 1), 2 * a + b[(depth - 1) - d]))),
           zeros, lambdaBytes);
    a = a * 2 + b[(depth - 1) - d];
  }
}

// TODO: use pdec here !!!
void vector_reconstruction(const uint8_t* pdec, const uint8_t* com_j, const uint8_t* b,
                           uint32_t lambda, uint32_t lambdaBytes, uint32_t NumVoleInstances,
                           vec_com_rec_t* VecComRec) {

  uint8_t depth      = ceil_log2(NumVoleInstances);
  uint64_t leafIndex = NumRec(depth, b);

  VecComRec->h   = malloc(lambdaBytes * 2);
  VecComRec->com = malloc(NumVoleInstances * lambdaBytes * 2);
  VecComRec->m   = malloc(NumVoleInstances * lambdaBytes);

  // TODO: actually reconstruct the k^d_0...2^d

  memcpy(VecComRec->com + (lambdaBytes * 2 * leafIndex), com_j, lambdaBytes * 2);
  for (uint32_t j = 0; j < NumVoleInstances; j++) {
    /* Reconstruct the coms and the m from the ks while keeping k_j* secret */
    if (j != leafIndex) {
      H0_context_t h0_ctx;
      switch (lambda) {
      case 128:
        H0_init(&h0_ctx, 128);
        break;
      default:
        H0_init(&h0_ctx, 256);
        break;
      }
      H0_update(&h0_ctx, VecComRec->k + (getNodeIndex(depth, j) * lambdaBytes), lambdaBytes);
      H0_final(&h0_ctx, VecComRec->m + (lambdaBytes * j), lambdaBytes,
               VecComRec->com + (lambdaBytes * 2 * j), lambdaBytes * 2);
    }
  }
  /* Compute the final com with all the computed coms and com_j */
  H1_context_t h1_ctx;
  switch (lambda) {
  case 128:
    H1_init(&h1_ctx, 128);
    break;
  default:
    H1_init(&h1_ctx, 256);
    break;
  }
  H1_update(&h1_ctx, VecComRec->com, lambdaBytes * 2 * NumVoleInstances);
  H1_final(&h1_ctx, VecComRec->h, lambdaBytes * 2);
}

int vector_verify(const uint8_t* pdec, const uint8_t* com_j, const uint8_t* b, uint32_t lambda,
                  uint32_t lambdaBytes, uint32_t numVoleInstances, const vec_com_t* VecCom,
                  vec_com_rec_t* VecComRec) {
  vector_reconstruction(pdec, com_j, b, lambda, lambdaBytes, numVoleInstances, VecComRec);
  if (memcmp(VecCom->h, VecComRec->h, lambdaBytes * 2) == 0) {
    return 1;
  } else {
    return 0;
  }
}