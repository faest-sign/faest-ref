#include "vc.h"

/* Gets the bit string of a node according to its position in the binary tree */
/* idx -> 2 -> {0,1},, Little Endian */
int BitDec(uint32_t leafIndex, uint32_t depth, uint8_t* out) {
  uint32_t i = leafIndex;
  if (leafIndex >= (uint32_t)(2 << depth)) {
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

void vector_commitment(const uint8_t* rootKey, const faest_paramset_t* params, vec_com_t* vecCom,
                       tree_t* tree, uint32_t voleInstances) {

  uint32_t lambda  = params->faest_param.lambda / 8;
  uint32_t lambda2 = params->faest_param.lambda / 4;

  /* Generating the tree with a rootkey */
  *tree = *(generateSeeds(rootKey, params, voleInstances));

  vecCom->h             = malloc(lambda2);
  vecCom->k_uint_size   = lambda;
  vecCom->k             = malloc(tree->numNodes * vecCom->k_uint_size);
  vecCom->com_unit_size = lambda2;
  vecCom->com           = malloc(voleInstances * vecCom->com_unit_size);
  vecCom->sd_uint_size  = lambda;
  vecCom->sd            = malloc(voleInstances * vecCom->sd_uint_size);

  /* Saving the tree nodes in K */
  for (uint32_t i = 0; i < tree->numNodes; i++) {
    memcpy(vecCom->k + (i * vecCom->k_uint_size), tree->nodes[i], vecCom->k_uint_size);
  }

  uint8_t** leaves = getLeaves(tree);

  for (uint32_t i = 0; i < voleInstances; i++) {
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
    H0_update(&h0_ctx, leaves[i], vecCom->k_uint_size);
    H0_final(&h0_ctx, vecCom->sd + i * vecCom->sd_uint_size, vecCom->sd_uint_size,
             vecCom->com + i * vecCom->com_unit_size, vecCom->com_unit_size);
  }

  H1_context_t h1_ctx;
  switch (lambda) {
  case 128:
    H1_init(&h1_ctx, 128);
    break;
  default:
    H1_init(&h1_ctx, 256);
    break;
  }
  for (uint32_t i = 0; i < voleInstances; i++) {
    H1_update(&h1_ctx, vecCom->com + (i * vecCom->com_unit_size), vecCom->com_unit_size);
  }
  /* Generating final commitment from all the com commitments */
  H1_final(&h1_ctx, vecCom->h, lambda2);
}

void vector_open(const uint8_t* k, const uint8_t* com, const uint8_t* b, uint8_t* pdec,
                 uint8_t* com_j, uint32_t numVoleInstances, uint32_t lambdabits,
                 vec_com_rec_t* vecComRec, const vec_com_t* vecCom) {

  uint32_t lambda    = lambdabits / 8;
  uint32_t lambda2   = lambdabits / 4;
  uint32_t depth     = ceil_log2(numVoleInstances);
  uint64_t leafIndex = NumRec(depth, b);

  memcpy(com_j, com + (leafIndex * lambda2), lambda2);
  uint32_t a = 0;
  for (uint32_t d = 0; d < depth; d++) {
    memcpy(pdec + (lambda * d), k + (lambda * getNodeIndex(d + 1, (2 * a) + b[(depth - 1) - d])),
           lambda);
    a = (2 * a) + b[(depth - 1) - d];
  }

  // TODO
  /* Unsure but perhaps we open the nodes (k) of the tree exlcuding the pdec nodes here and not in
   * vec reconstruct */
  vecComRec->k = malloc(getBinaryTreeNodeCount(numVoleInstances) * lambda);
  memset(vecComRec->k, 0, lambda); // setting root node to 0
  memcpy(vecComRec->k + lambda, vecCom->k + lambda,
         (getBinaryTreeNodeCount(numVoleInstances) - 1) * lambda);
  a              = 0;
  uint8_t* zeros = malloc(lambda);
  memset(zeros, 0, lambda);
  for (uint32_t d = 0; d < depth; d++) {
    memcpy(vecComRec->k + (lambda * (getNodeIndex((d + 1), 2 * a + b[(depth - 1) - d]))), zeros,
           lambda);
    a = a * 2 + b[(depth - 1) - d];
  }
}

void vector_reconstruction(const faest_paramset_t* params, const uint8_t* pdec,
                           const uint8_t* com_j, const uint8_t* b, uint32_t lambdaBits,
                           uint32_t NumVoleInstances, const vec_com_t* VecCom,
                           vec_com_rec_t* VecComRec) {

  uint32_t lambda    = lambdaBits / 8;
  uint32_t lambda2   = lambdaBits / 4;
  uint8_t depth      = ceil_log2(NumVoleInstances);
  uint64_t leafIndex = NumRec(depth, b);

  VecComRec->h   = malloc(lambda2);
  VecComRec->com = malloc(NumVoleInstances * lambda2);
  VecComRec->m   = malloc(NumVoleInstances * lambda);

  memcpy(VecComRec->com + (lambda2 * leafIndex), com_j, lambda2);
  for (uint32_t j = 0; j < NumVoleInstances; j++) {
    /* Reconstruct the coms and the m from the ks while keeping k_j* secret */
    if (j != leafIndex) {
      H0_context_t h0_ctx;
      switch (lambdaBits) {
      case 128:
        H0_init(&h0_ctx, 128);
        break;
      default:
        H0_init(&h0_ctx, 256);
        break;
      }
      H0_update(&h0_ctx, VecComRec->k + (getNodeIndex(depth, j) * lambda), lambda);
      H0_final(&h0_ctx, VecComRec->m + (lambda * j), lambda, VecComRec->com + (lambda2 * j),
               lambda2);
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
  H1_update(&h1_ctx, VecComRec->com, lambda2 * NumVoleInstances);
  H1_final(&h1_ctx, VecComRec->h, lambda2);
}

int vector_verify(const faest_paramset_t* params, const uint8_t* pdec, const uint8_t* com_j,
                  const uint8_t* b, uint32_t lambdaBits, uint32_t numVoleInstances,
                  const vec_com_t* VecCom, vec_com_rec_t* VecComRec) {
  vector_reconstruction(params, pdec, com_j, b, lambdaBits, numVoleInstances, VecCom, VecComRec);
  if (memcmp(VecCom->com, VecComRec->com,
             params->faest_param.k0 * (params->faest_param.seclvl / 4)) == 0) {
    return 1;
  } else {
    return 0;
  }
}