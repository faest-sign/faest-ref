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

void NumRec(uint32_t depth, const uint8_t* bi, uint64_t* out) {
  for (uint32_t i = 0; i < depth; i++) {
    *out = *out + ((uint64_t)bi[i] * (1 << i));
  }
}

void vector_commitment(const uint8_t* rootKey, const faest_paramset_t* params, vec_com_t* vecCom) {

  uint32_t lambda         = params->faest_param.lambda;
  uint32_t lambda2        = lambda * 2;
  uint32_t lambdaBytes    = lambda / 8;
  uint32_t lambda2Bytes   = lambda2 / 8;
  uint32_t vole_instances = params->faest_param.t;

  tree_t* tree = generateSeeds(rootKey, params);

#if 1
  printTree("tree_128_t_11", tree);
  printTreeInfo("tree_128_t_11_info", tree);
#endif

  vecCom->h   = malloc(lambda2Bytes);
  vecCom->k   = malloc(tree->numNodes * lambdaBytes);
  vecCom->com = malloc(vole_instances * lambda2Bytes);
  vecCom->sd  = malloc(vole_instances * lambdaBytes);

  /* Saving the tree nodes in K */
  for (uint32_t i = 0; i < tree->numNodes; i++) {
    memcpy(vecCom->k + (i * lambdaBytes), tree->nodes[i], lambdaBytes);
  }

  uint8_t** leaves = getLeaves(tree);

  /* Doing H_0 */
  H0_context_t h0_ctx;
  /* If lambda=128, SHAKE128 else SHAKE256*/
  switch (lambda) {
  case 128:
    H0_init(&h0_ctx, 128);
    break;
  default:
    H0_init(&h0_ctx, 256);
    break;
  }
  for (uint32_t i = 0; i < vole_instances; i++) {
    H0_update(&h0_ctx, leaves[i], lambdaBytes);
    H0_final(&h0_ctx, vecCom->sd + i * lambdaBytes, lambdaBytes, vecCom->com + i * lambda2Bytes,
             lambda2Bytes);
  }

  /* Doing H_1 */
  /* If lambda=128, SHAKE128 else SHAKE256*/
  H1_context_t h1_ctx;
  switch (lambda) {
  case 128:
    H1_init(&h1_ctx, 128);
    break;
  default:
    H1_init(&h1_ctx, 256);
    break;
  }
  for (uint32_t i = 0; i < vole_instances; i++) {
    H1_update(&h1_ctx, vecCom->com + (i * lambda2Bytes), lambda2Bytes);
  }
  H1_final(&h1_ctx, vecCom->h, lambda2Bytes);
}

void vector_open(const faest_paramset_t* params, const uint8_t* k, const uint8_t* com,
                 const uint8_t* b, uint32_t leafIndex, uint8_t* pdec, uint8_t* com_j) {

  uint32_t lambda         = params->faest_param.lambda;
  uint32_t lambda2        = lambda * 2;
  uint32_t lambdaBytes    = lambda / 8;
  uint32_t lambda2Bytes   = lambda2 / 8;
  uint32_t vole_instances = params->faest_param.t;

  uint32_t depth = ceil_log2(params->faest_param.t);

  memcpy(com_j, com + (leafIndex * lambda2Bytes), lambda2Bytes);

  uint32_t a = 0;
  for (uint32_t d = 0; d < depth; d++) {
    memcpy(pdec + (lambdaBytes * d),
           k + (lambdaBytes * getNodeIndex(d, (2 * a) + b[(depth - 1) - d])), lambdaBytes);
    a = (2 * a) + b[(depth - 1) - d];
  }
}

void vector_reconstruction(const faest_paramset_t* params, const uint8_t* pdec, const uint8_t* b,
                           vec_com_rec_t* VecComRec) {

  uint32_t lambda         = params->faest_param.lambda;
  uint32_t lambda2        = lambda * 2;
  uint32_t lambdaBytes    = lambda / 8;
  uint32_t lambda2Bytes   = lambda2 / 8;
  uint32_t vole_instances = params->faest_param.t;

  uint8_t depth = sizeof(b);
  uint64_t j_;
  NumRec(depth, b, &j_);

  VecComRec->h   = malloc(lambda2Bytes);
  VecComRec->k   = malloc(getBinaryTreeNodeCountFromIndex(depth) * lambdaBytes);
  VecComRec->com = malloc((vole_instances * lambda2Bytes) - 1);
  VecComRec->m   = malloc((vole_instances * lambdaBytes) - 1);

  memcpy(VecComRec->com, pdec[0], lambda2);

  uint32_t a = 0;
  for (uint32_t d = 1; d <= depth; d++) {
    memcpy(VecComRec->k + getNodeIndex(d, 2 * a + b[d]), pdec + (lambda2Bytes * d), lambda2Bytes);
    a = a * 2 + b[d];
  }

  /* Doing H_0 */
  H0_context_t h0_ctx;
  /* If lambda=128, SHAKE128 else SHAKE256*/
  switch (lambda) {
  case 128:
    H0_init(&h0_ctx, 128);
    break;
  default:
    H0_init(&h0_ctx, 256);
    break;
  }
  /* Doing H_1 */
  /* If lambda=128, SHAKE128 else SHAKE256*/
  H1_context_t h1_ctx;
  switch (lambda) {
  case 128:
    H1_init(&h1_ctx, 128);
    break;
  default:
    H1_init(&h1_ctx, 256);
    break;
  }

  for (uint32_t j = 0; j < (1 << depth); j++) {
    if (j == j_) {
      continue;
    }
    for (uint32_t i = 0; i < vole_instances; i++) {
      H0_update(&h0_ctx, VecComRec->k + getNodeIndex(depth, i), lambda2Bytes);
    }
    H0_final(&h0_ctx, VecComRec->m, lambda, VecComRec->com, lambda2Bytes);

    for (uint32_t i = 0; i < vole_instances; i++) {
      H1_update(&h1_ctx, VecComRec->com + (i * lambda2Bytes), lambda2Bytes);
    }
    H1_final(&h1_ctx, VecComRec->h, lambda2Bytes);
  }
}

int vector_verify(const faest_paramset_t* params, const uint8_t* pdec, const uint8_t* b,
                  vec_com_rec_t* VecComRec, vec_com_t* VecCom) {
  vector_reconstruction(params, pdec, b, VecComRec);
  if (memcmp(VecCom->com, VecComRec->com, sizeof(VecComRec->com)) == 0) {
    return 1;
  } else {
    return 0;
  }
}