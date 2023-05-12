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
                       tree_t* tree) {

  uint32_t lambda         = params->faest_param.lambda;
  uint32_t lambda2        = lambda * 2;
  uint32_t lambdaBytes    = lambda / 8;
  uint32_t lambda2Bytes   = lambda2 / 8;
  uint32_t vole_instances = params->faest_param.t;

  *tree = *(generateSeeds(rootKey, params));

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
                 const uint8_t* b, uint8_t* pdec, uint8_t* com_j) {

  uint32_t lambda         = params->faest_param.lambda;
  uint32_t lambda2        = lambda * 2;
  uint32_t lambdaBytes    = lambda / 8;
  uint32_t lambda2Bytes   = lambda2 / 8;
  uint32_t vole_instances = params->faest_param.t;

  uint32_t depth = ceil_log2(params->faest_param.t);

  uint64_t leafIndex = NumRec(depth, b);

  memcpy(com_j, com + (leafIndex * lambda2Bytes), lambda2Bytes);

  uint32_t a = 0;
  for (uint32_t d = 0; d < depth; d++) {
    memcpy(pdec + (lambdaBytes * d),
           k + (lambdaBytes * getNodeIndex(d, (2 * a) + b[(depth - 1) - d])), lambdaBytes);
    a = (2 * a) + b[(depth - 1) - d];
  }
}

void vector_reconstruction(const faest_paramset_t* params, const uint8_t* pdec,
                           const uint8_t* com_j, const uint8_t* b, const vec_com_t* VecCom,
                           vec_com_rec_t* VecComRec) {

  uint32_t lambda         = params->faest_param.lambda;
  uint32_t lambda2        = lambda * 2;
  uint32_t lambdaBytes    = lambda / 8;
  uint32_t lambda2Bytes   = lambda2 / 8;
  uint32_t vole_instances = params->faest_param.t;

  uint8_t depth      = ceil_log2(params->faest_param.t);
  uint64_t leafIndex = NumRec(depth, b);

  VecComRec->h = malloc(lambda2Bytes);
  VecComRec->k =
      malloc((getBinaryTreeNodeCount(params) - 1) * lambdaBytes); // excluding the root k_0,0
  VecComRec->com = malloc(vole_instances * lambda2Bytes);
  VecComRec->m   = malloc(vole_instances * lambdaBytes);

  memcpy(VecComRec->com + (lambda2Bytes * leafIndex), com_j, lambda2Bytes);

  memcpy(VecComRec->k, VecCom->k + lambdaBytes, (getBinaryTreeNodeCount(params) - 1) * lambdaBytes);

  uint32_t a = 0;
  for (uint32_t d = 0; d < depth; d++) {
    memcpy(VecComRec->k + (lambdaBytes * (getNodeIndex(d, 2 * a + b[(depth - 1) - d]) - 1)),
           pdec + (lambdaBytes * d), lambdaBytes);
    a = a * 2 + b[(depth - 1) - d];
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

  for (uint32_t j = 0; j < vole_instances; j++) {
    if (j == leafIndex) {
      continue;
    }
    /* -1 because root node is not present here */
    H0_update(&h0_ctx, VecComRec->k + ((getNodeIndex(depth, j) - 1) * lambdaBytes), lambdaBytes);
    H0_final(&h0_ctx, VecComRec->m + (lambdaBytes * j), lambdaBytes,
             VecComRec->com + (lambda2Bytes * j), lambda2Bytes);
  }

  H1_update(&h1_ctx, VecComRec->com, lambda2Bytes * vole_instances);
  H1_final(&h1_ctx, VecComRec->h, lambda2Bytes);
}

int vector_verify(const faest_paramset_t* params, const uint8_t* pdec, const uint8_t* com_j,
                  const uint8_t* b, const vec_com_t* VecCom, vec_com_rec_t* VecComRec) {
  vector_reconstruction(params, pdec, com_j, b, VecCom, VecComRec);
  if (memcmp(VecCom->com, VecComRec->com,
             params->faest_param.t * (params->faest_param.seclvl * 2)) == 0) {
    return 1;
  } else {
    return 0;
  }
}