#include "../vole.h"

int test_ConvertToVoleProver() {

  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  faest_paramset_t params = faest_get_paramset(1); // Just using the FAEST-128s
  vec_com_t vecCom;
  // vec_com_rec_t vecComRec;
  uint32_t numVoleInstances = (1 << params.faest_param.k0);
  uint32_t lambda           = params.faest_param.lambda;
  uint32_t lambdaBytes      = params.faest_param.lambdaBytes;
  vector_commitment(rootKey, &params, lambda, lambdaBytes, &vecCom, numVoleInstances);

  uint32_t leafIndex = 7;
  uint32_t depth     = ceil_log2(numVoleInstances);
  uint8_t* b         = malloc(depth);
  BitDec(leafIndex, depth, b);

  uint8_t* pdec  = malloc(depth * lambdaBytes);
  uint8_t* com_j = malloc(lambdaBytes * 2);
  vector_open(vecCom.k, vecCom.com, b, pdec, com_j, numVoleInstances, lambdaBytes);

  uint32_t outlen = 16;
  uint8_t* u      = malloc(outlen);
  uint8_t* v      = malloc(outlen * depth);
  ConvertToVoleProver(lambda, lambdaBytes, vecCom.sd, numVoleInstances, depth, outlen, u, v);

// TODO: write better test cases : )
#if 0
  printf("r\n");
  for (uint32_t i = 0; i < getBinaryTreeNodeCount(numVoleInstances); i++) {
    printf("%d ", i);
    for (uint32_t j = 0; j < outlen; j++) {
      printf("%.2x", *(r + j + (i * outlen)));
    }
    printf("\n");
  }
  printf("\n");
  printf("u\n");
  for (uint32_t i = 0; i < outlen; i++) {
    printf("%.2x", *(u + i));
  }
  printf("\nv_0...d-1\n");
  for (uint32_t d = 0; d < depth; d++) {
    for (uint32_t i = 0; i < outlen; i++) {
      printf("%.2x", *(v + (d * outlen) + i));
    }
    printf("\n");
  }
#endif

  // TODO: make tests !!
  return 1;
}

int test_ConvertToVoleVerifier() {

  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  faest_paramset_t params = faest_get_paramset(1); // Just using the FAEST-128s
  vec_com_t vecCom;
  vec_com_rec_t vecComRec;
  uint32_t numVoleInstances = (1 << params.faest_param.k0);
  uint32_t lambda           = params.faest_param.lambda;
  uint32_t lambdaBytes      = params.faest_param.lambdaBytes;
  vector_commitment(rootKey, &params, lambda, lambdaBytes, &vecCom, numVoleInstances);

  uint32_t leafIndex = 7;
  uint32_t depth     = ceil_log2(numVoleInstances);
  uint8_t* b         = malloc(depth);
  BitDec(leafIndex, depth, b);

  uint8_t* pdec  = malloc(depth * lambdaBytes);
  uint8_t* com_j = malloc(lambdaBytes * 2);
  vector_open(vecCom.k, vecCom.com, b, pdec, com_j, numVoleInstances, lambdaBytes);

  vector_verify(pdec, com_j, b, lambda, lambdaBytes, numVoleInstances, &vecComRec, &vecCom.h);

  uint32_t outlen = 16;
  uint8_t* v      = malloc(outlen * depth);
  // TODO: we do not input veccomRec.m but instead something else defined in
  ConvertToVoleVerifier(lambda, lambdaBytes, vecComRec.m, numVoleInstances, depth, outlen, v);

// TODO: write better test cases : )
#if 0
  printf("r\n");
  for (uint32_t i = 0; i < getBinaryTreeNodeCount(numVoleInstances); i++) {
    printf("%d ", i);
    for (uint32_t j = 0; j < outlen; j++) {
      printf("%.2x", *(r + j + (i * outlen)));
    }
    printf("\n");
  }
  printf("\n");
  printf("\nv_0...d-1\n");
  for (uint32_t d = 0; d < depth; d++) {
    for (uint32_t i = 0; i < outlen; i++) {
      printf("%.2x", *(v + (d * outlen) + i));
    }
    printf("\n");
  }
#endif

  // TODO: make tests !!
  return 1;
}

int main(void) {
  if (test_ConvertToVoleProver() == 1 && test_ConvertToVoleVerifier() == 1) {
    return 0;
  } else {
    return 1;
  }
}