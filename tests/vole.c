#include "../vole.h"

int test_ChalDec() {
  faest_paramset_t params = faest_get_paramset(1);
  uint32_t k0             = params.faest_param.k0;
  uint32_t t0             = params.faest_param.t0;
  uint32_t k1             = params.faest_param.k1;
  uint32_t t1             = params.faest_param.t1;
  uint8_t* chal           = malloc((k0 * t0) + (k1 * t1));
  uint8_t* chalout        = malloc((k0 * t0) + (k1 * t1));
  uint32_t idx            = 5; // between [0..t0+t1)

  for (uint32_t i = 0; i < (k0 * t0) + (k1 * t1); i++) {
    *(chal + i) = rand() % 2;
  }
  ChalDec(chal, idx, k0, t0, k1, t1, chalout);

  if (memcmp(chalout, chal + 60, 11) == 0) {
    return 1;
  } else {
    return 0;
  }
}

int test_FAESTVoleCommit() {

  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  faest_paramset_t params = faest_get_paramset(1); // Just using the FAEST-128s
                                                   //   vec_com_t vecCom;
                                                   //   vec_com_rec_t vecComRec;
  uint32_t outlen = 16;

  uint8_t* hcom      = malloc(params.faest_param.lambda / 8);
  vec_com_t** vecCom = malloc(params.faest_param.t * (sizeof(vec_com_t*)));
  uint8_t** c        = malloc((params.faest_param.t * sizeof(uint8_t*)) - 1);
  uint8_t* u         = malloc(outlen);
  uint8_t** v        = malloc(params.faest_param.t * sizeof(uint8_t*));
  voleCommit(rootKey, params.faest_param.lambda, params.faest_param.lambdaBytes, outlen, &params,
             hcom, vecCom, c, u, v);

  // TODO: make tests !!
  return 1;
}

int test_FAESTVoleVerify() {

  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  faest_paramset_t params = faest_get_paramset(1); // Just using the FAEST-128s

  uint32_t outlen = 16;

  uint32_t lambda      = params.faest_param.lambda;
  uint32_t lambdaBytes = params.faest_param.lambdaBytes;
  uint8_t* hcom        = malloc(lambdaBytes * 2);
  vec_com_t** vecCom   = malloc(params.faest_param.t * (sizeof(vec_com_t*)));
  uint8_t** c          = malloc((params.faest_param.t * sizeof(uint8_t*)) - 1);
  uint8_t* u           = malloc(outlen);
  uint8_t** v          = malloc(params.faest_param.t * sizeof(uint8_t*));

  voleCommit(rootKey, lambda, lambdaBytes, outlen, &params, hcom, vecCom, c, u, v);

  // TODO: this shouldn't be here !!
  uint8_t** pdec            = malloc(params.faest_param.t * sizeof(uint8_t*));
  uint8_t** b               = malloc(params.faest_param.t * sizeof(uint8_t*));
  uint8_t** com_j           = malloc(params.faest_param.t * sizeof(uint8_t*));
  vec_com_rec_t** vecComRec = malloc(params.faest_param.t * sizeof(vec_com_rec_t*));

  uint32_t depth;
  uint32_t numVoleInstances;
  for (uint32_t i = 0; i < params.faest_param.t; i++) {
    if (i < params.faest_param.t0) {
      depth            = params.faest_param.k0;
      numVoleInstances = 1 << depth;
    } else {
      depth            = params.faest_param.k1;
      numVoleInstances = 1 << depth;
    }
    b[i] = malloc(depth);
    memset(b[i], 0, depth); // TODO: LOOKOUT !! always opening the first leaf for this test
    pdec[i]  = malloc(depth * lambdaBytes);
    com_j[i] = malloc(lambdaBytes * 2);

    vector_open(vecCom[i]->k, vecCom[i]->com, b[i], pdec[i], com_j[i], numVoleInstances,
                lambdaBytes);
  }

  uint8_t* chal = malloc((params.faest_param.k0 * params.faest_param.t0) +
                         (params.faest_param.k1 * params.faest_param.t1));
  // TODO: LOOKOUT !! always setting it to 0s for testing
  memset(chal, 0,
         (params.faest_param.k0 * params.faest_param.t0) +
             (params.faest_param.k1 * params.faest_param.t1));
  uint8_t** q = malloc(params.faest_param.t * sizeof(uint8_t*));

  uint8_t* hcomRec = malloc(lambdaBytes * 2);

  voleVerify(chal, pdec, com_j, lambda, lambdaBytes, outlen, params.faest_param.t,
             params.faest_param.k0, params.faest_param.k1, hcomRec, q, vecComRec);

#if 0
  for (uint z = 0; z < params.faest_param.t; z++) {
    printf("%d) printing vecCom.com / vecComRed.com\n", z);
    for (uint32_t i = 0; i < params.faest_param.k0; i++) {
      for (uint32_t j = 0; j < lambdaBytes * 2; j++) {
        printf("%.2x", *(vecCom[z]->com + j + (i * lambdaBytes * 2)));
      }
      printf(" ");
      for (uint32_t j = 0; j < lambdaBytes * 2; j++) {
        printf("%.2x", *(vecComRec[z]->com + j + (i * lambdaBytes * 2)));
      }
      if (i == 0) {
        printf("<-- com_j*");
      }
      printf("\n");
    }
  }
  printf("\n");
  for (uint32_t i = 0; i < lambdaBytes; i++) {
    printf("%.2x", *(hcom + i));
  }
  printf(" ");
  for (uint32_t i = 0; i < lambdaBytes; i++) {
    printf("%.2x", *(hcomRec + i));
  }
  printf("\n");
#endif
  if (memcmp(hcom, hcomRec, lambdaBytes) == 0) {
    return 1;
  }
}

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
  if (test_ConvertToVoleProver() == 1 && test_ConvertToVoleVerifier() == 1 && test_ChalDec() == 1 &&
      test_FAESTVoleCommit() == 1 && test_FAESTVoleVerify() == 1) {
    return 0;
  } else {
    return 1;
  }
}