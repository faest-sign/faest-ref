#include "../vc.c"
#include "../fields.h"
#include "../compat.h"

int test_numrec_bitdec() {
  uint8_t expect_out_1[2] = {0x00, 0x01};
  uint8_t* b_1            = malloc(2);
  BitDec(2, 2, b_1);
  uint64_t idx_1 = NumRec(2, b_1);

  uint8_t expect_out_2[4] = {0x01, 0x01, 0x01, 0x00};
  uint8_t* b_2            = malloc(4);
  BitDec(7, 4, b_2);
  uint64_t idx_2 = NumRec(4, b_2);

  uint8_t expect_out_3[4] = {0x00, 0x01, 0x00, 0x01};
  uint8_t* b_3            = malloc(4);
  BitDec(10, 4, b_3);
  uint64_t idx_3 = NumRec(4, b_3);

  uint8_t expect_out_4[4] = {0x01, 0x00, 0x01, 0x01};
  uint8_t* b_4            = malloc(4);
  BitDec(13, 4, b_4);
  uint64_t idx_4 = NumRec(4, b_4);

  if (memcmp(b_1, &expect_out_1, 2) == 0 && idx_1 == 2 && memcmp(b_2, &expect_out_2, 4) == 0 &&
      idx_2 == 7 && memcmp(b_3, &expect_out_3, 4) == 0 && idx_3 == 10 &&
      memcmp(b_4, &expect_out_4, 4) == 0 && idx_4 == 13) {
    return 0;
  } else {
    return 1;
  }
}

int test_vector_commitment() {

  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  faest_paramset_t params = faest_get_paramset(1); // Just using the FAEST-128s
  vec_com_t vecCom;
  uint32_t numVoleInstances = (1 << params.faest_param.k0);
  vector_commitment(rootKey, &params, params.faest_param.lambda, params.faest_param.lambda / 8,
                    &vecCom, numVoleInstances);

#if 0
  uint32_t treeDepth = ceil_log2(numVoleInstances) + 1;
  uint32_t numNodes  = ((1 << treeDepth) - 1) - ((1 << (treeDepth - 1)) - numVoleInstances);
  printf("Printing h \n");
  for (uint32_t i = 0; i < params.faest_param.seclvl / 4; i++) {
    printf("%.2x ", *(vecCom.h + i));
  }
  printf("\n");
  printf("Printing k \n");
  for (uint32_t j = 0; j < numNodes; j++) {
    for (uint32_t i = 0; i < params.faest_param.seclvl / 8; i++) {
      printf("%.2x ", *((vecCom.k + i) + ((params.faest_param.seclvl / 8) * j)));
    }
    printf("\n");
  }
  printf("Printing com \n");
  for (uint32_t j = 0; j < numVoleInstances; j++) {
    for (uint32_t i = 0; i < params.faest_param.seclvl / 4; i++) {
      printf("%.2x ", *((vecCom.com + i) + (params.faest_param.seclvl / 4 * j)));
    }
    printf("\n");
  }
  printf("Printing sd \n");
  for (uint32_t j = 0; j < numVoleInstances; j++) {
    for (uint32_t i = 0; i < params.faest_param.seclvl / 8; i++) {
      printf("%.2x ", *((vecCom.sd + i) + (params.faest_param.seclvl / 8 * j)));
    }
    printf("\n");
  }
#endif

  return 0;
}

int test_vector_open_128() {
  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  faest_paramset_t params   = faest_get_paramset(1); // Just using the FAEST-128s
  vec_com_t* vecCom         = malloc(sizeof(vecCom));
  uint32_t numVoleInstances = 16;
  uint32_t lambda           = params.faest_param.lambda;
  uint32_t lambdaBytes      = lambda / 8;
  vector_commitment(rootKey, &params, lambda, lambdaBytes, vecCom, numVoleInstances);

  uint32_t leafIndex = 7;
  uint32_t depth     = ceil_log2(numVoleInstances);
  uint8_t* b         = malloc(depth);
  BitDec(leafIndex, depth, b);

  uint8_t* pdec  = malloc(depth * lambdaBytes);
  uint8_t* com_j = malloc(lambdaBytes * 2);
  vector_open(vecCom->k, vecCom->com, b, pdec, com_j, numVoleInstances, lambdaBytes);

#if 0
  printf("Hidden leaf index : %d\n", leafIndex);
  printf("b[i]\n");
  for (uint32_t i = 0; i < depth; i++) {
    printf("%.2x ", b[i]);
  }
  printf("\n");

  printf("printing all coms \n");
  for (uint32_t j = 0; j < numVoleInstances; j++) {
    for (uint32_t i = 0; i < lambdaBytes * 2; i++) {
      printf("%.2x", *(vecCom->com + (j * lambdaBytes * 2) + i));
    }
    printf("\n");
  }
  printf("printing com_j\n");
  for (uint32_t i = 0; i < lambdaBytes * 2; i++) {
    printf("%.2x", *(com_j + i));
  }
  printf("\n");

  printf("printing pdec\n");
  for (uint32_t j = 0; j < depth; j++) {
    for (uint32_t i = 0; i < lambdaBytes; i++) {
      printf("%.2x", *(pdec + i + (j * lambdaBytes)));
    }
    printf("\n");
  }
  printf("printing vecCom.k\n");
  for (uint32_t i = 0; i < getBinaryTreeNodeCount(numVoleInstances); i++) {
    for (uint32_t j = 0; j < lambdaBytes; j++) {
      printf("%.2x", *(vecCom->k + j + (i * lambdaBytes)));
    }
    printf("\n");
  }
#endif

  if (memcmp(vecCom->com + (leafIndex * lambdaBytes * 2), com_j, lambdaBytes * 2) == 0 &&
      memcmp(pdec, vecCom->k + (lambdaBytes * getNodeIndex(1, 1)), lambdaBytes) == 0 &&
      memcmp(pdec + lambdaBytes, vecCom->k + (lambdaBytes * getNodeIndex(2, 0)), lambdaBytes) ==
          0 &&
      memcmp(pdec + (lambdaBytes * 2), vecCom->k + (lambdaBytes * getNodeIndex(3, 2)),
             lambdaBytes) == 0 &&
      memcmp(pdec + (lambdaBytes * 3), vecCom->k + (lambdaBytes * getNodeIndex(4, 6)),
             lambdaBytes) == 0) {
    return 0;
  } else {
    return 1;
  }
}

int test_vector_open_192() {
  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  faest_paramset_t params = faest_get_paramset(3); // Just using the FAEST-192s
  vec_com_t* vecCom       = malloc(sizeof(vec_com_t));
  // vec_com_rec_t vecComRec;
  uint32_t numVoleInstances = 16;
  uint32_t lambda           = params.faest_param.lambda;
  uint32_t lambdaBytes      = lambda / 8;
  vector_commitment(rootKey, &params, lambda, lambdaBytes, vecCom, numVoleInstances);

  uint32_t leafIndex = 10;
  uint32_t depth     = ceil_log2(numVoleInstances);
  uint8_t* b         = malloc(depth);
  BitDec(leafIndex, depth, b);

  uint8_t* pdec  = malloc(depth * lambdaBytes);
  uint8_t* com_j = malloc(lambdaBytes * 2);
  vector_open(vecCom->k, vecCom->com, b, pdec, com_j, numVoleInstances, lambdaBytes);

#if 0
  printf("Hidden leaf index : %d\n", leafIndex);
  printf("b[i]\n");
  for (uint32_t i = 0; i < depth; i++) {
    printf("%.2x ", b[i]);
  }
  printf("\n");

  printf("printing all coms \n");
  for (uint32_t j = 0; j < params.faest_param.k0; j++) {
    for (uint32_t i = 0; i < lambdaBytes * 2; i++) {
      printf("%.2x", *(vecCom->com + (j * lambdaBytes * 2) + i));
    }
    printf("\n");
  }
  printf("printing com_j\n");
  for (uint32_t i = 0; i < lambdaBytes * 2; i++) {
    printf("%.2x", *(com_j + i));
  }
  printf("\n");

  printf("printing pdec\n");
  for (uint32_t j = 0; j < depth; j++) {
    for (uint32_t i = 0; i < lambdaBytes; i++) {
      printf("%.2x", *(pdec + i + (j * lambdaBytes)));
    }
    printf("\n");
  }
  printf("printing vecCom.k \n");
  for (uint32_t i = 0; i < getBinaryTreeNodeCount(numVoleInstances); i++) {
    for (uint32_t j = 0; j < lambdaBytes; j++) {
      printf("%.2x", *(vecCom->k + j + (i * lambdaBytes)));
    }
    printf("\n");
  }
#endif

  if (memcmp(vecCom->com + (leafIndex * lambdaBytes * 2), com_j, lambdaBytes * 2) == 0 &&
      memcmp(pdec, vecCom->k + (lambdaBytes * getNodeIndex(1, 0)), lambdaBytes) == 0 &&
      memcmp(pdec + lambdaBytes, vecCom->k + (lambdaBytes * getNodeIndex(2, 3)), lambdaBytes) ==
          0 &&
      memcmp(pdec + (lambdaBytes * 2), vecCom->k + (lambdaBytes * getNodeIndex(3, 4)),
             lambdaBytes) == 0 &&
      memcmp(pdec + (lambdaBytes * 3), vecCom->k + (lambdaBytes * getNodeIndex(4, 11)),
             lambdaBytes) == 0) {
    return 0;
  } else {
    return 1;
  }
}

int test_vector_open_256() {
  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  faest_paramset_t params   = faest_get_paramset(5); // Just using the FAEST-256s
  vec_com_t* vecCom         = malloc(sizeof(vec_com_t));
  uint32_t numVoleInstances = 16;
  uint32_t lambda           = params.faest_param.lambda;
  uint32_t lambdaBytes      = lambda / 8;
  vector_commitment(rootKey, &params, lambda, lambdaBytes, vecCom, numVoleInstances);

  uint32_t leafIndex = 7;
  uint32_t depth     = ceil_log2(numVoleInstances);
  uint8_t* b         = malloc(depth);
  BitDec(leafIndex, depth, b);

  uint8_t* pdec  = malloc(depth * lambdaBytes);
  uint8_t* com_j = malloc(lambdaBytes * 2);
  vector_open(vecCom->k, vecCom->com, b, pdec, com_j, numVoleInstances, lambdaBytes);

#if 0
  printf("Hidden leaf index : %d\n", leafIndex);
  printf("b[i]\n");
  for (uint32_t i = 0; i < depth; i++) {
    printf("%.2x ", b[i]);
  }
  printf("\n");

  printf("printing all coms \n");
  for (uint32_t j = 0; j < params.faest_param.k0; j++) {
    for (uint32_t i = 0; i < lambdaBytes * 2; i++) {
      printf("%.2x", *(vecCom->com + (j * lambdaBytes * 2) + i));
    }
    printf("\n");
  }
  printf("printing com_j\n");
  for (uint32_t i = 0; i < lambdaBytes * 2; i++) {
    printf("%.2x", *(com_j + i));
  }
  printf("\n");

  printf("printing pdec\n");
  for (uint32_t j = 0; j < depth; j++) {
    for (uint32_t i = 0; i < lambdaBytes; i++) {
      printf("%.2x", *(pdec + i + (j * lambdaBytes)));
    }
    printf("\n");
  }
  printf("printing vecCom.k \n");
  for (uint32_t i = 0; i < getBinaryTreeNodeCount(numVoleInstances); i++) {
    for (uint32_t j = 0; j < lambdaBytes; j++) {
      printf("%.2x", *(vecCom->k + j + (i * lambdaBytes)));
    }
    printf("\n");
  }
#endif

  if (memcmp(vecCom->com + (leafIndex * (lambdaBytes * 2)), com_j, lambdaBytes * 2) == 0 &&
      memcmp(pdec, vecCom->k + (lambdaBytes * getNodeIndex(1, 1)), lambdaBytes) == 0 &&
      memcmp(pdec + lambdaBytes, vecCom->k + (lambdaBytes * getNodeIndex(2, 0)), lambdaBytes) ==
          0 &&
      memcmp(pdec + (lambdaBytes * 2), vecCom->k + (lambdaBytes * getNodeIndex(3, 2)),
             lambdaBytes) == 0 &&
      memcmp(pdec + (lambdaBytes * 3), vecCom->k + (lambdaBytes * getNodeIndex(4, 6)),
             lambdaBytes) == 0) {
    return 0;
  } else {
    return 1;
  }
}

int test_vector_reconstruct_and_verify() {

  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  faest_paramset_t params   = faest_get_paramset(1); // Just using the FAEST-128s
  vec_com_t* vecCom         = malloc(sizeof(vec_com_t));
  vec_com_rec_t* vecComRec  = malloc(sizeof(vec_com_rec_t));
  uint32_t numVoleInstances = 16;
  uint32_t lambda           = params.faest_param.lambda;
  uint32_t lambdaBytes      = lambda / 8;
  vector_commitment(rootKey, &params, lambda, lambdaBytes, vecCom, numVoleInstances);

  uint32_t leafIndex = 7;
  uint32_t depth     = ceil_log2(numVoleInstances);
  uint8_t* b         = malloc(depth);
  BitDec(leafIndex, depth, b);

  uint8_t* pdec  = malloc(depth * lambdaBytes);
  uint8_t* com_j = malloc(lambdaBytes * 2);
  vector_open(vecCom->k, vecCom->com, b, pdec, com_j, numVoleInstances, lambdaBytes);

  int verify_ret =
      vector_verify(pdec, com_j, b, lambda, lambdaBytes, numVoleInstances, vecComRec, vecCom->h);

#if 0
  printf("Hidden leaf index : %d\n", leafIndex);
  printf("b[i]\n");
  for (uint32_t i = 0; i < depth; i++) {
    printf("%.2x ", b[i]);
  }
  printf("\n");

  printf("printing pdec\n");
  for (uint32_t j = 0; j < depth; j++) {
    for (uint32_t i = 0; i < lambdaBytes; i++) {
      printf("%.2x", *(pdec + i + (j * lambdaBytes)));
    }
    printf("\n");
  }

  printf("printing vecCom.k / vecComRed.k\n");
  for (uint32_t i = 0; i < getBinaryTreeNodeCount(numVoleInstances); i++) {
    for (uint32_t j = 0; j < lambdaBytes; j++) {
      printf("%.2x", *(vecCom->k + j + (i * lambdaBytes)));
    }
    printf(" ");
    for (uint32_t j = 0; j < lambdaBytes; j++) {
      printf("%.2x", *(vecComRec->k + j + (i * lambdaBytes)));
    }
    uint8_t* zeros = malloc(lambdaBytes);
    memset(zeros, 0, lambdaBytes);
    if (memcmp(vecComRec->k + (i * lambdaBytes), zeros, lambdaBytes) == 0 &&
        memcmp(vecCom->k + (i * lambdaBytes), zeros, lambdaBytes) != 0) {
      printf("<--Hidden");
    }

    printf("\n");
  }
  printf("\n");

  printf("\n");
  printf("printing vecCom.sd / vecComRed.m\n");
  for (uint32_t i = 0; i < params.faest_param.k0; i++) {
    for (uint32_t j = 0; j < lambdaBytes; j++) {
      printf("%.2x", *(vecCom->sd + j + (i * lambdaBytes)));
    }
    printf(" ");
    for (uint32_t j = 0; j < lambdaBytes; j++) {
      printf("%.2x", *(vecComRec->m + j + (i * lambdaBytes)));
    }
    if (i == leafIndex) {
      printf("<-- Hidden");
    }
    printf("\n");
  }
  printf("\n");

  printf("\n");
  printf("printing vecCom.com / vecComRed.com\n");
  for (uint32_t i = 0; i < params.faest_param.k0; i++) {
    for (uint32_t j = 0; j < lambdaBytes * 2; j++) {
      printf("%.2x", *(vecCom->com + j + (i * lambdaBytes * 2)));
    }
    printf(" ");
    for (uint32_t j = 0; j < lambdaBytes * 2; j++) {
      printf("%.2x", *(vecComRec->com + j + (i * lambdaBytes * 2)));
    }
    if (i == leafIndex) {
      printf("<-- com_j*");
    }
    printf("\n");
  }
  printf("\n");

  printf("\n");
  printf("printing vecCom.h / vecComRed.h\n");
  for (uint32_t j = 0; j < lambdaBytes * 2; j++) {
    printf("%.2x", *(vecCom->h + j));
  }
  printf(" ");
  for (uint32_t j = 0; j < lambdaBytes * 2; j++) {
    printf("%.2x", *(vecComRec->h + j));
  }
  printf("\n");

#endif

  if (verify_ret == 1) {
    return 0;
  } else {
    return 1;
  }
}

int main(void) {

  if (test_numrec_bitdec() == 0 && test_vector_commitment() == 0 && test_vector_open_128() == 0 &&
      test_vector_open_192() == 0 && test_vector_open_256() == 0 &&
      test_vector_reconstruct_and_verify() == 0) {
    printf("vc.c : All tests pass !!\n");
    return 0;
  } else {
    return 1;
  }
}