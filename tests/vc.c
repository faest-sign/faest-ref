#include "../vc.c"
#include "../fields.h"
#include "../compat.h"

int test_numrec_bitdec() {
  uint8_t expect_out_1[2] = {0x00, 0x01};
  uint8_t* b_1            = malloc(2);
  BitDec(2, 2, b_1);
  uint64_t* idx_1 = malloc(sizeof(uint64_t));
  NumRec(2, b_1, idx_1);

  uint8_t expect_out_2[4] = {0x01, 0x01, 0x01, 0x00};
  uint8_t* b_2            = malloc(4);
  BitDec(7, 4, b_2);
  uint64_t* idx_2 = malloc(sizeof(uint64_t));
  NumRec(4, b_2, idx_2);

  uint8_t expect_out_3[4] = {0x00, 0x01, 0x00, 0x01};
  uint8_t* b_3            = malloc(4);
  BitDec(10, 4, b_3);
  uint64_t* idx_3 = malloc(sizeof(uint64_t));
  NumRec(4, b_3, idx_3);

  uint8_t expect_out_4[4] = {0x01, 0x00, 0x01, 0x01};
  uint8_t* b_4            = malloc(4);
  BitDec(13, 4, b_4);
  uint64_t* idx_4 = malloc(sizeof(uint64_t));
  NumRec(4, b_4, idx_4);

  if (memcmp(b_1, &expect_out_1, 2) == 0 && *idx_1 == 2 && memcmp(b_2, &expect_out_2, 4) == 0 &&
      *idx_2 == 7 && memcmp(b_3, &expect_out_3, 4) == 0 && *idx_3 == 10 &&
      memcmp(b_4, &expect_out_4, 4) == 0 && *idx_4 == 13) {
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
  tree_t* tree;
  vector_commitment(rootKey, &params, &vecCom, tree);

#if 0
  uint32_t treeDepth = ceil_log2(params.faest_param.t) + 1;
  uint32_t numNodes  = ((1 << treeDepth) - 1) - ((1 << (treeDepth - 1)) - params.faest_param.t);
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
  for (uint32_t j = 0; j < params.faest_param.t; j++) {
    for (uint32_t i = 0; i < params.faest_param.seclvl / 4; i++) {
      printf("%.2x ", *((vecCom.com + i) + (params.faest_param.seclvl / 4 * j)));
    }
    printf("\n");
  }
  printf("Printing sd \n");
  for (uint32_t j = 0; j < params.faest_param.t; j++) {
    for (uint32_t i = 0; i < params.faest_param.seclvl / 8; i++) {
      printf("%.2x ", *((vecCom.sd + i) + (params.faest_param.seclvl / 8 * j)));
    }
    printf("\n");
  }
#endif

  /* NEED TO CROSS REF THIS WITH THE OTHER IMPLEMENTATION USING THE SAME PARAMS
  Printing h
  c2 4b 32 61 53 b4 fd b5 1c 59 ff ba 82 b0 55 38 76 3d 80 19 46 e6 47 f9 17 11 f0 b2 94 0b c3 c3
  Printing k
  00 01 02 03 04 05 06 07 08 09 0a 0b 0c 0d 0e 0f
  c6 a1 3b 37 87 8f 5b 82 6f 4f 81 62 a1 c8 d8 79
  73 46 13 95 95 c0 b4 1e 49 7b bd e3 65 f4 2d 0a
  2c 57 8f 79 27 a9 49 d3 b5 11 ae 8f b6 91 45 c6
  b7 5b 1a 66 b8 a4 21 3a b3 f5 d7 3e 3b a9 8a 87
  cd bd 38 92 5b e0 eb d4 ed db 4a ea bc d4 ef 6a
  00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
  7f d3 3c 93 31 62 41 be 4b e3 3f a2 1e b6 64 1c
  66 80 4f a3 a1 3a 7e 39 1c a2 cd e3 7c 7c 9e cf
  b1 28 c1 c4 cb 33 03 a0 07 6e e3 6d 47 30 58 ab
  d2 0d 33 dd ea b9 d7 f8 21 5b d1 5d d7 34 4c ea
  93 2c ed ba 96 80 d9 40 41 d7 34 3b a8 5d 97 e0
  a2 64 06 0c 84 ac 85 1e 1f 58 ee 8b 00 cd 55 cb
  00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
  00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
  83 84 e6 cd 73 58 8b b3 ba 12 0f b0 86 fe 4c fc
  53 81 5c 98 70 fa bc dc e3 25 1a e9 ba a1 0d dd
  e7 10 19 b7 88 81 34 0c bf 8e 82 6c 6e d6 3b c5
  81 90 d9 7a 1e db 75 95 22 5a 77 00 2d 04 e3 21
  7c fe dd a0 6e 53 b0 8a f0 18 95 a7 89 cd 36 ff
  b1 51 c3 3f 98 01 1f 33 0b 0b 2a 94 60 3f 98 80
  10 e3 7c 54 5d 91 e9 d2 35 a1 45 88 de 4a 9d 3e
  6d d0 cf 97 00 51 33 e4 b8 4f 29 91 87 46 5c 36
  de ff b0 ae f4 46 e3 3c 6d 0b e1 aa cb 73 4d f2
  81 85 4e fd de ee 7f 59 bf a8 c8 06 c3 cb d4 45
  73 6d b9 83 a7 90 53 1a d4 e6 a1 7d cb 9c cb 98
  Printing com
  0f bc 54 0a c5 0e ec 23 c7 9d a9 3a b9 54 fa 30 e7 72 6b a8 a8 c5 fb 6c 7b 84 49 cf 5d f1 8d 5a
  c6 e8 7c 41 86 53 0b f5 ea 38 e1 34 3d b3 e7 2f bb 80 4f 1a a0 5c 32 ea 11 48 6f 11 4c 77 1f 11
  51 f6 22 2b 53 ff da 3d cd 4d 18 97 d8 da 36 37 0d a9 28 3e 23 90 c9 45 c6 64 73 46 54 e3 27 f6
  bf ab 9f 2b 00 d1 68 07 4c 08 1d 63 af 9b e9 7f 70 84 c2 f1 ba 9f e9 1b 2a 12 64 53 a0 28 9e 4a
  a2 30 84 3c d7 7d 62 67 12 52 91 70 a7 88 0e 6e 27 6b 0d 52 79 9f 9c 00 8f c8 a3 4d ee 2f f3 30
  0c 5c 00 3b 18 eb e7 0f 8a d4 04 c5 f6 74 c6 8b 1f e1 d1 d3 f1 58 e8 16 d6 a3 99 6d fc d2 15 d1
  82 52 0f 12 85 b2 43 3b aa cb d7 ab 77 7a 8f 5f 10 0a c7 f7 53 42 4f 23 e8 10 4b b1 1f 47 7d 20
  87 dd 6d e7 4e 5b 58 4f 3f dc 3e 65 7d 9f 6e 20 77 06 58 67 7a c7 3b b5 3b 89 ca 23 e5 bc 01 c8
  f9 f0 3c 50 04 c1 eb 4c bb ba 82 dd e6 b5 cd 39 40 89 80 b0 16 26 5a 03 11 88 92 15 aa 10 66 e1
  aa 36 1e ba cc 74 ac bf 67 d7 8d 02 09 31 88 8a d3 fc 28 6f 46 53 72 ff 83 28 2c 16 92 fb b3 01
  c1 86 88 f0 42 1c 28 41 6f 4c 23 71 a4 f8 0b 43 e8 85 e5 15 17 7b 5a cb 8d c1 e4 1e f1 8b 48 f8
  Printing sd
  be e0 8f 5b 0b 46 5d 71 0b 53 e5 cf 68 a5 20 a4
  1f 7d f2 78 a4 4c e9 40 f7 93 26 4a 4e 76 bf d4
  af be 00 ed 94 ec 86 7a 92 82 e6 ee 59 54 df cc
  ed d8 aa 21 92 3d af 23 4e 87 27 67 6f 46 86 db
  c7 ad d8 a2 c4 f8 ef b6 82 ba f3 3c 17 88 c9 41
  97 7b b7 8e b3 4b f6 18 65 e4 83 1c 7b 81 4c 4a
  bd 52 1b d5 a5 b2 22 e0 d6 33 44 cd 2c bd 76 0e
  a5 a1 36 85 3f 2a 3f 83 82 64 2e 8f 5b be 3b 2b
  99 71 47 39 09 65 94 a6 40 5d 5f 01 2c dc 64 7f
  53 82 4e 89 58 ef f2 97 5e 04 eb 45 86 bb b7 97
  33 1a 30 08 60 b2 d2 4c 3f f2 0a 16 80 34 3b 77
    */

  return 0;
}

int test_vector_open() {
  uint8_t rootKey[32] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a,
                         0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15,
                         0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

  faest_paramset_t params = faest_get_paramset(1); // Just using the FAEST-128s
  vec_com_t vecCom;
  tree_t* tree;
  vector_commitment(rootKey, &params, &vecCom, tree);

  uint32_t leafIndex = 7;
  uint32_t depth     = ceil_log2(params.faest_param.t);
  uint8_t* b         = malloc(depth);
  BitDec(leafIndex, depth, b);

  for (uint32_t i = 0; i < depth; i++) {
    printf("%.2x ", b[i]);
  }
  printf("\n");

  uint32_t lambda  = params.faest_param.seclvl / 8;
  uint32_t lambda2 = params.faest_param.seclvl / 8;
  uint8_t* pdec    = malloc(depth * lambda);
  uint8_t* com_j   = malloc(lambda2);
  vector_open(&params, vecCom.k, vecCom.com, b, leafIndex, pdec, com_j);

  for (uint32_t i = 0; i < lambda2; i++) {
    printf("%.2x ", *(com_j + i));
  }
  printf("\n");

  for (uint32_t j = 0; j < depth; j++) {
    for (uint32_t i = 0; i < lambda; i++) {
      printf("%.2x", *(pdec + i + (j * lambda)));
    }
    printf("\n");
  }

  // TODO Cross verify with the tree nodes : )

  return 0;
}

int main(void) {
  if (test_vector_open()) {
    return 1;
  } else {
    return 0;
  }
}