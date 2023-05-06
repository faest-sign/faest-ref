#include "vc.h"
#include "tree.c"
#include "aes.h"

#include <math.h>

vec_com_t vector_commitment(uint64_t treeDepth, uint32_t seclvl, key128_t key, bf8_t* iv) {
    uint64_t N = 2 << treeDepth;

    key128_t* k;
    k = malloc(sizeof(key128_t) * getBinaryTreeNodeCount(treeDepth));
    memcpy(k[0], key, sizeof(key128_t));

    for(uint64_t d = 1; d <= treeDepth; d++) {
        for (uint64_t j = 0; j < (2 << d-1); j++) {

            bf8_t* output = aes_ctr_prg(k[getNodeIndex(d-1,j)], iv, seclvl);

            memcpy(k[getNodeIndex(d,2*j)],&output,seclvl/8);

            memcpy(k[getNodeIndex(d,(2*j)+1)],&output[seclvl/8],seclvl/8);

            // for(uint64_t i = 0; i < 16; i++) {
            //     printf("%.2x ", k[getNodeIndex(d-1,j)][i]);
            // }
            // printf("\n");
        }
    }

    vec_com_t out;
    return out;
}