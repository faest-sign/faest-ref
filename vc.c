#include "vc.h"
#include "tree.c"
#include "aes.h"

#include <math.h>

void vector_commitment(const uint8_t* rootKey, const faest_paramset_t* params, const uint8_t* salt, 
                            vec_com_t* vecCom) {

    tree_t* tree; 
    tree = generateSeeds(rootKey, params);
    uint8_t** buffer;
    buffer = malloc(params->faest_param.t * sizeof(uint8_t*));

    vecCom->dk = malloc(tree->numNodes * sizeof(uint8_t*));
    vecCom->dcom = malloc(params->faest_param.t * sizeof(uint8_t*));
    vecCom->sd = malloc(params->faest_param.t * sizeof(uint8_t*));   


    /* Doing H_0 */
    uint8_t** leaves = getLeaves(tree);
    for(uint32_t i = 0; i < params->faest_param.t; i++) {
        hash_context ctx;
        hash_init(&ctx, 128);
        // TODO - Concat i to the message input of the hash
        hash_update(&ctx, leaves, params->faest_param.seedSizeBytes);
        hash_update(&ctx, salt, params->faest_param.saltSizeBytes);
        hash_final(&ctx);
        hash_squeeze(&ctx, buffer[i], params->faest_param.h0digestSizeBytes);

        memcpy(vecCom->sd[i], buffer[i], params->faest_param.seclvl);
        memcpy(vecCom->dcom[i], buffer[i] + params->faest_param.seclvl, 2*params->faest_param.seclvl);
    }

    /* Doing H_1 */


    // key128_t* k;
    // k = malloc(sizeof(key128_t) * getBinaryTreeNodeCount(treeDepth));
    // memcpy(k[0], rootKey, sizeof(key128_t));
    // for(uint64_t d = 1; d <= treeDepth; d++) {
    //     for (uint64_t j = 0; j < (2 << d-1); j++) {
    //         bf8_t* output = aes_ctr_prg(k[getNodeIndex(d-1,j)], iv, seclvl);
    //         memcpy(k[getNodeIndex(d,2*j)],&output,seclvl/8);
    //         memcpy(k[getNodeIndex(d,(2*j)+1)],&output[seclvl/8],seclvl/8);
    //     }
    // }

    vec_com_t out;
    return out;
}