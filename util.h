#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include "instances.h"
#include "vbb.h"

void print_chall1(uint8_t* chall1, int lambda_bytes){
    printf("### chall 1 ###\n");
    for(int i = 0; i < (5 * lambda_bytes) + 8; i++){
        printf("%02x", chall1[i]);
    }
    printf("\n");
}

void print_chall2(uint8_t* chall2, int lambda_bytes){
    printf("### chall 2 ###\n");
    for(int i = 0; i < 3 * lambda_bytes + 8; i++){
        printf("%02x", chall2[i]);
    }
    printf("\n");
}

void print_b_tilde(uint8_t* b_tilde, int lambda_bytes){
    printf("### b_tilde ###\n");
    for(int i = 0; i < lambda_bytes; i++){
        printf("%02x", b_tilde[i]);
    }
    printf("\n");
}

void print_hv(uint8_t* hv, int lambda_bytes){
    printf("### hv ###\n");
    for(int i = 0; i < lambda_bytes*2; i++){
        printf("%02x", hv[i]);
    }
    printf("\n");
}

void compare_OLEs_cmo(vbb_t* vbb_full, vbb_t* vbb_nonfull){
    prepare_hash_sign(vbb_full);
    printf("### Comparing OLEs CMO ###\n");
    const unsigned int l           = vbb_nonfull->params->faest_param.l;
    const unsigned int lambda      = vbb_nonfull->params->faest_param.lambda;
    const unsigned int ell_hat     = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = (ell_hat + 7) / 8;
    for(unsigned int i = 0; i < lambda; i++){
        const uint8_t* OLE_full = get_vole_v_hash(vbb_full, i);
        const uint8_t* OLE_nonfull = get_vole_v_hash(vbb_nonfull, i);
        for(unsigned int j = 0; j < ell_hat_bytes; j++){
            if(OLE_full[j] != OLE_nonfull[j]){
                printf("OLEs CMO differ at index %d\n", i);
                break;
            }
        }
    }
}

void compare_OLEs(vbb_t* vbb_full, vbb_t* vbb_nonfull){
    prepare_aes_sign(vbb_full);
    printf("### Comparing OLEs ###\n");
    const unsigned int l           = vbb_full->params->faest_param.l;
    const unsigned int lambda      = vbb_full->params->faest_param.lambda;
    const unsigned int ell_hat     = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
    for(unsigned int i = 0; i < ell_hat; i++){
        uint8_t* OLE_full, *OLE_nonfull;
        if(lambda == 128){
            OLE_full = (uint8_t*)get_vole_aes_128(vbb_full, i);
            OLE_nonfull = (uint8_t*)get_vole_aes_128(vbb_nonfull, i);
        } else if(lambda == 192){
            OLE_full = (uint8_t*)get_vole_aes_192(vbb_full, i);
            OLE_nonfull = (uint8_t*)get_vole_aes_192(vbb_nonfull, i);
        } else if(lambda == 256){
            OLE_full = (uint8_t*)get_vole_aes_256(vbb_full, i);
            OLE_nonfull = (uint8_t*)get_vole_aes_256(vbb_nonfull, i);
        }
        for(unsigned int j = 0; j < lambda/8; j++){
            if(OLE_full[j] != OLE_nonfull[j]){
                printf("OLEs differ at index %d\n", i);
                break;
            }
        }
    }
}

#endif // UTIL_H
