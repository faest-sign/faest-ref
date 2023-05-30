#ifndef FAEST_VOLE_H
#define FAEST_VOLE_H

#include "vc.h"
#include <stdbool.h>

FAEST_BEGIN_C_DECL

// k_b is at most 12, so chalout needs to point to an array of at most 12 bytes
int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout);

void voleCommit(const uint8_t* rootKey, const uint8_t* iv, uint32_t ellhat,
                const faest_paramset_t* params, uint8_t* hcom, vec_com_t* vecCom, uint8_t** c,
                uint8_t* u, uint8_t** v);

void voleReconstruct(const uint8_t* iv, const uint8_t* chal, uint8_t** pdec, uint8_t** com_j,
                     uint8_t* hcom, uint8_t** q, uint32_t ellhat, const faest_paramset_t* params);

void ConvertToVole(const uint8_t* iv, const uint8_t* sd, bool sd0_bot, uint32_t lambda,
                   uint32_t lambdaBytes, uint32_t numVoleInstances, uint32_t depth,
                   uint32_t outLenBytes, uint8_t* u, uint8_t* v);

FAEST_END_C_DECL

#endif
