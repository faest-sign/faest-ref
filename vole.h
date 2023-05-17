#include "vc.h"

void ConvertToVoleProver(uint32_t lambda, const uint8_t* sd, uint32_t numVoleInstances,
                         uint32_t depth, uint32_t outLenBytes, uint8_t* u, uint8_t* v);

void ConvertToVoleVerifier(uint32_t lambda, const uint8_t* sd, uint32_t numVoleInstances,
                           uint32_t depth, uint32_t outLenBytes, uint8_t* v);