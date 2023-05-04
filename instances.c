/*
 *  SPDX-License-Identifier: MIT
 */

#include "instances.h"

const char* FAEST_GET_PARAM_NAME(faest_paramid_t paramid) {
  switch (paramid) {
  case PARAMETER_SET_INVALID:
    return "PARAMETER_SET_INVALID";
  case FAEST_L1_S:
    return "FAEST_L1_S";
  case FAEST_L1_F:
    return "FAEST_L1_F";
  case FAEST_L3_S:
    return "FAEST_L3_S";
  case FAEST_L3_F:
    return "FAEST_L3_F";
  case FAEST_L5_S:
    return "FAEST_L5_S";
  case FAEST_L5_F:
    return "FAEST_L5_F";
  default:
    return "PARAMETER_SET_MAX_INDEX";
  }
}

int FAEST_CHECK_PARAMSET(faest_paramset_t* paramset) {
  uint32_t tmp = (paramset->faest_param.k0 * paramset->faest_param.t0) +
                 (paramset->faest_param.k1 * paramset->faest_param.t1);
  if (tmp == paramset->faest_param.lambda)
    return 1;
  return -1;
}

// keySize, blockSize, numRounds, numSboxes, stateSizeBits, stateSizeBytes, stateSizeWords
const cipher_param_t AES_128_PARAMS     = {16, 16, 10, 200, 128, 16, 4};
const cipher_param_t AES_192_PARAMS     = {24, 16, 12, 220, 128, 16, 4};
const cipher_param_t AES_256_PARAMS     = {32, 16, 14, 240, 128, 16, 4};
const cipher_param_t AES_INVALID_PARAMS = {0, 0, 0, 0, 0, 0, 0};
// TODO: Later
// const cipher_param_t AES_128_EM_PARAMS = {16,16,10,160,128,16,4};

// TODO: Adapt the last three params (seedSizeBytes, saltSizeBytes, digestSizeBytes)
// TODO: Number of open rounds, number of MPC rounds
// lambda, lenExpAESWitnessBits, t, k0, k1, t0, t1, pkSizeBytes, skSizeBytes, uint32_t
// numOpenRounds,
//                                                          seedSizeBytes, saltSizeBytes,
//                                                          digestSizeBytes
const faest_param_t FAEST_128_S_PARAMS   = {128, 1600, 11, 12, 11, 7, 4, 32, 16, 0, 16, 16, 16};
const faest_param_t FAEST_128_F_PARAMS   = {128, 1600, 16, 8, 8, 0, 16, 32, 16, 0, 16, 16, 16};
const faest_param_t FAEST_192_S_PARAMS   = {192, 3328, 16, 12, 12, 0, 16, 32, 16, 0, 16, 16, 16};
const faest_param_t FAEST_192_F_PARAMS   = {192, 3328, 24, 8, 8, 0, 24, 32, 16, 0, 16, 16, 16};
const faest_param_t FAEST_256_S_PARAMS   = {256, 4000, 22, 12, 11, 14, 8, 32, 16, 0, 16, 16, 16};
const faest_param_t FAEST_256_F_PARAMS   = {256, 4000, 32, 8, 8, 0, 32, 32, 16, 0, 16, 16, 16};
const faest_param_t FAEST_INVALID_PARAMS = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static const faest_paramset_t faestInstances[PARAMETER_SET_MAX_INDEX] = {
    {AES_INVALID_PARAMS, FAEST_INVALID_PARAMS, PARAMETER_SET_INVALID},
    {AES_128_PARAMS, FAEST_128_S_PARAMS, FAEST_L1_S},
    {AES_128_PARAMS, FAEST_128_F_PARAMS, FAEST_L1_F},
    {AES_192_PARAMS, FAEST_192_S_PARAMS, FAEST_L3_S},
    {AES_192_PARAMS, FAEST_192_F_PARAMS, FAEST_L3_F},
    {AES_256_PARAMS, FAEST_256_S_PARAMS, FAEST_L5_S},
    {AES_256_PARAMS, FAEST_256_F_PARAMS, FAEST_L5_F}};

faest_paramset_t FAEST_GET_PARAMSET(faest_paramid_t paramid) {

  switch (paramid) {
  case FAEST_L1_S:
    return faestInstances[1];
  case FAEST_L1_F:
    return faestInstances[2];
  case FAEST_L3_S:
    return faestInstances[3];
  case FAEST_L3_F:
    return faestInstances[4];
  case FAEST_L5_S:
    return faestInstances[5];
  case FAEST_L5_F:
    return faestInstances[6];
  default:
    return faestInstances[0];
  }
}
