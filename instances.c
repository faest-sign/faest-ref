/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(WITH_CONFIG_H)
#include <config.h>
#else
#error "FIXME: Instances are currently only listed in config.h"
#endif

#include "instances.h"

const char* faest_get_param_name(faest_paramid_t paramid) {
  switch (paramid) {
  case PARAMETER_SET_INVALID:
    return "PARAMETER_SET_INVALID";
  case FAEST_128S:
    return "FAEST_128S";
  case FEAST_128F:
    return "FEAST_128F";
  case FAEST_192S:
    return "FAEST_192S";
  case FAEST_192F:
    return "FAEST_192F";
  case FAEST_256S:
    return "FAEST_256S";
  case FAEST_256F:
    return "FAEST_256F";
  default:
    return "PARAMETER_SET_MAX_INDEX";
  }
}

int faest_check_paramset(faest_paramset_t* paramset) {
  uint32_t tmp = (paramset->faest_param.k0 * paramset->faest_param.t0) +
                 (paramset->faest_param.k1 * paramset->faest_param.t1);
  if (tmp == paramset->faest_param.lambda)
    return 1;
  return -1;
}

// keySize, blockSize, numRounds, numSboxes, stateSizeBits, stateSizeBytes, stateSizeWords
static const cipher_param_t AES_128_PARAMS     = {16, 16, 10, 200, 128, 16, 4};
static const cipher_param_t AES_192_PARAMS     = {24, 16, 12, 220, 128, 16, 4};
static const cipher_param_t AES_256_PARAMS     = {32, 16, 14, 240, 128, 16, 4};
static const cipher_param_t AES_INVALID_PARAMS = {0, 0, 0, 0, 0, 0, 0};
// TODO: Later
// const cipher_param_t AES_128_EM_PARAMS = {16,16,10,160,128,16,4};

// TODO: Adapt the last three params (seedSizeBytes, saltSizeBytes, digestSizeBytes)
// TODO: Number of open rounds, number of MPC rounds
// lambda, lenExpAESWitnessBits, t, k0, k1, t0, t1, pkSizeBytes, skSizeBytes, uint32_t
// numOpenRounds,
//                                                          seedSizeBytes, saltSizeBytes,
//                                                          digestSizeBytes
static const faest_param_t FAEST_128_S_PARAMS = {FAEST_128S_LAMBDA,
                                                 FAEST_128S_ELL,
                                                 FAEST_128S_TAU,
                                                 FAEST_128S_K_0,
                                                 FAEST_128S_K_1,
                                                 FAEST_128S_T_0,
                                                 FAEST_128S_T_1,
                                                 FAEST_128S_PK_SIZE,
                                                 FAEST_128S_SK_SIZE,
                                                 0,
                                                 16,
                                                 16,
                                                 FAEST_128S_SIG_SIZE};
static const faest_param_t FAEST_128_F_PARAMS = {FAEST_128F_LAMBDA,
                                                 FAEST_128F_ELL,
                                                 FAEST_128F_TAU,
                                                 FAEST_128F_K_0,
                                                 FAEST_128F_K_1,
                                                 FAEST_128F_T_0,
                                                 FAEST_128F_T_1,
                                                 FAEST_128F_PK_SIZE,
                                                 FAEST_128F_SK_SIZE,
                                                 0,
                                                 16,
                                                 16,
                                                 FAEST_128S_SIG_SIZE};
static const faest_param_t FAEST_192_S_PARAMS = {FAEST_192S_LAMBDA,
                                                 FAEST_192S_ELL,
                                                 FAEST_192S_TAU,
                                                 FAEST_192S_K_0,
                                                 FAEST_192S_K_1,
                                                 FAEST_192S_T_0,
                                                 FAEST_192S_T_1,
                                                 FAEST_192S_PK_SIZE,
                                                 FAEST_192S_SK_SIZE,
                                                 0,
                                                 16,
                                                 16,
                                                 FAEST_192S_SIG_SIZE};
static const faest_param_t FAEST_192_F_PARAMS = {FAEST_192F_LAMBDA,
                                                 FAEST_192F_ELL,
                                                 FAEST_192F_TAU,
                                                 FAEST_192F_K_0,
                                                 FAEST_192F_K_1,
                                                 FAEST_192F_T_0,
                                                 FAEST_192F_T_1,
                                                 FAEST_192F_PK_SIZE,
                                                 FAEST_192F_SK_SIZE,
                                                 0,
                                                 16,
                                                 16,
                                                 FAEST_192S_SIG_SIZE};
static const faest_param_t FAEST_256_S_PARAMS = {FAEST_256S_LAMBDA,
                                                 FAEST_256S_ELL,
                                                 FAEST_256S_TAU,
                                                 FAEST_256S_K_0,
                                                 FAEST_256S_K_1,
                                                 FAEST_256S_T_0,
                                                 FAEST_256S_T_1,
                                                 FAEST_256S_PK_SIZE,
                                                 FAEST_256S_SK_SIZE,
                                                 0,
                                                 16,
                                                 16,
                                                 FAEST_256S_SIG_SIZE};
static const faest_param_t FAEST_256_F_PARAMS = {FAEST_256F_LAMBDA,
                                                 FAEST_256F_ELL,
                                                 FAEST_256F_TAU,
                                                 FAEST_256F_K_0,
                                                 FAEST_256F_K_1,
                                                 FAEST_256F_T_0,
                                                 FAEST_256F_T_1,
                                                 FAEST_256F_PK_SIZE,
                                                 FAEST_256F_SK_SIZE,
                                                 0,
                                                 16,
                                                 16,
                                                 FAEST_256S_SIG_SIZE};

static const faest_param_t FAEST_INVALID_PARAMS = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static const faest_paramset_t faestInstances[PARAMETER_SET_MAX_INDEX] = {
    {AES_INVALID_PARAMS, FAEST_INVALID_PARAMS, PARAMETER_SET_INVALID},
    {AES_128_PARAMS, FAEST_128_S_PARAMS, FAEST_128S},
    {AES_128_PARAMS, FAEST_128_F_PARAMS, FEAST_128F},
    {AES_192_PARAMS, FAEST_192_S_PARAMS, FAEST_192S},
    {AES_192_PARAMS, FAEST_192_F_PARAMS, FAEST_192F},
    {AES_256_PARAMS, FAEST_256_S_PARAMS, FAEST_256S},
    {AES_256_PARAMS, FAEST_256_F_PARAMS, FAEST_256F}};

faest_paramset_t faest_get_paramset(faest_paramid_t paramid) {
  switch (paramid) {
  case FAEST_128S:
    return faestInstances[1];
  case FEAST_128F:
    return faestInstances[2];
  case FAEST_192S:
    return faestInstances[3];
  case FAEST_192F:
    return faestInstances[4];
  case FAEST_256S:
    return faestInstances[5];
  case FAEST_256F:
    return faestInstances[6];
  default:
    return faestInstances[0];
  }
}
