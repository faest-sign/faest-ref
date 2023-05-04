/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef INSTANCES_H
#define INSTANCES_H

#include "faest_defines.h"

#include <stddef.h>
#include <stdint.h>

// TODO - Add paramset struct faest_param_t here, containing faest_L1, faest_L3, faest_L5.... ?
typedef enum faest_paramid_t {
  PARAMETER_SET_INVALID   = 0,
  FAEST_L1_S              = 1,
  FAEST_L1_F              = 2,
  FAEST_L3_S              = 3,
  FAEST_L3_F              = 4,
  FAEST_L5_S              = 5,
  FAEST_L5_F              = 6,
  PARAMETER_SET_MAX_INDEX = 7
} faest_paramid_t;

// Param for the cipher that is used,, AES in this case
typedef struct cipher_param_t {
  uint32_t keySize;
  uint32_t blockSize;
  uint32_t numRounds;
  uint32_t numSboxes;
  uint32_t stateSizeBits;
  uint32_t stateSizeBytes;
  uint32_t stateSizeWords;
} cipher_param_t;

typedef struct faest_param_t {
  uint32_t lambda;               // lambda
  uint32_t lenExpAESWitnessBits; // l
  uint32_t t;                    // t, number of VOLE instances (number of repetition)
  uint32_t k0;                   // k0 = ceil(lambda/t)
  uint32_t k1;                   // k1 = floor(lambda/t)
  uint32_t t0;                   // t0 = lambda mod t; number if BOLE instances of size 2^(k0)
  uint32_t t1;                   // t1 = t - (lambda mod t); number if BOLE instances of size 2^(k1)
  uint32_t pkSizeBytes;
  uint32_t skSizeBytes;
  uint32_t numOpenRounds; // TODO: unsure what this would be in the vole case ?
  uint32_t seedSizeBytes;
  uint32_t saltSizeBytes;
  uint32_t digestSizeBytes;
} faest_param_t;

// TODO - Add paramset struct paramset_t here, containing numround, numsbox, seedSizeBytes,
// digestSizeBytes.... ?
typedef struct faest_paramset_t {
  cipher_param_t cipher_param;
  faest_param_t faest_param;
  faest_paramid_t faest_paramid;
} faest_paramset_t;

const char* faest_get_param_name(faest_paramid_t paramid);
int faest_check_paramset(faest_paramset_t* paramset);
faest_paramset_t faest_get_paramset(faest_paramid_t paramid);

#endif