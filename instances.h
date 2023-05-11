/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef INSTANCES_H
#define INSTANCES_H

#include "macros.h"

#include <stddef.h>
#include <stdint.h>
FAEST_BEGIN_C_DECL

// TODO - Add paramset struct faest_param_t here, containing faest_L1, faest_L3, faest_L5.... ?
typedef enum faest_paramid_t {
  PARAMETER_SET_INVALID   = 0,
  FAEST_128S              = 1,
  FAEST_128F              = 2,
  FAEST_192S              = 3,
  FAEST_192F              = 4,
  FAEST_256S              = 5,
  FAEST_256F              = 6,
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
  uint32_t pkSizeBytes;          // pk part of the signature
  uint32_t skSizeBytes;          // sk part of the signature
  uint32_t numOpenRounds;        // round opened
  uint32_t seclvl;               // security level L1, L3 or L5
  uint32_t seedSizeBytes;        // size of the seed in bytes
  uint32_t saltSizeBytes;     // size of the salt in bytes (maybe not required !!) TODO remove it ??
  uint32_t h0digestSizeBytes; // size of the digest size in bytes for H_0
  uint32_t h1digestSizeBytes; // size of the digest size in bytes for H_0
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

FAEST_END_C_DECL

#endif
