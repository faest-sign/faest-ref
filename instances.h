/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef INSTANCES_H
#define INSTANCES_H

#include "macros.h"

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#define MAX_LAMBDA 256
#define MAX_LAMBDA_BYTES (MAX_LAMBDA / 8)
#define MAX_DEPTH 12
#define MAX_TAU 32
#define UNIVERSAL_HASH_B_BITS 16
#define UNIVERSAL_HASH_B (UNIVERSAL_HASH_B_BITS / 8)
#define IV_SIZE 16

FAEST_BEGIN_C_DECL

typedef enum faest_paramid_t {
  PARAMETER_SET_INVALID   = 0,
  FAEST_128S              = 1,
  FAEST_128F              = 2,
  FAEST_192S              = 3,
  FAEST_192F              = 4,
  FAEST_256S              = 5,
  FAEST_256F              = 6,
  FAEST_EM_128S           = 7,
  FAEST_EM_128F           = 8,
  FAEST_EM_192S           = 9,
  FAEST_EM_192F           = 10,
  FAEST_EM_256S           = 11,
  FAEST_EM_256F           = 12,
  PARAMETER_SET_MAX_INDEX = 13
} faest_paramid_t;

typedef struct faest_paramset_t {
  // main parameters
  uint16_t lambda;
  uint8_t tau;
  uint8_t w_grind;
  uint16_t T_open;
  uint16_t l;

  // extra parameters
  uint16_t k;
  uint8_t tau0;
  uint8_t tau1;
  uint32_t L;

  // OWF parameters
  uint16_t Nst;
  uint16_t Ske;
  uint16_t R;
  uint16_t Senc;
  uint16_t Lke;
  uint16_t Lenc;

  // additional parameters
  uint16_t sig_size;
  uint8_t owf_input_size;
  uint8_t owf_output_size;
} faest_paramset_t;

const char* ATTR_CONST faest_get_param_name(faest_paramid_t paramid);
const faest_paramset_t* ATTR_CONST faest_get_paramset(faest_paramid_t paramid);

static inline bool ATTR_PURE faest_is_em(const faest_paramset_t* params) {
  // EM instances do not have key expansion constraints
  return params->Ske == 0;
}

FAEST_END_C_DECL

#endif
