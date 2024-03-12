/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_VOLE_H
#define FAEST_VOLE_H

#include "vc.h"
#include <stdbool.h>

FAEST_BEGIN_C_DECL

// k_b is at most 12, so chalout needs to point to an array of at most 12 bytes
int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout);

FAEST_END_C_DECL

#endif
