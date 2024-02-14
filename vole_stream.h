#ifndef FAEST_VOLE_STREAM_H
#define FAEST_VOLE_STREAM_H

#include "vc_stream.h"
#include <stdbool.h>
/*
 *  SPDX-License-Identifier: MIT
 */

FAEST_BEGIN_C_DECL

// k_b is at most 12, so chalout needs to point to an array of at most 12 bytes
/*
int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout);
*/

void stream_vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                 const faest_paramset_t* params, uint8_t* hcom, stream_vec_com_t* sVecCom, uint8_t* c,
                 uint8_t* u, uint8_t** v);

FAEST_END_C_DECL

#endif
