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

void partial_vole_commit_cmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                             unsigned int start, unsigned int end, uint8_t* v, uint8_t* u,
                             uint8_t* hcom, uint8_t* c, const faest_paramset_t* params);

void partial_vole_commit_rmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int start,
                             unsigned int len, const faest_paramset_t* params, uint8_t* v);

void partial_vole_reconstruct_cmo(const uint8_t* iv, const uint8_t* chall,
                                  const uint8_t* const* pdec, const uint8_t* const* com_j,
                                  uint8_t* hcom, uint8_t* q, unsigned int ellhat,
                                  const faest_paramset_t* params, unsigned int start,
                                  unsigned int len);
void vole_reconstruct_hcom(const uint8_t* iv, const uint8_t* chall, const uint8_t* const* pdec,
                           const uint8_t* const* com_j, uint8_t* hcom, unsigned int ellhat,
                           const faest_paramset_t* params);
void partial_vole_reconstruct_rmo(const uint8_t* iv, const uint8_t* chall,
                                  const uint8_t* const* pdec, const uint8_t* const* com_j,
                                  uint8_t* q, unsigned int ellhat, const faest_paramset_t* params,
                                  unsigned int start, unsigned int len);

FAEST_END_C_DECL

#endif
