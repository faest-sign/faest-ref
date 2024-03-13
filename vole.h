/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_VOLE_H
#define FAEST_VOLE_H

#include <stdbool.h>
#include <assert.h>

#include "vc.h"
#include "macros.h"

FAEST_BEGIN_C_DECL

typedef enum { EXCLUDE_U_HCOM_C, EXCLUDE_V, INCLUDE_ALL } vole_mode_of_operation_t;

typedef struct vole_mode_t {
  vole_mode_of_operation_t mode;
  uint8_t* v;
  uint8_t* u;
  uint8_t* hcom;
  uint8_t* c;
} vole_mode_t;

ATTR_CONST ATTR_ALWAYS_INLINE inline vole_mode_t vole_mode_all(uint8_t* v, uint8_t* u,
                                                               uint8_t* hcom, uint8_t* c) {
  assert(v != NULL);
  assert(u != NULL);
  assert(hcom != NULL);
  assert(c != NULL);
  return (vole_mode_t){.mode = INCLUDE_ALL, .v = v, .u = u, .hcom = hcom, .c = c};
}

ATTR_CONST ATTR_ALWAYS_INLINE inline vole_mode_t vole_mode_u_hcom_c(uint8_t* u, uint8_t* hcom,
                                                                    uint8_t* c) {
  assert(u != NULL);
  assert(hcom != NULL);
  assert(c != NULL);
  return (vole_mode_t){.mode = EXCLUDE_V, .v = NULL, .u = u, .hcom = hcom, .c = c};
}

ATTR_CONST ATTR_ALWAYS_INLINE inline vole_mode_t vole_mode_v(uint8_t* v) {
  assert(v != NULL);
  return (vole_mode_t){.mode = EXCLUDE_U_HCOM_C, .v = v, .u = NULL, .hcom = NULL, .c = NULL};
}

// k_b is at most 12, so chalout needs to point to an array of at most 12 bytes
int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout);

// Signer
void partial_vole_commit_cmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                             unsigned int chunk_start, unsigned int chunk_end,
                             vole_mode_t vole_mode, const faest_paramset_t* params);

void partial_vole_commit_rmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int start,
                             unsigned int len, const faest_paramset_t* params, uint8_t* v);

// Verifier
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
