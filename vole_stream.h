#ifndef FAEST_VOLE_STREAM_H
#define FAEST_VOLE_STREAM_H

#include "vc_stream.h"
#include <stdbool.h>

FAEST_BEGIN_C_DECL

void partial_vole_commit_cmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                             unsigned int start, unsigned int len,
                             uint8_t* v, uint8_t* u, uint8_t* hcom, uint8_t* c,
                             const faest_paramset_t* params);

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
