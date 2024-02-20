#ifndef FAEST_VOLE_STREAM_H
#define FAEST_VOLE_STREAM_H

#include "vc_stream.h"
#include <stdbool.h>

FAEST_BEGIN_C_DECL

void partial_vole_commit_cmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                             const faest_paramset_t* params, stream_vec_com_t* sVecCom, uint8_t* v,
                             unsigned int start, unsigned int len);

void vole_commit_u_hcom_c(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                          const faest_paramset_t* params, uint8_t* hcom, stream_vec_com_t* sVecCom,
                          uint8_t* c, uint8_t* u);

void partial_vole_commit_rmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int start,
                             unsigned int len, const faest_paramset_t* params,
                             stream_vec_com_t* sVecCom, uint8_t* v);

FAEST_END_C_DECL

#endif
