/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_VOLE_H
#define FAEST_VOLE_H

#include <stdbool.h>

#include "bavc.h"

FAEST_BEGIN_C_DECL

// void vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
//                  const faest_paramset_t* params, bavc_t* vecCom, uint8_t* c, uint8_t* u,
//                  uint8_t** v);

// TODO: For now putting this here
static inline uint8_t get_bit_from_pt(const uint8_t *p, size_t i) {
    return (uint8_t)((p[i / 8] >> (i % 8)) & 1u);
}
static inline void set_bit_to_pt(uint8_t *p, size_t i, uint8_t v) {
    size_t byte = i / 8;
    uint8_t bit = (uint8_t)(1u << (i % 8));
    p[byte] = (uint8_t)((p[byte] & ~bit) | (bit & (uint8_t)(0u - v)));
}

void vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                 const faest_paramset_t* params, bavc_t* bavc, uint8_t* c, uint8_t* c_mult,
                 uint8_t* u, uint8_t** v, uint8_t* u_dash_m, uint8_t* v_dash);

// bool vole_reconstruct(uint8_t* com, uint8_t** q, const uint8_t* iv, const uint8_t* chall_3,
//                       const uint8_t* decom_i, const uint8_t* c, unsigned int ellhat,
//                       const faest_paramset_t* params);

bool vole_reconstruct(uint8_t* com, uint8_t** q, const uint8_t* iv, const uint8_t* chall_3,
                      const uint8_t* decom_i, const uint8_t* c, uint8_t* c_mult,
                      uint8_t* q_bar, unsigned int ellhat, const faest_paramset_t* params);

#if defined(FAEST_TESTS)
unsigned int convert_to_vole(const uint8_t* iv, const uint8_t* sd, bool sd0_bot, unsigned int i,
                             unsigned int outLenBytes, uint8_t* u, uint8_t* v,
                             const faest_paramset_t* params);
#endif

FAEST_END_C_DECL

#endif
