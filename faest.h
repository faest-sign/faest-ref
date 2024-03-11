/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_H
#define FAEST_H

#include <stdint.h>
#include <stddef.h>

#include "instances.h"

void faest_sign(uint8_t* sig, const uint8_t* msg, size_t msglen, const uint8_t* owf_key,
                const uint8_t* owf_input, const uint8_t* owf_output, const uint8_t* rho,
                size_t rholen, const faest_paramset_t* params);

int faest_verify(const uint8_t* msg, size_t msglen, const uint8_t* sig, const uint8_t* owf_input,
                 const uint8_t* owf_output, const faest_paramset_t* params);

ATTR_PURE const uint8_t* dsignature_iv(const uint8_t* base_ptr, const faest_paramset_t* params);
ATTR_PURE const uint8_t* dsignature_chall_3(const uint8_t* base_ptr,
                                            const faest_paramset_t* params);
ATTR_PURE const uint8_t* dsignature_pdec(const uint8_t* base_ptr, unsigned int index,
                                         const faest_paramset_t* params);
ATTR_PURE const uint8_t* dsignature_com(const uint8_t* base_ptr, unsigned int index,
                                        const faest_paramset_t* params);
ATTR_PURE const uint8_t* dsignature_u_tilde(const uint8_t* base_ptr,
                                            const faest_paramset_t* params);
ATTR_PURE const uint8_t* dsignature_c(const uint8_t* base_ptr, unsigned int index,
                                      const faest_paramset_t* params);
ATTR_PURE const uint8_t* dsignature_d(const uint8_t* base_ptr, const faest_paramset_t* params);

#endif
