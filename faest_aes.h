#ifndef FAEST_FAEST_AES_H
#define FAEST_FAEST_AES_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "instances.h"

void aes_key_schedule_forward(uint32_t lambda, uint32_t R, uint32_t Nwd, uint32_t Lke, uint32_t m,
                              const uint8_t* x, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                              uint8_t* out);

uint8_t* aes_key_schedule_backward(uint32_t lambda, uint32_t R, uint32_t Nwd, uint32_t Ske,
                                   uint8_t Lke, uint32_t m, const uint8_t* x, const uint8_t* xk,
                                   uint8_t Mtag, uint8_t Mkey, const uint8_t* delta);

int aes_key_schedule_constraints(uint32_t lambda, uint32_t R, uint32_t Nwd, uint32_t Ske,
                                 uint32_t Lke, const uint8_t* w, const uint8_t* v,
                                 const uint8_t Mkey, const uint8_t* q, const uint8_t* delta,
                                 uint8_t* A0, uint8_t* A1, uint8_t* k, uint8_t* vk, uint8_t* B,
                                 uint8_t* qk);

int aes_enc_forward(uint32_t lambda, uint32_t R, uint32_t m, uint32_t Lenc, const uint8_t* x,
                    const uint8_t* xk, const uint8_t* in, uint8_t Mtag, uint8_t Mkey,
                    const uint8_t* delta, uint8_t* y_out);

int aes_enc_backward(uint32_t lambda, uint32_t R, uint32_t m, uint32_t Lenc, const uint8_t* x,
                     const uint8_t* xk, uint8_t Mtag, uint8_t Mkey, const uint8_t* delta,
                     const uint8_t* out, uint8_t* y_out);

void aes_enc_constraints(uint32_t lambda, uint32_t R, uint32_t Lenc, uint32_t Senc,
                         const uint8_t* in, const uint8_t* out, const uint8_t* w, const uint8_t* v,
                         const uint8_t* k, const uint8_t* vk, uint8_t Mkey, const uint8_t* q,
                         const uint8_t* qk, const uint8_t* delta, uint8_t* A0, uint8_t* A1,
                         uint8_t* B);

void aes_prove(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
               const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params);

uint8_t* aes_verify(uint8_t* d, uint8_t** Q, const uint8_t* chall_2, const uint8_t* chall_3,
                    const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out,
                    const faest_paramset_t* params);

#endif
