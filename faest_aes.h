#include <stdbool.h>

#include "vc.h"
#include "vole.h"
#include "fields.h"
#include "universal_hashing.h"

void aes_extend_witness(uint32_t lambda, const uint8_t* sk, const uint8_t* in, uint8_t* w);

int aes_key_schedule_forward(uint32_t lambda, uint32_t m, const uint8_t* x, uint8_t Mtag,
                             uint8_t Mkey, const uint8_t* delta, uint8_t* y);

int aes_key_schedule_backward(uint32_t lambda, uint32_t m, const uint8_t* x, const uint8_t* xk,
                              uint8_t Mtag, uint8_t Mkey, const uint8_t* delta, uint32_t S_ke_Byte,
                              uint8_t* y);

int aes_key_schedule_constraints(uint32_t lambda, const uint8_t* w, const uint8_t* v,
                                 const uint8_t Mkey, const uint8_t* q, const uint8_t* delta,
                                 uint8_t** A, uint8_t* k, uint8_t* vk, uint8_t* B, uint8_t* qk);

void aes_prove(uint8_t* w, uint8_t* u, uint8_t** v, uint8_t* in, uint8_t* out, uint8_t* chal_2,
               uint32_t lambda, uint32_t tau, uint32_t l, uint8_t* a_tilde, uint8_t* b_tilde);