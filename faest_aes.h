#include <stdbool.h>

#include "vc.h"
#include "vole.h"
#include "fields.h"
#include "universal_hashing.h"

void aes_extend_witness(uint32_t lambda, const uint8_t* sk, const uint8_t* in, uint8_t* w);

void aes_prove(uint8_t* w, uint8_t* u, uint8_t** v, uint8_t* in, uint8_t* out, uint8_t* chal_2,
               uint32_t lambda, uint32_t tau, uint32_t l, uint8_t* a_tilde, uint8_t* b_tilde);