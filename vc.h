#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include "fields.h"

typedef uint8_t com_t;           // Using Keccak 128,256
typedef uint8_t dk_t;         // Using AES in CTR mode, with output (2*lambda)
typedef uint8_t dcom_t;
typedef uint8_t sd_t;         // Size lambda for each message

typedef struct vec_com_t {
    com_t* com;
    dk_t** dk;
    dcom_t** dcom;
    sd_t** sd;
} vec_com_t;

void vector_commitment(const uint8_t* rootKey, const faest_param_t* params, const uint8_t* salt, vec_com_t* vecCom);