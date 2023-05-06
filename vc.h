#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include "fields.h"

uint32_t treeDepth;

typedef uint8_t com_t;           // Using Keccak 128,256
typedef uint8_t decom_t;         // Using AES in CTR mode, with output (2*lambda)
typedef uint8_t mvec_t;      // Size lambda for each message

typedef uint8_t key128_t[16];

typedef struct vec_com_t {
    com_t *com;
    decom_t *decom;
    mvec_t *messages;
} vec_com_t;

vec_com_t vector_commitment(uint64_t treeDepth, uint32_t seclvl, key128_t key, bf8_t* iv);