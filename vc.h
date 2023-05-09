#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include "tree.h"
#include "aes.h"
#include "random_oracle.h"
#include "instances.h"

typedef uint8_t com_t;
typedef uint8_t dk_t;
typedef uint8_t dcom_t;
typedef uint8_t sd_t;

typedef struct vec_com_t {
  com_t* h;
  dk_t* k;
  dcom_t* com;
  sd_t* sd;
} vec_com_t;

typedef struct vec_com_rec_t {
  com_t* h;
  dk_t* k;
  dcom_t* com;
  sd_t* m;
} vec_com_rec_t;

void NumRec(uint32_t depth, const uint8_t* bi, uint64_t* out);

void vector_commitment(const uint8_t* rootKey, const faest_paramset_t* params, vec_com_t* vecCom);

void vector_open(const faest_paramset_t* params, const vec_com_t* vecCom, const uint8_t* b,
                 uint8_t* pdec);