//
//  rng.h
//
//  Created by Bassham, Lawrence E (Fed) on 8/29/17.
//  Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
//

#ifndef rng_h
#define rng_h

#include "../macros.h"

FAEST_BEGIN_C_DECL

#define RNG_SUCCESS 0
#define RNG_BAD_MAXLEN -1
#define RNG_BAD_OUTBUF -2
#define RNG_BAD_REQ_LEN -3

typedef struct {
  unsigned char buffer[16];
  int buffer_pos;
  unsigned long length_remaining;
  unsigned char key[32];
  unsigned char ctr[16];
} AES_XOF_struct;

void AES256_CTR_DRBG_Update(unsigned char* provided_data, unsigned char* Key, unsigned char* V);

int seedexpander_init(AES_XOF_struct* ctx, unsigned char* seed, unsigned char* diversifier,
                      unsigned long maxlen);

int seedexpander(AES_XOF_struct* ctx, unsigned char* x, unsigned long xlen);

void randombytes_init(unsigned char* entropy_input, unsigned char* personalization_string,
                      int security_strength);

int randombytes(unsigned char* x, unsigned long long xlen);

FAEST_END_C_DECL

#endif /* rng_h */
