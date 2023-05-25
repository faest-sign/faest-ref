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

void randombytes_init(unsigned char* entropy_input, unsigned char* personalization_string,
                      int security_strength);

int randombytes(unsigned char* x, unsigned long long xlen);

FAEST_END_C_DECL

#endif /* rng_h */
