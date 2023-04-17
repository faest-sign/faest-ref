// AES CTR mode 128/192/256

#include <stdint.h>

typedef uint8_t bf8_t;
typedef uint64_t bf64_t;

typedef bf8_t state_t[4][4];

bf8_t round_key[16*11];
bf8_t iv[16];

void initialize(bf8_t* key, bf8_t* iv_);
void encrypt(bf8_t* buffer);
