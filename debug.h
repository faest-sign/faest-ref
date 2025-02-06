#ifndef DEBUG_H
#define DEBUG_H

#include <stdio.h>
#include "random_oracle.h"
#include "fields.h"

inline void debug_print_buf(const char* name, const void* buf, size_t n) {
    const unsigned char* data = buf;
    printf("%s = ", name);
    for (size_t i = 0; i < n; ++i) {
        printf("%02x", (unsigned) data[i]);
    }
    printf("\n");
}

inline void debug_print_bf128(const char* name, const bf128_t* fe) {
    size_t n = sizeof(bf128_t);
    const unsigned char* data = (unsigned char*)fe;
    printf("%s = 0x", name);
    for (size_t i = 1; i <= n; ++i) {
        printf("%02x", data[n-i]);
    }
    printf("\n");
}

inline void debug_print_buf_bits(const char* name, const void* buf, size_t n) {
    const unsigned char* data = buf;
    printf("%s = ", name);
    for (size_t i = 0; i < n; ++i) {
        printf("%x", (unsigned) data[i]);
    }
    printf("\n");
}

inline void H_intermediate(const char* name, const hash_context* ctx) {

  uint8_t buf[16];
  hash_context new_ctx;
  memcpy(&new_ctx, ctx, sizeof(new_ctx));
  hash_final(&new_ctx);
  hash_squeeze(&new_ctx, buf, sizeof(buf));
  debug_print_buf(name, buf, sizeof(buf));
}


#endif
