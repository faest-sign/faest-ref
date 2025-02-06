#ifndef DEBUG_H
#define DEBUG_H

#include <assert.h>
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
    assert(sizeof(bf192_t) >= 16);
    const size_t n = 16;
    const unsigned char* data = (unsigned char*)fe;
    printf("%s = 0x", name);
    for (size_t i = 1; i <= n; ++i) {
        printf("%02x", data[n-i]);
    }
    printf("\n");
}

inline void debug_print_bf192(const char* name, const bf192_t* fe) {
    assert(sizeof(bf192_t) >= 24);
    const size_t n = 24;
    const unsigned char* data = (unsigned char*)fe;
    printf("%s = 0x", name);
    for (size_t i = 1; i <= n; ++i) {
        printf("%02x", data[n-i]);
    }
    printf("\n");
}

inline void debug_print_bf256(const char* name, const bf256_t* fe) {
    assert(sizeof(bf192_t) >= 32);
    const size_t n = 32;
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
