#include <stdlib.h>
#include <stdint.h>

static int32_t nlz(uint32_t x);
uint32_t ceil_log2(uint32_t x);
void printHex(const char* s, const uint8_t* data, size_t len);