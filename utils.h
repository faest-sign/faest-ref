#include <stdlib.h>
#include <stdint.h>

#include "compat.h"

void printHex(const char* s, const uint8_t* data, size_t len);

void xorUint8Arr(const uint8_t* a, const uint8_t* b, uint8_t* out, uint32_t len);

int printUint8Arr(uint8_t* arr, uint32_t unitSize, uint32_t len);

uint8_t* getUint8ArrPtr(uint8_t* in, uint32_t unitSize, uint32_t idx);