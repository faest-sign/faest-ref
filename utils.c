#include <stdio.h>

#include "utils.h"

void printHex(const char* s, const uint8_t* data, size_t len) {
  printf("%s: ", s);
  for (size_t i = 0; i < len; i++) {
    printf("%02X", data[i]);
  }
  printf("\n");
}