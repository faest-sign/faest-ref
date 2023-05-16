#include "vole.h"

void ConvertToVoleProver(uint32_t lambda, const uint8_t* sd, uint32_t numVoleInstances,
                         uint32_t depth, uint32_t outLenBytes, uint8_t* u, uint8_t* v) {

  uint8_t iv[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
  uint8_t* r     = malloc(getBinaryTreeNodeCount(numVoleInstances) * outLenBytes);

  /* Here r_d,0 is the root and r_0,0 is the first leaf, in this reference implementation we will
  keep it consistent like k_0,0 being the root and k_d,0 being the first leaf by manipulating it
  with getNodeIndex() due to ease of understanding... */

  uint8_t* allZeros = malloc(lambda / 8);
  memset(allZeros, 0, lambda / 8);
  if (memcmp(sd, allZeros, lambda / 8) == 0) {
    memset(r + (getNodeIndex(depth, 0) * outLenBytes), 0, outLenBytes);
  } else {
    uint8_t* out = malloc(outLenBytes);
    aes_prg(sd, iv, out, lambda, outLenBytes * 8);
    memcpy(r + (getNodeIndex(depth, 0) * outLenBytes), out, outLenBytes);
  }

  for (uint32_t i = 1; i < numVoleInstances; i++) {
    uint8_t iv_[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
    uint8_t* out    = malloc(outLenBytes);
    aes_prg(sd + ((lambda / 8) * i), iv_, out, lambda, outLenBytes * 8);
    memcpy(r + (outLenBytes * (getNodeIndex(depth, 0) + i)), out, outLenBytes);
  }

  memset(v, 0, depth * outLenBytes);
  for (uint32_t d = 0; d < depth; d++) {
    uint32_t depthloop = (numVoleInstances + ((1 << (d + 1)) - 1)) / (1 << (d + 1));
    for (uint32_t i = 0; i < depthloop; i++) {
      for (uint8_t b = 0; b < outLenBytes; b++) {
        *(v + ((depth - 1 - d) * outLenBytes) + b) =
            *(v + ((depth - 1 - d) * outLenBytes) + b) ^
            *(r + (getNodeIndex(depth - d, 2 * i + 1) * outLenBytes) + b);

        *(r + (getNodeIndex(depth - (d + 1), i) * outLenBytes) + b) =
            *(r + (getNodeIndex(depth - d, 2 * i) * outLenBytes) + b) ^
            *(r + (getNodeIndex(depth - d, (2 * i) + 1) * outLenBytes) + b);
      }
    }
  }
  memcpy(u, r, outLenBytes);
}

void ConvertToVoleVerifier(uint32_t lambda, const uint8_t* sd, uint32_t numVoleInstances,
                           uint32_t depth, uint32_t outLenBytes, uint8_t* u, uint8_t* v) {
  uint8_t iv[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
  uint8_t* r     = malloc(getBinaryTreeNodeCount(numVoleInstances) * outLenBytes);

  /* Here r_d,0 is the root and r_0,0 is the first leaf, in this reference implementation we will
  keep it consistent like k_0,0 being the root and k_d,0 being the first leaf by manipulating it
  with getNodeIndex() due to ease of understanding... */

  // TODO: stupid way to check,, change it !!
  uint8_t* allZeros = malloc(outLenBytes);
  memset(allZeros, 0, outLenBytes);
  if (memcmp(sd, allZeros, lambda) == 0) {
    memset(r + (getNodeIndex(depth, 0) * outLenBytes), 0, outLenBytes);
  } else {
    uint8_t* out = malloc(outLenBytes);
    aes_prg(sd, iv, out, lambda, outLenBytes * 8);
    memcpy(r + (getNodeIndex(depth, 0) * outLenBytes), out, outLenBytes);
  }

  for (uint32_t i = 1; i < numVoleInstances; i++) {
    uint8_t iv_[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
    uint8_t* out    = malloc(outLenBytes);
    aes_prg(sd + ((lambda / 8) * i), iv_, out, lambda, outLenBytes * 8);
    memcpy(r + (outLenBytes * (getNodeIndex(depth, 0) + i)), out, outLenBytes);
  }

  memset(v, 0, depth * outLenBytes);
  for (uint32_t d = 0; d < depth; d++) {
    uint32_t depthloop = (numVoleInstances + ((1 << (d + 1)) - 1)) / (1 << (d + 1));
    for (uint32_t i = 0; i < depthloop; i++) {
      for (uint8_t b = 0; b < outLenBytes; b++) {
        *(v + ((depth - 1 - d) * outLenBytes) + b) =
            *(v + ((depth - 1 - d) * outLenBytes) + b) ^
            *(r + (getNodeIndex(depth - d, 2 * i + 1) * outLenBytes) + b);

        *(r + (getNodeIndex(depth - (d + 1), i) * outLenBytes) + b) =
            *(r + (getNodeIndex(depth - d, 2 * i) * outLenBytes) + b) ^
            *(r + (getNodeIndex(depth - d, (2 * i) + 1) * outLenBytes) + b);
      }
    }
  }
}