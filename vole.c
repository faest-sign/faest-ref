#include "vole.h"

void ConvertToVole(const faest_paramset_t* params, const vec_com_t* vecCom, uint32_t depth,
                   uint32_t outLen, uint8_t* u, uint8_t* v) {
  uint8_t* r;
  r = malloc(getBinaryTreeNodeCount(depth) * outLen);

  uint8_t iv[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

  uint8_t* allZeros = malloc(outLen);
  memset(allZeros, 0, outLen);
  if (memcmp(vecCom->sd, allZeros, params->faest_param.seclvl) == 0) {
    memset(r, 0, outLen);
  } else {
    uint8_t* out = malloc(outLen);
    aes_prg(vecCom->sd, iv, out, params->faest_param.seclvl, outLen);
    memcpy(r, out, outLen);
  }

  uint64_t N = (1 << depth);
  for (uint32_t i = 1; i < N; i++) {
    uint8_t iv_[16] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
    uint8_t* out    = malloc(outLen);
    aes_prg(vecCom->sd + (vecCom->sd_uint_size * i), iv_, out, params->faest_param.seclvl, outLen);
    memcpy(r + (outLen * i), out, outLen);
  }

  memset(v, 0, depth * outLen);

  for (uint32_t j = 0; j < depth; j++) {
    for (uint32_t i = 0; i < N / (2 << (j + 1)); i++) {
      uint8_t tmp = *(v + (j * outLen)) ^ *(r + getNodeIndex(j, 2 * i + 1));
      memcpy(v + (j * outLen), &tmp, outLen);

      tmp = *(r + getNodeIndex(j, 2 * i)) ^ *(r + getNodeIndex(j, (2 * i) + 1));
      memcpy(r + getNodeIndex(j + 1, i), &tmp, outLen);
    }
  }
  memcpy(u, r + getNodeIndex(depth, 0), outLen);
}