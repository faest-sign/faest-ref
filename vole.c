#include "vole.h"

void ConvertToVole(const faest_paramset_t* params, const vec_com_t* vecCom, uint32_t depth,
                   uint32_t outLen, uint8_t* u, uint8_t* v) {
  uint8_t* r;
  r = malloc(getBinaryTreeNodeCount(depth) * outLen);

  memset(r, 0, outLen);
}