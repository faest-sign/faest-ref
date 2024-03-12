/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vc.h"
#include "random_oracle.h"
#include "compat.h"
#include "aes.h"
#include "instances.h"

#include <assert.h>
#include <string.h>

unsigned int NumRec(unsigned int depth, const uint8_t* bi) {
  unsigned int out = 0;
  for (unsigned int i = 0; i < depth; i++) {
    out += ((unsigned int)bi[i]) << i;
  }
  return out;
}
