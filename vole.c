/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "vole.h"
#include "aes.h"
#include "utils.h"
#include "random_oracle.h"

#include <stdbool.h>
#include <string.h>

int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout) {
  if (i >= t0 + t1) {
    return 0;
  }

  unsigned int lo;
  unsigned int hi;
  if (i < t0) {
    lo = i * k0;
    hi = ((i + 1) * k0);
  } else {
    unsigned int t = i - t0;
    lo             = (t0 * k0) + (t * k1);
    hi             = (t0 * k0) + ((t + 1) * k1);
  }

  assert(hi - lo == k0 || hi - lo == k1);
  for (unsigned int j = lo; j < hi; ++j) {
    // set_bit(chalout, i - lo, get_bit(chal, i));
    chalout[j - lo] = ptr_get_bit(chal, j);
  }
  return 1;
}

