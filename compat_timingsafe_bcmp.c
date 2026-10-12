/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#include "macros.h"

/* define HAVE_* for more known good configurations */
#if !defined(HAVE_CONSTTIME_MEMEQUAL) && NETBSD_CHECK(7, 0)
/* consttime_memequal was introduced in NetBSD 7.0 */
#define HAVE_CONSTTIME_MEMEQUAL
#endif /* HAVE_CONSTTIME_MEMEQUAL */
#endif /* HAVE_CONFIG_H */

#include "compat.h"

#if !defined(HAVE_TIMINGSAFE_BCMP) && !defined(OQS)
int faest_timingsafe_bcmp(const void* a, const void* b, size_t len) {
#if defined(HAVE_CONSTTIME_MEMEQUAL)
  return !consttime_memequal(a, b, len);
#else
  const unsigned char* p1 = a;
  const unsigned char* p2 = b;

  unsigned int res = 0;
  for (; len; --len, ++p1, ++p2) {
    res |= *p1 ^ *p2;
  }
  return res;
#endif
}
#endif /* HAVE_TIMINGSAFE_BCMP */