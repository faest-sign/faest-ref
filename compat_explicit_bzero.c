/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "compat.h"

#if !defined(HAVE_EXPLICIT_BZERO) && !defined(OQS)
#if defined(_WIN32)
#include <windows.h>
#endif

void faest_explicit_bzero(void* a, size_t len) {
#if defined(_WIN32)
  SecureZeroMemory(a, len);
#else
  volatile char* p = a;
  for (; len; ++p, --len) {
    *p = 0;
  }
#endif
}
#endif /* HAVE_EXPLICIT_BZERO */
