/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#include "macros.h"

/* define HAVE_* for more known good configurations */
#if !defined(HAVE_POSIX_MEMALIGN) &&                                                               \
    ((defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L) || defined(__APPLE__))
/* defined in POSIX and available on OS X */
#define HAVE_POSIX_MEMALIGN
#endif /* HAVE_POSIX_MEMALIGN */

#if !defined(HAVE_MEMALIGN) && defined(__linux__)
/* always available on Linux */
#define HAVE_MEMALIGN
#endif /* HAVE_MEMALIGN */
#endif /* HAVE_CONFIG_H */

#include "compat.h"

#if !defined(HAVE_ALIGNED_ALLOC) && !defined(OQS)
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#if !defined(HAVE_POSIX_MEMALIGN) || defined(__MING32__) || defined(__MING64__) || defined(_MSC_VER)
#include <malloc.h>
#endif

#if defined(HAVE_POSIX_MEMALIGN)
static_assert((sizeof(void*) & (sizeof(void*) - 1)) == 0, "sizeof void* is not a power of 2");
#endif

/* The fallback implementation tries to be as generic as possible. While all callers in this code
 * base satisfy all requirements, we still check them. Thereby, the fallback implementation may be
 * of use for others as well. */
void* faest_aligned_alloc(size_t alignment, size_t size) {
  /* check alignment (power of 2) and size (multiple of alignment) */
  if (alignment & (alignment - 1) || size & (alignment - 1)) {
    errno = EINVAL;
    return NULL;
  }

#if defined(HAVE_POSIX_MEMALIGN)
  /* check alignment (needs to be multiple of sizeof(void*)); posix_memalign has this additional
   * requirement */
  if (alignment < sizeof(void*)) {
    alignment = sizeof(void*);
    /* fix up size; needs to be a multiple of alignment, i.e., sizeof(void*) */
    size = (size + sizeof(void*) - 1) & ~(sizeof(void*) - 1);
  }

  void* ptr     = NULL;
  const int err = posix_memalign(&ptr, alignment, size);
  if (err) {
    errno = err;
  }
  return ptr;
#elif defined(HAVE_MEMALIGN)
  return memalign(alignment, size);
#elif defined(__MINGW32__) || defined(__MINGW64__)
  return __mingw_aligned_malloc(size, alignment);
#elif defined(_MSC_VER)
  return _aligned_malloc(size, alignment);
#else
  /* generic fallback implementation: malloc memory and return aligned address; in the byte before
   * the aligned location, place the offset to the original pointer */
  if (!size) {
    return NULL;
  }
  const size_t offset = alignment - 1 + sizeof(uint8_t);
  uint8_t* buffer     = malloc(size + offset);
  if (!buffer) {
    return NULL;
  }

  uint8_t* ptr   = (uint8_t*)(((uintptr_t)(buffer) + offset) & ~(alignment - 1));
  ptrdiff_t diff = ptr - buffer;
  if (diff > UINT8_MAX) {
    /* this should never happen in our code, but just to be safe */
    free(buffer);
    errno = EINVAL;
    return NULL;
  }
  ptr[-1] = diff;
  return ptr;
#endif
}

void faest_aligned_free(void* ptr) {
#if defined(HAVE_POSIX_MEMALIGN) || defined(HAVE_MEMALIGN)
  free(ptr);
#elif defined(__MINGW32__) || defined(__MINGW64__)
  __mingw_aligned_free(ptr);
#elif defined(_MSC_VER)
  _aligned_free(ptr);
#else
  if (ptr) {
    uint8_t* u8ptr = ptr;
    free(u8ptr - u8ptr[-1]);
  }
#endif
}
#endif /* HAVE_ALIGNED_ALLOC */