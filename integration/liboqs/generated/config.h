/*
 *  SPDX-License-Identifier: MIT
 *
 *  Generated for liboqs integration.
 */

#ifndef CONFIG_H
#define CONFIG_H

#ifndef OQS
#define OQS
#endif

#include <oqs/aes.h>
#include <oqs/common.h>
#include <oqs/rand.h>
#include <oqs/sha3.h>

#ifndef FAEST_LIBOQS_BUILD
#define FAEST_LIBOQS_BUILD
#endif

#if defined(__APPLE__) || (defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L)
#define HAVE_POSIX_MEMALIGN
#endif

#if defined(__linux__)
#define HAVE_MEMALIGN
#endif

#endif
