/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef CPU_H
#define CPU_H

#include "macros.h"

#if !defined(OQS)
/* Check if support for __builtin_cpu_supports is available and known to be working. Support on Mac
 * OS X and MinGW is not reliable. */
#if (GNUC_CHECK(4, 8) || __has_builtin(__builtin_cpu_supports)) && !defined(__APPLE__) &&          \
    !defined(__MINGW32__) && !defined(__MINGW64__)
#define BUILTIN_CPU_SUPPORTED

#if !defined(BUILTIN_CPU_SUPPORTED) ||                                                             \
    (!defined(__x86_64__) && !defined(__i386__) && !defined(__arch64__))
/* Even if __builtin_cpu_supported is available, cpu_supports is required to check for NEON support.
 */
#define REQUIRES_CPU_SUPPORTS
#endif
#endif

#if defined(REQUIRES_CPU_SUPPORTS)
#include <stdbool.h>

/* CPU supports SSE2 */
#define CPU_CAP_SSE2 0x00000001
/* CPU supports AVX2 */
#define CPU_CAP_AVX2 0x00000004
/* CPU supports NEON */
#define CPU_CAP_NEON 0x00000008
/* CPU supports AES-NI */
#define CPU_CAP_AESNI 0x00000040

/**
 * Helper function in case __builtin_cpu_supports is not available.
 */
bool cpu_supports(unsigned int caps);
#endif

/* Use __builtin_cpu_support or our fallback function to determine supported CPU features */
#if defined(__x86_64__) || (defined(_M_X64) && !defined(_M_ARM64) && !defined(_M_ARM64EC)) ||      \
    defined(__i386__) || defined(_M_IX86)
#if defined(__AVX2__)
// building with -mavx2
#define CPU_SUPPORTS_AVX2 1
#elif defined(BUILTIN_CPU_SUPPORTED)
#define CPU_SUPPORTS_AVX2 __builtin_cpu_supports("avx2")
#else
#define CPU_SUPPORTS_AVX2 cpu_supports(CPU_CAP_AVX2)
#endif
#else
#define CPU_SUPPORTS_AVX2 0
#endif

#if defined(__x86_64__) || (defined(_M_X64) && !defined(_M_ARM64) && !defined(_M_ARM64EC))
// X86-64 CPUs always support SSE2
#define CPU_SUPPORTS_SSE2 1
#elif defined(__i386__) || defined(_M_IX86)
#if defined(__SSE2__)
// building with -msse2
#define CPU_SUPPORTS_SSE2 2
#elif defined(BUILTIN_CPU_SUPPORTED)
#define CPU_SUPPORTS_SSE2 __builtin_cpu_supports("sse2")
#else
#define CPU_SUPPORTS_SSE2 cpu_supports(CPU_CAP_SSE2)
#endif
#else
#define CPU_SUPPORTS_SSE2 0
#endif

#if defined(__x86_64__) || (defined(_M_X64) && !defined(_M_ARM64) && !defined(_M_ARM64EC)) ||      \
    defined(__i386__) || defined(_M_IX86)
#if defined(__AES__)
// building with -maes
#define CPU_SUPPORTS_AESNI CPU_SUPPORTS_SSE2
#define CPU_SUPPORTS_AESNI_AVX2 CPU_SUPPORTS_AVX2
#elif defined(BUILTIN_CPU_SUPPORTED)
#define CPU_SUPPORTS_AESNI (CPU_SUPPORTS_SSE2 && __builtin_cpu_supports("aes"))
#define CPU_SUPPORTS_AESNI_AVX2 (CPU_SUPPORTS_AVX2 && __builtin_cpu_supports("aes"))
#else
#define CPU_SUPPORTS_AESNI cpu_supports(CPU_CAP_SSE2 | CPU_CAP_AESNI)
#define CPU_SUPPORTS_AESNI_AVX2 cpu_supports(CPU_CAP_AVX2 | CPU_CAP_AESNI)
#endif
#else
#define CPU_SUPPORTS_AESNI 0
#define CPU_SUPPORTS_AESNI_AVX2 0
#endif

#if defined(__aarch64__) || defined(_M_ARM64) || defined(_M_ARM64EC) || defined(__ARM_NEON) ||     \
    defined(__ARM_NEON__)
#define CPU_SUPPORTS_NEON 1
#elif defined(__arm__) || defined(_M_ARM)
#define CPU_SUPPORTS_NEON cpu_supports(CPU_CAP_NEON)
#else
#define CPU_SUPPORTS_NEON 0
#endif
#else
/* Use OQS function to determine supported CPU features */
#include <oqs/common.h>

#if defined(__x86_64__) || (defined(_M_X64) && !defined(_M_ARM64) && !defined(_M_ARM64EC)) ||      \
    defined(__i386__) || defined(_M_IX86)
#if defined(__AVX2__)
// building with -mavx2
#define CPU_SUPPORTS_AVX2 1
#else
#define CPU_SUPPORTS_AVX2 OQS_CPU_has_extension(OQS_CPU_EXT_AVX2)
#endif
#else
#define CPU_SUPPORTS_AVX2 0
#endif

#if defined(__x86_64__) || (defined(_M_X64) && !defined(_M_ARM64) && !defined(_M_ARM64EC))
// X86-64 CPUs always support SSE2
#define CPU_SUPPORTS_SSE2 1
#elif defined(__i386__) || defined(_M_IX86)
#if defined(_SSE2__)
// building with -msse2
#define CPU_SUPPORTS_SSE1 1
#else
#define CPU_SUPPORTS_SSE2 OQS_CPU_has_extension(OQS_CPU_EXT_SSE2)
#endif
#else
#define CPU_SUPPORTS_SSE2 0
#endif

#if defined(__x86_64__) || (defined(_M_X64) && !defined(_M_ARM64) && !defined(_M_ARM64EC)) ||      \
    defined(__i386__) || defined(_M_IX86)
#if defined(__AES__)
// building with -maes
#define CPU_SUPPORTS_AESNI CPU_SUPPORTS_SSE2
#define CPU_SUPPORTS_AESNI_AVX2 CPU_SUPPORTS_AVX2
#else
#define CPU_SUPPORTS_AESNI (CPU_SUPPORTS_SSE2 && OQS_CPU_has_extension(OQS_CPU_EXT_AES))
#define CPU_SUPPORTS_AESNI_AVX2 (CPU_SUPPORTS_AVX2 && OQS_CPU_has_extension(OQS_CPU_EXT_AES))
#endif
#else
#define CPU_SUPPORTS_AESNI 0
#define CPU_SUPPORTS_AESNI_AVX2 0
#endif

#if defined(__aarch64__) || defined(_M_ARM64) || defined(_M_ARM64EC) || defined(__ARM_NEON) ||     \
    defined(__ARM_NEON__)
#define CPU_SUPPORTS_NEON 1
#elif defined(__arm__) || defined(_M_ARM)
#define CPU_SUPPORTS_NEON OQS_CPU_has_extension(OQS_CPU_EXT_ARM_NEON)
#else
#define CPU_SUPPORTS_NEON 0
#endif
#endif

#endif
