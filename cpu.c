/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#if !defined(HAVE_SYS_AUXV_H) && defined(__linux__)
#define HAVE_SYS_AUXV_H
#endif

#if !defined(HAVE_ASM_HWCAP_H) && defined(__linux__) && (defined(__arm__) || defined(__aarch64__))
#define HAVE_ASM_HWCAP_H
#endif
#endif

#include "cpu.h"

#if defined(REQUIRES_CPU_SUPPORTS)
#if (defined(__arm__) || defined(__aarch64__)) && defined(HAVE_SYS_AUXV_H) &&                      \
    defined(HAVE_ASM_HWCAP_H)
#include <asm/hwcap.h>
#include <sys/auxv.h>

static unsigned int init_caps(void) {
  unsigned int caps         = 0;
  const unsigned long hwcap = getauxval(AT_HWCAP);
#if defined(__aarch64__)
  if (hwcap & HWCAP_ASIMD) {
    caps |= CPU_CAP_NEON;
  }
  if (hwcap & HWCAP_AES) {
    caps |= CPU_CAP_ARM_AES;
  }
#else
  if (hwcap & HWCAP_NEON) {
    caps |= CPU_CAP_NEON;
  }
#if defined(HWCAP2_AES)
  if ((hwcap & HWCAP_NEON) && (getauxval(AT_HWCAP2) & HWCAP2_AES)) {
    caps |= CPU_CAP_ARM_AES;
  }
#endif
#endif
  return caps;
}

#elif defined(_M_ARM64) || defined(_M_ARM64EC) || defined(_M_ARM)
#include <windows.h>

static unsigned int init_caps(void) {
  unsigned int caps = 0;
  if (IsProcessorFeaturePresent(PF_ARM_NEON_INSTRUCTIONS_AVAILABLE)) {
    caps |= CPU_CAP_NEON;
  }
  if (IsProcessorFeaturePresent(PF_ARM_V8_CRYPTO_INSTRUCTIONS_AVAILABLE)) {
    caps |= CPU_CAP_ARM_AES;
  }
  return caps;
}

#elif (defined(__x86_64__) || defined(__i386__) || defined(_M_IX86) ||                             \
       (defined(_M_X64) && !defined(_M_ARM64) && !defined(_M_ARM64EC))) &&                         \
    (defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER))

#ifdef _MSC_VER
#include <intrin.h>

static unsigned init_caps(void) {
  unsigned int caps = 0;

  union {
    struct {
      unsigned int eax, ebx, ecx, edx;
    };
    int data[4];
  } regs = {0};

  __cpuid(regs.data, 0);
  unsigned int max = regs.eax;

  if (max >= 1) {
    __cpuid(regs.data, 0);
    if (regs.edx & (1 << 26)) {
      caps |= CPU_CAP_SSE2;
    }
    if (regs.ecx & (1 << 25)) {
      caps |= CPU_CAP_AESNI;
    }
  }

  if (max >= 7) {
    __cpuidex(regs.data, 7, 0);
    if (regs.ebx & (1 << 5)) {
      caps |= CPU_CAP_AVX2;
    }
  }

  return caps;
}
#else
#include <cpuid.h>

static unsigned init_caps(void) {
  unsigned int caps = 0;
  unsigned int eax, ebx, ecx, edx;

  if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
    if (edx & (1 << 26)) {
      caps |= CPU_CAP_SSE2;
    }
    if (ecx & (1 << 25)) {
      caps |= CPU_CAP_AESNI;
    }
  }

  if (__get_cpuid(7, &eax, &ebx, &ecx, &edx)) {
    if (ebx & (1 << 5)) {
      caps |= CPU_CAP_AVX2;
    }
  }

  return caps;
}
#endif

#else

static unsigned init_caps(void) {
  return 0;
}

#endif

#include <limits.h>

static unsigned int cpu_caps = UINT_MAX;

bool cpu_supports(unsigned int caps) {
  if (cpu_caps == UINT_MAX) {
    cpu_caps = init_caps();
  }

  return (cpu_caps & caps) == caps;
}
#endif
