/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_DEFINES_H
#define FAEST_DEFINES_H

#if !defined(FAEST_EXPORT)
#if !defined(FAEST_STATIC) && (defined(_WIN16) || defined(_WIN32) || defined(_WIN64))
#define FAEST_EXPORT __declspec(dllimport)
#else
#define FAEST_EXPORT
#endif
#endif

#if defined(_WIN16) || defined(_WIN32)
#define FAEST_CALLING_CONVENTION __stdcall
#else
#define FAEST_CALLING_CONVENTION
#endif

#endif
