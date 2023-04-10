/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../fields.h"
#include "../randomness.h"

#include <NTL/GF2X.h>

#define BOOST_TEST_MODULE fields_test
#include <boost/test/included/unit_test.hpp>

using namespace NTL;

namespace {
  GF2X compute_bf8_residue() {
    GF2X residue;
    SetCoeff(residue, 8);
    SetCoeff(residue, 4);
    SetCoeff(residue, 3);
    SetCoeff(residue, 1);
    SetCoeff(residue, 0);
    return residue;
  }

  GF2X compute_bf64_residue() {
    GF2X residue;
    SetCoeff(residue, 64);
    SetCoeff(residue, 4);
    SetCoeff(residue, 3);
    SetCoeff(residue, 1);
    SetCoeff(residue, 0);
    return residue;
  }

  const GF2X bf8_residue = compute_bf8_residue();
  const GF2X bf64_residue = compute_bf64_residue();

  GF2X from_bf(bf8_t v) {
    GF2X ret;
    for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
      SetCoeff(ret, i, v & 1);
    }
    return ret;
  }

  GF2X from_bf(bf64_t v) {
    GF2X ret;
    for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
      SetCoeff(ret, i, v & 1);
    }
    return ret;
  }

  void check_mul_bf8(bf8_t lhs, bf8_t rhs, bf8_t expected) {
    const auto result = bf8_mul(lhs, rhs);
    BOOST_TEST(result == expected);
    const auto ntl_result = MulMod(from_bf(lhs), from_bf(rhs), bf8_residue);
    BOOST_TEST(ntl_result == from_bf(result));
  }

  void check_mul_bf64(bf64_t lhs, bf64_t rhs, bf64_t expected) {
    const auto result = bf64_mul(lhs, rhs);
    BOOST_TEST(result == expected);
    const auto ntl_result = MulMod(from_bf(lhs), from_bf(rhs), bf64_residue);
    BOOST_TEST(ntl_result == from_bf(result));
  }
}

BOOST_AUTO_TEST_CASE(test_bf8_mul) {
  check_mul_bf8(0xFF, 0x00, 0x00);
  check_mul_bf8(0x00, 0xFF, 0x00);
  check_mul_bf8(0xFF, 0x01, 0xFF);
  check_mul_bf8(0x01, 0xFF, 0xFF);

  for (unsigned int i = 50; i; --i) {
    bf8_t lhs, rhs;
    rand_bytes(&lhs, sizeof(lhs));
    rand_bytes(&rhs, sizeof(rhs));

    const auto result = bf8_mul(lhs, rhs);
    const auto ntl_result = MulMod(from_bf(lhs), from_bf(rhs), bf8_residue);
    BOOST_TEST(ntl_result == from_bf(result));
  }
}

BOOST_AUTO_TEST_CASE(test_bf64_mul) {
  check_mul_bf64(0xFF, 0x00, 0x00);
  check_mul_bf64(0x00, 0xFF, 0x00);
  check_mul_bf64(0xFF, 0x01, 0xFF);
  check_mul_bf64(0x01, 0xFF, 0xFF);

  for (unsigned int i = 50; i; --i) {
    bf64_t lhs, rhs;
    rand_bytes(reinterpret_cast<uint8_t*>(&lhs), sizeof(lhs));
    rand_bytes(reinterpret_cast<uint8_t*>(&rhs), sizeof(rhs));

    const auto result = bf64_mul(lhs, rhs);
    const auto ntl_result = MulMod(from_bf(lhs), from_bf(rhs), bf64_residue);
    BOOST_TEST(ntl_result == from_bf(result));
  }
}

