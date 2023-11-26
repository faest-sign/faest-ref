/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_TESTS_FIELDS_HPP
#define FAEST_TESTS_FIELDS_HPP

#include "../fields.h"

#include <boost/format.hpp>
#include <algorithm>
#include <iterator>
#include <array>

#if defined(HAVE_NTL)
#include <NTL/GF2X.h>
using namespace NTL;
#endif

namespace {
#if defined(HAVE_NTL)
  static inline GF2X DivMod(const GF2X& lhs, const GF2X& rhs, const GF2XModulus& modulus) {
    const auto denom = InvMod(rhs, modulus);
    return MulMod(lhs, denom, modulus);
  }
#endif

  class bf8 {
    bf8_t value;

  public:
    bf8() : value{0} {}
    bf8(bf8_t v) : value{v} {}
    bf8(const bf8&) = default;

    bf8& operator=(const bf8&) = default;

    bf8& operator+=(bf8 other) {
      value = bf8_add(value, other.value);
      return *this;
    }

    bf8& operator-=(bf8 other) {
      return *this += other;
    }

    bf8& operator*=(bf8 other) {
      value = bf8_mul(value, other.value);
      return *this;
    }

    bf8 operator+(bf8 other) const {
      return {bf8_add(value, other.value)};
    }

    bf8 operator-(bf8 other) const {
      return *this + other;
    }

    bf8 operator*(bf8 other) const {
      return {bf8_mul(value, other.value)};
    }

    bf8 operator/(bf8 other) const {
      return {bf8_mul(value, bf8_inv(other.value))};
    }

    bool operator==(bf8 other) const {
      return value == other.value;
    }

#if defined(HAVE_NTL)
    GF2X as_ntl() const {
      GF2X ret;
      auto v = value;
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i, v & 1);
      }
      return ret;
    }
#endif

    bf8_t as_internal() const {
      return value;
    }

    std::array<uint8_t, 1> as_uint8() const {
      std::array<uint8_t, 1> ret;
      bf8_store(ret.data(), value);
      return ret;
    }

#if defined(HAVE_NTL)
    static GF2X ntl_residue() {
      GF2X residue;
      SetCoeff(residue, 8);
      SetCoeff(residue, 4);
      SetCoeff(residue, 3);
      SetCoeff(residue, 1);
      SetCoeff(residue, 0);
      return residue;
    }
#endif

    static bf8 random() {
      return bf8_rand();
    }

    static bf8 zero() {
      return bf8_zero();
    }

    static bf8 one() {
      return bf8_one();
    }
  };

  static inline std::ostream& operator<<(std::ostream& stream, bf8 v) {
    stream << boost::format("%02x") % v.as_internal();
    return stream;
  }

  class bf64 {
    bf64_t value;

  public:
    typedef std::array<uint8_t, 8> bytes;

    bf64() : value{0} {}
    bf64(bf64_t v) : value{v} {}
    bf64(const bytes& b) : value(bf64_load(b.data())) {}
    bf64(const bf64&) = default;

    bf64& operator=(const bf64&) = default;

    bf64& operator+=(bf64 other) {
      value = bf64_add(value, other.value);
      return *this;
    }

    bf64& operator-=(bf64 other) {
      return *this += other;
    }

    bf64& operator*=(bf64 other) {
      value = bf64_mul(value, other.value);
      return *this;
    }

    bf64 operator+(bf64 other) const {
      return {bf64_add(value, other.value)};
    }

    bf64 operator-(bf64 other) const {
      return *this + other;
    }

    bf64 operator*(bf64 other) const {
      return {bf64_mul(value, other.value)};
    }

    bool operator==(bf64 other) const {
      return value == other.value;
    }

#if defined(HAVE_NTL)
    GF2X as_ntl() const {
      GF2X ret;
      auto v = value;
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i, v & 1);
      }
      return ret;
    }
#endif

    bf64_t as_internal() const {
      return value;
    }

    bytes as_uint8() const {
      bytes ret;
      bf64_store(ret.data(), value);
      return ret;
    }

#if defined(HAVE_NTL)
    static GF2X ntl_residue() {
      GF2X residue;
      SetCoeff(residue, 64);
      SetCoeff(residue, 4);
      SetCoeff(residue, 3);
      SetCoeff(residue, 1);
      SetCoeff(residue, 0);
      return residue;
    }
#endif

    static bf64 random() {
      return bf64_rand();
    }

    static bf64 zero() {
      return bf64_zero();
    }

    static bf64 one() {
      return bf64_one();
    }
  };

  static inline std::ostream& operator<<(std::ostream& stream, bf64 v) {
    stream << boost::format("%08x") % v.as_internal();
    return stream;
  }

  class bf128 {
    bf128_t value;

  public:
    typedef std::array<uint8_t, BF128_NUM_BYTES> bytes;

    bf128() : value{0} {}
    bf128(uint64_t v) : value{bf128_from_bf64(v)} {}
    bf128(bf128_t v) : value{v} {}
    bf128(bf64 v) : value{bf128_from_bf64(v.as_internal())} {}
    bf128(const bytes& b) : value{bf128_load(b.data())} {}
    bf128(const bf128&) = default;

    bf128& operator=(const bf128&) = default;

    bf128& operator+=(bf128 other) {
      value = bf128_add(value, other.value);
      return *this;
    }

    bf128& operator-=(bf128 other) {
      return *this += other;
    }

    bf128& operator*=(bf128 other) {
      value = bf128_mul(value, other.value);
      return *this;
    }

    bf128 operator+(bf128 other) const {
      return {bf128_add(value, other.value)};
    }

    bf128 operator-(bf128 other) const {
      return *this + other;
    }

    bf128 operator*(bf128 other) const {
      return {bf128_mul(value, other.value)};
    }

    bf128 operator*(bf64 other) const {
      return {bf128_mul_64(value, other.as_internal())};
    }

    bool operator==(bf128 other) const {
      return BF_VALUE(value, 0) == BF_VALUE(other.value, 0) &&
             BF_VALUE(value, 1) == BF_VALUE(other.value, 1);
    }

#if defined(HAVE_NTL)
    GF2X as_ntl() const {
      GF2X ret;
      auto v = BF_VALUE(value, 0);
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i, v & 1);
      }
      v = BF_VALUE(value, 1);
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 64, v & 1);
      }
      return ret;
    }
#endif

    bf128_t as_internal() const {
      return value;
    }

    bytes as_uint8() const {
      bytes ret;
      bf128_store(ret.data(), value);
      return ret;
    }

#if defined(HAVE_NTL)
    static GF2X ntl_residue() {
      GF2X residue;
      SetCoeff(residue, 128);
      SetCoeff(residue, 7);
      SetCoeff(residue, 2);
      SetCoeff(residue, 1);
      SetCoeff(residue, 0);
      return residue;
    }
#endif

    static bf128 random() {
      return bf128_rand();
    }

    static bf128 zero() {
      return bf128_zero();
    }

    static bf128 one() {
      return bf128_one();
    }
  };

  static inline std::ostream& operator<<(std::ostream& stream, bf128 v) {
    const auto value = v.as_internal();
    stream << boost::format("%08x %08x") % BF_VALUE(value, 1) % BF_VALUE(value, 0);
    return stream;
  }

  class bf192 {
    bf192_t value;

  public:
    typedef std::array<uint8_t, BF192_NUM_BYTES> bytes;

    bf192() : value{0} {}
    bf192(uint64_t v) : value{bf192_from_bf64(v)} {}
    bf192(bf192_t v) : value{v} {}
    bf192(bf64 v) : value{bf192_from_bf64(v.as_internal())} {}
    bf192(const bytes& b) : value{bf192_load(b.data())} {}
    bf192(const bf192&) = default;

    bf192& operator=(const bf192&) = default;

    bf192& operator+=(bf192 other) {
      value = bf192_add(value, other.value);
      return *this;
    }

    bf192& operator-=(bf192 other) {
      return *this += other;
    }

    bf192& operator*=(bf192 other) {
      value = bf192_mul(value, other.value);
      return *this;
    }

    bf192 operator+(bf192 other) const {
      return {bf192_add(value, other.value)};
    }

    bf192 operator-(bf192 other) const {
      return *this + other;
    }

    bf192 operator*(bf192 other) const {
      return {bf192_mul(value, other.value)};
    }

    bf192 operator*(bf64 other) const {
      return {bf192_mul_64(value, other.as_internal())};
    }

    bool operator==(bf192 other) const {
      return BF_VALUE(value, 0) == BF_VALUE(other.value, 0) &&
             BF_VALUE(value, 1) == BF_VALUE(other.value, 1) &&
             BF_VALUE(value, 2) == BF_VALUE(other.value, 2);
    }

#if defined(HAVE_NTL)
    GF2X as_ntl() const {
      GF2X ret;
      auto v = BF_VALUE(value, 0);
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i, v & 1);
      }
      v = BF_VALUE(value, 1);
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 64, v & 1);
      }
      v = BF_VALUE(value, 2);
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 128, v & 1);
      }
      return ret;
    }
#endif

    bf192_t as_internal() const {
      return value;
    }

    bytes as_uint8() const {
      bytes ret;
      bf192_store(ret.data(), value);
      return ret;
    }

#if defined(HAVE_NTL)
    static GF2X ntl_residue() {
      GF2X residue;
      SetCoeff(residue, 192);
      SetCoeff(residue, 7);
      SetCoeff(residue, 2);
      SetCoeff(residue, 1);
      SetCoeff(residue, 0);
      return residue;
    }
#endif

    static bf192 random() {
      return bf192_rand();
    }

    static bf192 zero() {
      return bf192_zero();
    }

    static bf192 one() {
      return bf192_one();
    }
  };

  static inline std::ostream& operator<<(std::ostream& stream, bf192 v) {
    const auto value = v.as_internal();
    stream << boost::format("%08x %08x %08x") % BF_VALUE(value, 2) % BF_VALUE(value, 1) %
                  BF_VALUE(value, 0);
    return stream;
  }

  class bf256 {
    bf256_t value;

  public:
    typedef std::array<uint8_t, BF256_NUM_BYTES> bytes;

    bf256() : value{0} {}
    bf256(uint64_t v) : value{bf256_from_bf64(v)} {}
    bf256(bf256_t v) : value{v} {}
    bf256(bf64 v) : value{bf256_from_bf64(v.as_internal())} {}
    bf256(const bytes& b) : value{bf256_load(b.data())} {}
    bf256(const bf256&) = default;

    bf256& operator=(const bf256&) = default;

    bf256& operator+=(bf256 other) {
      value = bf256_add(value, other.value);
      return *this;
    }

    bf256& operator-=(bf256 other) {
      return *this += other;
    }

    bf256& operator*=(bf256 other) {
      value = bf256_mul(value, other.value);
      return *this;
    }

    bf256 operator+(bf256 other) const {
      return {bf256_add(value, other.value)};
    }

    bf256 operator-(bf256 other) const {
      return *this + other;
    }

    bf256 operator*(bf256 other) const {
      return {bf256_mul(value, other.value)};
    }

    bf256 operator*(bf64 other) const {
      return {bf256_mul_64(value, other.as_internal())};
    }

    bool operator==(bf256 other) const {
      return BF_VALUE(value, 0) == BF_VALUE(other.value, 0) &&
             BF_VALUE(value, 1) == BF_VALUE(other.value, 1) &&
             BF_VALUE(value, 2) == BF_VALUE(other.value, 2) &&
             BF_VALUE(value, 3) == BF_VALUE(other.value, 3);
    }

#if defined(HAVE_NTL)
    GF2X as_ntl() const {
      GF2X ret;
      auto v = BF_VALUE(value, 0);
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i, v & 1);
      }
      v = BF_VALUE(value, 1);
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 64, v & 1);
      }
      v = BF_VALUE(value, 2);
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 128, v & 1);
      }
      v = BF_VALUE(value, 3);
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 192, v & 1);
      }
      return ret;
    }
#endif

    bf256_t as_internal() const {
      return value;
    }

    bytes as_uint8() const {
      bytes ret;
      bf256_store(ret.data(), value);
      return ret;
    }

#if defined(HAVE_NTL)
    static GF2X ntl_residue() {
      GF2X residue;
      SetCoeff(residue, 256);
      SetCoeff(residue, 10);
      SetCoeff(residue, 5);
      SetCoeff(residue, 2);
      SetCoeff(residue, 0);
      return residue;
    }
#endif

    static bf256 random() {
      return bf256_rand();
    }

    static bf256 zero() {
      return bf256_zero();
    }

    static bf256 one() {
      return bf256_one();
    }
  };

  static inline std::ostream& operator<<(std::ostream& stream, bf256 v) {
    const auto value = v.as_internal();
    stream << boost::format("%08x %08x %08x %08x") % BF_VALUE(value, 3) % BF_VALUE(value, 2) %
                  BF_VALUE(value, 1) % BF_VALUE(value, 0);
    return stream;
  }
} // namespace

#endif
