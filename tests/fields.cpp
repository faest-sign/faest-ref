/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../fields.h"
#include "../randomness.h"

#define BOOST_TEST_MODULE fields_test
#include <NTL/GF2X.h>
#include <boost/format.hpp>
#include <boost/test/included/unit_test.hpp>

using namespace NTL;

namespace {
  GF2X DivMod(const GF2X& lhs, const GF2X& rhs, const GF2XModulus& modulus) {
    const auto denom = InvMod(rhs, modulus);
    return MulMod(lhs, denom, modulus);
  }

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

    GF2X as_ntl() const {
      GF2X ret;
      auto v = value;
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i, v & 1);
      }
      return ret;
    }

    bf8_t as_internal() const {
      return value;
    }

    static GF2X ntl_residue() {
      GF2X residue;
      SetCoeff(residue, 8);
      SetCoeff(residue, 4);
      SetCoeff(residue, 3);
      SetCoeff(residue, 1);
      SetCoeff(residue, 0);
      return residue;
    }

    static bf8 random() {
      return bf8_rand();
    }
  };

  std::ostream& operator<<(std::ostream& stream, bf8 v) {
    stream << boost::format("%02x") % v.as_internal();
    return stream;
  }

  class bf64 {
    bf64_t value;

  public:
    bf64() : value{0} {}
    bf64(bf64_t v) : value{v} {}
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

    bf64 operator/(bf64 other) const {
      return {bf64_mul(value, bf64_inv(other.value))};
    }

    bool operator==(bf64 other) const {
      return value == other.value;
    }

    GF2X as_ntl() const {
      GF2X ret;
      auto v = value;
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i, v & 1);
      }
      return ret;
    }

    bf64_t as_internal() const {
      return value;
    }

    static GF2X ntl_residue() {
      GF2X residue;
      SetCoeff(residue, 64);
      SetCoeff(residue, 4);
      SetCoeff(residue, 3);
      SetCoeff(residue, 1);
      SetCoeff(residue, 0);
      return residue;
    }

    static bf64 random() {
      return bf64_rand();
    }
  };

  std::ostream& operator<<(std::ostream& stream, bf64 v) {
    stream << boost::format("%08x") % v.as_internal();
    return stream;
  }

  class bf128 {
    bf128_t value;

  public:
    bf128() : value{0} {}
    bf128(uint64_t v) : value{{v, 0}} {}
    bf128(bf128_t v) : value{v} {}
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

    bf128 operator/(bf128 other) const {
      return {bf128_mul(value, bf128_inv(other.value))};
    }

    bool operator==(bf128 other) const {
      return value.values[0] == other.value.values[0] && value.values[1] == other.value.values[1];
    }

    GF2X as_ntl() const {
      GF2X ret;
      auto v = value.values[0];
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i, v & 1);
      }
      v = value.values[1];
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 64, v & 1);
      }
      return ret;
    }

    bf128_t as_internal() const {
      return value;
    }

    static GF2X ntl_residue() {
      GF2X residue;
      SetCoeff(residue, 128);
      SetCoeff(residue, 7);
      SetCoeff(residue, 2);
      SetCoeff(residue, 1);
      SetCoeff(residue, 0);
      return residue;
    }

    static bf128 random() {
      return bf128_rand();
    }
  };

  std::ostream& operator<<(std::ostream& stream, bf128 v) {
    auto value = v.as_internal();
    stream << boost::format("%08x %08x") % value.values[1] % value.values[0];
    return stream;
  }

  class bf192 {
    bf192_t value;

  public:
    bf192() : value{0} {}
    bf192(uint64_t v) : value{{v, 0, 0}} {}
    bf192(bf192_t v) : value{v} {}
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

    bf192 operator/(bf192 other) const {
      return {bf192_mul(value, bf192_inv(other.value))};
    }

    bool operator==(bf192 other) const {
      return value.values[0] == other.value.values[0] && value.values[1] == other.value.values[1] &&
             value.values[2] == other.value.values[2];
    }

    GF2X as_ntl() const {
      GF2X ret;
      auto v = value.values[0];
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i, v & 1);
      }
      v = value.values[1];
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 64, v & 1);
      }
      v = value.values[2];
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 128, v & 1);
      }
      return ret;
    }

    bf192_t as_internal() const {
      return value;
    }

    static GF2X ntl_residue() {
      GF2X residue;
      SetCoeff(residue, 192);
      SetCoeff(residue, 7);
      SetCoeff(residue, 2);
      SetCoeff(residue, 1);
      SetCoeff(residue, 0);
      return residue;
    }

    static bf192 random() {
      return bf192_rand();
    }
  };

  std::ostream& operator<<(std::ostream& stream, bf192 v) {
    auto value = v.as_internal();
    stream << boost::format("%08x %08x %08x") % value.values[2] % value.values[1] % value.values[0];
    return stream;
  }

  class bf256 {
    bf256_t value;

  public:
    bf256() : value{0} {}
    bf256(uint64_t v) : value{{v, 0, 0, 0}} {}
    bf256(bf256_t v) : value{v} {}
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

    bf256 operator/(bf256 other) const {
      return {bf256_mul(value, bf256_inv(other.value))};
    }

    bool operator==(bf256 other) const {
      return value.values[0] == other.value.values[0] && value.values[1] == other.value.values[1] &&
             value.values[2] == other.value.values[2] && value.values[3] == other.value.values[3];
    }

    GF2X as_ntl() const {
      GF2X ret;
      auto v = value.values[0];
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i, v & 1);
      }
      v = value.values[1];
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 64, v & 1);
      }
      v = value.values[2];
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 128, v & 1);
      }
      v = value.values[3];
      for (unsigned int i = 0; i != sizeof(v) * 8 && v; ++i, v >>= 1) {
        SetCoeff(ret, i + 192, v & 1);
      }
      return ret;
    }

    bf256_t as_internal() const {
      return value;
    }

    static GF2X ntl_residue() {
      GF2X residue;
      SetCoeff(residue, 256);
      SetCoeff(residue, 10);
      SetCoeff(residue, 5);
      SetCoeff(residue, 2);
      SetCoeff(residue, 0);
      return residue;
    }

    static bf256 random() {
      return bf256_rand();
    }
  };

  std::ostream& operator<<(std::ostream& stream, bf256 v) {
    auto value = v.as_internal();
    stream << boost::format("%08x %08x %08x %08x") % value.values[3] % value.values[2] %
                  value.values[1] % value.values[0];
    return stream;
  }

  template <class B>
  void check_add(B lhs, B rhs, B expected) {
    const auto result = lhs + rhs;
    BOOST_TEST(result == expected);
    const auto ntl_result = lhs.as_ntl() + rhs.as_ntl();
    BOOST_TEST(ntl_result == result.as_ntl());
  }

  template <class B>
  void check_add(B lhs, B rhs) {
    const auto result     = lhs + rhs;
    const auto ntl_result = lhs.as_ntl() + rhs.as_ntl();
    BOOST_TEST(ntl_result == result.as_ntl());
  }

  template <class B>
  void check_mul(B lhs, B rhs, B expected) {
    const auto result = lhs * rhs;
    BOOST_TEST(result == expected);
    const auto ntl_result = MulMod(lhs.as_ntl(), rhs.as_ntl(), B::ntl_residue());
    BOOST_TEST(ntl_result == result.as_ntl());
  }

  template <class B>
  void check_mul(B lhs, B rhs) {
    const auto result     = lhs * rhs;
    const auto ntl_result = MulMod(lhs.as_ntl(), rhs.as_ntl(), B::ntl_residue());
    BOOST_TEST(ntl_result == result.as_ntl());
  }

  template <class B>
  void check_div(B lhs, B rhs, B expected) {
    const auto result = lhs / rhs;
    BOOST_TEST(result == expected);
    const auto ntl_result = DivMod(lhs.as_ntl(), rhs.as_ntl(), B::ntl_residue());
    BOOST_TEST(ntl_result == result.as_ntl());
  }

  template <class B>
  void check_div(B lhs, B rhs) {
    const auto result     = lhs / rhs;
    const auto ntl_result = DivMod(lhs.as_ntl(), rhs.as_ntl(), B::ntl_residue());
    BOOST_TEST(ntl_result == result.as_ntl());
  }
} // namespace

BOOST_AUTO_TEST_CASE(test_bf8_add_invariants) {
  check_add<bf8>(0xFF, 0x00, 0xFF);
  check_add<bf8>(0x00, 0xFF, 0xFF);
  check_add<bf8>(0xFF, 0xFF, 0x00);
}

BOOST_AUTO_TEST_CASE(test_bf8_add_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf8::random();
    auto rhs = bf8::random();
    check_add(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf8_mul_invariants) {
  check_mul<bf8>(0xFF, 0x00, 0x00);
  check_mul<bf8>(0x00, 0xFF, 0x00);
  check_mul<bf8>(0xFF, 0x01, 0xFF);
  check_mul<bf8>(0x01, 0xFF, 0xFF);
  check_div<bf8>(0xFF, 0x01, 0xFF);
}

BOOST_AUTO_TEST_CASE(test_bf8_mul_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf8::random();
    auto rhs = bf8::random();
    check_mul(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf8_div_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf8::random();
    auto rhs = bf8::random();
    if (rhs == bf8{}) {
      continue;
    }
    check_div(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf64_add_invariants) {
  check_add<bf64>(0xFF, 0x00, 0xFF);
  check_add<bf64>(0x00, 0xFF, 0xFF);
  check_add<bf64>(0xFF, 0xFF, 0x00);
}

BOOST_AUTO_TEST_CASE(test_bf64_add_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf64::random();
    auto rhs = bf64::random();
    check_add(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf64_mul_invariants) {
  check_mul<bf64>(0xFF, 0x00, 0x00);
  check_mul<bf64>(0x00, 0xFF, 0x00);
  check_mul<bf64>(0xFF, 0x01, 0xFF);
  check_mul<bf64>(0x01, 0xFF, 0xFF);
  check_div<bf64>(0xFF, 0x01, 0xFF);
}

BOOST_AUTO_TEST_CASE(test_bf64_mul_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf64::random();
    auto rhs = bf64::random();
    check_mul(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf64_div_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf64::random();
    auto rhs = bf64::random();
    if (rhs == bf64{}) {
      continue;
    }
    check_div(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf128_add_invariants) {
  check_add<bf128>(0xFF, 0x00, 0xFF);
  check_add<bf128>(0x00, 0xFF, 0xFF);
  check_add<bf128>(0xFF, 0xFF, 0x00);
  check_div<bf64>(0xFF, 0x01, 0xFF);
}

BOOST_AUTO_TEST_CASE(test_bf128_add_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf128::random();
    auto rhs = bf128::random();
    check_add(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf128_mul_invariants) {
  check_mul<bf128>(0xFF, 0x00, 0x00);
  check_mul<bf128>(0x00, 0xFF, 0x00);
  check_mul<bf128>(0xFF, 0x01, 0xFF);
  check_mul<bf128>(0x01, 0xFF, 0xFF);
  check_div<bf128>(0xFF, 0x01, 0xFF);
}

BOOST_AUTO_TEST_CASE(test_bf128_mul_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf128::random();
    auto rhs = bf128::random();
    check_mul(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf128_div_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf128::random();
    auto rhs = bf128::random();
    if (rhs == bf128{}) {
      continue;
    }
    check_div(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf192_add_invariants) {
  check_add<bf192>(0xFF, 0x00, 0xFF);
  check_add<bf192>(0x00, 0xFF, 0xFF);
  check_add<bf192>(0xFF, 0xFF, 0x00);
}

BOOST_AUTO_TEST_CASE(test_bf192_add_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf192::random();
    auto rhs = bf192::random();
    check_add(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf192_mul_invariants) {
  check_mul<bf192>(0xFF, 0x00, 0x00);
  check_mul<bf192>(0x00, 0xFF, 0x00);
  check_mul<bf192>(0xFF, 0x01, 0xFF);
  check_mul<bf192>(0x01, 0xFF, 0xFF);
  check_div<bf192>(0xFF, 0x01, 0xFF);
}

BOOST_AUTO_TEST_CASE(test_bf192_mul_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf192::random();
    auto rhs = bf192::random();
    check_mul(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf192_div_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf192::random();
    auto rhs = bf192::random();
    if (rhs == bf192{}) {
      continue;
    }
    check_div(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf256_add_invariants) {
  check_add<bf256>(0xFF, 0x00, 0xFF);
  check_add<bf256>(0x00, 0xFF, 0xFF);
  check_add<bf256>(0xFF, 0xFF, 0x00);
}

BOOST_AUTO_TEST_CASE(test_bf256_add_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf256::random();
    auto rhs = bf256::random();
    check_add(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf256_mul_invariants) {
  check_mul<bf256>(0xFF, 0x00, 0x00);
  check_mul<bf256>(0x00, 0xFF, 0x00);
  check_mul<bf256>(0xFF, 0x01, 0xFF);
  check_mul<bf256>(0x01, 0xFF, 0xFF);
  check_div<bf256>(0xFF, 0x01, 0xFF);
}

BOOST_AUTO_TEST_CASE(test_bf256_mul_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf256::random();
    auto rhs = bf256::random();
    check_mul(lhs, rhs);
  }
}

BOOST_AUTO_TEST_CASE(test_bf256_div_random) {
  for (unsigned int i = 50; i; --i) {
    auto lhs = bf256::random();
    auto rhs = bf256::random();
    if (rhs == bf256{}) {
      continue;
    }
    check_div(lhs, rhs);
  }
}
