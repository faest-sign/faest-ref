/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "fields.hpp"

#include <boost/test/unit_test.hpp>

namespace {
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

  template <class B>
  void add_invariants() {
    check_add<B>(0xFF, B::zero(), 0xFF);
    check_add<B>(B::zero(), 0xFF, 0xFF);
    check_add<B>(0xFF, 0xFF, B::zero());
  }

  template <class B>
  void add_random() {
    for (unsigned int i = 50; i; --i) {
      auto lhs = B::random();
      auto rhs = B::random();
      check_add(lhs, rhs);
    }
  }

  template <class B>
  void mul_invariants() {
    check_mul<B>(0xFF, B::zero(), B::zero());
    check_mul<B>(B::zero(), 0xFF, B::zero());
    check_mul<B>(0xFF, B::one(), 0xFF);
    check_mul<B>(B::one(), 0xFF, 0xFF);
    check_div<B>(0xFF, B::one(), 0xFF);
  }

  template <class B>
  void mul_random() {
    for (unsigned int i = 50; i; --i) {
      auto lhs = B::random();
      auto rhs = B::random();
      check_mul(lhs, rhs);
    }
  }

  template <class B>
  void div_random() {
    for (unsigned int i = 50; i; --i) {
      auto lhs = B::random();
      auto rhs = B::random();
      if (rhs == B::zero()) {
        continue;
      }
      check_div(lhs, rhs);
    }
  }
} // namespace

BOOST_AUTO_TEST_SUITE(fields);

BOOST_AUTO_TEST_CASE(test_bf8_add_invariants) {
  add_invariants<bf8>();
}

BOOST_AUTO_TEST_CASE(test_bf8_add_random) {
  add_random<bf8>();
}

BOOST_AUTO_TEST_CASE(test_bf8_mul_invariants) {
  mul_invariants<bf8>();
}

BOOST_AUTO_TEST_CASE(test_bf8_mul_random) {
  mul_random<bf8>();
}

BOOST_AUTO_TEST_CASE(test_bf8_div_random) {
  div_random<bf8>();
}

BOOST_AUTO_TEST_CASE(test_bf64_add_invariants) {
  add_invariants<bf64>();
}

BOOST_AUTO_TEST_CASE(test_bf64_add_random) {
  add_random<bf64>();
}

BOOST_AUTO_TEST_CASE(test_bf64_mul_invariants) {
  mul_invariants<bf64>();
}

BOOST_AUTO_TEST_CASE(test_bf64_mul_random) {
  mul_random<bf64>();
}

BOOST_AUTO_TEST_CASE(test_bf64_div_random) {
  div_random<bf64>();
}

BOOST_AUTO_TEST_CASE(test_bf128_add_invariants) {
  add_invariants<bf128>();
}

BOOST_AUTO_TEST_CASE(test_bf128_add_random) {
  add_random<bf128>();
}

BOOST_AUTO_TEST_CASE(test_bf128_mul_invariants) {
  mul_invariants<bf128>();
}

BOOST_AUTO_TEST_CASE(test_bf128_mul_random) {
  mul_random<bf128>();
}

BOOST_AUTO_TEST_CASE(test_bf128_div_random) {
  div_random<bf128>();
}

BOOST_AUTO_TEST_CASE(test_bf192_add_invariants) {
  add_invariants<bf192>();
}

BOOST_AUTO_TEST_CASE(test_bf192_add_random) {
  add_random<bf192>();
}

BOOST_AUTO_TEST_CASE(test_bf192_mul_invariants) {
  mul_invariants<bf192>();
}

BOOST_AUTO_TEST_CASE(test_bf192_mul_random) {
  mul_random<bf192>();
}

BOOST_AUTO_TEST_CASE(test_bf192_div_random) {
  div_random<bf192>();
}

BOOST_AUTO_TEST_CASE(test_bf256_add_invariants) {
  add_invariants<bf256>();
}

BOOST_AUTO_TEST_CASE(test_bf256_add_random) {
  add_random<bf256>();
}

BOOST_AUTO_TEST_CASE(test_bf256_mul_invariants) {
  mul_invariants<bf256>();
}

BOOST_AUTO_TEST_CASE(test_bf256_mul_random) {
  mul_random<bf256>();
}

BOOST_AUTO_TEST_CASE(test_bf256_div_random) {
  div_random<bf256>();
}

BOOST_AUTO_TEST_CASE(test_bf64_test_vectors) {
  constexpr bf64::bytes lhs{0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef};
  constexpr bf64::bytes rhs{0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01};
  constexpr bf64::bytes expected{0x96, 0xd6, 0xc9, 0x8c, 0x53, 0x13, 0x0c, 0x49};

  bf64 l{lhs};
  bf64 r{rhs};

  const auto o = (l * r).as_uint8();
  BOOST_TEST(o == expected);
}

BOOST_AUTO_TEST_CASE(test_bf128_test_vectors) {
  constexpr bf128::bytes lhs{0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef,
                             0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef};
  constexpr bf128::bytes rhs{0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01,
                             0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01};
  constexpr bf128::bytes expected{0xe2, 0x3d, 0x64, 0xab, 0xb2, 0x4c, 0x15, 0xda,
                                  0x43, 0x9c, 0xc5, 0xa,  0x13, 0xed, 0xb4, 0x7b};

  bf128 l{lhs};
  bf128 r{rhs};

  const auto o = (l * r).as_uint8();
  BOOST_TEST(o == expected);
}

BOOST_AUTO_TEST_CASE(test_bf192_test_vectors) {
  constexpr bf192::bytes lhs{0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef,
                             0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef,
                             0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef};
  constexpr bf192::bytes rhs{0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01,
                             0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01,
                             0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01};
  constexpr bf192::bytes expected{
      0x22, 0xdc, 0x85, 0x4a, 0x53, 0xad, 0xf4, 0x3b, 0xa2, 0x7d, 0x24, 0xeb,
      0xf2, 0x0c, 0x55, 0x9a, 0x03, 0xdc, 0x85, 0x4a, 0x53, 0xad, 0xf4, 0x3b,
  };

  bf192 l{lhs};
  bf192 r{rhs};

  const auto o = (l * r).as_uint8();
  BOOST_TEST(o == expected);
}

BOOST_AUTO_TEST_CASE(test_bf256_test_vectors) {
  constexpr bf256::bytes lhs{0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0x01, 0x23, 0x45,
                             0x67, 0x89, 0xab, 0xcd, 0xef, 0x01, 0x23, 0x45, 0x67, 0x89, 0xab,
                             0xcd, 0xef, 0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef};
  constexpr bf256::bytes rhs{0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01, 0xef, 0xcd, 0xab,
                             0x89, 0x67, 0x45, 0x23, 0x01, 0xef, 0xcd, 0xab, 0x89, 0x67, 0x45,
                             0x23, 0x01, 0xef, 0xcd, 0xab, 0x89, 0x67, 0x45, 0x23, 0x01};
  constexpr bf256::bytes expected{
      0x33, 0x5e, 0x5c, 0xd5, 0x5b, 0x5f, 0x50, 0xd9, 0x57, 0x5a, 0x54,
      0xdd, 0x53, 0x57, 0x58, 0xd1, 0x5f, 0x52, 0x5c, 0xd5, 0x5b, 0x5f,
      0x50, 0xd9, 0x57, 0x5a, 0x54, 0xdd, 0x53, 0x57, 0x58, 0xd1,
  };

  bf256 l{lhs};
  bf256 r{rhs};

  const auto o = (l * r).as_uint8();
  BOOST_TEST(o == expected);
}

BOOST_AUTO_TEST_SUITE_END();