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

BOOST_AUTO_TEST_SUITE_END();