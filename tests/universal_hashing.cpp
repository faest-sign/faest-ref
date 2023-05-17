#include "../universal_hashing.h"

#include "fields.hpp"

#include <boost/test/unit_test.hpp>
#include <array>
#include <vector>

namespace {
  static constexpr size_t xs = 8;

} // namespace

BOOST_AUTO_TEST_SUITE(universal_hashing)

BOOST_AUTO_TEST_CASE(test_vole_hash_128) {
  auto r0 = bf128::random().as_uint8();
  auto r1 = bf128::random().as_uint8();
  auto s  = bf128::random().as_uint8();
  auto t  = bf64::random().as_uint8();

  std::vector<uint8_t> x;
  for (unsigned int i = 0; i != xs; ++i) {
    auto tmp = bf128::random().as_uint8();
    x.insert(x.end(), tmp.begin(), tmp.end());
  }

  std::array<uint8_t, 16 + UNIVERSAL_HASH_B> digest;
  vole_hash_128(digest.data(), r0.data(), r1.data(), s.data(), t.data(), x.data(), xs * 128);
  BOOST_TEST(digest != decltype(digest){});
}

BOOST_AUTO_TEST_CASE(test_vole_hash_192) {
  auto r0 = bf192::random().as_uint8();
  auto r1 = bf192::random().as_uint8();
  auto s  = bf192::random().as_uint8();
  auto t  = bf64::random().as_uint8();

  std::vector<uint8_t> x;
  for (unsigned int i = 0; i != xs; ++i) {
    auto tmp = bf192::random().as_uint8();
    x.insert(x.end(), tmp.begin(), tmp.end());
  }

  std::array<uint8_t, 24 + UNIVERSAL_HASH_B> digest;
  vole_hash_192(digest.data(), r0.data(), r1.data(), s.data(), t.data(), x.data(), xs * 192);
  BOOST_TEST(digest != decltype(digest){});
}

BOOST_AUTO_TEST_CASE(test_vole_hash_256) {
  auto r0 = bf256::random().as_uint8();
  auto r1 = bf256::random().as_uint8();
  auto s  = bf256::random().as_uint8();
  auto t  = bf64::random().as_uint8();

  std::vector<uint8_t> x;
  for (unsigned int i = 0; i != xs; ++i) {
    auto tmp = bf256::random().as_uint8();
    x.insert(x.end(), tmp.begin(), tmp.end());
  }

  std::array<uint8_t, 32 + UNIVERSAL_HASH_B> digest;
  vole_hash_256(digest.data(), r0.data(), r1.data(), s.data(), t.data(), x.data(), xs * 256);
  BOOST_TEST(digest != decltype(digest){});
}

BOOST_AUTO_TEST_CASE(test_zk_hash_128) {
  auto r = bf128::random().as_uint8();
  auto s = bf128::random().as_uint8();
  auto t = bf64::random().as_uint8();

  std::array<bf128_t, xs> x;
  for (unsigned int i = 0; i != xs; ++i) {
    x[i] = bf128::random().as_internal();
  }

  std::array<uint8_t, 16> digest;
  zk_hash_128(digest.data(), r.data(), s.data(), t.data(), x.data(), xs);
  BOOST_TEST(digest != decltype(digest){});
}

BOOST_AUTO_TEST_CASE(test_zk_hash_192) {
  auto r = bf192::random().as_uint8();
  auto s = bf192::random().as_uint8();
  auto t = bf64::random().as_uint8();

  std::array<bf192_t, xs> x;
  for (unsigned int i = 0; i != xs; ++i) {
    x[i] = bf192::random().as_internal();
  }

  std::array<uint8_t, 16> digest;
  zk_hash_192(digest.data(), r.data(), s.data(), t.data(), x.data(), xs);
  BOOST_TEST(digest != decltype(digest){});
}

BOOST_AUTO_TEST_CASE(test_zk_hash_256) {
  auto r = bf256::random().as_uint8();
  auto s = bf256::random().as_uint8();
  auto t = bf64::random().as_uint8();

  std::array<bf256_t, xs> x;
  for (unsigned int i = 0; i != xs; ++i) {
    x[i] = bf256::random().as_internal();
  }

  std::array<uint8_t, 16> digest;
  zk_hash_256(digest.data(), r.data(), s.data(), t.data(), x.data(), xs);
  BOOST_TEST(digest != decltype(digest){});
}

BOOST_AUTO_TEST_SUITE_END()
