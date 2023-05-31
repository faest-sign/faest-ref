#ifndef TEST_AES_HPP
#define TEST_AES_HPP

#include <array>
#include <cstdint>

namespace aes_ctr_128_tv {
  extern const std::array<uint8_t, 16> key;
  extern const std::array<uint8_t, 16> in;
  extern const std::array<uint8_t, 16> out;
  extern const std::array<uint8_t, 200> expected_extended_witness;
}

namespace aes_ctr_192_tv {
  extern const std::array<uint8_t, 24> key;
  extern const std::array<uint8_t, 32> in;
  extern const std::array<uint8_t, 32> out;
  extern const std::array<uint8_t, 408> expected_extended_witness;
}

namespace aes_ctr_256_tv {
  extern const std::array<uint8_t, 32> key;
  extern const std::array<uint8_t, 32> in;
  extern const std::array<uint8_t, 32> out;
  extern const std::array<uint8_t, 500> expected_extended_witness;
}

#endif // TEST_AES_HPP
