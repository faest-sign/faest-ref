#ifndef TEST_AES_HPP
#define TEST_AES_HPP

#include <array>
#include <cstdint>

namespace aes_ctr_128_tv {
  extern const uint8_t key[16];
  extern const uint8_t in[16];
  extern const uint8_t out[16];
  extern const std::array<uint8_t, 200> expected_extended_witness;
}

namespace aes_ctr_192_tv {
  extern const uint8_t key[24];
  extern const uint8_t in[32];
  extern const uint8_t out[32];
  extern const std::array<uint8_t, 408> expected_extended_witness;
}

namespace aes_ctr_256_tv {
  extern const uint8_t key[32];
  extern const uint8_t in[32];
  extern const uint8_t out[32];
  extern const std::array<uint8_t, 500> expected_extended_witness;
}

#endif // TEST_AES_HPP
