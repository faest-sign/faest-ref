#ifndef TESTS_PRG_TVS_HPP
#define TESTS_PRG_TVS_HPP

#include <array>
#include <cstdint>

namespace prg_tvs {
  constexpr std::array<uint8_t, 16> key_128{
      0xd0, 0x22, 0xe7, 0xd5, 0x20, 0xf8, 0xe9, 0x38,
      0xa1, 0x4e, 0x18, 0x8c, 0x47, 0x30, 0x8c, 0xfe,
  };
  constexpr std::array<uint8_t, 16> iv_128{
      0xf5, 0xff, 0xf7, 0xf7, 0x28, 0xb9, 0xf8, 0xfb,
      0xf5, 0x1c, 0x7c, 0xcc, 0xcc, 0x4c, 0x24, 0x01,
  };
  constexpr uint32_t tweak_128 = 1811450929;
  constexpr std::array<uint8_t, 240> expected_128{
      0xbb, 0xdb, 0xe5, 0xfd, 0x57, 0xbc, 0x0c, 0x84, 0x4b, 0xdd, 0x31, 0x71, 0x63, 0xbd, 0x89,
      0x73, 0x75, 0x6e, 0x01, 0x70, 0x17, 0x0e, 0xf2, 0xb0, 0x26, 0x4c, 0x4f, 0xd6, 0x82, 0xc9,
      0xec, 0x17, 0x51, 0x61, 0x80, 0xac, 0x06, 0x05, 0xb1, 0x5f, 0x68, 0xbe, 0x65, 0xcf, 0x1f,
      0xfc, 0x67, 0xf0, 0xae, 0x0c, 0x46, 0x79, 0x31, 0xe3, 0x68, 0xdb, 0x5f, 0x01, 0x04, 0xae,
      0x63, 0x2f, 0x4e, 0x70, 0x69, 0x68, 0x5c, 0x6a, 0x00, 0x80, 0x4b, 0xf3, 0x37, 0x9a, 0xec,
      0x73, 0x76, 0x46, 0x3f, 0x82, 0xa1, 0x7e, 0x40, 0xc4, 0x33, 0x32, 0xf9, 0xd0, 0x53, 0xac,
      0xe0, 0x90, 0x2d, 0xef, 0xd4, 0x5a, 0xa2, 0xac, 0x3c, 0x75, 0x56, 0xc5, 0x13, 0x69, 0xad,
      0x8d, 0xe4, 0xad, 0x29, 0x3a, 0x1b, 0x60, 0x7a, 0x0c, 0x0a, 0x90, 0xc6, 0xfe, 0x54, 0xc1,
      0xa5, 0x2d, 0xdb, 0xdf, 0xb3, 0x4a, 0x59, 0xb3, 0xc9, 0x9e, 0x40, 0x53, 0x5e, 0x87, 0xaa,
      0x07, 0xdb, 0x23, 0x0c, 0x17, 0x58, 0x6a, 0xdd, 0x7e, 0xeb, 0xa0, 0x5e, 0x1e, 0x86, 0xac,
      0xfd, 0x9f, 0x6d, 0x2c, 0x5c, 0xa4, 0x94, 0x6f, 0xa3, 0x05, 0xa7, 0x22, 0x10, 0x1f, 0x79,
      0x6e, 0xe7, 0x6e, 0xfc, 0xf6, 0x2f, 0x93, 0xf5, 0xf0, 0x36, 0x1f, 0xc2, 0x77, 0xf9, 0xaf,
      0x10, 0xf7, 0x91, 0x47, 0xc7, 0xce, 0x45, 0x07, 0x73, 0x70, 0xfa, 0x1b, 0x2a, 0x2f, 0xab,
      0x51, 0x52, 0x07, 0x2d, 0x2a, 0xd6, 0x79, 0xf7, 0x39, 0x38, 0x8f, 0xbb, 0x83, 0x7f, 0x19,
      0x15, 0xcf, 0x89, 0x70, 0xd5, 0x52, 0xc5, 0x33, 0x45, 0xd4, 0x92, 0xf1, 0x6a, 0x0f, 0x6b,
      0xf8, 0xde, 0x65, 0x25, 0xe8, 0x5d, 0x58, 0xf2, 0x63, 0x61, 0xcc, 0x49, 0xaa, 0xa3, 0xba,
  };
  constexpr std::array<uint8_t, 24> key_192{
      0x1c, 0xea, 0xa3, 0xca, 0xe0, 0xf5, 0x80, 0xa7, 0xcc, 0x09, 0x5c, 0xd9,
      0x36, 0xef, 0xae, 0xad, 0x66, 0xc1, 0xbd, 0xbe, 0x79, 0x64, 0x6c, 0xa7,
  };
  constexpr std::array<uint8_t, 16> iv_192{
      0x2c, 0x2b, 0x4d, 0xb4, 0xcc, 0x08, 0x51, 0x46,
      0xdf, 0x0b, 0x26, 0x18, 0xfe, 0xd2, 0xd2, 0xb1,
  };
  constexpr uint32_t tweak_192 = 537655879;
  constexpr std::array<uint8_t, 240> expected_192{
      0xc5, 0x7a, 0x34, 0x03, 0xc6, 0x3c, 0x48, 0x76, 0x3c, 0x8a, 0xfa, 0xe0, 0xff, 0x42, 0x86,
      0x02, 0xdb, 0x8c, 0x77, 0x51, 0x17, 0x38, 0x96, 0x08, 0x72, 0x23, 0xef, 0x73, 0xba, 0x3b,
      0xb1, 0x1e, 0x82, 0x28, 0xde, 0xab, 0xf3, 0xe8, 0x50, 0xfa, 0xfd, 0xc0, 0x1f, 0x6c, 0xbb,
      0xec, 0xac, 0x1b, 0x61, 0x8d, 0x7f, 0xfa, 0xf8, 0x92, 0xb1, 0xdb, 0xc7, 0xea, 0x52, 0xdf,
      0xd3, 0xc5, 0xb6, 0x34, 0x23, 0x70, 0x52, 0x4d, 0xac, 0x1f, 0xf9, 0xa7, 0x93, 0x57, 0x56,
      0x19, 0xc3, 0x03, 0x2f, 0xa6, 0x6c, 0xb9, 0x4d, 0x90, 0xe0, 0x89, 0x11, 0x43, 0x25, 0x2e,
      0xdc, 0x7f, 0x57, 0x7e, 0x1c, 0x9e, 0xd7, 0x18, 0xd6, 0x57, 0x37, 0x03, 0xe7, 0x47, 0xd3,
      0xfc, 0x22, 0xc0, 0x87, 0x93, 0x78, 0x65, 0x00, 0xac, 0x58, 0xaf, 0x22, 0xf8, 0x46, 0xab,
      0xda, 0xea, 0x6a, 0xd1, 0xe2, 0x94, 0x9c, 0x92, 0xb2, 0x9f, 0x68, 0xe6, 0x5f, 0x87, 0x54,
      0x98, 0x46, 0x70, 0x5c, 0x83, 0xab, 0x2a, 0x15, 0xea, 0x27, 0xb9, 0xc6, 0x32, 0x3b, 0x53,
      0x55, 0xe5, 0x09, 0x02, 0x05, 0xcb, 0x03, 0x4c, 0x0b, 0xe6, 0x26, 0x7a, 0x23, 0xf0, 0x6b,
      0x2e, 0xd2, 0x02, 0x4c, 0x0f, 0x97, 0xe1, 0x96, 0x4f, 0x32, 0xe1, 0xe0, 0x28, 0xf8, 0x72,
      0x25, 0xc7, 0x6f, 0xcb, 0x92, 0x56, 0x9b, 0x4b, 0xe9, 0xad, 0xa6, 0x8e, 0xbe, 0x60, 0xc5,
      0xf6, 0xf6, 0x5c, 0x83, 0x04, 0x90, 0x27, 0xcf, 0xdf, 0x12, 0x26, 0xe1, 0xd2, 0x47, 0xac,
      0xec, 0x9e, 0x70, 0x8f, 0x16, 0xa7, 0xa7, 0x4e, 0x2e, 0x13, 0x70, 0xe6, 0xde, 0x0b, 0x35,
      0xc4, 0x19, 0x8b, 0xcb, 0x41, 0x8f, 0x42, 0x85, 0x99, 0xfb, 0x44, 0x62, 0x78, 0xaf, 0x44,
  };
  constexpr std::array<uint8_t, 32> key_256{
      0x51, 0xc3, 0xf3, 0x7d, 0x08, 0xa9, 0x70, 0x20, 0x61, 0x35, 0xc3,
      0x0d, 0xcb, 0x09, 0x2f, 0x68, 0x7d, 0x75, 0x72, 0x7c, 0xa5, 0xcb,
      0xb5, 0xeb, 0xc1, 0xce, 0x46, 0xb4, 0xae, 0x00, 0xa7, 0xb5,
  };
  constexpr std::array<uint8_t, 16> iv_256{
      0x29, 0xa4, 0x1e, 0x74, 0x7f, 0xc6, 0xf5, 0x92,
      0x57, 0xe0, 0x95, 0xce, 0x39, 0x04, 0xc0, 0xd2,
  };
  constexpr uint32_t tweak_256 = 1095625157;
  constexpr std::array<uint8_t, 240> expected_256{
      0x07, 0x9c, 0x25, 0xa2, 0xa1, 0xe2, 0x47, 0x8e, 0x3e, 0x59, 0x5f, 0x9b, 0x1f, 0xbc, 0x11,
      0x5b, 0xc6, 0x4b, 0xe3, 0x40, 0x9e, 0xae, 0xb7, 0xe4, 0x2d, 0xa4, 0x20, 0x79, 0x9c, 0x96,
      0x09, 0x3a, 0x2a, 0xa8, 0x8b, 0x05, 0xd7, 0x9e, 0x95, 0xb4, 0x39, 0x6e, 0x23, 0x03, 0x31,
      0x69, 0x19, 0xc6, 0x07, 0xb4, 0xfc, 0x72, 0xa3, 0x57, 0x97, 0x44, 0x38, 0xec, 0x04, 0xae,
      0x35, 0x2d, 0x1d, 0x3d, 0x93, 0xb6, 0x78, 0xdc, 0x34, 0x64, 0x97, 0xf6, 0x79, 0xaf, 0x2a,
      0x6f, 0x6d, 0x8d, 0x10, 0xa9, 0x26, 0x8c, 0x24, 0x77, 0xcc, 0xb1, 0x3e, 0xaa, 0x74, 0x0b,
      0x74, 0x0b, 0x4f, 0xe3, 0x8d, 0xe9, 0x2a, 0x01, 0x7e, 0x7e, 0xd0, 0xad, 0x58, 0xeb, 0x92,
      0x43, 0x58, 0x44, 0xdb, 0x85, 0xd0, 0xba, 0x4e, 0x49, 0x48, 0x29, 0xf2, 0xc8, 0x1f, 0xa9,
      0xef, 0xf9, 0xbc, 0x33, 0xcc, 0xcc, 0x1c, 0xfb, 0x7d, 0x04, 0x62, 0x05, 0x27, 0x4f, 0x50,
      0xa8, 0xdf, 0x49, 0x67, 0xfc, 0xc1, 0x99, 0x5f, 0x2c, 0xe0, 0x0b, 0x8b, 0x56, 0x91, 0xf1,
      0x01, 0xeb, 0xd6, 0xe9, 0x8d, 0x0a, 0x93, 0x1f, 0xcb, 0x20, 0x88, 0x97, 0xba, 0xdf, 0x52,
      0xb3, 0x9a, 0xca, 0x17, 0x39, 0x25, 0x03, 0x68, 0x58, 0x98, 0x0d, 0xf8, 0x2a, 0xa3, 0x2d,
      0x7f, 0x12, 0x33, 0x88, 0x59, 0x62, 0x7d, 0xdf, 0xcc, 0x69, 0xd9, 0xe0, 0x2f, 0xd8, 0x6c,
      0x51, 0x4a, 0xad, 0x1e, 0xad, 0x04, 0x33, 0x41, 0x69, 0x13, 0xfa, 0x40, 0x65, 0x7f, 0x01,
      0x34, 0xfa, 0xc6, 0x5e, 0x4a, 0x85, 0x9b, 0xc9, 0x90, 0x53, 0x70, 0xa1, 0xf8, 0xe8, 0x2c,
      0xe4, 0x62, 0x96, 0x26, 0xba, 0xb2, 0xc6, 0x04, 0xf0, 0x51, 0xa4, 0x57, 0x79, 0x14, 0x30,
  };
} // namespace prg_tvs

#endif
