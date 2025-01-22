#include "bavc_tvs.hpp"

namespace bavc_tvs {
  namespace FAEST_EM_192F {
    const std::array<uint8_t, 48> h{
        0x88, 0x27, 0x90, 0xbf, 0x29, 0x50, 0x80, 0xc6, 0x2f, 0xa3, 0x5e, 0x8d,
        0x24, 0xd1, 0x9d, 0x48, 0xcc, 0x70, 0xe8, 0x38, 0xbe, 0x63, 0xf6, 0x51,
        0x49, 0xc1, 0xf3, 0x82, 0x3c, 0xfb, 0x3e, 0x2b, 0x9f, 0xa5, 0x29, 0xa5,
        0x54, 0x6b, 0xd8, 0xd9, 0x78, 0x4b, 0xa1, 0x74, 0xab, 0x23, 0x00, 0x60,
    };
    const std::array<uint8_t, 64> hashed_k{
        0xa6, 0x1b, 0x98, 0xd0, 0x1d, 0x03, 0x98, 0x71, 0x74, 0xf7, 0x8c, 0x35, 0x5d,
        0xfa, 0xf4, 0xe4, 0xe4, 0x8b, 0xbc, 0x04, 0x4b, 0x27, 0x32, 0x70, 0x30, 0x06,
        0xe6, 0x27, 0x76, 0xbc, 0xff, 0x71, 0xae, 0x77, 0xbc, 0x8f, 0x5f, 0x5c, 0x16,
        0xb8, 0xc9, 0x71, 0xc3, 0xa0, 0xdb, 0xa8, 0x8d, 0xa5, 0x37, 0x76, 0x14, 0x76,
        0x29, 0x67, 0xa8, 0x9c, 0xf4, 0xa2, 0x81, 0x71, 0xf4, 0x2e, 0xe1, 0x14,
    };
    const std::array<uint8_t, 64> hashed_sd{
        0x01, 0x97, 0x6e, 0xfa, 0xe0, 0x9f, 0x17, 0xd0, 0xad, 0xfb, 0x23, 0xfb, 0xb8,
        0x58, 0xbc, 0x56, 0x9d, 0x7a, 0x11, 0xb9, 0x45, 0xd1, 0x8d, 0x9f, 0x0c, 0x8e,
        0xda, 0x00, 0x42, 0xbb, 0xf5, 0x46, 0x59, 0x9c, 0x62, 0x86, 0x58, 0x52, 0x71,
        0xce, 0xea, 0x43, 0x81, 0xfb, 0x5a, 0xcf, 0xa3, 0x5b, 0x09, 0x39, 0x19, 0x48,
        0xcd, 0x2f, 0x41, 0x51, 0x06, 0x89, 0xd4, 0xee, 0x56, 0x2c, 0x72, 0x95,
    };
    const std::array<uint16_t, 24> i_delta{
        0x00c9, 0x0040, 0x00b5, 0x00f2, 0x0004, 0x0067, 0x0040, 0x0005,
        0x0085, 0x0058, 0x0046, 0x008f, 0x0023, 0x008b, 0x0085, 0x00db,
        0x003f, 0x0035, 0x005f, 0x001f, 0x001e, 0x0028, 0x0074, 0x0015,
    };
    const std::array<uint8_t, 64> hashed_decom_i{
        0xcf, 0xa7, 0xc1, 0x1b, 0x0d, 0xcd, 0x07, 0x71, 0xf7, 0xb5, 0x39, 0x55, 0x38,
        0x22, 0x04, 0x3b, 0xc0, 0xad, 0x11, 0x61, 0x42, 0x9a, 0x66, 0xec, 0xdc, 0x80,
        0x68, 0x40, 0x3a, 0xdf, 0x86, 0xad, 0x47, 0xe7, 0xf8, 0x0d, 0x1a, 0x61, 0x79,
        0xdc, 0x7c, 0xf0, 0xe9, 0xc2, 0xaa, 0x48, 0xae, 0xcd, 0x99, 0x86, 0xf7, 0x00,
        0xc7, 0x61, 0x81, 0xa4, 0x2c, 0x7c, 0xa5, 0x06, 0x0d, 0x6d, 0x48, 0xa7,
    };
    const std::array<uint8_t, 64> hashed_rec_sd{
        0xca, 0x6a, 0x7f, 0xd6, 0x1f, 0x56, 0xe7, 0x8f, 0xbf, 0xdb, 0xf4, 0xb5, 0xcd,
        0x24, 0xc7, 0x16, 0x52, 0xa4, 0x60, 0x4c, 0x8f, 0xb2, 0x45, 0x8a, 0x4d, 0xaf,
        0x8b, 0x2b, 0xa1, 0xa8, 0xa1, 0x1b, 0xcb, 0xb5, 0xf0, 0x66, 0x01, 0xb8, 0xbd,
        0xa3, 0xef, 0x11, 0x31, 0xf0, 0x43, 0x3f, 0xe0, 0xa9, 0xd3, 0x71, 0xd4, 0x80,
        0x64, 0x60, 0x7b, 0x8d, 0x13, 0x44, 0xee, 0x60, 0xc2, 0xeb, 0x0c, 0x26,
    };
  } // namespace FAEST_EM_192F
} // namespace bavc_tvs
