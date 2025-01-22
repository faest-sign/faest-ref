#include "bavc_tvs.hpp"

namespace bavc_tvs {
  namespace FAEST_256F {
    const std::array<uint8_t, 64> h{
        0x3c, 0xc4, 0x6e, 0x66, 0x2d, 0xc1, 0x34, 0x91, 0xf5, 0x04, 0xb0, 0xd4, 0x53,
        0x1a, 0xa3, 0x28, 0x37, 0xf3, 0x37, 0xfb, 0x1a, 0x6f, 0x3f, 0x00, 0x36, 0xb2,
        0x34, 0x8f, 0xdd, 0x45, 0x0f, 0x1d, 0xf7, 0x2c, 0xc2, 0x8e, 0x2a, 0x61, 0xae,
        0x52, 0x99, 0xed, 0x09, 0x50, 0xe6, 0xbb, 0xd6, 0xa7, 0xe6, 0xdd, 0x3c, 0x52,
        0x7a, 0x09, 0xf2, 0x91, 0x89, 0xe1, 0x8f, 0xc9, 0xe1, 0xdd, 0xce, 0x91,
    };
    const std::array<uint8_t, 64> hashed_k{
        0x74, 0xfd, 0x6a, 0x67, 0x37, 0xc4, 0x4d, 0x6f, 0x0e, 0x08, 0x8c, 0xe9, 0x57,
        0x6b, 0x74, 0x1f, 0x24, 0x9b, 0x52, 0xf1, 0x41, 0x02, 0x23, 0xd7, 0xc0, 0x6b,
        0xbf, 0x12, 0xd1, 0x5b, 0xaa, 0x03, 0xe1, 0x4d, 0x40, 0x23, 0xe6, 0x50, 0x79,
        0x8d, 0x21, 0xe5, 0x73, 0x9d, 0x16, 0xb1, 0xf8, 0x82, 0x6b, 0xb2, 0x17, 0x2f,
        0x21, 0x41, 0xae, 0x7b, 0xa6, 0x70, 0x2a, 0xd0, 0xfd, 0x51, 0xdf, 0x6a,
    };
    const std::array<uint8_t, 64> hashed_sd{
        0xa6, 0x08, 0xa7, 0xea, 0x9b, 0x73, 0x60, 0xf9, 0xac, 0xbd, 0xbd, 0x54, 0xd2,
        0x6f, 0x0e, 0x9b, 0xc2, 0xd2, 0xb2, 0xa0, 0x0b, 0x9a, 0x88, 0x28, 0xc1, 0xf8,
        0xe4, 0x7c, 0x1a, 0xd6, 0xb4, 0x3f, 0x4a, 0x1b, 0xeb, 0x8e, 0x7d, 0x65, 0x15,
        0x17, 0xf8, 0xab, 0x24, 0xfa, 0x34, 0xde, 0xf3, 0xcb, 0x4e, 0x57, 0xf4, 0x4e,
        0x18, 0xfa, 0x40, 0x38, 0x7c, 0x23, 0xd8, 0xa1, 0x30, 0xc9, 0x42, 0x5c,
    };
    const std::array<uint16_t, 32> i_delta{
        0x00c9, 0x0040, 0x00b5, 0x00f2, 0x0004, 0x0067, 0x0040, 0x0005, 0x0085, 0x0058, 0x0046,
        0x008f, 0x0023, 0x008b, 0x0085, 0x00db, 0x007f, 0x006b, 0x00be, 0x003f, 0x003d, 0x0051,
        0x00e9, 0x002a, 0x001f, 0x0019, 0x005b, 0x007b, 0x0062, 0x000a, 0x003a, 0x0020,
    };
    const std::array<uint8_t, 64> hashed_decom_i{
        0x23, 0xe8, 0x48, 0x9c, 0x3e, 0x60, 0x2a, 0x0d, 0x40, 0xff, 0x79, 0x80, 0xc8,
        0x3b, 0x3f, 0xc5, 0xdd, 0x90, 0x8d, 0x02, 0xc5, 0x32, 0xf0, 0x5b, 0xb0, 0x9d,
        0x09, 0x3c, 0x76, 0xbb, 0xe9, 0x36, 0xc9, 0xe6, 0x39, 0x38, 0x84, 0x0f, 0xd8,
        0x8f, 0x4c, 0x6e, 0x78, 0xff, 0xfd, 0x08, 0xa1, 0x10, 0xa1, 0xff, 0x7a, 0x23,
        0xec, 0x7e, 0x48, 0xe3, 0x43, 0x39, 0x27, 0x35, 0x5a, 0x4d, 0xab, 0x12,
    };
    const std::array<uint8_t, 64> hashed_rec_sd{
        0xe8, 0xdb, 0x9b, 0x38, 0x92, 0xa8, 0x81, 0xcf, 0x9e, 0x13, 0xff, 0x9e, 0xa9,
        0x77, 0xea, 0x03, 0x24, 0xb8, 0x51, 0x9c, 0xe6, 0xaf, 0x3b, 0x6d, 0xea, 0x07,
        0x36, 0x9b, 0xf3, 0x53, 0x41, 0x28, 0x7a, 0x75, 0xfb, 0x47, 0xdc, 0x18, 0x51,
        0x53, 0x1a, 0xcd, 0xae, 0xf1, 0xdc, 0x20, 0xf4, 0xfc, 0xbe, 0xe9, 0x14, 0x6e,
        0x96, 0x80, 0x4d, 0xfc, 0x2d, 0x3f, 0x25, 0xef, 0x83, 0x8d, 0xd6, 0x27,
    };
  } // namespace FAEST_256F
} // namespace bavc_tvs
