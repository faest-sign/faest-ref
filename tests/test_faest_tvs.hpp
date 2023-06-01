#ifndef TEST_FAEST_TVS_HPP
#define TEST_FAEST_TVS_HPP

#include <array>
#include <cstdint>
#include <string>

namespace faest_tvs {

extern const std::string message;

namespace faest_128s_tvs {
    extern const std::array<uint8_t, 16 + 16> packed_sk;
    extern const std::array<uint8_t, 16 + 16> packed_pk;
    extern const std::array<uint8_t, 16> randomness;
    extern const std::array<uint8_t, 5006> signature;
}

namespace faest_128f_tvs {
    extern const std::array<uint8_t, 16 + 16> packed_sk;
    extern const std::array<uint8_t, 16 + 16> packed_pk;
    extern const std::array<uint8_t, 16> randomness;
    extern const std::array<uint8_t, 6336> signature;
}

namespace faest_192s_tvs {
    extern const std::array<uint8_t, 32 + 24> packed_sk;
    extern const std::array<uint8_t, 32 + 32> packed_pk;
    extern const std::array<uint8_t, 24> randomness;
    extern const std::array<uint8_t, 12744> signature;
}

namespace faest_192f_tvs {
    extern const std::array<uint8_t, 32 + 24> packed_sk;
    extern const std::array<uint8_t, 32 + 32> packed_pk;
    extern const std::array<uint8_t, 24> randomness;
    extern const std::array<uint8_t, 16792> signature;
}

namespace faest_256s_tvs {
    extern const std::array<uint8_t, 32 + 32> packed_sk;
    extern const std::array<uint8_t, 32 + 32> packed_pk;
    extern const std::array<uint8_t, 32> randomness;
    extern const std::array<uint8_t, 22100> signature;
}

namespace faest_256f_tvs {
    extern const std::array<uint8_t, 32 + 32> packed_sk;
    extern const std::array<uint8_t, 32 + 32> packed_pk;
    extern const std::array<uint8_t, 32> randomness;
    extern const std::array<uint8_t, 28400> signature;
}

}

#endif // TEST_FAEST_TVS_HPP
