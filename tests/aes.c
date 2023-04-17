#include <string.h>

#include "../aes/aes.c"

int main(void) {
    bf8_t key[16] = {0x2b, 0x7e, 0x15, 0x16,
                    0x28, 0xae, 0xd2, 0xa6,
                    0xab, 0xf7, 0x15, 0x88,
                    0x09, 0xcf, 0x4f, 0x3c};
    bf8_t iv_[16] = {0x32, 0x43, 0xf6, 0xa8,
                    0x88, 0x5a, 0x30, 0x8d,
                    0x31, 0x31, 0x98, 0xa2,
                    0xe0, 0x37, 0x07, 0x34};
    bf8_t expected[16] = {0x39,0x25,0x84,0x1d,
                            0x02,0xdc,0x09,0xfb,
                            0xdc,0x11,0x85,0x97,
                            0x19,0x6a,0x0b,0x32};
    Initialize(key, iv_);
    bf8_t buffer[16];
    Encrypt(buffer);

    return (memcmp(expected, buffer, sizeof(expected)) == 0) ? 0 : 1;
}

// #define BOOST_TEST_MODULE aes_test
// #include <boost/format.hpp>
// #include <boost/test/included/unit_test.hpp>

// namespace {
//     void check_aes_128(bf8_t *key, bf8_t *iv, bf8_t *expected) {
//         Initialize(key, iv);
//         bf8_t buffer[16];
//         Encrypt(buffer);
//         BOOST_TEST(buffer == expected);
//     }
// }

// BOOST_AUTO_TEST_CASE (test_aes_128) {
//     bf8_t key[16] = {0x2b, 0x7e, 0x15, 0x16,
//                     0x28, 0xae, 0xd2, 0xa6,
//                     0xab, 0xf7, 0x15, 0x88,
//                     0x09, 0xcf, 0x4f, 0x3c};
//     bf8_t iv[16] = {0x32, 0x43, 0xf6, 0xa8,
//                     0x88, 0x5a, 0x30, 0x8d,
//                     0x31, 0x31, 0x98, 0xa2,
//                     0xe0, 0x37, 0x07, 0x34};
//     bf8_t expected[16] = {0x39,0x25,0x84,0x1d,
//                             0x02,0xdc,0x09,0xfb,
//                             0xdc,0x11,0x85,0x97,
//                             0x19,0x6a,0x0b,0x32};
//     check_aes_128(key, iv, expected);
// }