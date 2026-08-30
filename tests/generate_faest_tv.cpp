/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest.h"
#include "utils.hpp"

#include <array>
#include <vector>

namespace {
  constexpr std::array<uint8_t, 16> random_seed{
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
      0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
  };

  // "This document describes and specifies the FAEST digital signature algorithm.";
  constexpr std::array<uint8_t, 76> message{
      0x54, 0x68, 0x69, 0x73, 0x20, 0x64, 0x6f, 0x63, 0x75, 0x6d, 0x65, 0x6e, 0x74,
      0x20, 0x64, 0x65, 0x73, 0x63, 0x72, 0x69, 0x62, 0x65, 0x73, 0x20, 0x61, 0x6e,
      0x64, 0x20, 0x73, 0x70, 0x65, 0x63, 0x69, 0x66, 0x69, 0x65, 0x73, 0x20, 0x74,
      0x68, 0x65, 0x20, 0x46, 0x41, 0x45, 0x53, 0x54, 0x20, 0x64, 0x69, 0x67, 0x69,
      0x74, 0x61, 0x6c, 0x20, 0x73, 0x69, 0x67, 0x6e, 0x61, 0x74, 0x75, 0x72, 0x65,
      0x20, 0x61, 0x6c, 0x67, 0x6f, 0x72, 0x69, 0x74, 0x68, 0x6d, 0x2e,
  };
} // namespace

void generate_faest_128f_tv(const uint8_t* message, size_t message_length,
                            const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_128s_tv(const uint8_t* message, size_t message_length,
                            const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_192f_tv(const uint8_t* message, size_t message_length,
                            const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_192s_tv(const uint8_t* message, size_t message_length,
                            const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_256f_tv(const uint8_t* message, size_t message_length,
                            const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_256s_tv(const uint8_t* message, size_t message_length,
                            const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_em_128f_tv(const uint8_t* message, size_t message_length,
                               const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_em_128s_tv(const uint8_t* message, size_t message_length,
                               const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_em_192f_tv(const uint8_t* message, size_t message_length,
                               const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_em_192s_tv(const uint8_t* message, size_t message_length,
                               const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_em_256f_tv(const uint8_t* message, size_t message_length,
                               const uint8_t* random_seed, size_t random_seed_length);
void generate_faest_em_256s_tv(const uint8_t* message, size_t message_length,
                               const uint8_t* random_seed, size_t random_seed_length);

int main() {
  std::cout << "#ifndef TEST_FAEST_TVS_HPP\n";
  std::cout << "#define TEST_FAEST_TVS_HPP\n\n";
  std::cout << "#include <array>\n";
  std::cout << "#include <cstdint>\n\n";
  std::cout << "namespace faest_tvs {\n";
  print_named_array("random_seed", "uint8_t", random_seed);
  print_named_array("message", "uint8_t", message);
  std::cout << "\n";

  generate_faest_128f_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_128s_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_192f_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_192s_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_256f_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_256s_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_em_128f_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_em_128s_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_em_192f_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_em_192s_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_em_256f_tv(message.data(), message.size(), random_seed.data(), random_seed.size());
  generate_faest_em_256s_tv(message.data(), message.size(), random_seed.data(), random_seed.size());

  std::cout << "}\n\n#endif" << std::endl;
}