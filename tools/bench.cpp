/*
 *  SPDX-License-Identifier: MIT
 */

extern "C" {
#include "api.h"
}

#include <array>
#include <boost/program_options.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

namespace {
  struct timing_and_size_t {
    microseconds keygen, sign, verify;
  };

  unsigned int parse_args(int argc, char** argv) {
    using namespace boost::program_options;

    options_description options{"Options"};
    options.add_options()("help", "produce help message");
    options.add_options()("iter,i", value<unsigned int>()->default_value(100),
                          "set number of iterations");

    variables_map vm;
    try {
      store(parse_command_line(argc, argv, options), vm);
      notify(vm);

      if (vm.count("help")) {
        std::cout << options << std::endl;
        return 0;
      }

      return vm["iter"].as<unsigned int>();
    } catch (const boost::exception& e) {
      std::cout << options << std::endl;
      return 0;
    }
  }

  void print_timings(const std::vector<timing_and_size_t>& timings) {
    for (const auto& timing : timings) {
      std::cout << timing.keygen.count() << ',' << timing.sign.count() << ','
                << timing.verify.count() << std::endl;
    }
  }

  void bench_sign_and_verify(unsigned int iter) {
    std::vector<timing_and_size_t> timings;
    timings.reserve(iter);

    std::array<uint8_t, 32> msg;
    {
      std::uniform_int_distribution<unsigned int> dist{0, 255};
      std::random_device rnd;
      std::default_random_engine eng(rnd());
      std::generate(msg.begin(), msg.end(), [&dist, &eng] { return dist(eng); });
    }

    std::vector<uint8_t> sig;
    sig.resize(msg.size() + CRYPTO_BYTES);

    for (unsigned int i = 0; i != iter; ++i) {
      timing_and_size_t timing;

      std::array<uint8_t, CRYPTO_PUBLICKEYBYTES> pk;
      std::array<uint8_t, CRYPTO_SECRETKEYBYTES> sk;

      // Generate the public/private keypair
      auto start_time = high_resolution_clock::now();
      auto ret        = crypto_sign_keypair(pk.data(), sk.data());
      timing.keygen   = duration_cast<microseconds>(high_resolution_clock::now() - start_time);
      if (ret != 0) {
        std::cout << "crypto_sign_keypair failed: <" << ret << std::endl;
        break;
      }

      unsigned long long smlen = sig.size();

      start_time  = high_resolution_clock::now();
      ret         = crypto_sign(sig.data(), &smlen, msg.data(), msg.size(), sk.data());
      timing.sign = duration_cast<microseconds>(high_resolution_clock::now() - start_time);
      if (ret != 0) {
        std::cout << "crypto_sign failed: " << ret << std::endl;
        break;
      }

      std::array<uint8_t, 32> msg1;
      unsigned long long msg1len = msg1.size();

      start_time    = high_resolution_clock::now();
      ret           = crypto_sign_open(msg1.data(), &msg1len, sig.data(), smlen, pk.data());
      timing.verify = duration_cast<microseconds>(high_resolution_clock::now() - start_time);
      if (ret != 0) {
        std::cout << "crypto_sign_open failed:" << ret << std::endl;
        break;
      }

      timings.emplace_back(timing);
    }

    print_timings(timings);
  }
} // namespace

int main(int argc, char** argv) {
  unsigned int iter = parse_args(argc, argv);
  if (!iter) {
    return 1;
  }

  bench_sign_and_verify(iter);
  return 0;
}
