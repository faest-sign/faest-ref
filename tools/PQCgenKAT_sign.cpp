
#include <algorithm>
#include <iomanip>
#include <string>
#include <string_view>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <boost/io/ios_state.hpp>

#include "rng.h"
extern "C" {
#include "api.h"
}

#define KAT_SUCCESS 0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_DATA_ERROR -3
#define KAT_CRYPTO_FAILURE -4

namespace {
  constexpr std::string_view AlgName{CRYPTO_ALGNAME};

  template <typename Iterator>
  void print_hex(std::ostream& os, Iterator begin, Iterator end) {
    boost::io::ios_flags_saver ifs{os};
    os << std::hex << std::setfill('0');
    std::for_each(begin, end, [&](unsigned int c) { os << std::setw(2) << c; });
  }
} // namespace

int main() {
  std::ostringstream oss1;
  oss1 << "PQCsignKAT_" << AlgName << ".req";
  const auto fn_req = oss1.str();

  std::ostringstream oss2;
  oss2 << "PQCsignKAT_" << AlgName << ".rsp";
  const auto fn_rsp = oss2.str();

  std::ofstream fp_req{fn_req};
  if (!fp_req) {
    std::cout << "Couldn't open <" << fn_req << "> for write" << std::endl;
    return KAT_FILE_OPEN_ERROR;
  }
  std::ofstream fp_rsp{fn_rsp};
  if (!fp_rsp) {
    std::cout << "Couldn't open <" << fn_rsp << "> for write" << std::endl;
    return KAT_FILE_OPEN_ERROR;
  }

  {
    std::array<uint8_t, 48> entropy_input;
    for (size_t i = 0; i < sizeof(entropy_input); i++)
      entropy_input[i] = i;
    randombytes_init(entropy_input.data(), NULL, 256);
  }

  std::vector<std::array<uint8_t, 48>> seeds;
  std::vector<std::vector<uint8_t>> messages;
  seeds.reserve(100);
  messages.reserve(100);

  for (unsigned int i = 0; i < 100; i++) {
    std::array<uint8_t, 48> seed;
    randombytes(seed.data(), 48);
    seeds.emplace_back(seed);

    fp_req << "count = " << i << "\n";
    fp_req << "seed = ";
    print_hex(fp_req, seed.begin(), seed.end());
    fp_req << "\n";

    const size_t mlen = 33 * (i + 1);
    fp_req << "mlen = " << mlen << "\n";

    std::vector<uint8_t> msg;
    msg.resize(mlen);
    randombytes(msg.data(), mlen);
    fp_req << "msg = ";
    print_hex(fp_req, msg.begin(), msg.end());
    fp_req << "\n";
    messages.emplace_back(msg);

    fp_req << "pk =\nsk = \nsmlen =\nsm=\n" << std::endl;
  }

  fp_rsp << "# " << CRYPTO_ALGNAME << "\n\n";
  for (unsigned int i = 0; i != seeds.size(); ++i) {
    auto& seed      = seeds[i];
    const auto& msg = messages[i];

    fp_rsp << "count = " << i << "\n";
    fp_rsp << "seed = ";
    print_hex(fp_rsp, seed.begin(), seed.end());
    fp_rsp << "\n";
    randombytes_init(seed.data(), NULL, 256);

    fp_rsp << "mlen = " << msg.size() << "\n";
    fp_rsp << "msg = ";
    print_hex(fp_rsp, msg.begin(), msg.end());
    fp_rsp << "\n";

    std::array<uint8_t, CRYPTO_PUBLICKEYBYTES> pk;
    std::array<uint8_t, CRYPTO_SECRETKEYBYTES> sk;

    // Generate the public/private keypair
    auto ret = crypto_sign_keypair(pk.data(), sk.data());
    if (ret != 0) {
      std::cout << "crypto_sign_keypair returned <" << ret << ">" << std::endl;
      return KAT_CRYPTO_FAILURE;
    }
    fp_rsp << "pk = ";
    print_hex(fp_rsp, pk.begin(), pk.end());
    fp_rsp << "\n";
    fp_rsp << "sk = ";
    print_hex(fp_rsp, sk.begin(), sk.end());
    fp_rsp << "\n";

    std::vector<uint8_t> sm;
    sm.resize(msg.size() + CRYPTO_BYTES);
    unsigned long long smlen = sm.size();

    ret = crypto_sign(sm.data(), &smlen, msg.data(), msg.size(), sk.data());
    if (ret != 0) {
      std::cout << "crypto_sign returned <" << ret << ">" << std::endl;
      return KAT_CRYPTO_FAILURE;
    }
    sm.resize(smlen);
    fp_rsp << "smlen = " << smlen << "\n";
    fp_rsp << "sm = ";
    print_hex(fp_rsp, sm.begin(), sm.end());
    fp_rsp << "\n";

    std::vector<uint8_t> msg1;
    msg1.resize(msg.size());
    unsigned long long msg1len = msg1.size();

    ret = crypto_sign_open(msg1.data(), &msg1len, sm.data(), smlen, pk.data());
    if (ret != 0) {
      std::cout << "crypto_sign_open returned <" << ret << ">" << std::endl;
      return KAT_CRYPTO_FAILURE;
    }

    if (msg.size() != msg1len) {
      std::cout << "crypto_sign_open returned bad 'mlen': Got <" << msg1len << ">, expected <"
                << msg.size() << ">" << std::endl;
      return KAT_CRYPTO_FAILURE;
    }

    if (msg != msg1) {
      std::cout << "crypto_sign_open returned bad 'm' value" << std::endl;
      return KAT_CRYPTO_FAILURE;
    }
  }

  return KAT_SUCCESS;
}
