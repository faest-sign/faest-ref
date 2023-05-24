#include "../instances.c"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <cmath>

namespace {
  constexpr faest_paramid_t all_parameters[] = {
      FAEST_128S,    FAEST_128F,    FAEST_192S,    FAEST_192F,    FAEST_256S,    FAEST_256F,
      FAEST_EM_128S, FAEST_EM_128F, FAEST_EM_192S, FAEST_EM_192F, FAEST_EM_256S, FAEST_EM_256F,
  };
} // namespace

BOOST_AUTO_TEST_SUITE(instances)

BOOST_DATA_TEST_CASE(test_keys, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    const auto param        = faest_get_paramset(param_id);
    const auto& faest_param = param.faest_param;

    BOOST_TEST(faest_param.lambda ==
               faest_param.k0 * faest_param.t0 + faest_param.k1 * faest_param.t1);
    BOOST_TEST(faest_param.k0 == std::ceil((float)faest_param.lambda / faest_param.tau));
    BOOST_TEST(faest_param.k1 == std::floor((float)faest_param.lambda / faest_param.tau));
    BOOST_TEST(faest_param.t0 == faest_param.lambda % faest_param.tau);
    BOOST_TEST(faest_param.t1 == faest_param.tau - (faest_param.lambda % faest_param.tau));

    BOOST_TEST(
        (faest_param.lambda == 128 || faest_param.lambda == 192 || faest_param.lambda == 256));
  }
}

BOOST_AUTO_TEST_SUITE_END()