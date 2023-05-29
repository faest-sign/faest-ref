/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef INSTANCES_HPP
#define INSTANCES_HPP

#include "instances.h"
#include "parameters.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <cmath>

constexpr faest_paramid_t all_parameters[] = {
    FAEST_128S,    FAEST_128F,    FAEST_192S,    FAEST_192F,
    FAEST_256S,    FAEST_256F,    FAEST_EM_128S, FAEST_EM_128F,
#if defined(HAVE_EM_192)
    FAEST_EM_192S, FAEST_EM_192F,
#endif
#if defined(HAVE_EM_256)
    FAEST_EM_256S, FAEST_EM_256F,
#endif
};

#endif