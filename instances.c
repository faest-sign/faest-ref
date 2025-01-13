/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "instances.h"
#include "parameters.h"

const char* faest_get_param_name(faest_paramid_t paramid) {
  switch (paramid) {
  case PARAMETER_SET_INVALID:
    return "PARAMETER_SET_INVALID";
  case FAEST_128S:
    return "FAEST_128S";
  case FAEST_128F:
    return "FAEST_128F";
  case FAEST_192S:
    return "FAEST_192S";
  case FAEST_192F:
    return "FAEST_192F";
  case FAEST_256S:
    return "FAEST_256S";
  case FAEST_256F:
    return "FAEST_256F";
  case FAEST_EM_128S:
    return "FAEST_EM_128S";
  case FAEST_EM_128F:
    return "FAEST_EM_128F";
  case FAEST_EM_192S:
    return "FAEST_EM_192S";
  case FAEST_EM_192F:
    return "FAEST_EM_192F";
  case FAEST_EM_256S:
    return "FAEST_EM_256S";
  case FAEST_EM_256F:
    return "FAEST_EM_256F";
  default:
    return "PARAMETER_SET_MAX_INDEX";
  }
}

#define PARAMS(name)                                                                               \
  {                                                                                                \
    name##_LAMBDA, name##_TAU, name##_W_GRIND, name##_T_OPEN, name##_ELL, name##_K, name##_T0,     \
        name##_T1, name##_L, name##_NwD, name##_Ske, name##_R, name##_Senc, name##_Lke,            \
        name##_Lenc, name##_SIG_SIZE, name##_PK_SIZE,                                              \
  }

#define FAEST_128S_PARAMS PARAMS(FAEST_128S)
#define FAEST_128F_PARAMS PARAMS(FAEST_128F)
#define FAEST_192S_PARAMS PARAMS(FAEST_192S)
#define FAEST_192F_PARAMS PARAMS(FAEST_192F)
#define FAEST_256S_PARAMS PARAMS(FAEST_256S)
#define FAEST_256F_PARAMS PARAMS(FAEST_256F)
#define FAEST_EM_128S_PARAMS PARAMS(FAEST_EM_128S)
#define FAEST_EM_128F_PARAMS PARAMS(FAEST_EM_128F)
#define FAEST_EM_192S_PARAMS PARAMS(FAEST_EM_192S)
#define FAEST_EM_192F_PARAMS PARAMS(FAEST_EM_192F)
#define FAEST_EM_256S_PARAMS PARAMS(FAEST_EM_256S)
#define FAEST_EM_256F_PARAMS PARAMS(FAEST_EM_256F)
#define FAEST_INVALID_PARAMS                                                                       \
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }

faest_paramset_t faest_get_paramset(faest_paramid_t paramid) {
  switch (paramid) {
  case FAEST_128S:
    return (faest_paramset_t){FAEST_128S_PARAMS, FAEST_128S};
  case FAEST_128F:
    return (faest_paramset_t){FAEST_128F_PARAMS, FAEST_128F};
  case FAEST_EM_128S:
    return (faest_paramset_t){FAEST_EM_128S_PARAMS, FAEST_EM_128S};
  case FAEST_EM_128F:
    return (faest_paramset_t){FAEST_EM_128F_PARAMS, FAEST_EM_128F};
  case FAEST_192S:
    return (faest_paramset_t){FAEST_192S_PARAMS, FAEST_192S};
  case FAEST_192F:
    return (faest_paramset_t){FAEST_192F_PARAMS, FAEST_192F};
  case FAEST_EM_192S:
    return (faest_paramset_t){FAEST_EM_192S_PARAMS, FAEST_EM_192S};
  case FAEST_EM_192F:
    return (faest_paramset_t){FAEST_EM_192F_PARAMS, FAEST_EM_192F};
  case FAEST_256S:
    return (faest_paramset_t){FAEST_256S_PARAMS, FAEST_256S};
  case FAEST_256F:
    return (faest_paramset_t){FAEST_256F_PARAMS, FAEST_256F};
  case FAEST_EM_256S:
    return (faest_paramset_t){FAEST_EM_256S_PARAMS, FAEST_EM_256S};
  case FAEST_EM_256F:
    return (faest_paramset_t){FAEST_EM_256F_PARAMS, FAEST_EM_256F};
  default:
    return (faest_paramset_t){FAEST_INVALID_PARAMS, PARAMETER_SET_INVALID};
  }
}
