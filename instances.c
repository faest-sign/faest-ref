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

#define FAEST_128S_PARAMS                                                                          \
  {                                                                                                \
    FAEST_128S_LAMBDA, FAEST_128S_Nwd, FAEST_128S_Ske, FAEST_128S_R, FAEST_128S_Senc,              \
        FAEST_128S_L, FAEST_128S_Lke, FAEST_128S_Lenc, FAEST_128S_TAU, FAEST_128S_K0,              \
        FAEST_128S_K1, FAEST_128S_T0, FAEST_128S_T1, FAEST_128S_SIG_SIZE, FAEST_128S_PK_SIZE       \
  }
#define FAEST_128F_PARAMS                                                                          \
  {                                                                                                \
    FAEST_128F_LAMBDA, FAEST_128F_Nwd, FAEST_128F_Ske, FAEST_128F_R, FAEST_128F_Senc,              \
        FAEST_128F_L, FAEST_128F_Lke, FAEST_128F_Lenc, FAEST_128F_TAU, FAEST_128F_K0,              \
        FAEST_128F_K1, FAEST_128F_T0, FAEST_128F_T1, FAEST_128F_SIG_SIZE, FAEST_128F_PK_SIZE       \
  }
#define FAEST_192S_PARAMS                                                                          \
  {                                                                                                \
    FAEST_192S_LAMBDA, FAEST_192S_Nwd, FAEST_192S_Ske, FAEST_192S_R, FAEST_192S_Senc,              \
        FAEST_192S_L, FAEST_192S_Lke, FAEST_192S_Lenc, FAEST_192S_TAU, FAEST_192S_K0,              \
        FAEST_192S_K1, FAEST_192S_T0, FAEST_192S_T1, FAEST_192S_SIG_SIZE, FAEST_192S_PK_SIZE       \
  }
#define FAEST_192F_PARAMS                                                                          \
  {                                                                                                \
    FAEST_192F_LAMBDA, FAEST_192F_Nwd, FAEST_192F_Ske, FAEST_192F_R, FAEST_192F_Senc,              \
        FAEST_192F_L, FAEST_192F_Lke, FAEST_192F_Lenc, FAEST_192F_TAU, FAEST_192F_K0,              \
        FAEST_192F_K1, FAEST_192F_T0, FAEST_192F_T1, FAEST_192F_SIG_SIZE, FAEST_192F_PK_SIZE       \
  }
#define FAEST_256S_PARAMS                                                                          \
  {                                                                                                \
    FAEST_256S_LAMBDA, FAEST_256S_Nwd, FAEST_256S_Ske, FAEST_256S_R, FAEST_256S_Senc,              \
        FAEST_256S_L, FAEST_256S_Lke, FAEST_256S_Lenc, FAEST_256S_TAU, FAEST_256S_K0,              \
        FAEST_256S_K1, FAEST_256S_T0, FAEST_256S_T1, FAEST_256S_SIG_SIZE, FAEST_256S_PK_SIZE       \
  }
#define FAEST_256F_PARAMS                                                                          \
  {                                                                                                \
    FAEST_256F_LAMBDA, FAEST_256F_Nwd, FAEST_256F_Ske, FAEST_256F_R, FAEST_256F_Senc,              \
        FAEST_256F_L, FAEST_256F_Lke, FAEST_256F_Lenc, FAEST_256F_TAU, FAEST_256F_K0,              \
        FAEST_256F_K1, FAEST_256F_T0, FAEST_256F_T1, FAEST_256F_SIG_SIZE, FAEST_256F_PK_SIZE       \
  }
#define FAEST_EM_128S_PARAMS                                                                       \
  {                                                                                                \
    FAEST_EM_128S_LAMBDA, FAEST_EM_128S_Nwd, FAEST_EM_128S_Ske, FAEST_EM_128S_R,                   \
        FAEST_EM_128S_Senc, FAEST_EM_128S_L, FAEST_EM_128S_Lke, FAEST_EM_128S_Lenc,                \
        FAEST_EM_128S_TAU, FAEST_EM_128S_K0, FAEST_EM_128S_K1, FAEST_EM_128S_T0, FAEST_EM_128S_T1, \
        FAEST_EM_128S_SIG_SIZE, FAEST_EM_128S_PK_SIZE                                              \
  }
#define FAEST_EM_128F_PARAMS                                                                       \
  {                                                                                                \
    FAEST_EM_128F_LAMBDA, FAEST_EM_128F_Nwd, FAEST_EM_128F_Ske, FAEST_EM_128F_R,                   \
        FAEST_EM_128F_Senc, FAEST_EM_128F_L, FAEST_EM_128F_Lke, FAEST_EM_128F_Lenc,                \
        FAEST_EM_128F_TAU, FAEST_EM_128F_K0, FAEST_EM_128F_K1, FAEST_EM_128F_T0, FAEST_EM_128F_T1, \
        FAEST_EM_128F_SIG_SIZE, FAEST_EM_128F_PK_SIZE                                              \
  }
#define FAEST_EM_192S_PARAMS                                                                       \
  {                                                                                                \
    FAEST_EM_192S_LAMBDA, FAEST_EM_192S_Nwd, FAEST_EM_192S_Ske, FAEST_EM_192S_R,                   \
        FAEST_EM_192S_Senc, FAEST_EM_192S_L, FAEST_EM_192S_Lke, FAEST_EM_192S_Lenc,                \
        FAEST_EM_192S_TAU, FAEST_EM_192S_K0, FAEST_EM_192S_K1, FAEST_EM_192S_T0, FAEST_EM_192S_T1, \
        FAEST_EM_192S_SIG_SIZE, FAEST_EM_192S_PK_SIZE                                              \
  }
#define FAEST_EM_192F_PARAMS                                                                       \
  {                                                                                                \
    FAEST_EM_192F_LAMBDA, FAEST_EM_192F_Nwd, FAEST_EM_192F_Ske, FAEST_EM_192F_R,                   \
        FAEST_EM_192F_Senc, FAEST_EM_192F_L, FAEST_EM_192F_Lke, FAEST_EM_192F_Lenc,                \
        FAEST_EM_192F_TAU, FAEST_EM_192F_K0, FAEST_EM_192F_K1, FAEST_EM_192F_T0, FAEST_EM_192F_T1, \
        FAEST_EM_192F_SIG_SIZE, FAEST_EM_192F_PK_SIZE                                              \
  }
#define FAEST_EM_256S_PARAMS                                                                       \
  {                                                                                                \
    FAEST_EM_256S_LAMBDA, FAEST_EM_256S_Nwd, FAEST_EM_256S_Ske, FAEST_EM_256S_R,                   \
        FAEST_EM_256S_Senc, FAEST_EM_256S_L, FAEST_EM_256S_Lke, FAEST_EM_256S_Lenc,                \
        FAEST_EM_256S_TAU, FAEST_EM_256S_K0, FAEST_EM_256S_K1, FAEST_EM_256S_T0, FAEST_EM_256S_T1, \
        FAEST_EM_256S_SIG_SIZE, FAEST_EM_256S_PK_SIZE                                              \
  }
#define FAEST_EM_256F_PARAMS                                                                       \
  {                                                                                                \
    FAEST_EM_256F_LAMBDA, FAEST_EM_256F_Nwd, FAEST_EM_256F_Ske, FAEST_EM_256F_R,                   \
        FAEST_EM_256F_Senc, FAEST_EM_256F_L, FAEST_EM_256F_Lke, FAEST_EM_256F_Lenc,                \
        FAEST_EM_256F_TAU, FAEST_EM_256F_K0, FAEST_EM_256F_K1, FAEST_EM_256F_T0, FAEST_EM_256F_T1, \
        FAEST_EM_256F_SIG_SIZE, FAEST_EM_256F_PK_SIZE                                              \
  }
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
