#include "sign.h"

void sign(const uint8_t* msg, const uint8_t* sk, const uint8_t* pk, const faest_paramset_t* params,
          uint32_t voleOutLen) {

  uint32_t lambda      = params->faest_param.lambda;
  uint32_t lambdaBytes = params->faest_param.lambdaBytes;
  uint32_t tau         = params->faest_param.t;

  // Step: 1
  uint8_t* p = malloc(lambdaBytes);
  // TODO: Setting it to all zero for now, ask in the group and change it accordingly
  memset(p, 0, lambdaBytes);

  // Step: 2
  uint8_t* u             = malloc(lambdaBytes * 2);
  uint8_t* pk_msg_concat = malloc(sizeof(pk) + sizeof(msg));
  memcpy(pk_msg_concat, pk, sizeof(pk));
  memcpy(pk_msg_concat + sizeof(pk), msg, sizeof(msg));
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  H1_update(&h1_ctx, pk_msg_concat, sizeof(pk_msg_concat));
  H1_final(&h1_ctx, u, lambdaBytes * 2);

  uint8_t* rootkey       = malloc(lambdaBytes);
  uint8_t* sk_u_p_concat = malloc(sizeof(sk) + sizeof(u) + sizeof(p));
  memcpy(sk_u_p_concat, sk, sizeof(sk));
  memcpy(sk_u_p_concat + sizeof(sk), u, sizeof(u));
  memcpy(sk_u_p_concat + sizeof(sk) + sizeof(u), p, sizeof(p));
  H3_context_t h3_ctx;
  H3_init(&h3_ctx, lambda);
  H3_update(&h3_ctx, sk_u_p_concat, sizeof(sk_u_p_concat));
  H3_final(&h3_ctx, rootkey, lambdaBytes);

  // Step: 3..4
  uint8_t* hcom      = malloc(lambdaBytes * 2);
  vec_com_t** vecCom = malloc(tau * sizeof(vec_com_t*));
  uint8_t** c        = malloc((tau + sizeof(uint8_t*)) - 1);
  uint8_t* u         = malloc(voleOutLen);
  uint8_t** v        = malloc(tau * sizeof(uint8_t*));
  voleCommit(rootkey, lambda, lambdaBytes, voleOutLen, params, hcom, vecCom, c, u, v);

  // Step: 5
  uint8_t* chal_1 = malloc(lambdaBytes);
  H2_context_t h2_ctx;
  H2_init(&h2_ctx, lambda);
  uint8_t* u_hcom_c_concat = malloc(sizeof(u) + sizeof(hcom) + (sizeof(c) * sizeof(c[0])));
  memcpy(u_hcom_c_concat, u, sizeof(u));
  memcpy(u_hcom_c_concat + sizeof(u), hcom, sizeof(hcom));
  for (uint32_t i = 0; i < (tau - 1); i++) {
    memcpy(u_hcom_c_concat + sizeof(u) + sizeof(hcom) + (voleOutLen * i), c[i], voleOutLen);
  }
  H2_update(&h2_ctx, u_hcom_c_concat, sizeof(u_hcom_c_concat));
  H2_final(&h2_ctx, chal_1, (5 * lambdaBytes) + 64);

  // Step: 7
  // TODO: Find what are these ???
  uint8_t* u_tilde;
  uint8_t* r0;
  uint8_t* r1;
  uint8_t* s;
  uint8_t* t;
  uint8_t* x;
  size_t ell;
  vole_hash(u_tilde, r0, r1, s, t, x, ell, lambda);

  // Step: 8
  uint8_t* V_tilde;
  uint8_t* r0_1;
  uint8_t* r1_1;
  uint8_t* s_1;
  uint8_t* t_1;
  uint8_t* x_1;
  size_t ell_1;
  vole_hash(V_tilde, r0_1, r1_1, s_1, t_1, x_1, ell_1, lambda);

  // Step: 9
  uint8_t* h_v = malloc(lambdaBytes * 2);
  H1_context_t h1_ctx_1;
  H1_init(&h1_ctx_1, lambda);
  H1_update(&h1_ctx_1, V_tilde, sizeof(V_tilde));
  H1_final(&h1_ctx_1, h_v, lambdaBytes * 2);

  // Step: 10..11
  // TODO: From where does in, out and B come ?
  uint8_t B;
  uint8_t* in;
  uint8_t* out;
  if (B == 1) {
    // TODO: unclear what they mean by parse here ?
  }
  // Step: 12..13
  // TODO
  uint8_t* w;
  aesExtendWitness(lambda, sk, in, w);
  uint8_t* d;
  uint32_t l;
  xorUint8Arr(w, u, d, l);

  // Step: 14
  // TODO
  uint8_t* chal_2;
  H2_context_t h2_ctx_1;
  uint8_t* chal_1_u_tilde_h_v_d_concat;
  H2_init(&h2_ctx_1, lambda);
  H2_update(&h2_ctx_1, chal_1_u_tilde_h_v_d_concat, sizeof(chal_1_u_tilde_h_v_d_concat));
  H2_final(&h2_ctx_1, chal_2, (3 * lambdaBytes) + 64);

  // Step: 15..17
}
