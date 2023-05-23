#include "faest.h"

// TODO: TEST EVERYTHING HERE !!!

// TODO: Do not pass lambdaBytes everywhere, compute it in the function....

// TODO: change q to Q where applicable

// TODO: Using the simple rand(), set to some secure sampling
void keyGen(uint32_t lambda, uint32_t lambdaBytes, uint8_t* sk, uint8_t* pk) {

  uint8_t* x0    = malloc(16);
  uint8_t* out_0 = malloc(16);
  for (uint32_t i = 0; i < 16; i++) {
    x0[i] = rand() % UINT8_MAX;
  }

  uint8_t* x1    = malloc(16);
  uint8_t* out_1 = malloc(16);
  if (lambda == 192 || lambda == 256) {
    for (uint32_t i = 0; i < 16; i++) {
      x1[i] = rand() % UINT8_MAX;
    }
  }

  bool zeroInpSB = false;
  uint8_t* key   = malloc(lambdaBytes);
  uint8_t* y;
  switch (lambda) {
  case 256:
    while (zeroInpSB == false) {
      for (uint32_t i = 0; i < lambdaBytes; i++) {
        key[i] = rand() % UINT8_MAX;
      }
      if (owf_256(key, x0, out_0) == true && owf_256(key, x1, out_1) == true) {
        y = malloc(32);
        memcpy(y, out_0, 16);
        memcpy(y + 16, out_1, 16);
        zeroInpSB = true;
      }
    }
    break;
  case 192:
    while (zeroInpSB == false) {
      for (uint32_t i = 0; i < lambdaBytes; i++) {
        key[i] = rand() % UINT8_MAX;
      }
      if (owf_192(key, x0, out_0) == true && owf_192(key, x1, out_1) == true) {
        y = malloc(32);
        memcpy(y, out_0, 16);
        memcpy(y + 16, out_1, 16);
        zeroInpSB = true;
      }
    }
    break;
  default:
    while (zeroInpSB == false) {
      for (uint32_t i = 0; i < lambdaBytes; i++) {
        key[i] = rand() % UINT8_MAX;
      }
      if (owf_128(key, x0, out_0) == true) {
        y = malloc(16);
        memcpy(y, out_0, 16);
        zeroInpSB = true;
      }
    }
  }
  sk = malloc(lambdaBytes);
  memcpy(sk, key, lambdaBytes);

  if (lambda == 128) {
    pk = malloc(32);
    memcpy(pk, x0, 16);
    memcpy(pk, y, 16);
  } else {
    pk = malloc(64);
    memcpy(pk, x0, 16);
    memcpy(pk + 16, x1, 16);
    memcpy(pk + 32, y, 32);
  }
}

// TODO: l is in bits, change it everywhere
void sign(const uint8_t* msg, const uint8_t* sk, const uint8_t* pk, const faest_paramset_t* params,
          uint32_t l, signature_t* signature) {

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
  vec_com_t** vecCom = malloc(tau * sizeof(vec_com_t*));
  uint8_t* u         = malloc(l);
  uint8_t** v        = malloc(tau * sizeof(uint8_t*));
  voleCommit(rootkey, lambda, lambdaBytes, l, params, signature->hcom, vecCom, signature->c, u, v);

  // Step: 5
  uint8_t* chal_1 = malloc(lambdaBytes);
  H2_context_t h2_ctx;
  H2_init(&h2_ctx, lambda);
  uint8_t* u_hcom_c_concat = malloc(sizeof(u) + sizeof(signature->hcom) + (tau * l));
  memcpy(u_hcom_c_concat, u, sizeof(u));
  memcpy(u_hcom_c_concat + sizeof(u), signature->hcom, sizeof(signature->hcom));
  for (uint32_t i = 0; i < (tau - 1); i++) {
    memcpy(u_hcom_c_concat + sizeof(u) + sizeof(signature->hcom) + (l * i), signature->c[i], l);
  }
  H2_update(&h2_ctx, u_hcom_c_concat, sizeof(u_hcom_c_concat));
  H2_final(&h2_ctx, chal_1, (5 * lambdaBytes) + 8);

  // Step: 7
  // TODO: Find what are these ???
  uint8_t* r0;
  uint8_t* r1;
  uint8_t* s;
  uint8_t* t;
  uint8_t* x;
  size_t ell;
  vole_hash(signature->u_tilde, r0, r1, s, t, x, ell, lambda);

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
  H1_context_t h1_ctx_1;
  H1_init(&h1_ctx_1, lambda);
  H1_update(&h1_ctx_1, V_tilde, sizeof(V_tilde));
  H1_final(&h1_ctx_1, signature->h_v, lambdaBytes * 2);

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
  aes_extend_witness(lambda, sk, in, w);
  uint32_t l;
  xorUint8Arr(w, u, signature->d, l);

  // Step: 14
  // TODO
  uint8_t* chal_2;
  uint8_t* chal_1_u_tilde_h_v_d_concat;
  memcpy(chal_1_u_tilde_h_v_d_concat, chal_1, sizeof(chal_1));

  memcpy(chal_1_u_tilde_h_v_d_concat + sizeof(chal_1), signature->u_tilde,
         sizeof(signature->u_tilde));

  memcpy(chal_1_u_tilde_h_v_d_concat + sizeof(chal_1) + sizeof(signature->u_tilde), signature->h_v,
         sizeof(signature->h_v));

  for (uint32_t i = 0; i < tau; i++) {

    memcpy(chal_1_u_tilde_h_v_d_concat + sizeof(chal_1) + sizeof(signature->u_tilde) +
               sizeof(signature->h_v) + (i * l),
           signature->c[i], sizeof(l));
  }
  H2_context_t h2_ctx_1;
  H2_init(&h2_ctx_1, lambda);
  H2_update(&h2_ctx_1, chal_1_u_tilde_h_v_d_concat, sizeof(chal_1_u_tilde_h_v_d_concat));
  H2_final(&h2_ctx_1, chal_2, (3 * lambdaBytes) + 8);

  // Step: 15..17
  // TODO
  // u = u[0..l+lambda)
  // v = v[0...l+lambda]

  // Step: 18
  // TODO
  aes_prove(w, u, v, in, out, chal_2, lambda, tau, l, signature->a_tilde, signature->b_tilde);

  // Step: 19
  // TODO
  uint8_t* chal_2_a_tilde_b_tilde_concat =
      malloc(sizeof(chal_2) + sizeof(signature->a_tilde) + sizeof(signature->b_tilde));
  memcpy(chal_2_a_tilde_b_tilde_concat, chal_2, sizeof(chal_2));
  memcpy(chal_2_a_tilde_b_tilde_concat + sizeof(chal_2), signature->a_tilde,
         sizeof(signature->a_tilde));
  memcpy(chal_2_a_tilde_b_tilde_concat + sizeof(chal_2) + sizeof(signature->a_tilde),
         signature->b_tilde, sizeof(signature->b_tilde));
  uint8_t* chal_3 = malloc(lambdaBytes);
  H2_context_t h2_ctx_2;
  H2_init(&h2_ctx_2, lambda);
  H2_update(&h2_ctx_2, chal_2_a_tilde_b_tilde_concat, sizeof(chal_2_a_tilde_b_tilde_concat));
  H2_final(&h2_ctx_2, chal_3, lambdaBytes);

  // Step: 20..23
  uint32_t numVoleInstances = 0;
  uint32_t depth            = 0;
  uint8_t** s               = malloc(tau * sizeof(uint8_t*));
  for (uint32_t i = 0; i < tau; i++) {
    ChalDec(chal_3, i, params->faest_param.k0, params->faest_param.t0, params->faest_param.k1,
            params->faest_param.t1, s[i]);

    if (i < (lambda % tau)) { // Computing the num of vole and the depth here
      numVoleInstances = 1 << params->faest_param.k0;
      depth            = params->faest_param.k0;
    } else {
      numVoleInstances = 1 << params->faest_param.k1;
      depth            = params->faest_param.k1;
    }
    signature->pdec[i]  = malloc(depth);
    signature->com_j[i] = malloc(lambdaBytes * 2);
    vector_open(vecCom[i]->k, vecCom[i]->com, s[i], signature->pdec[i], signature->com_j[i],
                numVoleInstances, lambdaBytes);
  }
}

// TODO: l is in bits, change it everywhere
int verify(const uint8_t* msg, const uint8_t* pk, const faest_paramset_t* params, uint32_t l,
           const signature_t* signature) {

  uint32_t lambda      = params->faest_param.lambda;
  uint32_t lambdaBytes = params->faest_param.lambdaBytes;
  uint32_t tau         = params->faest_param.t;

  // Step: 2..3
  // TODO: Where does B come from ?
  uint8_t B;
  uint8_t* in;
  uint8_t* out;
  if (B == 1) {
    // TODO: unclear what they mean by parse here ?
  } else {
    // TODO
  }

  // Step: 4
  if (params->faest_param.t0 != lambda % tau) {
    return -1;
  }
  uint8_t* mu            = malloc(lambdaBytes * 2);
  uint8_t* pk_msg_concat = malloc(sizeof(pk) + sizeof(msg));
  memcpy(pk_msg_concat, pk, sizeof(pk));
  memcpy(pk_msg_concat + sizeof(pk), msg, sizeof(msg));
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  H1_update(&h1_ctx, pk_msg_concat, sizeof(pk_msg_concat));
  H1_final(&h1_ctx, mu, lambdaBytes * 2);

  // Step: 6
  uint8_t* chal_1           = malloc((5 * lambdaBytes) + 8);
  uint8_t* mu_hcom_c_concat = malloc(sizeof(mu) + sizeof(signature->hcom) + (tau * l));
  memcpy(mu_hcom_c_concat, mu, sizeof(mu));

  memcpy(mu_hcom_c_concat + sizeof(mu), signature->hcom, sizeof(signature->hcom));

  for (uint32_t i = 0; i < tau; i++) {
    memcpy(mu_hcom_c_concat + sizeof(mu) + sizeof(signature->hcom) + (i * l), signature->c[i],
           sizeof(l));
  }
  H2_context_t* h2_ctx;
  H2_init(&h2_ctx, lambda);
  H2_update(&h2_ctx, mu_hcom_c_concat, sizeof(mu_hcom_c_concat));
  H2_final(&h2_ctx, chal_1, (5 * lambdaBytes) + 8);

  // Step: 7
  // TODO
  uint8_t* chal_2;
  H2_context_t h2_ctx_1;
  uint8_t* chal_1_u_tilde_h_v_d_concat = malloc(sizeof(chal_1) + sizeof(signature->u_tilde) +
                                                sizeof(signature->h_v) + sizeof(signature->d));
  memcpy(chal_1_u_tilde_h_v_d_concat, chal_1, sizeof(chal_1));
  memcpy(chal_1_u_tilde_h_v_d_concat + sizeof(chal_1), signature->u_tilde,
         sizeof(signature->u_tilde));
  memcpy(chal_1_u_tilde_h_v_d_concat + sizeof(chal_1) + sizeof(signature->u_tilde), signature->h_v,
         sizeof(signature->h_v));
  memcpy(chal_1_u_tilde_h_v_d_concat + sizeof(chal_1) + sizeof(signature->u_tilde) +
             sizeof(signature->h_v),
         signature->d, sizeof(signature->d));
  H2_init(&h2_ctx_1, lambda);
  H2_update(&h2_ctx_1, chal_1_u_tilde_h_v_d_concat, sizeof(chal_1_u_tilde_h_v_d_concat));
  H2_final(&h2_ctx_1, chal_2, (3 * lambdaBytes) + 8);

  // Step: 8
  // TODO :: Check
  uint8_t* chal_2_a_tilde_b_tilde_concat =
      malloc(sizeof(chal_2) + sizeof(signature->a_tilde) + sizeof(signature->b_tilde));
  memcpy(chal_2_a_tilde_b_tilde_concat, chal_2, sizeof(chal_2));
  memcpy(chal_2_a_tilde_b_tilde_concat + sizeof(chal_2), signature->a_tilde,
         sizeof(signature->a_tilde));
  memcpy(chal_2_a_tilde_b_tilde_concat + sizeof(chal_2) + sizeof(signature->a_tilde),
         signature->b_tilde, sizeof(signature->b_tilde));
  uint8_t* chal_3 = malloc(lambdaBytes);
  H2_context_t h2_ctx_2;
  H2_init(&h2_ctx_2, lambda);
  H2_update(&h2_ctx_2, chal_2_a_tilde_b_tilde_concat, sizeof(chal_2_a_tilde_b_tilde_concat));
  H2_final(&h2_ctx_2, chal_3, lambdaBytes);

  // Step: 9..11
  uint8_t** q              = malloc(tau * sizeof(uint8_t*));
  vec_com_rec_t* vecComRec = malloc(sizeof(vec_com_rec_t));
  uint8_t* hcomverify      = malloc(lambdaBytes * 2);
  voleVerify(chal_3, signature->pdec, signature->com_j, lambda, lambdaBytes, l, tau,
             params->faest_param.k0, params->faest_param.k1, hcomverify, q, vecComRec);
  if (memcmp(hcomverify, signature->hcom, lambdaBytes * 2) != 0) {
    return 0;
  }

  // Step: 12
  // TODO: no clue what this will take... FIX Q !!
  uint8_t* r0;
  uint8_t** q_tilde = malloc(tau * sizeof(uint8_t*));
  uint32_t depth    = 0;
  size_t ell;
  q_tilde[0] = malloc(l * depth);
  memcpy(q_tilde[0], q[0], sizeof(q[0]));
  for (uint32_t i = 1; i < tau; i++) {
    if (i < (lambda % tau)) { // Computing the num of vole and the depth here
      depth = params->faest_param.k0;
    } else {
      depth = params->faest_param.k1;
    }
    q_tilde[i] = malloc(l * depth);
    for (uint32_t d = 0; d < depth; d++) {
      for (uint32_t k = 0; k < l; k++) {
        *(q_tilde[i] + (d * l) + k) = *(q[i] + (d * l) + k) ^ *(signature->c[i] + (d * l) + k);
      }
    }
  }
  // TODO: currently guessing the params... Also the params are wrong : ((((
  vole_hash(vecComRec->h, r0, rintf128, chal_1, tau, q_tilde, ell, lambda);

  // Step: 13..16
  uint8_t b         = 0;
  uint8_t** d_tilde = malloc(tau * sizeof(uint8_t*));
  for (uint32_t i = 0; i < tau; i++) {
    if (i < params->faest_param.t0) {
      b = 0;
    } else {
      b = 1;
    }
    // TODO: the bits are represented as uint8,, hopefully its not too bad
    uint8_t* chalout;
    ChalDec(chal_3, i, params->faest_param.k0, params->faest_param.t0, params->faest_param.k1,
            params->faest_param.t1, chalout);
    d_tilde[i] = malloc(8);
    for (uint32_t j = 0; j < 8; j++) {
      d_tilde[i][j] = chalout[j] ^ signature->u_tilde[j];
    }
  }

  // Step: 17..18
  uint8_t* q_tilde_tmp = malloc(tau * 8);
  for (uint32_t i = 0; i < tau; i++) {
    for (uint8_t j = 0; j < 8; j++) {
      *(q_tilde_tmp + (i * tau) + j) = *(q_tilde[i] + j) ^ *(d_tilde[i] + j);
    }
  }
  uint8_t* hv_verify = malloc(lambdaBytes * 2);
  H1_context_t h1_ctx_1;
  H1_init(&h1_ctx_1, lambda);
  H1_update(&h1_ctx_1, q_tilde_tmp, tau * 8);
  H1_final(&h1_ctx_1, hv_verify, lambdaBytes * 2);
  if (memcmp(signature->h_v, hv_verify, lambdaBytes * 2) != 0) {
    return 0;
  }

  uint32_t ret;
  aesVerify(signature->d, q, chal_2, chal_3, tau, signature->a_tilde, signature->b_tilde, in, out,
            ret);

  if (ret == 0) {
    return 0;
  }
  return 1;
}