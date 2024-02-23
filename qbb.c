#include "qbb.h"

void init_qbb(qbb_t* qbb, unsigned int len, const uint8_t* iv, uint8_t* c,
              uint8_t* pdec_sig, uint8_t* com_sig, uint8_t* chall3, uint8_t* u_tilde, const faest_paramset_t* params) {
  qbb->iv        = iv;
  qbb->row_count = len;
  qbb->params    = params;
  qbb->iv        = iv;
  qbb->c         = c;
  qbb->com_hash  = calloc(MAX_LAMBDA_BYTES * 2, sizeof(uint8_t));

  const unsigned int lambda        = params->faest_param.lambda;
  const unsigned int l             = params->faest_param.l;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int tau           = params->faest_param.tau;
  const unsigned int tau0          = params->faest_param.t0;
  const size_t lambda_bytes        = params->faest_param.lambda / 8;
  const unsigned int utilde_bytes  = lambda_bytes + UNIVERSAL_HASH_B;
  const unsigned int k0            = params->faest_param.k0;
  const unsigned int k1            = params->faest_param.k1;

  uint8_t** qprime = malloc(lambda * sizeof(uint8_t*));
  qprime[0]        = calloc(lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    qprime[i] = qprime[0] + i * ell_hat_bytes;
  }
  const uint8_t* pdec[MAX_TAU];
  const uint8_t* com[MAX_TAU];
  for (unsigned int i = 0; i < tau; ++i) {
    // TODO compute offset
    if (i < tau0) {
      pdec[i] = pdec_sig + i * (params->faest_param.k0 + 2) * lambda_bytes;
      com[i]  = com_sig + (i * (params->faest_param.k0 + 2) + params->faest_param.k0) * lambda_bytes;
    } else {
      pdec[i] =
          pdec_sig + ((i - tau0) * (params->faest_param.k1 + 2) + tau0 * (params->faest_param.k0 + 2)) *
                     lambda_bytes;
      com[i] = com_sig + ((i - tau0) * (params->faest_param.k1 + 2) + params->faest_param.k1 +
                      tau0 * (params->faest_param.k0 + 2)) *
                         lambda_bytes;
    }
  }
  vole_reconstruct(iv, chall3, pdec, com, qbb->com_hash, qprime, ell_hat, params);

  uint8_t** q = malloc(lambda * sizeof(uint8_t*));
  q[0]        = calloc(lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    q[i] = q[0] + i * ell_hat_bytes;
  }

  uint8_t** Dtilde = malloc(lambda * sizeof(uint8_t*));
  Dtilde[0]        = calloc(lambda, (lambda_bytes + UNIVERSAL_HASH_B));
  for (unsigned int i = 1; i < lambda; ++i) {
    Dtilde[i] = Dtilde[0] + i * (lambda_bytes + UNIVERSAL_HASH_B);
  }

  unsigned int Dtilde_idx = 0;
  unsigned int q_idx      = 0;
  for (unsigned int i = 0; i < tau; i++) {
    const unsigned int depth = i < tau0 ? k0 : k1;

    // Step 11
    uint8_t delta[MAX_DEPTH];
    ChalDec(chall3, i, params->faest_param.k0, params->faest_param.t0,
            params->faest_param.k1, params->faest_param.t1, delta);
    // Step 16
    for (unsigned int j = 0; j != depth; ++j, ++Dtilde_idx) {
      // for scan-build
      assert(Dtilde_idx < lambda);
      masked_xor_u8_array(Dtilde[Dtilde_idx], u_tilde, Dtilde[Dtilde_idx],
                          delta[j], utilde_bytes);
    }

    if (i == 0) {
      // Step 8
      memcpy(q[q_idx], qprime[q_idx], ell_hat_bytes * depth);
      q_idx += depth;
    } else {
      // Step 14
      for (unsigned int d = 0; d < depth; ++d, ++q_idx) {
        masked_xor_u8_array(qprime[q_idx], c + (i - 1)*ell_hat_bytes, q[q_idx], delta[d],
                            ell_hat_bytes);
      }
    }
  }
  free(qprime[0]);
  free(qprime);
  qprime = NULL;
  qbb->vole_Q_cache = q[0];
}

uint8_t* get_vole_q_hash(qbb_t* qbb, unsigned int idx) {
  const unsigned int ellhat =
      qbb->params->faest_param.l + qbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int offset       = idx - qbb->cache_idx;
  return qbb->vole_Q_cache + offset * ellhat_bytes;
}

uint8_t* get_vole_q_prove(qbb_t* qbb, unsigned int idx) {
  unsigned int offset = (idx - qbb->cache_idx) * (qbb->params->faest_param.lambda / 8);
  return qbb->vole_Q_cache + offset;
}

bf256_t* get_vole_q_prove_256(qbb_t* qbb, unsigned int idx) {
  return (bf256_t*)get_vole_q_prove(qbb, idx);
}