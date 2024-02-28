#include "qbb.h"
#include "parameters.h"
#include "faest_aes.h"
#include "faest.h"

void init_qbb(qbb_t* qbb, unsigned int len, const uint8_t* iv, uint8_t* c, uint8_t* pdec_sig,
              uint8_t* com_sig, uint8_t* chall3, uint8_t* u_tilde, const faest_paramset_t* params,
              const uint8_t* sig) {
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
    pdec[i] = dsignature_pdec(sig, i, params);
    com[i]  = dsignature_com(sig, i, params);
  }

  vole_reconstruct(iv, chall3, pdec, com, qbb->com_hash, qprime, ell_hat, params);

  uint8_t** q = malloc(lambda * sizeof(uint8_t*));
  q[0]        = calloc(lambda, ell_hat_bytes);
  for (unsigned int i = 1; i < lambda; ++i) {
    q[i] = q[0] + i * ell_hat_bytes;
  }

  qbb->Dtilde    = malloc(lambda * sizeof(uint8_t*));
  qbb->Dtilde[0] = calloc(lambda, (lambda_bytes + UNIVERSAL_HASH_B));
  for (unsigned int i = 1; i < lambda; ++i) {
    qbb->Dtilde[i] = qbb->Dtilde[0] + i * (lambda_bytes + UNIVERSAL_HASH_B);
  }

  unsigned int Dtilde_idx = 0;
  unsigned int q_idx      = 0;
  for (unsigned int i = 0; i < tau; i++) {
    const unsigned int depth = i < tau0 ? k0 : k1;

    // Step 11
    uint8_t delta[MAX_DEPTH];
    ChalDec(chall3, i, params->faest_param.k0, params->faest_param.t0, params->faest_param.k1,
            params->faest_param.t1, delta);
    // Step 16
    for (unsigned int j = 0; j != depth; ++j, ++Dtilde_idx) {
      // for scan-build
      assert(Dtilde_idx < lambda);
      masked_xor_u8_array(qbb->Dtilde[Dtilde_idx], u_tilde, qbb->Dtilde[Dtilde_idx], delta[j],
                          utilde_bytes);
    }

    if (i == 0) {
      // Step 8
      memcpy(q[q_idx], qprime[q_idx], ell_hat_bytes * depth);
      q_idx += depth;
    } else {
      // Step 14
      for (unsigned int d = 0; d < depth; ++d, ++q_idx) {
        masked_xor_u8_array(qprime[q_idx], c + (i - 1) * ell_hat_bytes, q[q_idx], delta[d],
                            ell_hat_bytes);
      }
    }
  }
  free(qprime[0]);
  free(qprime);
  qprime                  = NULL;
  qbb->vole_Q_cache       = q[0];
  qbb->vole_Q_cache_index = q;
}

uint8_t* get_vole_q_hash(qbb_t* qbb, unsigned int idx) {
  qbb->cache_idx = 0;
  const unsigned int ellhat =
      qbb->params->faest_param.l + qbb->params->faest_param.lambda * 2 + UNIVERSAL_HASH_B_BITS;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int offset       = idx - qbb->cache_idx;
  return qbb->vole_Q_cache + offset * ellhat_bytes;
}

void prepare_verify_qbb(qbb_t* qbb, const uint8_t* sig_d, const uint8_t* sig_chall_3) {
  qbb->cache_idx                   = 0;
  const unsigned int lambda        = qbb->params->faest_param.lambda;
  const unsigned int l             = qbb->params->faest_param.l;
  const unsigned int ell_hat       = l + lambda * 2 + UNIVERSAL_HASH_B_BITS;
  const unsigned int ell_hat_bytes = ell_hat / 8;
  const unsigned int tau = qbb->params->faest_param.tau;
  const unsigned int t0  = qbb->params->faest_param.t0;
  const unsigned int k0  = qbb->params->faest_param.k0;
  const unsigned int t1  = qbb->params->faest_param.t1;
  const unsigned int k1  = qbb->params->faest_param.k1;

  // TODO: Actually EM use Lenc, but Lenc == L for all EM..
  unsigned int size = qbb->params->faest_param.l;

  for (unsigned int i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t decoded_challenge[MAX_DEPTH];
    ChalDec(sig_chall_3, i, k0, t0, k1, t1, decoded_challenge);
    for (unsigned int j = 0; j < depth; j++, ++col) {
      if (decoded_challenge[j] == 1) {
        xor_u8_array(sig_d, qbb->vole_Q_cache_index[col], qbb->vole_Q_cache_index[col], (size + 7) / 8);
      }
    }
  }

  switch (lambda) {
  case 256:
    qbb->vole_Q_cache_RMO =
        (uint8_t*)column_to_row_major_and_shrink_V_256(qbb->vole_Q_cache_index, FAEST_256F_L);
    break;
  case 192:
    qbb->vole_Q_cache_RMO =
        (uint8_t*)column_to_row_major_and_shrink_V_192(qbb->vole_Q_cache_index, FAEST_192F_L);
    break;
  default:
    qbb->vole_Q_cache_RMO =
        (uint8_t*)column_to_row_major_and_shrink_V_128(qbb->vole_Q_cache_index, FAEST_128F_L);
    break;
  }
}

uint8_t* get_vole_q_verify(qbb_t* qbb, unsigned int idx) {
  unsigned int offset = (idx - qbb->cache_idx) * (qbb->params->faest_param.lambda / 8);
  return qbb->vole_Q_cache_RMO + offset;
}

bf256_t* get_vole_q_verify_256(qbb_t* qbb, unsigned int idx) {
  return (bf256_t*)get_vole_q_verify(qbb, idx);
}
