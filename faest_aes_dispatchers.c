#include "faest_aes.h"

void aes_128_prover(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde, const uint8_t* w,
                    const uint8_t* u, uint8_t** V, const uint8_t* owf_in, const uint8_t* owf_out,
                    const uint8_t* chall_2, const faest_paramset_t* params);

void aes_192_prover(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde, const uint8_t* w,
                    const uint8_t* u, uint8_t** V, const uint8_t* owf_in, const uint8_t* owf_out,
                    const uint8_t* chall_2, const faest_paramset_t* params);

void aes_256_prover(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde, const uint8_t* w,
                    const uint8_t* u, uint8_t** V, const uint8_t* owf_in, const uint8_t* owf_out,
                    const uint8_t* chall_2, const faest_paramset_t* params);

void aes_128_verifier(uint8_t* a0_tilde, const uint8_t* d, uint8_t** Q, const uint8_t* owf_in,
                      const uint8_t* owf_out, const uint8_t* chall_2, const uint8_t* chall_3,
                      const uint8_t* a1_tilde, const uint8_t* a2_tilde,
                      const faest_paramset_t* params);

void aes_192_verifier(uint8_t* a0_tilde, const uint8_t* d, uint8_t** Q, const uint8_t* owf_in,
                      const uint8_t* owf_out, const uint8_t* chall_2, const uint8_t* chall_3,
                      const uint8_t* a1_tilde, const uint8_t* a2_tilde,
                      const faest_paramset_t* params);

void aes_256_verifier(uint8_t* a0_tilde, const uint8_t* d, uint8_t** Q, const uint8_t* owf_in,
                      const uint8_t* owf_out, const uint8_t* chall_2, const uint8_t* chall_3,
                      const uint8_t* a1_tilde, const uint8_t* a2_tilde,
                      const faest_paramset_t* params);

// AES(-EM) OWF dispatchers
void aes_prove(uint8_t* a0_tilde, uint8_t* a1_tilde, uint8_t* a2_tilde, const uint8_t* w,
               const uint8_t* u, uint8_t** V, const uint8_t* owf_in, const uint8_t* owf_out,
               const uint8_t* chall_2, const faest_paramset_t* params) {
  switch (params->lambda) {
  case 256:
    aes_256_prover(a0_tilde, a1_tilde, a2_tilde, w, u, V, owf_in, owf_out, chall_2, params);
    break;
  case 192:
    aes_192_prover(a0_tilde, a1_tilde, a2_tilde, w, u, V, owf_in, owf_out, chall_2, params);
    break;
  default:
    aes_128_prover(a0_tilde, a1_tilde, a2_tilde, w, u, V, owf_in, owf_out, chall_2, params);
  }
}

void aes_verify(uint8_t* a0_tilde, const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                const uint8_t* chall_3, const uint8_t* a1_tilde, const uint8_t* a2_tilde,
                const uint8_t* owf_in, const uint8_t* owf_out, const faest_paramset_t* params) {
  switch (params->lambda) {
  case 256:
    aes_256_verifier(a0_tilde, d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
    break;
  case 192:
    aes_192_verifier(a0_tilde, d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
    break;
  default:
    aes_128_verifier(a0_tilde, d, Q, owf_in, owf_out, chall_2, chall_3, a1_tilde, a2_tilde, params);
  }
}
