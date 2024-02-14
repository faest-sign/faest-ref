#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "vole_stream.h"
#include "aes.h"
#include "utils.h"
#include "random_oracle.h"

#include <stdbool.h>
#include <string.h>

/*
int ChalDec(const uint8_t* chal, unsigned int i, unsigned int k0, unsigned int t0, unsigned int k1,
            unsigned int t1, uint8_t* chalout) {
  if (i >= t0 + t1) {
    return 0;
  }

  unsigned int lo;
  unsigned int hi;
  if (i < t0) {
    lo = i * k0;
    hi = ((i + 1) * k0);
  } else {
    unsigned int t = i - t0;
    lo             = (t0 * k0) + (t * k1);
    hi             = (t0 * k0) + ((t + 1) * k1);
  }

  assert(hi - lo == k0 || hi - lo == k1);
  for (unsigned int j = lo; j < hi; ++j) {
    // set_bit(chalout, i - lo, get_bit(chal, i));
    chalout[j - lo] = ptr_get_bit(chal, j);
  }
  return 1;
}
*/

void StreamConstructVole(const uint8_t* iv, stream_vec_com_t* sVecCom, unsigned int lambda, unsigned int outLenBytes, uint8_t* u, uint8_t* v, uint8_t* h) {
  unsigned int depth = sVecCom->depth;
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;

  uint8_t* stack = calloc(depth+1, outLenBytes);
  memset(v, 0, depth * outLenBytes);
  unsigned int stack_index = 1;
#define V(idx) (v + (idx)*outLenBytes)
#define STACK_PEAK(depth) (stack + (stack_index - 1 - depth) * outLenBytes)


  uint8_t* sd = alloca(lambda_bytes);
  uint8_t* com = alloca(lambda_bytes * 2);

  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);

  get_sd_com(sVecCom, iv, lambda, 0, sd, com);
  H1_update(&h1_ctx, com, lambda_bytes * 2);
  prg(sd, iv, stack, lambda, outLenBytes);

  unsigned int j;
  // Step: 3..4
  for (unsigned int i = 1; i < num_instances; i++) {
    get_sd_com(sVecCom, iv, lambda, i, sd, com);
    H1_update(&h1_ctx, com, lambda_bytes * 2);

    prg(sd, iv, STACK_PEAK(-1), lambda, outLenBytes);
    stack_index++;

    j = 0;
    for (unsigned int x = i + 1; !(x & 1); x >>= 1) {
      xor_u8_array(V(j), STACK_PEAK(0), V(j), outLenBytes);
      ++j;
      xor_u8_array(STACK_PEAK(0), STACK_PEAK(1), STACK_PEAK(1), outLenBytes);
      stack_index--;
    }
  }
  
  // Step: 10
  if (u != NULL) {
    memcpy(u, stack, outLenBytes);
  }
  free(stack);

  H1_final(&h1_ctx, h, lambda_bytes * 2);
}

void stream_vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                        const faest_paramset_t* params, uint8_t* hcom, stream_vec_com_t* sVecCom, uint8_t* c,
                        uint8_t* u, uint8_t** v) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;
  unsigned int max_depth = MAX(k0, k1);

  uint8_t* ui = alloca(tau * ellhat_bytes);

  // Step 1
  uint8_t* expanded_keys = alloca(tau * lambda_bytes);
  prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);

  // for Step 12
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  uint8_t* h = alloca(lambda_bytes * 2);
  uint8_t* path = alloca(lambda_bytes * max_depth);

  unsigned int v_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    // Step 4
    unsigned int depth = i < tau0 ? k0 : k1;
    // Step 5
    stream_vector_commitment(expanded_keys + i * lambda_bytes, lambda, &sVecCom[i], depth);
    // Step 6
    sVecCom[i].path = path;
    StreamConstructVole(iv, &sVecCom[i], lambda, ellhat_bytes, ui + i * ellhat_bytes, v[v_idx], h);
    sVecCom[i].path = NULL;
    // Step 7 (and parts of 8)
    v_idx += depth;
    // Step 12 (part)
    H1_update(&h1_ctx, h, lambda_bytes * 2);
  }
  // Step 9
  memcpy(u, ui, ellhat_bytes);
  for (unsigned int i = 1; i < tau; i++) {
    // Step 11
    xor_u8_array(u, ui + i * ellhat_bytes, c + (i - 1) * ellhat_bytes, ellhat_bytes);
  }

  // Step 12: Generating final commitment from all the com commitments
  H1_final(&h1_ctx, hcom, lambda_bytes * 2);
}
