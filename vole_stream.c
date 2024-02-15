#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "vole_stream.h"
#include "aes.h"
#include "utils.h"
#include "random_oracle.h"

#include <stdbool.h>
#include <string.h>

static void StreamConstructVole(const uint8_t* iv, stream_vec_com_t* sVecCom, unsigned int lambda,
                                unsigned int outLenBytes, uint8_t* u, uint8_t* v, uint8_t* h,
                                unsigned int begin, unsigned int end) {
  unsigned int depth               = sVecCom->depth;
  const unsigned int num_instances = 1 << depth;
  const unsigned int lambda_bytes  = lambda / 8;
  unsigned int len                 = end - begin;

#define V(idx) (v + ((idx)-begin) * outLenBytes)

  uint8_t* sd  = malloc(lambda_bytes);
  uint8_t* com = malloc(lambda_bytes * 2);
  H1_context_t* h1_ctx;
  if (h != NULL) {
    h1_ctx = malloc(sizeof(H1_context_t));
    H1_init(h1_ctx, lambda);
  }

  uint8_t* r = malloc(outLenBytes);

  // Clear initial memory
  if (u != NULL) {
    memset(u, 0, outLenBytes);
  }
  memset(v, 0, len * outLenBytes);

  for (unsigned int i = 0; i < num_instances; i++) {
    get_sd_com(sVecCom, iv, lambda, i, sd, com);
    if (h != NULL) {
      H1_update(h1_ctx, com, lambda_bytes * 2);
    }

    prg(sd, iv, r, lambda, outLenBytes);
    if (u != NULL) {
      xor_u8_array(u, r, u, outLenBytes);
    }
    for (unsigned int j = begin; j < end; j++) {
      // Apply r if the j'th bit is set
      if ((i >> j) & 1) {
        xor_u8_array(V(j), r, V(j), outLenBytes); // FIXME - arhgh!
      }
    }
  }

  if (h != NULL) {
    H1_final(h1_ctx, h, lambda_bytes * 2);
  }

  if (h != NULL){
    free(h1_ctx);
  }
  free(sd);
  free(com);
  free(r);
}

void partial_vole_commit_cmo(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                             const faest_paramset_t* params, stream_vec_com_t* sVecCom, uint8_t** v,
                             unsigned int start, unsigned int len) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;
  unsigned int max_depth    = MAX(k0, k1);

  uint8_t* expanded_keys = malloc(tau * lambda_bytes);
  prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);
  uint8_t* path = malloc(lambda_bytes * max_depth);

  unsigned int end        = start + len;
  unsigned int tree_start = 0;
  for (unsigned int i = 0; i < tau; i++) {
    unsigned int depth = i < tau0 ? k0 : k1;

    if (tree_start + depth > start) {
      stream_vector_commitment(expanded_keys + i * lambda_bytes, lambda, &sVecCom[i], depth);

      unsigned int lbegin = 0;
      if (start > tree_start) {
        assert(start >= tree_start);
        lbegin = start - tree_start;
      }

      assert(end > tree_start);
      unsigned int lend  = MIN(depth, end - tree_start);
      unsigned int v_idx = 0;
      if (tree_start > start) {
        v_idx = tree_start - start;
      }

      assert(lend <= depth);

      sVecCom[i].path = path;
      StreamConstructVole(iv, &sVecCom[i], lambda, ellhat_bytes, NULL, v[v_idx], NULL, lbegin,
                          lend);
      sVecCom[i].path = NULL;
    }

    tree_start += depth;
    if (tree_start >= end) {
      break;
    }
  }
  free(expanded_keys);
  free(path);
}

void stream_vole_commit(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                        const faest_paramset_t* params, uint8_t* hcom, stream_vec_com_t* sVecCom,
                        uint8_t* c, uint8_t* u, uint8_t** v) {
  unsigned int lambda       = params->faest_param.lambda;
  unsigned int lambda_bytes = lambda / 8;
  unsigned int ellhat_bytes = (ellhat + 7) / 8;
  unsigned int tau          = params->faest_param.tau;
  unsigned int tau0         = params->faest_param.t0;
  unsigned int k0           = params->faest_param.k0;
  unsigned int k1           = params->faest_param.k1;
  unsigned int max_depth    = MAX(k0, k1);

  uint8_t* ui = malloc(tau * ellhat_bytes);

  // Step 1
  uint8_t* expanded_keys = malloc(tau * lambda_bytes);
  prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);

  // for Step 12
  H1_context_t h1_ctx;
  H1_init(&h1_ctx, lambda);
  uint8_t* h    = malloc(lambda_bytes * 2);
  uint8_t* path = malloc(lambda_bytes * max_depth);

  unsigned int v_idx = 0;
  for (unsigned int i = 0; i < tau; i++) {
    // Step 4
    unsigned int depth = i < tau0 ? k0 : k1;
    // Step 5
    stream_vector_commitment(expanded_keys + i * lambda_bytes, lambda, &sVecCom[i], depth);
    // Step 6
    sVecCom[i].path = path;
    StreamConstructVole(iv, &sVecCom[i], lambda, ellhat_bytes, ui + i * ellhat_bytes, v[v_idx], h,
                        0, depth);
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
  free(ui);
  free(expanded_keys);
  free(h);
  free(path);
}
