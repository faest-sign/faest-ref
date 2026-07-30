#ifndef FAEST_TABLES_H
#define FAEST_TABLES_H

#include <stdint.h>
#include "instances.h"   // for faest_paramid_t

#if defined(__cplusplus)
extern "C" {
#endif

typedef struct {
  const uint64_t* F;
  const uint64_t* G;
  const uint64_t* W_TREE;
  const uint64_t* W_GATE;
  const uint64_t* W_CRT;
  const uint64_t* TREE_MODULI;
  const uint64_t* M_TREE;

  unsigned int ngates;
  unsigned int ndelta_bits;
  unsigned int ntree_bits;
  unsigned int tau;
  unsigned int f_words, g_words, w_tree_words, w_gate_words, w_crt_words, m_tree_words;
} faest_tables_t;

const faest_tables_t* faest_get_tables(faest_paramid_t paramid);

#if defined(__cplusplus)
}
#endif

#endif