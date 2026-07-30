#include "tables.h"
#include "tables_128s.h"
#include "tables_128f.h"

#define TABLES(name)                                                 \
  {                                                                  \
    &name##_F[0][0], &name##_G[0][0],                                \
    &name##_W_TREE[0][0], &name##_W_GATE[0][0], &name##_W_CRT[0][0], \
    &name##_TREE_MODULI[0], &name##_M_TREE[0],                       \
    name##_NGATES, name##_NDELTA_BITS, name##_NTREE_BITS, name##_TAU,\
    name##_F_WORDS, name##_G_WORDS, name##_W_TREE_WORDS,             \
    name##_W_GATE_WORDS, name##_W_CRT_WORDS, name##_M_TREE_WORDS,    \
  }

#define CASE_TABLES(P) \
  case P: { static const faest_tables_t t = TABLES(P); return &t; }

const faest_tables_t* faest_get_tables(faest_paramid_t paramid) {
  switch (paramid) {
    CASE_TABLES(FAEST_128S)
    CASE_TABLES(FAEST_128F)
    default: return NULL;
  }
}