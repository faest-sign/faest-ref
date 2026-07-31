#include "tables.h"
#include "tables_128s.h"
#include "tables_128f.h"
#include "tables_192s.h"
#include "tables_192f.h"
#include "tables_em_192s.h"
#include "tables_256s.h"
#include "tables_256f.h"

// NOTE: Claude generated ---
#define TABLES(name)                                                 \
  {                                                                  \
    &name##_F[0][0], &name##_G[0][0],                                \
    &name##_W_TREE[0][0], &name##_W_GATE[0][0], &name##_W_CRT[0][0], \
    &name##_TREE_MODULI[0], &name##_M_TREE[0],                       \
    name##_NGATES, name##_NDELTA_BITS, name##_NTREE_BITS, name##_TAU,\
    name##_F_WORDS, name##_G_WORDS, name##_W_TREE_WORDS,             \
    name##_W_GATE_WORDS, name##_W_CRT_WORDS, name##_M_TREE_WORDS,    \
  }

#define CASE_TABLES(P)          CASE_TABLES_AS(P, P)
#define CASE_TABLES_AS(P, NAME) \
  case P: { static const faest_tables_t t = TABLES(NAME); return &t; }

const faest_tables_t* faest_get_tables(faest_paramid_t paramid) {
  switch (paramid) {
    CASE_TABLES(FAEST_128F)
    CASE_TABLES(FAEST_128S)
    CASE_TABLES(FAEST_192F)
    CASE_TABLES(FAEST_192S)
    CASE_TABLES(FAEST_256F)
    CASE_TABLES(FAEST_256S)

    CASE_TABLES_AS(FAEST_EM_128F, FAEST_128F)
    CASE_TABLES_AS(FAEST_EM_128S, FAEST_128S)
    CASE_TABLES_AS(FAEST_EM_192S, FAEST_EM_192S)
    CASE_TABLES_AS(FAEST_EM_192F, FAEST_192F)
    CASE_TABLES_AS(FAEST_EM_256F, FAEST_256F)
    CASE_TABLES_AS(FAEST_EM_256S, FAEST_256S)    

    default: return NULL;
  }
}
// ----