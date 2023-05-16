#include "vole.h"

int ChalDec(const uint8_t* chal, uint32_t i, uint32_t k0, uint32_t t0, uint32_t k1, uint32_t t1,
            uint8_t* chalout);

void voleCommit(uint8_t* rootKey, uint32_t lambda, uint32_t outlen, uint32_t tau, uint32_t k0,
                uint32_t k1, const faest_paramset_t* params, uint8_t* hcom, vec_com_t** vecCom,
                uint8_t** c, uint8_t* u, uint8_t** v);

void voleVerify(const faest_paramset_t* params, const uint8_t* chal, const uint8_t** pdec,
                const uint8_t** com_j, uint32_t lambda, uint32_t outlen, uint32_t tau, uint32_t k0,
                uint32_t k1, uint8_t* hcom, uint8_t** u, uint8_t** q, vec_com_t** vecCom,
                vec_com_rec_t** vecComRec);