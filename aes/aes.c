// AES CTR mode 128/192/256
// Reference - https://nvlpubs.nist.gov/nistpubs/fips/nist.fips.197.pdf
// Tested against the given key, plaintext, ciphertext in Appendix B - Cipher Example

// TODO - 192 and 256 support

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "aes.h"

#define NC 4 // total coloumns
#define NR 4 // total rows
#define NK 4 // 8*4 bits key
#define NROUND 10 // total rounds

#define MODULUS_BF8 (UINT8_C((1 << 4) | (1 << 3) | (1 << 1) | 1)) // 0x1B

bf8_t mul_bf8(bf8_t lhs, bf8_t rhs) {
  bf8_t result = 0;
  for (unsigned int idx = 8; idx; --idx, rhs >>= 1) {
    result ^= (-(rhs & 1)) & lhs;
    const uint8_t mask = -((lhs >> 7) & 1);
    lhs                = (lhs << 1) ^ (mask & MODULUS_BF8);
  }
  return result;
}

static bf8_t get_bit(bf8_t in, uint8_t index) {
    return (in >> index) & 0x01;
}

static bf8_t set_bit(bf8_t in, uint8_t index) {
    return (in << index);
}

static bf8_t modulo_inverse(bf8_t in) {
    if(0x00 == in) {
        return 0x00;
    }
    uint16_t t1 = in;
    uint16_t t2 = t1;
    for(size_t i = 0; i < 6; i++) {
        t2 = mul_bf8(t2,t2);
        t1 = mul_bf8(t1, t2);
    }
    t1 = mul_bf8(t1,t1);
    return (bf8_t)t1;
}

static bf8_t compute_sbox(bf8_t in) {
    bf8_t t = modulo_inverse(in);
    bf8_t t0 = set_bit( get_bit(t,0) ^ get_bit(t,4) ^ get_bit(t,5) ^ get_bit(t,6) ^ get_bit(t,7) ^ 0x01, 0);
    bf8_t t1 = set_bit( get_bit(t,0) ^ get_bit(t,1) ^ get_bit(t,5) ^ get_bit(t,6) ^ get_bit(t,7) ^ 0x01, 1);
    bf8_t t2 = set_bit( get_bit(t,0) ^ get_bit(t,1) ^ get_bit(t,2) ^ get_bit(t,6) ^ get_bit(t,7), 2);
    bf8_t t3 = set_bit( get_bit(t,0) ^ get_bit(t,1) ^ get_bit(t,2) ^ get_bit(t,3) ^ get_bit(t,7), 3);
    bf8_t t4 = set_bit( get_bit(t,0) ^ get_bit(t,1) ^ get_bit(t,2) ^ get_bit(t,3) ^ get_bit(t,4), 4);
    bf8_t t5 = set_bit( get_bit(t,1) ^ get_bit(t,2) ^ get_bit(t,3) ^ get_bit(t,4) ^ get_bit(t,5) ^ 0x01, 5);
    bf8_t t6 = set_bit( get_bit(t,2) ^ get_bit(t,3) ^ get_bit(t,4) ^ get_bit(t,5) ^ get_bit(t,6) ^ 0x01, 6);
    bf8_t t7 = set_bit( get_bit(t,3) ^ get_bit(t,4) ^ get_bit(t,5) ^ get_bit(t,6) ^ get_bit(t,7), 7);
    return t0 ^ t1 ^ t2 ^ t3 ^ t4 ^ t5 ^ t6 ^ t7;
}

// Round Functions

static void add_round_key(bf8_t round, state_t* state) {
    for(uint8_t c = 0; c < NC; c++) {
        for(uint8_t r = 0; r < NR; r++) {
            (*state)[c][r] ^= round_key[(NC * round * 4) + (c * NC) + r];
        }
    }
}

static void sub_bytes(state_t* state) {
    for(uint8_t c = 0; c < NC; c++) {
        for(uint8_t r = 0; r < NR; r++) {
            (*state)[c][r] = compute_sbox((*state)[c][r]);
        }
    }
}

static void shift_row(state_t* state) {
    bf8_t tmp = (*state)[0][1];
    (*state)[0][1] = (*state)[1][1];
    (*state)[1][1] = (*state)[2][1];
    (*state)[2][1] = (*state)[3][1];
    (*state)[3][1] = tmp;

    bf8_t tmp_20 = (*state)[0][2];
    bf8_t tmp_21 = (*state)[1][2];
    (*state)[0][2] = (*state)[2][2];
    (*state)[1][2] = (*state)[3][2];
    (*state)[2][2] = tmp_20;
    (*state)[3][2] = tmp_21;

    tmp = (*state)[0][3];
    (*state)[0][3] = (*state)[3][3];
    (*state)[3][3] = (*state)[2][3];
    (*state)[2][3] = (*state)[1][3];
    (*state)[1][3] = tmp;
}

static void mix_coloumn(state_t* state) {
    for(uint8_t c = 0; c < NC; c++) {
        bf8_t tmp = mul_bf8((*state)[c][0],0x02) ^ mul_bf8((*state)[c][1],0x03) 
                        ^ (*state)[c][2] ^ (*state)[c][3];
        bf8_t tmp_1 = (*state)[c][0] ^ mul_bf8((*state)[c][1],0x02) 
                        ^ mul_bf8((*state)[c][2],0x03) ^ (*state)[c][3];
        bf8_t tmp_2 = (*state)[c][0] ^ (*state)[c][1] ^ mul_bf8((*state)[c][2],0x02) 
                        ^ mul_bf8((*state)[c][3],0x03);
        bf8_t tmp_3 = mul_bf8((*state)[c][0],0x03) ^ (*state)[c][1] ^ (*state)[c][2] 
                        ^ mul_bf8((*state)[c][3],0x02);

        (*state)[c][0] = tmp;
        (*state)[c][1] = tmp_1;
        (*state)[c][2] = tmp_2;
        (*state)[c][3] = tmp_3;
    }
}

// Key Expansion functions

static void sub_words(bf8_t* words) {
    words[0] = compute_sbox(words[0]);
    words[1] = compute_sbox(words[1]);
    words[2] = compute_sbox(words[2]);
    words[3] = compute_sbox(words[3]);
}

static void rot_word(bf8_t* words) {
    bf8_t tmp = words[0];
    words[0] = words[1];
    words[1] = words[2];
    words[2] = words[3];
    words[3] = tmp;
}

static bf8_t round_const(uint8_t in) {
    bf8_t tmp = 1;
    for(uint8_t i = 1; i < in; i++) {
        tmp = mul_bf8(tmp,0x02);
    }
    return tmp;
}

// Main Functions

static void key_expansion(bf8_t* key) {
    for(uint8_t k = 0; k < NK; k++) {
        round_key[4*k] = key[4*k];
        round_key[(4*k)+1] = key[(4*k)+1];
        round_key[(4*k)+2] = key[(4*k)+2];
        round_key[(4*k)+3] = key[(4*k)+3];
    }
    bf8_t tmp[4];
    for(uint8_t k = NK; k < NC * (NROUND+1); k++) {
        uint8_t j = (k-1)*4;
        tmp[0] = round_key[j];
        tmp[1] = round_key[j+1];
        tmp[2] = round_key[j+2];
        tmp[3] = round_key[j+3];

        if (k % NK == 0) {
            rot_word(tmp);
            sub_words(tmp);
            tmp[0] ^= round_const(k/NK);
        }

        j = k*4;
        uint8_t m = (k - NK) * 4;
        round_key[j] = round_key[m] ^ tmp[0];
        round_key[j+1] = round_key[m + 1] ^ tmp[1];
        round_key[j+2] = round_key[m + 2] ^ tmp[2];
        round_key[j+3] = round_key[m + 3] ^ tmp[3];
    }
}

static void cipher(state_t* state) {
    uint8_t round = 0;

    // first round
    add_round_key(round, state);

    for(round = 1; round < NROUND; round++) {
        sub_bytes(state);
        shift_row(state);
        mix_coloumn(state);
        add_round_key(round,state);
    }

    // last round
    sub_bytes(state);
    shift_row(state);
    add_round_key(round,state);
}

static void increment_iv() {
    for(uint8_t i = 0; i < 16; i++) {
        if(iv[i] == 255) {
            iv[i] = 0x00;
            continue;
        }
        iv[i] ^= 0x01;
        break;
    }
}

// Calling Functions

void initialize(bf8_t* key, bf8_t* iv_) {
    key_expansion(key);
    memcpy(iv, iv_, 16);
}

void encrypt(bf8_t* buffer) {
    memcpy(buffer,iv,16);
    cipher((state_t*)buffer);
    increment_iv();
}

