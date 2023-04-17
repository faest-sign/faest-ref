// Reference - https://nvlpubs.nist.gov/nistpubs/fips/nist.fips.197.pdf
// Tested against the given key, plaintext, ciphertext in Appendix B - Cipher Example

// TODO - Clean up the mess
// TODO - Integrate with the existing framework
// TODO - 192 and 256 support
// TODO - Benchmark and optimize implementation

// AES CTR mode 128/192/256

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "aes.h"

#define NC 4 // total coloumns
#define NR 4 // total rows
#define NK 4 // 8*4 bits key
#define NROUND 10 // total rounds

#define modulus_bf8 (UINT8_C((1 << 4) | (1 << 3) | (1 << 1) | 1)) // 0x1B

bf8_t MULbf8(bf8_t lhs, bf8_t rhs) {
  bf8_t result = 0;
  for (unsigned int idx = 8; idx; --idx, rhs >>= 1) {
    result ^= (-(rhs & 1)) & lhs;
    const uint8_t mask = -((lhs >> 7) & 1);
    lhs                = (lhs << 1) ^ (mask & modulus_bf8);
  }
  return result;
}

static bf8_t GetBit(bf8_t in, uint8_t index) {
    return (in >> index) & 0x01;
}

static bf8_t SetBit(bf8_t in, uint8_t index) {
    return (in << index);
}

// TODO: Use Intel intrincics ??
static void AddRoundKey(bf8_t round, state_t* state) {
    // printf("\n");
    for(uint8_t c = 0; c < NC; c++) {
        for(uint8_t r = 0; r < NR; r++) {
            // printf("%.2x ",(*state)[c][r]);
            (*state)[c][r] ^= round_key[(NC * round * 4) + (c * NC) + r];
        }
        // printf("\n");
    }
}

static const bf8_t sbox[256] = {
    //  0       1   2       3   4       5   6       7   8       9   A       B   C       D   E       F
        0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76, // 0
        0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0, // 1
        0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15, // 2
        0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75, // 3
        0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84, // 4
        0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf, // 5
        0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8, // 6
        0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2, // 7
        0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73, // 8
        0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb, // 9
        0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79, // a
        0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08, // b
        0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a, // c
        0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e, // d
        0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf, // e
        0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 // f
};

// I guess this is the Itoh-Tsujii inverse... found it from a previous implementation of mine
static bf8_t ModuloInverse(bf8_t in) {

    if(0x00 == in) {
        return 0x00;
    }

    uint16_t t1 = in;
    uint16_t t2 = t1;

    for(size_t i = 0; i < 6; i++) {
        t2 = MULbf8(t2,t2);
        t1 = MULbf8(t1, t2);
    }

    t1 = MULbf8(t1,t1);
    return (bf8_t)t1;
}


static void SubBytes(state_t* state) {

    for(uint8_t c = 0; c < NC; c++) {
        for(uint8_t r = 0; r < NR; r++) {

            bf8_t t = ModuloInverse((*state)[c][r]);
            bf8_t t0 = SetBit( GetBit(t,0) ^ GetBit(t,4) ^ GetBit(t,5) ^ GetBit(t,6) ^ GetBit(t,7) ^ 0x01, 0);
            bf8_t t1 = SetBit( GetBit(t,0) ^ GetBit(t,1) ^ GetBit(t,5) ^ GetBit(t,6) ^ GetBit(t,7) ^ 0x01, 1);
            bf8_t t2 = SetBit( GetBit(t,0) ^ GetBit(t,1) ^ GetBit(t,2) ^ GetBit(t,6) ^ GetBit(t,7), 2);
            bf8_t t3 = SetBit( GetBit(t,0) ^ GetBit(t,1) ^ GetBit(t,2) ^ GetBit(t,3) ^ GetBit(t,7), 3);
            bf8_t t4 = SetBit( GetBit(t,0) ^ GetBit(t,1) ^ GetBit(t,2) ^ GetBit(t,3) ^ GetBit(t,4), 4);
            bf8_t t5 = SetBit( GetBit(t,1) ^ GetBit(t,2) ^ GetBit(t,3) ^ GetBit(t,4) ^ GetBit(t,5) ^ 0x01, 5);
            bf8_t t6 = SetBit( GetBit(t,2) ^ GetBit(t,3) ^ GetBit(t,4) ^ GetBit(t,5) ^ GetBit(t,6) ^ 0x01, 6);
            bf8_t t7 = SetBit( GetBit(t,3) ^ GetBit(t,4) ^ GetBit(t,5) ^ GetBit(t,6) ^ GetBit(t,7), 7);
            (*state)[c][r] = t0 ^ t1 ^ t2 ^ t3 ^ t4 ^ t5 ^ t6 ^ t7;

            // (*state)[c][r] = sbox[(*state)[c][r]];
        }
    }
}

static void ShiftRow(state_t* state) {
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


static void MixColoumn(state_t* state) {
    for(uint8_t c = 0; c < NC; c++) {

        bf8_t tmp = MULbf8((*state)[c][0],0x02) ^ MULbf8((*state)[c][1],0x03) 
                        ^ (*state)[c][2] ^ (*state)[c][3];

        bf8_t tmp_1 = (*state)[c][0] ^ MULbf8((*state)[c][1],(bf8_t)0x02) 
                        ^ MULbf8((*state)[c][2],(bf8_t)0x03) ^ (*state)[c][3];

        bf8_t tmp_2 = (*state)[c][0] ^ (*state)[c][1] ^ MULbf8((*state)[c][2],(bf8_t)0x02) 
                        ^ MULbf8((*state)[c][3],(bf8_t)0x03);
                        
        bf8_t tmp_3 = MULbf8((*state)[c][0],(bf8_t)0x03) ^ (*state)[c][1] ^ (*state)[c][2] 
                        ^ MULbf8((*state)[c][3],(bf8_t)0x02);

        (*state)[c][0] = tmp;
        (*state)[c][1] = tmp_1;
        (*state)[c][2] = tmp_2;
        (*state)[c][3] = tmp_3;
    }
}

// Key Expansion functions
static void SubWords(bf8_t* words) {
    words[0] = sbox[words[0]];
    words[1] = sbox[words[1]];
    words[2] = sbox[words[2]];
    words[3] = sbox[words[3]];
}

static void RotWord(bf8_t* words) {
    bf8_t tmp = words[0];
    words[0] = words[1];
    words[1] = words[2];
    words[2] = words[3];
    words[3] = tmp;
}

static bf8_t Rcon(uint8_t in) {
    bf8_t tmp = 1;
    for(uint8_t i = 1; i < in; i++) {
        tmp = MULbf8(tmp,(bf8_t)0x02);
    }
    return tmp;
}

static void KeyExpansion(bf8_t* key) {
    
    for(uint8_t k = 0; k < NK; k++) {
        round_key[4*k] = key[4*k];
        round_key[(4*k)+1] = key[(4*k)+1];
        round_key[(4*k)+2] = key[(4*k)+2];
        round_key[(4*k)+3] = key[(4*k)+3];
    }

    // printf("Round 0 Key \n");
    // for(uint8_t i = 0; i < 16; i += 4) {
    //     printf("%.2x %.2x %.2x %.2x \n",round_key[i],round_key[i+1],round_key[i+2],round_key[i+3]);
    // }

    bf8_t tmp[4];
    for(uint8_t k = NK; k < NC * (NROUND+1); k++) {
        uint8_t j = (k-1)*4;
        tmp[0] = round_key[j];
        tmp[1] = round_key[j+1];
        tmp[2] = round_key[j+2];
        tmp[3] = round_key[j+3];

        if (k % NK == 0) {
            RotWord(tmp);
            SubWords(tmp);
            tmp[0] ^= Rcon(k/NK);
            // printf("Round %d Key \n", k/nc);
        }

        j = k*4;
        uint8_t m = (k - NK) * 4;
        round_key[j] = round_key[m] ^ tmp[0];
        round_key[j+1] = round_key[m + 1] ^ tmp[1];
        round_key[j+2] = round_key[m + 2] ^ tmp[2];
        round_key[j+3] = round_key[m + 3] ^ tmp[3];

        // printf("%.2x %.2x %.2x %.2x \n",round_key[j],round_key[j+1],round_key[j+2],round_key[j+3]);
        
    }

}

// Cipher
static void Cipher(state_t* state) {
    uint8_t round = 0;

    // first round
    AddRoundKey(round, state);

    for(round = 1; round < NROUND; round++) {
        SubBytes(state);
        ShiftRow(state);
        MixColoumn(state);
        AddRoundKey(round,state);
    }

    // // last round
    SubBytes(state);
    ShiftRow(state);
    AddRoundKey(round,state);

}

// Getting all round keys
void Initialize(bf8_t* key, bf8_t* iv_) {
    KeyExpansion(key);
    memcpy(iv, iv_, 16);
}

static void IncrementIV() {
    for(uint8_t i = 0; i < 16; i++) {
        if(iv[i] == 255) {
            iv[i] = (bf8_t)0x00;
            continue;
        }
        iv[i] ^= (bf8_t)0x01;
        break;
    }
}

// Performs en/decrypt in CTR mode
void Encrypt(bf8_t* buffer) {

    memcpy(buffer,iv,16);

    Cipher((state_t*)buffer);

    IncrementIV();
}

/*int main(int argc, char* argv[]) {
    bf8_t key[16] = {0x2b, 0x7e, 0x15, 0x16,
                    0x28, 0xae, 0xd2, 0xa6,
                    0xab, 0xf7, 0x15, 0x88,
                    0x09, 0xcf, 0x4f, 0x3c};
    bf8_t iv[16] = {0x32, 0x43, 0xf6, 0xa8,
                    0x88, 0x5a, 0x30, 0x8d,
                    0x31, 0x31, 0x98, 0xa2,
                    0xe0, 0x37, 0x07, 0x34};
    bf8_t expect_out[16] = {0x39,0x25,0x84,0x1d,
                            0x02,0xdc,0x09,0xfb,
                            0xdc,0x11,0x85,0x97,
                            0x19,0x6a,0x0b,0x32};

    Initialize(key, iv);

    bf8_t buffer[16];
    Encrypt(buffer);

    printf("Round Output \n");
    for(uint8_t i = 0; i < 16; i += 4) {
        printf("%.2x %.2x %.2x %.2x \n",buffer[i],buffer[i+1],buffer[i+2],buffer[i+3]);
    }

    uint8_t match = 1;
    for(uint8_t i = 0; i < 16; i += 1) {
        if(expect_out[i] != buffer[i]) {
            match = 0;
        }
    }
    if (match == 1) {
        printf("Expected Matches");
    } else {
        printf("Expected Didn't Match");
    }
}*/