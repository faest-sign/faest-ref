#include "faest_aes.h"
// TODO: Do not pass lambdaBytes everywhere, compute it in the function....
// TODO: change q to Q where applicable

// TODO: Unsure what is happeninig here...
void aes_extend_witness(uint32_t lambda, const uint8_t* sk, const uint8_t* in, uint8_t* w) {}

int aes_key_schedule_forward(uint32_t lambda, uint32_t m, const uint8_t* x, uint8_t Mtag,
                             uint8_t Mkey, const uint8_t* delta, uint8_t* y) {

  uint32_t lambdaByte = lambda / 8;
  uint32_t numround   = 0;
  // TODO: verify if they mean the word size N is the key size and not the standard block size
  uint32_t N = lambdaByte / 4;
  switch (lambda) {
  case 256:
    numround = 14;
    break;
  case 192:
    numround = 12;
  default:
    numround = 10;
    break;
  }

  // Step: 1
  if ((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)) {
    return 0;
  }

  // STep: 2..3
  // TODO, change to field type
  uint8_t* y = malloc((numround + 1) * 16 * 16);
  memcpy(y, x, 16);

  // Step: 4
  uint32_t i_wd = lambdaByte;

  // Step: 5..10
  // TODO, check if this is not N * (numround + 1) ?
  for (uint32_t j = N; j < 4 * (numround + 1); j++) {
    if (j % N == 0 || (N > 6 && j % N == 4)) {
      // TODO: doing as per the byte thingy, hope its no problem
      memcpy(y + (4 * N), x + i_wd, 4);
      i_wd += 4;
    } else {
      for (uint32_t i = 0; i < 4; i++) {
        y[(4 * j) - i] = y[4 * (j - N) + i] + y[4 * (j - 1) + i];
      }
    }
  }
  return 1;
}

// w (extended witness bits), u (vole masks), v (vole tags), in (faest public-key input block bits),
// out (faest public-key output block bits), chal (quicksilver challenge), a_tilde and b_tilde
// (quicksilver response)
void aes_prove(uint8_t* w, uint8_t* u, uint8_t** v, uint8_t* in, uint8_t* out, uint8_t* chal,
               uint32_t lambda, uint32_t tau, uint32_t l, uint8_t* a_tilde, uint8_t* b_tilde) {

  uint32_t lByte      = l / 8;
  uint32_t lambdaByte = lambda / 8;

  // Step: 1
  // TODO: Is to field, equivalent to the function bf_load ??
  // TODO: Are we really passing it bit by bit to ToField function ?
  uint8_t* w_field = malloc(lByte);
  //   bf8_t* bf_w = malloc(lByte);
  //   for (uint32_t i = 0; lByte < l; i++) {
  //     bf_w[i] = bf8_load(w[i]);
  //   }

  // Step: 2
  // TODO: Is to field, equivalent to the function bf_load ??
  // TODO: Are we really passing it bit by bit to ToField function ?
  // TODO: Unsure what to do here
  uint8_t v_field = malloc(lByte + lambdaByte);
  //   for (uint32_t i = 0; i < lByte + lambdaBytes; i++) {
  //   }

  // Step: 3
  uint8_t* in_0  = malloc(16);
  uint8_t* out_0 = malloc(16);
  memcpy(in_0, in, 16);
  memcpy(out_0, out, 16);

  // Step: 4
  // TODO: What is B ?
  uint32_t B     = 2;
  uint8_t* in_1  = malloc(16);
  uint8_t* out_1 = malloc(16);
  if (B == 2) {
    memcpy(in_1, in + 16, 16);
    memcpy(out_1, out + 16, 16);
  }

  // Step: 5
  // TODO: Sampling some element from the field ??
  uint8_t *a0, a1;

  // Step: 6
  // TODO: What is l_ke ??
  uint32_t lByte_ke      = 16;
  uint8_t* w_field_tilde = malloc(lByte_ke);
  memcpy(w_field_tilde, w_field, lByte_ke);
  uint8_t* v_field_tilde = malloc(lByte_ke);
  memcpy(v_field_tilde, v_field, lByte_ke);

  // Step: 7
  uint8_t MKey = 0;
  uint8_t* q;
  uint8_t* delta;
  uint8_t* a0_tilde;
  uint8_t* a1_tilde;
  uint8_t* k;
  uint8_t* vk;
  aes_key_schedule_constrains(lambda, w_field_tilde, v_field_tilde, MKey, q, delta, a0_tilde,
                              a1_tilde, k, vk);

  // Step: 8
  // TODO: resize it before memcpy
  memcpy(a0, a0_tilde, sizeof(a0_tilde));
  memcpy(a1, a1_tilde, sizeof(a1_tilde));

  // Step: 9
  // TODO
  uint32_t lByte_enc;
  memcpy(w_field_tilde, w_field + lByte_ke, lByte_enc);
  memcpy(v_field_tilde, v_field + lByte_ke, lByte_enc);

  // Step: 10
  // TODO
  uint8_t* qk;
  aes_cipher_constrains(lambda, in_0, out_0, w_field_tilde, v_field_tilde, k, vk, MKey, q, qk,
                        delta, a0_tilde, a1_tilde);

  // Step: 11
  // TODO
  memcpy(a0, a0_tilde, sizeof(a0_tilde));
  memcpy(a1, a1_tilde, sizeof(a1_tilde));

  // STep: 12..15
  // TODO
  if (B == 2) {
    memcpy(w_field_tilde, w_field + (lByte_ke + lByte_enc), (lByte) - (lByte_ke + lByte_enc));
    memcpy(v_field_tilde, v_field + (lByte_ke + lByte_enc), (lByte) - (lByte_ke + lByte_enc));
    aes_cipher_constrains(lambda, in_1, out_1, w_field_tilde, v_field_tilde, k, vk, MKey, q, qk,
                          delta, a0_tilde, a1_tilde);
  }

  // Step: 16..17
  for (uint32_t i = lByte; i < lByte + lambdaByte; i++) {
    // TODO: unsure what is happening here...
  }

  // STep: 18
  // TODO: unsure
  uint8_t* us;
  uint8_t* vs;

  // Step: 19..20
  // TODO
  uint8_t* a1_us_concat = malloc(sizeof(a1) + sizeof(us));
  uint8_t* a0_vs_concat = malloc(sizeof(a0) + sizeof(vs));
  switch (lambda) {
  case 256:
    zk_hash_256(a_tilde, r, chal, t, a1_us_concat, ell);
    zk_hash_256(b_tilde, r, chal, t, a0_vs_concat, ell);
    break;
  case 192:
    zk_hash_192(a_tilde, r, chal, t, a1_us_concat, ell);
    zk_hash_192(b_tilde, r, chal, t, a0_vs_concat, ell);
    break;
  default:
    zk_hash_128(a_tilde, r, chal, t, a1_us_concat, ell);
    zk_hash_128(b_tilde, r, chal, t, a0_vs_concat, ell);
    break;
  }
}

bool aes_verify(uint8_t* d, uint8_t* Q, uint8_t* chal_2, uint8_t* chal_3, uint8_t* a_tilde,
                uint8_t* b_tilde, uint8_t* in, uint8_t* out, uint32_t lambda, uint32_t tau,
                uint32_t l, uint32_t k0, uint32_t k1) {

  uint32_t lByte = l / 8;

  // Step: 1
  // TODO: confirm if this is a check or we actually assign it here
  uint32_t t0 = lambda % tau;
  uint32_t t1 = (lambda - (k0 * t0)) / k1;

  // STep: 2
  // TODO: field operation
  uint8_t* delta;

  // Step: 3
  uint8_t* in0 = malloc(16);
  memcpy(in0, in, 16);
  uint8_t* out0 = malloc(16);
  memcpy(out0, out, 16);

  // STep: 4
  // TODO: What is B ?
  uint32_t B = 2;
  uint8_t *in1, out1;
  if (B == 2) {
    in1  = malloc(16);
    out1 = malloc(16);
  }

  // Step: 5..9
  // TODO: from where comes the small b ?
  // TODO:
  uint8_t* b;
  uint32_t kb = 0;
  uint8_t** q;
  for (uint32_t i = 0; i < tau; i++) {
    if (i < t0) {
      b  = 0;
      kb = k0;
    } else {
      b  = 1;
      kb = k1;
    }
    uint8_t* chalout;
    ChalDec(chal_3, i, k0, t0, k1, t1, chalout);
    for (uint32_t j = 0; j < kb; j++) {
      // TODO: check it, very likely it is wrong, need to do bf multiplication I guess
      q[i][j] = q[i][j] ^ (chalout[j] * d[j]);
    }
  }
  // Step: 10..12
  for (uint32_t i = 0; i < lByte + tau; i++) {
    // TODO: ToField problem....
  }
  // TODO: sampling b

  // Step: 13
  // TODO: pass and check for null instead of 0 ?
  // TODO:
  uint32_t Mkey = 1;
  uint32_t lByte_ke;
  uint8_t* q_;
  uint8_t* b1;
  uint8_t* qk;
  memcpy(q_, q, lByte_ke);
  aes_key_schedule_constrains(lambda, NULL, NULL, Mkey, q_, delta, b1, qk);

  // Step: 14
  // TODO
  uint8_t* q_1;
  uint32_t lByte_enc;
  uint8_t* b2;
  memcpy(q_1, q + lByte_ke, lByte_enc);
  aes_cipher_constrains(lambda, in0, out0, NULL, NULL, NULL, NULL, Mkey, q_1, qk, delta, b2);

  // Step: 15..17
  // TODO
  if (B == 1) {
    memcpy(b, b1, sizeof(b1));
    memcpy(b + sizeof(b1), b2, sizeof(b2));
  } else {
    uint8_t* b3;
    uint8_t* q_2;
    memcpy(q_2, q + (lByte_ke + lByte_enc), lByte - (lByte_ke + lByte_enc));
    aes_cipher_constrains(lambda, in1, out1, NULL, NULL, NULL, NULL, Mkey, q_2, qk, delta, b3);
    memcpy(b, b1, sizeof(b1));
    memcpy(b + sizeof(b1), b2, sizeof(b2));
    memcpy(b + sizeof(b1) + sizeof(b2), b3, sizeof(b3));
  }

  // Step: 18
  // TODO : unsure what to do
  uint8_t* qs;

  // Step: 19
  // TODO
  uint8_t* q_tilde;
  uint8_t* b_qs_conat;
  uint8_t *r, t;
  uint32_t ell;
  switch (lambda) {
  case 256:
    zk_hash_256(q_tilde, r, chal_2, t, b_qs_conat, ell);
    break;
  case 192:
    zk_hash_192(q_tilde, r, chal_2, t, b_qs_conat, ell);
    break;
  default:
    zk_hash_128(q_tilde, r, chal_2, t, b_qs_conat, ell);
    break;
  }
  uint8_t* q_check = malloc(sizeof(a_tilde));
  // TODO: its bf operation i guess
  for (uint32_t i = 0; i < sizeof(q_check); i++) {
    q_check[i] = a_tilde[i] + (b_tilde[i] * delta[i]);
  }

  if (memcmp(q_tilde, q_check, sizeof(q_check)) == 0) {
    return true;
  } else {
    return false;
  }
}