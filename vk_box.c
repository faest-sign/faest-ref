#include "vk_box.h"
#include "stdio.h"

const bf128_t* get_vk_128(vbb_t* vbb, unsigned int idx) {
  if (idx < FAEST_128F_LAMBDA) {
    const bf128_t* vk = get_V_k_128(vbb, idx);
    memcpy(vbb->vk_buf, vk, sizeof(bf128_t));
    return (bf128_t*)vbb->vk_buf;
  } else {
    unsigned int j = idx / 32 + FAEST_128F_Nwd;
    if ((j % FAEST_128F_Nwd) == 0 || (FAEST_128F_Nwd > 6 && (j % FAEST_128F_Nwd) == 4)) {
      unsigned int i_wd       = FAEST_128F_LAMBDA;
      unsigned int factor_128 = (idx / 128) - 1;
      unsigned int offset_128 = idx % 128;
      unsigned int index      = i_wd + factor_128 * 32 + offset_128;
      const bf128_t* vk       = get_V_k_128(vbb, index);
      memcpy(vbb->vk_buf, vk, sizeof(bf128_t));
      return (bf128_t*)vbb->vk_buf;
    } else {
      // Lhs recursive call
      const bf128_t* lhs_ptr = get_vk_128(vbb, idx - FAEST_128F_Nwd * 32);
      bf128_t lhs            = *lhs_ptr;
      // Rhs recursive call
      const bf128_t* rhs_ptr = get_vk_128(vbb, idx - 32);
      bf128_t rhs            = *rhs_ptr;

      bf128_t vk = bf128_add(lhs, rhs);
      memcpy(vbb->vk_buf, &vk, sizeof(bf128_t));
      return (bf128_t*)vbb->vk_buf;
    }
  }
}

const bf192_t* get_vk_192(vbb_t* vbb, unsigned int idx) {
  if (idx < FAEST_192F_LAMBDA) {
    const bf192_t* vk = get_V_k_192(vbb, idx);
    memcpy(vbb->vk_buf, vk, sizeof(bf192_t));
    return (bf192_t*)vbb->vk_buf;
  } else {
    unsigned int j = idx / 32 + FAEST_192F_Nwd;
    if ((j % FAEST_192F_Nwd) == 0 || (FAEST_192F_Nwd > 6 && (j % FAEST_192F_Nwd) == 4)) {
      unsigned int i_wd       = FAEST_192F_LAMBDA;
      unsigned int factor_192 = (idx / 192) - 1;
      unsigned int offset_192 = idx % 192;
      unsigned int index      = i_wd + factor_192 * 32 + offset_192;
      const bf192_t* vk       = get_V_k_192(vbb, index);
      memcpy(vbb->vk_buf, vk, sizeof(bf192_t));
      return (bf192_t*)vbb->vk_buf;
    } else {
      // Lhs recursive call
      const bf192_t* lhs_ptr = get_vk_192(vbb, idx - FAEST_192F_Nwd * 32);
      bf192_t lhs            = *lhs_ptr;
      // Rhs recursive call
      const bf192_t* rhs_ptr = get_vk_192(vbb, idx - 32);
      bf192_t rhs            = *rhs_ptr;

      bf192_t vk = bf192_add(lhs, rhs);
      memcpy(vbb->vk_buf, &vk, 3 * sizeof(uint64_t));
      return (bf192_t*)vbb->vk_buf;
    }
  }
}

const bf256_t* get_vk_256(vbb_t* vbb, unsigned int idx) {
  if (idx < FAEST_256F_LAMBDA) {
    const bf256_t* vk = get_V_k_256(vbb, idx);
    memcpy(vbb->vk_buf, vk, sizeof(bf256_t));
    return (bf256_t*)vbb->vk_buf;
  } else {
    // We go from j=N_k to j=4(R+1)
    // In our case this is j=8 to j=4(14+1)=60
    // Based on each j we have 32 tags in vk
    unsigned int j = idx / 32 + FAEST_256F_Nwd;
    if ((j % FAEST_256F_Nwd) == 0 || (FAEST_256F_Nwd > 6 && (j % FAEST_256F_Nwd) == 4)) {
      unsigned int i_wd       = FAEST_256F_LAMBDA;
      unsigned int factor_128 = (idx / 128) - 2;
      unsigned int offset_128 = idx % 128;
      unsigned int index      = i_wd + factor_128 * 32 + offset_128;
      const bf256_t* vk       = get_V_k_256(vbb, index);
      memcpy(vbb->vk_buf, vk, sizeof(bf256_t));
      return (bf256_t*)vbb->vk_buf;
    } else {
      // Lhs recursive call
      const bf256_t* lhs_ptr = get_vk_256(vbb, idx - FAEST_256F_Nwd * 32);
      bf256_t lhs            = *lhs_ptr;
      // Rhs recursive call
      const bf256_t* rhs_ptr = get_vk_256(vbb, idx - 32);
      bf256_t rhs            = *rhs_ptr;

      bf256_t vk = bf256_add(lhs, rhs);
      memcpy(vbb->vk_buf, &vk, sizeof(bf256_t));
      return (bf256_t*)vbb->vk_buf;
    }
  }
}