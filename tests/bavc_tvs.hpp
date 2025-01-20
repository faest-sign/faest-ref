#ifndef VC_TVS_HPP
#define VC_TVS_HPP

#include <array>
#include <cstdint>

namespace bavc_tvs {
  namespace FAEST_128F {
    extern const std::array<uint8_t, 32> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 16> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_128F

  namespace FAEST_128S {
    extern const std::array<uint8_t, 32> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 11> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_128S

  namespace FAEST_192F {
    extern const std::array<uint8_t, 48> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 24> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_192F

  namespace FAEST_192S {
    extern const std::array<uint8_t, 48> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 16> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_192S

  namespace FAEST_256F {
    extern const std::array<uint8_t, 64> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 32> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_256F

  namespace FAEST_256S {
    extern const std::array<uint8_t, 64> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 22> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_256S

  namespace FAEST_EM_128F {
    extern const std::array<uint8_t, 32> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 16> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_EM_128F

  namespace FAEST_EM_128S {
    extern const std::array<uint8_t, 32> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 11> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_EM_128S

  namespace FAEST_EM_192F {
    extern const std::array<uint8_t, 48> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 24> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_EM_192F

  namespace FAEST_EM_192S {
    extern const std::array<uint8_t, 48> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 16> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_EM_192S

  namespace FAEST_EM_256F {
    extern const std::array<uint8_t, 64> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 32> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_EM_256F

  namespace FAEST_EM_256S {
    extern const std::array<uint8_t, 64> h;
    extern const std::array<uint8_t, 64> hashed_k;
    extern const std::array<uint8_t, 64> hashed_sd;
    extern const std::array<uint16_t, 22> i_delta;
    extern const std::array<uint8_t, 64> hashed_decom_i;
    extern const std::array<uint8_t, 64> hashed_rec_sd;
  } // namespace FAEST_EM_256S
} // namespace bavc_tvs

#endif