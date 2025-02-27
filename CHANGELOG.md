# Changelog

## Version 2.0.3 -- 2025-02-28

* Fix stack smashing during OWF computation with OpenSSL.

## Version 2.0.2 -- 2025-02-25

* Add AES-NI-based implementation for PRNG.
* Small efficiency improvements.

## Version 2.0.1 -- 2025-02-21

* Reduce memory usage.
* Remove conversion of bit-packed storage to 1 bit/byte.
* Precompute finite field constants.
* Generate prover and verifier implementations from template.

## Version 2.0.0 -- 2025-02-08

* Update to version 2.0 of FAEST specification.

## Version 1.0.2 -- 2024-12-19

* Extend test coverage.
* Use compiler's support (GCC's vector size attribute) for vectorized implementations.
* Reduce memory consumption.
* Refactor for more efficient code generation.

## Version 1.0.1 -- 2023-10-02

* Post NIST PQC submission clean up.
* Reduce memory consumption when processing signatures.
* Implement more bit-sliced algorithms.
* Replace manual copies with `memcpy`.

## Version 1.0.0 -- 2023-06-02

* Initial release.
* Version submitted to the NIST PQC project.
