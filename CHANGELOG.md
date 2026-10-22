# Changelog

## Version 3.0.3 -- 2026-XX-XX

* Add SSE2 / AVX2+PCLMUL / NEON implementation of the finite field arithmetic (compile-time selection).
* Add ARM AES implementation (runtime dispatch).
* Add NEON-based matrix transpose (runtime dispatch).
* Add AES-NI implementation of EM-{192,256}{S,F} OWF (runtime dispatch).
* Fix OWF test vectors.
* Reduce memory consumption.

Benchmarks on a 2026 Intel laptop (Intel Core Ultra 7 265U) show the following improvements from 3.0.0 to 3.0.3:

| Parameter set | Keygen speedup | Sign: 3.0.0 → 3.0.3 (ms) | Speedup | Verify: 3.0.0 → 3.0.3 (ms) | Speedup |
|-------------- |---------------:|-------------------------:|--------:|---------------------------:|--------:|
| 128S          | 1.04×          | 36.135 → 6.016           | 6.01×   | 27.059 → 4.567             | 5.92×   |
| 128F          | 1.01×          | 16.423 → 0.988           | 16.62×  | 9.419 → 0.692              | 13.61×  |
| 192S          | 1.00×          | 143.625 → 19.328         | 7.43×   | 100.067 → 11.694           | 8.56×   |
| 192F          | 0.99×          | 79.716 → 2.673           | 29.82×  | 43.269 → 2.005             | 21.58×  |
| 256S          | 1.00×          | 214.451 → 25.423         | 8.44×   | 169.022 → 23.956           | 7.06×   |
| 256F          | 0.98×          | 102.520 → 4.185          | 24.50×  | 57.700 → 3.405             | 16.95×  |
| EM-128S       | 1.00×          | 15.710 → 4.968           | 3.16×   | 9.143 → 4.027              | 2.27×   |
| EM-128F       | 1.00×          | 11.920 → 0.987           | 12.08×  | 6.001 → 0.622              | 9.65×   |
| EM-192S       | 114.23×        | 62.194 → 12.633          | 4.92×   | 34.697 → 10.788            | 3.22×   |
| EM-192F       | 106.03×        | 51.405 → 2.038           | 25.22×  | 25.794 → 1.339             | 19.26×  |
| EM-256S       | 124.60×        | 98.327 → 20.663          | 4.76×   | 55.884 → 18.718            | 2.99×   |
| EM-256F       | 128.58×        | 82.686 → 3.257           | 25.39×  | 41.933 → 2.568             | 16.33×  |

## Version 3.0.2 -- 2026-10-08

* Fix shake header.

## Version 3.0.1 -- 2026-10-08

* Clear secrets from memory (#25).
* Update OQS integration. Thanks to Norman Ashley.
* Update SHAKE implementation from XKCP (#22).
* Remove unused code.

## Version 3.0.0 -- 2026-08-31

* Update to version 3.0 of FAEST specification (#18).

## Version 2.0.5 -- 2026-07-02

* Add option to use SHAKE implementation from OpenSSL.
* Refactor finite field implementations.

## Version 2.0.4 -- 2025-06-06

* Fix usage of bcrypt for AES192 and AES256.
* Refactor AES providers.

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
