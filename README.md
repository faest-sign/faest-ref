# FAEST - Reference implementation

This implementation of FAEST is the reference implementation of FAEST
specification version 3.0.

## Dependencies

For building:
* `meson` version 1 or newer
* `ninja` (depending on the build system generator selected via `meson`)

For tests (all optional):
* `boost` (unit test framework) for additional tests
* `NTL` for additional tests

On Debian-based Linux distributions:
```sh
apt install meson ninja-build # for build dependencies
apt install pkgconf libboost-test-dev libntl-dev # for test dependencies
```

Both `meson` and `ninja` are also available via PyPI:
```sh
pip install meson ninja
```

## Building

```sh
mkdir build
cd build
meson setup ..
ninja
ninja test
```

## Notes on Benchmarking

This implementation represents the reference implementation of FAEST. While it
aims to be as efficient as possible, its main goal is to provide an
implementation that is mapped from the specification. A optimized implementation
is available [here](https://github.com/faest-sign/faest-arch-opt). For
benchmarking and comparisons we recommend the optimized implementation.

For the reference implementation, all benchmarks can be executed via `bench.sh`.

# Implementations of AES and SHA3

For AES, the following implementations will be used (in order of preference and
availability):

* AES-NI (on x86 and x86-64)
* OpenSSL
* bcrypt (on Windows)
* unoptimized fallback implementation

For SHA3/SHAKE, the following implementation are available:

* `opt64` (optimized 64 bit implementation from the XKCP code package)
* `plain32` (32 bit implementation from the XKCP code package)
* `avx2` (AVX2 implementation from the XKCP code package)
* `armv8a-neon` (NEON implementation from the XKCP code package)
* `openssl` (for version >= 3.2)

The default on 64 bit architectures is `opt64` and on 32 bit architectures it is
`openssl` if available (and `plain32` oherwise). On aarch64 (if building under
gcc), the default implementation is `armv8a-neon`. On s390x, the default
implementation is `openssl` if available.

The SHA3 implementation does not peform runtime dispatching based on the
available CPU features. If, however, the target is known to support AVX2, `avx2`
is the recommended implementation on x86-64.