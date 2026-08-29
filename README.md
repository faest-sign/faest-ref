# FAEST - Reference implementation

This implementation of FAEST is the reference implementation of FAEST
specification version 3.0.

## Dependencies

For building:
* `meson` version 0.64 or newer
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