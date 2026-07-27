# FAEST - Reference implementation

## Dependencies

For building:
* `meson` version 0.57 or newer
* `ninja` (depending on the build system generator selected via `meson`)

For tests:
* `boost` (unit test framework)
* `NTL`

On Debian-based Linux distributions:
```sh
apt install meson ninja-build # for build dependencies
apt install libboost-test-dev libntl-dev # for test dependencies
```

Both `meson` and `ninja` are also available via PyPI:
```sh
pip install meson ninja
```

## Building

```sh
mkdir build
cd build
meson ..
ninja
ninja test
```

```
meson configure -Dtv-generators=disabled
./tv_binary
```


## Notes on Benchmarking

This implementation represents the reference implementation of FAEST. While it aims to be as
efficient as possible, its main goal is to provide an implementation that is mapped from the
specification. A optimized implementation using is available
[here](https://github.com/faest-sign/faest-arch-opt). For benchmarking and comparisons we recommend
the optimized implementation.
