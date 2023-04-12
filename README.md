# FAEST - Reference implementation

## Dependencies

For building:
* `meson` version 0.56 or newer
* `ninja` (depending on the build system generator selected via `meson`)

For tests:
* `boost` (unit test framework)
* `NTL`

On Debian-based linux distributions:
```sh
apt install meson ninja-build # for build dependencies
apt install libboost-test-dev libntl-dev # for test dependencies
```

## Building

```sh
mkdir build
cd build
meson ..
ninja
ninja test
```
