# FAEST liboqs integration support statement

This directory contains the upstream integration glue for importing FAEST into
[liboqs](https://github.com/open-quantum-safe/liboqs) via
`scripts/copy_from_upstream`.

## Maintenance

The FAEST reference implementation is actively maintained. Integration artifacts
in this directory are regenerated with:

```sh
python3 integration/liboqs/generate.py
```

Run this script after changing parameter sets in `meson.build` or the upstream
`.in` templates.

## liboqs primitive routing

When built for liboqs (`-DOQS`), FAEST routes cryptographic primitives to
liboqs:

| Primitive | Routing |
|-----------|---------|
| Random bytes | `OQS_randombytes` (`randomness.c`) |
| SHAKE | `OQS_SHA3_shake*` (`hash_shake.h`) |
| Memory helpers | `OQS_MEM_*` (`compat.h`) |
| CPU features | `OQS_CPU_has_extension` (`cpu.h`) |
| AES-128/256 ECB | `OQS_AES128/256_ECB_*` (`aes.c`) |
| AES-192 ECB | Internal FAEST implementation (`aes.c`) |

## Known follow-ups before merging into liboqs

1. Add optimized implementations (`avx2`, `neon`) as additional entries in each
   scheme's `META.yml`.
2. Run liboqs constant-time tests and document results under
   `tests/constant_time/sig/`.
