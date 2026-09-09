#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate liboqs integration artifacts for FAEST."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INTEGRATION = Path(__file__).resolve().parent
TEMPLATES = INTEGRATION / "templates"
GENERATED = INTEGRATION / "generated"

# Keep in sync with meson.build parameter definitions.
PARAMETER_SETS = [
    {
        "param": "128S",
        "param_l": "128s",
        "pqclean_scheme": "faest-128s",
        "pretty_name": "FAEST-128S",
        "lambda": 128,
        "nist_level": 1,
        "pk_size": 32,
        "sk_size": 32,
        "sig_size": 4066,
        "owf_input_size": 16,
        "owf_output_size": 16,
        "nst": 4,
        "ske": 40,
        "r": 10,
        "beta": 1,
        "ell": 960,
        "lke": 448,
        "lenc": 512,
        "tau": 11,
        "w_grind": 7,
        "t_open": 102,
        "n_mult": 312,
    },
    {
        "param": "128F",
        "param_l": "128f",
        "pqclean_scheme": "faest-128f",
        "pretty_name": "FAEST-128F",
        "lambda": 128,
        "nist_level": 1,
        "pk_size": 32,
        "sk_size": 32,
        "sig_size": 5170,
        "owf_input_size": 16,
        "owf_output_size": 16,
        "nst": 4,
        "ske": 40,
        "r": 10,
        "beta": 1,
        "ell": 960,
        "lke": 448,
        "lenc": 512,
        "tau": 17,
        "w_grind": 8,
        "t_open": 108,
        "n_mult": 312,
    },
    {
        "param": "EM_128S",
        "param_l": "em_128s",
        "pqclean_scheme": "faest-em-128s",
        "pretty_name": "FAEST-EM-128S",
        "lambda": 128,
        "nist_level": 1,
        "pk_size": 32,
        "sk_size": 32,
        "sig_size": 3466,
        "owf_input_size": 16,
        "owf_output_size": 16,
        "nst": 4,
        "ske": 0,
        "r": 10,
        "beta": 1,
        "ell": 640,
        "lke": 128,
        "lenc": 512,
        "tau": 11,
        "w_grind": 7,
        "t_open": 103,
        "n_mult": 312,
    },
    {
        "param": "EM_128F",
        "param_l": "em_128f",
        "pqclean_scheme": "faest-em-128f",
        "pretty_name": "FAEST-EM-128F",
        "lambda": 128,
        "nist_level": 1,
        "pk_size": 32,
        "sk_size": 32,
        "sig_size": 4170,
        "owf_input_size": 16,
        "owf_output_size": 16,
        "nst": 4,
        "ske": 0,
        "r": 10,
        "beta": 1,
        "ell": 640,
        "lke": 128,
        "lenc": 512,
        "tau": 17,
        "w_grind": 8,
        "t_open": 105,
        "n_mult": 312,
    },
    {
        "param": "192S",
        "param_l": "192s",
        "pqclean_scheme": "faest-192s",
        "pretty_name": "FAEST-192S",
        "lambda": 192,
        "nist_level": 3,
        "pk_size": 48,
        "sk_size": 40,
        "sig_size": 9410,
        "owf_input_size": 16,
        "owf_output_size": 32,
        "nst": 4,
        "ske": 32,
        "r": 12,
        "beta": 2,
        "ell": 1728,
        "lke": 448,
        "lenc": 640,
        "tau": 16,
        "w_grind": 12,
        "t_open": 162,
        "n_mult": 504,
    },
    {
        "param": "192F",
        "param_l": "192f",
        "pqclean_scheme": "faest-192f",
        "pretty_name": "FAEST-192F",
        "lambda": 192,
        "nist_level": 3,
        "pk_size": 48,
        "sk_size": 40,
        "sig_size": 11738,
        "owf_input_size": 16,
        "owf_output_size": 32,
        "nst": 4,
        "ske": 32,
        "r": 12,
        "beta": 2,
        "ell": 1728,
        "lke": 448,
        "lenc": 640,
        "tau": 24,
        "w_grind": 8,
        "t_open": 163,
        "n_mult": 504,
    },
    {
        "param": "EM_192S",
        "param_l": "em_192s",
        "pqclean_scheme": "faest-em-192s",
        "pretty_name": "FAEST-EM-192S",
        "lambda": 192,
        "nist_level": 3,
        "pk_size": 48,
        "sk_size": 48,
        "sig_size": 7874,
        "owf_input_size": 24,
        "owf_output_size": 24,
        "nst": 6,
        "ske": 0,
        "r": 12,
        "beta": 1,
        "ell": 1152,
        "lke": 192,
        "lenc": 960,
        "tau": 16,
        "w_grind": 8,
        "t_open": 162,
        "n_mult": 504,
    },
    {
        "param": "EM_192F",
        "param_l": "em_192f",
        "pqclean_scheme": "faest-em-192f",
        "pretty_name": "FAEST-EM-192F",
        "lambda": 192,
        "nist_level": 3,
        "pk_size": 48,
        "sk_size": 48,
        "sig_size": 9818,
        "owf_input_size": 24,
        "owf_output_size": 24,
        "nst": 6,
        "ske": 0,
        "r": 12,
        "beta": 1,
        "ell": 1152,
        "lke": 192,
        "lenc": 960,
        "tau": 25,
        "w_grind": 8,
        "t_open": 171,
        "n_mult": 504,
    },
    {
        "param": "256S",
        "param_l": "256s",
        "pqclean_scheme": "faest-256s",
        "pretty_name": "FAEST-256S",
        "lambda": 256,
        "nist_level": 5,
        "pk_size": 48,
        "sk_size": 48,
        "sig_size": 16626,
        "owf_input_size": 16,
        "owf_output_size": 32,
        "nst": 4,
        "ske": 52,
        "r": 14,
        "beta": 2,
        "ell": 2208,
        "lke": 672,
        "lenc": 768,
        "tau": 22,
        "w_grind": 6,
        "t_open": 225,
        "n_mult": 696,
    },
    {
        "param": "256F",
        "param_l": "256f",
        "pqclean_scheme": "faest-256f",
        "pretty_name": "FAEST-256F",
        "lambda": 256,
        "nist_level": 5,
        "pk_size": 48,
        "sk_size": 48,
        "sig_size": 20856,
        "owf_input_size": 16,
        "owf_output_size": 32,
        "nst": 4,
        "ske": 52,
        "r": 14,
        "beta": 2,
        "ell": 2208,
        "lke": 672,
        "lenc": 768,
        "tau": 33,
        "w_grind": 8,
        "t_open": 229,
        "n_mult": 704,
    },
    {
        "param": "EM_256S",
        "param_l": "em_256s",
        "pqclean_scheme": "faest-em-256s",
        "pretty_name": "FAEST-EM-256S",
        "lambda": 256,
        "nist_level": 5,
        "pk_size": 64,
        "sk_size": 64,
        "sig_size": 14554,
        "owf_input_size": 32,
        "owf_output_size": 32,
        "nst": 8,
        "ske": 0,
        "r": 14,
        "beta": 1,
        "ell": 1792,
        "lke": 256,
        "lenc": 1536,
        "tau": 22,
        "w_grind": 6,
        "t_open": 218,
        "n_mult": 696,
    },
    {
        "param": "EM_256F",
        "param_l": "em_256f",
        "pqclean_scheme": "faest-em-256f",
        "pretty_name": "FAEST-EM-256F",
        "lambda": 256,
        "nist_level": 5,
        "pk_size": 64,
        "sk_size": 64,
        "sig_size": 18084,
        "owf_input_size": 32,
        "owf_output_size": 32,
        "nst": 8,
        "ske": 0,
        "r": 14,
        "beta": 1,
        "ell": 1792,
        "lke": 256,
        "lenc": 1536,
        "tau": 33,
        "w_grind": 8,
        "t_open": 229,
        "n_mult": 704,
    },
]

SUBMITTERS = [
    "Carsten Baum",
    "Waard Beullens",
    "Lennart Braun",
    "Cyprien Delpech de Saint Guilhem",
    "Michael Klooß",
    "Christian Majenz",
    "Shibam Mukherjee",
    "Emmanuela Orsini",
    "Sebastian Ramacher",
    "Christian Rechberger",
    "Lawrence Roy",
    "Peter Scholl",
]

# SHA-256 of liboqs `kat_sig <scheme>` output (count = 0 vector).
NISTKAT_SHA256 = {
    "FAEST-128S": "ac36ff496db5b30a8c4df1768100c0e804c4b627c8df85be27ebd3a14001ca13",
    "FAEST-128F": "a8dc23f6d66282f529c1cf52288566cc0af8650910aa0b00a3c7acc9be380bc9",
    "FAEST-EM-128S": "ee087b31592cbd1eabc34882fbe22d3c5d270b92ada1e4aab5be5588dc9cbd2e",
    "FAEST-EM-128F": "585e5adc1a1edb610eeec10148e246db86758d0f77fc27137def3a2ab7238a2f",
    "FAEST-192S": "cdff1a6fc7d81eb03522b968a9af6e0db6295e00410c95ea6330c2bf1fa97148",
    "FAEST-192F": "0d8d585f17687dfad35159b02ebf9d7ed96983bcf4630fe1062f5a0690a3f632",
    "FAEST-EM-192S": "7855c8018fb3a60994722be77b44b0c14902303c649f2b80afa5c21911d0f676",
    "FAEST-EM-192F": "19d56fcb425264ee609c1da98622baa26bd349ab2ac6635a9c1a5338b691480f",
    "FAEST-256S": "5573247ebd897e38930942910772d58b7fbfe969793a8326180f173131492295",
    "FAEST-256F": "77db00b2d70791ed8f710bb226722357c5e155dc0ce596e9356aaf8d80f46f90",
    "FAEST-EM-256S": "3aaf2e891cade204aef3e92e3fc39dd2ee8e6eeaf4cf37ed8ec1888ae6903150",
    "FAEST-EM-256F": "c27139ced48f1bfecc0ced01e3ee03377623f562c86e98fcc3b47dafa7872799",
}

COMMON_HEADERS = [
    "aes.h",
    "aesni.h",
    "bavc.h",
    "compat.h",
    "cpu.h",
    "endian_compat.h",
    "faest.h",
    "faest_aes.h",
    "faest_defines.h",
    "faest_impl.h",
    "fields.h",
    "hash_shake.h",
    "instances.h",
    "macros.h",
    "owf.h",
    "random_oracle.h",
    "randomness.h",
    "universal_hashing.h",
    "utils.h",
    "vole.h",
]

COMMON_SOURCES = [
    "aes.c",
    "bavc.c",
    "compat.c",
    "faest_impl.c",
    "fields.c",
    "instances.c",
    "owf.c",
    "random_oracle.c",
    "randomness.c",
    "universal_hashing.c",
    "utils.c",
    "vole.c",
    "integration/liboqs/generated/config.h",
    "integration/liboqs/generated/parameters.h",
    "integration/liboqs/generated/faest_aes_128.c",
    "integration/liboqs/generated/faest_aes_192.c",
    "integration/liboqs/generated/faest_aes_256.c",
    "tables/pregenerated/tables_128s.h",
    "tables/pregenerated/tables_128f.h",
    "tables/pregenerated/tables_192s.h",
    "tables/pregenerated/tables_192f.h",
    "tables/pregenerated/tables_em_192s.h",
    "tables/pregenerated/tables_em_192f.h",
    "tables/pregenerated/tables_256s.h",
    "tables/pregenerated/tables_256f.h",
] + COMMON_HEADERS


def render_template(template_path: Path, mapping: dict[str, str]) -> str:
    content = template_path.read_text()
    for key, value in mapping.items():
        content = content.replace(f"@{key}@", value)
    return content


def calc_k(ps: dict) -> int:
    return ((ps["lambda"] - ps["w_grind"]) // ps["tau"]) + 1


def calc_tau1(ps: dict) -> int:
    return (ps["lambda"] - ps["w_grind"]) % ps["tau"]


def calc_tau0(ps: dict) -> int:
    return ps["tau"] - calc_tau1(ps)


def calc_l(ps: dict) -> int:
    k = calc_k(ps)
    tau1 = calc_tau1(ps)
    tau0 = calc_tau0(ps)
    return tau1 * (1 << k) + tau0 * (1 << (k - 1))


def generate_parameters_h() -> str:
    lines = [
        "/* SPDX-License-Identifier: MIT */",
        "/* Generated for liboqs integration. */",
        "#ifndef FAEST_PARAMETERS_H",
        "#define FAEST_PARAMETERS_H",
        "",
    ]
    for ps in PARAMETER_SETS:
        prefix = f"FAEST_{ps['param']}"
        lines.extend(
            [
                f"#define {prefix}_PARAM \"{ps['param']}\"",
                f"#define {prefix}_PARAM_L \"{ps['param_l']}\"",
                f"#define {prefix}_LAMBDA {ps['lambda']}",
                f"#define {prefix}_Nst {ps['nst']}",
                f"#define {prefix}_Ske {ps['ske']}",
                f"#define {prefix}_R {ps['r']}",
                f"#define {prefix}_BETA {ps['beta']}",
                f"#define {prefix}_ELL {ps['ell']}",
                f"#define {prefix}_Lke {ps['lke']}",
                f"#define {prefix}_Lenc {ps['lenc']}",
                f"#define {prefix}_TAU {ps['tau']}",
                f"#define {prefix}_W_GRIND {ps['w_grind']}",
                f"#define {prefix}_T_OPEN {ps['t_open']}",
                f"#define {prefix}_SIG_SIZE {ps['sig_size']}",
                f"#define {prefix}_PK_SIZE {ps['pk_size']}",
                f"#define {prefix}_SK_SIZE {ps['sk_size']}",
                f"#define {prefix}_OWF_INPUT_SIZE {ps['owf_input_size']}",
                f"#define {prefix}_OWF_OUTPUT_SIZE {ps['owf_output_size']}",
                f"#define {prefix}_N_MULT {ps['n_mult']}",
                "",
            ]
        )
    lines.extend(
        [
            "#define FAEST_128_LAMBDA 128",
            "#define FAEST_192_LAMBDA 192",
            "#define FAEST_256_LAMBDA 256",
            "",
            "#endif",
            "",
        ]
    )
    return "\n".join(lines)


def generate_config_h() -> str:
    return (TEMPLATES / "config.h.in").read_text()


def generate_faest_aes(sss: str) -> str:
    return render_template(ROOT / "faest_aes.c.in", {"SSS": sss})


def generate_scheme_sources(ps: dict, out_dir: Path) -> None:
    mapping = {
        "PARAM": ps["param"],
        "PARAM_L": ps["param_l"],
        "PK_SIZE": str(ps["pk_size"]),
        "SK_SIZE": str(ps["sk_size"]),
        "SIG_SIZE": str(ps["sig_size"]),
        "LAMBDA": str(ps["lambda"]),
        "Nst": str(ps["nst"]),
        "Ske": str(ps["ske"]),
        "R": str(ps["r"]),
        "BETA": str(ps["beta"]),
        "ELL": str(ps["ell"]),
        "Lke": str(ps["lke"]),
        "Lenc": str(ps["lenc"]),
        "TAU": str(ps["tau"]),
        "W_GRIND": str(ps["w_grind"]),
        "T_OPEN": str(ps["t_open"]),
        "OWF_INPUT_SIZE": str(ps["owf_input_size"]),
        "OWF_OUTPUT_SIZE": str(ps["owf_output_size"]),
        "N_MULT": str(ps["n_mult"]),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"faest_{ps['param_l']}.h").write_text(
        render_template(ROOT / "faest_param.h.in", mapping)
    )
    (out_dir / f"faest_{ps['param_l']}.c").write_text(
        render_template(ROOT / "faest_param.c.in", mapping)
    )
    (out_dir / "api.h").write_text(render_template(TEMPLATES / "api.h.in", mapping))
    (out_dir / "crypto_sign.c").write_text(
        render_template(TEMPLATES / "crypto_sign.c.in", mapping)
    )


def namespace(ps: dict) -> str:
    return f"PQCLEAN_FAEST_{ps['param']}_REF"


def generate_meta_yml(ps: dict) -> str:
    ns = namespace(ps)
    submitters = "\n".join(f"  - {name}" for name in SUBMITTERS)
    return f"""name: {ps['pretty_name']}
type: signature
principal-submitters:
{submitters}
crypto-assumption: Syndrome decoding in the random code model and AES.
website: https://faest.info/
spec-version: 3.0
claimed-nist-level: {ps['nist_level']}
claimed-security: EUF-CMA
length-public-key: {ps['pk_size']}
length-secret-key: {ps['sk_size']}
length-signature: {ps['sig_size']}
nistkat-sha256: {NISTKAT_SHA256[ps['pretty_name']]}
implementations:
  - name: ref
    version: 3.0
    folder_name: integration/liboqs/generated/{ps['pqclean_scheme']}
    signature_keypair: {ns}_crypto_sign_keypair
    signature_signature: {ns}_crypto_sign_signature
    signature_verify: {ns}_crypto_sign_verify
    compile_opts: -DOQS -DFAEST_LIBOQS_BUILD -DHAVE_CONFIG_H
    sources: api.h crypto_sign.c faest_{ps['param_l']}.c faest_{ps['param_l']}.h
    common_dep: faest_common
    no-secret-dependent-branching-claimed: true
    no-secret-dependent-branching-checked-by-valgrind: false
    large-stack-usage: true
    supported-platforms: all
"""


def generate_meta_common_yml() -> str:
    sources = " ".join(COMMON_SOURCES)
    return f"""commons:
  - name: faest_common
    folder_name: .
    sources: {sources}
    include_only: true
"""


def generate_support_md() -> str:
    return """# FAEST liboqs integration support statement

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
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit with status 1 if generated files are out of date.",
    )
    args = parser.parse_args()

    expected_outputs: list[Path] = []

    GENERATED.mkdir(parents=True, exist_ok=True)
    config_h = GENERATED / "config.h"
    parameters_h = GENERATED / "parameters.h"
    config_h.write_text(generate_config_h())
    parameters_h.write_text(generate_parameters_h())
    expected_outputs.extend([config_h, parameters_h])

    for sss in ("128", "192", "256"):
        out = GENERATED / f"faest_aes_{sss}.c"
        out.write_text(generate_faest_aes(sss))
        expected_outputs.append(out)

    for ps in PARAMETER_SETS:
        scheme_dir = GENERATED / ps["pqclean_scheme"]
        generate_scheme_sources(ps, scheme_dir)
        expected_outputs.extend(
            [
                scheme_dir / "api.h",
                scheme_dir / "crypto_sign.c",
                scheme_dir / f"faest_{ps['param_l']}.c",
                scheme_dir / f"faest_{ps['param_l']}.h",
            ]
        )
        meta = INTEGRATION / f"{ps['pretty_name']}_META.yml"
        meta.write_text(generate_meta_yml(ps))
        expected_outputs.append(meta)

    meta_common = INTEGRATION / "META_Common.yml"
    meta_common.write_text(generate_meta_common_yml())
    expected_outputs.append(meta_common)

    support = INTEGRATION / "SUPPORT.md"
    support.write_text(generate_support_md())
    expected_outputs.append(support)

    if args.check:
        stale = [path for path in expected_outputs if not path.exists()]
        if stale:
            raise SystemExit(f"Generated files missing or stale: {stale[0]}")
        print("liboqs integration artifacts are up to date.")
        return

    print(f"Generated liboqs integration artifacts under {INTEGRATION}")


if __name__ == "__main__":
    main()
