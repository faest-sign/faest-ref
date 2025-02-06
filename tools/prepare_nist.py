# SPDX-License-Identifier: MIT

import sys
import shutil
import subprocess
import os
from pathlib import Path


def generate(
    project_root: Path, build_root: Path, target_root: Path, param_name: str
) -> None:
    target = (target_root / "Reference_Implementation" / param_name).absolute()
    target_kat = (target_root / "KAT" / param_name).absolute()
    print(
        f"Preparing {param_name}: root: {project_root}, build root: {build_root}, target: {target}"
    )
    target.mkdir(parents=True, exist_ok=True)
    target_kat.mkdir(parents=True, exist_ok=True)

    target_sha3 = target / "sha3"
    target_sha3.mkdir(parents=True, exist_ok=True)
    target_nist_kat = target / "NIST-KATs"
    target_nist_kat.mkdir(parents=True, exist_ok=True)
    target_tests = target / "tests"
    target_tests.mkdir(parents=True, exist_ok=True)

    sha3_sources = project_root / "sha3"
    test_sources = project_root / "tests"
    tools_sources = project_root / "tools"

    # copy FAEST implementation
    for source in project_root.glob("*.c"):
        shutil.copy(source, target)
    for header in project_root.glob("*.h"):
        shutil.copy(header, target)

    # copy generated files
    for build_source in (f"{param_name}.c", f"{param_name}.h", "parameters.h"):
        shutil.copy(build_root / build_source, target)
    build_param = build_root / f"{param_name}"
    for build_source in ("crypto_sign.c", "crypto_sign.h", "api.h"):
        shutil.copy(build_param / build_source, target)

    # copy sha3 sources
    for source in sha3_sources.glob("*.c"):
        shutil.copy(source, target_sha3)
    for header in sha3_sources.glob("*.h"):
        shutil.copy(header, target_sha3)
    for source in sha3_sources.glob("*.macros"):
        shutil.copy(source, target_sha3)
    for source in sha3_sources.glob("*.inc"):
        shutil.copy(source, target_sha3)
    sha3_sources = sha3_sources / "opt64"
    for source in sha3_sources.glob("*.c"):
        shutil.copy(source, target_sha3)
    for header in sha3_sources.glob("*.h"):
        shutil.copy(header, target_sha3)
    for source in sha3_sources.glob("*.macros"):
        shutil.copy(source, target_sha3)
    for source in sha3_sources.glob("*.inc"):
        shutil.copy(source, target_sha3)

    # copy tests
    for test_source in ("api_test.c",):
        shutil.copy(test_sources / test_source, target_tests)
    # copy NIST files
    for tool_source in ("rng.c", "rng.h", "PQCgenKAT_sign.c"):
        shutil.copy(tools_sources / tool_source, target_nist_kat)
    for tool_source in ("Makefile",):
        shutil.copy(tools_sources / tool_source, target)

    # build and create KATs
    print(f"Building {param_name}")
    cpu_count = os.cpu_count()
    subprocess.check_call(
        ["make"] if cpu_count is None else ["make", "-j", str(max(2, cpu_count - 1))],
        cwd=target,
    )
    print(f"Generating KATs for {param_name}")
    subprocess.check_call(target_nist_kat / "PQCgenKAT_sign", cwd=target_kat)
    subprocess.check_call(["make", "clean"], cwd=target)


def main():
    project_root = Path(sys.argv[1])
    build_root = Path(sys.argv[2])
    target_root = Path(sys.argv[3])
    param_names = sys.argv[4:]

    for param_name in param_names:
        generate(project_root, build_root, target_root, param_name)


if __name__ == "__main__":
    main()
