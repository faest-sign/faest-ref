import sys
from pathlib import Path


def main():
    param_name = sys.argv[1]
    if param_name == "all":
        for bits in (128, 192, 256):
            for t in ("s", "f"):
                generate(f"{bits}{t}")
                generate(f"em_{bits}{t}")
    else:
        generate(param_name)


def generate(param_name):
    output_file = Path(f"faest_{param_name}") / "meson.build"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as outf:
        outf.write(
            f"""sources = configure_file(
    input: '../crypto_sign.c.in',
    output: 'crypto_sign.c',
    configuration: param_{param_name}
)
headers = [
    configure_file(
        input: '../crypto_sign.h.in',
        output: 'crypto_sign.h',
        configuration: param_{param_name}
    ),
    configure_file(
        input: '../crypto_sign.h.in',
        output: 'api.h',
        configuration: param_{param_name}
    )
]

libfaest_{param_name} = static_library('faest_{param_name}',
  sources,
  dependencies: libfaest_dependency,
  include_directories: include_directories,
  c_args: defines + c_flags
)
install_headers(headers, subdir: 'faest_{param_name}')
libfaest_{param_name}_dependency = declare_dependency(
  link_with: libfaest_{param_name},
  include_directories: include_directories
)
if openssl.found()
  tv_sources = files(join_paths(meson.project_source_root(), 'tools', 'rng.c'), join_paths(meson.project_source_root(), 'tools', 'PQCgenKAT_sign.cpp'))
  test_vector_generator = executable('feast_{param_name}_test_vectors', [sources] + faest_sources + tv_sources,
    dependencies: [openssl],
    include_directories: include_directories,
    c_args: defines + c_flags + ['-DHAVE_RANDOMBYTES'],
    cpp_args: defines + cpp_flags + ['-DHAVE_RANDOMBYTES']
  )
endif
if boost_program_options.found()
  bench_sources = files(join_paths(meson.project_source_root(), 'tools', 'bench.cpp'))
  bench = executable('feast_{param_name}_bench', [sources] + faest_sources + bench_sources,
    dependencies: [openssl, boost_program_options],
    include_directories: include_directories,
    c_args: defines + c_flags,
    cpp_args: defines + cpp_flags
  )
endif
"""
        )


if __name__ == "__main__":
    main()
