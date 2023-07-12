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
  dependencies: libfaest_static_dependency,
  include_directories: include_directories,
  c_args: defines + c_flags
)
install_headers(headers, subdir: 'faest_{param_name}')
libfaest_{param_name}_dependency = declare_dependency(
  link_with: libfaest_{param_name},
  include_directories: include_directories
)
if openssl.found()
  tv_sources = files(
    join_paths(meson.project_source_root(), 'randomness.c'),
    join_paths(meson.project_source_root(), 'tools', 'rng.c'),
    join_paths(meson.project_source_root(), 'tools', 'PQCgenKAT_sign.cpp')
  )
  test_vector_generator = executable('faest_{param_name}_test_vectors', [sources] + tv_sources,
    dependencies: [libfaest_no_random_static_dependency, openssl],
    include_directories: include_directories,
    c_args: defines + c_flags + ['-DHAVE_RANDOMBYTES'],
    cpp_args: defines + cpp_flags + ['-DHAVE_RANDOMBYTES'],
    override_options: ['b_lto=false'],
  )
endif
if boost_program_options.found()
  bench_sources = files(join_paths(meson.project_source_root(), 'tools', 'bench.cpp'))
  bench = executable('faest_{param_name}_bench', bench_sources,
    dependencies: [libfaest_{param_name}_dependency, boost_program_options],
    include_directories: include_directories,
    c_args: defines + c_flags,
    cpp_args: defines + cpp_flags
  )
endif
test_sources = files(join_paths(meson.project_source_root(), 'tests', 'api_test.c'))
faest_{param_name}_test = executable('faest_{param_name}_api_test', test_sources,
  dependencies: [libfaest_{param_name}_dependency, valgrind],
  include_directories: include_directories,
  c_args: defines + c_flags + valgrind_defines,
  override_options: ['b_lto=false'],
)
test('faest_{param_name}_api_test', faest_{param_name}_test,
  timeout: 6000,
)
if valgrind.found() and valgrind_exec.found()
  test('faest_{param_name}_api_test_ct', valgrind_exec,
    args: ['-q', '--error-exitcode=1', '--track-origins=yes', faest_{param_name}_test],
    timeout: 6000,
  )
endif
if valgrind_exec.found()
  custom_target('faest_{param_name}_memory_usage',
    command: [valgrind_exec, '-q', '--error-exitcode=1', '--tool=massif', '--stacks=yes', '--massif-out-file=@OUTPUT@', faest_{param_name}_test],
    output: 'faest_{param_name}.massif',
    depends: [faest_{param_name}_test],
    install: false,
    build_always_stale: true,
    build_by_default: false,
  )
endif
if get_option('benchmarks').enabled()
  bench_sources = files(
    join_paths(meson.project_source_root(), 'tools', 'bench_c2.cpp'),
  )
  bench_catch = executable('faest_{param_name}_bench_c2', bench_sources,
    dependencies: [libfaest_{param_name}_dependency, boost_program_options, catch2],
    include_directories: include_directories,
    c_args: defines + c_flags,
    cpp_args: defines + cpp_flags
  )
endif
"""
        )


if __name__ == "__main__":
    main()
