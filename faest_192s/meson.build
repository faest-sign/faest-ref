sources = configure_file(
    input: '../crypto_sign.c.in',
    output: 'crypto_sign.c',
    configuration: param_192s
)
headers = [
    configure_file(
        input: '../crypto_sign.h.in',
        output: 'crypto_sign.h',
        configuration: param_192s
    ),
    configure_file(
        input: '../crypto_sign.h.in',
        output: 'api.h',
        configuration: param_192s
    )
]

libfaest_192s = static_library('faest_192s',
  sources,
  dependencies: libfaest_static_dependency,
  include_directories: include_directories,
  c_args: defines + c_flags
)
install_headers(headers, subdir: 'faest_192s')
libfaest_192s_dependency = declare_dependency(
  link_with: libfaest_192s,
  include_directories: include_directories
)
if openssl.found()
  tv_sources = files(join_paths(meson.project_source_root(), 'tools', 'rng.c'), join_paths(meson.project_source_root(), 'tools', 'PQCgenKAT_sign.cpp'))
  test_vector_generator = executable('feast_192s_test_vectors', [sources] + faest_sources + tv_sources,
    dependencies: [openssl],
    include_directories: include_directories,
    c_args: defines + c_flags + ['-DHAVE_RANDOMBYTES'],
    cpp_args: defines + cpp_flags + ['-DHAVE_RANDOMBYTES']
  )
endif
if boost_program_options.found()
  bench_sources = files(join_paths(meson.project_source_root(), 'tools', 'bench.cpp'))
  bench = executable('feast_192s_bench', bench_sources,
    dependencies: [libfaest_192s_dependency, boost_program_options],
    include_directories: include_directories,
    c_args: defines + c_flags,
    cpp_args: defines + cpp_flags
  )
endif
test_sources = files(join_paths(meson.project_source_root(), 'tests', 'api_test.c'))
faest_192s_test = executable('faest_192s_api_test', test_sources,
  dependencies: [libfaest_192s_dependency],
  include_directories: include_directories,
  c_args: defines + c_flags,
)
test('faest_192s_api_test', faest_192s_test)
