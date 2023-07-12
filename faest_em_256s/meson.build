sources = configure_file(
    input: '../crypto_sign.c.in',
    output: 'crypto_sign.c',
    configuration: param_em_256s
)
headers = [
    configure_file(
        input: '../crypto_sign.h.in',
        output: 'crypto_sign.h',
        configuration: param_em_256s
    ),
    configure_file(
        input: '../crypto_sign.h.in',
        output: 'api.h',
        configuration: param_em_256s
    )
]

libfaest_em_256s = static_library('faest_em_256s',
  sources,
  dependencies: libfaest_static_dependency,
  include_directories: include_directories,
  c_args: defines + c_flags
)
install_headers(headers, subdir: 'faest_em_256s')
libfaest_em_256s_dependency = declare_dependency(
  link_with: libfaest_em_256s,
  include_directories: include_directories
)
if openssl.found()
  tv_sources = files(
    join_paths(meson.project_source_root(), 'randomness.c'),
    join_paths(meson.project_source_root(), 'tools', 'rng.c'),
    join_paths(meson.project_source_root(), 'tools', 'PQCgenKAT_sign.cpp')
  )
  test_vector_generator = executable('faest_em_256s_test_vectors', [sources] + tv_sources,
    dependencies: [libfaest_no_random_static_dependency, openssl],
    include_directories: include_directories,
    c_args: defines + c_flags + ['-DHAVE_RANDOMBYTES'],
    cpp_args: defines + cpp_flags + ['-DHAVE_RANDOMBYTES'],
    override_options: ['b_lto=false'],
  )
endif
if boost_program_options.found()
  bench_sources = files(join_paths(meson.project_source_root(), 'tools', 'bench.cpp'))
  bench = executable('faest_em_256s_bench', bench_sources,
    dependencies: [libfaest_em_256s_dependency, boost_program_options],
    include_directories: include_directories,
    c_args: defines + c_flags,
    cpp_args: defines + cpp_flags
  )
endif
test_sources = files(join_paths(meson.project_source_root(), 'tests', 'api_test.c'))
faest_em_256s_test = executable('faest_em_256s_api_test', test_sources,
  dependencies: [libfaest_em_256s_dependency, valgrind],
  include_directories: include_directories,
  c_args: defines + c_flags + valgrind_defines,
  override_options: ['b_lto=false'],
)
test('faest_em_256s_api_test', faest_em_256s_test,
  timeout: 6000,
)
if valgrind.found() and valgrind_exec.found()
  test('faest_em_256s_api_test_ct', valgrind_exec,
    args: ['-q', '--error-exitcode=1', '--track-origins=yes', faest_em_256s_test],
    timeout: 6000,
  )
endif
if valgrind_exec.found()
  custom_target('faest_em_256s_memory_usage',
    command: [valgrind_exec, '-q', '--error-exitcode=1', '--tool=massif', '--stacks=yes', '--massif-out-file=@OUTPUT@', faest_em_256s_test],
    output: 'faest_em_256s.massif',
    depends: [faest_em_256s_test],
    install: false,
    build_always_stale: true,
    build_by_default: false,
  )
endif
if get_option('benchmarks').enabled()
  bench_sources = files(
    join_paths(meson.project_source_root(), 'tools', 'bench_c2.cpp'),
  )
  bench_catch = executable('faest_em_256s_bench_c2', bench_sources,
    dependencies: [libfaest_em_256s_dependency, boost_program_options, catch2],
    include_directories: include_directories,
    c_args: defines + c_flags,
    cpp_args: defines + cpp_flags
  )
endif
