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

boost_program_options = dependency('boost', required: get_option('benchmarks'), modules: ['program_options'])

libfaest_192s = static_library('faest_192s',
  sources,
  dependencies: libfaest_static_dependency,
  include_directories: include_directories,
  c_args: defines + c_flags
)
libfaest_192s_dependency = declare_dependency(
  link_with: libfaest_192s,
  include_directories: include_directories
)
if openssl.found()
  tv_sources = files(
    join_paths(meson.project_source_root(), 'randomness.c'),
    join_paths(meson.project_source_root(), 'tools', 'rng.c'),
    join_paths(meson.project_source_root(), 'tools', 'PQCgenKAT_sign.cpp')
  )
  test_vector_generator = executable('faest_192s_test_vectors', [sources] + tv_sources,
    dependencies: [libfaest_no_random_static_dependency, openssl],
    include_directories: include_directories,
    c_args: defines + c_flags + ['-DHAVE_RANDOMBYTES'],
    cpp_args: defines + cpp_flags + ['-DHAVE_RANDOMBYTES'],
    link_args: linker_flags,
    override_options: ['b_lto=false'],
  )
endif
if boost_program_options.found() and get_option('benchmarks').enabled()
  threads = dependency('threads')
  bench_sources = files(join_paths(meson.project_source_root(), 'tools', 'bench.cpp'))
  bench = executable('faest_192s_bench', bench_sources,
    dependencies: [libfaest_192s_dependency, boost_program_options, threads],
    include_directories: include_directories,
    c_args: defines + c_flags,
    cpp_args: defines + cpp_flags,
    link_args: linker_flags
  )
endif
test_sources = files(join_paths(meson.project_source_root(), 'tests', 'api_test.c'))
faest_192s_test = executable('faest_192s_api_test', test_sources,
  dependencies: [libfaest_192s_dependency, valgrind],
  include_directories: include_directories,
  c_args: defines + c_flags + valgrind_defines,
  link_args: linker_flags,
  override_options: ['b_lto=false'],
)
test('faest_192s_api_test', faest_192s_test,
  timeout: 6000,
)
if valgrind.found() and valgrind_exec.found()
  test('faest_192s_api_test_ct', valgrind_exec,
    args: ['-q', '--error-exitcode=1', '--track-origins=yes', faest_192s_test],
    timeout: 6000,
  )
endif
if valgrind_exec.found()
  custom_target('faest_192s_memory_usage',
    command: [valgrind_exec, '-q', '--error-exitcode=1', '--tool=massif', '--stacks=yes', '--massif-out-file=@OUTPUT@', faest_192s_test],
    output: 'faest_192s.massif',
    depends: [faest_192s_test],
    install: false,
    build_always_stale: true,
    build_by_default: false,
  )
endif
if get_option('catch2').enabled() and get_option('benchmarks').enabled()
  bench_sources = files(
    join_paths(meson.project_source_root(), 'tools', 'bench_c2.cpp'),
  )
  bench_catch = executable('faest_192s_bench_c2', bench_sources,
    dependencies: [libfaest_192s_dependency, boost_program_options, catch2],
    include_directories: include_directories,
    c_args: defines + c_flags,
    cpp_args: defines + cpp_flags,
    link_args: linker_flags
  )
endif
