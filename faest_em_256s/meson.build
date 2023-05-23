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
  dependencies: libfaest_dependency,
  include_directories: include_directories,
  c_args: defines + c_flags
)
install_headers(headers, subdir: 'faest_em_256s')
libfaest_em_256s_dependency = declare_dependency(
  link_with: libfaest_em_256s,
  include_directories: include_directories
)
if openssl.found()
  tv_sources = files(join_paths(meson.source_root(), 'tools', 'rng.c'), join_paths(meson.source_root(), 'tools', 'PQCgenKAT_sign.cpp'))
  test_vector_generator = executable('feast_em_256s_test_vectors', [sources] + faest_sources + tv_sources,
    dependencies: [openssl],
    include_directories: include_directories,
    c_args: defines + c_flags + ['-DHAVE_RANDOMBYTES'],
    cpp_args: defines + cpp_flags + ['-DHAVE_RANDOMBYTES']
  )
endif
