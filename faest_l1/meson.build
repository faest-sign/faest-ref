sources = configure_file(
    input: '../crypto_sign.c.in',
    output: 'crypto_sign.c',
    configuration: param_l1
)
headers = [
    configure_file(
        input: '../crypto_sign.h.in',
        output: 'crypto_sign.h',
        configuration: param_l1
    ),
    configure_file(
        input: '../crypto_sign.h.in',
        output: 'api.h',
        configuration: param_l1
    )
]

libfaest_l1 = static_library('faest_l1',
  sources,
  dependencies: libfaest_dependency,
  include_directories: include_directories,
  c_args: defines + c_flags
)
install_headers(headers, subdir: 'faest_l1')
libfaest_dependency = declare_dependency(
  link_with: libfaest,
  include_directories: include_directories
)
