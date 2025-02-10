executable(
  'generate_field_constants',
  files('generate_field_constants.cpp'),
  dependencies: [libfaest_static_dependency],
  include_directories: include_directories,
  cpp_args: defines + test_defines + cpp_flags,
  link_args: linker_flags,
  override_options: ['b_lto=false'],
)