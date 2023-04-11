boost_test = dependency('boost', required: false, modules: ['unit_test_framework'])
ntl = cc.find_library('ntl', required: false, has_headers: ['NTL/version.h'])

hash_shake = executable('test_hash_shake', files('hash_shake.c'),
  dependencies: [libfaest_dependency],
  include_directories: include_directories,
  c_args: defines + c_flags
)
test('hash_shake', hash_shake)

if boost_test.found() and ntl.found()
  fields = executable('test_fields', files('fields.cpp'),
    dependencies: [libfaest_dependency, boost_test, ntl],
    include_directories: include_directories,
    cpp_args: defines + cpp_flags
  )
  test('fields', fields)
endif
