boost_test = dependency('boost', required: false, modules: ['unit_test_framework'])
ntl = cc.find_library('ntl', required: false, has_headers: ['NTL/version.h'])

hash_shake = executable('test_hash_shake', files('hash_shake.c'),
  dependencies: [libfaest_static_dependency],
  include_directories: include_directories,
  c_args: defines + c_flags,
  override_options: ['b_lto=false']
)
test('hash_shake', hash_shake)

tree = executable('tree', files('tree.c'),
  dependencies: [libfaest_static_dependency],
  include_directories: include_directories,
  c_args: defines + c_flags,
  override_options: ['b_lto=false']
)
test('tree', tree)

if boost_test.found()
  if ntl.found()
  fields = executable('test_fields', files('fields.cpp'),
    dependencies: [libfaest_static_dependency, boost_test, ntl],
    include_directories: include_directories,
    cpp_args: defines + cpp_flags,
    override_options: ['b_lto=false']
  )
  test('fields', fields)
  endif

  aes_ctr = executable('test_aes_ctr', files('aes.cpp'),
    dependencies: [libfaest_static_dependency, boost_test],
    include_directories: include_directories,
    cpp_args: defines + cpp_flags,
    override_options: ['b_lto=false']
  )
  test('aes_ctr', aes_ctr)

  instances = executable('test_instances', files('instances.cpp'),
    dependencies: [libfaest_static_dependency, boost_test],
    include_directories: include_directories,
    cpp_args: defines + cpp_flags,
    override_options: ['b_lto=false']
  )
  test('instances', instances)
endif