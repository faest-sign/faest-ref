boost_test = dependency('boost', required: false, modules: ['unit_test_framework'])
ntl = cc.find_library('ntl', required: false, has_headers: ['NTL/version.h'])

hash_shake = executable('test_hash_shake', files('hash_shake.c'),
  dependencies: [libfaest_static_dependency],
  include_directories: include_directories,
  c_args: defines + c_flags,
)
test('hash_shake', hash_shake)

tree = executable('tree', files('tree.c'),
  dependencies: [libfaest_static_dependency],
  include_directories: include_directories,
  c_args: defines + c_flags,
)
test('tree', tree)

vc = executable('vc', files('vc.c'),
  dependencies: [libfaest_static_dependency],
  include_directories: include_directories,
  c_args: defines + c_flags,
)
test('vc', vc)

vole = executable('vole', files('vole.c'),
  dependencies: [libfaest_static_dependency],
  include_directories: include_directories,
  c_args: defines + c_flags,
)
test('vole', vole)

faest = executable('faest', files('faest.c'),
  dependencies: [libfaest_static_dependency],
  include_directories: include_directories,
  c_args: defines + c_flags,
)
test('faest', faest)

if boost_test.found()
  test_sources = files('aes.cpp', 'extended_tests.cpp', 'instances.cpp', 'universal_hashing.cpp')
  test_deps = []
  if ntl.found()
    test_sources += files('fields.cpp')
    test_deps += ntl
  endif

  foreach parameter_set : parameter_sets
    test_sources += configure_file(
      input: 'basic.cpp.in',
      output: 'faest_@0@_basic.cpp'.format(parameter_set.get('PARAM_L')),
      configuration: parameter_set
    )
  endforeach

  extended_tests = executable('test_extended_tests', test_sources,
      dependencies: [libfaest_static_dependency, boost_test] + test_deps,
      include_directories: include_directories,
      cpp_args: defines + cpp_flags + ['-DBOOST_TEST_DYN_LINK'],
    )
  test('extended_tests', extended_tests)
endif