boost_test = dependency('boost', required: false, modules: ['unit_test_framework'])
ntl = cc.find_library('ntl', required: false, has_headers: ['NTL/version.h'])

if boost_test.found()
  test_sources = files(
    'aes.cpp',
    'aes_prove.cpp',
    'extended_tests.cpp',
    'hash_shake.cpp',
    'instances.cpp',
    'test_faest_tvs.cpp',
    'universal_hashing.cpp',
    'vc.cpp',
    'vole.cpp',
    'vole_tvs.cpp',
  )
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
    cpp_args: defines + cpp_flags + ['-DBOOST_TEST_DYN_LINK', '-DFAEST_TESTS'],
    override_options: ['b_lto=false']
  )
  test('AES', extended_tests,
    args: ['-t', 'aes'],
    timeout: 6000,
  )
  test('AES prove', extended_tests,
    args: ['-t', 'test_aes_prove'],
    timeout: 6000,
  )
  test('SHAKE', extended_tests,
    args: ['-t', 'hash_shake'],
    timeout: 6000,
  )
  test('Instances', extended_tests,
    args: ['-t', 'instances'],
    timeout: 6000,
  )
  test('Universal hashing', extended_tests,
    args: ['-t', 'universal_hashing'],
    timeout: 6000,
  )
  test('Vector commitments', extended_tests,
    args: ['-t', 'vector_commitments'],
    timeout: 6000,
  )
  test('VOLE', extended_tests,
    args: ['-t', 'vole'],
    timeout: 6000,
  )
  if ntl.found()
    test('Fields', extended_tests,
      args: ['-t', 'fields'],
      timeout: 6000,
    )
  endif
  foreach parameter_set : parameter_sets
    test('FAEST_@0@'.format(parameter_set.get('PARAM')), extended_tests,
      args: ['-t', 'faest_@0@'.format(parameter_set.get('PARAM_L'))],
      timeout: 6000,
    )
  endforeach
endif
