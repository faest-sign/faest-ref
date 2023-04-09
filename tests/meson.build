hash_shake = executable('test_hash_shake', files('hash_shake.c'),
  dependencies: [libfaest_dependency],
  include_directories: include_directories,
  c_args: defines + flags
)
test('hash_shake', hash_shake)
