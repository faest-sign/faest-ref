project('faest-ref', ['c', 'cpp'],
  version: '0.0.0',
  meson_version: '>=0.57',
  default_options: [
    'c_std=c11',
    'cpp_std=c++17',
    'warning_level=3',
    'b_lto=true',
    'b_ndebug=if-release',
  ],
)

cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
conf_data = configuration_data()

# compile flags and defines
flags = []
c_flags = []
cpp_flags = []
defines = []

if cc.get_id() != 'msvc'
  c_flags += [
    '-Werror=implicit-function-declaration',
    '-Werror=vla',
  ]
  flags += '-Wshadow'
  defines += '-D_GNU_SOURCE'
else
  warnings = []
  # function inlining
  warnings += ['4710', '4711']
  # padding of structs
  warnings += ['4820', '4324']
  # .. in include paths
  warnings += '4464'
  # undefined macros evaluating to 0 in MSVC headers
  warnings += '4668'
  # Spectre mitigation
  warnings += '5045'
  # initialization of aggregate types with non-constant values is valid C99
  warnings += '4204'
  # nameless structs/unions are part of C11
  warnings += '4201'
  # conditional expression is constant
  warnings += '4127'
  # unary minus of unsigned types
  warnings += '4146'

  foreach w : warnings
    flags += '/wd' + w
  endforeach

  defines += '-D_CRT_SECURE_NO_WARNINGS'
endif
c_flags = cc.get_supported_arguments(flags) + cc.get_supported_arguments(c_flags)
cpp_flags = cpp.get_supported_arguments(flags) +  cc.get_supported_arguments(cpp_flags)

# check availability of some headers
conf_data.set('HAVE_SYS_RANDOM_H', cc.has_header('sys/random.h'))

# check availability of some functions
conf_data.set('HAVE_ALIGNED_ALLOC', cc.has_header_symbol('stdlib.h', 'aligned_alloc', args: defines))
conf_data.set('HAVE_POSIX_MEMALIGN', cc.has_header_symbol('stdlib.h', 'posix_memalign', args: defines))
conf_data.set('HAVE_MEMALIGN', cc.has_header_symbol('malloc.h', 'memalign', args: defines))
conf_data.set('HAVE_GETRANDOM', cc.has_header_symbol('sys/random.h', 'getrandom', args: defines))
conf_data.set('HAVE_ARC4RANDOM_BUF', cc.has_header_symbol('stdlib.h', 'arc4random_buf', args: defines))
conf_data.set('HAVE_EXPLICIT_BZERO', cc.has_header_symbol('string.h', 'explicit_bzero', args: defines))
conf_data.set('HAVE_CONSTTIME_MEMEQUAL', cc.has_header_symbol('string.h', 'consttime_memequal', args: defines))
conf_data.set('HAVE_TIMINGSAFE_BCMP', cc.has_header_symbol('string.h', 'timingsafe_bcmp', args: defines))

# parameter sets
# After changing the parameter sets, be sure to run crypto_sign_generator.py to create the corresponding meson files.
parameter_sets = []

# parameter set: FAEST-128s
param_128s = configuration_data()
param_128s.set('LAMBDA', 128)
param_128s.set('PARAM', '128S')
param_128s.set('PARAM_L', '128s')
param_128s.set('ELL', 1600)
param_128s.set('TAU', 11)
param_128s.set('K_0', 12)
param_128s.set('K_1', 11)
param_128s.set('T_0', 7)
param_128s.set('T_1', 4)
param_128s.set('PK_SIZE', 32)
param_128s.set('SK_SIZE', 16)
param_128s.set('SIG_SIZE', 5034)
parameter_sets += param_128s

# parameter set: FAEST-128f
param_128f = configuration_data()
param_128f.set('LAMBDA', 128)
param_128f.set('PARAM', '128F')
param_128f.set('PARAM_L', '128f')
param_128f.set('ELL', 1600)
param_128f.set('TAU', 16)
param_128f.set('K_0', 8)
param_128f.set('K_1', 8)
param_128f.set('T_0', 0)
param_128f.set('T_1', 16)
param_128f.set('PK_SIZE', 32)
param_128f.set('SK_SIZE', 16)
param_128f.set('SIG_SIZE', 6384)
parameter_sets += param_128f

# parameter set: FAEST-192s
param_192s = configuration_data()
param_192s.set('LAMBDA', 192)
param_192s.set('PARAM', '192S')
param_192s.set('PARAM_L', '192s')
param_192s.set('ELL', 3328)
param_192s.set('TAU', 16)
param_192s.set('K_0', 12)
param_192s.set('K_1', 12)
param_192s.set('T_0', 0)
param_192s.set('T_1', 16)
param_192s.set('PK_SIZE', 48)
param_192s.set('SK_SIZE', 24)
param_192s.set('SIG_SIZE', 12952)
parameter_sets += param_192s

# parameter set: FAEST-192f
param_192f = configuration_data()
param_192f.set('LAMBDA', 192)
param_192f.set('PARAM', '192F')
param_192f.set('PARAM_L', '192f')
param_192f.set('ELL', 3328)
param_192f.set('TAU', 24)
param_192f.set('K_0', 8)
param_192f.set('K_1', 8)
param_192f.set('T_0', 0)
param_192f.set('T_1', 24)
param_192f.set('PK_SIZE', 48)
param_192f.set('SK_SIZE', 24)
param_192f.set('SIG_SIZE', 17064)
parameter_sets += param_192f

# parameter set: FAEST-256s
param_256s = configuration_data()
param_256s.set('LAMBDA', 256)
param_256s.set('PARAM', '256S')
param_256s.set('PARAM_L', '256s')
param_256s.set('ELL', 4000)
param_256s.set('TAU', 22)
param_256s.set('K_0', 12)
param_256s.set('K_1', 11)
param_256s.set('T_0', 14)
param_256s.set('T_1', 8)
param_256s.set('PK_SIZE', 64)
param_256s.set('SK_SIZE', 32)
param_256s.set('SIG_SIZE', 22212)
parameter_sets += param_256s

# parameter set: FAEST-256f
param_256f = configuration_data()
param_256f.set('LAMBDA', 256)
param_256f.set('PARAM', '256F')
param_256f.set('PARAM_L', '256f')
param_256f.set('ELL', 4000)
param_256f.set('TAU', 32)
param_256f.set('K_0', 8)
param_256f.set('K_1', 8)
param_256f.set('T_0', 0)
param_256f.set('T_1', 32)
param_256f.set('PK_SIZE', 64)
param_256f.set('SK_SIZE', 32)
param_256f.set('SIG_SIZE', 28512)
parameter_sets += param_256f

# also make data from all parameter sets available "globally"
foreach parameter_set : parameter_sets
  parameter_set_prefix = 'FAEST_@0@'.format(parameter_set.get('PARAM'))
  foreach key : parameter_set.keys()
    conf_data.set('@0@_@1@'.format(parameter_set_prefix, key), parameter_set.get(key))
  endforeach
endforeach

# dependencies
build_dependencies = []

if host_machine.system() == 'windows'
  # require new enough Windows for bcrypt to be available
  defines += '-D_WIN32_WINNT=0x0601'
  build_dependencies += cc.find_library('bcrypt')
endif

# generate config.h
config_header = configure_file(
  output: 'config.h',
  output_format: 'c',
  configuration: conf_data
)
defines += '-DWITH_CONFIG_H'

include_directories = [include_directories('.')]

# source files
sources = files(
  'aes.c',
  'compat.c',
  'fields.c',
  'instances.c',
  'randomness.c',
  'random_oracle.c',
  'tree.c'
)

# header files to install
headers = files(
  'faest_defines.h'
)

# SHA3 implementation
sha3 = get_option('SHA3')
if sha3 == 'auto'
  if cc.sizeof('void*') == 4
    sha3 = 'plain32'
  else
    sha3 = 'opt64'
  endif
endif
include_directories += [include_directories('sha3'), include_directories('sha3/' + sha3)]

if sha3 == 's390-cpacf'
  defines += '-DWITH_SHAKE_S390_CPACF'
else
  defines += '-DWITH_KECCAK_X4'
  sources += files(
    'sha3/KeccakHash.c',
    'sha3/KeccakSponge.c',
    'sha3/KeccakSpongetimes4.c',
    'sha3/KeccakHashtimes4.c',
  )

  if sha3 == 'avx2'
    sources += files(
      'sha3/avx2/KeccakP-1600-AVX2.s',
      'sha3/avx2/KeccakP-1600-times4-SIMD256.c',
    )
  elif sha3 == 'opt64'
    sources += files(
      'sha3/opt64/KeccakP-1600-opt64.c',
      'sha3/opt64/KeccakP-1600-times4-on1.c',
    )
  elif sha3 == 'plain32'
    sources += files(
      'sha3/plain32/KeccakP-1600-inplace32BI.c',
      'sha3/plain32/KeccakP-1600-times4-on1.c',
    )
  elif sha3 == 'armv8a-neon'
    sources += files(
      'sha3/armv8a-neon/KeccakP-1600-armv8a-neon.s',
      'sha3/opt64/KeccakP-1600-times4-on1.c',
    )
  endif

  # 32 bit ARM is in general unhappy with unaligned access, hence disable it in this case
  if host_machine.cpu_family() == 'arm'
    defines += '-DDISABLE_MISALIGNED_ACCESS=1'
  endif
endif

# generate files for parameter sets
foreach parameter_set : parameter_sets
  headers += configure_file(
    input: 'faest_param.h.in',
    output: 'faest_@0@.h'.format(parameter_set.get('PARAM_L')),
    configuration: parameter_set
  )
  sources += configure_file(
    input: 'faest_param.c.in',
    output: 'faest_@0@.c'.format(parameter_set.get('PARAM_L')),
    configuration: parameter_set
  )
endforeach

# static faest library
libfaest_static = static_library('faest',
  sources,
  # dependencies: build_dependencies,
  include_directories: include_directories,
  c_args: defines + c_flags + ['-DFAEST_STATIC', '-DFAEST_EXPORT='],
  install: true
)
libfaest_static_dependency = declare_dependency(
  link_with: libfaest_static,
  include_directories: include_directories
)

# shared library
visibility_define = []
if host_machine.system() == 'windows'
  visibility_define += '-DFAEST_EXPORT=__declspec(dllexport)'
elif cc.has_function_attribute('visibility:default')
  visibility_define += '-DFAEST_EXPORT=__attribute__((visibility("default")))'
endif
libfaest = shared_library('faest',
  sources,
  # dependencies: build_dependencies,
  include_directories: include_directories,
  c_args: defines + c_flags + visibility_define,
  install: true,
  version: '0.0.0',
  gnu_symbol_visibility: 'hidden'
)
libfaest_dependency = declare_dependency(
  link_with: libfaest,
  include_directories: include_directories
)

# install headers
install_headers(headers)

# install pkg config file
pkg = import('pkgconfig')
pkg.generate(
  name: 'faest',
  description: 'FAEST signature scheme',
  libraries: [libfaest]
)

# SUPERCOP / NIST files
foreach parameter_set : parameter_sets
  subdir('faest_@0@'.format(parameter_set.get('PARAM_L')))
endforeach

subdir('tests')
