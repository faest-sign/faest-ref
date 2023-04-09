project('faest-ref', 'c',
  version: '0.0.0',
  meson_version: '>=0.56',
  default_options: [
    'c_std=c11',
    'warning_level=3',
    'b_lto=true',
    'b_ndebug=if-release',
  ],
)

cc = meson.get_compiler('c')
conf_data = configuration_data()

# compile flags
flags = [
  '-Werror=implicit-function-declaration',
  '-Werror=vla',
]
flags = cc.get_supported_arguments(flags)

# defines
defines = ['-D_GNU_SOURCE']

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
  'compat.c',
  'randomness.c',
)

# header files to install
headers = files(
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

# faest library
libfaest = static_library('faest',
  sources,
  # dependencies: build_dependencies,
  include_directories: include_directories,
  c_args: defines + flags
)
libfaest_dependency = declare_dependency(
  link_with: libfaest,
  include_directories: include_directories
)

