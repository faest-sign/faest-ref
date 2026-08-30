presets = ['128s', '128f', '192s', '192f', '256s', '256f', 'em_192s', 'em_192f']
tables = []

boost = dependency('boost', required: false)
if boost.found()
  generator = executable(
    'vole_mult_tables',
    files('vole_mult_tables.cpp') + [parameters_header],
    native: true,
    cpp_args: cpp_flags,
    dependencies: boost,
    override_options: ['b_lto=false'],
  )

  foreach preset : presets
    tables += custom_target(
      'generate_tables_@0@'.format(preset),
      command: [
        generator,
        '--preset', preset,
        '--emit-c', join_paths('@OUTDIR@', 'tables_@0@.h'.format(preset)),
      ],
      output: [
        'tables_@0@.h'.format(preset),
      ],
      console: false,
    )
  endforeach
else
  fs = import('fs')

  foreach preset : presets
    header = 'tables_@0@.h'.format(preset)
    tables += fs.copyfile(
      join_paths(meson.current_source_dir(), 'pregenerated', header),
      header,
    )
  endforeach
endif