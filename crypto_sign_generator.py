import sys
from pathlib import Path


def main():
    param_name = sys.argv[1]
    output_file = Path(f"faest_{param_name}") / "meson.build"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as outf:
        outf.write(
            f"""sources = configure_file(
    input: '../crypto_sign.c.in',
    output: 'crypto_sign.c',
    configuration: param_{param_name}
)
headers = [
    configure_file(
        input: '../crypto_sign.h.in',
        output: 'crypto_sign.h',
        configuration: param_{param_name}
    ),
    configure_file(
        input: '../crypto_sign.h.in',
        output: 'api.h',
        configuration: param_{param_name}
    )
]

libfaest_{param_name} = static_library('faest_{param_name}',
  sources,
  dependencies: libfaest_dependency,
  include_directories: include_directories,
  c_args: defines + c_flags
)
install_headers(headers, subdir: 'faest_{param_name}')
libfaest_dependency = declare_dependency(
  link_with: libfaest,
  include_directories: include_directories
)
"""
        )


if __name__ == "__main__":
    main()
