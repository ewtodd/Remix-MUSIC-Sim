#!/bin/sh
# Compatibility launcher; the converter itself is compiled by `make`.
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
converter="$script_dir/../legacy-msc-to-toml"
if [ ! -x "$converter" ]; then
  echo "legacy_msc_to_toml.py: $converter is missing; run 'make legacy-msc-to-toml' first" >&2
  exit 127
fi
exec "$converter" "$@"
