#!/bin/bash
#
# Create a compressed tarball of this replay working directory.
#
# The archive is written next to the replay directory by default, e.g.
#   <parent>/lad_replay.tar.gz
# An existing archive is never overwritten: a numbered suffix is used
# instead (lad_replay.tar.gz.1, .2, ...).
#
# Usage:
#   ./make_tar.sh [output_directory]
#
# The output directory may also be set with the TAR_DIR environment variable.

set -euo pipefail

# Directory this script lives in, with symlinks resolved. Everything else is
# derived from it, so the script works from any checkout, run from any cwd.
SRC_DIR="$(cd -- "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")" && pwd)"
SRC_NAME="$(basename -- "$SRC_DIR")"

# Where the archive goes: argument, else $TAR_DIR, else the parent directory.
OUT_DIR="${1:-${TAR_DIR:-$(dirname -- "$SRC_DIR")}}"
if [ ! -d "$OUT_DIR" ]; then
  echo "make_tar.sh: output directory does not exist: $OUT_DIR" >&2
  exit 1
fi
OUT_DIR="$(cd -- "$OUT_DIR" && pwd)"

TAR_PATH="$OUT_DIR/$SRC_NAME.tar.gz"
if [ -e "$TAR_PATH" ]; then
  i=1
  while [ -e "$OUT_DIR/$SRC_NAME.tar.gz.$i" ]; do
    i=$((i+1))
  done
  TAR_PATH="$OUT_DIR/$SRC_NAME.tar.gz.$i"
fi

# Note: no -h/--dereference, deliberately. ROOTfiles, REPORT_OUTPUT and raw
# are symlinks into multi-TB data areas; they are archived as links, never
# followed. Adding -h here would try to pack the entire raw data set.
echo "Creating $TAR_PATH ..."
tar --exclude='*.root' -czf "$TAR_PATH" -C "$SRC_DIR" .

echo "Done: $TAR_PATH ($(du -h -- "$TAR_PATH" | cut -f1))"
