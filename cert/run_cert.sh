#!/usr/bin/env bash
# Usage (from a fresh clone, in targets/hart-rowley-ancestor/cert/):  bash run_cert.sh <script.py> <outfile> [args...]
# Writes <outfile> with command, date, commit, versions, script sha256, full output and exit status.
set -u
script="$1"; out="$2"; shift 2
{
  echo "command: python $script $*"
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "commit: $(git rev-parse HEAD) (fresh git clone of the local repository at this HEAD; the uncommitted cert/ scripts were copied in, see sha256)"
  echo "uncommitted cert files sha256:"; sha256sum hr_tables.py hr_tables_x.py hr_lean_gen.py run_cert.sh | sed 's/^/  /'
  echo "python: $(python -c 'import sys; print(sys.version.split()[0])')  platform: $(python -c 'import platform; print(platform.platform())')"
  echo "libraries: Python standard library only (itertools, collections); no numpy/sympy/flint used"
  echo "----"
  start=$(date +%s)
  python "$script" "$@"
  status=$?
  end=$(date +%s)
  echo "----"
  echo "exit status: $status ; runtime: $((end-start)) s"
} > "$out" 2>&1
tail -3 "$out"
