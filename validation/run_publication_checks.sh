#!/usr/bin/env bash
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python="${PYTHON:-python}"
out="${1:?Usage: bash validation/run_publication_checks.sh OUTPUT_DIRECTORY}"
mkdir -p "$out"
out="$(cd "$out" && pwd)"
cd "$repo"
git rev-parse HEAD > "$out/commit.txt"
"$python" -m pip freeze > "$out/installed_versions.txt"
export PYTHONPATH="$repo/src:$repo"
export MPLBACKEND=Agg
failed=0
for suite in tests validation/real_data; do
  label="${suite//\//_}"
  if "$python" -m unittest discover -s "$suite" -p 'test_*.py' -v > "$out/$label.log" 2>&1; then
    echo "$label PASS"
  else
    echo "$label FAIL: $out/$label.log"
    failed=1
  fi
done
exit "$failed"
