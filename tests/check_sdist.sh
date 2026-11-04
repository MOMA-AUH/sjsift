#!/usr/bin/env bash
# Verify a built source archive without importing from the original checkout.
set -euo pipefail

archive="$(realpath "$1")"
work="$(mktemp -d "${RUNNER_TEMP:-${TMPDIR:-/tmp}}/sjsift-sdist.XXXXXX")"
trap 'rm -rf "$work"' EXIT
mkdir "$work/source"
tar -xzf "$archive" --strip-components=1 -C "$work/source"
cd "$work"

unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
python -m venv environment
export PATH="$work/environment/bin:$PATH"
python -m pip install build
python -m build --wheel --outdir "$work/wheels" "$work/source"
wheel=("$work"/wheels/*.whl)
python -m pip install "${wheel[0]}[dev]"

cd "$work/source"
python - <<'PY'
from pathlib import Path
import sysconfig
import sjsift

installed = Path(sysconfig.get_path("purelib")).resolve()
assert Path(sjsift.__file__).resolve().is_relative_to(installed), sjsift.__file__
PY
python -m pytest
sjsift --help
sjsift --version
touch "$work/empty.SJ.out.tab"
sjsift --sample-name smoke-test --junctions "$work/empty.SJ.out.tab" \
  --definitions definitions/grch38.toml \
  --output "$work/smoke.tsv" --html-output "$work/smoke.html"
python tests/installed_alignment_smoke.py
