#!/usr/bin/env bash
# Move source paths away, then exercise copied reports with network access denied.
set -euo pipefail
root="${1:?Pass the absolute task directory}"
offline="${2:?Pass a new absolute directory for standalone reports}"
python="${PYTHON:-python}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd -- "$script_dir/../.." && pwd)"
[[ "$root" = /* && "$offline" = /* ]] || { echo 'Absolute directories required' >&2; exit 2; }
[[ "$offline/" != "$root/"* ]] || { echo 'Offline directory must be outside the input directory' >&2; exit 2; }
[[ ! -e "$root.unavailable" ]] || { echo 'Temporary input destination already exists' >&2; exit 2; }
mkdir "$offline"
for sample in public controls; do
  for format in bam cram; do
    cp "$root/$sample-$format-report/report.html" "$offline/$sample-$format.html"
  done
done
mv "$root" "$root.unavailable"
trap 'mv "$root.unavailable" "$root"' EXIT
[[ ! -e "$root" ]]
cd "$repo"
for sample in public controls; do
  for format in bam cram; do
    "$python" benchmarks/readiness/measure.py --seconds 90 --memory-mib 4096 \
      "$offline/$sample-$format-stage.json" -- node benchmarks/readiness/browser.cjs \
      "$offline/$sample-$format.html" "$offline/$sample-$format.json"
  done
done
