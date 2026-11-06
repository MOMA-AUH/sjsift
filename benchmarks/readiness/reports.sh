#!/usr/bin/env bash
# Profile and independently audit the installed candidate, one stage at a time.
set -euo pipefail
root="${1:?Pass the absolute task directory}"
python="$(command -v -- "${PYTHON:-python}")"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd -- "$script_dir/../.." && pwd)"
[[ "$root" = /* ]] || { echo 'Task directory must be absolute' >&2; exit 2; }
cd "$repo"
samples=("${@:2}")
if [[ ${#samples[@]} -eq 0 ]]; then samples=(public controls); fi
for sample in "${samples[@]}"; do
  [[ "$sample" == public || "$sample" == controls ]] || { echo 'Unknown input set' >&2; exit 2; }
  inputs="$root/$sample-inputs"
  audit_seconds=120
  if [[ "$sample" == public ]]; then audit_seconds=300; fi
  "$python" benchmarks/readiness/measure.py --seconds 120 --memory-mib 2048 \
    "$root/$sample-reference-stage.json" -- "$python" benchmarks/readiness/profile_workflow.py \
    reference "$inputs" "$root/$sample-reference-only.json"
  for format in bam cram; do
    output="$root/$sample-$format-report"
    "$python" benchmarks/readiness/measure.py --seconds 120 --memory-mib 2048 \
      "$root/$sample-$format-stage.json" -- "$python" benchmarks/readiness/profile_workflow.py \
      report "$inputs" "$output" --format "$format"
    "$python" benchmarks/readiness/measure.py --seconds "$audit_seconds" --memory-mib 2048 \
      "$root/$sample-$format-audit-stage.json" -- "$python" benchmarks/readiness/audit_evidence.py \
      "$inputs" "$output" "$root/$sample-$format-inventory.json" --format "$format"
  done
  cmp "$root/$sample-bam-report/main.tsv" "$root/$sample-cram-report/main.tsv"
  cmp "$root/$sample-bam-report/context.tsv" "$root/$sample-cram-report/context.tsv"
  # Run the public module entry point without alignment to check additive outputs.
  "$python" -m sjsift --sample-name readiness --definitions "$inputs/catalog.toml" \
    --junctions "$inputs/sample.SJ.out.tab" --output "$root/$sample-main-only.tsv" \
    --context-output "$root/$sample-context-only.tsv" --provenance-output "$root/$sample-only.json"
  cmp "$root/$sample-main-only.tsv" "$root/$sample-bam-report/main.tsv"
  cmp "$root/$sample-context-only.tsv" "$root/$sample-bam-report/context.tsv"
  # A second installed-console run checks byte-deterministic provenance.
  "$python" benchmarks/readiness/measure.py --seconds 120 --memory-mib 2048 \
    "$root/$sample-repeat-stage.json" -- "$(dirname -- "$python")/sjsift" --definitions "$inputs/catalog.toml" \
    --junctions "$inputs/sample.SJ.out.tab" --alignment "$inputs/reads.bam" \
    --html-output "$root/$sample-repeated.html" --output "$root/$sample-repeated.tsv" \
    --provenance-output "$root/$sample-repeated.json"
  cmp "$root/$sample-bam-report/provenance.json" "$root/$sample-repeated.json"
done
