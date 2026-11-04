#!/usr/bin/env bash
# Compare the published channel package with the separately measured candidate.
set -euo pipefail
root="${1:?Pass the absolute task directory}"
channel="${2:?Pass the absolute clean Conda environment}"
python="${PYTHON:-python}"
[[ "$root" = /* && "$channel" = /* ]] || exit 2
inputs="$root/public-inputs"
for format in bam cram; do
  output="$root/channel-$format"
  mkdir "$output"
  reference=()
  if [[ "$format" == cram ]]; then reference=(--reference "$inputs/reference.fa"); fi
  "$python" benchmarks/readiness/measure.py --seconds 120 --memory-mib 2048 \
    "$root/channel-$format-stage.json" -- "$channel/bin/sjsift" \
    --definitions "$inputs/catalog.toml" --junctions "$inputs/sample.SJ.out.tab" \
    --alignment "$inputs/reads.$format" "${reference[@]}" \
    --output "$output/main.tsv" --context-output "$output/context.tsv" \
    --html-output "$output/report.html"
  cmp "$output/main.tsv" "$root/public-$format-report/main.tsv"
  cmp "$output/context.tsv" "$root/public-$format-report/context.tsv"
done
"$python" - "$root" <<'PY'
from pathlib import Path
import json
import sys
sys.path.insert(0, 'benchmarks/readiness')
from audit_evidence import ReportInventory
root = Path(sys.argv[1])
result = {}
for fmt in ['bam', 'cram']:
    candidate = ReportInventory(root / f'public-{fmt}-report/report.html')
    channel = ReportInventory(root / f'channel-{fmt}/report.html')
    assert candidate.records == channel.records
    assert candidate.statuses == channel.statuses
    assert candidate.pairs == channel.pairs
    result[fmt] = {'tsv_comparison': 'passed', 'physical_inventory': 'passed',
                   'pair_inventory': 'passed', 'embedded_rows': sum(channel.records.values())}
(root / 'channel-comparison.json').write_text(json.dumps(result, indent=2) + '\n')
PY
