#!/usr/bin/env bash
set -euo pipefail
root="${1:?Pass the absolute task directory}"
python="${PYTHON:-python}"
star="${STAR:-STAR}"
samtools="${SAMTOOLS:-samtools}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd -- "$script_dir/../.." && pwd)"
[[ "$root" = /* ]] || { echo 'Task directory must be absolute' >&2; exit 2; }
[[ "$("$star" --version)" = 2.7.11b ]] || { echo 'STAR 2.7.11b required' >&2; exit 2; }
cd "$repo"
samples=("${@:2}")
if [[ ${#samples[@]} -eq 0 ]]; then samples=(public controls); fi
for sample in "${samples[@]}"; do
  [[ "$sample" == public || "$sample" == controls ]] || { echo 'Unknown input set' >&2; exit 2; }
  alignment_seconds=1800
  sort_seconds=120
  conversion_seconds=180
  if [[ "$sample" == public ]]; then
    alignment_seconds=3600
    sort_seconds=600
    conversion_seconds=900
  fi
  mkdir "$root/$sample-star"
  extra=()
  if [[ "$sample" == controls ]]; then
    extra=(--sjdbFileChrStartEnd "$root/controls-fastq/junctions.tab" --sjdbOverhang 99 --limitSjdbInsertNsj 100)
  fi
  "$python" benchmarks/readiness/measure.py --seconds "$alignment_seconds" --memory-mib 14336 \
    "$root/$sample-alignment.json" -- "$star" --runThreadN 2 \
    --genomeDir "$root/star-index" --genomeLoad NoSharedMemory \
    --readFilesIn "$root/$sample-fastq/R1.fastq" "$root/$sample-fastq/R2.fastq" \
    --outFileNamePrefix "$root/$sample-star/" --outSAMtype BAM Unsorted \
    --outSAMattrRGline ID:readiness SM:readiness \
    --outSAMattributes NH HI AS nM jM jI --outFilterMultimapNmax 20 \
    "${extra[@]}"
  mkdir "$root/$sample-inputs"
  "$python" benchmarks/readiness/measure.py --seconds "$sort_seconds" --memory-mib 2048 \
    "$root/$sample-sort.json" -- "$samtools" sort -@ 1 -m 256M \
    -o "$root/$sample-inputs/reads.bam" "$root/$sample-star/Aligned.out.bam"
  "$samtools" index "$root/$sample-inputs/reads.bam"
  cp definitions/grch38.toml "$root/$sample-inputs/catalog.toml"
  cp "$root/$sample-star/SJ.out.tab" "$root/$sample-inputs/sample.SJ.out.tab"
  ln -s "$root/data/reference.fa" "$root/$sample-inputs/reference.fa"
  ln -s "$root/data/reference.fa.fai" "$root/$sample-inputs/reference.fa.fai"
  "$python" benchmarks/readiness/measure.py --seconds "$conversion_seconds" --memory-mib 2048 \
    "$root/$sample-cram-conversion.json" -- "$samtools" view -@ 1 -C \
    -T "$root/data/reference.fa" -o "$root/$sample-inputs/reads.cram" \
    "$root/$sample-inputs/reads.bam"
  "$samtools" index "$root/$sample-inputs/reads.cram"
done
