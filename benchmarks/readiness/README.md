# Candidate workflow audit

These scripts reproduce the local workflow evidence in
[the readiness record](../../docs/RELEASE_READINESS.md).
They use an installed candidate package and the explicit reference catalog.
They do not change application code or define a supported Python library API.

The independent inspector reads STAR columns and BAM CIGAR operations directly.
It compares every catalog row and every embedded selected alignment membership.
Physical identities include coordinates, CIGAR, flags, and mapping metadata.
Displayed pairs must have reciprocal raw mate coordinates and compatible flags.
Both alignment sequence dictionaries must exactly match the complete FASTA index.
It never imports the sjsift matcher or evidence selector.
The output contains aggregate counts and representative alignment geometry without read names.
HTML reports still contain local read aliases. Keep reports outside the repository.

## Environment and resource limits

Use Linux, Python 3.11 or newer, pysam, STAR 2.7.11b, samtools, and the pinned
Playwright tools from `tests/package-lock.json`.
Install the candidate wheel into a clean environment.
Select that environment's Python for `profile_workflow.py` and `audit_evidence.py`.
Use a fresh task directory with at least 40 GiB free space.

Before preparing data, inspect the effective memory limit of every ancestor cgroup.
Host memory totals do not establish the account's allocation.
The recorded host exposed 431 GiB, but the account limit was 16 GiB.
Run one index, alignment, report, or browser phase at a time.
Reserve capacity for other processes within the same account.

`measure.py` reuses the [capacity watchdog](../capacity/README.md).
It writes limits before execution and retains failures.
It requires real `/proc` RSS counters and kills the command's process group at a limit.
RSS is resident memory. Process-tree RSS can count shared pages more than once.
The watchdog samples memory and cannot prevent every short overshoot.
The outer command exits nonzero when the measured command fails.

Place options before the result path and the command after `--`:

```bash
python benchmarks/readiness/measure.py --seconds 120 --memory-mib 2048 \
  /tmp/audit/report-stage.json -- \
  /tmp/candidate/bin/python benchmarks/readiness/profile_workflow.py report \
  /tmp/audit/inputs /tmp/audit/bam-report --format bam
```

## Data preparation

The exact acquisition URLs, hashes, STAR commands, and failed attempts are in the
readiness record. Download the complete files before creating the subset.
`prepare.py subset` verifies both published ENA MD5 values and paired names.
It replaces source read names with ordinal aliases.
Use `--pairs 22935521` for the complete recorded run. This also checks both file ends.
The default one-million-pair prefix is a quick pilot, not a random sample.
Neither mode establishes disease prevalence or sensitivity to rare variants.

```bash
python benchmarks/readiness/prepare.py subset /tmp/audit/data \
  /tmp/audit/public-fastq --pairs 22935521
python benchmarks/readiness/prepare.py controls /tmp/audit/data \
  /tmp/audit/controls-fastq --catalog definitions/grch38.toml
```

The `data` directory contains `R1.fastq.gz`, `R2.fastq.gz`, and an indexed
`reference.fa`. The reference contains all 455 UCSC GRCh38 contigs.
The controls create four artificial overlapping mate pairs per distinct catalog
junction from the full reference sequence. They also emit a STAR junction annotation.
Keep their alignment and reporting separate from the public observations.
These controls establish software behavior only.

After the full index finishes, run `align.sh` with an absolute task directory.
Set `PYTHON`, `STAR`, and `SAMTOOLS` to explicit executable paths.
The script expects the index at `star-index`, and the two prepared read directories.
It aligns the public data and controls sequentially, then sorts and converts each result.

## Reports and independent comparisons

Set `PYTHON` to the absolute candidate-environment interpreter and run:

```bash
PYTHON=/tmp/candidate/bin/python bash benchmarks/readiness/reports.sh /tmp/audit
```

This runs the bounded reference, report, and independent inventory stages below.
It also compares both TSV formats with TSV-only output and repeats installed-console
reporting to check byte-deterministic provenance.

For each input set, prepare `catalog.toml`, `sample.SJ.out.tab`, indexed `reads.bam`,
indexed `reads.cram`, and indexed `reference.fa`.
Preserve the full sequence dictionary when converting BAM to CRAM.
Give the alignment one neutral read-group sample label.
Use symlinks for large local inputs when appropriate.

`profile_workflow.py report` calls the installed CLI implementation under `cProfile`.
The process records whole-report peak RSS and separates reference validation,
extraction, and HTML rendering time. Profiling adds overhead.
`profile_workflow.py reference` runs only full-reference validation in a fresh process.
Its RSS includes Python, imports, and header loading. It is not a subtraction of peaks.
The function it instruments is private and can change after the pinned candidate.

```bash
python benchmarks/readiness/measure.py --seconds 120 --memory-mib 2048 \
  /tmp/audit/reference-stage.json -- \
  /tmp/candidate/bin/python benchmarks/readiness/profile_workflow.py reference \
  /tmp/audit/inputs /tmp/audit/reference-only.json
python benchmarks/readiness/audit_evidence.py /tmp/audit/inputs \
  /tmp/audit/bam-report /tmp/audit/bam-inventory.json
```

Repeat report profiling and inspection with `--format cram` and new destinations.
Pass the format to both scripts; the inspector checks it against provenance.
If the report used `--sample-name`, pass the same option to the inspector and
sensitivity check. Otherwise they independently infer the one header sample.
Readiness inputs must contain at most one distinct header sample; this inventory
does not implement selection among multiple samples.
The inspector scans the complete original BAM for an independent inventory for both formats.
Its full-public-run scan has a 300-second budget. This audit budget is separate from
the unchanged 120-second report-generation budget.
Compare both TSV files byte for byte between formats and against a TSV-only run.
The inspector also checks provenance schema, input hashes, membership counts,
sample-selection facts, assembly-label comparisons, and CRAM reference checks.
Expected metadata comes from the invocation, catalog, alignment header, and local
FASTA. CRAM M5 checksums are independently compared with complete local sequences;
the expected checked count is the number of available header M5 values.
BAM must report no reference object. These facts do not establish the anonymous
STAR input's sample identity or assembly. The inspector does not compare alignment counts with STAR support.
Those counts have different meanings.

## Offline browser validation

Use `offline.sh` to copy all four reports, move the complete task directory away,
and restore it on exit. The browser sees the original source paths as unavailable.
Select a watchdog interpreter outside the task directory because that directory moves.

```bash
PYTHON=/path/outside/task/bin/python bash benchmarks/readiness/offline.sh \
  /tmp/audit /tmp/offline
```
Run `browser.cjs` under `measure.py` with 90 seconds and 4096 MiB.
Set `NODE_PATH` and `PLAYWRIGHT_BROWSERS_PATH` for the pinned browser tools.

```bash
python benchmarks/readiness/measure.py --seconds 90 --memory-mib 4096 \
  /tmp/offline/browser-stage.json -- node benchmarks/readiness/browser.cjs \
  /tmp/offline/report.html /tmp/offline/browser.json
```

The browser uses `file:` URLs, offline mode, and blocked network routes.
It requires all 17 views, selected reads, and a safely associated mate pair.
It checks mapping and support filters, paging when present, pairing, and read inspection.
It compares the complete retained evidence before and after interaction.
It also opens the report with JavaScript disabled and inspects a native read preview.
Opening must finish within 30 seconds. Each interaction must finish within 10 seconds.
An absent required observation fails the audit instead of becoming a silent skip.

## Audit sensitivity

After preparing the controlled report, run this small mutation check:

```bash
python benchmarks/readiness/check_audit.py /tmp/audit/controls-inputs \
  /tmp/audit/controls-bam-report /tmp/audit/audit-sensitivity.json
```

It first requires the unmodified report to pass. It then alters a STAR count,
removes a physical record, changes a pair association, substitutes a same-name mate,
changes the provenance schema, and corrupts the header sample count in structured
sample-selection metadata.
Each altered report must fail an audit assertion. The restored report must pass.
Mutations use temporary copies and preserve the original reports.
The mate substitution keeps the old name/class/membership key unchanged.
The audit must distinguish the mates through coordinates, CIGAR, flags, and mapping metadata.

The small installed-CLI regression uses the committed STAR output and temporary
CRAM conversions. It requires neither STAR nor downloaded data:

```bash
python -m pytest -q tests/test_readiness_audit.py tests/test_star_interop.py
```

On 2026-11-06, the new regression first reproduced `KeyError: 'star_assurance'`
for BAM and CRAM. After the schema-1 repair, all six targeted cases passed,
including explicit sample names, inferred header samples, missing/partial assembly
labels, six rejected mutations, and unmodified/restored reports. This is tooling
regression evidence; the historical public-data readiness record is unchanged.

## Published channel comparison

Create a fresh environment with the recommended channel command in the readiness
record, then run `channel.sh /tmp/audit /tmp/audit/channel-env` from the repository.
Set `PYTHON` to the candidate interpreter. This compares both TSVs, physical HTML
inventories, strand statuses, and pair associations with the candidate reports.
It deliberately omits provenance: the pinned published 0.2.2 package predates that
option. The candidate wheel and candidate Conda CI are separate evidence.
