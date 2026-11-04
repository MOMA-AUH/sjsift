# Report capacity benchmark

This opt-in Linux benchmark measures retained alignment evidence through the real
`sjsift` CLI and Chromium's local `file:` report interface. It changes no production
code and performs no downsampling. It is separate from `pytest`'s fast `tests/`
suite and from required CI. Use the same command and settings for release comparisons.

## Run

Install the package and development dependencies as described in
[DEVELOPMENT.md](../../docs/DEVELOPMENT.md), then install the pinned browser tools:

```bash
npm ci --prefix tests
tests/node_modules/.bin/playwright install --with-deps chromium
NODE_PATH="$PWD/tests/node_modules" python benchmarks/capacity/run.py /tmp/capacity-run
```

Use a new output directory for every run. An existing separate Playwright install
can be selected with `NODE_PATH` and `PLAYWRIGHT_BROWSERS_PATH`, as in the required
browser checks. The runner requires real Linux `/proc` counters and permission to
launch Chromium. A sandbox that masks process memory is rejected before execution.
The benchmark performs no network requests after dependencies are installed.

Run its separate contract checks with:

```bash
python -m pytest -q benchmarks/capacity/test_capacity.py
```

`--sizes 1000` gives a shorter complete measurement. `--skip-browser` measures only
input/report generation and retained HTML validation; skipped browser cases never
establish a completed capacity. To vary reference-validation cost, repeat with
`--reference-bases 100000000` (100 million bases), keeping the other settings fixed.
These synthetic reference sizes are not a measurement of a full human assembly.

## Limits declared before execution

The runner writes settings and environment to `results.json` before any workload:

| Limit | Default |
| --- | ---: |
| Each fixture, report, or HTML validation worker | 120 seconds |
| Each browser process tree | 90 seconds |
| Local report opening | 30 seconds |
| Individual interaction | 10 seconds |
| Each Python worker process tree RSS | 2,048 MiB |
| Each Node/Chromium process tree RSS | 4,096 MiB |
| Reference sequence | 1,000,000 bases |
| Records per workload | 1,000; 10,000; 50,000; 100,000 |

All limits have corresponding command options (`--help`). RSS is sampled every
10 ms, summed across descendants, and includes shared pages more than once. This
is a watchdog, not a hard kernel memory allocation limit: a short overshoot between
samples is possible. It terminates the entire process group at a limit. Each stage
records elapsed time, sampled peak RSS, exit code, and a log. Python also reports
its kernel high-water RSS. A host OOM kill or interpreter/browser crash remains a
failure; it is never silently excluded. A missing downstream measurement is
`incomplete`, and an explicit omission is `skipped`. Browser progress is saved
between checks, including the active check when a timeout interrupts it. Completed
and failed workloads remain together in the final JSON. Any failed/incomplete
requested check makes a full run exit nonzero. The no-browser mode returns zero
only when all generation and retention checks pass.

Budget disk space for both formats' HTML and profiling output. The largest mixed
fixture embeds almost 200,000 rows in each HTML file. Outputs are deliberately not
committed; keep the small results JSON for comparison and regenerate the inputs.

## Workloads and independent evidence

Every size runs two scenarios, each through indexed BAM and explicit-reference
indexed CRAM. Both use one synthetic contig, deterministic sequence and names,
coordinate sorting, valid reciprocal read1/read2 flags and mate coordinates, and a
single sample read group. Inputs are semantically deterministic; CRAM headers
include local reference paths and encoding bytes can vary with HTSlib versions.
No human data or clinical validation is implied.

- **single:** all N records support one abundant defining splice junction.
- **multiple:** N−2 records support the abundant defining junction and the same
  junction as a reference for a second variant. The first N/2 records also support
  a second shared reference junction. The last two records form the rare defining
  pair for the second variant. There are N physical eligible records, 2N−2 embedded
  rows across two figures, and 3N−2 memberships across their groups.

Pair index modulo five cycles through `20M100N20M`, `5S20M100N20M5S`,
`10M1I10M100N20M`, `10M1D9M100N20M`, and `10=10X100N20M`. The shared subset appends
`100N20M`. This covers clipping, insertion, deletion, explicit match/mismatch,
and multiple skipped introns. Paired names share the same CIGAR. Mapping tags
cycle through three unique pairs, one multimapping pair, and one unknown pair.
All records carry agreeing `XS:+` evidence. STAR inputs deliberately report fixed
counts independent of retained alignment counts (abundant 800/200, shared 300/100,
rare 1/1); changing view controls must preserve them.

An independent streaming HTML inventory checks every expected read-name
multiplicity per figure, unique per-figure record IDs, exact row/membership counts,
and exact defining STAR support in TSV. It also compares BAM/CRAM defining and
reference TSVs. This detects omissions, duplicates, and loss of rare support.
Chromium verifies counts again, snapshots STAR support and a compact
order-independent alignment identity/membership fingerprint, exercises controls,
and compares the snapshot afterwards. The fingerprint is a regression detector,
not a cryptographic proof; the complete static name inventory is checked separately.

Browser checks open the local overview, select a variant, expand reads, paginate
forward/back, apply/reset mapping and support filters, pair/unpair reads, and
select the rare variant and its defining support. Assertions check the resulting
visible records, not merely dispatch of an event. Latency includes the action,
its observable result, and two animation frames. Each operation is measured once;
these are baseline observations, not percentile estimates or a service guarantee.
Pagination hides rows that are already embedded and loaded.

## Metrics and interpretation

`results.json` includes OS/kernel, CPU model/count/affinity, host memory at start,
Python/package/pysam/HTSlib versions, Git revision/dirty state, settings, and every
case. Browser records add Node, Playwright, Chromium, viewport-fixed timings,
element count, CDP DOM counters, and JavaScript heap sizes before/after interactions.
Process-tree RSS includes native Chromium memory; JS heap alone does not. RSS is
not a private native-heap metric. Host totals are not a reservation of resources.

Fixture creation records BAM/input generation and CRAM conversion separately,
including reference base/file size and BAM/CRAM bytes. Report workers profile the
actual CLI with standard-library `cProfile`. `extraction_seconds` includes input
validation and regional extraction; `reference_validation_seconds` is its nested
full-reference validation cost. The difference includes selected-region decoding,
classification and serialization plus other extraction overhead. It is not a
second extraction. `render_seconds` measures HTML rendering; outer wall time
includes startup, TSV writing, and profiling overhead. These profiled timings are
suitable for comparisons with the same harness, not uninstrumented CLI estimates.
The `.prof` files retain the detailed breakdown. No supported library API is added.

`largest_complete_records` requires both scenarios and both formats to finish
all checks at that size. It says nothing about untested sizes, longer reads,
additional variants/references, full human reference cost, different browsers,
or smaller machines. Read [the recorded baseline](../../docs/CAPACITY.md) before
choosing workloads for an interactive session.
