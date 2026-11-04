# Tested report capacity

The largest workload that completed every check in both scenarios and both
formats was **10,000 eligible records**, in the environment below.
This is an environment-specific tested envelope, not a universal record limit.
Every completed report passed the full retained-record inventory check.
Paging changes visibility and does not
reduce HTML size or initial browser memory.

This baseline was recorded on 2026-11-04 against sjsift 0.2.2 at production
commit `810d89e32a846d280f0cd619b8222ba7368c8d27` with the benchmark added in this change.
The source tree was dirty because the benchmark and this documentation were
being added. The committed harness accompanies the
[complete measured results](../benchmarks/capacity/baseline-2026-11-04.json).
The [benchmark instructions](../benchmarks/capacity/README.md) describe deterministic
inputs, independent retention checks, exact limits, and reproduction commands.

## Environment and limits

- Linux: `Linux-5.14.0-687.57.2.el9_8.x86_64-x86_64-with-glibc2.34`.
- CPU: AMD EPYC-Rome Processor; Python-reported CPU count 2; process affinity 16 CPUs.
- Python 3.14.7; pysam 0.24.1; HTSlib 1.24.
- Node/Playwright/Chromium versions are recorded per browser case in the raw JSON.
- Viewport: 1400 × 1000, headless Chromium, local `file:` opening, JavaScript enabled.
- Reference: 1,000,000 bases (1,012,506-byte FASTA), one synthetic contig.
- Limits declared before execution: 120 seconds and 2,048 MiB RSS per Python
  worker; 90 seconds and 4,096 MiB RSS per browser process tree; 30 seconds per
  opening and 10 seconds per interaction. Limits and every incomplete check are
  retained in the JSON. RSS is sampled, with possible brief overshoot.

The host is shared, with no reserved CPU or memory. Each operation has one sample;
there are no percentile or statistical confidence claims. Earlier exploratory runs
were superseded after fixture-base and process-memory accounting corrections.
The final baseline uses the corrected harness throughout. Available host memory
and dependency versions are retained in the raw environment record. The full
human-reference validation cost and representative biological samples are not
measured here.

## Measurements

Times below are seconds; sizes are MiB (1,048,576 bytes). Extraction includes
reference validation; the `ref` column isolates that nested cost. Render/extraction
times include `cProfile` overhead. Python RSS is the sampled process-tree peak.
Browser RSS includes Node and all Chromium descendants, with shared pages counted
more than once; it is not JavaScript heap size. Raw CDP heap and DOM counters are
also retained. A dash means the check did not complete, never zero cost.

| Records | Scenario | Format | Extract | Ref | Render | Python RSS | HTML | Browser RSS | Load | Slowest interaction | Browser result |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1,000 | single | BAM | 0.06 | 0.0000 | 0.21 | 44.36 | 1.88 | 662.30 | 0.40 | 0.19 | passed |
| 1,000 | single | CRAM | 0.07 | 0.0057 | 0.23 | 41.46 | 1.83 | 659.09 | 0.32 | 0.21 | passed |
| 1,000 | multiple | BAM | 0.12 | 0.0000 | 0.50 | 56.96 | 4.30 | 686.93 | 0.49 | 0.28 | passed |
| 1,000 | multiple | CRAM | 0.16 | 0.0061 | 0.53 | 53.09 | 4.21 | 678.72 | 0.41 | 0.19 | passed |
| 10,000 | single | BAM | 0.58 | 0.0000 | 2.02 | 172.09 | 18.54 | 961.68 | 2.03 | 1.61 | passed |
| 10,000 | single | CRAM | 0.65 | 0.0077 | 1.92 | 164.80 | 18.09 | 978.96 | 2.00 | 1.73 | passed |
| 10,000 | multiple | BAM | 1.48 | 0.0000 | 5.36 | 243.89 | 42.69 | 1221.43 | 3.50 | 2.41 | passed |
| 10,000 | multiple | CRAM | 1.53 | 0.0059 | 5.27 | 235.36 | 41.79 | 1303.22 | 3.87 | 2.19 | passed |
| 50,000 | single | BAM | 3.13 | 0.0000 | 10.91 | 692.83 | 92.68 | 2378.79 | 11.28 | 8.67 | passed |
| 50,000 | single | CRAM | 3.81 | 0.0064 | 9.67 | 675.44 | 90.44 | 2199.15 | 10.53 | 6.66 | passed |
| 50,000 | multiple | BAM | 6.57 | 0.0000 | 29.69 | 1091.55 | 213.48 | 3297.64 | 25.18 | 10.65 | failed |
| 50,000 | multiple | CRAM | 6.66 | 0.0041 | 29.44 | 1068.54 | 208.95 | 3284.21 | 25.28 | 12.53 | failed |
| 100,000 | single | BAM | 6.22 | 0.0000 | 20.48 | 1407.58 | 185.35 | 3889.10 | 24.66 | 14.69 | failed |
| 100,000 | single | CRAM | 5.88 | 0.0040 | 19.99 | 1345.77 | 180.87 | 3856.45 | 24.57 | 14.32 | failed |
| 100,000 | multiple | BAM | — | — | — | 2068.28 | — | — | — | — | incomplete |
| 100,000 | multiple | CRAM | — | — | — | 2058.10 | — | — | — | — | incomplete |

The single scenario has N embedded rows and N memberships. The multiple scenario
has 2N−2 embedded rows and 3N−2 memberships across two variant figures. This
includes exactly two rare defining records. Every completed static validation
checked the complete expected read-name inventory and unique per-figure record
IDs, defining STAR counts, and BAM/CRAM TSV equivalence. Every passing browser
case also preserved alignment identity/membership and STAR-support snapshots
through all viewing operations, including access to the rare defining pair.

## Observed constraints and reproduction

- 50,000 multiple BAM, browser: **failed**. AssertionError [ERR_ASSERTION]: mapping_all exceeded 10000 ms.
- 50,000 multiple CRAM, browser: **failed**. AssertionError [ERR_ASSERTION]: mapping_all exceeded 10000 ms.
- 100,000 single BAM, browser: **failed**. AssertionError [ERR_ASSERTION]: mapping_all exceeded 10000 ms.
- 100,000 single CRAM, browser: **failed**. AssertionError [ERR_ASSERTION]: mapping_all exceeded 10000 ms.
- 100,000 multiple BAM, report: **memory_limit**. No completed measurement.
- 100,000 multiple BAM, validation: **incomplete**. No completed measurement.
- 100,000 multiple BAM, browser: **incomplete**. No completed measurement.
- 100,000 multiple CRAM, report: **memory_limit**. No completed measurement.
- 100,000 multiple CRAM, validation: **incomplete**. No completed measurement.
- 100,000 multiple CRAM, browser: **incomplete**. No completed measurement.

Reproduce any row by running `python benchmarks/capacity/run.py /tmp/new-run
--sizes N` with the documented browser environment, replacing N with its record
count. This generates both scenarios and both formats with the same limits. The
JSON names each active browser check and preserves completed timings before an
interruption, so a failed operation can be distinguished from an unattempted one.
Inspect its per-stage logs and `.prof` output for local diagnosis. Those generated
HTML files, logs, and profiles are temporary outputs; the committed JSON is the
retained release-comparison baseline.

A timeout or memory-limit result demonstrates failure to complete this workload
under the declared budget. It does not establish data loss, a browser crash, or an
absolute maximum capacity. No downsampling or rendering redesign was introduced.
Workloads beyond 100,000 records, more figures/memberships, longer or more complex
reads, other browsers, and full human references remain unvalidated. Use a
representative local benchmark before treating a larger report as interactive.
