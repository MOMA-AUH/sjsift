# Candidate release-readiness record

This record validates candidate commit
[`484ce57891da9446bfb5a177d6d375b9be0c2e99`](https://github.com/MOMA-AUH/sjsift/commit/484ce57891da9446bfb5a177d6d375b9be0c2e99)
on 2026-11-04. It addresses [issue #56](https://github.com/MOMA-AUH/sjsift/issues/56).
The application still identifies itself as 0.2.2. This work changes no version,
publishes no package, and creates no tag or release.

The [compatibility contract](COMPATIBILITY.md) takes effect with v1.0.0.
The [capacity contract](CAPACITY.md) describes the previously measured synthetic envelope.
The audit applies both to the candidate without extending their claims.

## Candidate and artifact identity

The runtime Git tree is `dfdcf4f0a15e4924e2525241583b211587084c64`.
The complete schema-3 GRCh38 catalog contains 17 variants and 32 reference rows.
These 49 rows represent 39 distinct junctions. The catalog is an explicit input.

| Artifact | SHA-256 |
| --- | --- |
| `definitions/grch38.toml` | `a532f7bf56618991d403d7b3c87fc2d398d4203695b127caab68e8874c557ba3` |
| Candidate wheel, 41,890 bytes | `c7dba00f85cef06cad5c630a2c50e452149297db9cbc71cc6467143343f6525f` |
| Candidate source archive, 122,524 bytes | `72e78c039aa789f5c413b2fc907597cbbca7a4ca188a1e9dc731e18979380c21` |
| Published Conda `sjsift-0.2.2-py_0.conda` | `3fa5cd265e51b63502cd6dfa1423f67ba08e0b883c83b6114c9e13df3ac56192` |

The wheel and source archive were built from the candidate before this audit added
scripts and records. Their runtime and catalog inputs remain unchanged.
Archive bytes can differ across rebuilds because build metadata includes timestamps.
Retain the measured artifacts when comparing their exact hashes.

The published Conda package and the candidate both report 0.2.2, but their contents differ.
The published package does not provide `--provenance-output`.
Its clean installation establishes channel availability and its own behavior.
Candidate-wheel tests and candidate Conda CI establish behavior of the proposed code.
Do not use the shared version string to treat those artifacts as identical.

## Runtime and public data

The contract supports Python 3.11 or newer. Linux CI tests Python 3.11 through 3.14,
including minimum pysam 0.24.0 on Python 3.11.
The local candidate uses Python 3.14.7 and pysam 0.24.1.
The clean published Conda environment resolved Python 3.13.18 and pysam 0.24.1.
No additional platform or browser matrix is established by this audit.

The public RNA-seq source is [ENA run SRR1039508](https://www.ebi.ac.uk/ena/browser/view/SRR1039508),
experiment [SRX384345](https://www.ebi.ac.uk/ena/browser/view/SRX384345), study PRJNA229998.
ENA describes paired transcriptomic cDNA sequencing on Illumina HiSeq 2000.
The source contains 22,935,521 paired spots and 2,889,875,646 bases.
The initial pilot used the first 1,000,000 pairs. It supplied only 15 selected
alignment records and no complete selected mate group.
The audit therefore expanded to the complete public run. One sample is not a cohort.
Source read names are replaced by ordinal aliases before local alignment.
No private or clinical dataset was accessed.

| Full source file | Bytes | Published and observed MD5 |
| --- | ---: | --- |
| `SRR1039508_1.fastq.gz` | 1,243,786,980 | `d57df747bc142e9850074d512ab9d6db` |
| `SRR1039508_2.fastq.gz` | 1,230,997,455 | `3331c6a9e0183ff9d398a3292dd45f66` |
| UCSC `hg38.fa.gz` | 983,659,424 | `1c9dcaddfa41027f17cd8f7a82c7293b` |

The reference is the [complete UCSC initial GRCh38 assembly](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/),
accession GCA_000001405.15. It contains 455 sequences and 3,209,286,105 bases.
It includes primary chromosomes, unplaced sequences, and alternate loci.
It is not the GENCODE GRCh38.p14 FASTA. UCSC documents that patch additions do not
change the existing chromosome sequences used by this catalog.
The same complete indexed FASTA is used for alignment, CRAM encoding, and validation.

## Reproduction

Use the [audit scripts](../benchmarks/readiness/README.md) and pinned candidate wheel.
The commands below use a task-owned directory named `$audit` and explicit tool paths
named `$python`, `$star`, and `$samtools`. Set these for the local environment.
Never place raw genomic inputs or generated HTML in the repository.
The retained aggregate evidence contains no raw read names or sample labels.

```bash
mkdir -p "$audit/data"
curl -fL 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR1039508&result=read_run&fields=run_accession,library_layout,fastq_ftp,fastq_bytes,fastq_md5&format=tsv' \
  -o "$audit/data/ena.tsv"
curl -fL http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz \
  -o "$audit/data/R1.fastq.gz"
curl -fL http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_2.fastq.gz \
  -o "$audit/data/R2.fastq.gz"
curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
  -o "$audit/data/hg38.fa.gz"
curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt \
  -o "$audit/data/ucsc-md5sum.txt"
# Compare the three full-file MD5 hashes with the values above before proceeding.
gzip -dc "$audit/data/hg38.fa.gz" > "$audit/data/reference.fa"
"$samtools" faidx "$audit/data/reference.fa"
"$python" benchmarks/readiness/prepare.py subset "$audit/data" "$audit/public-fastq" --pairs 22935521
"$python" benchmarks/readiness/prepare.py controls "$audit/data" "$audit/controls-fastq"
mkdir "$audit/star-index"
"$python" benchmarks/readiness/measure.py --seconds 21600 --memory-mib 14336 \
  "$audit/index.json" -- "$star" --runMode genomeGenerate --runThreadN 2 \
  --genomeDir "$audit/star-index" --genomeFastaFiles "$audit/data/reference.fa" \
  --genomeSAsparseD 4 --genomeSAindexNbases 13 --limitGenomeGenerateRAM 12000000000 \
  --outFileNamePrefix "$audit/star-index/"
```

The sparse suffix array and shorter search index reduce index memory at a cost in
alignment speed. They are supported STAR settings, not a truncated reference.
The audit pins these settings and does not claim equivalence to every STAR configuration.
Run the sequential alignment and format-conversion script after the index finishes:

```bash
PYTHON="$python" STAR="$star" SAMTOOLS="$samtools" \
  bash benchmarks/readiness/align.sh "$audit"
```

It uses two STAR threads, unshared genome memory, and unsorted BAM output.
Samtools sorts with one worker and a 256 MiB buffer, then indexes both formats.
The public alignment uses no junction annotation. The separate controls supply
only the 39 catalog junctions and a 99-base annotation overhang.
Both alignments emit `NH HI AS nM jM jI` and allow up to 20 mapping loci.
The full public alignment permits 60 minutes and 14 GiB.
Its sorting permits 600 seconds and 2 GiB. Its CRAM conversion permits 900 seconds and 2 GiB.
The small control setup retains 30-minute alignment, 120-second sorting, and
180-second conversion limits, with the same respective memory limits.
These are input-preparation budgets. The report and browser limits remain unchanged.

STAR 2.7.11b generates the junction and BAM files. The independent comparison uses
the junction output from that exact alignment, not another alignment's support counts.

## Failures and recovery

The initial reference request to GENCODE's HTTPS FTP endpoint timed out.
The complete UCSC reference was reachable and passed its published MD5.
ENA HTTPS FASTQ transfers timed out or returned HTTP 503.
HTTP transfers succeeded, after one interrupted mate-2 transfer resumed with `curl -C -`.
Both final files passed the published hashes before use.
An official SRA S3 fallback archive was acquired but was not used for the retained workflow.

The first STAR index process exited 137 after 102.73 seconds.
Its sampled peak RSS was 15,773,290,496 bytes.
The account's ancestor cgroup limited total memory to 17,179,869,184 bytes.
Its `oom_kill` counter increased to one. Host RAM was not an available allocation.
The original 48 GiB watchdog therefore could not protect that smaller shared limit.
The retry used the smaller index settings above, a 12 GB STAR budget, and a 14 GiB watchdog.
Alignment, report, and browser measurements then ran sequentially.

The first custom browser preflight tried to use a detail selector from the overview.
The selector was correctly hidden there. The audit was fixed to follow the visible
overview link first. The corrected preflight passed in 3.26 seconds.
This was an audit-harness error, not an application change.

## Controlled observations

STAR aligned 156 artificial mate pairs across all 39 distinct catalog junctions.
All 17 defining junctions and all configured reference junctions have four unique
STAR contributions. The BAM contains 312 selected physical records.
Shared catalog references produce 392 memberships and embedded rows.
This illustrates why alignment-record counts must not replace STAR support.
All these controlled records are unique mappings. They establish no biological prevalence.

The independent inspector matched all 49 TSV rows and all embedded memberships in
both BAM and CRAM reports. Both TSVs also matched TSV-only output byte for byte.
Repeated installed-console runs produced identical provenance JSON bytes.
The reference-only process verified all 455 M5 checksums and read 3,209,286,105 bases.
Its validation took 8.81 seconds with 34,816,000 bytes of process peak RSS.
The controlled BAM and CRAM report profiles took 0.22 and 11.42 seconds, respectively.
Their process peaks were 35,127,296 and 195,268,608 bytes.

These are single observations with filesystem caches warmed by input preparation.
They are not cold-disk timings, latency percentiles, or additional capacity promises.
The separate process peak for reference validation includes interpreter and import costs.
It is measured directly and is never calculated by subtracting process peaks.

## Public observations

The complete run produced 193,080 STAR junction rows. All 17 defining junctions
had zero support. Reference junctions had positive support in each catalog gene
(AR, BRAF, EGFR, FGFR2, and MET). Two configured reference rows were zero:
EGFRvIIIb same-donor and EGFRvIIb same-acceptor. The aggregate
[measurement record](../benchmarks/readiness/record-2026-11-04.json) retains every
expected and observed catalog row for both formats.

The FGFR2-E18-C3 same-donor reference junction had two unique and one multimapping
STAR contributions. Its independently inspected alignments included a real
secondary read-2 record with flag 419, NH=2, and CIGAR `15M3676N48M`.
That record skips chr10:121480022–121483697, with negative-strand STAR motif 2.
This is reference-context evidence, not a positive defining-variant observation.
No custom auxiliary catalog was needed.

Both format inventories contained 343 selected physical records: 342 unique and
one multimapping. They produced 500 displayed memberships because a record may
support several catalog rows. The inspector verified seven displayed mate groups
against raw reciprocal mate coordinates. Eight complete selected mate groups existed
globally; a mate group need not appear together within one variant view.
All 500 memberships had agreeing transcript-strand evidence.
The audit compares these identities and memberships directly; it never substitutes
alignment-record counts for STAR support.

An independent method review exposed an insufficient early comparison based only
on name, mapping class, and membership. The final inspector also compares chromosome,
coordinates, CIGAR, MAPQ, NH, flags, read group, and hit index. It validates each
shown pair against raw mate fields and compares the complete sequence dictionaries.
All public and controlled inventories were rerun after this correction.
Five deliberate corruptions are rejected: wrong STAR count, missing physical record,
wrong pair association, same-name mate substitution, and wrong provenance schema.
Unmodified and restored controls pass.

## Measured capacity and offline behavior

The following values are separate measured stages. MiB means 1,048,576 bytes.
Report and reference peaks are kernel-reported process RSS. Setup and browser peaks
are sampled process-tree RSS; shared pages can be counted more than once.
Neither metric is total allocated account memory.

| Stage | Wall time (s) | Peak RSS (MiB) | Outcome |
| --- | ---: | ---: | --- |
| Complete reference index, successful retry | 1,739.58 | 12,637.12 | Passed 6-hour / 14 GiB setup limits |
| Complete public alignment | 2,174.51 | 9,852.69 | Passed 60-minute / 14 GiB setup limits |
| Public BAM sorting | 372.61 | 404.71 | Passed 600-second / 2 GiB setup limits |
| Public CRAM conversion | 166.93 | 283.28 | Passed 900-second / 2 GiB setup limits |
| Public reference-only process | 9.07 | 35.09 | All 455 sequences and M5 checks passed |
| Public BAM report | 1.89 | 39.50 | Passed 120-second / 2 GiB report limits |
| Public CRAM report | 10.95 | 52.50 | Passed 120-second / 2 GiB report limits |
| Controlled BAM report | 0.22 | 33.50 | Passed report limits |
| Controlled CRAM report | 11.42 | 186.22 | Passed report limits |
| Public BAM offline browser | 4.06 | 839.59 | Passed 90-second / 4 GiB browser limits |
| Public CRAM offline browser | 3.89 | 841.13 | Passed browser limits |
| Controlled BAM offline browser | 3.25 | 841.18 | Passed browser limits |
| Controlled CRAM offline browser | 3.58 | 843.45 | Passed browser limits |

Full-reference validation itself took 9.02 seconds in the public reference-only
process. Inside the public CRAM report, validation took 8.67 seconds; other extraction
took 0.51 seconds and rendering took 0.12 seconds. Public BAM extraction and rendering
took 0.07 and 0.12 seconds. Complete public HTML files were 1,033,151 and 1,033,152 bytes.
The separate exhaustive inventory scans took about 89 seconds each and passed their
300-second audit limits. These scan costs are not report-generation costs.

All four offline reports were copied out before the entire source directory was
moved away. Chromium ran with networking disabled and external requests blocked.
All 17 navigation destinations, mapping and junction-support filters, pair layout,
read inspection, unchanged retained evidence, and JavaScript-disabled inspection passed.
The public multimapping record was visible with All, excluded by Unique, restored
with All, and inspected with NH=2. No external request or browser error occurred.
The slowest opening was 0.49 seconds; the slowest interaction was 0.13 seconds,
within the published 30-second and 10-second limits. Selected views fit one page;
pagination remains covered by the required standard browser suite and capacity audit.

These measurements do not extend the published complete-evidence envelope of
10,000 records. They establish one public workflow and separate controlled cases
on Linux with Chromium 156.0.8078.4, Playwright 1.64.1, and Node 24.18.0.
The machine used AMD EPYC-Rome processors and two STAR threads. Setup, reporting,
and browser phases were serialized under the shared 16 GiB account limit.

## Installed artifacts and compatibility checks

The recommended channel installation succeeded in a fresh environment:

```bash
conda create --yes --prefix "$audit/channel-env" \
  --channel conda-forge --channel bioconda MOMA-AUH::sjsift
PYTHON="$python" bash benchmarks/readiness/channel.sh "$audit" "$audit/channel-env"
```

The published package's public BAM and CRAM TSV files matched the candidate byte
for byte. Physical HTML inventories, strand statuses, and pair associations also
matched. Their CLI wall times were 0.71 and 10.29 seconds, within the report limits.
The unavailable provenance option is a known artifact difference described above.

The candidate build and `twine check` passed. The local suite passed all 477 tests.
The extracted source archive passed all 477 tests and installed BAM/CRAM smoke checks.
The required standard browser runner passed against the installed candidate wheel.
These checks cover the documented CLI, exact TSV layout, schema-3 catalog handling,
error exits, additive outputs, membership/strand rules, and self-contained HTML.
The representative workflows additionally verify both TSVs against independent STAR
joins, exact physical HTML evidence, and compatible full-reference CRAM behavior.

Optional schema-1 provenance was included. Its catalog and junction fingerprints,
record/membership counts, and all 455 reference checks passed independent inspection.
Repeated installed-console runs produced byte-identical JSON. TSV-only provenance
and optional report outputs preserved the same TSV values.

Required CI checks for the audit pull request are `tests (3.11)`, `tests (3.12)`,
`tests (3.13)`, `tests (3.14)`, `package`, and `conda`. They must all pass on the exact
reviewed head before merge. The pull request linked from
[issue #56](https://github.com/MOMA-AUH/sjsift/issues/56) is the final-head CI evidence;
this local record does not predeclare an unobserved CI result. The package gate
includes the installed browser and extracted-source checks. The Conda gate builds
and tests the candidate, separately from the existing published package.

## Readiness decision

**Ready for v1 release preparation within the documented compatibility and capacity
limits, subject to the required final-head CI merge gates above.** No local workflow
or compatibility failure remains unresolved. The index resource failure was diagnosed,
its successful bounded retry is retained, and corrected audit checks were rerun.
The audit-only change does not modify the candidate runtime or catalog.

This decision does not establish clinical sensitivity, disease prevalence, or cohort
coverage. The public sample supplies real zero observations, reference support,
paired evidence, and one selected multimapping record. All positive defining-variant
observations come from explicitly separate artificial controls. Other platforms,
browsers, and workloads beyond the published envelope remain outside this evidence.

A later release change must select the version, rebuild and identify release artifacts,
review the contract, and pass its own required gates. This readiness change does not
publish v1, enable retrospective v1 guarantees for 0.2.2, or authorize hidden changes
to the measured application.
