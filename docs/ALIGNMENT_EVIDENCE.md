# Alignment evidence in offline reports

The current development version accepts a local coordinate-sorted indexed BAM
alongside the required STAR file and schema-3 catalog:

```bash
sjsift --junctions sample.SJ.out.tab --definitions catalog.toml \
  --alignments sample.bam --html-output sample.html --output sample.tsv
```

`--alignment-index PATH` overrides discovery of `sample.bam.bai`, `sample.bai`,
`sample.bam.csi`, or `sample.csi`. `--alignment-limit N` is a positive integer,
default 200, **per junction per mapping class**. Alignment options require
`--alignments`, which requires `--html-output`. Indexes are never created by
sjsift. This initial BAM slice shows defining-junction previews; reference reads
and expansion controls follow in the subsequent implementation milestones.
CRAM is rejected until explicit-reference support is implemented.

## What qualifies

An alignment must have the exact chromosome identifier and a CIGAR `N` matching
the selected intron's first and last intronic bases. With a zero-based reference
cursor `r` and skip length `L`, the STAR interval is `r+1` through `r+L`.
`M`, `=`, `X`, `D`, and `N` advance that cursor; insertions and clipping do not.
A deletion or the space between paired mates is not a splice. There is no MAPQ
or aligned-anchor threshold. Only selected-junction records are retained.

Transcript strand and forward/reverse alignment orientation are distinct:

- `TS:A:+/-` is the standard transcript-strand tag. `XS:A:+/-` follows the
  established RNA-aligner convention; numeric XS scores and other types are
  unusable as strand evidence.
- STAR `jM` is usable only with STAR program provenance in `@PG`, an integer
  array with one valid motif per CIGAR N, and (if supplied) `jI` exactly matching
  the CIGAR introns in order. Motifs 1/3/5 imply plus and 2/4/6 minus; subtract
  the annotation offset of 20 first. Motifs 0/20 give no strand evidence.
  Each selected junction uses its own motif, not a record-wide majority.
- Minimap2 `ts:A:+/-` is interpreted relative to alignment orientation only
  when minimap2 program provenance is available. It is distinct from `TS`.
- A record's `PG` selects its program/ancestor chain. Without a record `PG`,
  local tag semantics require the relevant program in every possible terminal
  header chain. Missing or ambiguous provenance leaves that source unusable.
- All usable sources must agree with the catalog strand. Opposite and mutually
  conflicting evidence have separate exclusion counts; excluded records are
  not embedded. Missing or unusable sources produce **strand unverified** rows,
  visibly distinguished from agreeing evidence. Details describe the sources.

NH alone determines mapping class: an integer 1 is **unique**, an integer above
1 is **multimapping**, and absent, nonpositive, or mistyped values are **unknown**.
Primary status and MAPQ cannot establish uniqueness. Qualifying secondary,
supplementary, duplicate-marked and QC-failure records are included and labeled.

## Counts, identity, and sampling

One row is one physical alignment record. Both mates may appear separately;
read name, read group and read-1/read-2 labels help recognize related records.
No unrelated mate is fetched. Separate records with identical fields remain
separate. A physical record retrieved at several regional queries is processed
only in the first requested region it overlaps; no name-based deduplication is
performed.

Every eligible record is counted before sampling. Deterministic SHA-256 bottom-k
sampling independently caps each junction/class subset, retaining all records
below its cap. Extraction memory for retained records is bounded by those caps;
it does not keep a set of every observed read. Caps limit records, not HTML bytes.
Eligible, embedded and currently shown counts have separate labels, as do strand
totals and exclusion reasons. The initial preview displays up to ten records.

**STAR support remains unchanged** in both TSVs and HTML. STAR's upstream
junction filtering, paired-template semantics, alignment output choices and the
explicit display policy can all produce different alignment-record totals. A
difference alone is not an input error. Do not add reference opportunities or
interpret these records as independent molecules or full-length transcripts.
Evidence not requested is distinct from an observed zero eligible count.

## Validation and portability

A usable local index, coordinate sort declaration, required literal contigs and
coordinate bounds are checked. Multiple distinct `@RG SM` values are rejected;
missing metadata is disclosed. Available `@SQ AS` labels must equal the catalog
assembly. Header labels do not establish sequence equivalence, and anonymous
STAR data cannot prove sample identity or genome assembly. The caller must
supply compatible inputs. Regional decoding errors stop the whole invocation.
All evidence extraction finishes before any output is created. Existing output
refusal, cleanup of newly created files on ordinary errors and stderr/exit-status
conventions are unchanged.

The single HTML contains structural operations, coordinates, CIGAR, MAPQ,
selected strand evidence, flags, read names, read-group/sample identifiers,
annotation and software provenance, and input basenames. It contains no read
sequence, base qualities, reference sequence, original alignment file, or
absolute source paths. Consider the retained raw identifiers before sharing.
All text is escaped and all assets are inline; no browser network access or
source BAM is needed. Counts and provenance remain readable without JavaScript.
The geometry is schematic; nucleotide-level inspection is outside this release.

## Dependency and validation basis

Pysam is a normal runtime dependency (`>=0.24.0,<0.25`) in both Python metadata
and the noarch conda recipe. There is no reads extra or external samtools
requirement. The 0.24 series provides Python 3.11–3.14 wheels; CI exercises the
minimum 0.24.0 on Python 3.11 and current releases on the other supported Python
versions. CI builds a wheel with dependency resolution in a clean environment
and builds/tests the conda recipe using conda-forge and bioconda. Installed checks
generate reports from real synthetic indexed alignments. These are software
fixtures, not clinical validation data.

Primary format/API references checked 2026-10-01:

- [SAM optional tags](https://samtools.github.io/hts-specs/SAMtags.pdf).
- [STAR output attributes](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault),
  and [BAM CIGAR/motif emission](https://github.com/alexdobin/STAR/blob/master/source/ReadAlign_alignBAM.cpp).
- [Minimap2 tag and splice orientation implementation](https://github.com/lh3/minimap2/blob/master/format.c).
- [Pysam 0.24 releases](https://pysam.readthedocs.io/en/latest/release.html)
  and [alignment APIs](https://pysam.readthedocs.io/en/latest/api.html).
