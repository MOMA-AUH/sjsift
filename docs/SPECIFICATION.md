# sjsift specification

## Purpose and scope

`sjsift` quantifies a predefined set of known RNA splice variants in one sample by exact matching against a STAR `SJ.out.tab` file. Each known splice variant has exactly one defining splice junction. The tool copies STAR's unique and multimapping support counts, calculates their sum, and produces one deterministic TSV row per catalog entry.

The domain terms in [`CONTEXT.md`](../CONTEXT.md) are normative. External format behavior described below comes from STAR; all matching, validation, CLI, and reporting rules are sjsift-specific design decisions.

This document describes the current source checkout, including catalog schemas
1 and 2, reference-junction context, and the optional HTML report. It consolidates
the original MVP specification and design decisions. Released v0.1.3 supports
gzip input and schema 1; reference context, schema 2, and HTML output require a
checkout containing these extensions until a new release includes them.

The [reference-context guide](REFERENCE_JUNCTION_CONTEXT.md) explains comparator
selection, the [curation record](REFERENCE_JUNCTION_CURATION.md) documents the
GRCh38 evidence, and the [HTML guide](HTML_REPORT.md) explains report navigation
and interpretation.

## Supported inputs

Exactly two inputs are required:

1. One STAR `SJ.out.tab` file for one sample.
2. One TOML variant-definition catalog.

Both inputs contain UTF-8 text addressed by filesystem paths. As of v0.1.3, STAR junction input may be plain text or gzip-compressed; gzip is detected by its file header regardless of filename extension and decompressed as a stream. Catalog input remains plain text. Invalid or truncated gzip data is an input error. Standard input, BAM, CRAM, SAM, other compression formats, directories, URLs, and multiple-sample input are not supported.

### STAR junction file

STAR describes `SJ.out.tab` as a tab-delimited file of collapsed splice junctions with these nine columns ([STAR manual](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/extras/doc-latex/STARmanual.tex#L441-L457)):

| Column | Meaning | Accepted values |
|---:|---|---|
| 1 | Chromosome identifier | Non-empty string |
| 2 | First intronic base | Positive integer, 1-based |
| 3 | Last intronic base | Positive integer, 1-based and not less than column 2 |
| 4 | Strand | `0` undefined, `1` plus, `2` minus |
| 5 | Intron motif | Integer `0` through `6` |
| 6 | Annotation status | `0` or `1` |
| 7 | Uniquely mapping reads crossing the junction | Non-negative integer |
| 8 | Multimapping reads crossing the junction | Non-negative integer |
| 9 | Maximum spliced-alignment overhang | Non-negative integer |

Every non-empty input line must contain exactly nine tab-separated fields. An empty file is valid. More than one STAR row matching the same defining or reference junction is an error because a standard `SJ.out.tab` is already collapsed. Non-target duplicate rows are not rejected, but all rows undergo field validation.

Columns 5, 6, and 9 are validated but do not affect matching or output. In two-pass STAR output, annotation status can include junctions found during the first pass and must not be interpreted as “present in the supplied GTF” ([STAR manual](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/extras/doc-latex/STARmanual.tex#L451-L457)).

`SJ.out.tab` is filtered upstream by STAR's `outSJfilter*` settings. Consequently, a missing row means that the junction was not present in this already-filtered file; it does not prove that no alignment in the original data supported the junction ([STAR parameter definitions](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/parametersDefault#L500-L533)).

### Variant-definition catalog

The catalog uses a closed, versioned TOML schema. TOML groups catalog metadata
and typed variant entries in a format that can be reviewed by humans and parsed
by Python's standard library; see [ADR 0001](adr/0001-use-toml-for-variant-definition-catalogs.md).

A minimal schema-1 catalog is:

```toml
schema_version = 1
genome_assembly = "GRCh38"

[[variants]]
id = "METex14"
chromosome = "chr7"
intron_start = 116771655
intron_end = 116774880
strand = "+"
```

Top-level fields:

| Field | Requirement |
|---|---|
| `schema_version` | Required integer; `1` or `2` |
| `genome_assembly` | Required, non-empty string without tabs or line breaks; applies to every entry |
| `variants` | Required, non-empty array of variant tables |

Every variant has these fields:

| Field | Requirement |
|---|---|
| `id` | Non-empty string without tabs or line breaks |
| `chromosome` | Non-empty string without tabs or line breaks, matched literally |
| `intron_start` | Positive integer using STAR's 1-based inclusive convention |
| `intron_end` | Positive integer, not less than `intron_start`, using the same convention |
| `strand` | `+` or `-` |

Schema 1 permits exactly those five fields and provides no reference context.
Schema 2 additionally requires `reference_junctions`, an array of tables that
may be empty (`reference_junctions = []`). Each reference has exactly `role`,
`chromosome`, `intron_start`, `intron_end`, and `strand`. Its junction fields
follow the same rules as a defining junction; `role` is a non-empty string
without tabs or line breaks. See the [schema-2 example](REFERENCE_JUNCTION_CONTEXT.md#catalog-schema-version-2).

Unknown and missing fields are errors, preventing silent typos. Boolean values
are not accepted as integers. Variant IDs and defining-junction tuples
`(chromosome, intron_start, intron_end, strand)` must each be unique across the
catalog. Duplicate definitions are rejected rather than merged or treated as
aliases.

Reference roles must be unique within each variant. A reference may not duplicate
its own variant's defining junction. References may otherwise share coordinates,
including across variants or with another variant's defining junction; the same
observed counts are reported for each use. Roles such as `same_donor` and
`same_acceptor` express curation intent in transcript orientation. The loader
accepts custom roles and does not enforce shared splice sites, chromosome, or
strand between a reference and its variant.

The required assembly label makes the catalog self-describing, but `SJ.out.tab` contains no corresponding assembly metadata. sjsift therefore cannot verify assembly compatibility. Users are responsible for ensuring that the STAR genome and catalog share the same reference coordinate namespace.

## Versioned GRCh38 reference catalog

[`definitions/grch38.toml`](../definitions/grch38.toml) uses schema 2 and contains
17 defining junctions across AR, BRAF, EGFR, FGFR2, and MET:

| Gene | Variant IDs |
|---|---|
| AR | `ARv7`, `ARv567es` |
| BRAF | `BRAFdel2-10`, `BRAFdel2-8`, `BRAFdel3-10`, `BRAFdel3-8`, `BRAFdel4-10`, `BRAFdel4-8` |
| EGFR | `EGFRvIVa`, `EGFRvIII`, `EGFRvIIIb`, `EGFRvIVb`, `EGFRvII`, `EGFRvIIb` |
| FGFR2 | `FGFR2-E18-C3` |
| MET | `METex7-8`, `METex14` |

It provides 32 reference entries representing 22 distinct reference junctions.
Coordinates use GRCh38 with literal `chr`-prefixed identifiers, including both
plus- and minus-strand events. The [curation record](REFERENCE_JUNCTION_CURATION.md)
pins transcript choices, explains project-specific EGFR suffixes, and records
annotation checksums and exon evidence. Versioned RefSeq `NM_` names are catalog
comments and documentation; transcript IDs and exon names are not schema fields.
`METex14` remains the project identifier for MET exon 14 skipping.

The catalog is explicit input, never an application default or an automatically
updated download. Users may copy or replace it with another schema-compatible
catalog. The [v0.1.3 catalog](https://github.com/MOMA-AUH/sjsift/blob/v0.1.3/definitions/grch38.toml)
remains available for the released schema-1-only application.

## Coordinate, chromosome, and strand rules

- `intron_start` and `intron_end` are the first and last intronic bases, both 1-based and inclusive. They correspond exactly to STAR columns 2 and 3; no conversion is performed.
- Chromosome identifiers match exactly. For example, `chr7` does not match `7`.
- Catalog `+` matches STAR strand `1`; catalog `-` matches STAR strand `2`.
- STAR strand `0` never matches a catalog entry.
- Motif, annotation status, and overhang are not part of the match key.
- Coordinate liftover, chromosome aliases, reference-sequence normalization, and tolerance windows are not supported.

STAR derives the junction strand from the intron motif; noncanonical motif `0` produces undefined strand `0` ([STAR source](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/ReadAlign_outputTranscriptSJ.cpp#L42-L50)). A catalog therefore supports only definitions expected to have a STAR-defined strand. All entries in the GRCh38 reference catalog meet that requirement on the standard GRCh38 reference.

## Matching and counting

For every defining and reference junction, sjsift:

1. Forms the exact key `(chromosome, intron_start, intron_end, STAR strand code)`.
2. Looks for the identical key in the STAR file.
3. If found, copies column 7 to `unique_support` and column 8 to `multimapping_support`.
4. Calculates `total_support = unique_support + multimapping_support`.
5. If not found, assigns zero to all three support fields.

No read is opened, recounted, deduplicated, filtered, fractionally weighted, or reclassified. STAR increments multimapping support without `1/N` weighting ([STAR source](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/ReadAlign_outputTranscriptSJ.cpp#L45-L51)); sjsift preserves that meaning. STAR treats paired-end input as one read for its mapping statistics, so sjsift reports STAR-defined support rather than inventing its own fragment semantics ([STAR manual](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/extras/doc-latex/STARmanual.tex#L286-L291)).

TSV results follow catalog order, regardless of STAR row order. Configured
references remain separate observations, ordered as declared within each
variant. A junction shared by several references is looked up once and its
counts are copied to each corresponding report row; these rows must not be
summed as independent evidence. Reference support does not change defining-junction
counts. A missing reference configuration differs from a configured junction
with zero support.

There is no threshold, score, positive/negative call, multi-junction aggregation,
or inferred whole-transcript isoform fraction. In particular, summing donor and
acceptor reference counts would count overlapping splice opportunities. Input
genome-index composition and STAR filtering can affect the observed counts;
sjsift preserves the input evidence.

## CLI contract

```text
sjsift --junctions PATH --definitions PATH [-o PATH] [--context-output PATH] [--html-output PATH]
sjsift --help
sjsift --version
```

- `--junctions PATH` is required and names one `SJ.out.tab` file.
- `--definitions PATH` is required and names one TOML catalog.
- `-o PATH` and `--output PATH` are equivalent and optional.
- Without `--output`, TSV is written to standard output.
- Diagnostics and warnings are written to standard error; standard output contains only the result table.
- `--context-output PATH` adds a reference-junction TSV.
- `--html-output PATH` adds an offline HTML report independently of `--context-output`.
- Both optional reports work with either a main TSV file or main TSV on stdout.
- All output paths must differ and must not already exist; parent directories must exist.
- Inputs cannot be provided via standard input.
- There are no subcommands.

## Output schemas

### Main TSV

The main output is UTF-8, tab-delimited text with one header and exactly one row per catalog entry:

| Column | Meaning |
|---|---|
| `variant_id` | Catalog `id` |
| `genome_assembly` | Catalog assembly label |
| `chromosome` | Exact catalog chromosome identifier |
| `intron_start` | Catalog STAR-native intron start |
| `intron_end` | Catalog STAR-native intron end |
| `strand` | Catalog `+` or `-` |
| `unique_support` | STAR column 7, or zero when unmatched |
| `multimapping_support` | STAR column 8, or zero when unmatched |
| `total_support` | Arithmetic sum of the preceding two fields |

No sample identifier is inferred from the input filename. Numeric values are emitted as base-10 integers without quoting.

### Reference-context TSV

The context TSV has one header and one row per configured reference junction,
including zero-support junctions. Its columns, in order, are:

```text
variant_id	genome_assembly	context_role	chromosome	intron_start	intron_end	strand	unique_support	multimapping_support	total_support
```

`variant_id` identifies the owning variant, `context_role` is the declared role,
and the coordinates and counts describe that reference junction. A catalog
without references produces a header-only context TSV. Schema 1 remains usable
with either optional report.

### HTML report

The HTML report presents the same results in a comparison table and per-variant
details, with exact counts and separate reference roles. It defaults to sorting
by defining-junction unique support, then total support, with catalog order
breaking ties. Interactive search, filters, and navigation change the view;
they do not alter either TSV. All assets are embedded for offline use, and all
results remain readable without JavaScript. See [HTML report](HTML_REPORT.md)
for evidence descriptions, chart scales, provenance, and browser behavior.

### Output failure handling

All input is validated before report creation. All requested output files are
opened exclusively before writing; file reports are flushed and closed before
the main TSV is emitted to stdout. If an output fails, newly created files are
removed where possible. Existing files, including inputs supplied as output
paths, are never overwritten. Bytes already emitted to stdout cannot be
retracted, and cleanup is not transactional against process termination or
power loss.

## Validation, warnings, and exit status

Expected user-facing errors include unreadable inputs, invalid TOML, unsupported schema version, missing or unknown catalog fields, invalid catalog values, duplicate IDs or junctions, malformed STAR rows, invalid STAR field values, duplicate matching STAR records, and refusal to overwrite an existing output file. These failures must identify the file and, where applicable, the catalog entry or STAR line number; they must not emit a Python traceback by default.

If the STAR file is empty, or if none of the defining or reference junctions' chromosome identifiers occurs anywhere in it, sjsift emits one compatibility warning to standard error and still writes zero-support rows. The HTML also displays this warning when requested. There is no warning for individual zero-support variants, and chromosome-name overlap alone does not prove assembly compatibility.

Exit statuses are:

| Status | Meaning |
|---:|---|
| `0` | Successful report, including reports containing warnings or only zero counts |
| `2` | Expected CLI, input, validation, or output-path error |
| `1` | Unexpected internal failure |

## Design and maintenance

The supported interface is one CLI with explicit named inputs and additive
outputs. Internal Python modules do not promise a public library API. Python
3.11 or newer provides standard-library TOML parsing; the package uses a `src/`
layout, setuptools, and no runtime dependencies. HTML assets are packaged with
the application rather than fetched from a CDN.

The [development guide](DEVELOPMENT.md) owns contributor setup, repository
policy, and release procedures. CI tests Python 3.11–3.14 and builds and installs
distributions, including an HTML rendering smoke test from the installed wheel.
Conda releases build and test a `noarch` package before publication to MOMA-AUH;
manual workflow dispatch builds and tests without publishing. Required checks
and squash-merge policy are documented in that guide.

## Verification

Automated checks cover the following contracts:

| Contract | Verification |
|---|---|
| Both catalog schemas, field validation, uniqueness, and reference constraints | `tests/test_catalog.py` |
| Committed catalog and reproducible curated reference coordinates | `tests/test_reference_catalog.py`, `docs/curation/gencode-v49-exons.tsv` |
| STAR validation, exact matching, gzip, duplicate targets, shared references, and zero counts | `tests/test_quantify.py` |
| TSV columns and order, CLI streams and statuses, compatibility warnings, and output refusal | `tests/test_cli.py` |
| HTML evidence, escaping, deterministic rendering, ordering, and report failure cleanup | `tests/test_html_report.py` |
| Version agreement and release configuration | `tests/test_release.py` |
| Wheel/sdist build, metadata, installed CLI, and packaged HTML assets | CI and publish workflows |

Fixtures and synthetic report examples exercise software behavior; they are
not clinical validation samples.

## Current boundaries

The application does not implement:

- BAM, CRAM, or SAM processing, alignment recounting, deduplication, or filtering;
- novel splice-event discovery or multiple defining junctions per variant;
- standard-input ingestion or multiple samples per invocation;
- exon-coordinate conversion, liftover, chromosome aliases, or fuzzy matching;
- scores, thresholds, clinical calls, or variant/reference fractions;
- inferred sample identity or independent verification of the input assembly;
- automatic catalog discovery, annotation download, or catalog updates;
- transcript-aware splice diagrams or structured transcript/exon metadata;
- workflow-engine integration, plugin interfaces, or a public Python library API.
