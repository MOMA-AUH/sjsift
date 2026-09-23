# sjsift v0.1 MVP specification

## Purpose and scope

`sjsift` quantifies a predefined set of known RNA splice variants in one sample by exact matching against a STAR `SJ.out.tab` file. Each known splice variant has exactly one defining splice junction. The tool copies STAR's unique and multimapping support counts, calculates their sum, and produces one deterministic TSV row per catalog entry.

The domain terms in [`CONTEXT.md`](../CONTEXT.md) are normative. External format behavior described below comes from STAR; all matching, validation, CLI, and reporting rules are sjsift-specific design decisions.

## Supported inputs

Exactly two inputs are required:

1. One STAR `SJ.out.tab` file for one sample.
2. One TOML variant-definition catalog.

Both inputs are UTF-8 text files addressed by filesystem paths. Standard input, BAM, CRAM, SAM, compressed inputs, directories, URLs, and multiple-sample input are not supported in v0.1.

### STAR junction file

STAR describes `SJ.out.tab` as a tab-delimited file of collapsed splice junctions with these nine columns ([STAR manual](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/extras/doc-latex/STARmanual.tex#L441-L457)):

| Column | Meaning | Accepted v0.1 values |
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

Every non-empty input line must contain exactly nine tab-separated fields. An empty file is valid. More than one STAR row matching the same catalog junction is an error because a standard `SJ.out.tab` is already collapsed.

Columns 5, 6, and 9 are validated but do not affect matching or output. In two-pass STAR output, annotation status can include junctions found during the first pass and must not be interpreted as “present in the supplied GTF” ([STAR manual](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/extras/doc-latex/STARmanual.tex#L451-L457)).

`SJ.out.tab` is filtered upstream by STAR's `outSJfilter*` settings. Consequently, a missing row means that the junction was not present in this already-filtered file; it does not prove that no alignment in the original data supported the junction ([STAR parameter definitions](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/parametersDefault#L500-L533)).

### Variant-definition catalog

The catalog uses TOML with this closed schema:

```toml
schema_version = 1
genome_assembly = "GRCh38"

[[variants]]
id = "EGFRvIII"
chromosome = "chr7"
intron_start = 55019366
intron_end = 55155829
strand = "+"

[[variants]]
id = "METx14skip"
chromosome = "chr7"
intron_start = 116771655
intron_end = 116774880
strand = "+"

[[variants]]
id = "ARv7"
chromosome = "chrX"
intron_start = 67686127
intron_end = 67694672
strand = "+"
```

Top-level fields:

| Field | Requirement |
|---|---|
| `schema_version` | Required integer; must equal `1` |
| `genome_assembly` | Required, non-empty string without tabs or line breaks; applies to every entry |
| `variants` | Required, non-empty array of variant tables |

Each variant table has exactly these required fields:

| Field | Requirement |
|---|---|
| `id` | Non-empty string without tabs or line breaks |
| `chromosome` | Non-empty string without tabs or line breaks, matched literally |
| `intron_start` | Positive integer using STAR's 1-based inclusive convention |
| `intron_end` | Positive integer, not less than `intron_start`, using the same convention |
| `strand` | `+` or `-` |

Unknown fields are errors so misspellings cannot be ignored silently. Variant IDs must be unique. The tuple `(chromosome, intron_start, intron_end, strand)` must also be unique; aliases for the same junction are not supported.

The required assembly label makes the catalog self-describing, but `SJ.out.tab` contains no corresponding assembly metadata. sjsift therefore cannot verify assembly compatibility. Users are responsible for ensuring that the STAR genome and catalog share the same reference coordinate namespace.

## Versioned GRCh38 reference catalog

[`definitions/grch38.toml`](../definitions/grch38.toml) is the v0.1 reference catalog. It contains the agreed GRCh38 definitions for six known EGFR variants, `METx14skip`, and `ARv7`. The GRCh38 STAR-style junctions for `EGFRvIII`, `METx14skip`, and `ARv7` are also documented by Illumina ([DRAGEN 4.4 documentation](https://help.connected.illumina.com/dragen/dragen-v4.4/product-guide/dragen-v4.4/dragen-rna-pipeline/splice-variant-caller#knowns-list)).

| Variant ID | Chromosome | Intron start | Intron end | Strand |
|---|---|---:|---:|:---:|
| `EGFRvIVa` | `chr7` | 55200414 | 55205255 | `+` |
| `EGFRvIII` | `chr7` | 55019366 | 55155829 | `+` |
| `EGFRvIIIb` | `chr7` | 55109959 | 55155829 | `+` |
| `EGFRvIVb` | `chr7` | 55200414 | 55202516 | `+` |
| `EGFRvII` | `chr7` | 55161632 | 55171174 | `+` |
| `EGFRvIIb` | `chr7` | 55161632 | 55170306 | `+` |
| `METx14skip` | `chr7` | 116771655 | 116774880 | `+` |
| `ARv7` | `chrX` | 67686127 | 67694672 | `+` |

The supported EGFR variants are `EGFRvII`, `EGFRvIIb`, `EGFRvIII`, `EGFRvIIIb`, `EGFRvIVa`, and `EGFRvIVb`. `METx14skip` is the stable project identifier for MET exon 14 skipping.

The catalog is explicit input, not a hidden default. Users may copy or replace it with another schema-compatible catalog.

## Coordinate, chromosome, and strand rules

- `intron_start` and `intron_end` are the first and last intronic bases, both 1-based and inclusive. They correspond exactly to STAR columns 2 and 3; no conversion is performed.
- Chromosome identifiers match exactly. For example, `chr7` does not match `7`.
- Catalog `+` matches STAR strand `1`; catalog `-` matches STAR strand `2`.
- STAR strand `0` never matches a v0.1 catalog entry.
- Motif, annotation status, and overhang are not part of the match key.
- Coordinate liftover, chromosome aliases, reference-sequence normalization, and tolerance windows are not supported.

STAR derives the junction strand from the intron motif; noncanonical motif `0` produces undefined strand `0` ([STAR source](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/ReadAlign_outputTranscriptSJ.cpp#L42-L50)). A v0.1 catalog therefore supports only definitions expected to have a STAR-defined strand. All entries in the GRCh38 reference catalog meet that requirement on the standard GRCh38 reference.

## Matching and counting

For each catalog entry, sjsift:

1. Forms the exact key `(chromosome, intron_start, intron_end, STAR strand code)`.
2. Looks for the identical key in the STAR file.
3. If found, copies column 7 to `unique_support` and column 8 to `multimapping_support`.
4. Calculates `total_support = unique_support + multimapping_support`.
5. If not found, assigns zero to all three support fields.

No read is opened, recounted, deduplicated, filtered, fractionally weighted, or reclassified. STAR increments multimapping support without `1/N` weighting ([STAR source](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/source/ReadAlign_outputTranscriptSJ.cpp#L45-L51)); sjsift preserves that meaning. STAR treats paired-end input as one read for its mapping statistics, so sjsift reports STAR-defined support rather than inventing its own fragment semantics ([STAR manual](https://github.com/alexdobin/STAR/blob/b1edc1208d91a53bf40ebae8669f71d50b994851/extras/doc-latex/STARmanual.tex#L286-L291)).

The result order is the catalog order, regardless of STAR row order. There is no threshold, score, positive/negative call, or multi-junction aggregation.

## CLI contract

```text
sjsift --junctions PATH --definitions PATH [-o PATH]
sjsift --help
sjsift --version
```

- `--junctions PATH` is required and names one `SJ.out.tab` file.
- `--definitions PATH` is required and names one TOML catalog.
- `-o PATH` and `--output PATH` are equivalent and optional.
- Without `--output`, TSV is written to standard output.
- Diagnostics and warnings are written to standard error; standard output contains only the result table.
- An existing output path is not overwritten.
- Inputs cannot be provided via standard input in v0.1.
- There are no subcommands.

## Output schema

The output is UTF-8, tab-delimited text with one header and exactly one row per catalog entry:

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

## Validation, warnings, and exit status

Expected user-facing errors include unreadable inputs, invalid TOML, unsupported schema version, missing or unknown catalog fields, invalid catalog values, duplicate IDs or junctions, malformed STAR rows, invalid STAR field values, duplicate matching STAR records, and refusal to overwrite an existing output file. These failures must identify the file and, where applicable, the catalog entry or STAR line number; they must not emit a Python traceback by default.

If the STAR file is empty, or if none of the catalog's chromosome identifiers occurs anywhere in it, sjsift emits one compatibility warning to standard error and still writes zero-support rows. It does not warn for individual zero-support variants.

Exit statuses are:

| Status | Meaning |
|---:|---|
| `0` | Successful report, including reports containing warnings or only zero counts |
| `2` | Expected CLI, input, validation, or output-path error |
| `1` | Unexpected internal failure |

## Acceptance criteria

Version 0.1 is acceptable when automated tests demonstrate that:

- all valid schema fields are parsed and every invalid or duplicate definition described above is rejected;
- all nine STAR fields are validated with useful line-numbered errors;
- exact chromosome, coordinate, and strand matches copy unique and multimapping counts correctly;
- chromosome aliases, off-by-one coordinates, opposite strand, and strand `0` do not match;
- total support is the unweighted arithmetic sum;
- unmatched variants, an empty STAR file, and an all-zero sample retain every catalog row;
- catalog order determines output order;
- the full reference catalog, including all six supported EGFR variants, `METx14skip`, and `ARv7`, is accepted;
- malformed and duplicate matching STAR records fail rather than being silently combined;
- stdout, stderr, output-file refusal, and exit statuses follow the CLI contract;
- a built wheel installs cleanly, imports, and runs `sjsift --help` and `sjsift --version`.

## Explicit non-goals for v0.1

- BAM, CRAM, or SAM processing
- recounting alignments or read-pair/fragment-level processing
- MAPQ, duplicate, SAM-flag, motif, annotation, or overhang filtering
- novel splice-event discovery
- more than one defining junction per variant
- exon-coordinate conversion, liftover, chromosome normalization, or fuzzy matching
- scores, thresholds, classifications, clinical interpretation, or sample-identity checks
- visualization, databases, downloadable genome bundles, or automatic catalog updates
- extensive provenance capture
- workflow-engine integration, plugin systems, or a public Python library interface
- speculative performance optimization
