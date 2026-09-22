# sjsift v0.1 decisions

This log records the decisions made during the design interview. It is intentionally lighter than a full architectural-decision-record system; the one decision that met the ADR threshold is linked below.

## Product and biological model

1. **Narrow quantifier, not a discovery workflow.** sjsift reads a predefined catalog and one STAR junction file. It does not discover events or inspect alignments.

2. **One variant equals one defining junction.** Every v0.1 variant has exactly one assembly-specific chromosome/start/end/strand identity. The initial cases do not require multi-junction aggregation.

3. **Catalog assembly is declared but not verified.** Every catalog has one required `genome_assembly` label. `SJ.out.tab` does not carry assembly metadata, so compatibility remains a documented user responsibility.

4. **Chromosome identifiers are literal.** `chr7` and `7` are different. Silent normalization was rejected because aliases become ambiguous for non-primary contigs.

5. **Coordinates are STAR-native.** Definitions store the first and last intronic bases using 1-based inclusive coordinates. Exonic and BED-style alternatives were rejected to avoid conversion and off-by-one errors.

6. **Strand is exact.** Catalog `+`/`-` matches STAR `1`/`2`; STAR `0` does not match. This deliberately excludes definitions whose motifs make STAR report undefined strand.

7. **STAR counts are preserved.** Output contains unique support, multimapping support, and their unweighted arithmetic sum. sjsift does not recount, filter, deduplicate, or fractionally weight evidence.

8. **Every catalog variant is reported.** A missing exact match produces three zero counts without a per-variant warning. Omitting unsupported entries and adding a redundant status column were rejected.

## Catalog and supported definitions

9. **TOML is the public definition format.** TOML cleanly separates catalog metadata from repeated typed variants and is available through Python 3.11's standard library. See [ADR 0001](adr/0001-use-toml-for-variant-definition-catalogs.md).

10. **The schema is closed and versioned.** `schema_version = 1`, one assembly label, and exactly five fields per variant allow typo detection and future explicit migrations.

11. **IDs and defining junctions are each unique.** Aliases, merging, and last-entry-wins behavior were rejected. Ambiguity is a validation error.

12. **A reference catalog is committed but never implicit.** [`definitions/grch38.toml`](../definitions/grch38.toml) remains a required explicit input. Building definitions into the application or making `--definitions` optional was rejected.

13. **Six EGFR variants are supported.** The v0.1 reference catalog includes `EGFRvII`, `EGFRvIIb`, `EGFRvIII`, `EGFRvIIIb`, `EGFRvIVa`, and `EGFRvIVb` as first-party sjsift definitions.

14. **`METx14del` is the stable identifier.** It denotes the MET exon 13-to-exon 15 RNA junction and does not imply genomic deletion of exon 14.

15. **Reference coordinates use GRCh38 and `chr` prefixes.** The committed catalog also contains `METx14del` and `ARv7`. All eight definitions are expected to produce STAR plus strand on the standard GRCh38 sequence.

## CLI and report

16. **One command, no subcommands.** The public form is `sjsift --junctions PATH --definitions PATH [-o PATH]`. A one-member subcommand hierarchy and positional arguments were rejected.

17. **TSV goes to stdout by default.** `--output` writes a new file, existing files are not overwritten, and diagnostics remain on stderr.

18. **Output is fixed and ordered by the catalog.** Each row carries the variant identity, assembly and junction fields, followed by unique, multimapping, and total support. A sample name is not guessed from the input filename.

19. **Expected user/input errors exit 2.** Success, including warnings and all-zero results, exits 0. Unexpected internal failures exit 1.

20. **Compatibility warnings stay narrow.** An empty STAR file or complete absence of catalog chromosomes produces one warning. Individual zero-support variants do not warn.

## Engineering and delivery

21. **Python 3.11+, setuptools, and `src/` layout.** These retain the useful packaging shape of skua and permit standard-library TOML parsing.

22. **No runtime dependencies.** `argparse`, `tomllib`, and the remaining required facilities are in the standard library. Frameworks are added only in response to a demonstrated need.

23. **The CLI is the supported external seam.** Internal catalog, quantification, and report modules remain small and cohesive; v0.1 does not promise a public Python library interface.

24. **CI covers every currently supported Python minor version.** Pytest runs on Python 3.11, 3.12, 3.13, and 3.14. A separate required package check builds, validates, installs, imports, and smoke-tests the CLI. No formatter, linter, type checker, OS matrix, or coverage gate is initially required.

25. **All changes use pull requests into protected `master`.** Required CI must pass; direct pushes, force pushes, and branch deletion are blocked. No reviewer approval is required.

26. **Squash is the only merge method.** Merge commits and rebase merges are disabled, merged head branches are deleted automatically, and `master` remains linear.

27. **Conda is the initial distribution route.** Tags build and test a pure-Python `noarch` package before publishing it to the `MOMA-AUH` Anaconda organization. Manual workflow runs do not publish.

## Recorded assumptions

- Input is standard nine-column STAR `SJ.out.tab`, already subject to the user's STAR filtering and genome-index choices.
- Reference-catalog coordinates assume the standard GRCh38 sequence and `chr`-prefixed contigs.
- Changing genome-index composition can change unique versus multimapping support even when coordinates are unchanged.
- The initial GitHub repository and `master` branch already exist; repository rules and merge settings still require external configuration.
- The `MOMA-AUH` organization can supply an Anaconda upload token and permit the documented repository settings.

## Unresolved implementation blockers

None. The matching, counting, schema, CLI, output, validation, test, packaging, and repository-workflow contracts are sufficiently defined for implementation.

## Later, not v0.1

- Multi-junction variant definitions and aggregation semantics
- Definition aliases or multiple labels for one junction
- Additional assemblies or maintained catalogs
- Explicit chromosome-alias maps or coordinate liftover
- Opt-in handling of STAR strand `0`
- Standard-input or compressed-file support
- Multiple samples per invocation
- BAM/CRAM processing or alignment-level filtering
- Reference-versus-alternative junction ratios, thresholds, or classifications
- Catalog discovery, download, update, or built-in catalog selection
- Workflow-engine adapters
- Visualization
- Public Python library guarantees
- Additional CI platforms, linting, formatting, typing, or coverage gates when justified by actual maintenance needs
