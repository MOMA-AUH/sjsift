# Stable v1 compatibility contract

This policy takes effect with sjsift v1.0.0 and covers stable v1 releases.
It is a release requirement prepared before v1, not an announcement of a release.
Existing v0.x releases retain their pre-1.0 policy, including possible breaking
changes in patch releases. They acquire no retrospective guarantees.

Compatibility means that an accepted invocation and its pinned inputs retain
their documented meaning after an upgrade within v1. Pin the application version
and catalog bytes when exact reproduction is required. Software versions,
dependency versions, diagnostics, and HTML presentation can differ after upgrades.

The [specification](SPECIFICATION.md) defines accepted inputs and output values.
The [alignment guide](ALIGNMENT_EVIDENCE.md), [HTML guide](HTML_REPORT.md), and
[provenance schema](PROVENANCE.md) supply their detailed behavior. This contract
defines which changes to those interfaces require a major release.

## Command and input guarantees

The supported program interfaces are `sjsift` and `python -m sjsift`.
The following option spellings, argument meanings, and dependencies remain
supported throughout v1:

| Option | Stable behavior |
| --- | --- |
| `--junctions PATH` | Required local UTF-8, nine-column STAR `SJ.out.tab`, plain text or gzip detected by file header. |
| `--definitions PATH` | Required explicit local TOML variant-definition catalog using schema 3. |
| `--sample-name NAME` | Nonblank sample name. Required unless one distinct nonblank sample name is available from the alignment header. |
| `-o PATH`, `--output PATH` | Main TSV destination. Without either option, the main TSV goes to stdout. |
| `--context-output PATH` | Additional reference-context TSV file. |
| `--html-output PATH` | Additional self-contained HTML file, independent of the context TSV option. |
| `--provenance-output PATH` | Additional versioned run provenance JSON file, including for TSV-only runs. |
| `--alignment PATH` | Optional local indexed, coordinate-sorted BAM or CRAM. Requires `--html-output`. |
| `--alignment-index PATH` | Explicit local alignment index instead of discovery. Requires `--alignment`. |
| `--reference PATH` | Explicit local indexed FASTA, required for every CRAM, including self-contained CRAM. Requires `--alignment` and applies only to CRAM. |
| `-h`, `--help` | Print help to stdout and exit successfully without inputs. Help layout and prose can change. |
| `--version` | Print `sjsift VERSION` and a newline to stdout and exit successfully without inputs. |

If the alignment header contains sample names, a supplied name must match one
exactly. Without header sample names, an explicit sample name remains valid.
Multiple alignment samples require explicit selection and retain only that
sample's assigned read groups. This selection never changes STAR support.
Alignment indexes and FASTA indexes must already exist. sjsift creates no indexes
and downloads no references. Standard input and remote input URLs remain outside
the supported input interface. Abbreviated long options are not stable spellings.

Reports remain additive. Requesting context, HTML, alignment evidence, or JSON
does not change either TSV's values or format. During report generation, stdout
contains only the main TSV, or is empty when `--output` selects a file.
Diagnostics and warnings go to stderr. Expected errors identify the input or
destination and relevant entry or line without a Python traceback.

| Exit status | Stable meaning |
| ---: | --- |
| `0` | Success, including warnings, zero-support reports, help, and version. |
| `2` | Expected CLI, input, validation, or output error, including a closed stdout pipe. |
| `1` | Unexpected internal failure. |

All output paths must differ, their parent directories must exist, and existing
files are never overwritten. Input validation precedes report creation.
File outputs finish before the main TSV reaches stdout. On an ordinary failure,
sjsift removes files created by that run where possible. Emitted stdout cannot
be recalled. Sudden termination, power loss, and failed filesystem cleanup are
outside this guarantee.

An empty STAR file or no catalog chromosome identifiers in STAR produces one
compatibility warning and a successful zero-support report. HTML also explains
the condition. Individual zero-support variants produce no such warning.
Warnings do not establish an assembly mismatch or change counts.
Additional warning conditions can be introduced without changing success status
or report semantics. Stderr wording is human-readable, not a parsing interface.
Use the provenance warning codes for machine decisions.

## Catalog and quantification guarantees

Schema 3 remains accepted throughout v1 with the field types, validation rules,
and meanings in the [catalog specification](SPECIFICATION.md#variant-definition-catalog).
This includes required annotations, optional exon extents, custom reference roles,
empty reference arrays, and shared reference junctions. Unknown fields remain
errors unless explicitly introduced as optional schema fields in a later release.
Such additions cannot invalidate existing catalogs or change omitted-field behavior.
An older binary is not required to accept fields introduced by a newer binary.
Schemas 1 and 2 remain unsupported by the v1 baseline.

Matching uses the exact chromosome identifier, 1-based inclusive intron start
and end, and transcript strand. `chr7` differs from `7`.
Catalog `+` matches STAR `1`, catalog `-` matches STAR `2`, and STAR `0` never
matches. There is no coordinate conversion, aliasing, or tolerance window.
Motif, annotation status, and overhang are validated but do not affect matching.

Unique and multimapping support retain STAR columns 7 and 8 without weighting.
Total support remains their arithmetic sum. Unmatched junctions retain zero
counts. Reference support remains separate from defining-junction support.
Shared reference rows are not independent evidence. Alignment evidence never
recounts or replaces STAR support. Zero support is not a clinical negative.

## TSV guarantees

The main TSV retains these columns in this exact order:

```text
variant_id	genome_assembly	chromosome	intron_start	intron_end	strand	unique_support	multimapping_support	total_support
```

The context TSV retains these columns in this exact order:

```text
variant_id	genome_assembly	context_role	chromosome	intron_start	intron_end	strand	unique_support	multimapping_support	total_support
```

Both formats retain UTF-8 encoding, tab separators, LF line endings, a header,
and base-10 integer values. The main TSV contains exactly one row per variant
in catalog order. Context rows follow variant order, then declared reference
order, including zero-support references. No references means a header-only
context file. STAR row order does not affect report order.

Adding, removing, renaming, reordering, or reinterpreting columns in either
existing TSV requires a major release. New opt-in output formats can be added
in a minor release while the existing formats remain unchanged.

## Versioned JSON and HTML boundaries

The opt-in output introduced before v1 by `--provenance-output` is covered from
v1.0.0. [Provenance schema 1](PROVENANCE.md) defines its fields, types, meanings,
warning codes, metadata limits, and deterministic serialization.
Consumers must inspect `schema_version`, tolerate unknown fields, and preserve
unknown warnings. Additive fields and warning codes are compatible.
Existing code meanings remain stable. Readable messages can
change while their structured facts retain their meaning.

Field removal or a change in type or meaning requires a new JSON schema version.
A schema increment alone does not permit a breaking change within sjsift v1.
Changing the schema emitted by the existing invocation incompatibly requires a
major application release. A minor release can add an explicit opt-in newer
schema while preserving schema 1 for existing invocations.
The same rule applies to any other opt-in versioned machine-readable output
introduced before v1. Its schema and compatibility checks must be documented
before the stable release.

HTML guarantees cover observable behavior: accurate counts, separate reference
roles, sample and input information, warnings, navigation, filters, and read
inspection as described in the guides. Alignment selection, classification,
membership counts, and their distinction from STAR support retain their meaning.
All eligible alignment evidence remains embedded. The combined variant figure
shows each physical record once and retains its qualifying junction memberships.
Support filtering and accessible labels expose shared memberships, as described
in the [reference-context guide](REFERENCE_JUNCTION_CONTEXT.md#alignment-evidence-and-shared-membership).
Paired viewing preserves individual records and does not reconstruct STAR support.
Reports open through `file:`
URLs without a server, external assets, network access, or original input files.
Without JavaScript, counts, variant details, and native preview inspection remain
readable. Interactive filtering and paging require JavaScript.

HTML is a human interface. DOM structure, element IDs, CSS classes, colors,
spacing, script organization, and embedded internal data are not stable machine
interfaces. Presentation can change while preserving the documented information
and capabilities. HTML bytes need not remain identical across software versions.
Use TSV and versioned JSON for automation.
Internal Python modules, functions, classes, and imports are not a supported
public library API and can change within v1.

## Catalog revisions and release versions

The application version, catalog schema number, and catalog content revision
identify different things. Schema 3 defines the format, not the accuracy or
membership of a particular catalog. The reference catalog's content is versioned
with the repository release tag and commit. A file checksum identifies exact bytes.

Reference catalog additions belong in a minor release. Corrections to coordinates,
annotations, or reference junctions can accompany a patch release when the schema
and application behavior remain compatible. Release notes must list affected
variant IDs, old and new values, the evidence, and the effect on counts or display.
Additions, removals, renames, and reordering also require explicit release notes.
Removing or renaming an existing variant ID requires a major release.
New catalogs can produce different results when users deliberately adopt them.

sjsift always reads the explicit `--definitions` file. An application upgrade
never replaces, edits, downloads, or silently substitutes a user's pinned catalog.
Release-tagged catalog copies remain available. Pin a tag or commit and retain
the file bytes instead of downloading from moving `master` for each run.
Review the catalog diff separately from the application upgrade.

A patch release fixes defects while preserving this contract. A minor release
adds compatible capabilities. A major release is required for these changes:

- Remove or rename an accepted CLI option, change its meaning or default, or add a required input to an existing invocation.
- Reject or incompatibly reinterpret a previously valid STAR, schema-3 catalog, alignment, or reference input.
- Change exact matching, STAR support semantics, existing TSV columns or ordering, exit statuses, output destinations, or overwrite guarantees.
- Break the existing versioned JSON interface, remove documented HTML capabilities, or require online resources for reports.
- Drop a supported runtime environment, including a supported Python minor version or supported platform.

The baseline requires Python 3.11 or newer, with Python 3.11–3.14 tested in Linux
CI and installed Python and Conda packages tested there. Metadata does not impose
an upper Python bound. Future Python releases require validation before support
is claimed. No additional operating-system or browser-version matrix is promised.
Dependency updates must preserve these supported environments and behaviors.
Changes that restore behavior to this contract are fixes. An incompatible change
to documented accepted behavior requires a major release even when called a bug fix.

## Deprecation and migration

Before removing or incompatibly changing a CLI option or catalog schema, announce
the replacement in at least one stable minor release of the current major series.
Keep the old behavior through the rest of that major series. Deprecation warnings
go to stderr and do not change successful exit statuses or report data.

The announcement and major-release notes must identify the last compatible
version, the first incompatible version, affected inputs, and migration steps.
Include before-and-after invocations or catalog examples and explain changes in
results. For a new catalog schema, document field mappings and information that
requires manual curation. Do not silently rewrite catalogs or invent missing
annotations. Users can retain the old application and pinned catalog during
migration. This policy does not promise indefinite maintenance of an older major.

## Automated protection and release review

Existing checks exercise the stable boundaries without treating internal Python
helpers or HTML markup as public APIs. Preserve their behavioral assertions when
refactoring tests. Update expected results only with a documented contract review.

| Promise | Automated evidence |
| --- | --- |
| CLI spellings, streams, statuses, gzip, both output aliases, refusal and cleanup | `tests/test_cli.py`, including `test_exact_defining_junction_is_reported`, `test_output_option_writes_the_stdout_report_byte_for_byte`, and closed-pipe tests. |
| Schema 3, annotations, optional exon extents, ordering, rejected fields | `tests/test_catalog.py`, `tests/test_annotations.py`. |
| Exact matching, unweighted counts, shared references, zero support, STAR row order | `tests/test_quantify.py`. |
| Main/context headers and ordered bytes through the installed CLI | `tests/test_cli.py`, `tests/test_star_interop.py`, and independent `tests/fixtures/star/expected-main.tsv` / `expected-context.tsv`. |
| Sample selection, alignment option dependencies, BAM/CRAM equivalence, explicit offline reference access | `tests/test_sample_name.py`, `tests/test_alignment_cli.py`, `tests/test_alignment_evidence.py`, `tests/test_cram_evidence.py`. |
| Provenance schema 1, deterministic bytes, exact fingerprints, warning codes, input checks, output cleanup | `tests/test_provenance.py`, `tests/test_star_interop.py`, `tests/installed_alignment_smoke.py`. |
| HTML counts, navigation, filters, inspection, offline and JavaScript-disabled behavior | `tests/test_html_report.py`, read/pairing tests, and `tests/paired_view_browser.py` with its Chromium runners. |
| Runtime and installed artifact compatibility | Python 3.11–3.14 CI, `tests/check_sdist.sh`, wheel/browser checks, Conda build/install checks, and `tests/test_release.py`. |

The [release review](DEVELOPMENT.md#compatibility-review) requires a contract diff,
version classification, catalog notes, migration instructions when applicable,
and passing tests against the candidate artifacts. Automated tests support the
review. They do not replace review of undocumented compatibility gaps.
