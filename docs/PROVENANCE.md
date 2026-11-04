# Run provenance JSON

Run provenance records the inputs and available assurance for one successful run.
Add `--provenance-output` to save this information as a separate JSON file.
The option works with TSV-only runs and with context, HTML, BAM, or CRAM reports.
It does not change either TSV format or add columns to those files.

```bash
sjsift --sample-name sample-01 \
  --junctions sample.SJ.out.tab \
  --definitions definitions/grch38.toml \
  --output sample.tsv \
  --provenance-output sample.provenance.json
```

Keep input files unchanged throughout the run. sjsift reads the catalog and STAR
file again to compute SHA-256 fingerprints. A fingerprint identifies exact file
bytes, including whitespace and comments. For gzip-compressed STAR input, it
identifies the compressed bytes, including the gzip header. Compression detection
uses the file header, regardless of the filename extension.

All output files use exclusive creation. Existing files remain unchanged.
If an ordinary input or output error occurs, sjsift removes files created by that
run, including the JSON file. File writes finish before TSV output reaches stdout.
Already emitted stdout cannot be recalled. Sudden process termination and failed
filesystem cleanup are outside this guarantee.

## Schema version 1

The root is a JSON object. All fields below are present unless their description
states a condition. Names refer to basenames: filenames without directory paths.
Integers are JSON numbers with no fractional part. SHA-256 values contain 64
lowercase hexadecimal characters.

| Root field | Type | Meaning |
| --- | --- | --- |
| `schema_version` | integer | `1`, the provenance format version. |
| `sample` | object | `name` is the resolved sample name. `source` is `explicit` for `--sample-name`, or `alignment_header` for inference from one distinct header sample. |
| `software` | object | `name` is `sjsift`. `version` is the installed sjsift version string. |
| `genome_assembly` | string | The catalog's declared assembly. This does not establish the STAR file's assembly. |
| `catalog` | object | `basename`, `sha256`, `fingerprint_scope: "exact_file_bytes"`, and integer `schema_version` for the catalog format, currently `3`. |
| `star_input` | object | `basename`, `sha256`, `fingerprint_scope: "exact_file_bytes"`, and `compression`, either `none` or `gzip`. |
| `star_assurance` | object | `sample_identity` and `genome_assembly` are both `unverified`. |
| `warnings` | array | Zero or more objects with string `code` and `message` fields, in documented code order. |
| `alignment` | object | Either `{"status": "not_requested"}` or the requested-evidence object described below. |

The schema has one compatibility warning code, `NO_CATALOG_CHROMOSOMES`.
It means that no STAR row contains a catalog chromosome identifier.
An empty file also produces this warning. Its message is
`no catalog chromosome identifiers found; check genome assembly and chromosome naming compatibility`.
This is the stderr message without its input path and prefix.
HTML explains the same condition. A warning does not alter support counts or
establish a sample or assembly mismatch.

## Requested alignment evidence

`--alignment` still requires `--html-output`. A requested alignment object has
`status: "requested"`, even when no eligible alignment records exist.
Successful validation and zero eligible records differ from absent alignment input.
The remaining fields describe the selected evidence and available assurance.

| Alignment field | Type | Meaning |
| --- | --- | --- |
| `basename` | string | Alignment input basename. |
| `format` | string | `BAM` or `CRAM`, detected from the input. |
| `index_basename` | string | The discovered or explicitly selected index basename. |
| `dependencies` | object | String versions `pysam` and `htslib` from the installed alignment dependency. |
| `eligible_junction_memberships` | integer | Sum of eligible record counts across selected junction groups. Zero means no eligible evidence. A record shared by several groups contributes once to each group. This is neither STAR support nor a molecule count. |
| `sample_selection` | object | Selection fields described below. |
| `assurance` | object | Readable string fields `sample` and `assembly` describing available assurance and its limits. |
| `reference` | object or null | `null` for BAM. For CRAM, the reference information described below. |

`sample_selection.header_sample_count` is the number of distinct, nonblank header
sample names. `selected_sample_in_header` is a boolean that records whether the
resolved sample matches a header sample. `policy` has two possible values:

| Policy | Meaning |
| --- | --- |
| `all_read_groups` | Zero or one header sample. No sample-based read-group filter applies, and eligible unassigned records remain included. |
| `selected_sample_read_groups` | Multiple header samples. Only groups for the explicit selected sample contribute. Unassigned records are excluded. |

Available assembly labels on required contigs must agree with the catalog.
Missing labels leave assembly assurance unavailable. Contig names and bounds
alone do not establish sequence equivalence. Header sample agreement does not
independently establish the identity of the anonymous STAR input.

## CRAM reference information

The reference object identifies the explicit local FASTA and its indexes.
`basename` is the FASTA basename. `index_basenames` contains the `.fai` basename,
followed by the `.gzi` basename when the FASTA is compressed.
`access` is `explicit_local_only`. `assurance` is a readable string that describes
the performed checks and missing-header-checksum limits.

The `checks` object has these fields:

| Check field | Type | Meaning |
| --- | --- | --- |
| `contig_names` | boolean | `true`: every alignment header contig exists in the FASTA. |
| `contig_lengths` | boolean | `true`: every header contig length agrees with the FASTA index. |
| `full_sequences_read` | boolean | `true`: full reference sequences for all header contigs were read in bounded chunks. |
| `header_contigs` | integer | Number of contigs in the alignment sequence dictionary. |
| `header_m5_checked` | integer | Number of available header M5 checksums compared successfully with the local sequences. Missing M5 values limit sequence assurance. |
| `decoded_slice_checks` | boolean | `true`: HTSlib checksum checks remain enabled when decoding requested regions. This does not claim decoding of the whole CRAM file. |

M5 is the alignment header's MD5 sequence checksum. These checks establish only
the available CRAM/reference compatibility. They do not establish the anonymous
STAR file's sample identity or assembly.

## Retained metadata and compatibility

Schema 1 has no content fingerprints for alignment files, indexes, or reference
files. Basenames and successful checks do not uniquely identify those files.
Archive those inputs separately when exact replay requires them.

The JSON contains no read names, biological sequences, qualities, full headers,
directory paths, or copied command lines. It retains the resolved sample name,
assembly label, and input basenames. These are input metadata, not anonymized labels.

With identical input bytes, basenames, configuration, and software/dependency
versions, sjsift produces identical JSON bytes. Output destinations do not affect
the JSON. Serialization uses sorted object keys, two-space indentation, ASCII
escapes for non-ASCII characters, LF line endings, and one final newline.
No run timestamp or output path appears in the document.

Consumers must inspect `schema_version` before reading fields. Within version 1,
new fields and warning codes can be added. Consumers must tolerate unknown fields
and preserve unknown warnings. Removing fields or changing their types or meanings
requires a new schema version. Readable messages can change without a schema
change. Use warning codes and structured fields for machine decisions.
