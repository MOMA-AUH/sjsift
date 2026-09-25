# sjsift

[![Conda Version](https://img.shields.io/conda/vn/MOMA-AUH/sjsift?style=for-the-badge&cacheSeconds=300)](https://anaconda.org/MOMA-AUH/sjsift) [![Conda Downloads](https://img.shields.io/conda/dn/MOMA-AUH/sjsift?style=for-the-badge&cacheSeconds=300)](https://anaconda.org/MOMA-AUH/sjsift)

**sjsift** quantifies predefined RNA splice variants from a STAR
`SJ.out.tab` file. It matches exact, assembly-specific splice junctions from a
human-readable TOML catalog and reports STAR's unique and multimapping support
counts in a deterministic TSV table.

sjsift does not reopen alignments, discover novel splice events, or assign
biological or clinical significance. Version 0.1 supports one defining
junction per variant.

## Installation

The recommended way to install **sjsift** is via
[conda](https://docs.conda.io/), using the `MOMA-AUH` channel:

```bash
conda install MOMA-AUH::sjsift
```

sjsift requires Python 3.11 or newer and has no third-party runtime
dependencies.

## Usage

sjsift needs one STAR `SJ.out.tab` file and one variant-definition catalog. To
use the reference GRCh38 catalog from this repository:

```bash
curl -L \
  https://raw.githubusercontent.com/MOMA-AUH/sjsift/master/definitions/grch38.toml \
  -o grch38.toml

sjsift \
  --junctions sample.SJ.out.tab \
  --definitions grch38.toml \
  --output sample.sjsift.tsv
```

`--junctions` accepts plain text or gzip-compressed data, detected from the
file contents regardless of filename extension. For example:

```bash
sjsift --junctions sample.SJ.out.tab.gz --definitions grch38.toml
```

Omit `--output` to write the report to standard output:

```bash
sjsift --junctions sample.SJ.out.tab --definitions grch38.toml
```

Run `sjsift --help` for the complete command-line reference. Existing output
files are never overwritten.

## Reference catalog

[`definitions/grch38.toml`](definitions/grch38.toml) contains GRCh38
definitions for six EGFR variants (`EGFRvII`, `EGFRvIIb`, `EGFRvIII`,
`EGFRvIIIb`, `EGFRvIVa`, and `EGFRvIVb`), `METex14`, `METex7-8`, `ARv7`,
`ARv567es`, six BRAF exon-deletion junctions, and `FGFR2-E18-C3`.

Catalog coordinates use STAR's 1-based, inclusive intron convention and are
matched exactly. Chromosome names are also literal: for example, `chr7` does
not match `7`. Make sure the catalog and the STAR genome use the same assembly
and chromosome naming scheme.

You can copy the reference catalog or provide your own schema-compatible TOML
file:

```toml
schema_version = 1
genome_assembly = "GRCh38"

[[variants]]
id = "EGFRvIII"
chromosome = "chr7"
intron_start = 55019366
intron_end = 55155829
strand = "+"
```

Each variant must have a unique `id` and a unique combination of chromosome,
intron coordinates, and strand. See the [MVP
specification](docs/MVP_SPEC.md#variant-definition-catalog) for the complete
catalog schema.

## Output

The report contains one row per catalog entry, in catalog order:

```text
variant_id	genome_assembly	chromosome	intron_start	intron_end	strand	unique_support	multimapping_support	total_support
EGFRvIII	GRCh38	chr7	55019366	55155829	+	14	2	16
```

`unique_support` and `multimapping_support` are copied from STAR columns 7 and
8. `total_support` is their sum. Variants without a matching junction remain in
the report with zero counts; zero support is a quantification result, not a
clinical interpretation.

An empty STAR file, or one with no chromosome identifiers in common with the
catalog, still produces the complete zero-support report and emits a
compatibility warning on standard error. The input filename is not included in
the report.

## Scope

sjsift accepts one plain-text or gzip-compressed, nine-column STAR
`SJ.out.tab` file at a time.
It does not process BAM, CRAM, or SAM files; normalize chromosome names; lift
coordinates between assemblies; apply thresholds; or combine multiple
junctions into a call. For the full input, matching, validation, and exit-status
contract, see the [MVP specification](docs/MVP_SPEC.md).

## Development

See the [development and release guide](docs/DEVELOPMENT.md) for local setup,
testing, packaging, and release instructions. Additional design context is in
the [design decisions](docs/DECISIONS.md) and [domain
language](CONTEXT.md).

## License

MIT. See [LICENSE](LICENSE) for details.
