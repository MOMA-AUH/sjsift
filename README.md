# sjsift

`sjsift` is a deliberately small command-line tool for quantifying predefined RNA splice variants from one STAR `SJ.out.tab` file. It matches exact, assembly-specific splice junctions from a human-editable TOML catalog and reports STAR's unique and multimapping support counts without reopening alignments or applying a biological classification.

Version 0.1 is limited to one defining junction per known variant. It does not process BAM or CRAM files, discover novel events, normalize chromosome names, score clinical significance, or filter alignments.

## Installation

The intended release channel is the `MOMA-AUH` Anaconda channel:

```bash
conda install -c MOMA-AUH sjsift
```

This command describes the v0.1 distribution target; it will become usable after the first package is published.

## Example

The repository includes [`definitions/grch38.toml`](definitions/grch38.toml), a versioned GRCh38 catalog containing six known EGFR variants plus `METx14del` and `ARv7`.

```bash
sjsift \
  --junctions sample.SJ.out.tab \
  --definitions definitions/grch38.toml \
  --output sample.sjsift.tsv
```

Representative output:

```text
variant_id	genome_assembly	chromosome	intron_start	intron_end	strand	unique_support	multimapping_support	total_support
EGFRvIII	GRCh38	chr7	55019366	55155829	+	14	2	16
METx14del	GRCh38	chr7	116771655	116774880	+	0	0	0
ARv7	GRCh38	chrX	67686127	67694672	+	3	0	3
```

Unsupported variants remain in the output with zero counts. A zero is a junction-quantification result, not a clinical interpretation.
`METx14del` names the MET exon 14-skipping RNA splice event; it does not assert a deletion in genomic DNA.

## Documentation

- [MVP specification](docs/MVP_SPEC.md)
- [Development and release workflow](docs/DEVELOPMENT.md)
- [Design decisions](docs/DECISIONS.md)
- [Implementation plan](docs/IMPLEMENTATION_PLAN.md)
- [Domain language](CONTEXT.md)
