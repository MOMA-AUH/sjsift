# Reference-junction context

This design extends sjsift with local, catalogued reference-junction evidence.
It does not alter defining-junction matching, make a clinical call, or infer a
whole-transcript isoform fraction.

## Why the context is explicit

The useful comparator for a skipped-exon junction usually shares one splice
site with it. For example, a `13→15` variant junction can be interpreted with
the canonical `13→14` junction (same donor) and `14→15` junction (same
acceptor). Those two comparators answer different questions and must remain
separate; adding them would double-count an opportunity to splice.

Alternative terminal exons need an event-specific comparator. `ARv7`, for
example, can compare `3→CE3` with `3→4`; `FGFR2-E18-C3` can compare its
alternative terminal exon junction with the canonical terminal-exon junction.

## Catalog schema version 2

Version 1 catalogs remain supported and have no reference context. A version 2
variant has a required `reference_junctions` array. Each entry has a unique
`role` within that variant and uses the same STAR-native coordinate convention
as a defining junction.

```toml
schema_version = 2
genome_assembly = "GRCh38"

[[variants]]
id = "METex14"
chromosome = "chr7"
intron_start = 116771655
intron_end = 116774880
strand = "+"

[[variants.reference_junctions]]
role = "same_donor"
chromosome = "chr7"
intron_start = 116771655
intron_end = 116771848
strand = "+"

[[variants.reference_junctions]]
role = "same_acceptor"
chromosome = "chr7"
intron_start = 116771990
intron_end = 116774880
strand = "+"
```

These are the curated MET exon 13→14 and 14→15 reference junctions from
RefSeq transcript NM_000245.4 (GENCODE v49 coordinate evidence). The repository's
[`definitions/grch38.toml`](../definitions/grch38.toml) uses schema version 2
and provides 32 reference entries across all 17 variants. See the
[curation record](REFERENCE_JUNCTION_CURATION.md) for pinned transcripts,
annotation checksum, exon evidence, and every derived coordinate.

`same_donor` and `same_acceptor` describe sites in transcript orientation.
On the plus strand they share `intron_start` and `intron_end`, respectively;
on the minus strand they share `intron_end` and `intron_start`. These are
curated roles, not constraints imposed on arbitrary user-defined role names.

Reference junctions may be shared by multiple variants. A reference junction
may not duplicate its own variant's defining junction. STAR rows matching any
catalogue target must remain unique, because `SJ.out.tab` is already collapsed.

## Command-line output

`load_catalog` returns ordered reference-junction definitions and `quantify`
returns their STAR support with each `VariantSupport`. The existing main TSV is
unchanged. Pass `--context-output PATH` to write an additional long-form TSV:

```text
variant_id	genome_assembly	context_role	chromosome	intron_start	intron_end	strand	unique_support	multimapping_support	total_support
METex14	GRCh38	same_donor	chr7	…	…	+	34	2	36
METex14	GRCh38	same_acceptor	chr7	…	…	+	29	1	30
```

The file has a header even when the catalog contains no reference junctions.
Rows are ordered first by catalog variant and then by the variant's declared
reference-junction order. A long-form table is used because variants may have
different context roles.

For a readable overview and per-variant evidence view, add `--html-output PATH`.
The [HTML report](HTML_REPORT.md) displays the same counts alongside the defining
junction and keeps the main and context TSV formats unchanged.
