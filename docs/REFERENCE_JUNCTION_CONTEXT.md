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
intron_end = 116771656
strand = "+"

[[variants.reference_junctions]]
role = "same_acceptor"
chromosome = "chr7"
intron_start = 116774879
intron_end = 116774880
strand = "+"
```

The coordinates above illustrate the schema only; they are deliberately not a
curated MET definition. Curation must pin a transcript accession and annotation
release, derive coordinates from its adjacent exon boundaries, and record its
source before the reference catalog moves to version 2.

Reference junctions may be shared by multiple variants. A reference junction
may not duplicate its own variant's defining junction. STAR rows matching any
catalogue target must remain unique, because `SJ.out.tab` is already collapsed.

## First implementation slice

`load_catalog` returns ordered reference-junction definitions and `quantify`
returns their STAR support with each `VariantSupport`. The stable command-line
TSV is intentionally unchanged in this slice. Before exposing the values, the
report design needs a decision between wide role-specific columns and a
separate long-form context table. A long-form table is the more extensible
choice because variants may have differing context roles.
