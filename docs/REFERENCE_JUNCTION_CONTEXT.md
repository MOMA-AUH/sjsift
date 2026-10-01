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

## Catalog schema version 3

Only schema 3 is supported; older schemas are rejected. A current
variant has a required `reference_junctions` array. Each entry has a unique
`role` within that variant and uses the same STAR-native coordinate convention
as a defining junction.

```toml
schema_version = 3
genome_assembly = "GRCh38"

[[variants]]
id = "METex14"
chromosome = "chr7"
intron_start = 116771655
intron_end = 116774880
strand = "+"
reference_transcript = "NM_000245.4"
donor_exon = "13"
acceptor_exon = "15"
annotation_source = "GENCODE v49 / ENST00000397752.8; pinned exon evidence"

[[variants.reference_junctions]]
role = "same_donor"
chromosome = "chr7"
intron_start = 116771655
intron_end = 116771848
strand = "+"
reference_transcript = "NM_000245.4"
donor_exon = "13"
acceptor_exon = "14"
annotation_source = "GENCODE v49 / ENST00000397752.8; pinned exon evidence"

[[variants.reference_junctions]]
role = "same_acceptor"
chromosome = "chr7"
intron_start = 116771990
intron_end = 116774880
strand = "+"
reference_transcript = "NM_000245.4"
donor_exon = "14"
acceptor_exon = "15"
annotation_source = "GENCODE v49 / ENST00000397752.8; pinned exon evidence"
```

These are the curated MET exon 13→14 and 14→15 reference junctions from
RefSeq transcript NM_000245.4 (GENCODE v49 coordinate evidence). The repository's
[`definitions/grch38.toml`](../definitions/grch38.toml) uses schema version 3
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

## Alignment evidence and shared membership

Optional indexed alignment input supplies separate previews for every reference
role, including custom roles and references on different contigs or strands.
Each group reports its own STAR support, eligible and embedded alignment counts,
strand totals and exclusions. Sampling is independent per junction and mapping
class. A configured zero is different from no reference configuration or evidence
not requested.

The same physical alignment may support several selected junctions. It appears
in every qualifying group with the same record ID and an **Also matches** list;
read name and mate labels identify related records without collapsing them.
Shared references and overlapping regional fetches do not count a physical
record twice within one group. Separate records with identical visible fields
are retained. Shared rows and reference counts must not be summed into molecule,
template, transcript or isoform counts. See [alignment evidence](ALIGNMENT_EVIDENCE.md).
