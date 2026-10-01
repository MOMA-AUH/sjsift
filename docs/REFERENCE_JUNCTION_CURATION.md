# GRCh38 reference-junction curation

Curation date: 2026-09-30.

The references in `definitions/grch38.toml` provide local evidence for the
alternative splice sites used by each defining junction. The role describes the
shared site in transcript orientation. On the plus strand, `same_donor` shares
`intron_start` and `same_acceptor` shares `intron_end`. On the minus strand,
`same_donor` shares `intron_end` and `same_acceptor` shares `intron_start`.

For a junction skipping internal exons, the two selected references are the
adjacent annotated junction leaving the same donor exon and the adjacent
annotated junction entering the same acceptor exon in the selected reference
transcript. They remain separate observations: their sum is not a count of
full-length transcripts.

## Alternative terminal exons

### ARv7

AR-V7 joins canonical exon 3 to the cryptic terminal exon CE3. Hu et al.
identified the cryptic exons and experimentally characterized AR-V7; the study's
Figure 1 uses hg18 genomic coordinates, so those positions must not be copied
into a GRCh38 catalog. The appropriate local comparator is canonical exon
`3→4`, assigned `same_donor`. This comparator choice follows from the shared
exon-3 donor; it is a curation decision, not an assay threshold from the paper.
[Hu et al., Cancer Research (2009), Figure 1](https://pubmed.ncbi.nlm.nih.gov/19117982/).

Illumina's versioned DRAGEN 4.4 documentation independently identifies
`chrX:67686127-67694672` on the plus strand as the hg38 ARv7 junction, using the
same first/last intron-base convention as the catalog.
[DRAGEN 4.4, Splice Variant Caller, Knowns List](https://help.connected.illumina.com/dragen/dragen-v4.4/product-guide/dragen-v4.4/dragen-rna-pipeline/splice-variant-caller#knowns-list).

There is no second canonical reference entering CE3 in the full-length
transcript: CE3 is the alternative terminal exon. Accordingly, the catalog uses
one reference for this event.

### FGFR2-E18-C3

Zingg et al. describe E18-C3 as an alternative terminal exon within FGFR2 intron
17. Their full-length terminal exon is E18-C1, and Extended Data Figure 14c
compares junctions from E17 to E18-C1, alternative E18-C2/C3/C4, and other
partners. The selected comparator is therefore `E17→E18-C1`, assigned
`same_donor`. These exon names follow the paper's biological nomenclature;
annotation-specific exon ranks must be recorded separately if they differ.
[Zingg et al., Nature (2022), Expression of FGFR2ΔE18 in human cancer and Extended Data Figure 14](https://pmc.ncbi.nlm.nih.gov/articles/PMC9436779/).

FGFR2 is on the minus strand. For the catalog's defining junction
`chr10:121482178-121483697`, the shared donor is therefore the higher genomic
endpoint, `121483697`. The reference must join that donor to the canonical
terminal-exon acceptor. A junction sharing only the lower endpoint would be a
`same_acceptor` reference and would not implement the intended comparison.

This curation selects one canonical terminal-exon comparator. The paper also
reports C2/C4 and rearranged partners, so the variant and its selected reference
do not exhaust all junctions leaving E17. They cannot by themselves establish a
whole-transcript isoform fraction.

## EGFR event labels

The established EGFRvIII event removes exons 2–7, and EGFRvII removes exons
14–15. These mappings are described in primary tumor sequencing work.
[Francis et al., Cancer Discovery (2014)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4125473/).

The exact names `EGFRvIIIb` and `EGFRvIIb` did not yield a supporting primary
publication in this curation's searches. They are retained as existing project
identifiers. Their comparator choices must be derived from the actual catalog
endpoints and the pinned transcript's exon boundaries, rather than inferred
from similarity to the EGFRvIII or EGFRvII name. An endpoint mapping establishes
the local splice comparison; it does not establish that the project suffix is a
standard literature name or that the junction has the same biological effect as
the similarly named event.

## Pinned annotation and transcript choices

Coordinates were derived from the chromosome-only comprehensive GTF for
[GENCODE human release 49 (GRCh38.p14)](https://www.gencodegenes.org/human/release_49.html),
whose header identifies Ensembl 115 and annotation date 2025-07-08. The exact
source is
[`gencode.v49.annotation.gtf.gz`](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz),
retrieved on 2026-09-30. Its compressed-file SHA-256 is:

```text
d6e6fe0515c95b2a8cd36a853c1989cee9115c736c60237c56ae92b9daaaf7c4
```

The five baseline transcripts carry the `MANE_Select` tag in that release.
Two additional EGFR transcripts are necessary to match the existing catalog's
alternative endpoints. Transcript selection here is a documented local
comparison, not a claim that only this transcript is expressed.

| Gene / entries | Reference transcript | GENCODE coordinate source | Selection in GENCODE v49 |
|---|---|---|---|
| AR | NM_000044.6 | ENST00000374690.9 | MANE Select |
| BRAF | NM_004333.6 | ENST00000646891.2 | MANE Select |
| EGFR, except the two entries below | NM_005228.5 | ENST00000275493.7 | MANE Select |
| EGFRvIIIb | NM_001346900.2 | ENST00000450046.2 | Alternative first exon; `basic`, `CCDS`, `RNA_Seq_supported_only` |
| EGFRvIIb | NM_201284.2 | ENST00000344576.7 | Alternative terminal exon; `basic`, `CCDS` |
| FGFR2 | NM_000141.5 | ENST00000358487.10 | MANE Select |
| MET | NM_000245.4 | ENST00000397752.8 | MANE Select |

The catalog labels use versioned RefSeq `NM_` accessions. The GENCODE IDs remain
in the source evidence so the coordinate extraction is reproducible. The five
MANE pairs were verified against NCBI's
[MANE GRCh38 v1.4 summary](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz).
[MANE pairs have identical sequence and exon structure](https://www.ncbi.nlm.nih.gov/refseq/MANE/).
All seven mappings also appear in the pinned
[GENCODE v49 RefSeq metadata](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.metadata.RefSeq.gz).
Both files were retrieved on 2026-09-30; their compressed-file SHA-256 values are:

```text
MANE.GRCh38.v1.4.summary.txt.gz
4b3992457556e302a5e47e18a305bb4763718377696ca37d5a6b34df7db630d2
gencode.v49.metadata.RefSeq.gz
fee9360f88af80833fa957775b1aab968b57d9a5fd662fcec5fab87ba390d9ba
```

For the two non-MANE EGFR transcripts, the exon alignments were checked directly
against the [UCSC hg38 NCBI RefSeq track](https://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeq;chrom=chr7;start=55000000;end=55220000)
on 2026-09-30. NM_201284.2 has the same 16 exon intervals as ENST00000344576.7.
NM_001346900.2 has the same 28-exon splice chain as ENST00000450046.2, but its
outer transcript boundaries differ: RefSeq exon 1 begins at 55109839 instead of
55109723, and its final exon ends at 55211628 instead of 55211536. All splice
boundaries used here, including the exon-1 donor at 55109958, agree. Therefore
its `NM_` name is valid for these local junctions, without implying that the two
complete transcript sequences are identical.

The small [exon evidence snapshot](curation/gencode-v49-exons.tsv) contains all
137 exon records of these seven transcripts, extracted from that GTF. It records
gene, versioned transcript and exon accessions, chromosome, strand, exon rank,
and 1-based inclusive exon coordinates. The tests use these exon boundaries to
check every curated reference without downloading annotations.

GENCODE defines exon rank from the transcript's 5′ end, including for minus-strand
transcripts. For two adjacent exons ordered by genomic coordinate, the STAR
intron interval is `lower exon end + 1` through `higher exon start - 1`.
See the [GENCODE GTF field definitions](https://www.gencodegenes.org/pages/data_format.html)
and the [STAR junction coordinate convention](SPECIFICATION.md#star-junction-file).

### EGFR alternative endpoints

The annotations map `EGFRvIIIb` to exon 1→8 of NM_001346900.2. Its alternative
first exon ends at 55109958; exon 2 starts at 55142286. The `same_donor`
reference is therefore 55109959–55142285, not the MANE transcript's exon-1
junction. Its exon 7→8 acceptor reference is shared with `EGFRvIII`.
The alternative first exon is untranslated in this annotation, with coding
sequence beginning in exon 2; this endpoint mapping does not establish a
protein effect equivalent to EGFRvIII.

The annotations map `EGFRvIIb` to exon 13→16 of NM_201284.2. Its terminal
exon 16 starts at 55170307, distinct from the MANE exon 16 at 55171175.
Within this selected transcript, the reference path includes exons 14 and 15,
so both exon 13→14 (`same_donor`) and exon 15→16 (`same_acceptor`) are
appropriate local references. These are annotation-derived mappings, not
literature-based reinterpretations of the existing identifiers.

## Curated junctions

The table below uses exon ranks in the selected transcript. Each coordinate
interval is a 1-based inclusive intron. `—` means that no acceptor-side reference
is selected for that alternative terminal-exon event. The catalog contains
32 context entries representing 22 distinct junctions across 17 variants.
All existing defining-junction coordinates and variant order are preserved.

| Variant | Chromosome / strand | Donor reference (exons; intron) | Acceptor reference (exons; intron) |
|---|---|---|---|
| ARv7 | chrX / + | 3→4; 67686127–67711401 | — |
| ARv567es | chrX / + | 4→5; 67711690–67717477 | 7→8; 67722985–67723685 |
| BRAFdel2-10 | chr7 / - | 1→2; 140850213–140924565 | 10→11; 140781694–140783020 |
| BRAFdel2-8 | chr7 / - | 1→2; 140850213–140924565 | 8→9; 140787585–140794307 |
| BRAFdel3-10 | chr7 / - | 2→3; 140834873–140850110 | 10→11; 140781694–140783020 |
| BRAFdel3-8 | chr7 / - | 2→3; 140834873–140850110 | 8→9; 140787585–140794307 |
| BRAFdel4-10 | chr7 / - | 3→4; 140808996–140834608 | 10→11; 140781694–140783020 |
| BRAFdel4-8 | chr7 / - | 3→4; 140808996–140834608 | 8→9; 140787585–140794307 |
| EGFRvIVa | chr7 / + | 24→25; 55200414–55201187 | 27→28; 55202626–55205255 |
| EGFRvIII | chr7 / + | 1→2; 55019366–55142285 | 7→8; 55154153–55155829 |
| EGFRvIIIb | chr7 / + | 1→2; 55109959–55142285 | 7→8; 55154153–55155829 |
| EGFRvIVb | chr7 / + | 24→25; 55200414–55201187 | 26→27; 55201783–55202516 |
| EGFRvII | chr7 / + | 13→14; 55161632–55163732 | 15→16; 55165438–55171174 |
| EGFRvIIb | chr7 / + | 13→14; 55161632–55163732 | 15→16; 55165438–55170306 |
| FGFR2-E18-C3 | chr10 / - | 17→18; 121480022–121483697 | — |
| METex7-8 | chr7 / + | 6→7; 116755516–116757436 | 8→9; 116757775–116758458 |
| METex14 | chr7 / + | 13→14; 116771655–116771848 | 14→15; 116771990–116774880 |

Shared references intentionally appear under each relevant variant. They are
counted once per STAR junction and reported under each configured role; the
context rows must not be summed into gene expression or transcript fractions.
The roles are descriptive catalog labels; the generic schema does not enforce
shared splice sites, so the curated-catalog regression tests check that relation
on both strands.

### Reference-sequence cross-check

All 22 distinct reference junctions were also checked against GRCh38 sequence
retrieved through UCSC's `getData/sequence` API on 2026-09-30. The two intronic
bases at each end form GT–AG in transcript orientation for every reference
(genomic CT–AC on the minus strand). This is a sequence check, not a claim about
observed support in a sample. The queried regions are linked below; API interval
starts are 0-based and ends are exclusive, unlike the catalog intron coordinates.

| Gene | hg38 sequence query |
|---|---|
| AR | [chrX:67544000–67730700](https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrX;start=67544000;end=67730700) |
| BRAF | [chr7:140730000–140925000](https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr7;start=140730000;end=140925000) |
| EGFR | [chr7:55019000–55211700](https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr7;start=55019000;end=55211700) |
| FGFR2 | [chr10:121478000–121599000](https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr10;start=121478000;end=121599000) |
| MET | [chr7:116672000–116799000](https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr7;start=116672000;end=116799000) |

## Reproducing the exon evidence

Download the pinned GTF above, verify its SHA-256, and run:

```bash
python docs/curation/extract_gencode_exons.py gencode.v49.annotation.gtf.gz > exons.tsv
cmp exons.tsv docs/curation/gencode-v49-exons.tsv
```

This is an offline curation helper, not a runtime download or automatic catalog
update. Re-curation against a different release requires reviewing the transcript
choices, exon evidence, and reference junctions together.

## Annotated schema 3 (2026-10-01)

The current catalog promotes the selected RefSeq accession, donor/acceptor
labels, and annotation provenance to required fields on every defining and
reference junction. All 17 defining and 32 reference entries retain their
coordinates, strand, identifiers, roles, and order. Labels are tested against
the pinned exon evidence above. ARv7 explicitly uses `3→CE3`; FGFR2 uses the
paper's `E17→E18-C3` and comparator `E17→E18-C1`, with annotation ranks 17/18
recorded in provenance. Neither canonical accession is asserted to contain the
alternative terminal exon; its full extent is not curated here. The two alternate
EGFR accessions remain the local annotation basis. Schematics show boundaries
only and do not invent a full transcript.
