# HTML report

`--html-output PATH` adds a self-contained HTML report to the existing main TSV
and optional context TSV. It consumes the same quantification results, including
zero-support entries. Pysam is installed as a normal runtime dependency for optional alignment
evidence; there are no remote resources or browser network requests.

## Creating and opening a report

Run the current development version with a schema-3 catalog with a new destination path:

```bash
sjsift --sample-name sample --junctions sample.SJ.out.tab \
  --definitions definitions/grch38.toml \
  --output sample.tsv \
  --context-output sample.context.tsv \
  --html-output sample.html
```

Open `sample.html` directly in a browser. The file can be moved or shared on its
own. `--context-output` is optional; omitting `--output` sends the main TSV to
stdout while still creating the HTML. Only annotated schema 3 is supported; empty reference arrays mean no reference context. See the [specification](SPECIFICATION.md)
for the complete input, schema, and CLI contract.

## Overview and variant details

The overview contains one row per catalog variant. Defining-junction support
and each distinct reference role get separate columns. Standard `same_donor`
and `same_acceptor` roles have readable headings; custom role names are shown
literally. A role missing from a variant is labeled **Not configured**, whereas
an unmatched configured junction retains its numeric zero counts.

The default order is decreasing variant unique support, then decreasing variant
total support, with ties preserving catalog order. Users can restore catalog
order, search variant IDs, or show only variants with any unique or multimapping
defining-junction support. This filter applies no minimum support threshold.

Selecting a variant opens its detail view, which shows the defining junction
and every configured reference in declared order. Each junction has unique,
multimapping, and total counts, plus its chromosome, 1-based inclusive intron
interval, and strand. The variant selector and Previous/Next buttons follow the
current filters and order. The browser's back button and links to individual
variant details work locally.

Bars stack unique and multimapping counts without weighting. Overview bars use
one linear scale across the entire catalog, including when filtered. Detail
bars use one scale within the selected variant. Tiny bars can be visually hard
to distinguish from zero; exact numeric counts are always displayed.

## Reading the evidence

Evidence descriptions depend only on whether total support is zero or nonzero:

| Description | Meaning |
|---|---|
| Variant and reference support | The defining junction and at least one selected reference have support. |
| Variant support only | The defining junction has support; all selected references have zero reported support. |
| Reference support only | At least one selected reference has support; the defining junction has zero reported support. |
| No selected junction support | All selected junctions have zero reported support. |
| No reference context configured | This entry provides no reference junctions; this is different from observed zero reference support. |

Unique and multimapping counts remain separate. An entry supported only by
multimapping counts is included by the “any support” filter, even if its unique
count is zero. No description constitutes a positive/negative or clinical call.
Zero refers to the already-filtered STAR file, not to a recount of the original
alignments.

Compare the donor and acceptor references separately. Shared reference junctions
may appear for several variants. Summing these reference counts would not yield
gene expression or a whole-transcript isoform fraction. The report does not
calculate a single variant percentage.

An empty input or lack of chromosome-name overlap produces the same stderr
compatibility warning as before and embeds a corresponding warning in the HTML.

## Provenance and file behavior

The report opens with **SJSIFT · JUNCTION REPORT** and the sample name as a subtitle;
the browser title also includes the sample name. Supply `--sample-name NAME`, or
let sjsift infer the single distinct `@RG SM` name from a BAM/CRAM input. Missing
sample metadata requires an explicit name. Multiple names require an explicit
selection matching the header; alignment evidence then includes only that sample's
read groups. A supplied name must also match a single named alignment sample.

**Inputs and configuration** appears at the bottom as two-column key-value rows:
Alignment file, STAR junctions, Assembly, and sjsift version. The alignment file
is labeled “Not supplied” when absent.
Basenames identify the
supplied files; neither a supplied name nor alignment metadata verifies the anonymous
STAR file's sample identity. The catalog's assembly cannot independently establish the assembly
used to create the STAR file.

The file embeds all styles and scripts. With JavaScript disabled, the complete
static overview and all detail sections remain visible and linked. Search,
sorting, and the single-variant view require a modern JavaScript-enabled browser.
Labels from catalogs and filenames are escaped as text, never inserted into
executable JavaScript. Output is deterministic for identical inputs, basenames,
and software version; no generation timestamp is added.

All requested output paths must differ and must not already exist. sjsift opens
every destination exclusively before writing any report, flushes and closes file
outputs, then writes the main TSV to stdout if requested. On a write failure,
newly created report files are removed where possible; existing files are never
overwritten. Bytes already emitted to stdout cannot be retracted. This cleanup
does not provide transactional guarantees against process termination or power
loss.

## Optional alignment previews

Supply `--alignment sample.bam` to embed selected structural evidence. Read
[alignment evidence](ALIGNMENT_EVIDENCE.md) for exact inclusion rules, independent
counts, input checks and the raw identifiers retained when sharing.
For CRAM, use `--alignment sample.cram --reference genome.fa` with local existing
alignment and FASTA indexes. The explicit FASTA is always required and remote
reference retrieval is disabled. The generated report contains no reference bases
or absolute paths and opens independently of both source CRAM and FASTA.

Defining and configured reference junctions share one pileup per variant, with
colored junction guides labeled by exon connection and STAR support. Variant reads
are blue, same-donor reads teal, and same-acceptor reads purple. Custom reference
roles are gray. Reads use solid colors at full opacity regardless of mapping class
or transcript-strand evidence. Shared reads appear once, colored by their first
matching group in catalog order (variant first), with a **Supports** list.
Mapping class and strand evidence remain available in read details.
**Also matches** identifies shared physical records across
qualifying junction groups; record IDs are stable within the report and identical
input/settings. Empty reference arrays, zero eligible records and alignment
input not requested remain distinct states.

All eligible alignment records are embedded without downsampling. The initial
figure shows up to ten alignment rows per support group. **Expand reads**
pages through the embedded evidence, up to 50 rows per group at a time;
**Collapse reads** restores the preview. Shared reads belong to the first matching
group in catalog order for pagination. **All / Unique only** and **Support** change
the displayed subset, leaving STAR support unchanged. Unique only means NH=1, regardless of primary,
duplicate, supplementary or QC flags. The page status labels the last page and
filter-empty state. Paging affects visibility only; every eligible read remains
available in the report.

Reads form a compact pileup sorted by support group, then genomic start, beneath a shared coordinate
ruler. Matching positions align vertically; thin splice lines connect aligned
blocks. Gaps longer than 80 nt without aligned blocks or displayed soft clips are compressed
at shared `//` marks. The scale remains fixed when filtering or paging.
Deletions, insertions and clipping have distinct marks with exact lengths and
coordinate tooltips. Soft clips are gold blocks extending from the alignment ends
at the same width per base as aligned blocks, with space for their full lengths.
They have no assigned reference coordinates. Hard clips remain gold ticks.
References on other chromosomes have separate labeled axes
within the figure. Click a read, or focus it and press Enter, to open
**Read details**, which includes identity, flags,
strand sources and the explicitly defined genomic left/right aligned anchors;
these are not STAR maximum overhang. The alignment guide specifies the anchor
calculation. Without JavaScript, ten previews per group, native read details,
STAR support counts and labeled junction guides remain readable; filtering and paging require JavaScript.
