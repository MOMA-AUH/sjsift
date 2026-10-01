# HTML report

`--html-output PATH` adds a self-contained HTML report to the existing main TSV
and optional context TSV. It consumes the same quantification results, including
zero-support entries. Pysam is installed as a normal runtime dependency for optional alignment
evidence; there are no remote resources or browser network requests.

## Creating and opening a report

Run the current development version with a schema-3 catalog with a new destination path:

```bash
sjsift --junctions sample.SJ.out.tab \
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

The header records the input and catalog basenames, the catalog assembly label,
and the sjsift version. Basenames identify the supplied files; no sample identity
is inferred. The catalog's assembly cannot independently establish the assembly
used to create the STAR file.

Each junction includes a local schematic with donor and acceptor exon labels,
exact intronic boundary coordinates, transcript strand, and intron length.
The schematic reads in transcript direction (decreasing coordinates on minus
strand), is not to scale, and does not imply full exon extents or a complete
transcript model. Dashed outer edges mark the partial exon views.

Reference accessions identify the annotation or comparison basis, not necessarily
a transcript containing the aberrant junction. Each group's annotation provenance
is visible, including explicit ARv7 CE3 and FGFR2 alternative terminal-exon caveats.
All labels come from required schema fields, never guesses from comments or IDs.

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

Supply `--alignments sample.bam` to embed selected structural evidence. Read
[alignment evidence](ALIGNMENT_EVIDENCE.md) for exact inclusion rules, independent
counts, sampling, input checks and the raw identifiers retained when sharing.

Every configured reference has an immediate read preview and retains its own
annotation and counts. **Also matches** identifies shared physical records across
qualifying junction groups; record IDs are stable within the report and identical
input/settings. Empty reference arrays, zero eligible records and alignment
input not requested remain distinct states.
