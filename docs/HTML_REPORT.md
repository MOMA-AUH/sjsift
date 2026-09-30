# HTML report

`--html-output PATH` adds a self-contained HTML report to the existing main TSV
and optional context TSV. It consumes the same quantification results, including
zero-support entries. There are no additional runtime dependencies, remote
resources, or browser network requests.

## Creating and opening a report

Run the current source checkout with a new destination path:

```bash
sjsift --junctions sample.SJ.out.tab \
  --definitions definitions/grch38.toml \
  --output sample.tsv \
  --context-output sample.context.tsv \
  --html-output sample.html
```

Open `sample.html` directly in a browser. The file can be moved or shared on its
own. `--context-output` is optional; omitting `--output` sends the main TSV to
stdout while still creating the HTML. Both catalog schemas are supported;
schema 1 has no reference context. These report extensions are not included in
the released v0.1.3 package. See the [current specification](SPECIFICATION.md)
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

The current catalog schema does not expose its transcript comments or exon
labels to the loader. The HTML therefore displays junction coordinates, rather
than extracting NM accessions from comments or guessing exon names. The curated
catalog's [transcript provenance](REFERENCE_JUNCTION_CURATION.md) remains in its
documentation. Splice diagrams are deferred until the necessary metadata is
explicit in the schema.

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
