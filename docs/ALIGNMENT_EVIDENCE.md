# Alignment evidence in offline reports

sjsift v0.2.1 accepts a local coordinate-sorted indexed BAM or CRAM
alongside the required STAR file and schema-3 catalog:

```bash
sjsift --junctions sample.SJ.out.tab --definitions catalog.toml \
  --alignment sample.bam --html-output sample.html --output sample.tsv
```

`--alignment-index PATH` overrides discovery of `sample.bam.bai`, `sample.bai`,
`sample.bam.csi`, or `sample.csi`. All eligible records are embedded without
downsampling; `--alignment-limit` is no longer supported. Alignment options require
`--alignment`, which requires `--html-output`. Indexes are never created by
sjsift. Defining and reference junctions share one pileup figure per variant.
Custom roles, references on other contigs or strands, and shared references use
the same inclusion rules. Different chromosomes have separate axes within
the figure. Previews include up to ten records per support group; expansion shows
up to 50 per group per page.

## CRAM and its explicit local FASTA

For CRAM, add `--reference genome.fa`. This is mandatory even when the CRAM
embeds reference bases or uses reference-free compression. `--reference` applies
only to CRAM. Index discovery checks `sample.cram.crai` then `sample.crai`, or
uses the explicit `--alignment-index`. The local FASTA needs an existing `.fai`;
BGZF-compressed FASTA also needs `.gzi`. Plain gzip is not a random-access FASTA
format. No alignment or reference index is built automatically.

All alignment-header contigs must exist at exactly the declared length in the
FASTA. Their sequence is read in bounded chunks and checked against any available
header M5 checksum before regional extraction; this can require reading the full
reference once. Missing M5 checksums limit sequence assurance.
HTSlib's CRAM slice checks remain enabled. Pysam's bundled header/EOF check
rejects a missing CRAM end marker, which regional decoding alone can merely warn
about; this does not require an external samtools executable. A missing or invalid reference/index,
detectable mismatch, or decoding failure aborts before any report output. Input
files must remain unchanged during the invocation.

The explicit FASTA is passed directly to HTSlib. For the single CLI process,
inherited `REF_PATH` and `REF_CACHE` are temporarily replaced with an isolated
empty local directory, then restored even on failure. Supported HTSlib versions
also reject remote header UR fallbacks. This combination prevents remote
reference retrieval; a wrong FASTA is never repaired by downloading another one.
Tests include a real HTTP listener with a successful unguarded-HTSlib positive
control, followed by zero requests from sjsift with valid, wrong or absent input.

BAM and CRAM use the same record-classification, count, retention and rendering
path. Reference validation establishes available CRAM/FASTA compatibility only;
the anonymous STAR input still cannot establish sample identity or assembly.
Internal extraction provenance records FASTA/index basenames and checks performed;
the compact HTML displays the alignment basename. Reference
sequence, read bases/qualities and absolute paths are never embedded. The HTML
remains usable after both CRAM and FASTA have been removed.

## What qualifies

An alignment must have the exact chromosome identifier and a CIGAR `N` matching
the selected intron's first and last intronic bases. With a zero-based reference
cursor `r` and skip length `L`, the STAR interval is `r+1` through `r+L`.
`M`, `=`, `X`, `D`, and `N` advance that cursor; insertions and clipping do not.
A deletion or the space between paired mates is not a splice. There is no MAPQ
or aligned-anchor threshold. Only selected-junction records are retained.

Transcript strand and forward/reverse alignment orientation are distinct:

- `TS:A:+/-` is the standard transcript-strand tag. `XS:A:+/-` follows the
  established RNA-aligner convention; numeric XS scores and other types are
  unusable as strand evidence.
- STAR `jM` is usable only with STAR program provenance in `@PG`, an integer
  array with one valid motif per CIGAR N, and (if supplied) `jI` exactly matching
  the CIGAR introns in order. Motifs 1/3/5 imply plus and 2/4/6 minus; subtract
  the annotation offset of 20 first. Motifs 0/20 give no strand evidence.
  Each selected junction uses its own motif, not a record-wide majority.
- Minimap2 `ts:A:+/-` is interpreted relative to alignment orientation only
  when minimap2 program provenance is available. It is distinct from `TS`.
- A record's `PG` selects its program/ancestor chain. Without a record `PG`,
  local tag semantics require the relevant program in every possible terminal
  header chain. Missing or ambiguous provenance leaves that source unusable.
- All usable sources must agree with the catalog strand. Opposite and mutually
  conflicting evidence have separate exclusion counts; excluded records are
  not embedded. Missing or unusable sources produce **strand unverified** rows,
  identified in accessible read labels. The compact read details omit strand-source diagnostics.

NH alone determines mapping class: an integer 1 is **unique**, an integer above
1 is **multimapping**, and absent, nonpositive, or mistyped values are **unknown**.
Primary status and MAPQ cannot establish uniqueness. Qualifying secondary,
supplementary, duplicate-marked and QC-failure records are included and labeled.

## Counts, identity, and read retention

One row is one physical alignment record. Both mates may appear separately;
the read name and numeric flags help recognize related records.
No unrelated mate is fetched. Separate records with identical fields remain
separate. A physical record retrieved at several regional queries is processed
only in the first requested region it overlaps; no name-based deduplication is
performed. A record matching several selected junctions retains membership in each
qualifying group but appears only once per variant figure, identified by its stable
record ID in the embedded HTML. Highlighted splice lines, the support filter,
and accessible read labels retain its local memberships.
Membership indicates exact supporting groups. Excluded opposite/conflicting groups
are not listed as supporting membership. Every eligible defining and reference
record is retained, so abundant references cannot displace rare defining reads.

Every eligible record is counted and embedded. There is no per-junction or
per-mapping-class cap. Extraction memory and HTML size grow with the number of
eligible records.
The figure displays each physical alignment only once, even when it has multiple
memberships. Its preview displays up to ten records per support
group, sorted by support group, genomic start and end. A shared record is assigned
to its first matching group in catalog order (variant first) for paging and block
color; its highlighted splice lines and support filter still reflect all local memberships.
**All** is the initial filter. **Unique only** shows exactly NH=1 records,
including eligible secondary, supplementary, duplicate and QC-flagged records.
It excludes unknown multiplicity. The **Support** filter selects any matching
membership, including shared reads. Changing either filter resets to the first page;
collapsing restores the preview. Paging and filtering change only which records
are visible. The last page is labeled and its Next button disabled.

## Structural geometry and anchors

Reads form a compact pileup beneath a shared genomic coordinate ruler and colored
sashimi-style junction arcs. Arc tooltips identify support types,
donor/acceptor exons, and STAR support.
Arcs share a baseline and the same coordinate transform as the reads. Longer
connections arch over shorter ones. Transcript strand remains in coordinate tooltips.
Line thickness uses total STAR support on a logarithmic scale, shared across
chromosomes within the variant. Zero-support junctions remain visible as dashed
arcs. Exact unique and multimapping counts appear in the preceding support table and do not change
with read filters. Arcs mark splice boundaries, not complete exon extents;
no coverage track is inferred from the selected alignment records. Blue identifies variant support,
teal same-donor support, and purple same-acceptor support. Custom reference roles
receive distinct colors, consistently shared by the tables, arcs, and alignment
tracks throughout the report. All records use solid block colors at full opacity. Shared records use the
color of their first matching group in catalog order (variant first).
Each row displays the complete alignment CIGAR in
increasing genomic coordinates, independent of alignment orientation and
transcript strand. Matching coordinates line up across all rows. M, = and X are
aligned blocks; thin N lines are splice gaps; red D lines are deletions; purple I
marks are insertions; gold S blocks show soft clipping at actual length, gold H
ticks mark hard clipping, and gray P is padding.
Mapping class and transcript-strand status are included in accessible read labels;
the compact read details show NH without strand-source diagnostics.
Colored splice lines highlight selected junctions. Hover over an
operation for its exact length and 1-based inclusive coordinates or interbase
position. Soft clips extend left from the alignment start or right from its end,
using the same width per base as aligned blocks. These visual extensions do not
assign reference coordinates to clipped bases or change the aligned coordinates.
Extra space at the figure edges accommodates the longest clips; the genomic ruler
labels only the reference span. Hard clips and insertions remain fixed-width
boundary markers.

The scale is nonuniform. For each annotated exon, the plot retains 120 nt of context
around aligned blocks, displayed soft clips, and selected splice boundaries,
clipped to the full annotated extent. Short outer remainders of at most 80 nt
remain visible; longer uncovered tails are trimmed. Gaps longer than 80 nt outside
that context, aligned blocks, and displayed soft clips are compressed, including
uncovered sections within long exons. Compressed gaps have straight connections
with exact intervals in tooltips. A zig-zag in an intron connection indicates
omitted numbered exons between displayed neighbors from the same transcript and
strand, using known extents to place the mark. Its tooltip lists the omitted exon
numbers. Nonnumeric labels such as CE3 do not establish omitted exon numbers.
Neutral exon boxes show the retained segments on the
aligned-base scale. Zig-zag edges mark trimmed ends and both sides of compressed
exon sections. Aligned bases are never compressed or omitted.
Labels follow genomic orientation on either strand. Unknown extents use vertical
ticks at the known splice boundary, with labels above the intron line. Shared exon
extents and shared unknown boundaries appear once, retaining all annotation labels,
transcript accessions, and sources in tooltips. Known extents include inclusive
coordinates and length, plus the displayed segment when truncated. Read spans
determine the view, never the annotated exon extent.
All other regions
share one linear scale. Selected splice boundaries are preserved. The transform
uses all embedded records across the variant's groups on that chromosome and stays
fixed when filtering or paging. No nucleotide sequence or mismatch bases are
shown. Click a read (or focus it and press Enter) to disclose its name and
metadata. With JavaScript, opening a read closes the previous inspection in
that figure. Native read inspection also works without JavaScript.

Read details show the read name; one genomic interval · CIGAR · MAPQ · NH line;
a **Read** line with mate number (or unspecified), forward/reverse orientation;
flags with their numeric value; and one compact line with both **genomic left/right aligned anchors**.
Anchors refer to the extraction junction: for a shared read, this is the first
group in catalog order that the record supports.
Each anchor is the sum of consecutive M, = and X operations immediately abutting
that N. Any I, D, N, S, H or P terminates the sum; a nonaligned immediate neighbor
gives zero. Left/right do not swap on the minus strand. For example,
`10M2I5=2D3X100N7M2I10M` has anchors 3 and 7 at its skip. These measurements
are not STAR maximum overhang, inferred exon extents or quality thresholds.

**STAR support remains unchanged** in both TSVs and HTML. STAR's upstream
junction filtering, paired-template semantics, alignment output choices and the
explicit display policy can all produce different alignment-record totals. A
difference alone is not an input error. Do not add reference opportunities or
interpret these records as independent molecules or full-length transcripts.
Evidence not requested is distinct from an observed zero eligible count.

In v0.2.1, alignment records are not reconciled to the reads contributing to
STAR's reported support. Both mates can appear separately, and alignment output
can include records that did not contribute to STAR's junction counts. A report
can therefore contain multimapping alignments beside zero STAR multimapping
support, or more unique alignment records than STAR unique support. The TSV and
support tables retain STAR's original values. Reconciliation is tracked separately
in [issue #40](https://github.com/MOMA-AUH/sjsift/issues/40).

## Validation and portability

A usable local index, coordinate sort declaration, required literal contigs and
coordinate bounds are checked. The report sample name is inferred from a single
distinct nonblank `@RG SM` value unless `--sample-name NAME` is supplied. Explicit
names must match a header sample when sample metadata exists. No usable `SM` values
requires an explicit name. Multiple distinct names require an explicit matching
selection: only records with an `RG` belonging to that sample contribute counts
and embedded evidence. Records without a sample-assigned read group are excluded
in this multi-sample case. For single-sample or unnamed files, records without
read groups remain eligible as before. Unknown or invalid read-group tags remain
errors. Available `@SQ AS` labels must equal the catalog
assembly. Header labels do not establish sequence equivalence, and anonymous
STAR data cannot prove sample identity or genome assembly. The caller must
supply compatible inputs. Regional decoding errors stop the whole invocation.
All evidence extraction finishes before any output is created. Existing output
refusal, cleanup of newly created files on ordinary errors and stderr/exit-status
conventions are unchanged.

The single HTML contains structural operations, coordinates, CIGAR, MAPQ, NH,
transcript-strand status labels, flags, read names, the selected sample name,
junction exon labels and annotation tooltips, the software version, and STAR and
alignment input basenames. Detailed extraction diagnostics and reference-validation
metadata are not displayed in the compact report. It contains no read
sequence, base qualities, reference sequence, original alignment file, or
absolute source paths. Consider the retained raw identifiers before sharing.
All text is escaped and all assets are inline; no browser network access or
source BAM/CRAM or FASTA is needed. Counts and provenance remain readable without JavaScript.
The document's content security policy also blocks external resources and connections.
The geometry is schematic; nucleotide-level inspection is outside this release.

## Dependency and validation basis

Pysam is a normal runtime dependency (`>=0.24.0,<0.25`) in both Python metadata
and the noarch conda recipe. There is no reads extra or external samtools
requirement. The 0.24 series provides Python 3.11–3.14 wheels; CI exercises the
minimum 0.24.0 on Python 3.11 and current releases on the other supported Python
versions. CI builds a wheel with dependency resolution in a clean environment
and builds/tests the conda recipe using conda-forge and bioconda. Installed checks
generate reports from real synthetic indexed alignments. These are software
fixtures, not clinical validation data.

Primary format/API references checked 2026-10-01:

- [SAM optional tags](https://samtools.github.io/hts-specs/SAMtags.pdf).
- [STAR output attributes](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault),
  and [BAM CIGAR/motif emission](https://github.com/alexdobin/STAR/blob/master/source/ReadAlign_alignBAM.cpp).
- [Minimap2 tag and splice orientation implementation](https://github.com/lh3/minimap2/blob/master/format.c).
- [Pysam 0.24 releases](https://pysam.readthedocs.io/en/latest/release.html)
  and [alignment APIs](https://pysam.readthedocs.io/en/latest/api.html).
- [HTSlib reference resolution and validation](https://www.htslib.org/doc/reference_seqs.html)
  and [HTSlib 1.24 reference fallback implementation](https://github.com/samtools/htslib/blob/1.24/cram/cram_io.c).
- [Bundled samtools header/EOF check semantics](https://www.htslib.org/doc/samtools-quickcheck.html).
