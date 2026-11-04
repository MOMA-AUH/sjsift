# STAR interoperability fixture

These artificial sequences test interoperability with STAR output. They do not
provide biological or clinical validation. All inputs fit in this directory.
Routine tests use the committed files and never run STAR or access the network.

## Inputs and expected results

`reference.fa` contains four 1,000-base contigs. `regenerate.py` specifies every
base with SHA-256 digests and explicit splice motifs. `multiA` and `multiB` are
identical. No external reference or annotation is used.

Each contig has exons at 101–200, 301–340, and 441–540, in 1-based inclusive
coordinates. The two introns are 201–300 and 341–440. Plus-strand introns have
GT/AG boundaries. The `minus` contig has CT/AC boundaries and reversed transcript
orientation. Exon numbers in the catalog follow transcript orientation.

The transcript concatenates the three exons. For `minus`, it is the reverse
complement of that concatenation. Read 1 contains transcript bases 61–180.
Read 2 is the reverse complement of bases 121–240. Both reads have 120 bases
and constant `I` qualities. `reads_1.fastq` and `reads_2.fastq` contain four pairs:
`plus-1`, `plus-2`, `minus-1`, and `multi-1`. The last pair maps to both duplicates.

The expected junction inventory applies to both introns on each contig:

| Contig | Strand code | Motif code | Unique support | Multimapping support |
| --- | --- | --- | --- | --- |
| plus | 1 (+) | 1 (GT/AG) | 2 | 0 |
| minus | 2 (-) | 2 (CT/AC) | 1 | 0 |
| multiA | 1 (+) | 1 (GT/AG) | 0 | 1 |
| multiB | 1 (+) | 1 (GT/AG) | 0 | 1 |

All eight STAR rows have annotation flag 0 and maximum overhang 40.
`expected-main.tsv` and `expected-context.tsv` are literal expectations derived
from this design, independent of sjsift. Catalog order deliberately differs from
STAR contig order and support order. The unsupported variant uses plus:201–440,
which neither the reads nor STAR support. Every supported variant uses 201–300
as its defining junction and 341–440 as its `adjacent` reference junction.

`expected-alignments.tsv` inventories all ten alignment records. Its positions
are 1-based. It includes flags, CIGAR, MAPQ, mate positions, template lengths,
NH, HI, jM, and jI. NH is the number of mapping loci. HI identifies a mapping.
jM records intron motifs and jI records intron boundaries.
The duplicate mapping on `multiB` is secondary. Its NH=2 records remain eligible.

Read 1 spans both junctions and has shared alignment membership. On `plus`,
the two defining records and four reference records come from two pairs.
On `minus`, the defining junction has two records and the reference has one.
Each duplicate has one defining record and two reference records.
These record counts deliberately differ from STAR support. STAR counts are
retained independently, without forcing agreement with paired alignment totals.

## Regeneration

Use Python 3.11 or newer, STAR **2.7.11b**, and pysam **0.24.1** with bundled
samtools/HTSlib **1.24**. The script refuses other STAR or conversion-tool versions.
The original generation used Python 3.14.7 and the Bioconda STAR 2.7.11b executable.
Install the pinned tools in a separate environment if they are unavailable.
No STAR index is committed because the script builds it from the tiny FASTA.

From the repository root, run:

```bash
python tests/fixtures/star/regenerate.py /tmp/star-regenerated --star /path/to/STAR
python -m pytest -q tests/test_star_interop.py
```

The output directory must not exist. The script contains both complete STAR
commands, including index configuration and every nondefault alignment parameter.
It uses one thread, seed 777, no mismatches, and explicit NH/HI/jM/jI attributes.
The small-junction filters retain singly supported and multimapping junctions.
All other parameters use the pinned release defaults. See the
[STAR 2.7.11b parameter definitions](https://github.com/alexdobin/STAR/blob/2.7.11b/source/parametersDefault).

STAR writes SAM and `SJ.out.tab`. The script copies `SJ.out.tab` unchanged,
then uses `pysam.sort --no-PG` and `pysam.index` to make the committed indexed BAM.
`pysam.faidx` indexes the local FASTA. `Log.final.out` retains STAR's summary.
The BAM retains STAR's original program header and complete records.

Regeneration compares input files and STAR junction rows byte for byte, then
compares complete SAM records with the committed BAM. BAM headers contain
invocation paths, and STAR logs contain times, so their bytes can differ.
The expected report files and alignment inventory are maintained independently.
Do not replace them with sjsift output during regeneration.

The installed CLI test compares complete main/context TSV files and static HTML
counts, memberships, strand evidence, and offline resource declarations.
It also converts the BAM to CRAM with `pysam.AlignmentFile` in `wc` mode and the
local FASTA, then indexes the CRAM. The same inventory and report expectations
apply to both formats. CRAM files remain temporary. The source archive, wheel,
and Conda checks exercise these fixtures without installing STAR.
