# sjsift

This context defines the language used to quantify predefined RNA splice variants from splice-junction evidence.

## Language

**Known splice variant**:
A predefined RNA splicing event represented by exactly one defining splice junction.
_Avoid_: Candidate variant, discovered variant

**Genome assembly**:
The named reference coordinate system to which every junction in one variant-definition catalog belongs. The label documents an assumption about the input STAR file; it does not prove that the file was generated from that assembly.
_Avoid_: Genome version, reference version

**Reference coordinate namespace**:
The combination of a genome assembly and its exact chromosome identifiers. Documentation examples use GRCh38 with `chr`-prefixed identifiers.
_Avoid_: Assembly, when chromosome naming is also intended

**Chromosome identifier**:
The exact contig name shared by a catalog junction and column 1 of its compatible STAR junction file. Identifiers such as `1` and `chr1` are distinct.
_Avoid_: Normalized chromosome name

**Junction coordinates**:
The first and last intronic bases, expressed as 1-based inclusive coordinates exactly as in columns 2 and 3 of STAR `SJ.out.tab`.
_Avoid_: Exon-boundary coordinates, BED coordinates

**Defining splice junction**:
The assembly-specific combination of chromosome identifier, junction coordinates, and transcript strand that identifies one known splice variant. sjsift supports only definitions expected to have a STAR-defined strand; a STAR junction with strand `0` does not match.
_Avoid_: Coordinate-only junction, signature region

**Junction support**:
The counts already reported by STAR for one defining or reference splice junction, retained as unique support and multimapping support. Total support is their arithmetic sum, not an independent recount, weighted value, or classification score.
_Avoid_: Variant score, filtered support, recounted support

**Unsupported variant**:
A known splice variant whose defining splice junction has no exact match in the sample's STAR junction file. It is reported with zero support and does not imply a clinical negative or that the variant was not examined.
_Avoid_: Negative variant, absent variant, untested variant

**Variant-definition catalog**:
A human-maintained collection of known splice variants whose defining splice junctions all use one declared genome assembly. Variant identifiers and defining splice junctions are each unique, preserving a one-to-one relationship.
_Avoid_: Variant database, genome library

**GRCh38 reference catalog**:
The versioned sjsift catalog containing 17 variants across AR, BRAF, EGFR, FGFR2, and MET, with curated reference junctions. It remains an explicit input rather than a hidden application default.
_Avoid_: Built-in database, comprehensive splice-variant catalog

**Reference junction**:
A catalogued splice junction whose support provides local context for a variant's defining junction. Each has a role unique within that variant, such as `same_donor` or `same_acceptor` in transcript orientation. References remain separate observations and may be shared across variants.
_Avoid_: Gene expression, full-length transcript count

## Selected variant identifiers

**EGFRvIII**:
The EGFR splice variant defined by the exon 1-to-exon 8 junction.
_Avoid_: EGFR variant III

**Supported EGFR variants**:
The six EGFR variants included in the reference catalog: `EGFRvII`, `EGFRvIIb`, `EGFRvIII`, `EGFRvIIIb`, `EGFRvIVa`, and `EGFRvIVb`.
_Avoid_: EGFR variant set

**METex14**:
The project identifier for MET exon 14 skipping, defined by the exon 13-to-exon 15 junction.
_Avoid_: MET exon 14 deletion, when describing the RNA event

**ARv7**:
The AR splice variant defined by the exon 3-to-cryptic exon 3 junction.
_Avoid_: AR variant 7

**Junction annotation**:
Explicit reference-transcript accession, donor/acceptor exon labels, and source
provenance for a local splice connection. An accession is an annotation or
comparison basis, not a claim that it contains the aberrant junction.
_Avoid_: Variant transcript accession, inferred exon label

**Local junction schematic**:
Labeled splice boundaries in transcript orientation, with exact intronic
coordinates and no claim about full exon extents or complete transcript structure.
_Avoid_: Full transcript model, exon coverage

**Alignment evidence**:
Individual selected-junction alignment records from optional indexed input.
Their eligible, embedded and shown counts are separate from STAR support.
_Avoid_: Recounted STAR support, molecule count

**Mapping class**:
Unique (valid NH=1), multimapping (valid NH>1), or unknown multiplicity. Alignment
flags and MAPQ do not establish uniqueness.
_Avoid_: Primary means unique

**Transcript-strand evidence**:
Usable transcript-strand tags or validated per-junction STAR motifs. Agreeing and
unverified records are displayed separately; opposite/conflicting evidence is
excluded with counts. Forward/reverse alignment orientation is a different fact.
_Avoid_: Read orientation proves transcript strand

**Shared alignment membership**:
One physical alignment record supporting several selected junction groups. It
has the same record ID in each and an Also matches indication. This is distinct
from related mates sharing a read name, and does not imply independent molecules.
_Avoid_: Deduplicated template, independent reference molecules

**CRAM reference FASTA**:
The explicit local reference sequence used to validate and decode CRAM
alignments, distinct from selected reference junctions used as local splice
comparators. Available sequence checks establish CRAM/FASTA compatibility,
not the sample identity or assembly of the anonymous STAR input.
_Avoid_: Reference junction file, automatically downloaded genome, verified STAR assembly
