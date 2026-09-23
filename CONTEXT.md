# sjsift

This context defines the language used to quantify predefined RNA splice variants from splice-junction evidence.

## Language

**Known splice variant**:
A predefined RNA splicing event represented in version 0.1 by exactly one defining splice junction.
_Avoid_: Candidate variant, discovered variant

**Genome assembly**:
The named reference coordinate system to which every junction in one variant-definition catalog belongs. The label documents an assumption about the input STAR file; it does not prove that the file was generated from that assembly.
_Avoid_: Genome version, reference version

**Reference coordinate namespace**:
The combination of a genome assembly and its exact chromosome identifiers. Version 0.1 documentation examples use GRCh38 with `chr`-prefixed identifiers.
_Avoid_: Assembly, when chromosome naming is also intended

**Chromosome identifier**:
The exact contig name shared by a catalog junction and column 1 of its compatible STAR junction file. Identifiers such as `1` and `chr1` are distinct.
_Avoid_: Normalized chromosome name

**Junction coordinates**:
The first and last intronic bases, expressed as 1-based inclusive coordinates exactly as in columns 2 and 3 of STAR `SJ.out.tab`.
_Avoid_: Exon-boundary coordinates, BED coordinates

**Defining splice junction**:
The assembly-specific combination of chromosome identifier, junction coordinates, and transcript strand that identifies one known splice variant. Version 0.1 supports only definitions expected to have a STAR-defined strand; a STAR junction with strand `0` does not match.
_Avoid_: Coordinate-only junction, signature region

**Junction support**:
The counts already reported by STAR for one defining splice junction, retained as unique support and multimapping support. Total support is their arithmetic sum, not an independent recount, weighted value, or classification score.
_Avoid_: Variant score, filtered support, recounted support

**Unsupported variant**:
A known splice variant whose defining splice junction has no exact match in the sample's STAR junction file. It is reported with zero support and does not imply a clinical negative or that the variant was not examined.
_Avoid_: Negative variant, absent variant, untested variant

**Variant-definition catalog**:
A human-maintained collection of known splice variants whose defining splice junctions all use one declared genome assembly. Variant identifiers and defining splice junctions are each unique, preserving a one-to-one relationship.
_Avoid_: Variant database, genome library

**GRCh38 reference catalog**:
The versioned sjsift catalog containing six known EGFR variants together with `METex14` and `ARv7`. It remains an explicit input rather than a hidden application default.
_Avoid_: Built-in database, comprehensive splice-variant catalog

## Initial variant identifiers

**EGFRvIII**:
The EGFR splice variant defined by the exon 1-to-exon 8 junction.
_Avoid_: EGFR variant III

**Supported EGFR variants**:
The six EGFR variants included in the v0.1 reference catalog: `EGFRvII`, `EGFRvIIb`, `EGFRvIII`, `EGFRvIIIb`, `EGFRvIVa`, and `EGFRvIVb`.
_Avoid_: EGFR variant set

**METex14**:
The project identifier for MET exon 14 skipping, defined by the exon 13-to-exon 15 junction.
_Avoid_: MET exon 14 deletion, when describing the RNA event

**ARv7**:
The AR splice variant defined by the exon 3-to-cryptic exon 3 junction.
_Avoid_: AR variant 7
