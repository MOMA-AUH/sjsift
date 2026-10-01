"""Check curated junctions against pinned, independently extracted exon evidence."""

import csv
from pathlib import Path

import pytest

from sjsift.catalog import ReferenceJunction, load_catalog


ROOT = Path(__file__).parents[1]

# Transcript and adjacent exon ranks selected in REFERENCE_JUNCTION_CURATION.md.
# Order is donor, then acceptor; the two terminal-exon events have only a donor.
CURATION = {
    "ARv7": ("ENST00000374690.9", (3,)),
    "ARv567es": ("ENST00000374690.9", (4, 7)),
    "BRAFdel2-10": ("ENST00000646891.2", (1, 10)),
    "BRAFdel2-8": ("ENST00000646891.2", (1, 8)),
    "BRAFdel3-10": ("ENST00000646891.2", (2, 10)),
    "BRAFdel3-8": ("ENST00000646891.2", (2, 8)),
    "BRAFdel4-10": ("ENST00000646891.2", (3, 10)),
    "BRAFdel4-8": ("ENST00000646891.2", (3, 8)),
    "EGFRvIVa": ("ENST00000275493.7", (24, 27)),
    "EGFRvIII": ("ENST00000275493.7", (1, 7)),
    "EGFRvIIIb": ("ENST00000450046.2", (1, 7)),
    "EGFRvIVb": ("ENST00000275493.7", (24, 26)),
    "EGFRvII": ("ENST00000275493.7", (13, 15)),
    "EGFRvIIb": ("ENST00000344576.7", (13, 15)),
    "FGFR2-E18-C3": ("ENST00000358487.10", (17,)),
    "METex7-8": ("ENST00000397752.8", (6, 8)),
    "METex14": ("ENST00000397752.8", (13, 14)),
}


def test_every_reference_catalog_variant_has_curated_context() -> None:
    catalog = load_catalog(ROOT / "definitions" / "grch38.toml")

    assert catalog.genome_assembly == "GRCh38"
    assert tuple(variant.identifier for variant in catalog.variants) == tuple(CURATION)
    assert sum(len(variant.reference_junctions) for variant in catalog.variants) == 32


@pytest.mark.parametrize("identifier", CURATION)
def test_curated_references_match_adjacent_exons_and_shared_splice_sites(
    identifier: str,
) -> None:
    catalog = load_catalog(ROOT / "definitions" / "grch38.toml")
    variant = next(v for v in catalog.variants if v.identifier == identifier)
    transcript, ranks = CURATION[identifier]
    with (ROOT / "docs" / "curation" / "gencode-v49-exons.tsv").open(
        encoding="utf-8", newline=""
    ) as stream:
        exons = {
            int(row["exon_number"]): row
            for row in csv.DictReader(stream, delimiter="\t")
            if row["transcript_id"] == transcript
        }

    expected = []
    for role, rank in zip(("same_donor", "same_acceptor"), ranks):
        upstream, downstream = exons[rank], exons[rank + 1]
        assert upstream["chromosome"] == downstream["chromosome"] == variant.chromosome
        assert upstream["strand"] == downstream["strand"] == variant.strand
        if variant.strand == "+":
            start = int(upstream["exon_end"]) + 1
            end = int(downstream["exon_start"]) - 1
        else:
            start = int(downstream["exon_end"]) + 1
            end = int(upstream["exon_start"]) - 1
        assert start <= end
        if (role == "same_donor") == (variant.strand == "+"):
            assert start == variant.intron_start
            assert end != variant.intron_end
        else:
            assert end == variant.intron_end
            assert start != variant.intron_start
        expected.append((role, variant.chromosome, start, end, variant.strand))

    assert tuple((r.role, r.chromosome, r.intron_start, r.intron_end, r.strand) for r in variant.reference_junctions) == tuple(expected)


@pytest.mark.parametrize("identifier", CURATION)
def test_annotations_match_pinned_exon_boundaries_and_explicit_terminal_curation(identifier):
    variant = next(v for v in load_catalog(ROOT / "definitions/grch38.toml").variants if v.identifier == identifier)
    enst, ranks = CURATION[identifier]
    transcripts = {
        "ENST00000374690.9": "NM_000044.6", "ENST00000646891.2": "NM_004333.6",
        "ENST00000275493.7": "NM_005228.5", "ENST00000450046.2": "NM_001346900.2",
        "ENST00000344576.7": "NM_201284.2", "ENST00000358487.10": "NM_000141.5",
        "ENST00000397752.8": "NM_000245.4",
    }
    with (ROOT / "docs/curation/gencode-v49-exons.tsv").open() as stream:
        exons = {row["exon_number"]: row for row in csv.DictReader(stream, delimiter="\t") if row["transcript_id"] == enst}
    for junction in (variant, *variant.reference_junctions):
        annotation = junction.annotation
        assert annotation.reference_transcript == transcripts[enst]
        assert enst in annotation.annotation_source and "GENCODE v49" in annotation.annotation_source
        donor, acceptor = annotation.donor_exon, annotation.acceptor_exon
        if identifier == "FGFR2-E18-C3":
            assert donor == "E17"
            assert acceptor == ("E18-C3" if junction is variant else "E18-C1")
            assert "not present in NM_000141.5" in annotation.annotation_source
            donor, acceptor = "17", (None if junction is variant else "18")
        if identifier == "ARv7" and junction is variant:
            assert acceptor == "CE3"
            assert "not present in NM_000044.6" in annotation.annotation_source
            acceptor = None
        if junction.strand == "+":
            assert int(exons[donor]["exon_end"]) + 1 == junction.intron_start
            if acceptor:
                assert int(exons[acceptor]["exon_start"]) - 1 == junction.intron_end
        else:
            assert int(exons[donor]["exon_start"]) - 1 == junction.intron_end
            if acceptor:
                assert int(exons[acceptor]["exon_end"]) + 1 == junction.intron_start
