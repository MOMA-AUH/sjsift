"""Behavior at the focused alignment-evidence entry point, using real BAM."""
import pytest

from alignment_fixtures import CATALOG, bam, record
from sjsift.alignment_evidence import extract_evidence


def test_exact_cigar_splice_counts_each_record_and_keeps_star_independent(tmp_path):
    path = bam(tmp_path, [record("exact", tags={"NH": 1, "XS": "+"}),
                          record("near", cigar="20M99N20M"),
                          record("deletion", cigar="20M100D20M")])
    result = extract_evidence(CATALOG, path)
    group = result.groups["v0-j0"]
    assert group.eligible == {"unique": 1, "multimapping": 0, "unknown": 0}
    assert [r.read_name for r in group.records] == ["exact"]
    assert group.records[0].strand_status == "agreeing"
    assert group.records[0].geometry[1].start == 100
    assert group.records[0].geometry[1].end == 200


@pytest.mark.parametrize("tags,flag,status", [
    ({"NH": 1}, 16, "unverified"),
    ({"XS": "+"}, 16, "agreeing"),
    ({"TS": "+"}, 0, "agreeing"),
    ({"XS": "-"}, 0, "opposite"),
    ({"XS": "+", "TS": "-"}, 16, "conflicting"),
    ({"jM": [21], "jI": [101, 200]}, 16, "agreeing"),
    ({"jM": [22]}, 0, "opposite"),
    ({"jM": [1], "XS": "-"}, 0, "conflicting"),
    ({"jM": [0]}, 0, "unverified"),
    ({"jM": [20]}, 0, "unverified"),
    ({"jM": [7]}, 0, "unverified"),
    ({"jM": [1, 2]}, 0, "unverified"),
    ({"jM": [1], "jI": [100, 200]}, 0, "unverified"),
    ({"XS": 1, "TS": "?"}, 0, "unverified"),
])
def test_transcript_strand_sources_are_validated_and_orientation_is_independent(tmp_path, tags, flag, status):
    result = extract_evidence(CATALOG, bam(tmp_path, [record(tags=tags, flag=flag)]))
    group = result.groups["v0-j0"]
    if status in ("opposite", "conflicting"):
        assert group.exclusions[status] == 1
        assert not group.records and sum(group.eligible.values()) == 0
    else:
        assert group.records[0].strand_status == status
        assert group.strand_totals[status] == 1


def test_star_motifs_require_program_provenance(tmp_path):
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}]}
    path = bam(tmp_path, [record(tags={"jM": [1]})], header=header)
    assert extract_evidence(CATALOG, path).groups["v0-j0"].records[0].strand_status == "unverified"


def test_motif_strand_is_interpreted_per_junction(tmp_path):
    from dataclasses import replace
    catalog = replace(CATALOG, variants=(CATALOG.variants[0], replace(CATALOG.variants[0], identifier="second", intron_start=221, intron_end=320, strand="-")))
    path = bam(tmp_path, [record(cigar="20M100N20M100N20M", tags={"jM": [21, 22], "jI": [101, 200, 221, 320]})])
    result = extract_evidence(catalog, path)
    assert all(g.strand_totals == {"agreeing": 1, "unverified": 0} for g in result.groups.values())
    assert all(len(g.records) == 1 for g in result.groups.values())


@pytest.mark.parametrize("nh,expected", [(1, "unique"), (2, "multimapping"), (0, "unknown"), (-1, "unknown"), ("1", "unknown"), (1.0, "unknown"), (None, "unknown")])
def test_nh_alone_defines_mapping_class_even_with_low_mapq_and_flags(tmp_path, nh, expected):
    path = bam(tmp_path, [record(flag=256|512|1024|2048, mapq=0, tags={} if nh is None else {"NH": nh})])
    row = extract_evidence(CATALOG, path).groups["v0-j0"].records[0]
    assert row.mapping_class == expected
    assert row.flags == ("secondary", "QC failure", "duplicate", "supplementary")
    assert row.mapq == 0


def test_sampling_is_bounded_per_class_but_counts_and_identical_records_are_preserved(tmp_path):
    reads = [record("identical", tags={"NH": 1}) for _ in range(7)]
    reads += [record("multi", tags={"NH": 2}), record("unknown")]
    path = bam(tmp_path, reads)
    result = extract_evidence(CATALOG, path, limit=3)
    assert result == extract_evidence(CATALOG, path, limit=3)
    group = result.groups["v0-j0"]
    assert group.eligible == {"unique": 7, "multimapping": 1, "unknown": 1}
    assert len(group.records) == 5
    assert len({r.record_id for r in group.records}) == 5
    assert sum(r.mapping_class == "unique" for r in group.records) == 3
    assert len(extract_evidence(CATALOG, path, limit=10).groups["v0-j0"].records) == 9


def test_structural_operations_retain_coordinates_without_sequences(tmp_path):
    read = record(start=80, cigar="3H4S10M2I5=2D3X100N7M2S", tags={"NH": 1})
    path = bam(tmp_path, [read])
    row = extract_evidence(CATALOG, path).groups["v0-j0"].records[0]
    assert [(o.op, o.start, o.end, o.length) for o in row.geometry] == [
        ("H", 80, 80, 3), ("S", 80, 80, 4), ("M", 80, 90, 10),
        ("I", 90, 90, 2), ("=", 90, 95, 5), ("D", 95, 97, 2),
        ("X", 97, 100, 3), ("N", 100, 200, 100), ("M", 200, 207, 7), ("S", 207, 207, 2),
    ]
    assert "ACGTACGT" not in repr(row) and "query_sequence" not in repr(row)


@pytest.mark.parametrize("reverse,value,status", [(False, "+", "agreeing"), (True, "-", "agreeing"), (True, "+", "opposite")])
def test_minimap_relative_transcript_strand_has_explicit_program_provenance(tmp_path, reverse, value, status):
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}], "PG": [{"ID": "mm2", "PN": "minimap2"}]}
    group = extract_evidence(CATALOG, bam(tmp_path, [record(flag=16 if reverse else 0, tags={"ts": value})], header=header)).groups["v0-j0"]
    if status == "opposite": assert group.exclusions[status] == 1
    else: assert group.records[0].strand_status == status
