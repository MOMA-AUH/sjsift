"""Paired inspection and styling through real alignments and generated reports."""
from html.parser import HTMLParser

import pytest

from alignment_fixtures import CATALOG, bam, cram, fasta, record
from sjsift.alignment_evidence import extract_evidence
from test_read_inspection import render


class Reads(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.rows = []
        self.summaries = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        classes = attrs.get("class", "").split()
        if "read-row" in classes:
            self.rows.append(attrs)
        if "read-hit" in classes:
            self.summaries.append(attrs)


@pytest.mark.parametrize("nh,flag,style,reason", [
    (1, 0, "solid", "unique"),
    (2, 0, "hatched", "multimapping"),
    (None, 0, "muted", "unknown"),
    (1, 256, "muted", "secondary"),
    (2, 256, "muted", "secondary"),
    (1, 2048, "muted", "supplementary"),
    (2, 1024, "muted", "duplicate"),
    (1, 512, "muted", "QC failure"),
])
def test_record_style_and_hover_explain_mapping_and_flags(tmp_path, nh, flag, style, reason):
    path = bam(tmp_path, [record("inspect", flag=flag, tags={} if nh is None else {"NH": nh})])
    evidence = extract_evidence(CATALOG, path)
    html = render(CATALOG, evidence)
    parsed = Reads(html)
    assert len(parsed.rows) == 1
    assert f"read-{style}" in parsed.rows[0]["class"].split()
    assert reason in parsed.summaries[0]["title"]
    assert reason in parsed.summaries[0]["aria-label"]
    assert "NH: " + (str(nh) if nh else "unavailable/unusable") in parsed.summaries[0]["title"]
    assert sum(evidence.groups["v0-j0"].eligible.values()) == 1
    assert "Displayed alignments may include records not counted in STAR junction support" in html


def mates(name="pair", tags=None, start=70, second_start=80):
    first = record(name, start=start, cigar=f"{100 - start}M100N20M", flag=99,
                   tags={"NH": 1} if tags is None else tags)
    second = record(name, start=second_start, cigar=f"{100 - second_start}M100N30M", flag=147,
                    tags={"NH": 1} if tags is None else tags)
    for read, mate in ((first, second), (second, first)):
        read.next_reference_id = mate.reference_id
        read.next_reference_start = mate.reference_start
    return first, second


def test_reciprocal_mates_are_associated_without_collapsing_records_or_support(tmp_path):
    evidence = extract_evidence(CATALOG, bam(tmp_path, mates()))
    html = render(CATALOG, evidence)
    rows = Reads(html).rows
    assert len(rows) == 2
    assert rows[0]["data-pair"] == rows[1]["data-pair"] != ""
    assert rows[0]["data-record-id"] != rows[1]["data-record-id"]
    assert evidence.groups["v0-j0"].eligible == {"unique": 2, "multimapping": 0, "unknown": 0}
    assert 'class="read-pairs"' in html and "View as pairs" in html
    assert "8 unique + 2 multi (STAR)" in html


@pytest.mark.parametrize("problem", [
    "missing_mate", "nonreciprocal", "different_group", "ambiguous", "missing_hi",
    "different_hi", "invalid_hi", "wrong_orientation", "secondary_mismatch", "supplementary",
])
def test_uncertain_mate_associations_remain_individual_records(tmp_path, problem):
    first, second = mates(tags={"NH": 2, "HI": 1})
    reads = [first, second]
    if problem == "missing_mate": reads = [first]
    if problem == "nonreciprocal": second.next_reference_start = 60
    if problem == "different_group": second.set_tag("RG", "other")
    if problem == "ambiguous": reads.append(mates(tags={"NH": 2, "HI": 1})[1])
    if problem == "missing_hi": second.set_tag("HI", None)
    if problem == "different_hi": second.set_tag("HI", 2)
    if problem == "invalid_hi":
        for read in reads: read.set_tag("HI", "1")
    if problem == "wrong_orientation": second.flag &= ~16
    if problem == "secondary_mismatch": second.flag |= 256
    if problem == "supplementary":
        for read in reads: read.flag |= 2048
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}],
              "RG": [{"ID": "other"}]}
    evidence = extract_evidence(CATALOG, bam(tmp_path, reads, header=header))
    rows = Reads(render(CATALOG, evidence)).rows
    assert len(rows) == len(reads)
    assert all(row["data-pair"] == "" for row in rows)


def test_chimeric_tag_is_inspectable_without_inferring_star_exclusion_or_muting(tmp_path):
    read = record('chimeric <"read">', tags={"NH": 1})
    read.set_tag("ch", "1", "A")
    html = render(CATALOG, extract_evidence(CATALOG, bam(tmp_path, [read])))
    parsed = Reads(html)
    assert "read-solid" in parsed.rows[0]["class"].split()
    assert "ch:A:1" in parsed.summaries[0]["title"]
    assert "ch:A:1" in parsed.summaries[0]["aria-label"]
    assert 'chimeric <"read">' not in html


@pytest.mark.parametrize("format", ["BAM", "CRAM"])
def test_placements_and_read_groups_keep_distinct_pairs_and_shared_memberships(tmp_path, format):
    from dataclasses import replace
    from test_reference_alignments import combined_report, reference

    catalog = replace(CATALOG, variants=(replace(CATALOG.variants[0], reference_junctions=(reference(),)),))
    reads = []
    # Identical coordinates/names are distinguishable by read group and HI.
    for group in ("a", "b"):
        for hi in (0, 1):
            pair = mates(tags={"NH": 2, "HI": hi, "RG": group})
            for read in pair:
                read.cigarstring += "100N20M"
                read.query_sequence += "A" * 20
                if hi == 1: read.flag |= 256
            reads.extend(pair)
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}],
              "RG": [{"ID": "a"}, {"ID": "b"}]}
    path = bam(tmp_path, reads, header=header)
    options = {}
    if format == "CRAM":
        options["reference_path"] = fasta(tmp_path)
        path = cram(tmp_path, path, options["reference_path"])
    evidence = extract_evidence(catalog, path, **options)
    rows = Reads(combined_report(catalog, evidence)).rows
    assert len(rows) == 8
    pairs = [row["data-pair"] for row in rows]
    assert "" not in pairs and len(set(pairs)) == 4
    assert all(pairs.count(pair) == 2 for pair in pairs)
    assert sum(row["data-support"] == "0 1" for row in rows) == 4
    records = {record.record_id: record for record in evidence.groups["v0-j0"].records}
    for pair in set(pairs):
        members = [records[row["data-record-id"]] for row in rows if row["data-pair"] == pair]
        assert len({(member.read_group, member.hit_index) for member in members}) == 1
    assert evidence.groups["v0-j0"].eligible["multimapping"] == 8
    assert evidence.groups["v0-j1"].eligible["multimapping"] == 4
