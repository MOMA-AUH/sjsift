"""Structural inspection through real alignments and the public report boundary."""
from dataclasses import replace
from html.parser import HTMLParser
import io

import pytest

from alignment_fixtures import CATALOG, bam, record
from sjsift.alignment_evidence import extract_evidence
from sjsift.html_report import write_html
from sjsift.quantify import VariantSupport


def render(catalog, evidence):
    output = io.StringIO()
    write_html("synthetic", tuple(VariantSupport(v, 8, 2) for v in catalog.variants),
               output, alignment_evidence=evidence)
    return output.getvalue()


class PileupParser(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.rows = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        classes = attrs.get("class", "").split()
        if tag == "li" and "read-row" in classes:
            self.rows.append([])
        if "cigar-op" in classes:
            style = dict(item.split(":") for item in attrs["style"].split(";"))
            self.rows[-1].append((attrs["title"], float(style["left"].strip("%")),
                                 float(style["width"].strip("%"))))


def test_pileup_aligns_shared_coordinates_and_preserves_relative_anchor_lengths(tmp_path):
    catalog = replace(CATALOG, variants=(replace(CATALOG.variants[0], intron_end=1200),))
    path = bam(tmp_path, [
        record("early", start=60, cigar="4S40M1100N20M3S", tags={"NH": 1}),
        record("late", start=80, cigar="10M2I10M1100N40M", tags={"NH": 2}),
    ])
    html = render(catalog, extract_evidence(catalog, path))
    early, late = PileupParser(html).rows
    left_anchor = early[1]
    late_anchor = late[0]
    assert left_anchor[1] < late_anchor[1]
    assert left_anchor[2] == pytest.approx(late_anchor[2] * 4, abs=0.00001)
    # The same junction and downstream genomic start occupy identical positions.
    assert early[2][1:] == late[3][1:]
    assert early[3][1] == late[4][1]
    assert early[3][2] * 2 == pytest.approx(late[4][2], abs=0.00001)
    # Clips/insertions mark boundaries without consuming reference width.
    assert early[0][2] == late[1][2] == early[-1][2] == 0
    assert 'aria-label="Shared genomic coordinate ruler"' in html
    assert 'Compressed gap: 101–1,200 (1,100 nt)' in html
    # Names and metadata are inside a closed native disclosure, below its summary.
    assert html.index('</summary><div class="read-metadata">') < html.index('<strong>early</strong>')
    assert '<details class="read-inspection" open' not in html


def test_pileup_never_compresses_across_another_embedded_reads_aligned_bases(tmp_path):
    catalog = replace(CATALOG, variants=(replace(CATALOG.variants[0], intron_start=501, intron_end=1600),))
    path = bam(tmp_path, [
        record("extra-splice", start=100, cigar="20M350N30M1100N20M"),
        record("continuous", start=100, cigar="400M1100N40M"),
    ])
    html = render(catalog, extract_evidence(catalog, path))
    split, continuous = PileupParser(html).rows
    assert split[1][2] == pytest.approx(continuous[0][2] * 350 / 400, abs=0.00001)
    assert split[3][1:] == continuous[1][1:]
    assert 'Compressed gap: 121–470' not in html
    assert 'Compressed gap: 501–1,600' in html


@pytest.mark.parametrize("cigar,start,anchors", [
    ("10M5=5X100N7=3X10M", 80, (20, 20)),
    ("10M2I5=2D3X100N7M2I10M", 80, (3, 7)),
    ("4S20M100N2D20M3S", 80, (20, 0)),
    ("20M100N100N20M", 80, (20, 0)),
])
def test_anchors_are_only_contiguous_aligned_bases_abutting_selected_skip(tmp_path, cigar, start, anchors):
    path = bam(tmp_path, [record(cigar=cigar, start=start)])
    evidence = extract_evidence(CATALOG, path)
    row = evidence.groups["v0-j0"].records[0]
    assert row.anchors == anchors
    html = render(CATALOG, evidence)
    assert f"Aligned anchors: genomic left {anchors[0]} nt; genomic right {anchors[1]} nt" in html
    assert "not STAR maximum overhang" in html


def test_long_gap_geometry_coordinates_flags_and_safe_offline_report(tmp_path):
    catalog = replace(CATALOG, variants=(replace(CATALOG.variants[0], intron_end=1200),))
    hostile = '</script><img/src=https://example.invalid/onerror=alert(1)>'
    path = bam(tmp_path, [record(hostile, cigar="2H4S10M2I8=2X1100N3M2D7M3S1H", flag=3840, tags={"NH": 1})])
    evidence = extract_evidence(catalog, path)
    html = render(catalog, evidence)
    assert html == render(catalog, evidence)
    assert 'class="gap-break"' in html and "1100 nt; 101–1200" in html
    assert "Nonuniform scale" in html
    for operation in ("M", "equal", "X", "N", "I", "D", "S", "H"):
        assert f"op-{operation}" in html
    for flag in ("secondary", "supplementary", "duplicate", "QC failure"):
        assert flag in html
    assert hostile not in html and "&lt;/script&gt;" in html
    assert "ACGTACGTACGT" not in html and "IIIIIIIIIIII" not in html
    assert "Content-Security-Policy" in html and "connect-src 'none'" in html


def test_sampled_counts_initial_preview_and_controls_are_separate(tmp_path):
    path = bam(tmp_path, [record(f"read-{i}", tags={"NH": 1}) for i in range(70)]
               + [record("unknown")])
    html = render(CATALOG, extract_evidence(CATALOG, path, limit=55))
    assert "Eligible: 71 · Embedded: 56" in html and "Currently shown: 10" in html
    assert "Sampled evidence" in html and "55 records per junction per mapping class" in html
    assert '<option value="all" selected>All</option>' in html
    assert '<option value="unique">Unique only</option>' in html
    assert "Expand reads" in html and "Previous reads" in html and "Next reads" in html
    assert html.count('<li class="read-row') == 56
    # Static HTML itself provides only ten immediate previews without JavaScript.
    rows = html.split('<li class="read-row')[1:]
    assert sum(' hidden' not in row.split('>', 1)[0] for row in rows) == 10
