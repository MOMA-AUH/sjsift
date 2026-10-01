"""Structural inspection through real alignments and the public report boundary."""
from dataclasses import replace
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
