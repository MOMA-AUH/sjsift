"""Structural inspection through real alignments and the public report boundary."""
from dataclasses import replace
from html.parser import HTMLParser
import io
from pathlib import Path
import re

import pytest

from alignment_fixtures import CATALOG, bam, record
from sjsift.alignment_evidence import extract_evidence
from sjsift.catalog import ReferenceJunction, load_catalog
from sjsift.html_report import write_html
from sjsift.quantify import ReferenceJunctionSupport, VariantSupport


def render(catalog, evidence):
    output = io.StringIO()
    write_html("synthetic", tuple(VariantSupport(v, 8, 2) for v in catalog.variants),
               output, alignment_evidence=evidence)
    return output.getvalue()


class PileupParser(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.rows = []
        self.arcs = []
        self.exon_edges = []
        self.exon_boundaries = []
        self.axis_paths = []
        self.omitted_exons = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        classes = attrs.get("class", "").split()
        if tag == "path" and "junction-arc" in classes:
            self.arcs.append(attrs)
        if tag in {"rect", "path"} and "exon-edge" in classes:
            self.exon_edges.append(attrs)
        if tag == "path" and "exon-boundary" in classes:
            self.exon_boundaries.append(attrs)
        if tag == "path" and set(classes) & {"axis-line", "gap-break"}:
            self.axis_paths.append(attrs)
        if tag == "path" and "omitted-exons" in classes:
            self.omitted_exons.append(attrs)
        if tag == "li" and "read-row" in classes:
            self.rows.append([])
        if "cigar-op" in classes:
            style = dict(item.split(":") for item in attrs["style"].split(";"))
            self.rows[-1].append((attrs["title"], float(style["left"].strip("%")),
                                 float(style["width"].strip("%"))))


def test_ar_exon_three_matches_a_read_spanning_its_full_annotated_extent(tmp_path):
    catalog = load_catalog(Path(__file__).parents[1] / "definitions/grch38.toml")
    catalog = replace(catalog, variants=(catalog.variants[0],))
    path = bam(tmp_path, [record(start=67686009, cigar="117M8546N20M")],
               header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chrX", "LN": 100000000}]})
    html = render(catalog, extract_evidence(catalog, path))
    parsed = PileupParser(html)
    exon = parsed.exon_edges[0]
    aligned = parsed.rows[0][0]
    assert float(exon["x"]) / 10 == pytest.approx(aligned[1], abs=0.00001)
    assert float(exon["width"]) / 10 == pytest.approx(aligned[2], abs=0.00001)
    assert "67,686,010–67,686,126 (117 nt)" in html
    assert len(parsed.exon_boundaries) == 1  # CE3 has no asserted full extent.


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("empty", [False, True])
def test_full_exons_preserve_uncovered_bases_on_both_strands(tmp_path, strand, empty):
    # Only 20 nt of each exon is covered; the rest must retain the same base scale.
    extents = ((1, 100), (201, 400))
    donor, acceptor = extents if strand == "+" else extents[::-1]
    annotation = replace(CATALOG.variants[0].annotation, donor_extent=donor, acceptor_extent=acceptor)
    variant = replace(CATALOG.variants[0], strand=strand, annotation=annotation)
    catalog = replace(CATALOG, variants=(variant,))
    path = bam(tmp_path, [] if empty else [record(tags={"TS": strand})])
    html = render(catalog, extract_evidence(catalog, path))
    parsed = PileupParser(html)
    left, right = parsed.exon_edges
    assert not parsed.exon_boundaries
    assert float(right["width"]) == pytest.approx(2 * float(left["width"]), abs=0.00001)
    assert float(left["x"]) == pytest.approx(10)
    assert float(right["x"]) + float(right["width"]) == pytest.approx(990)
    assert html.count('<path class="gap-break"') == 1  # Only the intron is compressed.
    assert "Compressed gap: 101–200 (100 nt)" in html
    arc = list(map(float, re.findall(r"-?\d+(?:\.\d+)?", parsed.arcs[0]["d"])))
    assert float(left["x"]) + float(left["width"]) == pytest.approx(arc[0], abs=0.00001)
    assert float(right["x"]) == pytest.approx(arc[-2], abs=0.00001)
    if not empty:
        aligned_left, _, aligned_right = parsed.rows[0]
        assert float(left["width"]) / 10 == pytest.approx(5 * aligned_left[2], abs=0.00001)
        assert float(right["width"]) / 10 == pytest.approx(10 * aligned_right[2], abs=0.00001)


def test_known_exon_shared_by_two_junctions_is_drawn_once(tmp_path):
    annotation = CATALOG.variants[0].annotation
    first = replace(CATALOG.variants[0], annotation=replace(annotation, acceptor_extent=(201, 400)))
    reference = ReferenceJunction("next", "chr1", 401, 500, "+",
                                 annotation=replace(annotation, donor_exon="3", donor_extent=(201, 400)))
    variant = replace(first, reference_junctions=(reference,))
    catalog = replace(CATALOG, variants=(variant,))
    path = bam(tmp_path, [])
    support = VariantSupport(variant, 0, 0, (ReferenceJunctionSupport(reference, 0, 0),))
    stream = io.StringIO()
    write_html("synthetic", (support,), stream, alignment_evidence=extract_evidence(catalog, path))
    parsed = PileupParser(stream.getvalue())
    assert len(parsed.exon_edges) == 1
    assert len(parsed.exon_boundaries) == 2


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("empty", [False, True])
def test_long_outer_exons_are_trimmed_to_local_evidence(tmp_path, strand, empty):
    left, right = (1, 1000), (1101, 2000)
    donor, acceptor = (left, right) if strand == "+" else (right, left)
    annotation = replace(CATALOG.variants[0].annotation, donor_extent=donor, acceptor_extent=acceptor)
    variant = replace(CATALOG.variants[0], intron_start=1001, intron_end=1100,
                      strand=strand, annotation=annotation)
    catalog = replace(CATALOG, variants=(variant,))
    path = bam(tmp_path, [] if empty else [record(start=980, tags={"TS": strand})])
    html = render(catalog, extract_evidence(catalog, path))
    parsed = PileupParser(html)
    assert len(parsed.exon_edges) == 2
    assert all("exon-truncated" in exon["class"] for exon in parsed.exon_edges)
    ticks = [int(label.replace(",", "")) for label in re.findall(
        r'<span class="pileup-tick[^>]*>([^<]+)</span>', html)]
    assert ticks[0] > 700 and ticks[-1] < 1400
    assert "1–1,000 (1,000 nt); annotated exon extent" in html
    assert "1,101–2,000 (900 nt); annotated exon extent" in html
    if not empty:
        aligned_left, _, aligned_right = parsed.rows[0]
        assert aligned_left[2] > 5  # Read anchors remain legible after trimming tails.
        assert aligned_left[2] == pytest.approx(aligned_right[2], abs=0.00001)


@pytest.mark.parametrize("covered_middle", [False, True])
def test_long_internal_exon_cuts_only_uncovered_sequence(tmp_path, covered_middle):
    annotation = CATALOG.variants[0].annotation
    first = replace(CATALOG.variants[0], annotation=replace(annotation, acceptor_extent=(201, 1600)))
    reference = ReferenceJunction("next", "chr1", 1601, 1700, "+",
                                 annotation=replace(annotation, donor_exon="3", donor_extent=(201, 1600)))
    variant = replace(first, reference_junctions=(reference,))
    catalog = replace(CATALOG, variants=(variant,))
    reads = ([record("continuous", cigar="20M100N1400M100N20M")] if covered_middle else
             [record("left"), record("right", start=1580)])
    path = bam(tmp_path, reads)
    support = VariantSupport(variant, 8, 0, (ReferenceJunctionSupport(reference, 8, 0),))
    stream = io.StringIO()
    write_html("synthetic", (support,), stream, alignment_evidence=extract_evidence(catalog, path))
    html = stream.getvalue()
    parsed = PileupParser(html)
    assert len(parsed.exon_edges) == (1 if covered_middle else 2)
    assert all(("exon-truncated" in exon["class"]) != covered_middle for exon in parsed.exon_edges)
    assert ("Compressed gap: 341–1,460 (1,120 nt)" in html) != covered_middle
    assert "201–1,600 (1,400 nt); annotated exon extent" in html
    if not covered_middle:
        assert parsed.rows[0][0][2] == pytest.approx(parsed.rows[1][0][2], abs=0.00001)


def test_ar_exon_eight_tail_does_not_squeeze_reference_reads(tmp_path):
    catalog = load_catalog(Path(__file__).parents[1] / "definitions/grch38.toml")
    variant = next(variant for variant in catalog.variants if variant.identifier == "ARv567es")
    catalog = replace(catalog, variants=(variant,))
    reads = [record(ref.role, start=ref.intron_start - 21,
                    cigar=f"20M{ref.intron_end - ref.intron_start + 1}N20M")
             for ref in variant.reference_junctions]
    path = bam(tmp_path, reads,
               header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chrX", "LN": 100000000}]})
    stream = io.StringIO()
    support = VariantSupport(variant, 0, 0, tuple(
        ReferenceJunctionSupport(ref, 20, 0) for ref in variant.reference_junctions))
    write_html("GRCh38", (support,), stream, alignment_evidence=extract_evidence(catalog, path))
    html = stream.getvalue()
    parsed = PileupParser(html)
    assert any("exon-truncated" in exon["class"] for exon in parsed.exon_edges)
    assert "67,723,686–67,730,619 (6,934 nt)" in html
    ticks = [int(label.replace(",", "")) for label in re.findall(
        r'<span class="pileup-tick[^>]*>([^<]+)</span>', html)]
    assert ticks[-1] == 67723825  # 20 aligned bases + 120 nt of context into exon 8.
    for row in parsed.rows:
        anchors = [op for op in row if op[0].startswith("M:")]
        assert all(anchor[2] > 2 for anchor in anchors)


@pytest.mark.parametrize("identifier,omitted", [("BRAFdel2-10", "3–9"), ("ARv567es", "6")])
def test_intron_marks_numbered_exons_omitted_between_displayed_neighbors(tmp_path, identifier, omitted):
    catalog = load_catalog(Path(__file__).parents[1] / "definitions/grch38.toml")
    variant = next(v for v in catalog.variants if v.identifier == identifier)
    catalog = replace(catalog, variants=(variant,))
    path = bam(tmp_path, [], header={"HD": {"SO": "coordinate"},
                                   "SQ": [{"SN": variant.chromosome, "LN": 250000000}]})
    support = VariantSupport(variant, 0, 0, tuple(
        ReferenceJunctionSupport(ref, 0, 0) for ref in variant.reference_junctions))
    stream = io.StringIO()
    write_html("GRCh38", (support,), stream, alignment_evidence=extract_evidence(catalog, path))
    html = stream.getvalue()
    markers = PileupParser(html).omitted_exons
    assert len(markers) == 1  # The adjacent reference connections stay straight.
    assert f"Exons omitted from illustration: {omitted}" in html
    transcript = "NM_004333.6" if identifier.startswith("BRAF") else "NM_000044.6"
    assert f"Exons omitted from illustration: {omitted} ({transcript})" in html
    assert re.findall(r"[A-Za-z]", markers[0]["d"]) == ["M", "L", "L", "L", "L"]


@pytest.mark.parametrize("different_transcript", [False, True])
def test_omission_marks_require_numbered_exons_from_the_same_transcript(tmp_path, different_transcript):
    annotation = replace(CATALOG.variants[0].annotation, donor_exon="1", acceptor_exon="CE3",
                         donor_extent=(1, 100), acceptor_extent=None)
    variant = replace(CATALOG.variants[0], annotation=annotation)
    if different_transcript:
        other = replace(annotation, reference_transcript="other", donor_exon="3", donor_extent=(201, 300))
        variant = replace(variant, reference_junctions=(
            ReferenceJunction("next", "chr1", 301, 400, "+", annotation=other),))
    catalog = replace(CATALOG, variants=(variant,))
    path = bam(tmp_path, [])
    support = VariantSupport(variant, 0, 0, tuple(
        ReferenceJunctionSupport(ref, 0, 0) for ref in variant.reference_junctions))
    stream = io.StringIO()
    write_html("synthetic", (support,), stream, alignment_evidence=extract_evidence(catalog, path))
    assert not PileupParser(stream.getvalue()).omitted_exons


@pytest.mark.parametrize("strand", ["+", "-"])
def test_sashimi_arcs_align_with_reads_and_keep_zero_star_support_visible(tmp_path, strand):
    annotation = CATALOG.variants[0].annotation
    donor = ReferenceJunction("same_donor", "chr1", 101, 200, strand, annotation=annotation)
    acceptor = ReferenceJunction("same_acceptor", "chr1", 201, 300, strand, annotation=annotation)
    variant = replace(CATALOG.variants[0], intron_end=300, strand=strand,
                      reference_junctions=(donor, acceptor))
    catalog = replace(CATALOG, variants=(variant,))
    path = bam(tmp_path, [record("variant", cigar="20M200N20M"),
                          record("donor"), record("acceptor", start=180)])
    support = VariantSupport(variant, 8, 2, (ReferenceJunctionSupport(donor, 20, 3),
                                           ReferenceJunctionSupport(acceptor, 0, 0)))
    stream = io.StringIO()
    write_html("synthetic", (support,), stream, alignment_evidence=extract_evidence(catalog, path))
    html = stream.getvalue()
    parsed = PileupParser(html)
    arcs = {arc["stroke"]: arc for arc in parsed.arcs}
    assert len(arcs) == 3
    # The defining junction shares an edge with each reference junction.
    assert len(parsed.exon_boundaries) == 4
    assert not parsed.exon_edges
    baselines = set()
    for color, row in zip(("#3159ad", "#168375", "#9950a5"), parsed.rows):
        coordinates = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", arcs[color]["d"])]
        _, left, width = next(op for op in row if op[0].startswith("N:"))
        assert coordinates[0] / 10 == pytest.approx(left, abs=0.00001)
        assert coordinates[-2] / 10 == pytest.approx(left + width, abs=0.00001)
        baselines.add(coordinates[-1])
    assert len(baselines) == 1
    assert float(arcs["#168375"]["stroke-width"]) > float(arcs["#3159ad"]["stroke-width"])
    assert 'zero-support' in arcs["#9950a5"]["class"]
    assert '0 unique + 0 multi (STAR)' in html
    assert 'arc-direction' not in html
    assert 'aria-label="Transcript direction"' not in html
    assert html.index('class="junction-arcs"') < html.index('class="pileup-ruler"') < html.index('<ol class="read-rows">')


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
    # Soft clips use the aligned-base scale; insertions remain boundary markers.
    assert early[0][2] == pytest.approx(left_anchor[2] * 4 / 40, abs=0.00001)
    assert early[-1][2] == pytest.approx(early[3][2] * 3 / 20, abs=0.00001)
    assert early[0][1] + early[0][2] == pytest.approx(left_anchor[1], abs=0.00001)
    assert early[-1][1] == pytest.approx(early[3][1] + early[3][2], abs=0.00001)
    assert late[1][2] == 0
    assert 'aria-label="Shared genomic coordinate ruler"' in html
    assert 'Compressed gap: 101–1,200 (1,100 nt)' in html
    # Names and metadata are inside a closed native disclosure, below its summary.
    assert html.index('</summary><div class="read-metadata">') < html.index('<strong>early</strong>')
    assert '<details class="read-inspection" open' not in html


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("empty", [False, True])
def test_exon_edges_follow_genomic_sides_and_gaps_form_a_straight_connected_axis(tmp_path, strand, empty):
    variant = replace(CATALOG.variants[0], strand=strand)
    catalog = replace(CATALOG, variants=(variant,))
    path = bam(tmp_path, [] if empty else [record(tags={"TS": strand})])
    html = render(catalog, extract_evidence(catalog, path))
    parsed = PileupParser(html)
    left, right = (list(map(float, re.findall(r"-?\d+(?:\.\d+)?", marker["d"])))
                   for marker in parsed.exon_boundaries)
    arc = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", parsed.arcs[0]["d"])]
    # The exon-facing sides end/start exactly where the junction arc lands.
    assert left[0] == pytest.approx(arc[0], abs=0.00001)
    assert right[0] == pytest.approx(arc[-2], abs=0.00001)
    assert 0 <= left[0] < right[0] <= 1000
    assert all(re.findall(r"[A-Za-z]", marker["d"]) == ["M", "V"] for marker in parsed.exon_boundaries)
    assert not parsed.exon_edges
    labels = re.findall(r'<span class="exon-edge-label exon-boundary-label"[^>]*>([^<]+)</span>', html)
    assert labels == (["1", "3"] if strand == "+" else ["3", "1"])
    assert "full exon extent unknown" in html
    assert len(parsed.axis_paths) == 3
    first, gap, last = parsed.axis_paths
    start, _, stop = map(float, re.findall(r"-?\d+(?:\.\d+)?", first["d"]))
    gap_start, gap_y, gap_end = map(float, re.findall(r"-?\d+(?:\.\d+)?", gap["d"]))
    next_start, _, end = map(float, re.findall(r"-?\d+(?:\.\d+)?", last["d"]))
    assert start <= stop == gap_start < gap_end == next_start <= end
    assert gap_y == 12
    assert gap["class"] == "gap-break"
    assert all(re.findall(r"[A-Za-z]", path["d"]) == ["M", "H"] for path in parsed.axis_paths)
    assert "Compressed gap: 101–200 (100 nt)" in html
    assert "//</span>" not in html


@pytest.mark.parametrize("flag", [0, 16])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_long_soft_clips_fit_at_both_ends_and_keep_cigar_order_on_either_strand(tmp_path, flag, strand):
    catalog = replace(CATALOG, variants=(replace(CATALOG.variants[0], strand=strand),))
    path = bam(tmp_path, [record(cigar="3H120S20M100N20M240S2H", flag=flag, tags={"TS": strand})])
    evidence = extract_evidence(catalog, path)
    html = render(catalog, evidence)
    hard_left, soft_left, aligned_left, gap, aligned_right, soft_right, hard_right = PileupParser(html).rows[0]
    assert soft_left[2] == pytest.approx(aligned_left[2] * 6, abs=0.00001)
    assert soft_right[2] == pytest.approx(aligned_right[2] * 12, abs=0.00001)
    assert soft_left[1] + soft_left[2] == pytest.approx(aligned_left[1], abs=0.00001)
    assert soft_right[1] == pytest.approx(aligned_right[1] + aligned_right[2], abs=0.00001)
    assert soft_left[1] >= 1 and soft_right[1] + soft_right[2] <= 99.00001
    assert hard_left[2] == hard_right[2] == 0
    assert "S: 120 nt; boundary 80; unaligned soft clip" in html
    assert "S: 240 nt; boundary 220; unaligned soft clip" in html
    # Even a clip extending past chromosome start does not create genomic labels.
    ticks = [int(label.replace(",", "")) for label in re.findall(
        r'<span class="pileup-tick[^>]*>([^<]+)</span>', html)]
    assert ticks[0] == 81 and ticks[-1] == 220
    assert ticks == sorted(ticks) and all(81 <= tick <= 220 for tick in ticks)
    assert 'Compressed gap: 101–200' in html
    row = evidence.groups["v0-j0"].records[0]
    assert row.start == 80 and row.end == 220
    assert all(op.start == op.end for op in row.geometry if op.op == "S")


def test_soft_clips_keep_the_base_scale_when_extending_into_another_reads_intron(tmp_path):
    path = bam(tmp_path, [
        record("clipped", cigar="20M100N20M120S"),
        record("extra-splice", cigar="20M100N20M500N20M"),
    ])
    html = render(CATALOG, extract_evidence(CATALOG, path))
    clipped, extended = PileupParser(html).rows
    assert clipped[-1][2] == pytest.approx(clipped[-2][2] * 6, abs=0.00001)
    assert clipped[-1][1] == extended[-2][1]
    assert clipped[0][1:] == extended[0][1:]
    assert clipped[2][1:] == extended[2][1:]
    assert 'Compressed gap: 341–720' in html
    assert 'Compressed gap: 221–720' not in html


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
    assert f"Aligned anchors: left {anchors[0]} nt, right {anchors[1]} nt" in html


def test_long_gap_geometry_coordinates_flags_and_safe_offline_report(tmp_path):
    catalog = replace(CATALOG, variants=(replace(CATALOG.variants[0], intron_end=1200),))
    hostile = '</script><img/src=https://example.invalid/onerror=alert(1)>'
    path = bam(tmp_path, [record(hostile, cigar="2H4S10M2I8=2X1100N3M2D7M3S1H", flag=3840, tags={"NH": 1})])
    evidence = extract_evidence(catalog, path)
    html = render(catalog, evidence)
    assert html == render(catalog, evidence)
    assert 'class="gap-break"' in html and "1100 nt; 101–1200" in html
    for operation in ("M", "equal", "X", "N", "I", "D", "S", "H"):
        assert f"op-{operation}" in html
    for flag in ("secondary", "supplementary", "duplicate", "QC failure"):
        assert flag in html
    assert hostile not in html and "&lt;/script&gt;" in html
    assert "ACGTACGTACGT" not in html and "IIIIIIIIIIII" not in html
    assert "Content-Security-Policy" in html and "connect-src 'none'" in html


def test_all_reads_are_embedded_with_preview_controls_and_without_sampling_text(tmp_path):
    path = bam(tmp_path, [record(f"read-{i}", tags={"NH": 1}) for i in range(270)]
               + [record("unknown")])
    html = render(CATALOG, extract_evidence(CATALOG, path))
    for removed in ("distinct embedded alignment records", "Currently shown:", "Sampling limit:",
                    "Sampled evidence:", "embedded samples are deterministic", "shown-count"):
        assert removed not in html
    assert '<option value="all" selected>All</option>' in html
    assert '<option value="unique">Unique only</option>' in html
    assert "Expand reads" in html and "Previous reads" in html and "Next reads" in html
    assert html.count('<li class="read-row') == 271
    # Static HTML itself provides only ten immediate previews without JavaScript.
    rows = html.split('<li class="read-row')[1:]
    assert sum(' hidden' not in row.split('>', 1)[0] for row in rows) == 10
