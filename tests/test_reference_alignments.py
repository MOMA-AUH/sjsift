"""Reference groups are independent observations, including shared physical records."""
from dataclasses import replace
from html.parser import HTMLParser
import io

from alignment_fixtures import ANNOTATION, CATALOG, bam, record
from sjsift.alignment_evidence import extract_evidence
from sjsift.catalog import ReferenceJunction
from sjsift.html_report import write_html
from sjsift.quantify import ReferenceJunctionSupport, VariantSupport


def reference(role="same_acceptor", start=221, end=320, chromosome="chr1", strand="+"):
    return ReferenceJunction(role, chromosome, start, end, strand, annotation=ANNOTATION)


def test_shared_junctions_and_overlapping_fetches_preserve_physical_identity(tmp_path):
    ref = reference()
    first = replace(CATALOG.variants[0], reference_junctions=(ref,))
    second = replace(first, identifier="another", intron_start=401, intron_end=500)
    catalog = replace(CATALOG, variants=(first, second))
    path = bam(tmp_path, [record("identical", cigar="20M100N20M100N20M", tags={"NH": 1}) for _ in range(2)])
    result = extract_evidence(catalog, path)
    assert set(result.groups) == {"v0-j0", "v0-j1", "v1-j0", "v1-j1"}
    for key in ("v0-j0", "v0-j1", "v1-j1"):
        group = result.groups[key]
        assert group.eligible == {"unique": 2, "multimapping": 0, "unknown": 0}
        assert len(group.records) == 2
        assert len({r.record_id for r in group.records}) == 2
        assert all(r.matches == ("v0-j0", "v0-j1", "v1-j1") for r in group.records)
    assert result.groups["v0-j0"].records[0].record_id == result.groups["v1-j1"].records[0].record_id
    assert result.groups["v1-j0"].eligible == {"unique": 0, "multimapping": 0, "unknown": 0}


def test_reference_previews_show_shared_membership_and_escape_custom_text(tmp_path):
    hostile = '</script><img/src=x/onerror=alert(1)>'
    ref = reference(role=hostile)
    variant = replace(CATALOG.variants[0], reference_junctions=(ref,))
    catalog = replace(CATALOG, variants=(variant,))
    path = bam(tmp_path, [record(hostile, cigar="20M100N20M100N20M")])
    evidence = extract_evidence(catalog, path)
    support = VariantSupport(variant, 8, 2, (ReferenceJunctionSupport(ref, 20, 3),))
    def render():
        stream = io.StringIO()
        write_html("synthetic", (support,), stream, alignment_evidence=evidence)
        return stream.getvalue()
    html = render()
    assert html == render()
    assert html.count('<li class="read-row') == 1
    assert 'data-support="0 1"' in html
    assert 'style="--read-color:#3159ad"' in html
    assert 'repeating-linear-gradient' not in html
    assert "Supports: Variant; &lt;/script&gt;" in html
    assert "Also matches: skip · &lt;/script&gt;" in html
    assert hostile not in html and '<img' not in html
    assert "Alignment evidence not requested" not in html


def test_rare_defining_and_abundant_reference_evidence_are_both_fully_retained(tmp_path):
    variant = replace(CATALOG.variants[0], reference_junctions=(reference(),))
    catalog = replace(CATALOG, variants=(variant,))
    path = bam(tmp_path, [record("rare", tags={"NH": 1})] + [record(f"abundant-{i}", start=200, tags={"NH": 1}) for i in range(250)])
    groups = extract_evidence(catalog, path).groups
    assert groups["v0-j0"].eligible["unique"] == 1
    assert [r.read_name for r in groups["v0-j0"].records] == ["rare"]
    assert groups["v0-j1"].eligible["unique"] == 250
    assert len(groups["v0-j1"].records) == 250


def test_reference_on_other_contig_and_strand_uses_its_own_policy(tmp_path):
    variant = replace(CATALOG.variants[0], reference_junctions=(reference(chromosome="chr2", strand="-"),))
    catalog = replace(CATALOG, variants=(variant,))
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}, {"SN": "chr2", "LN": 2000}]}
    path = bam(tmp_path, [record("minus", start=200, contig=1, flag=512, tags={"NH": 1, "TS": "-"}),
                          record("opposite", start=200, contig=1, tags={"TS": "+"})], header=header)
    groups = extract_evidence(catalog, path).groups
    assert sum(groups["v0-j0"].eligible.values()) == 0
    assert groups["v0-j1"].eligible["unique"] == 1
    assert groups["v0-j1"].exclusions["opposite"] == 1
    assert groups["v0-j1"].records[0].orientation == "forward"
    assert groups["v0-j1"].records[0].strand_status == "agreeing"


def test_reference_bam_to_cli_html_preserves_both_tsv_outputs(tmp_path):
    from alignment_fixtures import cli_inputs
    from test_alignment_cli import invoke
    args = cli_inputs(tmp_path)
    definitions = tmp_path / "catalog.toml"
    definitions.write_text(definitions.read_text().replace("reference_junctions = []\n", "") + '''
[[variants.reference_junctions]]
role = "custom comparator"
chromosome = "chr1"
intron_start = 221
intron_end = 320
strand = "+"
reference_transcript = "NM_synthetic.1"
donor_exon = "3"
acceptor_exon = "4"
annotation_source = "Synthetic comparator"
''')
    path = bam(tmp_path, [record("shared", cigar="20M100N20M100N20M")])
    before = tmp_path / "before.context.tsv"
    after = tmp_path / "after.context.tsv"
    base = invoke(*args, "--context-output", before)
    result = invoke(*args, "--alignment", path, "--context-output", after, "--html-output", tmp_path / "report.html")
    assert result.returncode == 0, result.stderr
    assert base.stdout == result.stdout
    assert before.read_bytes() == after.read_bytes()
    html = (tmp_path / "report.html").read_text()
    assert "Also matches: skip · custom comparator" in html


class PileupRows(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.rows = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "li" and "read-row" in attrs.get("class", "").split():
            self.rows.append(attrs)


def combined_report(catalog, evidence):
    stream = io.StringIO()
    write_html("synthetic", tuple(
        VariantSupport(variant, 8, 2, tuple(ReferenceJunctionSupport(ref, 20, 3)
                                          for ref in variant.reference_junctions))
        for variant in catalog.variants
    ), stream, alignment_evidence=evidence)
    return stream.getvalue()


def test_single_figure_colors_three_roles_and_previews_each_independently(tmp_path):
    donor = reference("same_donor", 101, 200)
    acceptor = reference("same_acceptor", 201, 300)
    variant = replace(CATALOG.variants[0], intron_end=300, reference_junctions=(donor, acceptor))
    catalog = replace(CATALOG, variants=(variant,))
    path = bam(tmp_path, [record("variant", cigar="20M200N20M", tags={"NH": 1})]
               + [record(f"donor-{i}", tags={"NH": 2}) for i in range(60)]
               + [record("acceptor", start=180)])
    html = combined_report(catalog, extract_evidence(catalog, path))
    rows = PileupRows(html).rows
    assert html.count('class="alignment-group"') == 1
    assert html.count('aria-label="Shared genomic coordinate ruler"') == 1
    assert len(rows) == 62
    assert len([row for row in rows if 'hidden' not in row]) == 12
    variant_rows = [row for row in rows if row['data-support'] == '0']
    donor_rows = [row for row in rows if row['data-support'] == '1']
    acceptor_rows = [row for row in rows if row['data-support'] == '2']
    assert len(variant_rows) == len(acceptor_rows) == 1 and len(donor_rows) == 60
    assert '--read-color:#3159ad' in variant_rows[0]['style']
    assert '--read-color:#168375' in donor_rows[0]['style']
    assert '--read-color:#9950a5' in acceptor_rows[0]['style']
    assert variant_rows[0]['data-mapping'] == 'unique'
    assert donor_rows[0]['data-mapping'] == 'multimapping'
    assert acceptor_rows[0]['data-mapping'] == 'unknown'
    assert 'Variant · exon 1 → 3' in html
    assert 'Same donor · exon 1 → 3' in html
    assert 'Same acceptor · exon 1 → 3' in html
    assert '8 unique + 2 multi (STAR)' in html


def test_combined_figure_keeps_other_contigs_separate_and_empty_guides_visible(tmp_path):
    ref = reference(chromosome="chr2", strand="-")
    variant = replace(CATALOG.variants[0], reference_junctions=(ref,))
    catalog = replace(CATALOG, variants=(variant,))
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}, {"SN": "chr2", "LN": 2000}]}
    path = bam(tmp_path, [record("minus", start=200, contig=1, tags={"NH": 1, "TS": "-"})], header=header)
    html = combined_report(catalog, extract_evidence(catalog, path))
    assert html.count('class="alignment-group"') == 1
    assert html.count('class="pileup-track"') == 2
    assert 'chr1 · genomic coordinates →' in html and 'chr2 · genomic coordinates →' in html
    assert html.count('<li class="read-row') == 1
    assert 'data-support="1"' in html
    assert '>←</span>' in html
