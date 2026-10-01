"""Reference groups are independent observations, including shared physical records."""
from dataclasses import replace
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
    assert html.count("Eligible: 1 · Embedded: 1") == 2
    assert "Also matches: skip · Defining junction" in html
    assert "Also matches: skip · &lt;/script&gt;" in html
    assert hostile not in html and '<img' not in html
    assert "Alignment evidence not requested" not in html


def test_rare_defining_evidence_is_independent_from_abundant_reference_samples(tmp_path):
    variant = replace(CATALOG.variants[0], reference_junctions=(reference(),))
    catalog = replace(CATALOG, variants=(variant,))
    path = bam(tmp_path, [record("rare", tags={"NH": 1})] + [record(f"abundant-{i}", start=200, tags={"NH": 1}) for i in range(150)])
    groups = extract_evidence(catalog, path, limit=3).groups
    assert groups["v0-j0"].eligible["unique"] == 1
    assert [r.read_name for r in groups["v0-j0"].records] == ["rare"]
    assert groups["v0-j1"].eligible["unique"] == 150
    assert len(groups["v0-j1"].records) == 3


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
    result = invoke(*args, "--alignments", path, "--context-output", after, "--html-output", tmp_path / "report.html")
    assert result.returncode == 0, result.stderr
    assert base.stdout == result.stdout
    assert before.read_bytes() == after.read_bytes()
    html = (tmp_path / "report.html").read_text()
    assert 'id="reads-v0-j1"' in html and "Also matches: skip · custom comparator" in html
    assert "Eligible: 1 · Embedded: 1" in html
