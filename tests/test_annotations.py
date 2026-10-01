"""Current catalog annotation and local schematic behavior."""
import io
from pathlib import Path

import pytest

from sjsift.catalog import CatalogError, load_catalog
from sjsift.html_report import write_html
from sjsift.quantify import VariantSupport


ANNOTATED = '''schema_version = 3
genome_assembly = "synthetic"
[[variants]]
id = "terminal"
chromosome = "chr1"
intron_start = 101
intron_end = 200
strand = "-"
reference_junctions = []
reference_transcript = "NM_test.1"
donor_exon = "17"
acceptor_exon = "E18-C3"
annotation_source = "Curated terminal exon; accession is comparison basis"
'''


def test_schema_three_annotation_loads_and_old_versions_are_rejected(tmp_path):
    path = tmp_path / "catalog.toml"
    path.write_text(ANNOTATED)
    catalog = load_catalog(path)
    assert catalog.variants[0].annotation.donor_exon == "17"
    assert catalog.variants[0].annotation.acceptor_exon == "E18-C3"
    for version in (1, 2):
        path.write_text(ANNOTATED.replace("schema_version = 3", f"schema_version = {version}"))
        with pytest.raises(CatalogError, match=f"unsupported catalog schema version {version}; expected 3"):
            load_catalog(path)


@pytest.mark.parametrize("junction_kind", ["defining", "reference"])
@pytest.mark.parametrize("field", ["reference_transcript", "donor_exon", "acceptor_exon", "annotation_source"])
@pytest.mark.parametrize("value", ['""', '" "', '7', 'true', '"bad\\nline"'])
def test_annotations_are_required_nonempty_text(tmp_path, field, value, junction_kind):
    path = tmp_path / "catalog.toml"
    prefix = ""
    target = ANNOTATED
    if junction_kind == "reference":
        prefix = ANNOTATED.replace("reference_junctions = []\n", "")
        target = ANNOTATED.split("[[variants]]", 1)[1].replace('id = "terminal"', 'role = "comparator"').replace("intron_start = 101", "intron_start = 151").replace("reference_junctions = []\n", "")
        prefix += "\n[[variants.reference_junctions]]\n"
    lines = target.splitlines()
    path.write_text(prefix + "\n".join(f"{field} = {value}" if line.startswith(field + " =") else line for line in lines))
    with pytest.raises(CatalogError, match=field):
        load_catalog(path)
    path.write_text(prefix + "\n".join(line for line in lines if not line.startswith(field + " =")))
    with pytest.raises(CatalogError, match=field):
        load_catalog(path)


def test_minus_strand_schematic_labels_local_boundaries_and_comparison_basis(tmp_path):
    path = tmp_path / "catalog.toml"
    path.write_text(ANNOTATED)
    catalog = load_catalog(path)
    stream = io.StringIO()
    write_html(catalog.genome_assembly, (VariantSupport(catalog.variants[0], 2, 3),), stream)
    html = stream.getvalue()
    assert "Donor exon 17" in html and "Acceptor exon E18-C3" in html
    assert "Donor boundary: 200" in html and "Acceptor boundary: 101" in html
    assert "Annotation / comparison basis: NM_test.1" in html
    assert "Local splice boundaries" in html and "not to scale" in html
