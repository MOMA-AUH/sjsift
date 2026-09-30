"""HTML evidence, safe rendering, and coordinated report output behavior."""

from html.parser import HTMLParser
import io
from pathlib import Path

import pytest

from sjsift.catalog import ReferenceJunction, VariantDefinition
from sjsift.html_report import write_html
from sjsift.quantify import ReferenceJunctionSupport, VariantSupport
from sjsift import report


class Document(HTMLParser):
    def __init__(self, text: str) -> None:
        super().__init__(convert_charrefs=True)
        self.tags = []
        self.text = []
        self.scripts = []
        self.in_script = False
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))
        if tag == "script":
            self.in_script = True

    def handle_endtag(self, tag):
        if tag == "script":
            self.in_script = False

    def handle_data(self, data):
        (self.scripts if self.in_script else self.text).append(data)


def support(identifier="skip", unique=18, multi=2, reference_counts=((90, 5), (84, 4))):
    references = (
        ReferenceJunction("same_donor", "chr7", 10, 20, "+"),
        ReferenceJunction("same_acceptor", "chr7", 21, 30, "+"),
    )[:len(reference_counts)]
    variant = VariantDefinition(identifier, "chr7", 10, 30, "+", references)
    return VariantSupport(variant, unique, multi, tuple(
        ReferenceJunctionSupport(ref, u, m)
        for ref, (u, m) in zip(references, reference_counts)
    ))


def render(results, **kwargs):
    stream = io.StringIO()
    write_html("GRCh38", results, stream, **kwargs)
    return stream.getvalue()


def test_html_has_exact_counts_both_comparators_and_coordinates_without_javascript():
    html = render((support(),), junctions_name="sample.SJ.out.tab", definitions_name="catalog.toml")
    doc = Document(html)
    text = " ".join(doc.text)
    for expected in ("skip", "Same donor", "Same acceptor", "sample.SJ.out.tab", "catalog.toml",
                     "GRCh38", "18", "90", "84", "20", "95", "88", "chr7:10–20 (+)"):
        assert expected in text
    assert "Variant and reference support" in text
    assert ("a", {"href": "#variant-0"}) in doc.tags
    assert len([tag for tag, _ in doc.tags if tag == "script"]) == 1
    assert not any("src" in attrs for _, attrs in doc.tags)
    assert not any(tag == "link" for tag, _ in doc.tags)


@pytest.mark.parametrize(("unique", "multi", "references", "label"), [
    (0, 4, ((0, 0),), "Variant support only"),
    (0, 0, ((3, 0),), "Reference support only"),
    (0, 0, ((0, 0),), "No selected junction support"),
    (0, 0, (), "No variant support; no reference context configured"),
    (1, 0, (), "Variant support observed; no reference context configured"),
])
def test_evidence_descriptions_distinguish_zero_multimapping_and_unconfigured_context(
    unique, multi, references, label,
):
    html = render((support(unique=unique, multi=multi, reference_counts=references),))
    assert label in " ".join(Document(html).text)


def test_arbitrary_roles_and_catalog_text_cannot_inject_html_or_javascript():
    hostile = '</script><img src=x onerror="alert(1)"> & $script'
    ref = ReferenceJunction(hostile, hostile, 10, 20, "-")
    variant = VariantDefinition(hostile, hostile, 10, 30, "-", (ref,))
    result = VariantSupport(variant, 3, 1, (ReferenceJunctionSupport(ref, 9, 2),))
    stream = io.StringIO()
    write_html(hostile, (result,), stream, junctions_name=hostile, definitions_name=hostile)
    doc = Document(stream.getvalue())

    assert hostile in "".join(doc.text)
    assert hostile not in "".join(doc.scripts)
    assert not any(tag == "img" or "onerror" in attrs for tag, attrs in doc.tags)
    assert len([tag for tag, _ in doc.tags if tag == "script"]) == 1
    assert not any("$script" in value for _, attrs in doc.tags for value in attrs.values() if value)


def test_report_is_deterministic_and_preserves_large_integer_counts():
    result = support(unique=9007199254740993, reference_counts=())
    first = render((result,), compatibility_warning=True)
    assert first == render((result,), compatibility_warning=True)
    assert 'data-unique="9007199254740993"' in first
    assert "No catalog chromosome identifiers were found" in first


def test_overview_ranks_unique_support_then_total_with_catalog_order_as_tiebreaker():
    rows = (support("first", 0, 8), support("second", 1, 0), support("third", 0, 8))
    doc = Document(render(rows))
    indices = [attrs["data-index"] for tag, attrs in doc.tags if tag == "tr" and "data-index" in attrs]
    assert indices == ["1", "0", "2"]


def test_role_columns_include_custom_roles_and_distinguish_missing_from_zero():
    ref = ReferenceJunction("terminal comparator", "chr10", 10, 20, "-")
    v = VariantDefinition("terminal", "chr10", 15, 20, "-", (ref,))
    rows = (support(), VariantSupport(v, 0, 0, (ReferenceJunctionSupport(ref, 0, 0),)))
    html = render(rows)
    assert "Reference · terminal comparator" in html
    assert "Not configured" in html
    assert "chr10:10–20 (-)" in html


@pytest.mark.parametrize("failed_writer", ["write_tsv", "write_context_tsv", "write_html"])
def test_write_failure_removes_all_new_files_and_does_not_emit_stdout(tmp_path, monkeypatch, failed_writer):
    def fail(*args, **kwargs):
        args[2].write("partial")
        raise OSError("simulated disk failure")

    monkeypatch.setattr(report, failed_writer, fail)
    paths = [tmp_path / name for name in ("main.tsv", "context.tsv", "report.html")]
    stdout = io.StringIO()
    with pytest.raises(report.ReportError, match="simulated disk failure") as caught:
        report.write_reports("GRCh38", (support(),), *paths[:2], stdout, html_output_path=paths[2])
    failed_index = ["write_tsv", "write_context_tsv", "write_html"].index(failed_writer)
    assert str(paths[failed_index]) in str(caught.value)
    assert not any(path.exists() for path in paths)
    assert stdout.getvalue() == ""


def test_html_failure_precedes_stdout_and_removes_context(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise OSError("HTML destination failed")

    monkeypatch.setattr(report, "write_html", fail)
    context, html = tmp_path / "context.tsv", tmp_path / "report.html"
    stdout = io.StringIO()
    with pytest.raises(report.ReportError, match="HTML destination failed"):
        report.write_reports("GRCh38", (support(),), None, context, stdout, html_output_path=html)
    assert stdout.getvalue() == ""
    assert not context.exists() and not html.exists()


def test_flush_failure_removes_outputs(tmp_path, monkeypatch):
    original_open = Path.open

    class FailingFlush:
        def __init__(self, stream):
            self.stream = stream
        def __enter__(self):
            return self
        def __exit__(self, *args):
            self.stream.close()
        def write(self, value):
            return self.stream.write(value)
        def flush(self):
            raise OSError("simulated flush failure")

    html = tmp_path / "report.html"
    context = tmp_path / "context.tsv"
    def open_path(path, *args, **kwargs):
        stream = original_open(path, *args, **kwargs)
        return FailingFlush(stream) if path == html else stream
    monkeypatch.setattr(Path, "open", open_path)
    with pytest.raises(report.ReportError, match="simulated flush failure"):
        report.write_reports("GRCh38", (support(),), None, context, io.StringIO(), html_output_path=html)
    assert not html.exists() and not context.exists()
