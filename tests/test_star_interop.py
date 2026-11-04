"""Installed CLI interoperability with unmodified, pinned STAR output."""
import csv
from html.parser import HTMLParser
from pathlib import Path
import subprocess
import sysconfig

import pysam
import pytest


FIXTURES = Path(__file__).parent / "fixtures" / "star"


class Report(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.overview = []
        self.reads = {}
        self.variant = None
        self.counts = {}
        self.index = None
        self.count_depth = 0
        self.count_text = []
        self.labels = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "article" and attrs.get("id", "").startswith("variant-"):
            self.variant = attrs["id"]
            self.reads[self.variant] = []
        if tag == "span" and (self.count_depth or attrs.get("class") == "counts"):
            self.count_depth += 1
        if tag == "summary" and attrs.get("class") == "read-hit":
            self.labels.append(attrs["title"])
        if tag == "tr" and "data-index" in attrs:
            self.index = attrs["data-index"]
            self.counts[self.index] = []
            self.overview.append((attrs["data-index"], attrs["data-unique"], attrs["data-total"]))
        if tag == "li" and "read-row" in attrs.get("class", "").split():
            self.reads[self.variant].append((attrs["data-mapping"], attrs["data-support"]))

    def handle_data(self, text):
        if self.count_depth and self.index is not None:
            self.count_text.append(text)

    def handle_endtag(self, tag):
        if tag == "span" and self.count_depth:
            self.count_depth -= 1
            if not self.count_depth and self.index is not None:
                self.counts[self.index].append("".join(self.count_text))
                self.count_text = []
        if tag == "tr":
            self.index = None


@pytest.mark.parametrize("format", ["bam", "cram"])
def test_star_output_through_installed_cli_preserves_support_and_memberships(tmp_path, format):
    alignment = FIXTURES / "alignments.bam"
    reference = FIXTURES / "reference.fa"
    options = []
    if format == "cram":
        converted = tmp_path / "alignments.cram"
        with pysam.AlignmentFile(str(alignment), "rb") as source:
            with pysam.AlignmentFile(str(converted), "wc", header=source.header,
                                     reference_filename=str(reference)) as destination:
                for record in source:
                    destination.write(record)
        pysam.index(str(converted))
        alignment = converted
        options = ["--reference", str(reference)]
    # The inventory is specified from read construction and SAM conventions,
    # never obtained by asking sjsift to interpret its own output.
    with pysam.AlignmentFile(str(alignment), reference_filename=str(reference)) as source:
        actual = [list(map(str, [
            r.query_name, r.reference_name, r.reference_start + 1, r.flag, r.mapping_quality,
            r.cigarstring, r.next_reference_start + 1, r.template_length, r.get_tag("NH"), r.get_tag("HI"),
            ",".join(map(str, r.get_tag("jM"))), ",".join(map(str, r.get_tag("jI"))),
        ])) for r in source]
    with (FIXTURES / "expected-alignments.tsv").open() as stream:
        inventory = list(csv.reader(stream, delimiter="\t"))
    assert actual == inventory[1:]
    main = tmp_path / "main.tsv"
    context = tmp_path / "context.tsv"
    html = tmp_path / "report.html"
    command = Path(sysconfig.get_path("scripts")) / "sjsift"
    result = subprocess.run([
        str(command), "--sample-name", "STAR interoperability", "--junctions", str(FIXTURES / "SJ.out.tab"),
        "--definitions", str(FIXTURES / "catalog.toml"), "--alignment", str(alignment),
        "--output", str(main), "--context-output", str(context), "--html-output", str(html), *options,
    ], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert result.stdout == "" and result.stderr == ""
    assert main.read_bytes() == (FIXTURES / "expected-main.tsv").read_bytes()
    assert context.read_bytes() == (FIXTURES / "expected-context.tsv").read_bytes()
    text = html.read_text()
    parsed = Report(text)
    assert sorted(parsed.overview) == [("0", "1", "1"), ("1", "0", "0"), ("2", "2", "2"), ("3", "0", "1"), ("4", "0", "1")]
    assert {key: sorted(value) for key, value in parsed.reads.items()} == {
        "variant-0": [("unique", "0"), ("unique", "0 1")],
        "variant-1": [],
        "variant-2": [("unique", "0 1"), ("unique", "0 1"), ("unique", "1"), ("unique", "1")],
        "variant-3": [("multimapping", "0 1"), ("multimapping", "1")],
        "variant-4": [("multimapping", "0 1"), ("multimapping", "1")],
    }
    assert parsed.counts == {
        "0": ["1 unique + 0 multi"] * 2, "1": ["0 unique + 0 multi"],
        "2": ["2 unique + 0 multi"] * 2,
        "3": ["0 unique + 1 multi"] * 2, "4": ["0 unique + 1 multi"] * 2,
    }
    assert len(parsed.labels) == 10
    assert all("strand evidence agrees" in label for label in parsed.labels)
    assert any("read 1, reverse" in label for label in parsed.labels)
    assert any("secondary" in label and "NH: 2" in label for label in parsed.labels)
    assert "connect-src 'none'" in text
    assert "<script src=" not in text and 'rel="stylesheet"' not in text
