"""Alignment input, output protection, and portable report behavior through CLI."""
from pathlib import Path
import subprocess
import sys

import pytest

from alignment_fixtures import bam, cli_inputs, record


def invoke(*args):
    return subprocess.run([sys.executable, "-m", "sjsift", *map(str, args)], capture_output=True, text=True)


def test_bam_report_embeds_previews_and_preserves_exact_star_tsv(tmp_path):
    args = cli_inputs(tmp_path)
    reads = bam(tmp_path, [record("raw-name", tags={"NH": 1, "XS": "+"}, flag=65),
                           record("raw-name", flag=129, tags={"NH": 2})])
    html = tmp_path / "report.html"
    original = invoke(*args)
    result = invoke(*args, "--alignments", reads, "--html-output", html, "--alignment-limit", "1")
    assert result.returncode == 0, result.stderr
    assert result.stdout == original.stdout
    text = html.read_text()
    for expected in ("raw-name", "read 1", "read 2", "Alignment records", "Eligible: 2", "Embedded: 2", "20M100N20M", "strand unverified"):
        assert expected in text
    assert str(tmp_path) not in text
    assert "ACGTACGTACGT" not in text and "IIIIIIIIIIII" not in text
    assert "STAR support is separate" in text


@pytest.mark.parametrize("options,message", [
    (["--alignments", "missing.bam"], "requires --html-output"),
    (["--alignment-index", "missing.bai"], "require --alignments"),
    (["--alignment-limit", "1"], "require --alignments"),
    (["--alignment-limit", "0"], "positive integer"),
    (["--alignment-limit", "-1"], "positive integer"),
])
def test_alignment_option_dependencies_fail_before_output(tmp_path, options, message):
    result = invoke(*cli_inputs(tmp_path), *options)
    assert result.returncode == 2 and message in result.stderr
    assert not result.stdout


@pytest.mark.parametrize("problem", ["missing_index", "corrupt_index", "truncated", "sort", "contig", "bounds", "samples", "sam", "missing_file"])
def test_invalid_alignment_inputs_stop_before_any_report_output(tmp_path, problem):
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}]}
    if problem == "sort": header["HD"]["SO"] = "queryname"
    if problem == "contig": header["SQ"][0]["SN"] = "1"
    if problem == "bounds": header["SQ"][0]["LN"] = 150
    if problem == "samples": header["RG"] = [{"ID": "a", "SM": "one"}, {"ID": "b", "SM": "two"}]
    reads = bam(tmp_path, [record()], header=header)
    if problem == "missing_index": Path(str(reads) + ".bai").unlink()
    if problem == "corrupt_index": Path(str(reads) + ".bai").write_bytes(b"broken")
    if problem == "truncated": reads.write_bytes(reads.read_bytes()[:-40])
    if problem == "sam": reads.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n")
    if problem == "missing_file": reads.unlink()
    outputs = [tmp_path / f"out.{ext}" for ext in ("tsv", "context.tsv", "html")]
    result = invoke(*cli_inputs(tmp_path), "--alignments", reads, "--output", outputs[0], "--context-output", outputs[1], "--html-output", outputs[2])
    assert result.returncode == 2, result.stderr
    assert not result.stdout and "Traceback" not in result.stderr
    assert not any(p.exists() for p in outputs)


def test_explicit_index_and_limits_preserve_tsv_and_deterministic_safe_html(tmp_path):
    args = cli_inputs(tmp_path)
    hostile = '</script><img/src=x/onerror=alert(1)>'
    reads = bam(tmp_path, [record(hostile, tags={"NH": 1}) for _ in range(12)])
    index = tmp_path / "custom.index"
    Path(str(reads) + ".bai").rename(index)
    outputs = [tmp_path / f"report{i}.html" for i in range(3)]
    results = [invoke(*args, "--alignments", reads, "--alignment-index", index, "--alignment-limit", cap, "--html-output", output)
               for cap, output in zip(("2", "2", "20"), outputs)]
    assert all(r.returncode == 0 for r in results)
    assert len({r.stdout for r in results}) == 1
    assert outputs[0].read_bytes() == outputs[1].read_bytes()
    html = outputs[2].read_text()
    assert "Eligible: 12" in html and "Embedded: 12" in html and "Currently shown: 10" in html
    assert hostile not in html and "&lt;/script&gt;" in html
    assert html.count('<li class="read-row') == 12
    assert html.count(' hidden>') >= 2
    assert "No eligible alignment records (0)" in invoke_empty_report(tmp_path, args)


def invoke_empty_report(tmp_path, args):
    reads = bam(tmp_path, [], name="empty.bam")
    output = tmp_path / "empty.html"
    result = invoke(*args, "--alignments", reads, "--html-output", output)
    assert result.returncode == 0
    return output.read_text()


def test_alignment_report_respects_existing_output_and_cleans_new_files(tmp_path):
    args = cli_inputs(tmp_path)
    reads = bam(tmp_path, [record()])
    output = tmp_path / "keep.html"
    output.write_text("keep")
    new = tmp_path / "new.tsv"
    result = invoke(*args, "--alignments", reads, "--html-output", output, "--output", new)
    assert result.returncode == 2 and not result.stdout
    assert output.read_text() == "keep" and not new.exists()


@pytest.mark.parametrize("assembly,success", [("synthetic", True), ("different", False)])
def test_available_alignment_assembly_labels_must_agree_without_claiming_star_identity(tmp_path, assembly, success):
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000, "AS": assembly}], "RG": [{"ID": "lane1", "SM": "sample"}]}
    path = bam(tmp_path, [record(tags={"RG": "lane1"})], header=header)
    output = tmp_path / "report.html"
    result = invoke(*cli_inputs(tmp_path), "--alignments", path, "--html-output", output)
    if success:
        assert result.returncode == 0
        assert "STAR sample identity cannot be verified" in output.read_text()
        assert "sequence equivalence" in output.read_text()
    else:
        assert result.returncode == 2 and "assembly metadata disagrees" in result.stderr
        assert not output.exists() and not result.stdout
