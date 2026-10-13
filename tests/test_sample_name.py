"""Sample identity through real BAM/CRAM inputs and the public CLI."""

from html import escape

import pytest

from alignment_fixtures import bam, cli_inputs, cram, fasta, record
from test_alignment_cli import invoke


def alignment(tmp_path, kind, read_groups, reads):
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}], "RG": read_groups}
    path = bam(tmp_path, reads, header=header)
    if kind == "CRAM":
        reference = fasta(tmp_path)
        return ["--alignments", cram(tmp_path, path, reference), "--reference", reference]
    return ["--alignments", path]


@pytest.mark.parametrize("html", [False, True])
@pytest.mark.parametrize("name", [None, "", "   "])
def test_sample_name_required_without_alignments_before_any_output(tmp_path, html, name):
    args = cli_inputs(tmp_path, sample_name=name)
    tsv, report = tmp_path / "out.tsv", tmp_path / "out.html"
    result = invoke(*args, "--output", tsv, *(["--html-output", report] if html else []))
    assert result.returncode == 2 and "sample" in result.stderr
    assert not result.stdout and not tsv.exists() and not report.exists()


@pytest.mark.parametrize("kind", ["BAM", "CRAM"])
def test_infers_one_distinct_sample_across_read_groups(tmp_path, kind):
    name = 'Patient <A> & "B"'
    options = alignment(tmp_path, kind, [{"ID": "a", "SM": name}, {"ID": "b", "SM": name}],
                        [record("lane-a", tags={"RG": "a"}), record("lane-b", tags={"RG": "b"})])
    report = tmp_path / "out.html"
    result = invoke(*cli_inputs(tmp_path, sample_name=None), *options, "--html-output", report)
    assert result.returncode == 0, result.stderr
    html = report.read_text()
    assert f'<title>{escape(name)} · SJSIFT junction report</title>' in html
    assert f'<p class="sample-name">{escape(name)}</p>' in html
    assert "lane-a" in html and "lane-b" in html
    header = html.split('<header>')[1].split('</header>')[0]
    assert 'Inputs' not in header and 'sample.SJ.out.tab' not in header
    footer = html.split('<footer>')[1]
    assert 'Inputs and configuration' in footer and 'sample.SJ.out.tab' in footer


@pytest.mark.parametrize("kind", ["BAM", "CRAM"])
@pytest.mark.parametrize("name", [None, "patient"])
def test_missing_alignment_sample_metadata_requires_explicit_name(tmp_path, kind, name):
    options = alignment(tmp_path, kind, [], [record()])
    report = tmp_path / "out.html"
    result = invoke(*cli_inputs(tmp_path, sample_name=name), *options, "--html-output", report)
    if name:
        assert result.returncode == 0, result.stderr
        assert '<p class="sample-name">patient</p>' in report.read_text()
    else:
        assert result.returncode == 2 and '--sample-name is required' in result.stderr
        assert not report.exists() and not result.stdout


@pytest.mark.parametrize("kind", ["BAM", "CRAM"])
@pytest.mark.parametrize("name", [None, "missing", "chosen"])
def test_multiple_samples_require_valid_selection_and_only_embed_that_sample(tmp_path, kind, name):
    options = alignment(tmp_path, kind, [
        {"ID": "a", "SM": "chosen"}, {"ID": "b", "SM": "other"}, {"ID": "c", "SM": "chosen"}, {"ID": "unknown"},
    ], [record("chosen-a", tags={"RG": "a"}), record("other-read", tags={"RG": "b"}),
        record("chosen-c", tags={"RG": "c"}), record("no-sample", tags={"RG": "unknown"}), record("no-rg")])
    args = cli_inputs(tmp_path, sample_name=name)
    report, tsv = tmp_path / "out.html", tmp_path / "out.tsv"
    result = invoke(*args, *options, "--html-output", report, "--output", tsv)
    if name != "chosen":
        assert result.returncode == 2 and '--sample-name' in result.stderr
        assert not report.exists() and not tsv.exists() and not result.stdout
        return
    assert result.returncode == 0, result.stderr
    html = report.read_text()
    assert '<p class="sample-name">chosen</p>' in html
    assert html.count('<li class="read-row') == 2
    assert 'chosen-a' in html and 'chosen-c' in html
    assert 'other-read' not in html and 'no-sample' not in html and 'no-rg' not in html
    baseline = invoke(*args)
    assert baseline.returncode == 0 and tsv.read_text() == baseline.stdout


@pytest.mark.parametrize("kind", ["BAM", "CRAM"])
def test_explicit_name_must_match_single_alignment_sample(tmp_path, kind):
    options = alignment(tmp_path, kind, [{"ID": "a", "SM": "actual"}], [record(tags={"RG": "a"})])
    report = tmp_path / "out.html"
    result = invoke(*cli_inputs(tmp_path, sample_name="wrong"), *options, "--html-output", report)
    assert result.returncode == 2 and 'does not exist in alignment header' in result.stderr
    assert not report.exists() and not result.stdout
