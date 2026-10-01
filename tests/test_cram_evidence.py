"""Real CRAM/reference validation and equivalence with BAM evidence."""
from dataclasses import replace
from pathlib import Path

import pysam
import pytest

from alignment_fixtures import CATALOG, bam, cli_inputs, cram, fasta, record
from sjsift.alignment_evidence import AlignmentError, extract_evidence
from test_alignment_cli import invoke
from test_reference_alignments import reference


def test_equivalent_cram_preserves_bam_evidence_and_sampling(tmp_path):
    variant = replace(CATALOG.variants[0], reference_junctions=(reference(),))
    catalog = replace(CATALOG, variants=(variant,))
    reads = [record(f"read-{i}", cigar="4S20M100N20M100N20M3S", flag=flag,
                    tags=tags, mapq=0) for i, (flag, tags) in enumerate([
        (65, {"NH": 1, "XS": "+"}), (129, {"NH": 2}), (256, {"NH": 1}),
        (512, {"NH": 1}), (1024, {"NH": 1}), (2048, {"NH": 1}),
        (16, {"NH": 1, "jM": [21, 22]}), (0, {"TS": "-"}),
        (0, {"TS": "+", "XS": "-"}), (0, {}),
    ])]
    reads += [record("identical", tags={"NH": 1}) for _ in range(2)]
    source = bam(tmp_path, reads)
    ref = fasta(tmp_path)
    path = cram(tmp_path, source, ref)
    expected = extract_evidence(catalog, source, limit=3)
    actual = extract_evidence(catalog, path, reference_path=ref, limit=3)
    assert expected.groups == actual.groups
    assert actual.provenance["Format"] == "CRAM"
    assert actual.provenance["Reference FASTA"] == "reference.fa"
    assert "M5" in actual.provenance["Reference assurance"]


@pytest.mark.parametrize("mode", ["external", "embedded", "no_ref"])
def test_every_cram_requires_explicit_reference_even_when_self_contained(tmp_path, mode):
    ref = fasta(tmp_path)
    options = {"external": None, "embedded": ["embed_ref=1"], "no_ref": ["no_ref=1"]}[mode]
    path = cram(tmp_path, bam(tmp_path, [record()]), ref, options)
    with pytest.raises(AlignmentError, match="requires --reference"):
        extract_evidence(CATALOG, path)


@pytest.mark.parametrize("problem", [
    "missing_reference", "missing_fai", "invalid_fai", "wrong_bases", "short_reference",
    "missing_contig", "truncated_reference", "missing_crai", "invalid_crai", "truncated_cram", "missing_cram_eof", "corrupt_cram",
])
def test_bad_cram_inputs_create_no_outputs_or_implicit_indexes(tmp_path, problem):
    ref = fasta(tmp_path)
    path = cram(tmp_path, bam(tmp_path, [record()]), ref)
    if problem == "missing_reference": ref.unlink()
    if problem == "missing_fai": Path(str(ref) + ".fai").unlink()
    if problem == "invalid_fai": Path(str(ref) + ".fai").write_text("broken\n")
    if problem == "wrong_bases": ref = fasta(tmp_path, "C" * 2000, name="wrong.fa")
    if problem == "short_reference": ref = fasta(tmp_path, "A" * 100, name="short.fa")
    if problem == "missing_contig": ref = fasta(tmp_path, contig="1", name="other.fa")
    if problem == "truncated_reference": ref.write_text(">chr1\nAAAA\n")
    if problem == "missing_crai": Path(str(path) + ".crai").unlink()
    if problem == "invalid_crai": Path(str(path) + ".crai").write_bytes(b"broken")
    if problem == "truncated_cram": path.write_bytes(path.read_bytes()[:-50])
    if problem == "missing_cram_eof": path.write_bytes(path.read_bytes()[:-38])
    if problem == "corrupt_cram":
        content = bytearray(path.read_bytes())
        content[-70] ^= 0xff
        path.write_bytes(content)
    args = cli_inputs(tmp_path)
    before = {p: p.read_bytes() for p in tmp_path.iterdir()}
    outputs = [tmp_path / f"out.{ext}" for ext in ("tsv", "context.tsv", "html")]
    result = invoke(*args, "--alignments", path, "--reference", ref,
                    "--output", outputs[0], "--context-output", outputs[1], "--html-output", outputs[2])
    assert result.returncode == 2, result.stderr
    assert "cannot extract alignment evidence" in result.stderr
    assert not result.stdout and "Traceback" not in result.stderr
    assert not any(p.exists() for p in outputs)
    assert {p: p.read_bytes() for p in tmp_path.iterdir()} == before


@pytest.mark.parametrize("index_name", ["reads.crai", "custom.index"])
def test_cram_cli_explicit_and_conventional_indexes_preserve_tsv_and_privacy(tmp_path, index_name):
    args = cli_inputs(tmp_path)
    ref = fasta(tmp_path)
    path = cram(tmp_path, bam(tmp_path, [record("raw-read", tags={"NH": 1})]), ref)
    index = tmp_path / index_name
    Path(str(path) + ".crai").rename(index)
    options = ["--alignment-index", index] if index_name == "custom.index" else []
    output = tmp_path / "report.html"
    baseline = invoke(*args)
    result = invoke(*args, "--alignments", path, "--reference", ref, "--html-output", output, *options)
    assert result.returncode == 0, result.stderr
    assert result.stdout == baseline.stdout
    html = output.read_text()
    assert "raw-read" in html and "CRAM" in html and "reference.fa" in html
    assert "ACGTACGTACGT" not in html and "IIIIIIIIIIII" not in html and "AAAAAAAAAAAA" not in html
    assert str(tmp_path) not in html
    assert "STAR sample identity cannot be verified" in html


def test_bgzip_reference_requires_existing_gzi_and_supports_real_decoding(tmp_path):
    ref = fasta(tmp_path)
    path = cram(tmp_path, bam(tmp_path, [record()]), ref)
    compressed = tmp_path / "reference.fa.gz"
    pysam.tabix_compress(str(ref), str(compressed))
    pysam.faidx(str(compressed))
    assert sum(extract_evidence(CATALOG, path, reference_path=compressed).groups["v0-j0"].eligible.values()) == 1
    Path(str(compressed) + ".gzi").unlink()
    with pytest.raises(AlignmentError, match="gzi"):
        extract_evidence(CATALOG, path, reference_path=compressed)
    assert not Path(str(compressed) + ".gzi").exists()


def test_missing_m5_discloses_limited_assurance_and_still_requires_local_reference(tmp_path):
    ref = fasta(tmp_path)
    source = bam(tmp_path, [record()])
    path = tmp_path / "no-reference.cram"
    with pysam.AlignmentFile(str(source), "rb") as stream:
        with pysam.AlignmentFile(str(path), "wc", header=stream.header, format_options=["no_ref=1"]) as output:
            for read in stream:
                output.write(read)
    pysam.index(str(path))
    result = extract_evidence(CATALOG, path, reference_path=ref)
    assert "0 of 1 header M5" in result.provenance["Reference assurance"]
    assert "Missing M5 limits" in result.provenance["Reference assurance"]


def test_remote_reference_environment_is_isolated_and_restored(tmp_path, monkeypatch):
    import os
    ref = fasta(tmp_path)
    path = cram(tmp_path, bam(tmp_path, [record()]), ref)
    settings = {"REF_PATH": "https://example.invalid/%s", "REF_CACHE": "inherited-cache/%s"}
    for key, value in settings.items(): monkeypatch.setenv(key, value)
    extract_evidence(CATALOG, path, reference_path=ref)
    assert {key: os.environ[key] for key in settings} == settings
    wrong = fasta(tmp_path, "C" * 2000, name="wrong.fa")
    with pytest.raises(AlignmentError, match="M5"):
        extract_evidence(CATALOG, path, reference_path=wrong)
    assert {key: os.environ[key] for key in settings} == settings


def test_cram_preserves_existing_output_and_removes_new_destinations(tmp_path):
    ref = fasta(tmp_path)
    path = cram(tmp_path, bam(tmp_path, [record()]), ref)
    keep = tmp_path / "keep.html"
    keep.write_text("existing report")
    new = tmp_path / "new.tsv"
    result = invoke(*cli_inputs(tmp_path), "--alignments", path, "--reference", ref,
                    "--html-output", keep, "--output", new)
    assert result.returncode == 2 and not result.stdout
    assert keep.read_text() == "existing report" and not new.exists()


def test_cram_never_contacts_inherited_refget_with_real_network_detector(tmp_path, monkeypatch):
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    import subprocess
    import sys
    from threading import Thread

    requests = []

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            requests.append(self.path)
            self.send_response(200)
            self.send_header("Content-Length", "2000")
            self.end_headers()
            self.wfile.write(b"A" * 2000)

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        original = fasta(tmp_path, name="original.fa")
        path = cram(tmp_path, bam(tmp_path, [record()]), original)
        # Replace the writer's local UR metadata with a remote URI, too.
        with pysam.AlignmentFile(str(path), "rc") as stream:
            header = str(stream.header)
        url = f"http://127.0.0.1:{server.server_port}/%s"
        header_path = tmp_path / "remote-header.sam"
        header_path.write_text(header.replace(str(original.resolve()), url))
        replaced = path
        pysam.reheader("-i", "-P", str(header_path), str(replaced))
        pysam.index(str(replaced))
        original.unlink()
        Path(str(original) + ".fai").unlink()
        monkeypatch.setenv("REF_PATH", url)
        monkeypatch.setenv("REF_CACHE", str(tmp_path / "cache" / "%s"))
        # Positive control: unguarded HTSlib really reaches this HTTP detector.
        code = "import pysam,sys; f=pysam.AlignmentFile(sys.argv[1],'rc'); assert next(f).query_sequence"
        probe = subprocess.run([sys.executable, "-c", code, str(replaced)], capture_output=True, text=True, timeout=15)
        assert probe.returncode == 0, probe.stderr
        assert requests, "positive control must detect a real HTSlib network request"
        requests.clear()
        monkeypatch.setenv("REF_CACHE", str(tmp_path / "unused-cache" / "%s"))
        ref = fasta(tmp_path)
        args = cli_inputs(tmp_path)
        good = invoke(*args, "--alignments", replaced, "--reference", ref, "--html-output", tmp_path / "good.html")
        assert good.returncode == 0, good.stderr
        wrong = fasta(tmp_path, "C" * 2000, name="wrong.fa")
        bad = invoke(*args, "--alignments", replaced, "--reference", wrong, "--html-output", tmp_path / "bad.html")
        missing = invoke(*args, "--alignments", replaced, "--html-output", tmp_path / "missing.html")
        assert bad.returncode == missing.returncode == 2
        assert not requests, "sjsift must never attempt remote reference resolution"
        assert not (tmp_path / "bad.html").exists() and not (tmp_path / "missing.html").exists()
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
