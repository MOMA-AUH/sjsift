"""Run provenance through the installed console entry point."""

import json
import gzip
import hashlib
import os
from pathlib import Path
import signal
import subprocess
import sysconfig

import pysam
import pytest

from alignment_fixtures import bam, cli_inputs, cram, fasta, record
from test_cli import run_console_script


def invoke(*args):
    return run_console_script(*map(str, args))


def test_tsv_only_provenance_identifies_inputs_and_is_deterministic(tmp_path):
    args = cli_inputs(tmp_path, sample_name="patient α")
    baseline = invoke(*args)
    documents = []
    for name in ("first.json", "second.json"):
        output = tmp_path / name
        result = invoke(*args, "--provenance-output", output)
        assert result.returncode == 0, result.stderr
        assert result.stdout == baseline.stdout and result.stderr == baseline.stderr == ""
        documents.append(output.read_bytes())
    assert documents[0] == documents[1]
    data = json.loads(documents[0])
    assert data["schema_version"] == 1
    assert next(iter(data)) == "schema_version"
    assert next(iter(data["catalog"])) == "schema_version"
    assert data["sample"] == {"name": "patient α", "source": "explicit"}
    assert data["software"]["name"] == "sjsift"
    assert data["software"]["version"] == invoke("--version").stdout.strip().split()[1]
    assert data["genome_assembly"] == "synthetic"
    assert data["catalog"]["schema_version"] == 3
    assert data["catalog"]["basename"] == "catalog.toml"
    assert len(data["catalog"]["sha256"]) == 64
    assert data["star_input"]["basename"] == "sample.SJ.out.tab"
    assert data["star_input"]["compression"] == "none"
    assert data["star_input"]["sha256"] == "bb5f835a08c52516b17efb1d9e03b5f509801943b7b4521961beb8335d84d34e"
    assert data["star_input"]["fingerprint_scope"] == "exact_file_bytes"
    assert data["warnings"] == []
    assert data["alignment"] == {"status": "not_requested"}
    assert b"assurance" not in documents[0]
    assert str(tmp_path).encode() not in documents[0]


@pytest.mark.parametrize("content", [b"", b"other\t101\t200\t1\t1\t0\t8\t2\t20\n"])
def test_compressed_fingerprints_and_warnings_preserve_all_existing_reports(tmp_path, content):
    args = cli_inputs(tmp_path)
    catalog = tmp_path / "catalog.toml"
    # A harmless comment is still part of the exact catalog identity.
    catalog.write_bytes(catalog.read_bytes() + b"\n# tracked catalog revision\n")
    sj = tmp_path / "sample.SJ.out.tab"
    compressed = gzip.compress(content, mtime=0)
    sj.write_bytes(compressed)  # Detection must follow bytes, not the suffix.
    outputs = []
    for provenance in (False, True):
        directory = tmp_path / str(provenance)
        directory.mkdir()
        main, context, html, sidecar = [directory / name for name in ("main.tsv", "context.tsv", "out.html", "run.json")]
        result = invoke(*args, "--output", main, "--context-output", context, "--html-output", html,
                        *(["--provenance-output", sidecar] if provenance else []))
        assert result.returncode == 0, result.stderr
        outputs.append([p.read_bytes() for p in (main, context, html)])
    assert outputs[0] == outputs[1]
    data = json.loads(sidecar.read_text())
    assert data["catalog"]["sha256"] == hashlib.sha256(catalog.read_bytes()).hexdigest()
    assert data["star_input"]["sha256"] == hashlib.sha256(compressed).hexdigest()
    assert data["star_input"]["sha256"] != hashlib.sha256(content).hexdigest()
    assert data["star_input"]["compression"] == "gzip"
    assert data["star_input"]["fingerprint_scope"] == "exact_file_bytes"
    assert len(data["warnings"]) == 1
    warning = data["warnings"][0]
    assert warning["code"] == "NO_CATALOG_CHROMOSOMES"
    assert warning["message"] in result.stderr
    assert b"No catalog chromosome identifiers were found" in outputs[1][2]


@pytest.mark.parametrize("kind", ["BAM", "CRAM"])
@pytest.mark.parametrize("eligible", [False, True])
def test_alignment_provenance_reports_checks_without_record_or_header_payloads(tmp_path, kind, eligible):
    args = cli_inputs(tmp_path, sample_name=None)
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000, "AS": "synthetic"}],
              "RG": [{"ID": "lane", "SM": "inferred patient"}],
              "PG": [{"ID": "aligner", "CL": "SECRET_COMMAND /private/raw.fastq"}],
              "CO": ["SECRET_COMMENT"]}
    reads = bam(tmp_path, [record("SECRET_READ_NAME", tags={"RG": "lane", "NH": 1})] if eligible else [], header=header)
    options = []
    if kind == "CRAM":
        reference = fasta(tmp_path)
        reads = cram(tmp_path, reads, reference)
        options = ["--reference", reference]
    index = tmp_path / "custom.index"
    reads.with_name(reads.name + (".bai" if kind == "BAM" else ".crai")).rename(index)
    documents = []
    for number in (1, 2):
        sidecar, html = tmp_path / f"run{number}.json", tmp_path / f"run{number}.html"
        result = invoke(*args, "--alignment", reads, "--alignment-index", index, *options,
                        "--html-output", html, "--provenance-output", sidecar)
        assert result.returncode == 0, result.stderr
        documents.append(sidecar.read_bytes())
    assert documents[0] == documents[1]
    assert b"assurance" not in documents[0]
    data = json.loads(documents[0])
    assert next(iter(data)) == "schema_version"
    assert data["sample"] == {"name": "inferred patient", "source": "alignment_header"}
    alignment = data["alignment"]
    assert alignment["status"] == "requested"
    assert alignment["eligible_junction_memberships"] == int(eligible)
    assert alignment["format"] == kind
    assert alignment["basename"] == reads.name and alignment["index_basename"] == "custom.index"
    assert alignment["sample_selection"] == {"policy": "all_read_groups", "header_sample_count": 1,
                                             "selected_sample_in_header": True,
                                             "header_read_group_count": 1, "read_groups_without_sample": 0}
    assert set(alignment["dependencies"]) == {"pysam", "htslib"}
    assert all(alignment["dependencies"].values())
    assert alignment["assembly_checks"] == {"required_contigs": 1, "labeled_contigs": 1,
                                            "labels_match_catalog": True}
    if kind == "BAM":
        assert alignment["reference"] is None
    else:
        reference = alignment["reference"]
        assert reference["basename"] == "reference.fa"
        assert reference["index_basenames"] == ["reference.fa.fai"]
        assert reference["checks"] == {"contig_names": True, "contig_lengths": True, "full_sequences_read": True,
                                      "header_m5_checked": 1, "header_contigs": 1, "decoded_slice_checks": True}
        assert reference["access"] == "explicit_local_only"
    for private in (str(tmp_path), "SECRET", "ACGTACGTACGT", "IIIIIIIIIIII", "AAAAAAAAAAAA", "/private/"):
        assert private.encode() not in documents[0]
    assert "sha256" not in json.dumps(alignment)


@pytest.mark.parametrize("samples", [[], ["chosen"], ["chosen", "other"]])
def test_explicit_sample_selection_and_missing_assembly_labels(tmp_path, samples):
    args = cli_inputs(tmp_path, sample_name="chosen")
    groups = [{"ID": name, "SM": name} for name in samples]
    reads = [record("unassigned")]
    reads.extend(record(name, tags={"RG": name}) for name in samples)
    path = bam(tmp_path, reads, header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}], "RG": groups})
    sidecar = tmp_path / "run.json"
    result = invoke(*args, "--alignment", path, "--html-output", tmp_path / "report.html", "--provenance-output", sidecar)
    assert result.returncode == 0, result.stderr
    data = json.loads(sidecar.read_text())
    assert data["sample"] == {"name": "chosen", "source": "explicit"}
    selection = data["alignment"]["sample_selection"]
    assert selection["header_sample_count"] == len(samples)
    assert selection["selected_sample_in_header"] == bool(samples)
    assert selection["header_read_group_count"] == len(samples)
    assert selection["read_groups_without_sample"] == 0
    assert selection["policy"] == ("selected_sample_read_groups" if len(samples) > 1 else "all_read_groups")
    assert data["alignment"]["eligible_junction_memberships"] == (1 if len(samples) > 1 else 1 + len(samples))
    assert data["alignment"]["assembly_checks"] == {"required_contigs": 1, "labeled_contigs": 0,
                                                   "labels_match_catalog": None}


def test_cram_provenance_records_compressed_reference_indexes_and_missing_m5(tmp_path):
    args = cli_inputs(tmp_path)
    source = bam(tmp_path, [record()])
    reference = fasta(tmp_path)
    path = tmp_path / "no-reference.cram"
    with pysam.AlignmentFile(str(source), "rb") as stream:
        with pysam.AlignmentFile(str(path), "wc", header=stream.header, format_options=["no_ref=1"]) as output:
            for read in stream:
                output.write(read)
    pysam.index(str(path))
    compressed = tmp_path / "reference.fa.gz"
    pysam.tabix_compress(str(reference), str(compressed))
    pysam.faidx(str(compressed))
    sidecar = tmp_path / "run.json"
    result = invoke(*args, "--alignment", path, "--reference", compressed,
                    "--html-output", tmp_path / "report.html", "--provenance-output", sidecar)
    assert result.returncode == 0, result.stderr
    metadata = json.loads(sidecar.read_text())["alignment"]["reference"]
    assert metadata["index_basenames"] == ["reference.fa.gz.fai", "reference.fa.gz.gzi"]
    assert metadata["checks"]["header_m5_checked"] == 0
    assert metadata["checks"]["header_contigs"] == 1


@pytest.mark.parametrize("failure", ["existing", "missing_parent", "input_alias", "directory"])
def test_provenance_creation_failure_protects_existing_files_and_cleans_reports(tmp_path, failure):
    args = cli_inputs(tmp_path)
    sidecar = tmp_path / "run.json"
    if failure == "existing":
        sidecar.write_bytes(b"keep me\n")
    elif failure == "missing_parent":
        sidecar = tmp_path / "missing" / "run.json"
    elif failure == "input_alias":
        sidecar = tmp_path / "catalog.toml"
    else:
        sidecar.mkdir()
    original = sidecar.read_bytes() if sidecar.is_file() else None
    outputs = [tmp_path / name for name in ("main.tsv", "context.tsv", "report.html")]
    result = invoke(*args, "--output", outputs[0], "--context-output", outputs[1], "--html-output", outputs[2],
                    "--provenance-output", sidecar)
    assert result.returncode == 2 and not result.stdout
    assert "Traceback" not in result.stderr and "internal error" not in result.stderr
    assert not any(path.exists() for path in outputs)
    if original is not None:
        assert sidecar.read_bytes() == original


@pytest.mark.parametrize("other", ["--output", "--context-output", "--html-output"])
def test_provenance_cannot_alias_another_output(tmp_path, other):
    args = cli_inputs(tmp_path)
    output = tmp_path / "result"
    directory = tmp_path / "subdir"
    directory.mkdir()
    result = invoke(*args, other, output, "--provenance-output", directory / ".." / "result")
    assert result.returncode == 2 and not result.stdout
    assert "output paths must differ" in result.stderr and not output.exists()


def test_provenance_write_failure_removes_partial_json_and_tsvs(tmp_path):
    resource = pytest.importorskip("resource")
    args = cli_inputs(tmp_path)
    outputs = [tmp_path / name for name in ("main.tsv", "context.tsv", "run.json")]

    def limit_file_size():
        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_FSIZE, (512, 512))

    # Real OS write/flush error in the installed process, after TSV writes fit.
    result = subprocess.run([str(Path(sysconfig.get_path("scripts")) / "sjsift"), *args,
                             "--output", str(outputs[0]), "--context-output", str(outputs[1]),
                             "--provenance-output", str(outputs[2])],
                            capture_output=True, text=True, preexec_fn=limit_file_size, check=False)
    assert result.returncode == 2 and not result.stdout
    assert "run.json: cannot write report" in result.stderr
    assert "Traceback" not in result.stderr
    assert not any(path.exists() for path in outputs)


def test_stdout_failure_removes_completed_provenance(tmp_path):
    args = cli_inputs(tmp_path)
    sidecar = tmp_path / "run.json"
    read_fd, write_fd = os.pipe()
    os.close(read_fd)
    with os.fdopen(write_fd, "wb") as output:
        result = subprocess.run([str(Path(sysconfig.get_path("scripts")) / "sjsift"), *args,
                                 "--provenance-output", str(sidecar)],
                                stdout=output, stderr=subprocess.PIPE, text=True, check=False)
    assert result.returncode == 2
    assert "standard output: cannot write report: Broken pipe" in result.stderr
    assert "Traceback" not in result.stderr and "Exception ignored" not in result.stderr
    assert not sidecar.exists()
