"""The independent readiness audit accepts actual installed CLI reports."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig

import pysam
import pytest


FIXTURES = Path(__file__).parent / "fixtures" / "star"
AUDIT = Path(__file__).resolve().parents[1] / "benchmarks" / "readiness"


@pytest.mark.parametrize("format", ["bam", "cram"])
@pytest.mark.parametrize("header_metadata", [False, True])
def test_current_cli_reports_pass_independent_readiness_audit(tmp_path, format, header_metadata):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for source, target in [
        ("catalog.toml", "catalog.toml"), ("SJ.out.tab", "sample.SJ.out.tab"),
        ("alignments.bam", "reads.bam"), ("alignments.bam.bai", "reads.bam.bai"),
        ("reference.fa", "reference.fa"), ("reference.fa.fai", "reference.fa.fai"),
    ]:
        shutil.copyfile(FIXTURES / source, inputs / target)
    sample_options = ["--sample-name", "STAR interoperability"]
    if header_metadata:
        # Keep the committed STAR records unchanged; add independently known
        # header facts to exercise inference, blank SM, and partial AS coverage.
        with pysam.AlignmentFile(str(FIXTURES / "alignments.bam"), "rb") as source:
            header = source.header.to_dict()
            header["RG"] = [{"ID": "sample", "SM": "STAR interoperability"},
                            {"ID": "blank", "SM": " "}]
            header["SQ"][0]["AS"] = "STAR-synthetic"
            with pysam.AlignmentFile(str(inputs / "reads.bam"), "wb", header=header) as destination:
                for record in source:
                    destination.write(record)
        pysam.index(str(inputs / "reads.bam"))
        sample_options = []
    reference = inputs / "reference.fa"
    with pysam.AlignmentFile(str(inputs / "reads.bam"), "rb") as source:
        with pysam.AlignmentFile(str(inputs / "reads.cram"), "wc", header=source.header,
                                 reference_filename=str(reference)) as destination:
            for record in source:
                destination.write(record)
    pysam.index(str(inputs / "reads.cram"))
    report = tmp_path / "report"
    report.mkdir()
    result = subprocess.run([
        str(Path(sysconfig.get_path("scripts")) / "sjsift"),
        *sample_options,
        "--definitions", str(inputs / "catalog.toml"),
        "--junctions", str(inputs / "sample.SJ.out.tab"),
        "--alignment", str(inputs / f"reads.{format}"),
        "--output", str(report / "main.tsv"),
        "--context-output", str(report / "context.tsv"),
        "--html-output", str(report / "report.html"),
        "--provenance-output", str(report / "provenance.json"),
        *(["--reference", str(reference)] if format == "cram" else []),
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    inventory_path = tmp_path / "inventory.json"
    result = subprocess.run([
        sys.executable, str(AUDIT / "audit_evidence.py"),
        str(inputs), str(report), str(inventory_path),
        "--format", format, *sample_options,
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    inventory = json.loads(inventory_path.read_text())
    assert inventory["status"] == "passed"
    assert inventory["embedded_rows"] == 10
    assert inventory["eligible_memberships"] == 15
    sensitivity_path = tmp_path / "sensitivity.json"
    result = subprocess.run([
        sys.executable, str(AUDIT / "check_audit.py"),
        str(inputs), str(report), str(sensitivity_path),
        "--format", format, *sample_options,
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(sensitivity_path.read_text()) == {
        "unmodified_report": "passed", "restored_report": "passed",
        "wrong_star_count": "rejected", "missing_physical_record": "rejected",
        "wrong_pair_association": "rejected", "same_name_mate_substitution": "rejected",
        "wrong_provenance_schema": "rejected", "wrong_sample_selection": "rejected",
    }
