"""Tests for the installed sjsift command-line interface."""

from __future__ import annotations

import gzip
import io
import os
from importlib.metadata import version
from pathlib import Path
import subprocess
import sys
import sysconfig

import pytest

import sjsift
from sjsift import cli


REPOSITORY_ROOT = Path(__file__).parents[1]


def write_single_variant_inputs(tmp_path: Path) -> tuple[Path, Path]:
    definitions = tmp_path / "definitions.toml"
    definitions.write_text(
        """\
schema_version = 1
genome_assembly = "GRCh38"

[[variants]]
id = "EGFRvIII"
chromosome = "chr7"
intron_start = 55019366
intron_end = 55155829
strand = "+"
""",
        encoding="utf-8",
    )
    junctions = tmp_path / "sample.SJ.out.tab"
    junctions.write_text(
        "chr7\t55019366\t55155829\t1\t1\t1\t23\t4\t71\n",
        encoding="utf-8",
    )
    return junctions, definitions


def run_module(*arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "sjsift", *arguments],
        check=False,
        capture_output=True,
        text=True,
    )


def run_console_script(*arguments: str) -> subprocess.CompletedProcess[str]:
    suffix = ".exe" if os.name == "nt" else ""
    executable = Path(sysconfig.get_path("scripts")) / f"sjsift{suffix}"
    assert executable.is_file(), "tests must run against an installed package"
    return subprocess.run(
        [str(executable), *arguments],
        check=False,
        capture_output=True,
        text=True,
    )


def test_package_version_matches_installed_metadata() -> None:
    assert sjsift.__version__ == version("sjsift")


@pytest.mark.parametrize("invoke", [run_module, run_console_script])
def test_help_describes_single_command_interface(invoke) -> None:
    result = invoke("--help")

    assert result.returncode == 0
    assert result.stderr == ""
    assert "--junctions PATH" in result.stdout
    assert "--definitions PATH" in result.stdout
    assert "-o" in result.stdout
    assert "--output PATH" in result.stdout
    assert "{command}" not in result.stdout


@pytest.mark.parametrize("invoke", [run_module, run_console_script])
def test_version_succeeds(invoke) -> None:
    result = invoke("--version")

    assert result.returncode == 0
    assert result.stdout == f"sjsift {sjsift.__version__}\n"
    assert result.stderr == ""


def test_required_arguments_are_enforced() -> None:
    result = run_console_script()

    assert result.returncode == 2
    assert "--junctions" in result.stderr
    assert "--definitions" in result.stderr


@pytest.mark.parametrize("compressed", [False, True], ids=["plain", "gzip"])
@pytest.mark.parametrize("suffix", [".tab", ".tab.gz"])
def test_exact_defining_junction_is_reported(
    tmp_path: Path, compressed: bool, suffix: str
) -> None:
    definitions = tmp_path / "definitions.toml"
    definitions.write_text(
        """\
schema_version = 1
genome_assembly = "GRCh38"

[[variants]]
id = "EGFRvIII"
chromosome = "chr7"
intron_start = 55019366
intron_end = 55155829
strand = "+"
""",
        encoding="utf-8",
    )
    junctions = tmp_path / f"sample.SJ.out{suffix}"
    junctions.write_text(
        "chr7\t55019366\t55155829\t1\t1\t1\t23\t4\t71\n",
        encoding="utf-8",
    )

    if compressed:
        junctions.write_bytes(gzip.compress(junctions.read_bytes()))

    result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
    )

    assert result.returncode == 0
    assert result.stdout == (
        "variant_id\tgenome_assembly\tchromosome\tintron_start\tintron_end\tstrand\t"
        "unique_support\tmultimapping_support\ttotal_support\n"
        "EGFRvIII\tGRCh38\tchr7\t55019366\t55155829\t+\t23\t4\t27\n"
    )
    assert result.stderr == ""
    assert junctions.name not in result.stdout


@pytest.mark.parametrize("option", ["-o", "--output"])
def test_output_option_writes_the_stdout_report_byte_for_byte(
    tmp_path: Path, option: str
) -> None:
    junctions, definitions = write_single_variant_inputs(tmp_path)
    stdout_result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
    )
    output = tmp_path / f"report-{option.lstrip('-')}.tsv"

    file_result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
        option,
        str(output),
    )

    assert file_result.returncode == 0
    assert file_result.stdout == ""
    assert file_result.stderr == ""
    assert output.read_bytes() == stdout_result.stdout.encode("utf-8")


def test_existing_output_is_refused_without_modification(tmp_path: Path) -> None:
    junctions, definitions = write_single_variant_inputs(tmp_path)
    output = tmp_path / "report.tsv"
    original = b"do not replace me\n"
    output.write_bytes(original)

    result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
        "--output",
        str(output),
    )

    assert result.returncode == 2
    assert result.stdout == ""
    assert str(output) in result.stderr
    assert "already exists" in result.stderr
    assert "Traceback" not in result.stderr
    assert output.read_bytes() == original


@pytest.mark.parametrize(
    "junction_text",
    ["", "chr1\t1\t2\t0\t0\t0\t5\t6\t7\n"],
    ids=["empty", "no-catalog-chromosomes"],
)
def test_incompatible_star_file_reports_zeros_with_one_compatibility_warning(
    tmp_path: Path, junction_text: str
) -> None:
    junctions, definitions = write_single_variant_inputs(tmp_path)
    junctions.write_text(junction_text, encoding="utf-8")

    result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
    )

    assert result.returncode == 0
    assert result.stdout == (
        "variant_id\tgenome_assembly\tchromosome\tintron_start\tintron_end\tstrand\t"
        "unique_support\tmultimapping_support\ttotal_support\n"
        "EGFRvIII\tGRCh38\tchr7\t55019366\t55155829\t+\t0\t0\t0\n"
    )
    assert result.stderr == (
        f"sjsift: warning: {junctions}: no catalog chromosome identifiers found; "
        "check genome assembly and chromosome naming compatibility\n"
    )


def test_reference_catalog_reports_every_variant_in_catalog_order(
    tmp_path: Path,
) -> None:
    junctions = tmp_path / "sample.SJ.out.tab"
    junctions.write_text(
        """\
chrX\t67686127\t67694672\t1\t1\t0\t7\t1\t48
chr7\t116771655\t116774880\t1\t1\t1\t19\t0\t62
chr7\t116755516\t116758458\t1\t1\t0\t5\t1\t40
chr7\t140787585\t140834608\t2\t1\t0\t13\t2\t37
chr10\t121482178\t121483697\t2\t1\t0\t4\t0\t35
chr7\t55161632\t55171174\t1\t2\t1\t11\t2\t35
chr7\t55200414\t55202516\t1\t1\t0\t3\t4\t29
""",
        encoding="utf-8",
    )

    result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(REPOSITORY_ROOT / "definitions" / "grch38.toml"),
    )

    assert result.returncode == 0
    assert result.stdout == (
        "variant_id\tgenome_assembly\tchromosome\tintron_start\tintron_end\tstrand\t"
        "unique_support\tmultimapping_support\ttotal_support\n"
        "ARv7\tGRCh38\tchrX\t67686127\t67694672\t+\t7\t1\t8\n"
        "ARv567es\tGRCh38\tchrX\t67711690\t67723685\t+\t0\t0\t0\n"
        "BRAFdel2-10\tGRCh38\tchr7\t140781694\t140924565\t-\t0\t0\t0\n"
        "BRAFdel2-8\tGRCh38\tchr7\t140787585\t140924565\t-\t0\t0\t0\n"
        "BRAFdel3-10\tGRCh38\tchr7\t140781694\t140850110\t-\t0\t0\t0\n"
        "BRAFdel3-8\tGRCh38\tchr7\t140787585\t140850110\t-\t0\t0\t0\n"
        "BRAFdel4-10\tGRCh38\tchr7\t140781694\t140834608\t-\t0\t0\t0\n"
        "BRAFdel4-8\tGRCh38\tchr7\t140787585\t140834608\t-\t13\t2\t15\n"
        "EGFRvIVa\tGRCh38\tchr7\t55200414\t55205255\t+\t0\t0\t0\n"
        "EGFRvIII\tGRCh38\tchr7\t55019366\t55155829\t+\t0\t0\t0\n"
        "EGFRvIIIb\tGRCh38\tchr7\t55109959\t55155829\t+\t0\t0\t0\n"
        "EGFRvIVb\tGRCh38\tchr7\t55200414\t55202516\t+\t3\t4\t7\n"
        "EGFRvII\tGRCh38\tchr7\t55161632\t55171174\t+\t11\t2\t13\n"
        "EGFRvIIb\tGRCh38\tchr7\t55161632\t55170306\t+\t0\t0\t0\n"
        "FGFR2-E18-C3\tGRCh38\tchr10\t121482178\t121483697\t-\t4\t0\t4\n"
        "METex7-8\tGRCh38\tchr7\t116755516\t116758458\t+\t5\t1\t6\n"
        "METex14\tGRCh38\tchr7\t116771655\t116774880\t+\t19\t0\t19\n"
    )
    assert result.stderr == ""


@pytest.mark.parametrize(
    ("catalog_text", "expected_message"),
    [
        (
            'schema_version = "one"\n'
            'genome_assembly = "GRCh38"\n'
            "[[variants]]\n"
            'id = "test"\n'
            'chromosome = "chr1"\n'
            "intron_start = 1\n"
            "intron_end = 2\n"
            'strand = "+"\n',
            "schema_version",
        ),
        (
            'schema_version = 1\n'
            'genome_assembly = "GRCh38"\n'
            "[[variants]]\n"
            'id = "test"\n'
            "intron_start = 1\n"
            "intron_end = 2\n"
            'strand = "+"\n',
            "variant 1 ('test'): missing field(s) 'chromosome'",
        ),
        ('schema_version = 1\ngenome_assembly = "GRCh38"\nvariants = [', "invalid TOML"),
    ],
)
def test_invalid_catalog_is_a_concise_cli_error(
    tmp_path: Path, catalog_text: str, expected_message: str
) -> None:
    definitions = tmp_path / "invalid.toml"
    definitions.write_text(catalog_text, encoding="utf-8")
    junctions = tmp_path / "sample.SJ.out.tab"
    junctions.write_text("", encoding="utf-8")

    result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
    )

    assert result.returncode == 2
    assert result.stdout == ""
    assert str(definitions) in result.stderr
    assert expected_message in result.stderr
    assert "Traceback" not in result.stderr
    assert "variant_id" not in result.stderr


def test_unreadable_catalog_is_a_concise_cli_error(tmp_path: Path) -> None:
    definitions = tmp_path / "missing.toml"
    junctions = tmp_path / "sample.SJ.out.tab"
    junctions.write_text("", encoding="utf-8")

    result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
    )

    assert result.returncode == 2
    assert result.stdout == ""
    assert str(definitions) in result.stderr
    assert "cannot read catalog" in result.stderr
    assert "Traceback" not in result.stderr


def test_unreadable_star_file_is_a_concise_cli_error(tmp_path: Path) -> None:
    _, definitions = write_single_variant_inputs(tmp_path)
    junctions = tmp_path / "missing.SJ.out.tab"

    result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
    )

    assert result.returncode == 2
    assert result.stdout == ""
    assert str(junctions) in result.stderr
    assert "cannot read STAR junction file" in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize(
    ("data", "message"),
    [
        (gzip.compress(b"\ninvalid\n"), "line 2: expected exactly 9"),
        (gzip.compress(b"\xff\n"), "invalid UTF-8"),
        (gzip.compress(b"")[:-4], "cannot read STAR junction file"),
        (gzip.compress(b"\n")[:-8] + b"\x00" * 8, "cannot read STAR junction file"),
        (b"\x1f\x8b\x08\x00" + b"\x00" * 6 + b"\x07", "cannot read STAR junction file"),
    ],
    ids=["invalid-row", "invalid-utf8", "truncated", "bad-crc", "bad-deflate"],
)
def test_invalid_gzip_input_is_a_concise_error_without_output(
    tmp_path: Path, data: bytes, message: str
) -> None:
    junctions, definitions = write_single_variant_inputs(tmp_path)
    junctions.write_bytes(data)
    output = tmp_path / "report.tsv"

    result = run_console_script(
        "--junctions", str(junctions), "--definitions", str(definitions),
        "--output", str(output),
    )

    assert result.returncode == 2
    assert result.stdout == ""
    assert str(junctions) in result.stderr
    assert message in result.stderr
    assert "Traceback" not in result.stderr
    assert "internal error" not in result.stderr
    assert not output.exists()


def test_output_creation_failure_is_a_concise_cli_error(tmp_path: Path) -> None:
    junctions, definitions = write_single_variant_inputs(tmp_path)
    output = tmp_path / "missing-directory" / "report.tsv"

    result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
        "--output",
        str(output),
    )

    assert result.returncode == 2
    assert result.stdout == ""
    assert str(output) in result.stderr
    assert "cannot write report" in result.stderr
    assert "Traceback" not in result.stderr
    assert not output.exists()


def test_input_failure_does_not_create_a_partial_output(tmp_path: Path) -> None:
    junctions, definitions = write_single_variant_inputs(tmp_path)
    junctions.write_text("malformed\n", encoding="utf-8")
    output = tmp_path / "report.tsv"

    result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
        "--output",
        str(output),
    )

    assert result.returncode == 2
    assert result.stdout == ""
    assert "Traceback" not in result.stderr
    assert not output.exists()


def test_unexpected_top_level_failure_returns_status_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    junctions, definitions = write_single_variant_inputs(tmp_path)
    errors = io.StringIO()

    class FailingStandardOutput:
        def write(self, value: str) -> int:
            raise RuntimeError("simulated internal failure")

    monkeypatch.setattr(cli.sys, "stdout", FailingStandardOutput())
    monkeypatch.setattr(cli.sys, "stderr", errors)

    status = cli.main(
        [
            "--junctions",
            str(junctions),
            "--definitions",
            str(definitions),
        ]
    )

    assert status == 1
    assert errors.getvalue() == "sjsift: internal error: simulated internal failure\n"


@pytest.mark.parametrize(
    ("junction_text", "expected_message"),
    [
        ("chr7\t10\t20\t1\n", "exactly 9 tab-separated fields"),
        ("chr7\tbad\t20\t1\t0\t0\t0\t0\t0\n", "intron start must be an integer"),
        (
            "chr7\t10\t20\t1\t0\t0\t1\t2\t3\n"
            "chr7\t10\t20\t1\t0\t0\t4\t5\t6\n",
            "duplicate STAR row",
        ),
    ],
)
def test_invalid_star_input_is_a_concise_cli_error(
    tmp_path: Path, junction_text: str, expected_message: str
) -> None:
    definitions = tmp_path / "definitions.toml"
    definitions.write_text(
        """\
schema_version = 1
genome_assembly = "GRCh38"

[[variants]]
id = "test"
chromosome = "chr7"
intron_start = 10
intron_end = 20
strand = "+"
""",
        encoding="utf-8",
    )
    junctions = tmp_path / "sample.SJ.out.tab"
    junctions.write_text(junction_text, encoding="utf-8")

    result = run_console_script(
        "--junctions",
        str(junctions),
        "--definitions",
        str(definitions),
    )

    assert result.returncode == 2
    assert result.stdout == ""
    assert str(junctions) in result.stderr
    assert "line " in result.stderr
    assert expected_message in result.stderr
    assert "Traceback" not in result.stderr
    assert "variant_id" not in result.stderr
