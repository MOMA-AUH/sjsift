"""Tests for the installed sjsift command-line interface."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import sysconfig

import pytest
from importlib.metadata import version

import sjsift


REPOSITORY_ROOT = Path(__file__).parents[1]


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


def test_exact_defining_junction_is_reported(tmp_path: Path) -> None:
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


def test_reference_catalog_reports_every_variant_in_catalog_order(
    tmp_path: Path,
) -> None:
    junctions = tmp_path / "sample.SJ.out.tab"
    junctions.write_text(
        """\
chrX\t67686127\t67694672\t1\t1\t0\t7\t1\t48
chr7\t116771655\t116774880\t1\t1\t1\t19\t0\t62
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
        "EGFRvIVa\tGRCh38\tchr7\t55200414\t55205255\t+\t0\t0\t0\n"
        "EGFRvIII\tGRCh38\tchr7\t55019366\t55155829\t+\t0\t0\t0\n"
        "EGFRvIIIb\tGRCh38\tchr7\t55109959\t55155829\t+\t0\t0\t0\n"
        "EGFRvIVb\tGRCh38\tchr7\t55200414\t55202516\t+\t3\t4\t7\n"
        "EGFRvII\tGRCh38\tchr7\t55161632\t55171174\t+\t11\t2\t13\n"
        "EGFRvIIb\tGRCh38\tchr7\t55161632\t55170306\t+\t0\t0\t0\n"
        "METx14del\tGRCh38\tchr7\t116771655\t116774880\t+\t19\t0\t19\n"
        "ARv7\tGRCh38\tchrX\t67686127\t67694672\t+\t7\t1\t8\n"
    )
    assert result.stderr == ""
