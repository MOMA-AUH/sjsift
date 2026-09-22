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
    assert "-o, --output PATH" in result.stdout
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


def test_quantification_is_explicitly_deferred() -> None:
    result = run_console_script(
        "--junctions",
        "sample.SJ.out.tab",
        "--definitions",
        "definitions.toml",
    )

    assert result.returncode == 2
    assert result.stdout == ""
    assert result.stderr == "sjsift: error: quantification is not implemented yet\n"
