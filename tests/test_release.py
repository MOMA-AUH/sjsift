"""Tests for release-package metadata and automation."""

from __future__ import annotations

from pathlib import Path
import re
import tomllib


REPOSITORY_ROOT = Path(__file__).parents[1]


def test_conda_recipe_version_matches_python_package() -> None:
    pyproject = tomllib.loads(
        (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )
    recipe = (REPOSITORY_ROOT / "conda-recipe" / "meta.yaml").read_text(
        encoding="utf-8"
    )
    recipe_version = re.search(
        r'^\{% set version = "(?P<version>[^\"]+)" %\}$', recipe, re.MULTILINE
    )

    assert recipe_version is not None, "conda recipe must declare a literal version"
    assert recipe_version.group("version") == pyproject["project"]["version"]


def test_manual_release_build_cannot_publish() -> None:
    workflow = (
        REPOSITORY_ROOT / ".github" / "workflows" / "publish.yml"
    ).read_text(encoding="utf-8")

    assert "workflow_dispatch:" in workflow
    assert '"v[0-9]+.[0-9]+.[0-9]+"' in workflow
    assert "needs: [python-distributions, conda-package]" in workflow
    assert "if: github.event_name == 'push' && startsWith(github.ref, 'refs/tags/v')" in workflow
    assert "actions/upload-artifact@v7" in workflow
    assert "ANACONDA_API_TOKEN: ${{ secrets.ANACONDA_API_TOKEN }}" in workflow
    assert "anaconda upload --skip-existing --user MOMA-AUH" in workflow
