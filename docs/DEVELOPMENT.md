# Development and release guide

## Development baseline

sjsift uses Python 3.11 or newer, a `src/` package layout, setuptools as its PEP 517 build backend, and `pyproject.toml` for project metadata. This follows useful conventions observed in [`MOMA-AUH/skua`](https://github.com/MOMA-AUH/skua/tree/4b2162c4d2df11830ba030420899af7cb6165d04) while omitting its BAM-oriented dependencies, subcommand hierarchy, and advanced domain machinery.

The Python Packaging User Guide recommends declaring the build backend in `[build-system]` and new-project metadata in `[project]` ([PyPA guide](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/)). Python 3.11 is the minimum because it supplies `tomllib` in the standard library.

### Local setup

Once the implementation files exist:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
pytest
```

Build and inspect distributions with:

```bash
python -m build
python -m twine check dist/*
```

Generated environments, build directories, coverage data, and distributions remain untracked.

## Proposed repository structure

```text
.
├── .github/
│   ├── dependabot.yml
│   └── workflows/
│       ├── ci.yml
│       └── publish.yml
├── conda-recipe/
│   ├── build.sh
│   └── meta.yaml
├── definitions/
│   └── grch38.toml
├── docs/
│   ├── adr/
│   ├── DECISIONS.md
│   ├── DEVELOPMENT.md
│   ├── IMPLEMENTATION_PLAN.md
│   └── MVP_SPEC.md
├── src/
│   └── sjsift/
│       ├── __init__.py
│       ├── __main__.py
│       ├── catalog.py
│       ├── cli.py
│       ├── quantify.py
│       └── report.py
├── tests/
│   ├── fixtures/
│   ├── test_catalog.py
│   ├── test_cli.py
│   ├── test_quantify.py
│   ├── test_reference_catalog.py
│   └── test_report.py
├── CONTEXT.md
├── LICENSE
├── README.md
└── pyproject.toml
```

The intended external seam is the `sjsift` command. `cli.py` only parses arguments, translates expected failures into messages and exit statuses, and coordinates deeper modules:

- `catalog.py` owns TOML decoding and all catalog invariants.
- `quantify.py` owns STAR validation, exact matching, and result construction behind one small callable interface.
- `report.py` owns the fixed TSV schema and destination handling.

Avoid generic `utils.py`, adapter hierarchies, plugin seams, or a separately supported Python library interface. New modules should exist only when they hide a coherent body of behavior behind a smaller interface.

## Dependency management

Runtime dependencies: none. Use `argparse`, `csv`, `dataclasses`, `pathlib`, and `tomllib` from the standard library as needed.

Development dependencies belong in a `dev` optional dependency group and initially comprise:

- `pytest` for automated tests;
- `build` for isolated sdist and wheel creation;
- `twine` for distribution metadata validation.

Do not add a formatter, linter, type checker, CLI framework, schema framework, dataframe library, or bioinformatics library in v0.1 unless a concrete defect demonstrates the need. `git diff --check` is the only initial whitespace check.

Keep the version static and explicit in `pyproject.toml` and the conda recipe. A small test must fail if they diverge; a custom version-management framework is unnecessary.

## Tests and fixtures

Tests are organized by behavior rather than internal helper functions. Use temporary paths for CLI and output-file tests and keep hand-written fixtures small enough to inspect directly.

Minimum fixture set:

- a valid catalog containing the eight entries in `definitions/grch38.toml`;
- a small STAR file with supported junctions carrying both unique and multimapping counts;
- a STAR file with near misses for chromosome spelling, off-by-one coordinates, opposite strand, and strand `0`;
- an empty STAR file;
- a STAR file whose chromosomes have no overlap with the catalog;
- malformed STAR rows covering field count, integer parsing, range checks, and duplicate target records;
- invalid TOML catalogs covering missing, unknown, mistyped, duplicate, and out-of-range values;
- exact expected TSV output, including zero-support rows and catalog ordering.

Unit tests cover catalog validation, STAR parsing/matching, arithmetic, and report serialization. CLI integration tests invoke `main(argv)` with real temporary files and assert stdout, stderr, output files, refusal to overwrite, and exit statuses. At least one installed-package smoke test runs the console entry point.

Every normative rule in [`MVP_SPEC.md`](MVP_SPEC.md) must have a named test or a direct build/install check. Test fixtures are synthetic evidence, not clinical validation samples.

## Continuous integration

GitHub Actions supports Python version matrices and standard build/test steps ([GitHub's Python Actions guide](https://docs.github.com/en/actions/tutorials/build-and-test-code/python)). Keep the initial workflow modest:

1. `tests (3.11)` through `tests (3.14)` each install the package with development dependencies and run pytest.
2. `package` runs `git diff --check`, builds the sdist and wheel, runs `twine check`, installs the wheel into a clean environment, imports `sjsift`, and smoke-tests `sjsift --help` and `sjsift --version`.

Run CI for pull requests targeting `master` and for pushes to `master`. The four test checks and `package` become required checks. Do not add operating-system matrices, coverage gates, formatting gates, linting, or static typing to v0.1.

Dependabot may check GitHub Actions weekly, matching the lightweight convention used by skua.

## Pull-request workflow

The repository uses `master` as its protected default branch.

1. Synchronize local `master` with `origin/master` using fast-forward only.
2. Create a descriptive branch before editing.
3. Keep each pull request limited to one implementation milestone or one coherent maintenance change.
4. Open a pull request into `master` and wait for every required CI check.
5. Do not merge while a required check is failing or pending.
6. Squash-merge through GitHub; do not merge or rebase-merge.
7. Let GitHub delete the merged head branch automatically.
8. Fast-forward the local `master` after merge and remove the local topic branch.

No reviewer approval is required for v0.1. A pull request and passing required checks are mandatory.

## Repository bootstrap and GitHub configuration

This repository's bootstrap push has already occurred: `master` exists on `origin` at the initial commit. Configure the ruleset and merge settings before any subsequent change is merged. For a new replacement repository, the one-time exception is to create and push the initial `master` commit, set it as the default branch, and immediately apply the settings below; all later changes use pull requests.

GitHub requires a candidate status check to have completed successfully in the repository during the preceding seven days before it can be selected as required ([GitHub troubleshooting documentation](https://docs.github.com/en/pull-requests/how-tos/merge-and-close-pull-requests/troubleshooting-required-status-checks)). Bootstrap the checks in this order:

1. Activate the ruleset's pull-request, linear-history, force-push, deletion, and no-bypass rules, and configure squash-only merging.
2. Add `.github/workflows/ci.yml` on a branch and merge that bootstrap pull request only after all of its checks pass, even though their names cannot yet be selected as required.
3. Immediately add `tests (3.11)`, `tests (3.12)`, `tests (3.13)`, `tests (3.14)`, and `package` to the ruleset's required checks.
4. Do not merge any later pull request while any required check is failing or pending.

Step 2 is the sole required-check bootstrap exception; it is not an exception to branch, pull-request, or passing-CI requirements.

### Configuration committed to the repository

- `.github/workflows/ci.yml` defines the checks.
- `.github/workflows/publish.yml` defines build and publication mechanics.
- `.github/dependabot.yml` schedules Actions updates.
- `pyproject.toml`, tests, and the conda recipe define build and validation behavior.
- This guide documents the intended policy.

Committed files create check runs, but they cannot make those checks mandatory or prohibit direct pushes by themselves.

### GitHub branch ruleset

Create an active branch ruleset targeting the default branch `master`. GitHub documents that multiple applicable rulesets accumulate and provides visibility into which rules apply ([ruleset documentation](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-rulesets/about-rulesets)). Configure it to:

- require a pull request before merging;
- require `tests (3.11)`, `tests (3.12)`, `tests (3.13)`, `tests (3.14)`, and `package` to pass;
- require linear history;
- block force pushes;
- block deletion of `master`;
- require zero approving reviews;
- define no bypass actors, so administrators do not become an undocumented direct-push path.

Do not require conversation resolution, signed commits, merge queue, deployments, code-owner review, or an up-to-date branch unless a later concrete need justifies the added friction.

### GitHub merge settings

In the repository's pull-request settings:

- enable squash merging;
- disable merge commits;
- disable rebase merging;
- enable automatic deletion of head branches after merge.

GitHub's linear-history rule alone is insufficient for squash-only policy because rebased commits are also linear. The allowed merge methods must be configured separately ([merge-method documentation](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/configuring-pull-request-merges/about-merge-methods-on-github), [squash configuration](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/configuring-pull-request-merges/configuring-commit-squashing-for-pull-requests)). Automatic head-branch deletion is another repository setting and can still be limited by branch rules ([GitHub documentation](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/configuring-pull-request-merges/managing-the-automatic-deletion-of-branches)).

These ruleset and merge settings cannot be enforced solely by files committed to this repository. Organization policy and GitHub plan capabilities may also constrain available settings; record any enforced equivalent if the exact UI control is unavailable.

## Release process

Use semantic versions and `vMAJOR.MINOR.PATCH` Git tags, following skua's release convention.

1. On a release branch, update the version in `pyproject.toml` and `conda-recipe/meta.yaml`, update user-facing documentation, and add or update release tests.
2. Open and squash-merge the release pull request after all required checks pass.
3. From the resulting `master` commit, create and push an annotated version tag.
4. The tag-triggered publish workflow rebuilds and tests the Python distributions and conda package from that exact checkout.
5. Retain build artifacts in GitHub Actions.
6. Publish the conda package only after its build and recipe tests pass.
7. Create the corresponding GitHub Release and summarize user-visible changes.

Manual workflow dispatch may build and test release artifacts but must not publish them.

## MOMA-AUH Anaconda publishing outline

Use a conventional `conda-recipe/meta.yaml` and mark the package `noarch: python` while sjsift remains pure Python. The recipe should:

- build from the checked-out repository source;
- install with pip without dependency resolution;
- declare Python `>=3.11` in host and run requirements;
- test `import sjsift`;
- test `sjsift --help` and `sjsift --version`;
- carry license and source metadata consistent with `pyproject.toml`.

The publish job installs `conda-build` and `anaconda-client`, builds the recipe, and uploads the resulting artifact with `anaconda upload --skip-existing --user MOMA-AUH`. This matches conda-build's documented build/upload flow ([conda-build tutorial](https://docs.conda.io/projects/conda-build/en/stable/user-guide/tutorials/build-pkgs.html#optional-uploading-new-packages-to-anaconda-org)) and the established [`skua` publish workflow](https://github.com/MOMA-AUH/skua/blob/4b2162c4d2df11830ba030420899af7cb6165d04/.github/workflows/publish.yml).

Store `ANACONDA_API_TOKEN` as a GitHub Actions secret available only to the publish job. The token, repository ruleset, required-check selection, merge-method settings, and automatic branch deletion are external GitHub/Anaconda configuration and cannot be committed to the repository.
