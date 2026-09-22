# sjsift v0.1 implementation plan

This plan implements [`MVP_SPEC.md`](MVP_SPEC.md) in small pull requests. Each milestone is independently reviewable and must preserve all explicit v0.1 exclusions.

## Milestone 0: finish repository-policy bootstrap

### Scope

Protect `master` before implementation work and establish the one-time path for introducing required CI checks.

### Expected files or configuration

- No production files.
- GitHub branch ruleset and merge settings described in [`DEVELOPMENT.md`](DEVELOPMENT.md).

### Verification

- Confirm direct pushes, force pushes, and deletion of `master` are blocked.
- Confirm pull requests are required with zero approving reviews.
- Confirm squash is the only enabled merge method and merged branches are deleted automatically.

### Completion criteria

- The ruleset is active with no bypass actors.
- The initial CI workflow must be added by pull request.
- After that workflow produces named checks, `tests (3.11)`, `tests (3.12)`, `tests (3.13)`, `tests (3.14)`, and `package` are added as required checks before later pull requests merge.

### Explicitly out of scope

- Reviewer approval requirements
- Merge queue, signed commits, deployments, or code owners
- Any product implementation

## Milestone 1: package and CLI skeleton

### Scope

Create an installable Python package with a console entry point, version/help behavior, development dependencies, and modest CI. The quantification invocation may fail with a clear “not implemented” message until later milestones; no biological behavior is implemented here.

### Expected files or components

- `pyproject.toml`
- `src/sjsift/__init__.py`
- `src/sjsift/__main__.py`
- `src/sjsift/cli.py`
- `tests/test_cli.py`
- `.github/workflows/ci.yml`
- `.github/dependabot.yml`

### Tests to add

- Package version is importable.
- `python -m sjsift --help` and installed `sjsift --help` succeed.
- `--version` prints the package version.
- Required CLI arguments are named as specified.
- Package build, metadata check, wheel installation, import, and CLI smoke checks run in CI.

### Completion criteria

- Tests pass on Python 3.11, 3.12, 3.13, and 3.14.
- The `package` job validates both sdist and wheel and installs the wheel cleanly.
- The four test check names and `package` are configured as required checks in GitHub before the next milestone merges.

### Explicitly out of scope

- Catalog parsing
- STAR parsing or matching
- TSV result generation
- Conda packaging or publication
- Formatting, linting, typing, coverage, or OS matrices

## Milestone 2: catalog loading and validation

### Scope

Implement the complete TOML schema and return an ordered, immutable in-memory catalog behind one small module interface.

### Expected files or components

- `src/sjsift/catalog.py`
- `tests/test_catalog.py`
- `tests/test_reference_catalog.py`
- Focused fixtures under `tests/fixtures/catalogs/`

### Tests to add

- Valid schema version, assembly, and variant fields parse successfully.
- Entry order is preserved.
- Missing, unknown, mistyped, empty, and out-of-range fields fail with useful entry context.
- Invalid strand fails.
- Duplicate IDs fail.
- Duplicate defining junctions under different IDs fail.
- Unsupported schema versions fail explicitly.
- [`definitions/grch38.toml`](../definitions/grch38.toml) contains exactly the six supported EGFR variants, `METx14del`, and `ARv7` with the specified coordinates and strand.

### Completion criteria

- Every catalog rule in the MVP specification has a focused test.
- No third-party TOML or schema dependency is introduced.
- Error messages identify the catalog and relevant variant where possible.

### Explicitly out of scope

- STAR parsing
- Chromosome normalization or liftover
- Catalog aliases, downloads, built-in selection, or migrations beyond rejecting unsupported versions
- Multi-junction variants

## Milestone 3: STAR validation and exact quantification

### Scope

Implement streaming validation of one `SJ.out.tab`, exact lookup of catalog junctions, and ordered result construction.

### Expected files or components

- `src/sjsift/quantify.py`
- `tests/test_quantify.py`
- Focused fixtures under `tests/fixtures/star/`

### Tests to add

- Each of the nine STAR columns accepts valid documented values.
- Wrong field counts, invalid integers, invalid ranges, invalid strand/motif/annotation codes, and negative counts fail with line numbers.
- Exact chromosome/start/end/strand matches copy unique and multimapping support.
- `total_support` is their unweighted sum.
- `chr7` versus `7`, coordinate differences of one base, opposite strand, and STAR strand `0` do not match.
- Motif, annotation, overhang, and STAR row order do not change matching.
- Missing junctions and an empty STAR file produce ordered zero results.
- Duplicate records matching one catalog target fail rather than sum or override.
- No chromosome overlap produces exactly one warning condition.
- Mixed supported and unsupported variants retain catalog order.

### Completion criteria

- Matching and counting tests cover all acceptance cases without BAM, CRAM, or SAM fixtures.
- The quantification interface returns results and warning information without writing output or exiting the process.
- Memory use depends primarily on catalog size rather than retaining all STAR rows.

### Explicitly out of scope

- Alignment-level processing or filtering
- Novel junction discovery
- Approximate matching, chromosome aliases, strand fallback, or coordinate conversion
- Thresholds, ratios, classifications, or multi-junction aggregation
- Speculative parallelism or indexing

## Milestone 4: deterministic reporting and complete CLI

### Scope

Serialize the fixed TSV schema and connect catalog loading, quantification, warnings, files, and exit statuses through the command-line interface.

### Expected files or components

- `src/sjsift/report.py`
- Completed `src/sjsift/cli.py`
- `tests/test_report.py`
- Expanded `tests/test_cli.py`
- Exact expected-output fixtures under `tests/fixtures/expected/`

### Tests to add

- Header names and column order match the specification exactly.
- One row is emitted per variant in catalog order, including zero-support rows.
- Integer and strand formatting is deterministic.
- Stdout contains only TSV when no output path is supplied.
- `--output` writes the same bytes and normal success leaves stderr empty.
- Existing output paths are refused without modification.
- Expected input and output failures emit concise stderr messages without tracebacks and exit `2`.
- Compatibility warnings use stderr, still emit a complete report, and exit `0`.
- An injected unexpected failure exits `1` through the top-level entry point.
- A representative whole-file run covers all eight reference-catalog variants.

### Completion criteria

- Every CLI, output, warning, and exit-status acceptance criterion is covered.
- README commands and representative output agree with actual behavior.
- No subcommands or public Python-library guarantees have been added.

### Explicitly out of scope

- Input from stdin or compressed files
- Sample-name inference
- JSON, CSV, VCF, or visualization output
- Overwrite flags or interactive prompts
- Workflow-engine wrappers

## Milestone 5: conda package and release automation

### Scope

Add a tested pure-Python conda recipe and tag-triggered publication to the `MOMA-AUH` Anaconda organization.

### Expected files or components

- `conda-recipe/meta.yaml`
- `conda-recipe/build.sh`
- `.github/workflows/publish.yml`
- Version-synchronization test
- Any final release-command corrections in `docs/DEVELOPMENT.md`

### Tests to add

- Recipe and Python package versions agree.
- Conda package build succeeds.
- Recipe tests import `sjsift` and run `sjsift --help` and `sjsift --version`.
- Manual workflow dispatch builds and tests but does not reach the upload step.
- Tag-trigger conditions and publish-job dependencies are inspectable in workflow tests or review.

### Completion criteria

- A dry-run/manual Actions build produces retained Python and conda artifacts without publishing.
- `ANACONDA_API_TOKEN` is configured for the publish job.
- A protected-`master` tag can publish only after all build and recipe tests pass.
- Release instructions match the implemented workflow.

### Explicitly out of scope

- PyPI publication
- Multiple conda channels or platform builds
- Automatic GitHub Release generation unless it can be added without obscuring the package publication path
- Nightly builds or automatic catalog updates

## Requirement-to-test traceability

| Requirement area | Primary verification |
|---|---|
| Catalog schema and uniqueness | `test_catalog.py`, `test_reference_catalog.py` |
| STAR field validation | `test_quantify.py` malformed-field cases |
| Exact coordinates, chromosome, and strand | `test_quantify.py` exact and near-miss cases |
| Unique, multimapping, and total counts | `test_quantify.py` count cases |
| Zero-support and ordering behavior | `test_quantify.py`, `test_report.py` |
| Output columns and formatting | `test_report.py` golden TSV fixture |
| CLI streams, files, warnings, and exits | `test_cli.py` integration cases |
| Installation and import | CI `package` job |
| Supported Python versions | CI `tests` matrix endpoints |
| Conda construction and CLI smoke tests | Recipe tests and publish-workflow dry run |
| Pull-request and protected-branch policy | GitHub settings audit described in Milestone 0 |
