"""Match catalog definitions against STAR junction support."""

from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn

from .catalog import Catalog, VariantDefinition


class QuantifyError(ValueError):
    """A user-correctable problem with a STAR junction file."""


@dataclass(frozen=True)
class VariantSupport:
    """STAR support associated with a known splice variant."""

    definition: VariantDefinition
    unique: int
    multimapping: int

    @property
    def total(self) -> int:
        """Return unweighted total STAR support."""
        return self.unique + self.multimapping


@dataclass(frozen=True)
class Quantification:
    """Ordered support results and input-compatibility information."""

    results: tuple[VariantSupport, ...]
    compatibility_warning: bool


def _fail(path: Path, line_number: int | None, message: str) -> NoReturn:
    location = str(path)
    if line_number is not None:
        location += f": line {line_number}"
    raise QuantifyError(f"{location}: {message}")


def _parse_integer(
    path: Path, line_number: int, value: str, field_name: str
) -> int:
    if value != value.strip():
        _fail(path, line_number, f"{field_name} must be an integer")
    try:
        return int(value)
    except ValueError:
        _fail(path, line_number, f"{field_name} must be an integer")


def _parse_row(
    path: Path, line_number: int, line: str
) -> tuple[str, int, int, int, int, int]:
    fields = line.split("\t")
    if len(fields) != 9:
        _fail(path, line_number, "expected exactly 9 tab-separated fields")

    chromosome = fields[0]
    if not chromosome:
        _fail(path, line_number, "chromosome must not be empty")

    intron_start = _parse_integer(path, line_number, fields[1], "intron start")
    intron_end = _parse_integer(path, line_number, fields[2], "intron end")
    strand = _parse_integer(path, line_number, fields[3], "strand")
    motif = _parse_integer(path, line_number, fields[4], "motif")
    annotation = _parse_integer(path, line_number, fields[5], "annotation status")
    unique = _parse_integer(path, line_number, fields[6], "unique support")
    multimapping = _parse_integer(
        path, line_number, fields[7], "multimapping support"
    )
    overhang = _parse_integer(path, line_number, fields[8], "overhang")

    if intron_start < 1:
        _fail(path, line_number, "intron start must be positive")
    if intron_end < 1:
        _fail(path, line_number, "intron end must be positive")
    if intron_end < intron_start:
        _fail(path, line_number, "intron end must not be less than intron start")
    if strand not in {0, 1, 2}:
        _fail(path, line_number, "strand must be one of 0, 1, or 2")
    if not 0 <= motif <= 6:
        _fail(path, line_number, "motif must be between 0 and 6")
    if annotation not in {0, 1}:
        _fail(path, line_number, "annotation status must be 0 or 1")
    if unique < 0:
        _fail(path, line_number, "unique support must not be negative")
    if multimapping < 0:
        _fail(path, line_number, "multimapping support must not be negative")
    if overhang < 0:
        _fail(path, line_number, "overhang must not be negative")

    return chromosome, intron_start, intron_end, strand, unique, multimapping


def quantify(catalog: Catalog, junctions_path: Path) -> Quantification:
    """Validate STAR rows and copy support from exact defining-junction matches."""
    star_strands = {"+": 1, "-": 2}
    targets = {
        (
            variant.chromosome,
            variant.intron_start,
            variant.intron_end,
            star_strands[variant.strand],
        ): index
        for index, variant in enumerate(catalog.variants)
    }
    counts = [(0, 0) for _ in catalog.variants]
    matched_targets: set[int] = set()
    catalog_chromosomes = {variant.chromosome for variant in catalog.variants}
    compatible_chromosome_seen = False

    try:
        with junctions_path.open(encoding="utf-8", newline="") as stream:
            for line_number, raw_line in enumerate(stream, start=1):
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                chromosome, start, end, strand, unique, multimapping = _parse_row(
                    junctions_path, line_number, line
                )
                if chromosome in catalog_chromosomes:
                    compatible_chromosome_seen = True
                key = (chromosome, start, end, strand)
                if key not in targets:
                    continue

                target = targets[key]
                if target in matched_targets:
                    identifier = catalog.variants[target].identifier
                    _fail(
                        junctions_path,
                        line_number,
                        f"duplicate STAR row for catalog variant {identifier!r}",
                    )
                counts[target] = (unique, multimapping)
                matched_targets.add(target)
    except UnicodeDecodeError as error:
        _fail(junctions_path, None, f"invalid UTF-8: {error}")
    except OSError as error:
        detail = error.strerror or str(error)
        _fail(junctions_path, None, f"cannot read STAR junction file: {detail}")

    return Quantification(
        results=tuple(
            VariantSupport(variant, unique, multimapping)
            for variant, (unique, multimapping) in zip(catalog.variants, counts)
        ),
        compatibility_warning=not compatible_chromosome_seen,
    )
