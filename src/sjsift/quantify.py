"""Match catalog definitions against STAR junction support."""

from dataclasses import dataclass
from contextlib import contextmanager
from collections.abc import Iterator
import gzip
import io
from pathlib import Path
from typing import NoReturn, TextIO
import zlib

from .catalog import Catalog, ReferenceJunction, VariantDefinition


class QuantifyError(ValueError):
    """A user-correctable problem with a STAR junction file."""


@dataclass(frozen=True)
class VariantSupport:
    """STAR support associated with a known splice variant."""

    definition: VariantDefinition
    unique: int
    multimapping: int
    reference_junctions: tuple["ReferenceJunctionSupport", ...] = ()

    @property
    def total(self) -> int:
        """Return unweighted total STAR support."""
        return self.unique + self.multimapping


@dataclass(frozen=True)
class ReferenceJunctionSupport:
    """STAR support associated with one catalogued reference junction."""

    definition: ReferenceJunction
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


@contextmanager
def _open_junctions(path: Path) -> Iterator[TextIO]:
    """Stream UTF-8 junctions, detecting gzip by its header rather than suffix."""
    with path.open("rb") as raw:
        binary = gzip.GzipFile(fileobj=raw) if raw.peek(2)[:2] == b"\x1f\x8b" else raw
        with io.TextIOWrapper(binary, encoding="utf-8", newline="") as stream:
            yield stream


def quantify(catalog: Catalog, junctions_path: Path) -> Quantification:
    """Validate STAR rows and copy support from exact defining-junction matches."""
    star_strands = {"+": 1, "-": 2}
    defining_targets = {
        (
            variant.chromosome,
            variant.intron_start,
            variant.intron_end,
            star_strands[variant.strand],
        ): index
        for index, variant in enumerate(catalog.variants)
    }
    counts = [(0, 0) for _ in catalog.variants]
    reference_targets: dict[tuple[str, int, int, int], list[tuple[int, int]]] = {}
    reference_counts = [
        [(0, 0) for _ in variant.reference_junctions] for variant in catalog.variants
    ]
    for variant_index, variant in enumerate(catalog.variants):
        for reference_index, reference in enumerate(variant.reference_junctions):
            key = (
                reference.chromosome,
                reference.intron_start,
                reference.intron_end,
                star_strands[reference.strand],
            )
            reference_targets.setdefault(key, []).append((variant_index, reference_index))

    matched_junctions: set[tuple[str, int, int, int]] = set()
    catalog_chromosomes = {
        junction[0] for junction in defining_targets | reference_targets
    }
    compatible_chromosome_seen = False

    try:
        with _open_junctions(junctions_path) as stream:
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
                if key not in defining_targets and key not in reference_targets:
                    continue

                if key in matched_junctions:
                    if key in defining_targets:
                        identifier = catalog.variants[defining_targets[key]].identifier
                        detail = f" for catalog variant {identifier!r}"
                    else:
                        detail = " for a catalog reference junction"
                    _fail(
                        junctions_path,
                        line_number,
                        f"duplicate STAR row{detail}",
                    )
                matched_junctions.add(key)
                if key in defining_targets:
                    counts[defining_targets[key]] = (unique, multimapping)
                for variant_index, reference_index in reference_targets.get(key, []):
                    reference_counts[variant_index][reference_index] = (unique, multimapping)
    except UnicodeDecodeError as error:
        _fail(junctions_path, None, f"invalid UTF-8: {error}")
    except OSError as error:
        detail = error.strerror or str(error)
        _fail(junctions_path, None, f"cannot read STAR junction file: {detail}")
    except (EOFError, zlib.error) as error:
        _fail(junctions_path, None, f"cannot read STAR junction file: {error}")

    return Quantification(
        results=tuple(
            VariantSupport(
                variant,
                unique,
                multimapping,
                tuple(
                    ReferenceJunctionSupport(reference, reference_unique, reference_multi)
                    for reference, (reference_unique, reference_multi) in zip(
                        variant.reference_junctions, reference_counts[index]
                    )
                ),
            )
            for index, (variant, (unique, multimapping)) in enumerate(
                zip(catalog.variants, counts)
            )
        ),
        compatibility_warning=not compatible_chromosome_seen,
    )
