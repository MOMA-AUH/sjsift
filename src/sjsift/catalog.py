"""Load and validate variant-definition catalogs."""

from dataclasses import dataclass
from pathlib import Path
import tomllib
from typing import NoReturn


_TOP_LEVEL_FIELDS = frozenset({"schema_version", "genome_assembly", "variants"})
_VARIANT_FIELDS = frozenset(
    {"id", "chromosome", "intron_start", "intron_end", "strand"}
)


class CatalogError(ValueError):
    """A user-correctable problem with a variant-definition catalog."""


@dataclass(frozen=True)
class VariantDefinition:
    """One known splice variant and its defining junction."""

    identifier: str
    chromosome: str
    intron_start: int
    intron_end: int
    strand: str


@dataclass(frozen=True)
class Catalog:
    """An ordered set of definitions in one genome assembly."""

    genome_assembly: str
    variants: tuple[VariantDefinition, ...]


def _fail(path: Path, message: str) -> NoReturn:
    raise CatalogError(f"{path}: {message}")


def _describe_fields(fields: set[str]) -> str:
    return ", ".join(repr(field) for field in sorted(fields))


def _validate_fields(
    path: Path,
    value: dict[str, object],
    expected: frozenset[str],
    context: str,
) -> None:
    actual = set(value)
    missing = expected - actual
    unknown = actual - expected
    problems = []
    if missing:
        problems.append(f"missing field(s) {_describe_fields(missing)}")
    if unknown:
        problems.append(f"unknown field(s) {_describe_fields(unknown)}")
    if problems:
        _fail(path, f"{context}: {'; '.join(problems)}")


def _validate_string(path: Path, value: object, field: str, context: str) -> str:
    if not isinstance(value, str):
        _fail(path, f"{context}: field {field!r} must be a string")
    if not value:
        _fail(path, f"{context}: field {field!r} must not be empty")
    if "\t" in value or "\n" in value or "\r" in value:
        _fail(path, f"{context}: field {field!r} must not contain tabs or line breaks")
    return value


def _validate_positive_integer(
    path: Path, value: object, field: str, context: str
) -> int:
    # bool is a subclass of int, but TOML booleans are not schema integers.
    if type(value) is not int:
        _fail(path, f"{context}: field {field!r} must be an integer")
    if value < 1:
        _fail(path, f"{context}: field {field!r} must be a positive integer")
    return value


def _variant_context(index: int, value: dict[str, object]) -> str:
    identifier = value.get("id")
    if isinstance(identifier, str) and identifier:
        return f"variant {index} ({identifier!r})"
    return f"variant {index}"


def _parse_variant(path: Path, value: object, index: int) -> VariantDefinition:
    context = f"variant {index}"
    if not isinstance(value, dict):
        _fail(path, f"{context}: must be a table")

    context = _variant_context(index, value)
    _validate_fields(path, value, _VARIANT_FIELDS, context)
    identifier = _validate_string(path, value["id"], "id", context)
    context = f"variant {index} ({identifier!r})"
    chromosome = _validate_string(path, value["chromosome"], "chromosome", context)
    intron_start = _validate_positive_integer(
        path, value["intron_start"], "intron_start", context
    )
    intron_end = _validate_positive_integer(
        path, value["intron_end"], "intron_end", context
    )
    if intron_end < intron_start:
        _fail(path, f"{context}: field 'intron_end' must not be less than 'intron_start'")
    strand = value["strand"]
    if not isinstance(strand, str) or strand not in {"+", "-"}:
        _fail(path, f"{context}: field 'strand' must be '+' or '-'")

    return VariantDefinition(
        identifier=identifier,
        chromosome=chromosome,
        intron_start=intron_start,
        intron_end=intron_end,
        strand=strand,
    )


def _validate_unique(path: Path, variants: tuple[VariantDefinition, ...]) -> None:
    identifiers: dict[str, int] = {}
    junctions: dict[tuple[str, int, int, str], tuple[int, str]] = {}
    for index, variant in enumerate(variants, start=1):
        if variant.identifier in identifiers:
            first_index = identifiers[variant.identifier]
            _fail(
                path,
                f"variant {index} ({variant.identifier!r}): duplicate id; "
                f"first defined by variant {first_index}",
            )
        identifiers[variant.identifier] = index

        junction = (
            variant.chromosome,
            variant.intron_start,
            variant.intron_end,
            variant.strand,
        )
        if junction in junctions:
            first_index, first_identifier = junctions[junction]
            _fail(
                path,
                f"variant {index} ({variant.identifier!r}): duplicate defining junction; "
                f"first defined by variant {first_index} ({first_identifier!r})",
            )
        junctions[junction] = (index, variant.identifier)


def load_catalog(path: Path) -> Catalog:
    """Load and validate a schema-version-1 catalog from *path*."""
    try:
        with path.open("rb") as stream:
            document = tomllib.load(stream)
    except tomllib.TOMLDecodeError as error:
        _fail(path, f"invalid TOML: {error}")
    except UnicodeDecodeError as error:
        _fail(path, f"invalid UTF-8: {error}")
    except OSError as error:
        detail = error.strerror or str(error)
        _fail(path, f"cannot read catalog: {detail}")

    _validate_fields(path, document, _TOP_LEVEL_FIELDS, "catalog")

    schema_version = document["schema_version"]
    if type(schema_version) is not int:
        _fail(path, "catalog: field 'schema_version' must be an integer")
    if schema_version != 1:
        _fail(path, f"unsupported catalog schema version {schema_version}; expected 1")

    genome_assembly = _validate_string(
        path, document["genome_assembly"], "genome_assembly", "catalog"
    )
    raw_variants = document["variants"]
    if not isinstance(raw_variants, list):
        _fail(path, "catalog: field 'variants' must be an array of tables")
    if not raw_variants:
        _fail(path, "catalog: field 'variants' must not be empty")

    variants = tuple(
        _parse_variant(path, variant, index)
        for index, variant in enumerate(raw_variants, start=1)
    )
    _validate_unique(path, variants)
    return Catalog(genome_assembly=genome_assembly, variants=variants)
