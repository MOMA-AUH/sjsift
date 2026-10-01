"""Load and validate variant-definition catalogs."""

from dataclasses import dataclass, field
from pathlib import Path
import tomllib
from typing import NoReturn


_TOP_LEVEL_FIELDS = frozenset({"schema_version", "genome_assembly", "variants"})
_ANNOTATION_FIELDS = frozenset(
    {"reference_transcript", "donor_exon", "acceptor_exon", "annotation_source"}
)
_VARIANT_FIELDS = frozenset(
    {"id", "chromosome", "intron_start", "intron_end", "strand", "reference_junctions"}
) | _ANNOTATION_FIELDS
_REFERENCE_JUNCTION_FIELDS = frozenset(
    {"role", "chromosome", "intron_start", "intron_end", "strand"}
) | _ANNOTATION_FIELDS


class CatalogError(ValueError):
    """A user-correctable problem with a variant-definition catalog."""


@dataclass(frozen=True)
class JunctionAnnotation:
    """Explicit labels and their comparison basis, not a whole transcript model."""

    reference_transcript: str
    donor_exon: str
    acceptor_exon: str
    annotation_source: str


@dataclass(frozen=True)
class VariantDefinition:
    """One known splice variant and its defining junction."""

    identifier: str
    chromosome: str
    intron_start: int
    intron_end: int
    strand: str
    annotation: JunctionAnnotation = field(kw_only=True)
    reference_junctions: tuple["ReferenceJunction", ...] = ()


@dataclass(frozen=True)
class ReferenceJunction:
    """A named junction used to place a defining junction in local context."""

    role: str
    chromosome: str
    intron_start: int
    intron_end: int
    strand: str
    annotation: JunctionAnnotation = field(kw_only=True)


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
    if not value.strip():
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


def _parse_annotation(path: Path, value: dict[str, object], context: str) -> JunctionAnnotation:
    return JunctionAnnotation(**{
        name: _validate_string(path, value[name], name, context)
        for name in sorted(_ANNOTATION_FIELDS)
    })


def _parse_reference_junction(
    path: Path, value: object, variant_context: str, index: int
) -> ReferenceJunction:
    context = f"{variant_context}, reference junction {index}"
    if not isinstance(value, dict):
        _fail(path, f"{context}: must be a table")

    _validate_fields(path, value, _REFERENCE_JUNCTION_FIELDS, context)
    role = _validate_string(path, value["role"], "role", context)
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

    return ReferenceJunction(
        role, chromosome, intron_start, intron_end, strand,
        annotation=_parse_annotation(path, value, context),
    )


def _parse_reference_junctions(
    path: Path, value: object, variant_context: str
) -> tuple[ReferenceJunction, ...]:
    if not isinstance(value, list):
        _fail(
            path,
            f"{variant_context}: field 'reference_junctions' must be an array of tables",
        )

    junctions = tuple(
        _parse_reference_junction(path, junction, variant_context, index)
        for index, junction in enumerate(value, start=1)
    )
    roles: set[str] = set()
    for index, junction in enumerate(junctions, start=1):
        if junction.role in roles:
            _fail(
                path,
                f"{variant_context}, reference junction {index} ({junction.role!r}): "
                "duplicate role",
            )
        roles.add(junction.role)
    return junctions


def _parse_variant(
    path: Path, value: object, index: int
) -> VariantDefinition:
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

    reference_junctions = _parse_reference_junctions(path, value["reference_junctions"], context)
    defining_junction = (chromosome, intron_start, intron_end, strand)
    for reference_index, reference in enumerate(reference_junctions, start=1):
        reference_junction = (
            reference.chromosome,
            reference.intron_start,
            reference.intron_end,
            reference.strand,
        )
        if reference_junction == defining_junction:
            _fail(
                path,
                f"{context}, reference junction {reference_index} ({reference.role!r}): "
                "must not duplicate the defining junction",
            )

    return VariantDefinition(
        identifier=identifier,
        chromosome=chromosome,
        intron_start=intron_start,
        intron_end=intron_end,
        strand=strand,
        reference_junctions=reference_junctions,
        annotation=_parse_annotation(path, value, context),
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
    """Load and validate a schema-version-3 catalog."""
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
    if schema_version != 3:
        _fail(path, f"unsupported catalog schema version {schema_version}; expected 3")

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
