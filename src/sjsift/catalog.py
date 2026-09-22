"""Load variant-definition catalogs."""

from dataclasses import dataclass
from pathlib import Path
import tomllib


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


def load_catalog(path: Path) -> Catalog:
    """Load a schema-version-1 catalog from *path*."""
    with path.open("rb") as stream:
        document = tomllib.load(stream)

    if document["schema_version"] != 1:
        raise ValueError("unsupported catalog schema version")

    variants = tuple(
        VariantDefinition(
            identifier=variant["id"],
            chromosome=variant["chromosome"],
            intron_start=variant["intron_start"],
            intron_end=variant["intron_end"],
            strand=variant["strand"],
        )
        for variant in document["variants"]
    )
    return Catalog(genome_assembly=document["genome_assembly"], variants=variants)
