"""Match catalog definitions against STAR junction support."""

from dataclasses import dataclass
from pathlib import Path

from .catalog import Catalog, VariantDefinition


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


def quantify(catalog: Catalog, junctions_path: Path) -> tuple[VariantSupport, ...]:
    """Copy support from exact defining-junction matches in a STAR file."""
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

    with junctions_path.open(encoding="utf-8", newline="") as stream:
        for line in stream:
            fields = line.rstrip("\r\n").split("\t")
            key = (fields[0], int(fields[1]), int(fields[2]), int(fields[3]))
            if key in targets:
                counts[targets[key]] = (int(fields[6]), int(fields[7]))

    return tuple(
        VariantSupport(variant, unique, multimapping)
        for variant, (unique, multimapping) in zip(catalog.variants, counts)
    )
