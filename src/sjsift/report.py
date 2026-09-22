"""Serialize quantification results as the fixed TSV report."""

import csv
from collections.abc import Iterable
from typing import TextIO

from .quantify import VariantSupport


HEADER = (
    "variant_id",
    "genome_assembly",
    "chromosome",
    "intron_start",
    "intron_end",
    "strand",
    "unique_support",
    "multimapping_support",
    "total_support",
)


def write_tsv(
    genome_assembly: str,
    results: Iterable[VariantSupport],
    stream: TextIO,
) -> None:
    """Write *results* to *stream* using the report schema."""
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(HEADER)
    for result in results:
        definition = result.definition
        writer.writerow(
            (
                definition.identifier,
                genome_assembly,
                definition.chromosome,
                definition.intron_start,
                definition.intron_end,
                definition.strand,
                result.unique,
                result.multimapping,
                result.total,
            )
        )
