"""Serialize quantification results as the fixed TSV report."""

import csv
from collections.abc import Iterable
from pathlib import Path
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


class ReportError(OSError):
    """A user-correctable problem writing a report destination."""


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


def write_report(
    genome_assembly: str,
    results: Iterable[VariantSupport],
    output_path: Path | None,
    stdout: TextIO,
) -> None:
    """Write a report to standard output or a newly created file."""
    output_created = False
    try:
        if output_path is None:
            write_tsv(genome_assembly, results, stdout)
        else:
            with output_path.open("x", encoding="utf-8", newline="") as stream:
                output_created = True
                write_tsv(genome_assembly, results, stream)
    except FileExistsError:
        raise ReportError(f"{output_path}: output path already exists") from None
    except OSError as error:
        if output_created and output_path is not None:
            try:
                output_path.unlink()
            except OSError:
                pass
        destination = str(output_path) if output_path is not None else "standard output"
        detail = error.strerror or str(error)
        raise ReportError(f"{destination}: cannot write report: {detail}") from None
