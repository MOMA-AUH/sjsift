"""Serialize quantification results as the fixed TSV report."""

import csv
from collections.abc import Iterable
from contextlib import ExitStack
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

CONTEXT_HEADER = (
    "variant_id",
    "genome_assembly",
    "context_role",
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


def write_context_tsv(
    genome_assembly: str,
    results: Iterable[VariantSupport],
    stream: TextIO,
) -> None:
    """Write one context-evidence row for every configured reference junction."""
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(CONTEXT_HEADER)
    for result in results:
        for support in result.reference_junctions:
            definition = support.definition
            writer.writerow(
                (
                    result.definition.identifier,
                    genome_assembly,
                    definition.role,
                    definition.chromosome,
                    definition.intron_start,
                    definition.intron_end,
                    definition.strand,
                    support.unique,
                    support.multimapping,
                    support.total,
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


def write_reports(
    genome_assembly: str,
    results: Iterable[VariantSupport],
    output_path: Path | None,
    context_output_path: Path,
    stdout: TextIO,
) -> None:
    """Write the main report and a context report without leaving partial files."""
    if output_path is not None and output_path == context_output_path:
        raise ReportError("main and context output paths must differ")

    result_rows = tuple(results)
    created_paths: list[Path] = []
    completed = False
    try:
        with ExitStack() as stack:
            if output_path is None:
                main_stream = stdout
            else:
                main_stream = output_path.open("x", encoding="utf-8", newline="")
                created_paths.append(output_path)
                stack.enter_context(main_stream)
            context_stream = context_output_path.open("x", encoding="utf-8", newline="")
            created_paths.append(context_output_path)
            stack.enter_context(context_stream)
            write_tsv(genome_assembly, result_rows, main_stream)
            write_context_tsv(genome_assembly, result_rows, context_stream)
        completed = True
    except FileExistsError as error:
        destination = error.filename or context_output_path
        raise ReportError(f"{destination}: output path already exists") from None
    except OSError as error:
        destination = error.filename or context_output_path
        detail = error.strerror or str(error)
        raise ReportError(f"{destination}: cannot write report: {detail}") from None
    finally:
        if not completed:
            for path in created_paths:
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass
                except OSError:
                    pass
