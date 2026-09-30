"""Serialize TSV reports and manage report destinations."""

import csv
from collections.abc import Iterable
from contextlib import ExitStack
from functools import partial
from pathlib import Path
from typing import TextIO

from .quantify import VariantSupport
from .html_report import write_html


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


def write_reports(
    genome_assembly: str,
    results: Iterable[VariantSupport],
    output_path: Path | None,
    context_output_path: Path | None,
    stdout: TextIO,
    *,
    html_output_path: Path | None = None,
    junctions_name: str = "",
    definitions_name: str = "",
    compatibility_warning: bool = False,
) -> None:
    """Create requested reports exclusively; remove new files if any output fails.

    All file destinations are opened before writing, and standard output is
    deferred until files are flushed and closed. Emitted stdout cannot be undone.
    """
    writers = [(output_path, write_tsv), (context_output_path, write_context_tsv)]
    if html_output_path is not None:
        writers.append((html_output_path, partial(
            write_html,
            junctions_name=junctions_name,
            definitions_name=definitions_name,
            compatibility_warning=compatibility_warning,
        )))
    writers = [(path, writer) for path, writer in writers if path is not None]
    result_rows = tuple(results)
    created_paths: list[Path] = []
    completed = False
    destination = "output paths"
    try:
        paths = [path.resolve() for path, _ in writers]
        if len(set(paths)) != len(paths):
            raise ReportError("main, context, and HTML output paths must differ")
        with ExitStack() as stack:
            streams = []
            for path, writer in writers:
                destination = str(path)
                stream = stack.enter_context(path.open("x", encoding="utf-8", newline=""))
                created_paths.append(path)
                streams.append((path, writer, stream))
            for path, writer, stream in streams:
                destination = str(path)
                writer(genome_assembly, result_rows, stream)
                stream.flush()
        if output_path is None:
            destination = "standard output"
            write_tsv(genome_assembly, result_rows, stdout)
            stdout.flush()
        completed = True
    except ReportError:
        raise
    except FileExistsError as error:
        destination = error.filename or destination
        raise ReportError(f"{destination}: output path already exists") from None
    except OSError as error:
        destination = error.filename or destination
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
