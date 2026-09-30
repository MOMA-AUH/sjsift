"""Command-line entry point for sjsift."""

import argparse
from collections.abc import Sequence
from pathlib import Path
import sys

from . import __version__
from .catalog import CatalogError, load_catalog
from .quantify import QuantifyError, quantify
from .report import ReportError, write_reports


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sjsift",
        description="Quantify predefined RNA splice variants from a STAR junction file.",
    )
    parser.add_argument(
        "--junctions",
        metavar="PATH",
        required=True,
        help="STAR SJ.out.tab file (plain text or gzip-compressed)",
    )
    parser.add_argument(
        "--definitions",
        metavar="PATH",
        required=True,
        help="TOML variant-definition catalog",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="write TSV output to PATH instead of standard output",
    )
    parser.add_argument(
        "--context-output",
        metavar="PATH",
        help="write reference-junction context TSV to PATH",
    )
    parser.add_argument(
        "--html-output",
        metavar="PATH",
        help="write an offline HTML overview and variant details to PATH",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def _run(parser: argparse.ArgumentParser, arguments: argparse.Namespace) -> int:
    try:
        catalog = load_catalog(Path(arguments.definitions))
    except CatalogError as error:
        parser.error(str(error))
    try:
        quantification = quantify(catalog, Path(arguments.junctions))
    except QuantifyError as error:
        parser.error(str(error))
    try:
        write_reports(
            catalog.genome_assembly,
            quantification.results,
            Path(arguments.output) if arguments.output is not None else None,
            Path(arguments.context_output) if arguments.context_output is not None else None,
            sys.stdout,
            html_output_path=Path(arguments.html_output) if arguments.html_output else None,
            junctions_name=Path(arguments.junctions).name,
            definitions_name=Path(arguments.definitions).name,
            compatibility_warning=quantification.compatibility_warning,
        )
    except ReportError as error:
        parser.error(str(error))
    if quantification.compatibility_warning:
        print(
            f"sjsift: warning: {arguments.junctions}: no catalog chromosome "
            "identifiers found; check genome assembly and chromosome naming "
            "compatibility",
            file=sys.stderr,
        )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Run the sjsift command-line interface."""
    parser = _parser()
    arguments = parser.parse_args(argv)
    try:
        return _run(parser, arguments)
    except Exception as error:
        print(f"{parser.prog}: internal error: {error}", file=sys.stderr)
        return 1
