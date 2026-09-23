"""Command-line entry point for sjsift."""

import argparse
from collections.abc import Sequence
from pathlib import Path
import sys

from . import __version__
from .catalog import CatalogError, load_catalog
from .quantify import QuantifyError, quantify
from .report import ReportError, write_report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sjsift",
        description="Quantify predefined RNA splice variants from a STAR junction file.",
    )
    parser.add_argument(
        "--junctions",
        metavar="PATH",
        required=True,
        help="STAR SJ.out.tab file",
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
        write_report(
            catalog.genome_assembly,
            quantification.results,
            Path(arguments.output) if arguments.output is not None else None,
            sys.stdout,
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
