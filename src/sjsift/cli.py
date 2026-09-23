"""Command-line entry point for sjsift."""

import argparse
from collections.abc import Sequence
from pathlib import Path
import sys

from . import __version__
from .catalog import CatalogError, load_catalog
from .quantify import QuantifyError, quantify
from .report import write_tsv


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


def main(argv: Sequence[str] | None = None) -> int:
    """Run the sjsift command-line interface."""
    parser = _parser()
    arguments = parser.parse_args(argv)

    try:
        catalog = load_catalog(Path(arguments.definitions))
    except CatalogError as error:
        parser.error(str(error))
    try:
        results = quantify(catalog, Path(arguments.junctions))
    except QuantifyError as error:
        parser.error(str(error))
    write_tsv(catalog.genome_assembly, results, sys.stdout)
    return 0
