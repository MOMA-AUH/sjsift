"""Generate a real report and exercise it with the optional Playwright check.

Run with PYTHONPATH=src:tests and NODE_PATH pointing to an installation of
playwright. No browser dependencies are needed by the normal pytest suite.
"""
from dataclasses import replace
import io
from pathlib import Path
import subprocess
import sys

from alignment_fixtures import ANNOTATION, CATALOG, bam, record
from sjsift.alignment_evidence import extract_evidence
from sjsift.catalog import ReferenceJunction
from sjsift.html_report import write_html
from sjsift.quantify import ReferenceJunctionSupport, VariantSupport
from test_paired_view import mates


def report(destination):
    destination.mkdir(parents=True, exist_ok=True)
    reference = ReferenceJunction("next", "chr1", 301, 400, "+", annotation=ANNOTATION)
    variant = replace(CATALOG.variants[0], reference_junctions=(reference,))
    catalog = replace(CATALOG, variants=(variant,))
    # Enough pairs for two expanded pages, plus styling and filter examples.
    reads = [read for i in range(62) for read in mates(f"pair-{i:02}")]
    # One-base terminal blocks must still have legible direction arrows.
    reads[-2].cigarstring = "30M100N19M1X"
    reads[-1].cigarstring = "1=19M100N30M"
    reads += [record("multi", tags={"NH": 2}), record("unknown"),
              record("secondary", flag=256, tags={"NH": 2})]
    # One pair belongs to different support groups and has separated geometry.
    first = record("separated", flag=99, tags={"NH": 1})
    second = record("separated", start=223, cigar="77M100N20M", flag=147, tags={"NH": 1})
    first.next_reference_id = second.next_reference_id = 0
    first.next_reference_start, second.next_reference_start = 223, 80
    reads += [first, second]
    evidence = extract_evidence(catalog, bam(destination, reads))
    output = io.StringIO()
    write_html("synthetic", (VariantSupport(variant, 8, 2, (
        ReferenceJunctionSupport(reference, 20, 3),)),), output, alignment_evidence=evidence)
    path = destination / "paired-view.html"
    path.write_text(output.getvalue())
    return path


if __name__ == "__main__":
    path = report(Path(sys.argv[1]))
    subprocess.run(["node", str(Path(__file__).with_suffix(".cjs")), str(path)], check=True)
