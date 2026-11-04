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


def shared_strand_report(destination):
    references = (ReferenceJunction("same_donor", "chr1", 101, 200, "+", annotation=ANNOTATION),
                  ReferenceJunction("same_acceptor", "chr1", 221, 320, "+", annotation=ANNOTATION))
    variant = replace(CATALOG.variants[0], intron_end=320, reference_junctions=references)
    reordered = replace(variant, identifier="reordered", intron_start=401, intron_end=500,
                        reference_junctions=references[::-1])
    catalog = replace(CATALOG, variants=(variant, reordered))
    reads = mates("mixed-pair", tags={"NH": 2, "jM": [21, 20]})
    for read in reads:
        read.cigarstring = f"{100 - read.reference_start}M100N20M100N20M"
        read.query_sequence = "A" * (140 - read.reference_start)
    evidence = extract_evidence(catalog, bam(destination, reads, name="shared.bam"))
    output = io.StringIO()
    write_html("synthetic", tuple(VariantSupport(definition, 8, 2, tuple(
        ReferenceJunctionSupport(reference, 20, 3) for reference in definition.reference_junctions
    )) for definition in catalog.variants), output, alignment_evidence=evidence)
    path = destination / "shared-strand.html"
    path.write_text(output.getvalue())
    return path


def mapping_legend_report(destination):
    references = tuple(ReferenceJunction(role, "chr1", start, start + 99, "+", annotation=ANNOTATION)
                       for role, start in (("same_donor", 301), ("same_acceptor", 501),
                                           ("custom A", 701), ("custom B", 901)))
    variant = replace(CATALOG.variants[0], reference_junctions=references)
    catalog = replace(CATALOG, variants=(variant,))
    reads = [record(f"{start}-{name}", start=start - 21, cigar="3H10S20M100N20M15S2H",
                    flag=flag, tags=tags)
             for start in (101, 301, 501, 701, 901)
             for name, flag, tags in (("unique", 0, {"NH": 1}), ("multi", 0, {"NH": 2}),
                                      ("unknown", 0, {}), ("secondary", 256, {"NH": 2}),
                                      ("supplementary", 2048, {"NH": 1}),
                                      ("duplicate", 1024, {"NH": 2}), ("qc", 512, {"NH": 1}))]
    evidence = extract_evidence(catalog, bam(destination, reads, name="legend.bam"))
    output = io.StringIO()
    write_html("synthetic", (VariantSupport(variant, 8, 2, tuple(
        ReferenceJunctionSupport(reference, 20, 3) for reference in references
    )),), output, alignment_evidence=evidence)
    path = destination / "mapping-legend.html"
    path.write_text(output.getvalue())
    return path


if __name__ == "__main__":
    destination = Path(sys.argv[1])
    path = report(destination)
    shared = shared_strand_report(destination)
    legend = mapping_legend_report(destination)
    subprocess.run(["node", str(Path(__file__).with_suffix(".cjs")), str(path), str(shared), str(legend)], check=True)
