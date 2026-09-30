"""Reproduce the pinned exon evidence used for GRCh38 reference curation."""

import argparse
import csv
import gzip
import hashlib
from pathlib import Path
import re
import sys


SOURCE_SHA256 = "d6e6fe0515c95b2a8cd36a853c1989cee9115c736c60237c56ae92b9daaaf7c4"
TRANSCRIPTS = {
    "ENST00000374690.9": 8,
    "ENST00000646891.2": 18,
    "ENST00000275493.7": 28,
    "ENST00000450046.2": 28,
    "ENST00000344576.7": 16,
    "ENST00000358487.10": 18,
    "ENST00000397752.8": 21,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gtf", type=Path, help="gencode.v49.annotation.gtf.gz")
    path = parser.parse_args().gtf
    with path.open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != SOURCE_SHA256:
            parser.error("GTF checksum differs from the pinned GENCODE v49 source")

    rows = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "exon":
                continue
            transcript = re.search(r'transcript_id "([^"]+)"', fields[8])[1]
            if transcript not in TRANSCRIPTS:
                continue
            attributes = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
            rank = int(re.search(r"exon_number (\d+)", fields[8])[1])
            rows.append((
                attributes["gene_name"], transcript, fields[0], fields[6],
                rank, int(fields[3]), int(fields[4]), attributes["exon_id"],
            ))

    for transcript, count in TRANSCRIPTS.items():
        ranks = sorted(row[4] for row in rows if row[1] == transcript)
        if ranks != list(range(1, count + 1)):
            parser.error(f"unexpected exon ranks for {transcript}")

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow((
        "gene", "transcript_id", "chromosome", "strand", "exon_number",
        "exon_start", "exon_end", "exon_id",
    ))
    writer.writerows(sorted(rows, key=lambda row: (row[0], row[1], row[4])))


if __name__ == "__main__":
    main()
