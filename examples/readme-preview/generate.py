"""Generate the README's synthetic METex14 report through the real CLI."""

import argparse
from pathlib import Path
import subprocess
import sys
import tomllib

import pysam


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path, help="New directory for generated files")
    args = parser.parse_args()
    args.destination.mkdir(parents=True, exist_ok=False)

    catalog_path = Path(__file__).resolve().parents[2] / "definitions/grch38.toml"
    catalog = tomllib.loads(catalog_path.read_text(encoding="utf-8"))
    variant = next(item for item in catalog["variants"] if item["id"] == "METex14")
    junctions = [variant, *variant["reference_junctions"]]
    reads = []
    star_rows = []
    for junction, (unique_count, multi_count) in zip(
        junctions, ((6, 0), (59, 1), (240, 0)), strict=True
    ):
        start, end = junction["intron_start"], junction["intron_end"]
        role = junction.get("role", "variant")
        # These are invented STAR-format counts, not output from a STAR run.
        star_rows.append(
            f"chr7\t{start}\t{end}\t1\t1\t0\t{unique_count}\t{multi_count}\t75\n"
        )
        for index in range(unique_count + multi_count):
            is_multi = bool(multi_count and index == 8)
            if role == "variant":
                left = 45 + index * 5
                right = 150 - left
            elif index < 10:
                # Keep varied anchors in the initial ten-row preview. Later
                # records start closer to the junction and appear on paging.
                left = 110 - index * 6
                right = (65, 90, 45, 80, 60, 105, 70, 85, 55, 95)[index]
            else:
                left = 25 + (index * 7) % 30
                right = 35 + (index * 17) % 71
            if is_multi:
                right = 20
            length = left + right
            read = pysam.AlignedSegment()
            read.query_name = f"synthetic-{role}-{index + 1:02}"
            read.reference_id = 0
            read.reference_start = start - 1 - left
            read.flag = 16 if index % 2 else 0
            read.mapping_quality = 3 if is_multi else 255
            read.cigarstring = f"{left}M{end - start + 1}N{right}M"
            read.query_sequence = ("ACGT" * ((length + 3) // 4))[:length]
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            read.set_tag("NH", 2 if is_multi else 1)
            read.set_tag("XS", "+", "A")
            read.set_tag("RG", "synthetic")
            reads.append(read)

    bam_path = args.destination / "synthetic.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "chr7", "LN": 159345973},
            {"SN": "chrX", "LN": 156040895},
            {"SN": "chr10", "LN": 133797422},
        ],
        "RG": [{"ID": "synthetic", "SM": "Synthetic example"}],
    }
    with pysam.AlignmentFile(bam_path, "wb", header=header) as stream:
        for read in sorted(reads, key=lambda item: (item.reference_start, item.query_name)):
            stream.write(read)
    pysam.index(str(bam_path))
    star_path = args.destination / "synthetic.SJ.out.tab"
    star_path.write_text("".join(star_rows), encoding="utf-8")
    subprocess.run([
        sys.executable, "-m", "sjsift",
        "--sample-name", "Synthetic example",
        "--junctions", str(star_path),
        "--definitions", str(catalog_path),
        "--alignment", str(bam_path),
        "--output", str(args.destination / "synthetic.tsv"),
        "--context-output", str(args.destination / "synthetic.context.tsv"),
        "--html-output", str(args.destination / "synthetic.html"),
    ], check=True)


if __name__ == "__main__":
    main()
