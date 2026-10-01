"""Small real indexed alignments for behavior and installed-package checks."""
from array import array
from pathlib import Path

import pysam

from sjsift.catalog import Catalog, JunctionAnnotation, VariantDefinition


ANNOTATION = JunctionAnnotation("NM_synthetic.1", "1", "3", "Synthetic fixture")
CATALOG = Catalog("synthetic", (VariantDefinition(
    "skip", "chr1", 101, 200, "+", annotation=ANNOTATION,
),))


def record(name="read", start=80, cigar="20M100N20M", flag=0, mapq=60, tags=None, contig=0):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.flag = flag
    read.reference_id = contig
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigarstring = cigar
    length = sum(n for op, n in read.cigartuples if op in (0, 1, 4, 7, 8))
    read.query_sequence = "ACGT" * (length // 4) + "ACGT"[:length % 4]
    read.query_qualities = pysam.qualitystring_to_array("I" * length)
    for tag, value in (tags or {}).items():
        if tag in ("TS", "XS", "ts") and isinstance(value, str) and len(value) == 1:
            read.set_tag(tag, value, "A")
        elif tag in ("jM", "jI") and isinstance(value, list):
            read.set_tag(tag, array("i", value))
        else:
            read.set_tag(tag, value)
    return read


def bam(tmp_path, reads, header=None, name="reads.bam", index=True):
    path = tmp_path / name
    header = header or {"HD": {"VN": "1.6", "SO": "coordinate"},
                        "SQ": [{"SN": "chr1", "LN": 2000}],
                        "PG": [{"ID": "STAR", "PN": "STAR", "VN": "2.7.11b"}]}
    with pysam.AlignmentFile(path, "wb", header=header) as stream:
        for read in sorted(reads, key=lambda r: (r.reference_id, r.reference_start)):
            stream.write(read)
    if index:
        pysam.index(str(path))
    return path


def cli_inputs(tmp_path):
    catalog = tmp_path / "catalog.toml"
    catalog.write_text('''schema_version = 3
genome_assembly = "synthetic"
[[variants]]
id = "skip"
chromosome = "chr1"
intron_start = 101
intron_end = 200
strand = "+"
reference_junctions = []
reference_transcript = "NM_synthetic.1"
donor_exon = "1"
acceptor_exon = "3"
annotation_source = "Synthetic fixture"
''')
    sj = tmp_path / "sample.SJ.out.tab"
    sj.write_text("chr1\t101\t200\t1\t1\t0\t8\t2\t20\n")
    return ["--junctions", str(sj), "--definitions", str(catalog)]


def fasta(tmp_path, sequence="A" * 2000, name="reference.fa", contig="chr1"):
    path = tmp_path / name
    path.write_text(f">{contig}\n{sequence}\n")
    pysam.faidx(str(path))
    return path


def cram(tmp_path, bam_path, reference_path, options=None):
    path = tmp_path / "reads.cram"
    with pysam.AlignmentFile(str(bam_path), "rb") as source:
        with pysam.AlignmentFile(str(path), "wc", header=source.header,
                                 reference_filename=str(reference_path),
                                 format_options=options) as output:
            for read in source:
                output.write(read)
    pysam.index(str(path))
    return path
