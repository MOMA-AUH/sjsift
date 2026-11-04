"""Regenerate tiny STAR interoperability inputs/output; never invoked by CI."""
import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import tempfile

import pysam


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def sequence(label):
    # Each digest byte contributes one base; no platform-dependent PRNG.
    data = b"".join(hashlib.sha256(f"{label}:{i}".encode()).digest() for i in range(32))
    return "".join("ACGT"[byte % 4] for byte in data[:1000])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new output directory")
    parser.add_argument("--star", default="STAR")
    args = parser.parse_args()
    assert subprocess.check_output([args.star, "--version"], text=True).strip() == "2.7.11b"
    assert pysam.__version__ == "0.24.1" and pysam.__samtools_version__ == "1.24"
    args.output.mkdir(parents=True, exist_ok=False)
    root = args.output.resolve()
    sequences = {}
    for contig, label, minus in [("plus", "plus", False), ("minus", "minus", True),
                                  ("multiA", "multi", False), ("multiB", "multi", False)]:
        bases = list(sequence(label))
        for start, end in [(200, 300), (340, 440)]:
            bases[start:start + 2] = "CT" if minus else "GT"
            bases[end - 2:end] = "AC" if minus else "AG"
        sequences[contig] = "".join(bases)
    reference = root / "reference.fa"
    reference.write_text("".join(f">{name}\n{bases}\n" for name, bases in sequences.items()))
    pysam.faidx(str(reference))
    with (root / "reads_1.fastq").open("w") as first, (root / "reads_2.fastq").open("w") as second:
        for name, contig in [("plus-1", "plus"), ("plus-2", "plus"), ("minus-1", "minus"), ("multi-1", "multiA")]:
            bases = sequences[contig]
            transcript = bases[100:200] + bases[300:340] + bases[440:540]
            if contig == "minus":
                transcript = reverse_complement(transcript)
            for stream, read in [(first, transcript[60:180]), (second, reverse_complement(transcript[120:240]))]:
                stream.write(f"@{name}\n{read}\n+\n{'I' * len(read)}\n")
    with tempfile.TemporaryDirectory(prefix="sjsift-star-") as directory:
        work = Path(directory)
        index = work / "index"
        index.mkdir()
        subprocess.run([args.star, "--runMode", "genomeGenerate", "--runThreadN", "1",
                        "--genomeDir", str(index), "--genomeFastaFiles", str(reference),
                        "--genomeSAindexNbases", "3", "--genomeChrBinNbits", "10",
                        "--outFileNamePrefix", str(work / "index-")], check=True)
        subprocess.run([args.star, "--runThreadN", "1", "--runRNGseed", "777",
                        "--genomeDir", str(index), "--readFilesIn", str(root / "reads_1.fastq"), str(root / "reads_2.fastq"),
                        "--outFileNamePrefix", str(work / "aligned-"),
                        "--outSAMtype", "SAM", "--outSAMattributes", "NH", "HI", "AS", "nM", "jM", "jI",
                        "--alignIntronMin", "20", "--alignIntronMax", "500",
                        "--alignMatesGapMax", "500", "--outFilterMismatchNmax", "0",
                        "--outSJfilterCountUniqueMin", "1", "1", "1", "1",
                        "--outSJfilterCountTotalMin", "1", "1", "1", "1",
                        "--outSJfilterOverhangMin", "12", "12", "12", "12",
                        "--outSJfilterDistToOtherSJmin", "0", "0", "0", "0"], check=True)
        shutil.copyfile(work / "aligned-SJ.out.tab", root / "SJ.out.tab")
        shutil.copyfile(work / "aligned-Log.final.out", root / "Log.final.out")
        pysam.sort("--no-PG", "-o", str(root / "alignments.bam"), str(work / "aligned-Aligned.out.sam"))
        pysam.index(str(root / "alignments.bam"))
    baseline = Path(__file__).resolve().parent
    for name in ("reference.fa", "reference.fa.fai", "reads_1.fastq", "reads_2.fastq", "SJ.out.tab"):
        assert (root / name).read_bytes() == (baseline / name).read_bytes(), name
    with pysam.AlignmentFile(str(root / "alignments.bam")) as actual:
        with pysam.AlignmentFile(str(baseline / "alignments.bam")) as expected:
            assert [r.to_string() for r in actual] == [r.to_string() for r in expected]
    print("Regeneration matches committed inputs, STAR counts, and complete SAM records.")


if __name__ == "__main__":
    main()
