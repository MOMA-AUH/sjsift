"""Run against the installed wheel/conda package, without repository imports."""
from pathlib import Path
import subprocess
import tempfile

from alignment_fixtures import bam, cli_inputs, cram, fasta, record


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        path = bam(root, [record("installed-smoke", tags={"NH": 1, "XS": "+"})])
        reference = fasta(root)
        cram_path = cram(root, path, reference)
        results = []
        for alignment, options in ((path, []), (cram_path, ["--reference", str(reference)])):
            output = root / f"{alignment.suffix[1:]}.html"
            result = subprocess.run([
                "sjsift", *cli_inputs(root), "--alignments", str(alignment),
                "--html-output", str(output), *options,
            ], capture_output=True, text=True, check=True)
            assert "skip\tsynthetic\tchr1\t101\t200\t+\t8\t2\t10" in result.stdout
            results.append(result.stdout)
            text = output.read_text()
            assert "installed-smoke" in text and "Eligible: 1" in text
            assert "Currently shown: 1" in text
            assert "CIGAR: 20M100N20M" in text
            assert "ACGTACGTACGT" not in text and "AAAAAAAAAAAA" not in text
        assert results[0] == results[1]
        print("Installed BAM and explicit-reference CRAM reports passed")


if __name__ == "__main__":
    main()
