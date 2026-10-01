"""Run against the installed wheel/conda package, without repository imports."""
from pathlib import Path
import subprocess
import tempfile

from alignment_fixtures import bam, cli_inputs, record


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        path = bam(root, [record("installed-smoke", tags={"NH": 1, "XS": "+"})])
        output = root / "report.html"
        result = subprocess.run([
            "sjsift", *cli_inputs(root), "--alignments", str(path),
            "--html-output", str(output),
        ], capture_output=True, text=True, check=True)
        assert "skip\tsynthetic\tchr1\t101\t200\t+\t8\t2\t10" in result.stdout
        text = output.read_text()
        assert "installed-smoke" in text and "Eligible: 1" in text
        assert "Currently shown: 1" in text
        assert "CIGAR: 20M100N20M" in text
        assert "ACGTACGTACGT" not in text


if __name__ == "__main__":
    main()
