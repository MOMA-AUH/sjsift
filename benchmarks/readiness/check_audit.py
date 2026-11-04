"""Prove that the independent audit rejects wrong counts and substituted mates."""
import argparse
import csv
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path)
    parser.add_argument('report', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output already exists')
    original_html = (args.report / 'report.html').read_text()
    rows = re.findall(r'<li class="read-row.*?</li>', original_html, re.S)
    first = next(row for row in rows if '; read 1,' in row)
    name = re.search(r'title="([^;]+);', first)[1]
    second = next(row for row in rows if f'title="{name};' in row and '; read 2,' in row)
    # Both mates have the same legacy name/class/membership key. The physical
    # identity comparison must distinguish them even when that weak key cannot.
    for attribute in ['data-mapping', 'data-support']:
        pattern = rf'{attribute}="([^"]*)"'
        assert re.search(pattern, first)[1] == re.search(pattern, second)[1]
    pair_id = re.search(r'data-pair="([^"]+)"', first)[1]
    other_pair_id = next(re.search(r'data-pair="([^"]+)"', row)[1] for row in rows
                         if 'data-pair="' in row and re.search(r'data-pair="([^"]*)"', row)[1] not in ['', pair_id])
    wrong_pair = first.replace(f'data-pair="{pair_id}"', f'data-pair="{other_pair_id}"', 1)
    table = list(csv.reader(io.StringIO((args.report / 'main.tsv').read_text()), delimiter='\t'))
    for column in ['unique_support', 'total_support']:
        index = table[0].index(column)
        table[1][index] = str(int(table[1][index]) + 1)
    buffer = io.StringIO()
    csv.writer(buffer, delimiter='\t', lineterminator='\n').writerows(table)
    provenance = json.loads((args.report / 'provenance.json').read_text())
    provenance['schema_version'] = 2
    mutations = {
        'wrong_star_count': ('main.tsv', buffer.getvalue()),
        'missing_physical_record': ('report.html', original_html.replace(first, '', 1)),
        'wrong_pair_association': ('report.html', original_html.replace(first, wrong_pair, 1)),
        'same_name_mate_substitution': ('report.html', original_html.replace(second, first, 1)),
        'wrong_provenance_schema': ('provenance.json', json.dumps(provenance)),
    }
    results = {}
    with tempfile.TemporaryDirectory(prefix='sjsift-audit-check-') as temporary:
        root = Path(temporary)
        report = root / 'report'
        def audit():
            return subprocess.run([sys.executable, str(Path(__file__).with_name('audit_evidence.py')),
                str(args.inputs.resolve()), str(report), str(root / 'inventory.json')],
                capture_output=True, text=True)
        shutil.copytree(args.report, report)
        baseline = audit()
        assert baseline.returncode == 0, baseline.stderr
        results['unmodified_report'] = 'passed'
        for name, (filename, contents) in mutations.items():
            path = report / filename
            original = path.read_text()
            assert original != contents
            try:
                path.write_text(contents)
                result = audit()
                assert result.returncode != 0 and 'AssertionError' in result.stderr, name
                results[name] = 'rejected'
            finally:
                path.write_text(original)
        restored = audit()
        assert restored.returncode == 0, restored.stderr
        results['restored_report'] = 'passed'
    args.output.write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
