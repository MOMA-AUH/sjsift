"""Opt-in benchmark contract check; never collected by the fast test suite."""
import json
from pathlib import Path
import subprocess
import sys


def test_small_benchmark_retains_records_and_separates_measurements(tmp_path):
    subprocess.run([sys.executable, str(Path(__file__).with_name('run.py')),
                    str(tmp_path / 'run'), '--sizes', '20', '--skip-browser'], check=True)
    result = json.loads((tmp_path / 'run' / 'results.json').read_text())
    assert result['settings']['sizes'] == [20]
    assert len(result['cases']) == 4
    for case in result['cases']:
        assert case['report']['status'] == 'passed'
        assert case['validation']['status'] == 'passed'
        assert case['validation']['eligible_records'] == 20
        assert case['validation']['embedded_rows'] == (20 if case['scenario'] == 'single' else 38)
        assert case['validation']['memberships'] == (20 if case['scenario'] == 'single' else 58)
        assert case['browser']['status'] == 'skipped'
        assert case['report']['peak_tree_rss_bytes'] > 0
        assert case['report']['metrics']['extraction_seconds'] > 0
        assert case['report']['metrics']['render_seconds'] > 0
        assert case['report']['metrics']['reference_validation_seconds'] >= 0
    assert result['largest_complete_records'] is None


def test_timeout_is_retained_and_dependent_checks_are_incomplete(tmp_path):
    completed = subprocess.run([sys.executable, str(Path(__file__).with_name('run.py')),
        str(tmp_path / 'limited'), '--sizes', '20', '--skip-browser', '--stage-seconds', '0.001'])
    assert completed.returncode == 1
    result = json.loads((tmp_path / 'limited' / 'results.json').read_text())
    assert result['largest_complete_records'] is None
    assert len(result['cases']) == 4
    for case in result['cases']:
        assert case['generation']['status'] == 'timeout'
        assert case['report']['status'] == 'incomplete'
        assert case['validation']['status'] == 'incomplete'
        assert case['browser']['status'] == 'incomplete'
