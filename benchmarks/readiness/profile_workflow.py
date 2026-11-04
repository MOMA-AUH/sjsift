"""Measure installed candidate CLI or reference-only validation in a fresh process."""
import argparse
import cProfile
import json
from pathlib import Path
import pstats
import resource
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['report', 'reference'])
    parser.add_argument('inputs', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--format', choices=['bam', 'cram'], default='bam')
    args = parser.parse_args()
    root = args.inputs
    started = time.perf_counter()
    if args.action == 'reference':
        import pysam
        from sjsift.alignment_evidence import _validate_reference
        # Instrumentation only: this private function is not a supported API.
        with pysam.AlignmentFile(str(root / 'reads.cram'), 'rc',
                                 reference_filename=str(root / 'reference.fa')) as source:
            header = source.header.to_dict()
        start = time.perf_counter()
        _, metadata = _validate_reference(header, root / 'reference.fa')
        result = {'reference_validation_seconds': time.perf_counter() - start,
                  'checks': metadata['checks']}
    else:
        from sjsift.cli import main as cli
        args.output.mkdir(exist_ok=False)
        command = ['--definitions', str(root / 'catalog.toml'), '--junctions', str(root / 'sample.SJ.out.tab'),
                   '--alignment', str(root / f'reads.{args.format}'),
                   '--output', str(args.output / 'main.tsv'),
                   '--context-output', str(args.output / 'context.tsv'),
                   '--html-output', str(args.output / 'report.html'),
                   '--provenance-output', str(args.output / 'provenance.json')]
        if args.format == 'cram':
            command += ['--reference', str(root / 'reference.fa')]
        profile = cProfile.Profile()
        assert profile.runcall(cli, command) == 0
        profile.dump_stats(str(args.output / 'report.prof'))
        stats = pstats.Stats(profile).stats
        def elapsed(filename, name):
            return sum(v[3] for k, v in stats.items() if k[0].endswith(filename) and k[2] == name)
        reference = elapsed('alignment_evidence.py', '_validate_reference')
        extraction = elapsed('alignment_evidence.py', 'extract_evidence')
        result = {'reference_validation_seconds': reference, 'extraction_seconds': extraction,
                  'extraction_without_reference_validation_seconds': extraction - reference,
                  'render_seconds': elapsed('html_report.py', 'write_html'),
                  'html_bytes': (args.output / 'report.html').stat().st_size}
    result.update(wall_seconds=time.perf_counter() - started,
                  peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024)
    target = args.output if args.action == 'reference' else args.output / 'profile.json'
    target.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
