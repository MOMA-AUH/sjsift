"""Isolated fixture, CLI profiling, and retained-evidence verification workers."""
from collections import Counter
import cProfile
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import pstats
import resource
import sys
import time

import pysam


def generate(directory, count, scenario, reference_bases):
    started = time.perf_counter()
    reference = directory / 'reference.fa'
    with reference.open('w') as out:
        out.write('>chr1\n')
        for start in range(0, reference_bases, 80):
            out.write(('ACGT' * 20)[:min(80, reference_bases - start)] + '\n')
    pysam.faidx(str(reference))
    annotation = ('reference_transcript = "NM_synthetic.1"\ndonor_exon = "1"\n'
                  'acceptor_exon = "3"\nannotation_source = "Capacity fixture"\n')
    catalog = 'schema_version = 3\ngenome_assembly = "synthetic"\n'
    for identifier, start in [('abundant', 101)] + ([('rare', 501)] if scenario == 'multiple' else []):
        catalog += (f'[[variants]]\nid = "{identifier}"\nchromosome = "chr1"\n'
                    f'intron_start = {start}\nintron_end = {start + 99}\nstrand = "+"\n' + annotation)
        if scenario == 'multiple':
            for role, position in [('shared', 221)] + ([('abundant', 101)] if identifier == 'rare' else []):
                catalog += (f'[[variants.reference_junctions]]\nrole = "{role}"\nchromosome = "chr1"\n'
                            f'intron_start = {position}\nintron_end = {position + 99}\nstrand = "+"\n' + annotation)
        else:
            catalog += 'reference_junctions = []\n'
    (directory / 'catalog.toml').write_text(catalog)
    # Deliberately independent STAR counts: viewing must not recount alignments.
    (directory / 'sample.SJ.out.tab').write_text('chr1\t101\t200\t1\t1\t0\t800\t200\t20\n'
        + ('chr1\t221\t320\t1\t1\t0\t300\t100\t20\nchr1\t501\t600\t1\t1\t0\t1\t1\t20\n'
           if scenario == 'multiple' else ''))
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': 'chr1', 'LN': reference_bases}],
              'RG': [{'ID': 'benchmark', 'SM': 'synthetic'}]}
    cigars = ['20M100N20M', '5S20M100N20M5S', '10M1I10M100N20M',
              '10M1D9M100N20M', '10=10X100N20M']
    with pysam.AlignmentFile(str(directory / 'reads.bam'), 'wb', header=header) as output:
        for i in range(count):
            rare = scenario == 'multiple' and i >= count - 2
            read = pysam.AlignedSegment(output.header)
            read.query_name = 'rare' if rare else f'pair-{i // 2:06}'
            read.flag = 99 if i % 2 == 0 else 147
            read.reference_id = read.next_reference_id = 0
            read.reference_start = read.next_reference_start = 480 if rare else 80
            read.mapping_quality = 60
            read.cigarstring = cigars[(i // 2) % len(cigars)] + (
                '100N20M' if scenario == 'multiple' and i < count // 2 else '')
            length = sum(size for op, size in read.cigartuples if op in (0, 1, 4, 7, 8))
            sequence = ('ACGT' * ((length + 3) // 4))[:length]
            if (i // 2) % len(cigars) == 4:
                # '=' bases agree with reference; 'X' bases must actually differ.
                sequence = sequence[:10] + sequence[10:20].translate(str.maketrans('ACGT', 'TGCA')) + sequence[20:]
            read.query_sequence = sequence
            read.query_qualities = pysam.qualitystring_to_array('I' * length)
            read.template_length = read.reference_length * (1 if i % 2 == 0 else -1)
            read.set_tag('RG', 'benchmark')
            read.set_tag('XS', '+', 'A')
            if (i // 2) % 5 != 4:
                read.set_tag('NH', 2 if (i // 2) % 5 == 3 else 1)
            output.write(read)
    pysam.index(str(directory / 'reads.bam'))
    bam_seconds = time.perf_counter() - started
    started = time.perf_counter()
    with pysam.AlignmentFile(str(directory / 'reads.bam'), 'rb') as source:
        with pysam.AlignmentFile(str(directory / 'reads.cram'), 'wc', header=source.header,
                                 reference_filename=str(reference)) as output:
            for read in source:
                output.write(read)
    pysam.index(str(directory / 'reads.cram'))
    return {'bam_and_inputs_seconds': bam_seconds, 'cram_conversion_seconds': time.perf_counter() - started,
            'reference_bases': reference_bases, 'reference_bytes': reference.stat().st_size,
            'bam_bytes': (directory / 'reads.bam').stat().st_size,
            'cram_bytes': (directory / 'reads.cram').stat().st_size}


def report(directory, format_name):
    from sjsift.cli import main
    args = ['--definitions', str(directory / 'catalog.toml'), '--junctions', str(directory / 'sample.SJ.out.tab'),
            '--alignment', str(directory / f'reads.{format_name}'),
            '--output', str(directory / f'{format_name}.tsv'),
            '--context-output', str(directory / f'{format_name}.context.tsv'),
            '--html-output', str(directory / f'{format_name}.html')]
    if format_name == 'cram':
        args += ['--reference', str(directory / 'reference.fa')]
    profile = cProfile.Profile()
    code = profile.runcall(main, args)
    profile.dump_stats(str(directory / f'{format_name}.prof'))
    if code:
        raise RuntimeError(f'sjsift exited {code}')
    stats = pstats.Stats(profile).stats
    def elapsed(filename, name):
        return sum(value[3] for key, value in stats.items() if key[0].endswith(filename) and key[2] == name)
    extraction = elapsed('alignment_evidence.py', 'extract_evidence')
    validation = elapsed('alignment_evidence.py', '_validate_reference')
    return {'extraction_seconds': extraction, 'reference_validation_seconds': validation,
            'extraction_without_reference_validation_seconds': extraction - validation,
            'render_seconds': elapsed('html_report.py', 'write_html'),
            'html_bytes': (directory / f'{format_name}.html').stat().st_size,
            'peak_process_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024}


class Inventory(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = 0
        self.memberships = 0
        self.ids = set()
        self.labels = Counter()
        self.variant = None
        self.dom_elements = 0

    def handle_starttag(self, tag, attributes):
        self.dom_elements += 1
        attrs = dict(attributes)
        classes = attrs.get('class', '').split()
        if 'variant-detail' in classes:
            self.variant = attrs['id']
        if 'read-row' in classes:
            self.rows += 1
            self.memberships += len(attrs['data-support'].split())
            self.ids.add((self.variant, attrs['data-record-id']))
        if 'read-hit' in classes:
            self.labels[(self.variant, attrs['title'].split(';')[0])] += 1


def validate(directory, format_name, count, scenario):
    inventory = Inventory()
    path = directory / f'{format_name}.html'
    with path.open() as source:
        while chunk := source.read(1024 * 1024):
            inventory.feed(chunk)
    expected = Counter()
    for i in range(count):
        rare = scenario == 'multiple' and i >= count - 2
        name = 'rare' if rare else f'pair-{i // 2:06}'
        if not rare:
            expected[('variant-0', name)] += 1
        if scenario == 'multiple':
            expected[('variant-1', name)] += 1
    rows = count if scenario == 'single' else 2 * count - 2
    memberships = count if scenario == 'single' else 3 * count - 2
    assert inventory.rows == len(inventory.ids) == rows, (inventory.rows, rows)
    assert inventory.memberships == memberships, (inventory.memberships, memberships)
    assert inventory.labels == expected, 'Missing or repeated eligible record names'
    import csv
    with (directory / f'{format_name}.tsv').open() as source:
        results = list(csv.DictReader(source, delimiter='\t'))
    assert [(row['variant_id'], row['unique_support'], row['multimapping_support']) for row in results] == (
        [('abundant', '800', '200')] + ([('rare', '1', '1')] if scenario == 'multiple' else []))
    if format_name == 'cram':
        for suffix in ['tsv', 'context.tsv']:
            assert (directory / f'bam.{suffix}').read_bytes() == (directory / f'cram.{suffix}').read_bytes()
    with path.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    return {'eligible_records': count, 'embedded_rows': rows, 'memberships': memberships,
            'dom_elements': inventory.dom_elements,
            'html_sha256': digest}


if __name__ == '__main__':
    action, target, destination, *arguments = sys.argv[1:]
    directory = Path(destination)
    if action == 'generate':
        result = generate(directory, int(arguments[0]), arguments[1], int(arguments[2]))
    elif action == 'report':
        result = report(directory, arguments[0])
    else:
        result = validate(directory, arguments[0], int(arguments[1]), arguments[2])
    Path(target).write_text(json.dumps(result, indent=2) + '\n')
