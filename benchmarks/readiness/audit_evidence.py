"""Independently join raw STAR/CIGAR evidence to the full catalog and reports.

This audit imports no sjsift matcher or evidence selector. The aggregate output
contains no read names. Exact read inventories remain in memory only.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import tomllib
from typing import NamedTuple

import pysam


def junction(row):
    return row['chromosome'], int(row['intron_start']), int(row['intron_end']), row['strand']


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


class RecordIdentity(NamedTuple):
    name: str
    chromosome: str
    start: int
    end: int
    cigar: str
    mapq: int
    nh: int
    flag: int
    read_group: str
    hit_index: int


class ReportInventory(HTMLParser):
    def __init__(self, path):
        super().__init__()
        self.variant = None
        self.row = None
        self.records = Counter()
        self.statuses = Counter()
        self.metadata_depth = 0
        self.metadata_text = []
        self.name = None
        self.pair = None
        self.pairs = defaultdict(list)
        with path.open() as stream:
            while text := stream.read(1024 * 1024):
                self.feed(text)

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        classes = attrs.get('class', '').split()
        if 'variant-detail' in classes:
            self.variant = int(attrs['id'].split('-')[1])
        if 'read-row' in classes:
            self.row = (self.variant, attrs['data-mapping'], attrs['data-support'])
            self.pair = attrs.get('data-pair')
            self.statuses.update(set(classes) & {'agreeing', 'unverified', 'mixed'})
        if 'read-hit' in classes:
            self.name = attrs['title'].split(';')[0]
        if 'read-metadata' in classes:
            self.metadata_depth = 1
            self.metadata_text = []
        elif self.metadata_depth and tag == 'div':
            self.metadata_depth += 1

    def handle_data(self, text):
        if self.metadata_depth:
            self.metadata_text.append(text)

    def handle_endtag(self, tag):
        if tag != 'div' or not self.metadata_depth:
            return
        self.metadata_depth -= 1
        if self.metadata_depth:
            return
        text = ' '.join(self.metadata_text)
        location = re.search(r'(\S+):(\d+)–(\d+) · CIGAR: (\S+) · MAPQ: (\d+) · NH: (\d+)', text)
        flags = re.search(r'Flags:.*?\((\d+)\)', text)
        group = re.search(r'Read group: (.*?) · HI: (\d+)', text)
        assert location and flags and group, 'missing physical alignment metadata'
        chromosome, start, end, cigar, mapq, nh = location.groups()
        identity = RecordIdentity(self.name, chromosome, int(start), int(end), cigar, int(mapq),
                    int(nh), int(flags[1]), group[1], int(group[2]))
        self.records[(*self.row, identity)] += 1
        if self.pair:
            self.pairs[(self.variant, self.pair)].append(identity)


def inspect_metadata(provenance, inputs, catalog, header, format, sample_name):
    """Compare structured facts with the invocation, raw header, and local FASTA."""
    assert provenance['schema_version'] == 1
    assert provenance['genome_assembly'] == catalog['genome_assembly']
    assert provenance['catalog']['sha256'] == sha256(inputs / 'catalog.toml')
    assert provenance['star_input']['sha256'] == sha256(inputs / 'sample.SJ.out.tab')
    groups = header.get('RG', [])
    samples = {group['SM'] for group in groups if group.get('SM', '').strip()}
    # Readiness inputs have at most one sample; the raw inventory includes all RGs.
    assert len(samples) <= 1, 'readiness audit requires single-sample alignments'
    if sample_name is None:
        assert len(samples) == 1, 'pass the report invocation\'s --sample-name to the audit'
        sample_name = next(iter(samples))
        sample_source = 'alignment_header'
    else:
        assert not samples or sample_name in samples
        sample_source = 'explicit'
    assert provenance['sample']['name'] == sample_name
    assert provenance['sample']['source'] == sample_source
    alignment = provenance['alignment']
    assert alignment['status'] == 'requested'
    assert alignment['format'] == format.upper()
    assert alignment['basename'] == f'reads.{format}'
    expected_selection = {
        'policy': 'all_read_groups',
        'header_sample_count': len(samples),
        'selected_sample_in_header': sample_name in samples,
        'header_read_group_count': len(groups),
        'read_groups_without_sample': sum(not group.get('SM', '').strip() for group in groups),
    }
    for key, expected in expected_selection.items():
        assert alignment['sample_selection'][key] == expected, f'sample_selection.{key}'
    required = {row['chromosome'] for variant in catalog['variants']
                for row in [variant, *variant['reference_junctions']]}
    labels = [sq['AS'] for sq in header['SQ'] if sq['SN'] in required and sq.get('AS')]
    assert all(label == catalog['genome_assembly'] for label in labels)
    expected_assembly = {
        'required_contigs': len(required), 'labeled_contigs': len(labels),
        'labels_match_catalog': True if labels else None,
    }
    for key, expected in expected_assembly.items():
        assert alignment['assembly_checks'][key] == expected, f'assembly_checks.{key}'
    if format == 'bam':
        assert alignment['reference'] is None
        return
    reference = alignment['reference']
    assert reference['basename'] == 'reference.fa'
    assert reference['index_basenames'] == ['reference.fa.fai']
    assert reference['access'] == 'explicit_local_only'
    checked = 0
    with pysam.FastaFile(str(inputs / 'reference.fa')) as fasta:
        for sq in header['SQ']:
            assert fasta.get_reference_length(sq['SN']) == sq['LN']
            digest = hashlib.md5()
            for start in range(0, sq['LN'], 1024 * 1024):
                digest.update(fasta.fetch(sq['SN'], start, min(start + 1024 * 1024, sq['LN'])).upper().encode('ascii'))
            if sq.get('M5'):
                assert digest.hexdigest() == sq['M5'].lower(), 'CRAM M5 differs from local FASTA'
                checked += 1
    expected_checks = {
        'contig_names': True, 'contig_lengths': True, 'full_sequences_read': True,
        'header_contigs': len(header['SQ']), 'header_m5_checked': checked,
        'decoded_slice_checks': True,
    }
    for key, expected in expected_checks.items():
        assert reference['checks'][key] == expected, f'reference.checks.{key}'


def inspect(inputs, report, format='bam', sample_name=None):
    catalog_path = inputs / 'catalog.toml'
    catalog = tomllib.loads(catalog_path.read_text())
    variants = catalog['variants']
    targets = defaultdict(list)
    rows = []
    for vi, variant in enumerate(variants):
        for bucket, row in enumerate([variant, *variant['reference_junctions']]):
            key = junction(row)
            targets[key].append((vi, bucket))
            rows.append((variant['id'], 'defining' if bucket == 0 else row['role'], key))
    star = {}
    with (inputs / 'sample.SJ.out.tab').open() as stream:
        for columns in csv.reader(stream, delimiter='\t'):
            key = columns[0], int(columns[1]), int(columns[2]), {'1': '+', '2': '-', '0': '.'}[columns[3]]
            assert key not in star, 'duplicate STAR key'
            star[key] = int(columns[6]), int(columns[7])
    expected = {key: star.get(key, (0, 0)) for key in targets}
    observed = []
    for filename, defining in [('main.tsv', True), ('context.tsv', False)]:
        with (report / filename).open() as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            columns = ['variant_id', 'genome_assembly'] + ([] if defining else ['context_role']) + [
                'chromosome', 'intron_start', 'intron_end', 'strand',
                'unique_support', 'multimapping_support', 'total_support']
            assert reader.fieldnames == columns, 'TSV columns or order changed'
            actual = list(reader)
        relevant = [row for row in rows if (row[1] == 'defining') == defining]
        assert len(actual) == len(relevant)
        for actual_row, (identifier, role, key) in zip(actual, relevant):
            unique, multi = expected[key]
            assert (actual_row['variant_id'], junction(actual_row)) == (identifier, key)
            assert actual_row['genome_assembly'] == catalog['genome_assembly']
            if not defining:
                assert actual_row['context_role'] == role
            assert tuple(int(actual_row[k]) for k in ['unique_support', 'multimapping_support', 'total_support']) == (unique, multi, unique + multi)
            observed.append({'variant': identifier, 'role': role, 'junction': key,
                             'star_unique': unique, 'star_multi': multi})
    membership_counts = Counter()
    mapping_counts = Counter()
    strand_counts = Counter()
    cigar_counts = Counter()
    nh_counts = Counter()
    flags = Counter()
    inventory = Counter()
    row_statuses = Counter()
    opposite_memberships = 0
    pair_mates = defaultdict(set)
    mate_coordinates = defaultdict(list)
    inspected = []
    example_keys = set()
    coordinate_targets = {(c, start, end) for c, start, end, _ in targets}
    with pysam.AlignmentFile(str(inputs / 'reads.bam'), 'rb') as source:
        bam_header = source.header.to_dict()
        dictionary = [(r, n) for r, n in zip(source.references, source.lengths)]
        for read in source:
            if read.is_unmapped:
                continue
            cursor = read.reference_start
            introns = []
            for operation, size in read.cigartuples:
                if operation == 3:
                    introns.append((cursor + 1, cursor + size))
                if operation in [0, 2, 3, 7, 8]:
                    cursor += size
            if not any((read.reference_name, start, end) in coordinate_targets for start, end in introns):
                continue
            # STAR's emitted motif/coordinate tags are inspected independently.
            motifs = list(read.get_tag('jM'))
            assert len(motifs) == len(introns)
            assert all(m in [0, 1, 2, 3, 4, 5, 6, 20, 21, 22, 23, 24, 25, 26] for m in motifs)
            assert list(read.get_tag('jI')) == [x for pair in introns for x in pair]
            memberships = defaultdict(list)
            statuses = defaultdict(set)
            for index, (start, end) in enumerate(introns):
                motif = motifs[index] % 20
                direction = '.' if motif == 0 else '+' if motif % 2 else '-'
                for strand in ['+', '-']:
                    key = read.reference_name, start, end, strand
                    if key not in targets:
                        continue
                    if direction not in ['.', strand]:
                        opposite_memberships += len(targets[key])
                        continue
                    for vi, bucket in targets[key]:
                        memberships[vi].append(bucket)
                        statuses[vi].add('unverified' if direction == '.' else 'agreeing')
                        membership_counts[(vi, bucket)] += 1
                        strand_counts['unverified' if direction == '.' else 'agreeing'] += 1
            if not memberships:
                continue
            nh = read.get_tag('NH')
            assert isinstance(nh, int) and nh > 0
            mapping = 'unique' if nh == 1 else 'multimapping'
            mapping_counts[mapping] += 1
            nh_counts[str(nh)] += 1
            flags[str(read.flag)] += 1
            cigar_counts[read.cigarstring] += 1
            if read.is_paired:
                pair_mates[(read.query_name, read.get_tag('HI'))].add(1 if read.is_read1 else 2)
            identity = RecordIdentity(read.query_name, read.reference_name, read.reference_start + 1,
                        read.reference_end, read.cigarstring, read.mapping_quality, nh, read.flag,
                        read.get_tag('RG') if read.has_tag('RG') else 'unassigned', read.get_tag('HI'))
            mate_coordinates[identity].append((read.next_reference_name, read.next_reference_start + 1))
            for vi, buckets in memberships.items():
                row_statuses.update(statuses[vi] if len(statuses[vi]) == 1 else ['mixed'])
                inventory[(vi, mapping, ' '.join(map(str, sorted(buckets))), identity)] += 1
            example_key = (tuple((vi, tuple(buckets)) for vi, buckets in sorted(memberships.items())),
                           mapping, read.is_read1)
            if example_key not in example_keys:
                example_keys.add(example_key)
                inspected.append({'chromosome': read.reference_name, 'start_1based': read.reference_start + 1,
                                  'variant_memberships': {variants[vi]['id']: buckets for vi, buckets in memberships.items()},
                                  'flag': read.flag, 'cigar': read.cigarstring, 'NH': nh,
                                  'jM': motifs, 'jI': list(read.get_tag('jI'))})
    reference_dictionary = [(fields[0], int(fields[1]))
                            for line in (inputs / 'reference.fa.fai').read_text().splitlines()
                            for fields in [line.split('\t')]]
    assert dictionary == reference_dictionary, 'BAM does not retain the full reference dictionary'
    with pysam.AlignmentFile(str(inputs / 'reads.cram'), 'rc',
                             reference_filename=str(inputs / 'reference.fa')) as cram:
        assert list(zip(cram.references, cram.lengths)) == dictionary, 'CRAM sequence dictionary differs from BAM'
        cram_header = cram.header.to_dict()
    html = ReportInventory(report / 'report.html')
    assert html.records == inventory, 'embedded physical identities/classes/memberships differ from raw CIGAR inventory'
    assert html.statuses == row_statuses, 'per-membership strand labels disagree'
    for identities in html.pairs.values():
        assert len(identities) == 2, 'displayed pair does not contain two records'
        a, b = identities
        assert (a.name, a.chromosome, a.read_group, a.hit_index) == (
            b.name, b.chromosome, b.read_group, b.hit_index), 'pair name/contig/RG/HI mismatch'
        af, bf = a.flag, b.flag
        assert {af & 192, bf & 192} == {64, 128}, 'pair must contain read1 and read2'
        assert af & 1 and bf & 1 and not (af | bf) & (8 | 2048), 'ineligible pair flags'
        assert bool(af & 16) == bool(bf & 32) and bool(bf & 16) == bool(af & 32), 'mate orientation mismatch'
        assert bool(af & 256) == bool(bf & 256), 'primary/secondary pair mismatch'
        assert len(mate_coordinates[a]) == len(mate_coordinates[b]) == 1, 'ambiguous physical pair'
        assert mate_coordinates[a][0] == (b.chromosome, b.start) and mate_coordinates[b][0] == (a.chromosome, a.start), 'nonreciprocal mate coordinates'
    provenance = json.loads((report / 'provenance.json').read_text())
    inspect_metadata(provenance, inputs, catalog,
                     cram_header if format == 'cram' else bam_header, format, sample_name)
    assert provenance['alignment']['eligible_junction_memberships'] == sum(membership_counts.values())
    for row in observed:
        vi = next(i for i, v in enumerate(variants) if v['id'] == row['variant'])
        bucket = 0 if row['role'] == 'defining' else next(i + 1 for i, r in enumerate(variants[vi]['reference_junctions']) if r['role'] == row['role'])
        row['eligible_memberships'] = membership_counts[(vi, bucket)]
    return {'status': 'passed', 'variants': len(variants), 'catalog_rows': len(rows),
            'distinct_catalog_junctions': len(targets), 'star_junctions': len(star),
            'complete_bam_cram_reference_dictionary_equal': True,
            'header_contigs': len(dictionary), 'header_bases': sum(n for _, n in dictionary),
            'embedded_rows': sum(inventory.values()), 'eligible_physical_records': sum(mapping_counts.values()),
            'eligible_memberships': sum(membership_counts.values()),
            'validated_display_pair_groups': len(html.pairs),
            'selected_complete_mate_groups': sum(mates == {1, 2} for mates in pair_mates.values()),
            'selected_incomplete_mate_groups': sum(mates != {1, 2} for mates in pair_mates.values()),
            'excluded_opposite_memberships': opposite_memberships,
            'mapping_classes': dict(mapping_counts), 'NH': dict(nh_counts), 'flags': dict(flags),
            'strand_memberships': dict(strand_counts), 'html_row_strand_status': dict(html.statuses),
            'cigar_patterns': dict(cigar_counts), 'inspected_alignments_without_names': inspected,
            'junction_observations': observed, 'html_sha256': sha256(report / 'report.html')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path)
    parser.add_argument('report', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--format', choices=['bam', 'cram'], default='bam')
    parser.add_argument('--sample-name', help='the explicit sample name used for report generation')
    args = parser.parse_args()
    args.output.write_text(json.dumps(inspect(args.inputs, args.report, args.format, args.sample_name), indent=2) + '\n')


if __name__ == '__main__':
    main()
