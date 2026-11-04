"""Prepare documented public-read subsets and separately labeled sequence controls."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tomllib


def key(row):
    return row['chromosome'], row['intron_start'], row['intron_end'], row['strand']


def subset(data, output, pairs):
    expected = ['d57df747bc142e9850074d512ab9d6db', '3331c6a9e0183ff9d398a3292dd45f66']
    hashes = {}
    for mate, digest in enumerate(expected, 1):
        path = data / f'R{mate}.fastq.gz'
        with path.open('rb') as source:
            actual = hashlib.file_digest(source, 'md5').hexdigest()
        assert actual == digest, f'mate {mate}: source MD5 mismatch'
        hashes[path.name] = actual
    with gzip.open(data / 'R1.fastq.gz', 'rt') as first, gzip.open(data / 'R2.fastq.gz', 'rt') as second:
        with (output / 'R1.fastq').open('w') as out1, (output / 'R2.fastq').open('w') as out2:
            for number in range(1, pairs + 1):
                records = [[source.readline().rstrip('\n') for _ in range(4)] for source in [first, second]]
                names = [record[0].split()[0].removesuffix('/1').removesuffix('/2') for record in records]
                assert names[0] == names[1] and names[0].startswith('@'), f'unpaired spot {number}'
                for record, destination in zip(records, [out1, out2]):
                    assert record[2].startswith('+') and len(record[1]) == len(record[3]) > 0
                    # Alias names keep raw source read identifiers out of local report diagnostics.
                    destination.write(f'@public-{number:08}\n{record[1]}\n+\n{record[3]}\n')
            if pairs == 22935521:
                assert not first.read(1) and not second.read(1), 'expected complete run has extra records'
    result = {'pairs': pairs, 'sampling': 'complete paired run' if pairs == 22935521 else
              'first complete paired records; not randomized', 'source_md5': hashes}
    for mate in [1, 2]:
        with (output / f'R{mate}.fastq').open('rb') as source:
            result[f'mate_{mate}_sha256'] = hashlib.file_digest(source, 'sha256').hexdigest()
    return result


def controls(reference, catalog, output):
    import pysam
    variants = tomllib.loads(catalog.read_text())['variants']
    junctions = sorted({key(row) for variant in variants for row in [variant, *variant['reference_junctions']]})
    complement = str.maketrans('ACGTN', 'TGCAN')
    with pysam.FastaFile(str(reference)) as genome, (output / 'R1.fastq').open('w') as first, \
            (output / 'R2.fastq').open('w') as second, (output / 'junctions.tab').open('w') as annotations:
        for number, (chromosome, start, end, strand) in enumerate(junctions):
            sequence = (genome.fetch(chromosome, start - 51, start - 1) + genome.fetch(chromosome, end, end + 50)).upper()
            assert len(sequence) == 100 and set(sequence) <= set('ACGT'), (chromosome, start, end)
            annotations.write(f'{chromosome}\t{start}\t{end}\t{1 if strand == "+" else 2}\n')
            for replicate in range(4):
                name = f'control-{number:02}-{replicate}'
                first.write(f'@{name}\n{sequence}\n+\n{"I" * 100}\n')
                second.write(f'@{name}\n{sequence.translate(complement)[::-1]}\n+\n{"I" * 100}\n')
    return {'unique_junctions': len(junctions), 'pairs': len(junctions) * 4,
            'description': 'Artificial 100-base fragments across exact catalog junctions, both mates overlap the junction. Not biological observations.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['subset', 'controls'])
    parser.add_argument('data', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--pairs', type=int, default=1000000)
    parser.add_argument('--catalog', type=Path, default=Path('definitions/grch38.toml'))
    args = parser.parse_args()
    if args.pairs < 1:
        parser.error('pairs must be positive')
    args.output.mkdir(exist_ok=False)
    result = subset(args.data, args.output, args.pairs) if args.action == 'subset' else controls(
        args.data / 'reference.fa', args.catalog, args.output)
    (args.output / 'preparation.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
