"""Score de novo alignments against pinned Rfam seed sequence rows."""
import argparse
from collections import Counter
import gzip
import hashlib
from io import StringIO
from itertools import combinations
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from seqlib.alignment import align_records
from seqlib.stockholm import read_stockholm

FAMILIES = {'RF00001', 'RF00003', 'RF00004', 'RF00005', 'RF00015', 'RF00020', 'RF00026'}
RELEASE_SHA256 = 'b06ac74d772609e863c5e1906eb491d699e9d325c5f22f4fc78126a808ad0bed'
FIXTURE = HERE / 'rfam-15.0.sto.gz'


def extract(release):
    if hashlib.sha256(release.read_bytes()).hexdigest() != RELEASE_SHA256:
        raise ValueError('expected the pinned Rfam 15.0 release archive')
    selected = []
    with gzip.open(release, 'rt', encoding='latin1') as source:
        block = []
        accession = None
        for line in source:
            if line.startswith('#=GF AC'):
                accession = line.split()[2]
            if line.startswith(('# STOCKHOLM', '#=GF AC', '#=GF ID')) or (line.strip() and not line.startswith('#')):
                block.append(line)
            if line.strip() == '//':
                if accession in FAMILIES:
                    selected.extend(block)
                block, accession = [], None
    # No timestamp or platform-specific header in the fixture.
    with FIXTURE.open('wb') as raw, gzip.GzipFile(fileobj=raw, filename='', mode='wb', mtime=0) as out:
        out.write(''.join(selected).encode('ascii'))


def residue_pairs(a, b):
    """Ungapped coordinates from paired reference or predicted columns, mismatches included."""
    i = j = 0
    paired = set()
    for x, y in zip(a, b, strict=True):
        i += x != '-'
        j += y != '-'
        if x != '-' and y != '-':
            paired.add((i, j))
    return paired


def metrics(counts):
    tp, fp, fn = (counts[key] for key in ('tp', 'fp', 'fn'))
    return {'tp': tp, 'fp': fp, 'fn': fn, 'precision': round(tp / (tp + fp), 6),
            'recall': round(tp / (tp + fn), 6), 'f1': round(2 * tp / (2 * tp + fp + fn), 6)}


def evaluate():
    families = {}
    total = {key: Counter() for key in ('aligned', 'left_justified')}
    data = gzip.decompress(FIXTURE.read_bytes())
    for block in data.decode('ascii').split('//\n'):
        if not block.strip():
            continue
        accession = next(line.split()[2] for line in block.splitlines() if line.startswith('#=GF AC'))
        label = next(line.split()[2] for line in block.splitlines() if line.startswith('#=GF ID'))
        reference, _ = read_stockholm(StringIO(block + '//\n'))
        unique = {}
        for name in sorted(reference):
            unique.setdefault(reference[name].replace('-', '').upper(), name)
        ids = sorted(unique.values(), key=lambda name: hashlib.sha256(name.encode('ascii')).digest())
        groups = [ids[i:i + 16] for i in range(0, len(ids), 16)]
        if len(groups[-1]) == 1:
            groups[-2].extend(groups.pop())
        counts = {key: Counter() for key in total}
        pair_count = 0
        for group in groups:
            raw = {name: reference[name].replace('-', '').upper() for name in group}
            predicted = align_records(raw.items())
            longest = max(map(len, raw.values()))
            baseline = {key: seq.ljust(longest, '-') for key, seq in raw.items()}
            for a, b in combinations(group, 2):
                expected = residue_pairs(reference[a], reference[b])
                pair_count += 1
                for key, rows in [('aligned', predicted), ('left_justified', baseline)]:
                    actual = residue_pairs(rows[a], rows[b])
                    tp = len(expected & actual)
                    counts[key].update(tp=tp, fp=len(actual) - tp, fn=len(expected) - tp)
        for key in total:
            total[key].update(counts[key])
        families[accession] = {'name': label, 'source_rows': len(reference), 'unique_sequences': len(ids),
                               'groups': len(groups), 'sequence_pairs': pair_count,
                               **{k: metrics(v) for k, v in counts.items()}}
    if families.keys() != FAMILIES:
        raise ValueError('reference family set differs from pinned selection')
    return {'release': 'Rfam 15.0', 'reference_sha256': hashlib.sha256(data).hexdigest(),
            'families': families, 'micro': {k: metrics(v) for k, v in total.items()},
            **{key: sum(f[key] for f in families.values())
               for key in ('source_rows', 'unique_sequences', 'groups', 'sequence_pairs')}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extract', type=Path, metavar='Rfam.seed.gz', help='rebuild the fixture from release 15.0')
    parser.add_argument('--check', action='store_true', help='compare measured results to the recorded JSON')
    args = parser.parse_args()
    if args.extract:
        extract(args.extract)
    result = evaluate()
    print(json.dumps(result, sort_keys=True, indent=2))
    if args.check and result != json.loads((HERE / 'results.json').read_text()):
        sys.exit('measured alignment results differ from the recorded run')
