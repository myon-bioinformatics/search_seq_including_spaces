from io import StringIO
from itertools import product
import json
from pathlib import Path
import subprocess
import sys
import unittest

from seqlib.alignment import align_records
from seqlib.stockholm import annotated_matches, read_stockholm

ROOT = Path(__file__).resolve().parents[1]


def optimal_score(a, b):
    # Exhaustive alignment paths, independent of the aligner's traceback.
    if not a or not b:
        return -2 * (len(a) + len(b))
    return max((2 if a[0] == b[0] else -1) + optimal_score(a[1:], b[1:]),
               -2 + optimal_score(a, b[1:]), -2 + optimal_score(a[1:], b))


class TestDeNovoAlignment(unittest.TestCase):
    def test_pairwise_optimum(self):
        words = [''.join(word) for n in (1, 2, 3) for word in product('AC', repeat=n)]
        for a, b in product(words, repeat=2):
            rows = align_records([('a', a), ('b', b)])
            scores = [-2 if '-' in (x, y) else 2 if x == y else -1 for x, y in zip(*rows.values(), strict=True)]
            self.assertEqual(sum(scores), optimal_score(a, b))
            self.assertEqual([s.replace('-', '') for s in rows.values()], [a, b])

    def test_three_rows_ties_and_invalid_input(self):
        records = [('a', 'acgu'), ('b', 'aacgu'), ('c', 'acgugg')]
        expected = {'a': '-ACGU--', 'b': 'AACGU--', 'c': '-ACGUGG'}
        self.assertEqual(align_records(records), expected)
        self.assertEqual(align_records(reversed(records)), expected)
        self.assertEqual(align_records([('a', 'AA'), ('b', 'A')])['b'], '-A')
        self.assertEqual(align_records([('a', 'UNT'), ('b', 'TNU')]), {'a': 'UNT', 'b': 'TNU'})
        for records in ([], [('a', 'AC')], [('a', 'AC'), ('a', 'A')], [('a', 'A-C'), ('b', 'AC')],
                        [('a', ''), ('b', 'A')], [('a', 'A?'), ('b', 'A')],
                        [('a', 'A' * 3000), ('b', 'A' * 3000)]):
            with self.assertRaises(ValueError):
                align_records(records)

    def test_cli_round_trip(self):
        cmd = [sys.executable, '-S', str(ROOT / 'ssis.py'), 'align', '-']
        seqs = '>a\nACGU\n>b\nAGU\n'
        sto = subprocess.check_output(cmd, input=seqs, text=True)
        rows, tracks = read_stockholm(StringIO(sto))
        self.assertEqual(rows, {'a': 'ACGU', 'b': 'A-GU'})
        self.assertEqual(tracks['coverage'], '*.**')
        hit, = annotated_matches(StringIO(sto), track='coverage', left='AC', right='GU')['results']
        self.assertEqual(hit['columns'], [1, 4])
        self.assertEqual(subprocess.check_output(cmd + ['--format', 'fasta'], input=seqs, text=True),
                         '>a\nACGU\n>b\nA-GU\n')

    def test_rfam_reference_benchmark(self):
        report = json.loads(subprocess.check_output(
            [sys.executable, '-S', str(ROOT / 'experiments/alignment/evaluate.py'), '--check']))
        self.assertEqual(report['sequence_pairs'], 17800)
