"""Check all ambiguity pairs against an independent set-based oracle."""
import unittest

from seqlib.anchor_matcher import find_matches
from seqlib.streaming_matcher import match_chunks

BASES = dict(zip('ACGTRYSWKMBDHVN',
                 ('A', 'C', 'G', 'T', 'AG', 'CT', 'CG', 'AT', 'GT', 'AC',
                  'CGT', 'AGT', 'ACT', 'ACG', 'ACGT')))


class IupacMaskTests(unittest.TestCase):
    def test_all_symbol_pairs_in_both_anchor_positions_and_chunk_boundary(self):
        for observed, bases in BASES.items():
            for anchor, allowed in BASES.items():
                compatible = bool(set(bases) & set(allowed))
                for text, left, right in ((observed + 'T', anchor, 'T'),
                                          ('A' + observed, 'A', anchor)):
                    with self.subTest(text=text, left=left, right=right):
                        report = find_matches(text, left, right, iupac=True)
                        self.assertEqual(len(report['matches']), int(compatible))
                        events = list(match_chunks(iter(text), left, right, iupac=True))
                        self.assertEqual([e['match'] for e in events[:-1]], report['matches'])
                        self.assertEqual(events[-1]['match_count'], int(compatible))
                        literal = find_matches(text, left, right)
                        self.assertEqual(bool(literal['matches']), text == left + right)
