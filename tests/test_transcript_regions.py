import json
from pathlib import Path
import subprocess
import sys
import unittest

from seqlib.transcript_regions import transcript_regions
from sequence_service import call_tool

ROOT = Path(__file__).resolve().parents[1]
UTRS = ROOT / 'tests' / 'fixtures' / 'lin-utr.fa'


def search(*args):
    return json.loads(subprocess.check_output(
        [sys.executable, '-S', str(ROOT / 'ssis.py'), 'transcript-search', str(UTRS), *args]))


class TestTranscriptRegions(unittest.TestCase):
    def test_mirna_across_two_distinct_transcripts(self):
        report = search('--query', 'UCCCUGAGACCUCAAGUGUGA')
        self.assertEqual(report['records_searched'], 2)
        rows = {r['record']['id']: r for r in report['results']}
        for name, spans in [('NM_077515.5', [[832, 839], [880, 887], [1111, 1118]]),
                            ('NM_001025914.6', [[335, 342]])]:
            sites = rows[name]['seed_patterns']
            self.assertEqual(sites['counts'], {'8mer': len(spans)})
            self.assertEqual([s['target_span'] for s in sites['sites']], spans)

    def test_mrna_region_search_preserves_original_coordinates(self):
        report = search('--query-file', str(UTRS), '--query-record', 'NM_077515.5', '--query-span', '825', '846')
        rows = report['results']
        self.assertEqual(rows[0]['record']['id'], 'NM_077515.5')
        hit = rows[1]['similar']['hits'][0]
        self.assertEqual((hit['query_span'], hit['query_record_span'], hit['target_span']),
                         ([8, 15], [832, 839], [335, 342]))
        self.assertEqual(hit['matches'], 8)
        self.assertAlmostEqual(hit['query_coverage'], 8 / 22)

    def test_reverse_orientation_and_gaps(self):
        result = call_tool('transcript_regions', {'query': 'AUGCCUAGUA', 'target': 'TTTCTAGGCAT'})
        hit = result['reverse_complement']['hits'][0]
        self.assertEqual((hit['query_span'], hit['target_span'], hit['aligned_query']),
                         ([1, 8], [4, 11], 'CUAGGCAU'))
        for t, gaps, substitutions in [('AAAAAGCCCCC', 1, 0), ('AAAAAUCCCC', 0, 1)]:
            hit = transcript_regions('AAAAACCCCC', t)['similar']['hits'][0]
            self.assertEqual((hit['gap_columns'], hit['substitutions']), (gaps, substitutions))

    def test_ties_overlaps_and_ambiguous_bases(self):
        hits = transcript_regions('AAA', 'AAACAAA', top=1)['similar']
        self.assertEqual(hits['hits'][0]['target_span'], [1, 3])
        self.assertTrue(hits['more_hits'])
        self.assertEqual(transcript_regions('NNNN', 'NNNN')['similar']['hits'], [])
        self.assertIsNone(transcript_regions('ACGU', 'ACGU')['seed_patterns'])
        query = 'UCCCUGAGACCUCAAGUGUGA'
        sites = transcript_regions(query, 'CUCAGGGACUCAGGGA', top=1)['seed_patterns']
        self.assertEqual((sites['counts'], len(sites['sites']), sites['omitted']), ({'8mer': 2}, 1, 1))

    def test_site_classes_and_limits(self):
        q = 'UCCCUGAGACCUCAAGUGUGA'
        for target, kind in [('CUCAGGGA', '8mer'), ('CUCAGGG', '7mer-m8'),
                             ('UCAGGGA', '7mer-A1'), ('UCAGGG', '6mer')]:
            self.assertEqual(transcript_regions(q, target)['seed_patterns']['counts'], {kind: 1})
        for q, t in [('', 'A'), ('AC-', 'ACG'), ('A' * 3000, 'C' * 3000)]:
            with self.subTest(q=q[:8]), self.assertRaises(ValueError):
                transcript_regions(q, t)
