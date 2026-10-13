import csv
from pathlib import Path
import unittest

from comparison_s2 import sliding_identity
from result_summary import summarize_hits, summarize_csv, summarize_windows

ROOT = Path(__file__).resolve().parents[1]


class SummaryTests(unittest.TestCase):
    def test_real_csv_conserves_all_5966_hits_and_original(self):
        path = next((ROOT / 'output_folder/example').glob('semi*'))
        before = path.read_bytes()
        doc = summarize_csv(path, top=3)
        self.assertEqual(doc['total_hits'], 5966)
        self.assertEqual(doc['total_regions'], 1246)
        self.assertEqual([r['hit_count'] for r in doc['regions']], [49, 38, 35])
        self.assertEqual(sum(r['hit_count'] for r in doc['regions']) + doc['omitted_items'], 5966)
        self.assertEqual(doc['omitted_regions'], 1243)
        self.assertEqual(before, path.read_bytes())
        self.assertEqual(summarize_csv(path, gap_merge=0)['total_regions'], 5966)

    def test_boundaries_sort_duplicates_and_original_rows(self):
        rows = [dict(start_num=s, end_num=s+2, back_elements=b) for s,b in [(21,'B'),(1,'A'),(11,'A'),(11,'B'),(32,'C')]]
        doc = summarize_hits(rows, gap_merge=10, top=1)
        self.assertEqual(doc['regions'], [{'start':1,'end':23,'hit_count':4,'source_rows':[1,4]}])
        self.assertEqual(doc['omitted_items'], 1)
        self.assertEqual(doc['omitted_back_element_hits'], 3)
        self.assertEqual(summarize_hits(rows,gap_merge=0)['total_regions'],5)

    def test_density_coordinate_semantics_empty_and_bad_input(self):
        rows = [dict(start_num=1,end_num=1,back_elements='A')]
        self.assertEqual(summarize_hits(rows,coordinate_convention='1-based-inclusive')['density_per_mb'],1000000)
        self.assertEqual(summarize_hits(rows,coordinate_convention='1-based-inclusive',sequence_length=100)['density_per_mb'],10000)
        self.assertIsNone(summarize_hits([])['density_per_mb'])
        for options in ({'gap_merge':-1},{'top':False},{'sequence_length':0}):
            with self.assertRaises(ValueError): summarize_hits([],**options)
        with self.assertRaises(ValueError): summarize_hits(rows)
        with self.assertRaises(ValueError): summarize_hits([dict(start_num=0.5,end_num=2,back_elements='A')])

    def test_window_summary_preserves_unknown_and_counts(self):
        doc = sliding_identity('AACCGG--','AATTGG--',alphabet='dna',window=2,step=2)
        report = summarize_windows(doc,threshold=100,top=2)
        self.assertEqual(report['classification_counts'], {'conserved':2,'divergent':1,'unknown':1})
        self.assertEqual(report['total_regions'],4)
        self.assertEqual(sum(r['window_count'] for r in report['regions'])+report['omitted_items'],4)
        self.assertEqual(report['regions'][0]['source_rows'],[1,1])

    def test_window_gaps_not_bridged_overlap_disclosed_and_coverage(self):
        doc = sliding_identity('AAAAAAAA','AAAAAAAA',alphabet='dna',window=2,step=3)
        self.assertEqual(summarize_windows(doc)['total_regions'],3)
        doc = sliding_identity('A-','AA',alphabet='dna',window=2)
        self.assertEqual(summarize_windows(doc,min_coverage=60)['regions'][0]['classification'],'unknown')
        for threshold in (float('nan'),-1,101,True):
            with self.assertRaises(ValueError): summarize_windows(doc,threshold=threshold)
