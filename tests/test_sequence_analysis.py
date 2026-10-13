from functools import lru_cache
from itertools import product
import unittest
from sequence_analysis import nucleotide_features, nucleotide_pairwise, fold_rna, marker_clues


class AnalysisTests(unittest.TestCase):
    def test_features_denominators_and_ambiguity(self):
        r=nucleotide_features('AACGN',k=2)
        self.assertEqual((r['canonical_count'],r['ambiguous_count'],r['gc_fraction']),(4,1,0.5))
        self.assertEqual(r['kmers'],{'AA':1,'AC':1,'CG':1})
        self.assertEqual(r['excluded_kmers'],1)
        self.assertEqual(r['homopolymers'],[{'symbol':'A','start':1,'end':2,'length':2}])
        self.assertIsNone(nucleotide_features('NN')['gc_fraction'])
        self.assertEqual(nucleotide_features('AU',alphabet='rna')['at_or_au_fraction'],1)
        with self.assertRaises(ValueError): nucleotide_features('AC',k=7)

    def test_transitions_transversions_gaps_and_ambiguity(self):
        r=nucleotide_pairwise('ACGTN-','GTCTAA')
        self.assertEqual(r['substitutions']['transitions'],2)
        self.assertEqual(r['substitutions']['transversions'],1)

    def test_marker_boundaries_and_candidates(self):
        self.assertEqual(marker_clues('CAAAAAAT',kind='poly-a')['results'][0]['start'],2)
        r=marker_clues('CATATATG',kind='str',motif='AT',min_length=3)
        self.assertEqual(r['results'][0]['copies'],3)
        r=marker_clues('AN',kind='snp',position=2,allele='A')['results'][0]
        self.assertFalse(r['exact_match']);self.assertTrue(r['compatible'])
        r=marker_clues('ACGGAC',kind='tsd',repeat_length=2,min_gap=2,max_gap=2)
        self.assertEqual(r['results'][0]['right_span'],[5,6])
        with self.assertRaises(ValueError): marker_clues('AAAA',kind='tsd',repeat_length=1,max_results=1)

    def test_nussinov_exhaustive_short_oracle_and_structure(self):
        # Independent recurrence enumerates whether LEFT endpoint is paired.
        @lru_cache(None)
        def oracle(sequence):
            if not sequence: return 0
            best=oracle(sequence[1:])
            for j in range(2,len(sequence)):
                if sequence[0]+sequence[j] in ('AU','UA','GC','CG','GU','UG'):
                    best=max(best,1+oracle(sequence[1:j])+oracle(sequence[j+1:]))
            return best
        for letters in product('ACGU',repeat=5):
            sequence=''.join(letters)
            r=fold_rna(sequence,min_loop=1)
            self.assertEqual(r['pair_count'],oracle(sequence))
            self.assertEqual(r['dot_bracket'].count('('),r['pair_count'])
            positions=[x for pair in r['pairs'] for x in pair]
            self.assertEqual(len(positions),len(set(positions)))
            for a,b in r['pairs']:
                self.assertGreater(b-a,1)
                for c,d in r['pairs']: self.assertFalse(a<c<b<d)
        r=fold_rna('GGGAAACCC')
        self.assertEqual((r['pair_count'],r['dot_bracket']),(3,'(((...)))'))
        self.assertEqual(r['loops'][0]['kind'],'hairpin')
        for text,opts in (('AN',{}),('AT',{}),('A'*301,{}),('AU',{'min_loop':-1})):
            with self.assertRaises(ValueError): fold_rna(text,**opts)

    def test_explicit_codon_frame_counts_and_partial_tail(self):
        result=nucleotide_features('AACGTNNNA',codon_frame=1)
        self.assertEqual(result['codons'],{'frame':1,'counts':{'ACG':1},'excluded_triplets':1,'trailing_bases':2})
