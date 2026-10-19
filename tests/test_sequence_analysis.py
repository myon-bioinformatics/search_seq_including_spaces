from functools import lru_cache
from itertools import product
import unittest
from sequence_analysis import nucleotide_features, nucleotide_pairwise, fold_rna, marker_clues, translate_nucleotide


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

    def test_translate_nucleotide_functionality(self):
        r = translate_nucleotide('ATGGCC')
        self.assertEqual(r['schema'], 'nucleotide-translation/1')
        self.assertEqual(r['translations'][0]['protein_sequence'], 'MA')
        self.assertEqual(r['translations'][0]['nucleotide_span'], [1, 6])

        # Test RNA input
        r_rna = translate_nucleotide('AUGGCC')
        self.assertEqual(r_rna['translations'][0]['protein_sequence'], 'MA')

        # Test IUPAC degenerate codons
        # ACN -> T, CCN -> P, TTR -> L, TTY -> F, RTY -> X
        r_iupac = translate_nucleotide('ACNCCNTTRTTYRTY')
        self.assertEqual(r_iupac['translations'][0]['protein_sequence'], 'TPLFX')

        # Test reading frames and both strands
        r_all = translate_nucleotide('ATGGCC', frame='all', strand='both')
        self.assertEqual(len(r_all['translations']), 6) # 3 frames * 2 strands

        # Reverse strand of ATGGCC (complement GGCCAT -> GGCCAT rev = GGCCAT): GGC(G) CAT(H)
        r_rev = translate_nucleotide('ATGGCC', strand='-')
        self.assertEqual(r_rev['translations'][0]['protein_sequence'], 'GH')
        self.assertEqual(r_rev['translations'][0]['nucleotide_span'], [1, 6])

        # Test to_stop
        r_stop = translate_nucleotide('ATGTAAGCC', to_stop=True)
        self.assertEqual(r_stop['translations'][0]['protein_sequence'], 'M')
        self.assertFalse(r_stop['translations'][0]['has_internal_stop'])

        r_nostop = translate_nucleotide('ATGTAAGCC', to_stop=False)
        self.assertEqual(r_nostop['translations'][0]['protein_sequence'], 'M*A')
        self.assertTrue(r_nostop['translations'][0]['has_internal_stop'])

        r_term = translate_nucleotide('ATGGCCTAA', to_stop=False)
        self.assertEqual(r_term['translations'][0]['protein_sequence'], 'MA*')
        self.assertFalse(r_term['translations'][0]['has_internal_stop'])

        # Error cases
        with self.assertRaises(ValueError): translate_nucleotide('')
        with self.assertRaises(ValueError): translate_nucleotide('ATGGCC', table=2)
        with self.assertRaises(ValueError): translate_nucleotide('ATGGCC', frame=5)
        with self.assertRaises(ValueError): translate_nucleotide('ATGGCC', strand='invalid')
        with self.assertRaises(ValueError): translate_nucleotide('1234')


class TranslationStopSymbolTests(unittest.TestCase):
    def test_custom_symbol_rendering_counts_and_truncation(self):
        for symbol in ('#', 'M', 'X'):
            with self.subTest(symbol=symbol):
                full = translate_nucleotide('ATGTAAGCC', stop_symbol=symbol)['translations'][0]
                self.assertEqual(full['protein_sequence'], 'M' + symbol + 'A')
                self.assertEqual(full['stop_count'], 1)
                self.assertTrue(full['has_internal_stop'])
                self.assertEqual(full['length'], 3)
                self.assertEqual(full['nucleotide_span'], [1, 9])
                stopped = translate_nucleotide('ATGTAAGCC', stop_symbol=symbol, to_stop=True)['translations'][0]
                self.assertEqual(stopped['protein_sequence'], 'M')
                self.assertEqual(stopped['stop_count'], 0)
                self.assertFalse(stopped['has_internal_stop'])
                self.assertEqual(stopped['nucleotide_span'], [1, 3])

    def test_terminal_and_degenerate_stops_and_reverse_coordinates(self):
        terminal = translate_nucleotide('ATGTAR', stop_symbol='#')['translations'][0]
        self.assertEqual(terminal['protein_sequence'], 'M#')
        self.assertEqual(terminal['stop_count'], 1)
        self.assertFalse(terminal['has_internal_stop'])
        reverse = translate_nucleotide('GGCTTACAT', strand='-', stop_symbol='#', to_stop=True)['translations'][0]
        self.assertEqual(reverse['protein_sequence'], 'M')
        self.assertEqual(reverse['nucleotide_span'], [7, 9])

    def test_stop_symbol_must_be_one_character(self):
        for symbol in ('', '**', None, 1):
            with self.subTest(symbol=symbol), self.assertRaises(ValueError):
                translate_nucleotide('ATGTAA', stop_symbol=symbol)
