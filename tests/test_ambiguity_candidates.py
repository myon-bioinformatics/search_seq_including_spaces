from itertools import combinations, product
import re
import unittest

from ambiguity_candidates import generate_candidates, select_representatives, exact_hit_matrix
from anchor_matcher import _BASES


class CandidateTests(unittest.TestCase):
    def test_pairwise_obligations_independently_enumerated(self):
        for text in ('A', 'N', 'NN', 'RYSWKMBDHVN', 'NNNNNNNNNN'):
            doc = generate_candidates(text,mode='pairwise',seed=17)
            candidates = doc['candidates']
            self.assertEqual(len(set(candidates)),len(candidates))
            for candidate in candidates:
                self.assertTrue(all(base in _BASES[code] for base,code in zip(candidate,text)))
            positions = [i for i,c in enumerate(text) if len(_BASES[c])>1]
            for i,j in combinations(positions,2):
                self.assertEqual({(row[i],row[j]) for row in candidates},set(product(_BASES[text[i]],_BASES[text[j]])))
            self.assertEqual(doc,generate_candidates(text,mode='pairwise',seed=17))
        self.assertLess(len(generate_candidates('N'*10,mode='pairwise')['candidates']),4**10)

    def test_sample_large_space_unique_and_reproducible(self):
        doc = generate_candidates('N'*100,mode='sample',seed=3,limit=17)
        self.assertEqual(doc['possible_candidates'],4**100)
        self.assertEqual(len(set(doc['candidates'])),17)
        self.assertEqual(doc,generate_candidates('N'*100,mode='sample',seed=3,limit=17))

    def test_modes_and_exhaustive_boundary(self):
        self.assertEqual(generate_candidates('AR',mode='exhaustive',limit=2)['candidates'],['AA','AG'])
        self.assertEqual(generate_candidates('AR')['sequence'],'AR')
        pattern=generate_candidates('AR',mode='regex')['pattern']
        self.assertIsNotNone(re.fullmatch(pattern,'AG'))
        self.assertIsNone(re.fullmatch(pattern,'AC'))
        for options in ({'mode':'exhaustive','limit':1},{'mode':'pairwise','limit':1},{'seed':-1},{'mode':'bad'},{'limit':True}):
            with self.assertRaises(ValueError): generate_candidates('NN',**options)
        with self.assertRaises(ValueError): generate_candidates('U')
        with self.assertRaises(ValueError): generate_candidates('N'*65,mode='pairwise')

    def test_exact_matrix_and_approximation_loss_are_auditable(self):
        candidates=['AA','CC','GG','TT']
        targets=['AACC','CCGG','TT']
        self.assertEqual(exact_hit_matrix(candidates,targets),[1,3,2,4])
        report=select_representatives(candidates,targets,rank=1,seed=0)
        selected=report['selected_indices']
        all_hits=set(j for i,row in enumerate(candidates) for j,t in enumerate(targets) if row in t)
        retained=set(j for i in selected for j,t in enumerate(targets) if candidates[i] in t)
        self.assertEqual(report['lost_target_hits'],len(all_hits-retained))
        self.assertEqual(report['omitted_candidates'],len(candidates)-len(selected))
        self.assertEqual(report,select_representatives(candidates,targets,rank=1,seed=0))
        self.assertFalse(report['pairwise_coverage_preserved'])
