import io
import unittest
from protein_ionization import ionization
from sequence_preprocess import preprocess, export_records
from sequence_readers import read_sequences
from sequence_router import route_item


class PreprocessRouterTests(unittest.TestCase):
    def test_normalization_name_collision_dedup_and_receipts(self):
        records=read_sequences(io.StringIO('>x?\nacgu\n>x!\nACGT\n>x_\nNNNN\n'),format='fasta')
        out,report=preprocess(records,normalize_ut=True,deduplicate='sequence',max_ambiguity=0.5)
        self.assertEqual([r.sequence for r in out],['ACGT'])
        self.assertEqual([r['reason'] for r in report['receipts']],[None,'duplicate','ambiguity'])
        self.assertEqual(report['receipts'][0]['input_id'],'x?')
        self.assertEqual(export_records(out),'>x_\nACGT\n')
        records=read_sequences(io.StringIO('>x?\nAC\n>x!\nAT\n>x_2\nAG\n'),format='fasta')
        out,_=preprocess(records)
        self.assertEqual([r.id for r in out],['x_','x__2','x_2'])

    def test_alignment_alphabet_export_validation(self):
        records=list(read_sequences(io.StringIO('>x\nAC\n>y\nACT\n'),format='fasta'))
        with self.assertRaises(ValueError): preprocess(records,aligned=True)
        out,_=preprocess(records)
        with self.assertRaises(ValueError): export_records(out,format='mega')
        with self.assertRaises(ValueError): export_records(out,format='tsv')
        with self.assertRaises(ValueError): preprocess(records,alphabet='protein',normalize_ut=True)

    def test_ionization_reference_and_monotonicity(self):
        # Published Bio.SeqUtils.IsoelectricPoint example values, rounded there.
        self.assertAlmostEqual(ionization('INGAR')['estimated_pi'],9.75,places=2)
        self.assertAlmostEqual(ionization('INGAR')['net_charge'],0.76,places=2)
        pi=ionization('PETER')['estimated_pi']
        self.assertAlmostEqual(ionization('PETER',ph=pi)['net_charge'],0,places=10)
        self.assertGreater(ionization('DKR',ph=3)['net_charge'],ionization('DKR',ph=10)['net_charge'])
        with self.assertRaises(ValueError): ionization('X')

    def test_routing_failure_uncertainty_no_side_effect(self):
        item={'kind':'analysis','features':{'score':4}}
        self.assertEqual(route_item(item)['route'],'local')
        def fail(request): raise TimeoutError('private transport details')
        r=route_item(item,classifier=fail,threshold=0.8)
        self.assertEqual(r['reason'],'classifier-failed')
        self.assertNotIn('private',str(r))
        good=lambda request:{'decision':'no','confidence':0.9}
        self.assertEqual(route_item(item,classifier=good,threshold=0.8)['route'],'exclude')
        self.assertEqual(route_item(item,classifier=good,threshold=0.95)['route'],'local')
        self.assertEqual(route_item(item,classifier=fail,threshold=0.8,important=True)['route'],'pro')
        with self.assertRaises(ValueError): route_item(item,classifier=good)
        r=route_item({'kind':'sequence','sequence':'ACGT','alphabet':'dna'},question='sequence_type',classifier=fail,threshold=0.8)
        self.assertEqual(r['route'],'deterministic')
