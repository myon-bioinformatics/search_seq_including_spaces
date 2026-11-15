"""Preprocessing with transform receipts and record export (seqlib/sequence_preprocess.py)."""
import io, unittest
from seqlib.sequence_preprocess import preprocess, export_records
from seqlib.sequence_readers import read_sequences


class PreprocessTests(unittest.TestCase):
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
