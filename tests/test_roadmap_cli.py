import json
from pathlib import Path
import tempfile
import unittest
from test_cli import run


class RoadmapCliTests(unittest.TestCase):
    def test_match_jsonl_pack_stdin_and_errors(self):
        code,out,err=run('match','-','--input-format','fasta','--left','AC','--right','GT','--min-gap',1,'--jsonl',stdin=b'>x\nACNGT\n>y\nTTTT\n')
        self.assertEqual((code,err),(0,''))
        docs=[json.loads(line) for line in out.splitlines()]
        self.assertEqual(docs[0]['result']['matches'][0]['middle'],'N')
        self.assertEqual(docs[1]['result']['matches'],[])
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'patterns.json'
            p.write_text(json.dumps({'schema':'sequence-patterns/1','patterns':[{'id':'p','left':'A','right':'T'}]}))
            code,out,err=run('match','-','--input-format','jsonl','--patterns',p,'--sequence-path','/payload/seq',stdin=b'{"payload":{"seq":"AT"}}\n')
            self.assertEqual((code,err),(0,''))
            self.assertEqual(json.loads(out)['pattern_id'],'p')
            code,out,err=run('match','-','--input-format','fasta','--patterns',p,'--left','A',stdin=b'')
            self.assertEqual((code,out),(2,b''))
        for opts in (('--min-gap','-1'),('--max-matches','0')):
            code,out,err=run('match','-','--input-format','fasta','--left','A','--right','T',*opts,stdin=b'')
            self.assertEqual((code,out),(2,b''))
            self.assertNotIn('Traceback',err)

    def test_chunk_jsonl_and_summary_cli(self):
        code,out,err=run('match','-','--input-format','fasta','--left','A','--right','T','--chunk-size','4','--jsonl',stdin=b'>x\nCCCCAT\n')
        self.assertEqual((code,err),(0,''))
        self.assertEqual(json.loads(out.splitlines()[-1])['match_count'],1)
        with tempfile.TemporaryDirectory() as d:
            a,b=Path(d)/'a.fa',Path(d)/'b.fa'
            a.write_text('>a\nAACC');b.write_text('>b\nAATT')
            code,out,err=run('compare',a,b,'--alphabet','dna','--window','2','--step','2','--summary','--top','1')
            self.assertEqual((code,err),(0,''))
            self.assertEqual(json.loads(out)['omitted_items'],1)

    def test_analysis_candidates_and_preprocess_receipt(self):
        for args,data,schema in (
            (('candidates','--mode','pairwise'),b'>x\nNNN\n','sequence-candidates/1'),
            (('nucleotide-features',),b'>x\nACGTN\n','nucleotide-features/1'),
            (('rna-fold',),b'>x\nGGGAAACCC\n','rna-nussinov/1'),
            (('protein-ionization',),b'>x\nINGAR\n','protein-ionization/1'),
            (('preset','--kind','poly-a'),b'>x\nAAAAAA\n','sequence-marker-clues/1')):
            code,out,err=run(*args,'-','--input-format','fasta',stdin=data)
            self.assertEqual((code,err),(0,''))
            self.assertEqual(json.loads(out)['schema'],schema)
        with tempfile.TemporaryDirectory() as d:
            receipt=Path(d)/'receipt.json'
            args=('preprocess','-','--input-format','fasta','--export-format','mega','--receipt',receipt)
            code,out,err=run(*args,stdin=b'>a\nACGT\n>b\nAC-T\n')
            self.assertEqual((code,err),(0,''))
            self.assertTrue(out.startswith(b'#MEGA'))
            self.assertEqual(json.loads(receipt.read_text())['output_records'],2)
            before=receipt.read_bytes()
            code,out,err=run(*args,stdin=b'>a\nAAAA\n')
            self.assertEqual((code,out),(2,b''))
            self.assertEqual(receipt.read_bytes(),before)
