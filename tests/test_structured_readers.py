import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from sequence_readers import read_sequences, SequenceReadError

GB='''LOCUS       DEMO 8 bp DNA
DEFINITION  demo sequence
            continuation.
ACCESSION   DEMO
VERSION     DEMO.1
FEATURES             Location/Qualifiers
ORIGIN
        1 acgt nnac
//
'''


class StructuredReadersTests(unittest.TestCase):
    def read(self,text,format='json',**options):
        return list(read_sequences(io.StringIO(text),format=format,**options))

    def test_json_pointer_and_bounding_span(self):
        text='[\n{"id":"α", "payload":{"a/b":["a c\\nT"]}}\n]'
        record=self.read(text,sequence_path='/payload/a~1b/0',id_path='id')[0]
        self.assertEqual((record.id,record.sequence),('α','a c\nT'))
        self.assertEqual((record.provenance.start_line,record.provenance.end_line),(1,3))

    def test_jsonl_laziness_line_spans_ownership(self):
        handle=io.StringIO('\ufeff\n{"sequence":"AC"}\nBAD\n')
        records=read_sequences(handle,format='jsonl')
        record=next(records)
        self.assertEqual(record.provenance.start_line,2)
        with self.assertRaises(SequenceReadError): next(records)
        self.assertFalse(handle.closed)

    def test_bad_json_shapes_types_duplicates_and_limits(self):
        for text in ('null','[null]','{"sequence":4}','{"sequence":""}',
                     '{"sequence":"A","sequence":"T"}','{"sequence":"A","v":NaN}'):
            with self.assertRaises(SequenceReadError): self.read(text)
        for fmt in ('json','jsonl'):
            with patch('structured_readers.MAX_JSON_CHARS',20):
                with self.assertRaises(SequenceReadError): self.read('{"sequence":"'+'A'*30+'"}',format=fmt)
        with self.assertRaises(SequenceReadError): self.read('{"sequence":"AC"}',id_path='missing')
        with self.assertRaises(SequenceReadError): self.read('>x\nAC\n',format='fasta',sequence_path='sequence')
        with self.assertRaises(SequenceReadError): self.read('{"sequence":"AC"}',sequence_column='sequence')

    def test_genbank_metadata_multirecord_and_errors(self):
        records=self.read(GB+GB,format='genbank')
        self.assertEqual(len(records),2)
        self.assertEqual(records[0].sequence,'acgtnnac')
        self.assertEqual(records[0].metadata,{'locus':'DEMO','accession':'DEMO','version':'DEMO.1'})
        self.assertEqual(records[0].description,'demo sequence continuation.')
        self.assertEqual(records[1].provenance.record_index,2)
        for text in (GB.replace('//',''),GB.replace('ORIGIN','CONTIG'),GB.replace('1 acgt','2 acgt'),GB.replace('acgt','acgt!')):
            with self.assertRaises(SequenceReadError): self.read(text,format='genbank')

    def test_compression_and_extension_detection(self):
        with tempfile.TemporaryDirectory() as directory:
            for suffix,text in (('jsonl','{"sequence":"AC"}\n'),('gb',GB)):
                path=Path(directory)/('sequence.'+suffix+'.gz')
                with gzip.open(path,'wt',encoding='utf-8') as handle: handle.write(text)
                record=list(read_sequences(path))[0]
                self.assertEqual(record.provenance.compression,'gzip')
