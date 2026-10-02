import bz2
from dataclasses import asdict, FrozenInstanceError
import gzip
import io
import json
import lzma
from pathlib import Path
import tempfile
import unittest

from anchor_matcher import find_matches
from sequence_readers import read_sequences, SequenceReadError
from sequence_records import SCHEMA


class NonSeek(io.StringIO):
    def seekable(self):
        return False

    def read(self, *args):
        raise AssertionError('no sniff/read-all on non-seek stream')

    def seek(self, *args):
        raise AssertionError('no seek on non-seek stream')


class ReaderTests(unittest.TestCase):
    def test_fasta_metadata_and_line_endings(self):
        for newline in ('\n', '\r\n'):
            text = newline.join(['\ufeff;comment', '>one description here', 'ac N',
                                 'T-*', '', '>two', 'GG'])
            records = list(read_sequences(io.StringIO(text), format='fasta', source='fixture'))
            self.assertEqual([(r.id, r.description, r.sequence, r.quality) for r in records],
                             [('one', 'description here', 'acNT-*', None), ('two', '', 'GG', None)])
            self.assertEqual(asdict(records[0].provenance), dict(format='fasta', compression=None,
                             record_index=1, start_line=2, end_line=5))
            self.assertEqual(records[1].provenance.record_index, 2)
            self.assertEqual(records[0].source, 'fixture')
            json.dumps(dict(schema=SCHEMA, **asdict(records[0])))
            with self.assertRaises(FrozenInstanceError):
                records[0].id = 'changed'

    def test_fastq_wrapped_quality_and_line_endings(self):
        for newline in ('\n', '\r\n'):
            text = newline.join(['@a desc', 'ac', 'GT', '+a desc', '@+', '!!',
                                 '@b', 'N', '+', '~'])
            rows = list(read_sequences(io.StringIO(text), format='fastq'))
            self.assertEqual([(r.sequence, r.quality) for r in rows], [('acGT', '@+!!'), ('N', '~')])
            self.assertEqual((rows[0].provenance.start_line, rows[0].provenance.end_line), (1, 6))
            self.assertEqual(rows[1].provenance.record_index, 2)

    def test_all_plain_and_compressed_paths(self):
        codecs = [('', lambda b: b, None), ('.gz', gzip.compress, 'gzip'),
                  ('.bz2', bz2.compress, 'bz2'), ('.xz', lzma.compress, 'xz'),
                  ('.lzma', lambda b: lzma.compress(b, format=lzma.FORMAT_ALONE), 'lzma')]
        with tempfile.TemporaryDirectory() as directory:
            for ext, text in [('.fasta', '>id desc\r\nAC\r\nGT\r\n'),
                              ('.fastq', '@id desc\r\nACGT\r\n+\r\n!!!!\r\n')]:
                for suffix, compress, name in codecs:
                    with self.subTest(ext=ext, suffix=suffix):
                        path = Path(directory) / ('input' + ext + suffix)
                        path.write_bytes(compress(text.encode()))
                        row, = read_sequences(path)
                        self.assertEqual(row.sequence, 'ACGT')
                        self.assertEqual(row.provenance.compression, name)
                        self.assertEqual(row.source, str(path))
                        alias, = read_sequences(path, source='logical-source')
                        self.assertEqual(alias.source, 'logical-source')

    def test_detection_priority(self):
        stream = io.StringIO('>a\nAC\n')
        list(read_sequences(stream, format='fasta', source='wrong.fastq'))
        with self.assertRaisesRegex(SequenceReadError, '@ header'):
            list(read_sequences(io.StringIO('>a\nAC\n'), source='wrong.fastq'))
        for text, expected in [('>a\nAC\n', 'fasta'), ('@a\nAC\n+\n!!\n', 'fastq')]:
            row, = read_sequences(io.StringIO(text))
            self.assertEqual(row.provenance.format, expected)
        for text in ['ACGT', '@a\nAC\n+\n!', '>\nAC', 'x' * 9000 + '>a\nAC']:
            with self.assertRaisesRegex(SequenceReadError, 'conservatively detect'):
                list(read_sequences(io.StringIO(text)))
        with self.assertRaisesRegex(SequenceReadError, 'format must'):
            list(read_sequences(io.StringIO('>a\nAC'), format='genbank', source='a.fa'))

    def test_sniff_restores_current_position(self):
        stream = io.StringIO('ignored\n>a\nAC\n')
        stream.seek(len('ignored\n'))
        self.assertEqual(next(read_sequences(stream)).id, 'a')

    def test_nonseek_and_stream_ownership(self):
        stream = NonSeek('>a\nAC\n>b\nGT\n')
        with self.assertRaisesRegex(SequenceReadError, 'non-seek'):
            next(read_sequences(stream))
        self.assertEqual(stream.tell(), 0)
        rows = read_sequences(stream, format='fasta')
        self.assertEqual(next(rows).id, 'a')
        rows.close()
        self.assertFalse(stream.closed)
        named = NonSeek('@a\nAC\n+\n!!\n')
        self.assertEqual(next(read_sequences(named, source='a.fq')).quality, '!!')

    def test_malformed_fastq(self):
        cases = [('@a\nAC\n+\n!', 'quality-length mismatch'),
                 ('@a\nAC\n+\n!!!\n', 'quality-length mismatch'),
                 ('@a\nAC\n+', 'quality-length mismatch'),
                 ('@a\nAC\n', 'missing FASTQ +'),
                 ('@a\n+\n', 'empty FASTQ sequence'),
                 ('@a\nAC\n+b\n!!', 'does not match'),
                 ('a\nAC\n+\n!!', '@ header'),
                 ('@\nAC\n+\n!!', '@ header'),
                 ('@a\nA C\n+\n!!!', 'whitespace'),
                 ('@a\nAC\n+\n! ', 'invalid FASTQ quality')]
        for text, reason in cases:
            with self.subTest(text=text), self.assertRaisesRegex(SequenceReadError, reason) as caught:
                list(read_sequences(io.StringIO(text), format='fastq', source='bad'))
            self.assertEqual(caught.exception.source, 'bad')
            self.assertIsNotNone(caught.exception.line)

    def test_malformed_fasta_and_empty_file(self):
        for text in ['AC\n>a\nGG', '>\nAC', '>a\n>b\nAC', '>a\n']:
            with self.assertRaises(SequenceReadError):
                list(read_sequences(io.StringIO(text), format='fasta'))
        self.assertEqual(list(read_sequences(io.StringIO(''), format='fasta')), [])
        self.assertEqual(list(read_sequences(io.StringIO(''), format='fastq')), [])

    def test_records_are_yielded_before_later_error(self):
        rows = read_sequences(NonSeek('@ok\nAC\n+\n!!\nwrong\n'), format='fastq')
        self.assertEqual(next(rows).id, 'ok')
        with self.assertRaises(SequenceReadError):
            next(rows)

    def test_matcher_consumes_only_sequence(self):
        row, = read_sequences(io.StringIO('>a\nacN\nGT\n'), format='fasta')
        result = find_matches(row.sequence, 'AC', 'GT', min_gap=1)
        self.assertEqual(result['schema'], 'anchor-gap/1')
        self.assertEqual(result['matches'][0]['middle'], 'N')
        self.assertEqual(result['matches'][0]['source_span'], [1, 5])


if __name__ == '__main__':
    unittest.main()
