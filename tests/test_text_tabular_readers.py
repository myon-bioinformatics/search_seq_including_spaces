import bz2
import csv
import gzip
import io
import lzma
from pathlib import Path
import tempfile
import unittest

from sequence_readers import read_sequences, SequenceReadError
from test_sequence_readers import NonSeek


class TextTabularTests(unittest.TestCase):
    def rows(self, text, format='csv', **kwargs):
        return list(read_sequences(io.StringIO(text), format=format, **kwargs))

    def test_plain_text_one_record_and_span(self):
        row, = self.rows('\ufeff\n ac N\n\nT-*\n \n', 'text')
        self.assertEqual((row.id, row.description, row.sequence, row.quality),
                         ('record-1', '', 'acNT-*', None))
        self.assertEqual((row.provenance.start_line, row.provenance.end_line), (2, 5))
        self.assertEqual(self.rows(' \n\t', 'text'), [])
        self.assertEqual(self.rows('', 'text'), [])
        self.assertEqual(self.rows('éΩ\n', 'text')[0].sequence, 'éΩ')

    def test_explicit_extension_and_no_generic_sniff(self):
        for text in ('ACGT\n', 'id,sequence\na,AC\n', 'name\tsequence\na\tAC\n'):
            with self.assertRaisesRegex(SequenceReadError, 'conservatively detect'):
                list(read_sequences(io.StringIO(text), sequence_column='sequence'))
        row, = read_sequences(NonSeek('AC GT\n'), source='input.txt.gz')
        self.assertEqual(row.sequence, 'ACGT')
        self.assertIsNone(row.provenance.compression)
        row, = read_sequences(NonSeek('sequence\nAC\n'), source='a.csv', sequence_column='sequence')
        self.assertEqual(row.id, 'record-1')
        row, = read_sequences(io.StringIO('>id\nAC\n'), format='text', source='a.fa')
        self.assertEqual(row.sequence, '>idAC')
        with self.assertRaisesRegex(SequenceReadError, 'column options'):
            self.rows('AC', 'text', sequence_column='sequence')

    def test_csv_tsv_columns_unicode_and_multiline(self):
        for delimiter, format in [(',', 'csv'), ('\t', 'tsv')]:
            for newline in ('\n', '\r\n', '\r'):
                text = newline.join(['\ufeffname' + delimiter + 'sequence' + delimiter + 'note',
                                     '"名' + delimiter + '前"' + delimiter + '"ac',
                                     'N T"' + delimiter + '"説明', '二行"', '',
                                     '次' + delimiter + 'GG' + delimiter + 'end', ''])
                stream = io.TextIOWrapper(io.BytesIO(text.encode()), encoding='utf-8', newline='')
                with stream:
                    rows = list(read_sequences(stream, format=format, sequence_column='sequence',
                                               name_column='name', description_column='note'))
                self.assertEqual([(r.id, r.sequence) for r in rows], [('名' + delimiter + '前', 'acNT'), ('次', 'GG')])
                self.assertEqual(rows[0].description, '説明' + newline + '二行')
                self.assertEqual([(r.provenance.record_index, r.provenance.start_line,
                                   r.provenance.end_line) for r in rows], [(1, 2, 4), (2, 6, 6)])
                self.assertTrue(all(r.quality is None for r in rows))

    def test_generated_ids_and_quoted_escape(self):
        rows = self.rows('sequence,unused\n"A C","say ""hello"""\nGT,x\n', sequence_column='sequence')
        self.assertEqual([(r.id, r.sequence) for r in rows], [('record-1', 'AC'), ('record-2', 'GT')])
        self.assertEqual(self.rows('sequence\n', sequence_column='sequence'), [])
        self.assertEqual(self.rows('', sequence_column='sequence'), [])

    def test_invalid_tables_have_source_and_line(self):
        cases = [('id,other\na,AC\n', 'missing column', 1),
                 ('sequence,sequence\nAC,GT\n', 'duplicate', 1),
                 ('sequence,\nAC,x\n', 'empty or duplicate', 1),
                 ('id,sequence\na\n', 'column count', 2),
                 ('id,sequence\na,AC,x\n', 'column count', 2),
                 ('id,sequence\na," \n "\n', 'empty tabular sequence', 2),
                 ('id,sequence\n,AC\n', 'empty tabular id', 2),
                 ('id,sequence\na,"AC\n', 'malformed csv', 2)]
        for text, reason, line in cases:
            with self.subTest(reason=reason), self.assertRaisesRegex(SequenceReadError, reason) as caught:
                self.rows(text, sequence_column='sequence', id_column='id', source='fixture')
            self.assertEqual((caught.exception.source, caught.exception.line), ('fixture', line))
        for options in ({}, {'sequence_column': ''}, {'sequence_column': 1},
                        {'sequence_column': 'sequence', 'id_column': 'id', 'name_column': 'name'},
                        {'sequence_column': 'sequence', 'id_column': ''},
                        {'sequence_column': 'sequence', 'description_column': 'missing'}):
            with self.subTest(options=options), self.assertRaises(SequenceReadError):
                self.rows('sequence\nAC\n', **options)

    def test_lazy_error_and_stream_lifetime(self):
        stream = NonSeek('sequence\nAC\n""\n')
        rows = read_sequences(stream, format='csv', sequence_column='sequence')
        self.assertEqual(next(rows).sequence, 'AC')
        with self.assertRaisesRegex(SequenceReadError, 'empty tabular sequence'):
            next(rows)
        self.assertFalse(stream.closed)

    def test_all_compression_formats_and_physical_extension(self):
        codecs = [('', lambda b: b, None), ('.gz', gzip.compress, 'gzip'),
                  ('.bz2', bz2.compress, 'bz2'), ('.xz', lzma.compress, 'xz'),
                  ('.lzma', lambda b: lzma.compress(b, format=lzma.FORMAT_ALONE), 'lzma')]
        with tempfile.TemporaryDirectory() as directory:
            for ext, text, options in [('txt', 'ac N\r\nGT\r\n', {}),
                                        ('csv', 'name,sequence\r\n名,"ac N\r\nGT"\r\n',
                                         {'sequence_column': 'sequence', 'name_column': 'name'}),
                                        ('tsv', 'name\tsequence\r\n名\t"ac N\r\nGT"\r\n',
                                         {'sequence_column': 'sequence', 'name_column': 'name'})]:
                for suffix, compress, codec in codecs:
                    with self.subTest(ext=ext, suffix=suffix):
                        path = Path(directory) / ('input.' + ext.upper() + suffix)
                        path.write_bytes(compress(text.encode('utf-8-sig')))
                        row, = read_sequences(path, source='logical.fa', **options)
                        self.assertEqual((row.sequence, row.source, row.provenance.format,
                                          row.provenance.compression), ('acNGT', 'logical.fa', ext if ext != 'txt' else 'text', codec))

    def test_tabular_field_limit_is_contextual_and_unchanged(self):
        limit = csv.field_size_limit()
        with self.assertRaisesRegex(SequenceReadError, 'malformed csv') as caught:
            self.rows('sequence\n' + 'A' * (limit + 1) + '\n',
                      sequence_column='sequence', source='large')
        self.assertEqual(caught.exception.source, 'large')
        self.assertEqual(caught.exception.line, 2)
        self.assertEqual(csv.field_size_limit(), limit)


if __name__ == '__main__':
    unittest.main()
