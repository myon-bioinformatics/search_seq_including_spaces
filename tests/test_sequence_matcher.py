import bz2
from contextlib import closing
from dataclasses import asdict
import gc
import gzip
import io
import json
import lzma
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import weakref

from anchor_matcher import find_matches
from sequence_matcher import find_sequence_matches
from sequence_readers import read_sequences, SequenceReadError
from sequence_records import ReadProvenance, SequenceRecord


class GuardedStream:
    """Non-seek fixture permitting only explicitly requested physical lines."""
    def __init__(self, text):
        self.lines = iter(text.splitlines(keepends=True))
        self.allowed = self.consumed = 0
        self.closed = False

    def __iter__(self):
        return self

    def __next__(self):
        if self.consumed >= self.allowed:
            raise AssertionError('adapter read ahead of the requested record')
        line = next(self.lines)
        self.consumed += 1
        return line

    def seekable(self):
        return False

    def read(self, *args):
        raise AssertionError('no read-all/sniff on explicit non-seek input')

    readline = seek = read

    def close(self):
        self.closed = True


class SequenceMatcherTests(unittest.TestCase):
    def test_fasta_envelope_coordinates_and_duplicate_ids(self):
        text = '\ufeff;leading\r\n>same 説明\r\n ac N\r\nGT\r\n;tail\r\n\r\n>same next\r\nCC\r\n'
        envelopes = list(find_sequence_matches(io.StringIO(text), 'AC', 'GT',
                         format='fasta', source='logical', min_gap=1))
        records = list(read_sequences(io.StringIO(text), format='fasta', source='logical'))
        for envelope, record in zip(envelopes, records):
            self.assertEqual(set(envelope), {'schema', 'record', 'result'})
            self.assertEqual(envelope['schema'], 'sequence-matches/1')
            self.assertEqual(envelope['record'], {
                'id': record.id, 'description': record.description,
                'source': record.source, 'provenance': asdict(record.provenance)})
            self.assertEqual(envelope['result'], find_matches(record.sequence, 'AC', 'GT', min_gap=1))
            self.assertEqual(json.loads(json.dumps(envelope)), envelope)
        self.assertEqual(len(envelopes), 2)
        first, second = envelopes
        self.assertEqual(first['record']['provenance'], dict(format='fasta', compression=None,
                         record_index=1, start_line=2, end_line=6))
        self.assertEqual(second['record']['provenance']['record_index'], 2)
        self.assertEqual(first['result']['matches'][0]['source_span'], [1, 5])
        self.assertEqual(first['result']['matches'][0]['middle_source_span'], [3, 3])
        self.assertEqual(second['result']['matches'], [])
        self.assertTrue(second['result']['complete'])
        first['record']['provenance']['record_index'] = 999
        self.assertEqual(second['record']['provenance']['record_index'], 2)

    def test_fastq_wrapped_quality_is_not_matched_or_exported(self):
        text = '@one desc\nac\nNGT\n+one desc\n@@\n+!!\n@two\nacNGT\n+\n!!!!!\n'
        first, second = find_sequence_matches(io.StringIO(text), 'AC', 'GT',
                                             format='fastq', min_gap=1)
        self.assertEqual(first['result'], find_matches('acNGT', 'AC', 'GT', min_gap=1))
        self.assertEqual(first['result'], second['result'])
        self.assertIsNone(first['record']['source'])
        self.assertEqual(first['record']['provenance']['end_line'], 6)
        self.assertEqual(second['record']['provenance']['start_line'], 7)
        self.assertNotIn('sequence', first['record'])
        self.assertNotIn('quality', first['record'])

    def test_plain_and_compressed_paths_all_supported_formats(self):
        codecs = [('', lambda b: b, None), ('.gz', gzip.compress, 'gzip'),
                  ('.bz2', bz2.compress, 'bz2'), ('.xz', lzma.compress, 'xz'),
                  ('.lzma', lambda b: lzma.compress(b, format=lzma.FORMAT_ALONE), 'lzma')]
        cases = [('fasta', '>id\r\nacN\r\nGT\r\n', {}),
                 ('fastq', '@id\r\nacNGT\r\n+\r\n!!!!!\r\n', {}),
                 ('txt', ' acN\r\nGT\r\n', {}),
                 ('csv', 'seq\r\nacNGT\r\n', {'sequence_column': 'seq'}),
                 ('tsv', 'seq\r\nacNGT\r\n', {'sequence_column': 'seq'})]
        with tempfile.TemporaryDirectory() as directory:
            for extension, text, options in cases:
                for suffix, compress, compression in codecs:
                    with self.subTest(extension=extension, suffix=suffix):
                        path = Path(directory) / ('input.' + extension + suffix)
                        path.write_bytes(compress(text.encode('utf-8')))
                        envelope, = find_sequence_matches(path, 'AC', 'GT', min_gap=1, **options)
                        self.assertEqual(envelope['record']['source'], str(path))
                        self.assertEqual(envelope['record']['provenance']['compression'], compression)
                        self.assertEqual(envelope['result'], find_matches('acNGT', 'AC', 'GT', min_gap=1))
                        alias, = find_sequence_matches(path, 'AC', 'GT', min_gap=1,
                                                       source='logical.wrong', **options)
                        self.assertEqual(alias['record']['source'], 'logical.wrong')
                        self.assertEqual(alias['result'], envelope['result'])

    def test_csv_tsv_multiline_unicode_and_identity_alias(self):
        for format, delimiter in [('csv', ','), ('tsv', '\t')]:
            text = delimiter.join(['name', 'seq', 'notes']) + '\n'
            text += delimiter.join(['試料', '"ac\nNGT"', '"説明\n続き"']) + '\n'
            text += delimiter.join(['次', 'CC', 'none']) + '\n'
            for identity_option in ('id_column', 'name_column'):
                with self.subTest(format=format, identity_option=identity_option):
                    rows = list(find_sequence_matches(io.StringIO(text), 'AC', 'GT',
                                format=format, sequence_column='seq', description_column='notes',
                                min_gap=1, **{identity_option: 'name'}))
                    self.assertEqual([r['record']['id'] for r in rows], ['試料', '次'])
                    self.assertEqual(rows[0]['record']['description'], '説明\n続き')
                    self.assertEqual(rows[0]['record']['provenance']['start_line'], 2)
                    self.assertEqual(rows[0]['record']['provenance']['end_line'], 4)
                    self.assertEqual(rows[0]['result']['matches'][0]['source_span'], [1, 5])
                    self.assertEqual(rows[1]['result']['matches'], [])

    def test_records_are_independent_not_concatenated(self):
        rows = list(find_sequence_matches(io.StringIO('>left\nAC\n>right\nGT\n'),
                                         'AC', 'GT', format='fasta'))
        self.assertEqual([r['result']['matches'] for r in rows], [[], []])

    def test_all_matcher_options_and_reverse_coordinates_are_unchanged(self):
        cases = [('TTACG', 'CG', 'AA', {'min_gap': 1, 'strand': '-'}),
                 ('AT', 'A', 'T', {'strand': 'both'}),
                 ('ART', 'AG', 'T', {'iupac': True}),
                 ('AAAA', 'A', 'A', {'min_gap': 0, 'max_gap': 1, 'max_matches': 5}),
                 ('MQKLV', 'MQ', 'LV', {'min_gap': 1})]
        for sequence, left, right, options in cases:
            with self.subTest(sequence=sequence, options=options):
                row, = find_sequence_matches(io.StringIO(sequence), left, right,
                                             format='text', **options)
                self.assertEqual(row['result'], find_matches(sequence, left, right, **options))
                self.assertEqual(row['result']['schema'], 'anchor-gap/1')
                self.assertTrue(row['result']['matches'])
        self.assertEqual(find_matches('TTACG', 'CG', 'AA', min_gap=1,
                                     strand='-')['matches'][0]['left_span'], [4, 5])

    def test_lazy_nonseek_input_no_next_record_prefetch(self):
        stream = GuardedStream('@one\nAC\n+\n!!\n@two\nAC\n+\n!!\n')
        rows = find_sequence_matches(stream, 'A', 'C', format='fastq')
        self.assertEqual(stream.consumed, 0)
        stream.allowed = 4
        self.assertEqual(next(rows)['record']['id'], 'one')
        self.assertEqual(stream.consumed, 4)
        stream.allowed = 8
        self.assertEqual(next(rows)['record']['id'], 'two')
        self.assertEqual(stream.consumed, 8)
        rows.close()
        self.assertFalse(stream.closed)

    def test_reader_error_after_valid_record_keeps_source_and_line(self):
        cases = [('fastq', '@ok\nAC\n+\n!!\nbroken\n', {}, 5),
                 ('fasta', '>ok\nAC\n>\nGT\n', {}, 3),
                 ('csv', 'seq\nAC\n""\n', {'sequence_column': 'seq'}, 3)]
        for format, text, options, line in cases:
            with self.subTest(format=format):
                stream = io.StringIO(text)
                rows = find_sequence_matches(stream, 'A', 'C', format=format,
                                             source='bad-input', **options)
                first = next(rows)
                self.assertTrue(first['result']['complete'])
                self.assertNotIn('complete', first)  # no whole-file success claim
                with self.assertRaises(SequenceReadError) as caught:
                    next(rows)
                self.assertEqual(caught.exception.source, 'bad-input')
                self.assertEqual(caught.exception.line, line)
                with self.assertRaises(StopIteration):
                    next(rows)
                self.assertFalse(stream.closed)

    def test_limit_resets_for_each_record_and_exact_limit_succeeds(self):
        rows = list(find_sequence_matches(io.StringIO('>a\nAAA\n>b\nAAA\n'),
                                         'A', 'A', format='fasta', max_matches=2))
        self.assertEqual([len(r['result']['matches']) for r in rows], [2, 2])

    def test_limit_overflow_yields_no_partial_record(self):
        rows = find_sequence_matches(io.StringIO('>a\nAAAA\n'), 'A', 'A',
                                     format='fasta', max_matches=2)
        with self.assertRaisesRegex(ValueError, 'match limit exceeded'):
            next(rows)
        with self.assertRaises(StopIteration):
            next(rows)

    def test_later_limit_error_preserves_prior_record_only(self):
        stream = GuardedStream('@ok\nAA\n+\n!!\n@many\nAAA\n+\n!!!\nNEVER READ\n')
        stream.allowed = 8
        rows = find_sequence_matches(stream, 'A', 'A', format='fastq', max_matches=1)
        first = next(rows)
        with self.assertRaisesRegex(ValueError, 'match limit exceeded'):
            next(rows)
        self.assertEqual(first['record']['id'], 'ok')
        self.assertEqual(len(first['result']['matches']), 1)
        self.assertEqual(stream.consumed, 8)
        self.assertFalse(stream.closed)

    def test_invalid_matcher_options_still_raise(self):
        for options in [{'max_matches': 0}, {'max_matches': True}, {'min_gap': -1},
                        {'min_gap': 2, 'max_gap': 1}, {'iupac': 1}, {'strand': 'unknown'}]:
            with self.subTest(options=options), self.assertRaises(ValueError):
                next(find_sequence_matches(io.StringIO('AC'), 'A', 'C', format='text', **options))
        with self.assertRaisesRegex(ValueError, 'anchors must not be empty'):
            next(find_sequence_matches(io.StringIO('AC'), '', 'C', format='text'))

    def test_invalid_sequence_is_not_silently_normalized_by_adapter(self):
        with self.assertRaisesRegex(ValueError, 'ASCII letters'):
            next(find_sequence_matches(io.StringIO('>id\nA-C\n'), 'A', 'C', format='fasta'))

    def test_empty_inputs_yield_no_envelopes(self):
        for format in ('fasta', 'fastq', 'text', 'csv', 'tsv'):
            options = {'sequence_column': 'seq'} if format in ('csv', 'tsv') else {}
            with self.subTest(format=format):
                self.assertEqual(list(find_sequence_matches(io.StringIO(''), 'A', 'C',
                                     format=format, **options)), [])

    def test_reader_detection_and_errors_are_delegated(self):
        row, = find_sequence_matches(io.StringIO('>id\nAC\n'), 'A', 'C')
        self.assertEqual(row['record']['provenance']['format'], 'fasta')
        row, = find_sequence_matches(io.StringIO('>id\nAC\n'), 'A', 'C',
                                     format='fasta', source='wrong.fastq')
        self.assertEqual(row['record']['source'], 'wrong.fastq')
        for options in ({'format': 'json'}, {'format': 'jsonl'},
                        {'format': 'csv'}, {'format': 'fasta', 'sequence_column': 'seq'}):
            with self.subTest(options=options), self.assertRaises(SequenceReadError):
                next(find_sequence_matches(io.StringIO('>id\nAC\n'), 'A', 'C', **options))

    def test_stream_compression_label_does_not_invent_compression(self):
        row, = find_sequence_matches(io.StringIO('>id\nAC\n'), 'A', 'C', source='caller.fa.gz')
        self.assertIsNone(row['record']['provenance']['compression'])

    def test_path_is_opened_lazily_and_io_errors_propagate(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'missing.fa'
            rows = find_sequence_matches(path, 'A', 'C')
            rows.close()  # closing an unstarted iterator must not open the path
            rows = find_sequence_matches(path, 'A', 'C')
            with self.assertRaises(FileNotFoundError):
                next(rows)

    def test_owned_file_closes_on_exhaustion_close_and_errors(self):
        real_open = open
        cases = [('exhaust', '@ok\nAC\n+\n!!\n'),
                 ('close', '@ok\nAC\n+\n!!\n@later\nAC\n+\n!!\n'),
                 ('reader-error', '@ok\nAC\n+\n!!\nbroken\n'),
                 ('matcher-error', '@ok\nAC\n+\n!!\n@bad\nA-C\n+\n!!!\n')]
        with tempfile.TemporaryDirectory() as directory:
            for action, text in cases:
                with self.subTest(action=action):
                    path = Path(directory) / 'input.fq'
                    path.write_text(text, encoding='utf-8')
                    handles = []
                    def tracked_open(*args, **kwargs):
                        handle = real_open(*args, **kwargs)
                        handles.append(handle)
                        return handle
                    with patch('sequence_readers.open', tracked_open, create=True):
                        rows = find_sequence_matches(path, 'A', 'C')
                        self.assertEqual(handles, [])
                        next(rows)
                        self.assertFalse(handles[0].closed)
                        if action == 'close':
                            rows.close()
                        elif action == 'exhaust':
                            with self.assertRaises(StopIteration):
                                next(rows)
                        else:
                            with self.assertRaises(ValueError):
                                next(rows)
                        self.assertTrue(handles[0].closed)

    def test_caller_stream_remains_open_after_exhaustion_or_early_exit(self):
        stream = io.StringIO('>id\nAC\n')
        list(find_sequence_matches(stream, 'A', 'C', format='fasta'))
        self.assertFalse(stream.closed)
        stream.seek(0)
        with self.assertRaisesRegex(RuntimeError, 'consumer'):
            with closing(find_sequence_matches(stream, 'A', 'C', format='fasta')) as rows:
                next(rows)
                raise RuntimeError('consumer')
        self.assertFalse(stream.closed)

    def test_previous_record_and_result_released_before_next_read(self):
        record_refs, result_refs = [], []
        class WeakResult(dict):
            pass
        def reader(*args, **kwargs):
            for index in range(1, 4):
                gc.collect()
                self.assertTrue(all(ref() is None for ref in record_refs))
                self.assertTrue(all(ref() is None for ref in result_refs))
                record = SequenceRecord(str(index), '', 'AC', None,
                                        ReadProvenance('fasta', None, index, index, index))
                record_refs.append(weakref.ref(record))
                yield record
                del record
        def matcher(*args, **kwargs):
            result = WeakResult(find_matches(*args, **kwargs))
            result_refs.append(weakref.ref(result))
            return result
        with patch('sequence_matcher.read_sequences', reader), patch('sequence_matcher.find_matches', matcher):
            rows = find_sequence_matches(None, 'A', 'C')
            for index in range(1, 4):
                envelope = next(rows)
                self.assertIs(envelope['result'], result_refs[-1]())
                self.assertEqual(envelope['record']['id'], str(index))
                del envelope
            with self.assertRaises(StopIteration):
                next(rows)
        gc.collect()
        self.assertTrue(all(ref() is None for ref in record_refs + result_refs))


if __name__ == '__main__':
    unittest.main()
