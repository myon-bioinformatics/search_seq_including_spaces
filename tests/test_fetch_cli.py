"""Exercise CLI parsing, real acquisition dispatch and output with offline I/O."""
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import io
import json
import unittest
from unittest.mock import Mock, patch

from acquisition import Client
from sequence_tool import main


class FetchCliTests(unittest.TestCase):
    def invoke(self, args, payload=b'', error=None):
        transport = Mock(return_value=payload, side_effect=error)
        def client(**kwargs):
            return Client(transport=transport, retries=0, **kwargs)
        buffer = io.BytesIO()
        stdout = io.TextIOWrapper(buffer, encoding='utf-8')
        stderr = io.StringIO()
        with patch('sequence_commands.Client', side_effect=client), \
                redirect_stdout(stdout), redirect_stderr(stderr):
            try:
                code = main(['fetch', *args])
            except SystemExit as exc:
                code = exc.code
        stdout.flush()
        return code, buffer.getvalue(), stderr.getvalue(), transport

    def test_structure_cif_pdb_and_default(self):
        for fmt in (None, 'cif', 'pdb'):
            with self.subTest(format=fmt):
                expected = fmt or 'cif'
                payload = b'data_1CRN\n#\n' if expected == 'cif' else b'HEADER    TEST\nEND\n'
                args = ['--operation', 'records', '--database', 'structure',
                        '--ids', '1crn', '--timeout', '2', '--max-bytes', '1000']
                if fmt:
                    args += ['--format', fmt]
                code, out, err, transport = self.invoke(args, payload)
                self.assertEqual((code, err), (0, ''))
                result = json.loads(out)
                url = f'https://files.rcsb.org/download/1CRN.{expected}'
                transport.assert_called_once_with(url, timeout=2, max_bytes=1000)
                self.assertEqual(result['schema'], 'structure-download/1')
                self.assertEqual((result['id'], result['format']), ('1CRN', expected))
                self.assertEqual(result['text'], payload.decode())
                self.assertEqual(result['source'], url)
                self.assertEqual(result['sha256'], hashlib.sha256(payload).hexdigest())
                self.assertTrue(result['mmdb_lookup'].endswith('1CRN'))

    def test_invalid_records_fail_before_transport(self):
        cases = [
            ('structure', []),
            ('structure', ['--ids', '1CRN', '2ABC']),
            ('structure', ['--ids', '../secret']),
            ('structure', ['--ids', '1CRN', '--query', 'x']),
            *[('structure', ['--ids', '1CRN', '--format', fmt])
              for fmt in ('fasta', 'genbank', 'xml')],
            *[(db, ['--ids', 'NM_1.2', '--format', fmt])
              for db in ('nuccore', 'protein') for fmt in ('cif', 'pdb')],
            ('pubmed', ['--ids', '123', '--format', 'fasta']),
        ]
        for db, extra in cases:
            with self.subTest(database=db, args=extra):
                code, out, err, transport = self.invoke(
                    ['--operation', 'records', '--database', db, *extra])
                self.assertEqual((code, out), (2, b''))
                self.assertIn('error:', err)
                self.assertNotIn('Traceback', err)
                transport.assert_not_called()

    def test_sequence_and_pubmed_records_unchanged(self):
        for db in ('nuccore', 'protein'):
            for fmt in (None, 'fasta', 'genbank'):
                with self.subTest(database=db, format=fmt):
                    payload = (b'LOCUS       TEST 4 bp DNA\nACCESSION   TEST\n'
                               b'ORIGIN\n        1 acgt\n//\n'
                               if fmt == 'genbank' else b'>TEST\nACGT\n')
                    args = ['--operation', 'records', '--database', db, '--ids', 'TEST']
                    if fmt:
                        args += ['--format', fmt]
                    code, out, err, transport = self.invoke(args, payload)
                    self.assertEqual((code, err), (0, ''))
                    result = json.loads(out)
                    self.assertEqual(result['schema'], 'ncbi-sequences/1')
                    self.assertEqual(result['records'][0]['sequence'],
                                     'acgt' if fmt == 'genbank' else 'ACGT')
                    rettype = ('gb' if db == 'nuccore' else 'gp') if fmt == 'genbank' else 'fasta'
                    self.assertIn('rettype=' + rettype, transport.call_args.args[0])
        code, out, err, transport = self.invoke(
            ['--operation', 'records', '--database', 'pubmed', '--ids', '123'],
            b'<PubmedArticleSet/>')
        self.assertEqual((code, err), (0, ''))
        self.assertEqual(json.loads(out)['articles'], [])
        self.assertIn('efetch.fcgi?db=pubmed', transport.call_args.args[0])

    def test_structure_search_summary_keep_ncbi_dispatch(self):
        for operation, extra, payload in (
            ('search', ['--query', 'crambin'],
             b'{"esearchresult":{"count":"1","idlist":["123"]}}'),
            ('summary', ['--ids', '123'], b'{"result":{"uids":["123"]}}'),
        ):
            with self.subTest(operation=operation):
                args = ['--operation', operation, '--database', 'structure', *extra]
                code, out, err, transport = self.invoke(args, payload)
                self.assertEqual((code, err), (0, ''))
                self.assertEqual(json.loads(out)['database'], 'structure')
                self.assertIn(f'/e{operation}.fcgi?db=structure', transport.call_args.args[0])
                code, out, err, transport = self.invoke(args + ['--format', 'cif'])
                self.assertEqual((code, out), (2, b''))
                transport.assert_not_called()

    def test_structure_transport_failure_and_size_limit(self):
        args = ['--operation', 'records', '--database', 'structure', '--ids', '1CRN']
        code, out, err, transport = self.invoke(args, error=TimeoutError('private detail'))
        self.assertEqual((code, out), (2, b''))
        self.assertIn('fetch failed (TimeoutError)', err)
        self.assertNotIn('private detail', err)
        self.assertEqual(transport.call_count, 1)
        code, out, err, transport = self.invoke(args + ['--max-bytes', '2'], b'abc')
        self.assertEqual((code, out), (2, b''))
        self.assertIn('oversized', err)
