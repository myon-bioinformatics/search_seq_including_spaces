"""Native discovery contracts run in both pytest and stdlib-only environments."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from sequence_readers import discover_sequences

ROOT = Path(__file__).resolve().parents[1]


class SequenceDiscoveryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def inputs(self):
        for name in ('sample.fa', 'sample.FASTQ.GZ', 'other.fasta.bz2', 'data.jsonl.xz',
                     'unrelated.png', '.hidden.fa', 'nested/sample.gbk'):
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            # Invalid compressed contents: discovery must not open them.
            path.write_bytes(b'not a readable sequence container')

    def test_formats_compression_and_injected_names(self):
        self.inputs()
        report = discover_sequences(self.root, formats=('fasta', 'fastq'), patterns=('sample*', 'other*'))
        self.assertEqual(report, {'files': [
            {'path': 'other.fasta.bz2', 'format': 'fasta', 'compression': 'bz2'},
            {'path': 'sample.FASTQ.GZ', 'format': 'fastq', 'compression': 'gzip'},
            {'path': 'sample.fa', 'format': 'fasta', 'compression': None},
        ], 'truncated': False})
        self.assertEqual(discover_sequences(self.root, formats=('genbank',))['files'][0]['path'], 'nested/sample.gbk')
        self.assertEqual(discover_sequences(self.root, formats=('jsonl',))['files'][0]['compression'], 'xz')
        self.assertEqual(discover_sequences(self.root, formats=('genbank',), recursive=False)['files'], [])
        self.assertEqual(len(discover_sequences(self.root, formats=('fasta',), include_hidden=True)['files']), 3)

    def test_limit_is_global_exact_and_filters_before_counting(self):
        self.inputs()
        report = discover_sequences(self.root, formats=('fasta',), max_files=2)
        self.assertEqual(len(report['files']), 2)
        self.assertFalse(report['truncated'])
        self.assertTrue(discover_sequences(self.root, max_files=2)['truncated'])
        self.assertEqual(discover_sequences(self.root, patterns=('nested/*',), max_files=1)['files'][0]['format'], 'genbank')
        self.assertEqual(discover_sequences(self.root, patterns=())['files'], [])

    def test_bad_rules_rejected_on_empty_root(self):
        for options, error in [
            ({'formats': ('unknown',)}, ValueError), ({'formats': 'fasta'}, TypeError),
            ({'patterns': '*.fa'}, TypeError), ({'patterns': ('',)}, ValueError),
            ({'max_files': 0}, ValueError), ({'max_files': True}, ValueError),
        ]:
            with self.subTest(options=options), self.assertRaises(error):
                discover_sequences(self.root, **options)

    def test_links_and_missing_root(self):
        self.inputs()
        with self.assertRaises(ValueError):
            discover_sequences(self.root / 'missing')
        link = self.root / 'link.fa'
        try:
            link.symlink_to(self.root / 'sample.fa')
            (self.root / 'loop').symlink_to(self.root, target_is_directory=True)
        except (OSError, NotImplementedError):
            self.skipTest('symlinks unavailable')
        self.assertNotIn('link.fa', [row['path'] for row in discover_sequences(self.root)['files']])
        with self.assertRaises(ValueError):
            discover_sequences(self.root / 'loop')

    def test_native_cli_without_site_packages_and_truncated_status(self):
        self.inputs()
        args = [sys.executable, '-S', str(ROOT / 'sequence_tool.py'), 'discover', str(self.root),
                '--format', 'fasta', '--pattern', 'sample*', '--pattern', 'other*']
        result = subprocess.run(args, capture_output=True, text=True, encoding='utf-8', timeout=20)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(json.loads(result.stdout)['files']), 2)
        result = subprocess.run(args + ['--max-files', '1'], capture_output=True, text=True, encoding='utf-8', timeout=20)
        self.assertEqual(result.returncode, 2)
        self.assertTrue(json.loads(result.stdout)['truncated'])

    def test_walk_errors_are_not_an_empty_success(self):
        def failed(root, **options):
            options['onerror'](PermissionError('denied'))
            return iter(())
        with patch('sequence_readers.os.walk', failed), self.assertRaises(PermissionError):
            discover_sequences(self.root)
