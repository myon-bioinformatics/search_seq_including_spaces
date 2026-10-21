"""Native sequence discovery uses reader format/compression rules, not imports."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

from sequence_readers import discover_sequences

ROOT = Path(__file__).resolve().parents[1]


def inputs(root):
    for name in ('sample.fa', 'sample.FASTQ.GZ', 'other.fasta.bz2', 'data.jsonl.xz',
                 'unrelated.png', '.hidden.fa', 'nested/sample.gbk'):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        # Deliberately invalid compressed contents: discovery must not open them.
        path.write_bytes(b'not a readable sequence container')


def test_formats_compression_and_injected_names(tmp_path):
    inputs(tmp_path)
    report = discover_sequences(tmp_path, formats=('fasta', 'fastq'), patterns=('sample*', 'other*'))
    assert report == {'files': [
        {'path': 'other.fasta.bz2', 'format': 'fasta', 'compression': 'bz2'},
        {'path': 'sample.FASTQ.GZ', 'format': 'fastq', 'compression': 'gzip'},
        {'path': 'sample.fa', 'format': 'fasta', 'compression': None},
    ], 'truncated': False}
    assert discover_sequences(tmp_path, formats=('genbank',))['files'][0]['path'] == 'nested/sample.gbk'
    assert discover_sequences(tmp_path, formats=('jsonl',))['files'][0]['compression'] == 'xz'
    assert discover_sequences(tmp_path, formats=('genbank',), recursive=False)['files'] == []
    assert len(discover_sequences(tmp_path, formats=('fasta',), include_hidden=True)['files']) == 3


def test_limit_is_global_exact_and_filters_before_counting(tmp_path):
    inputs(tmp_path)
    report = discover_sequences(tmp_path, formats=('fasta',), max_files=2)
    assert len(report['files']) == 2 and not report['truncated']
    assert discover_sequences(tmp_path, max_files=2)['truncated']
    assert discover_sequences(tmp_path, patterns=('nested/*',), max_files=1)['files'][0]['format'] == 'genbank'
    assert discover_sequences(tmp_path, patterns=())['files'] == []


@pytest.mark.parametrize('options,error', [
    ({'formats': ('unknown',)}, ValueError), ({'formats': 'fasta'}, TypeError),
    ({'patterns': '*.fa'}, TypeError), ({'patterns': ('',)}, ValueError),
    ({'max_files': 0}, ValueError), ({'max_files': True}, ValueError),
])
def test_bad_rules_rejected_on_empty_root(tmp_path, options, error):
    with pytest.raises(error):
        discover_sequences(tmp_path, **options)


def test_links_and_missing_root(tmp_path):
    inputs(tmp_path)
    with pytest.raises(ValueError):
        discover_sequences(tmp_path / 'missing')
    link = tmp_path / 'link.fa'
    try:
        link.symlink_to(tmp_path / 'sample.fa')
        (tmp_path / 'loop').symlink_to(tmp_path, target_is_directory=True)
    except (OSError, NotImplementedError):
        pytest.skip('symlinks unavailable')
    assert 'link.fa' not in [row['path'] for row in discover_sequences(tmp_path)['files']]
    with pytest.raises(ValueError):
        discover_sequences(tmp_path / 'loop')


def test_native_cli_without_site_packages_and_truncated_status(tmp_path):
    inputs(tmp_path)
    args = [sys.executable, '-S', str(ROOT / 'sequence_tool.py'), 'discover', str(tmp_path),
            '--format', 'fasta', '--pattern', 'sample*', '--pattern', 'other*']
    result = subprocess.run(args, capture_output=True, text=True, encoding='utf-8', timeout=20)
    assert result.returncode == 0, result.stderr
    assert len(json.loads(result.stdout)['files']) == 2
    result = subprocess.run(args + ['--max-files', '1'], capture_output=True, text=True, encoding='utf-8', timeout=20)
    assert result.returncode == 2 and json.loads(result.stdout)['truncated']


def test_walk_errors_are_not_an_empty_success(tmp_path, monkeypatch):
    def failed(root, **options):
        options['onerror'](PermissionError('denied'))
        return iter(())
    monkeypatch.setattr('sequence_readers.os.walk', failed)
    with pytest.raises(PermissionError):
        discover_sequences(tmp_path)
