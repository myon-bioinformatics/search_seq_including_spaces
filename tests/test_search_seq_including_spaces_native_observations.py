"""Downstream integration of the pinned shared pytest adapter."""
import os, subprocess
from pathlib import Path

from vendor.xprobe import pytest_receipt

ROOT = Path(__file__).resolve().parents[1]


def test_native_failure_evidence(tmp_path):
    suite = tmp_path / 'test_sample.py'
    # Keep only the reader -> pytest -> shared-validator integration here.
    # Generic xfail/XPASS/setup/teardown and corrupt receipts are upstream tests.
    suite.write_text('''import io
from seqlib.sequence_readers import read_sequences
def test_failure():
    record, = read_sequences(io.StringIO(">fixture\\nAC\\n"), format="fasta")
    assert record.sequence == "AG", "private diagnostic"
''', encoding='utf-8')
    destination = tmp_path / 'events.jsonl'
    env = dict(os.environ, PYTHONPATH=str(ROOT), PYTEST_DISABLE_PLUGIN_AUTOLOAD='1')
    env.pop('PYTEST_ADDOPTS', None)
    result = subprocess.run(
        ['uv', 'run', '--locked', '--project', str(ROOT), 'pytest', '-c', os.devnull,
         '-p', 'vendor.xprobe_pytest',
         '--xprobe-jsonl=' + str(destination), '--xprobe-repository=myon-bioinformatics/search_seq_including_spaces',
         '--xprobe-run-id=reader-integration',
         str(suite)], cwd=tmp_path, env=env, capture_output=True, text=True, timeout=60)
    assert result.returncode == 1, result.stdout + result.stderr
    text = destination.read_text(encoding='utf-8')
    receipt = pytest_receipt(text, expected_context={
        'repository': 'myon-bioinformatics/search_seq_including_spaces',
        'commit_sha': None, 'report_id': 'reader-integration'})
    assert receipt['exitstatus'] == result.returncode
    assert receipt['phase_counts']['call'] == {'failed': 1}
    assert 'private' not in text
