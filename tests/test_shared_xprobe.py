"""Consumer integrity check; helper behavior tests belong to upstream xprobe."""
import hashlib
import json
from pathlib import Path
import unittest


class SharedXprobeTests(unittest.TestCase):
    def test_pinned_runtime_helper(self):
        root = Path(__file__).resolve().parents[1]
        data = (root / 'vendor/xprobe.py').read_bytes()
        meta = json.loads((root / 'vendor/xprobe.provenance.json').read_text())
        self.assertEqual(meta['repository'], 'myon-bioinformatics/xprobe')
        self.assertEqual(meta['commit'], '326acd667e13b21bf53ccc1590af960edf8cbf6c')
        self.assertEqual(meta['path'], 'xprobe.py')
        self.assertEqual(meta['blob_sha'], 'dbc5b7d55005d6288c072a7612584d6170c216f4')
        self.assertEqual(hashlib.sha256(data).hexdigest(), meta['sha256'])
        self.assertEqual(hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest(), meta['blob_sha'])
        self.assertTrue((root / 'vendor' / meta['license']).is_file())
        from vendor import xprobe
        self.assertTrue(callable(xprobe.corpus_from_json))
