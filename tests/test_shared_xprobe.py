"""Consumer integrity check; helper behavior tests belong to upstream xprobe."""
import json
from pathlib import Path
import unittest


def _lock_entry(root: Path, destination: str) -> dict:
    lock = json.loads((root / "vendor.lock.json").read_text(encoding="utf-8"))
    for entry in lock["files"]:
        if entry["destination"] == destination:
            return entry
    raise AssertionError(f"missing lock entry for {destination}")


class SharedXprobeTests(unittest.TestCase):
    def test_pinned_runtime_helper(self):
        root = Path(__file__).resolve().parents[1]
        meta = _lock_entry(root, "vendor/xprobe.py")
        self.assertEqual(meta["repository"], "myon-bioinformatics/xprobe")
        self.assertEqual(meta["commit"], "642999cea4185a68bffa7f7ccc46bd78dde03e5a")
        self.assertEqual(meta["source"], "xprobe.py")
        self.assertTrue((root / "vendor" / "xprobe-LICENSE").is_file())
        from vendor import xprobe
        self.assertTrue(callable(xprobe.pytest_receipt))
