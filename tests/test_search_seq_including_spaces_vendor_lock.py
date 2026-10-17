"""Checked-in vendor bytes must match vendor.lock.json."""
import hashlib
import json
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class VendorLockTests(unittest.TestCase):
    def test_vendor_lock_matches_checked_in_bytes(self) -> None:
        lock = json.loads((ROOT / "vendor.lock.json").read_text(encoding="utf-8"))
        self.assertEqual(lock["schema"], "vendor-lock/1")
        for entry in lock["files"]:
            self.assertRegex(entry["ref"], r"^[0-9a-f]{40}$", entry["destination"])
            self.assertEqual(entry["ref"], entry["commit"], entry["destination"])
            data = (ROOT / entry["destination"]).read_bytes()
            blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
            self.assertEqual(blob, entry["blob_sha"], entry["destination"])
            self.assertEqual(hashlib.sha256(data).hexdigest(), entry["sha256"], entry["destination"])
