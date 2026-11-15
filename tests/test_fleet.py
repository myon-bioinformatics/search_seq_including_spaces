"""Fleet contracts through the parent's shared checks (vendor/fleet_check.py): vendor lock and
xprobe provenance, parent pins, CI Python setup, benchmark guard; stdlib unittest, also under pytest.

Helper behavior tests belong to upstream xprobe; this repository checks provenance and consumer integration.
"""
import re
import unittest

from vendor import fleet_check as fleet

ROOT = fleet.repo_root(__file__)
PARENT_WORKFLOW_PIN = "9e3da7ae1846e9604ebc24bb5178f7de094a4e77"
# Node sources and the shared checks stay on accepted upstream commits with exact bytes.
ACCEPTED = {
    "vendor/xprobe.py": "1ffbb7b984075c683690bbe4e9deb06199cdfe94",
    "vendor/xprobe_pytest.py": "326acd667e13b21bf53ccc1590af960edf8cbf6c",
    "vendor/fleet_check.py": "6c29503cf17eda856766ddfdf4d66a5b0965a8b1",
}


class FleetTests(unittest.TestCase):
    def test_vendor_lock_matches_checked_in_bytes_and_accepted_commits(self):
        fleet.check_lock(ROOT, commits=ACCEPTED)

    def test_pinned_runtime_helper(self):
        row = fleet.check_lock(ROOT, commits=ACCEPTED)["vendor/xprobe.py"]
        self.assertEqual((row["repository"], row["source"]), ("myon-bioinformatics/xprobe", "xprobe.py"))
        self.assertTrue((ROOT / "vendor" / "xprobe-LICENSE").is_file())
        from vendor import xprobe
        self.assertTrue(callable(xprobe.pytest_receipt))

    def test_adapter_provenance(self):
        # check_lock verifies ref == commit and the placed bytes against blob_sha / sha256.
        row = fleet.check_lock(ROOT, commits=ACCEPTED)["vendor/xprobe_pytest.py"]
        self.assertEqual((row["repository"], row["source"]), ("myon-bioinformatics/xprobe", "scripts/xprobe_pytest.py"))
        self.assertEqual((row["sha256"], row["blob_sha"]),
                         ("c0ea71f9d971bf47a9f5f1794ef8a7f641dd17694ea8e7f29795d326c5229285",
                          "70fac53151e97f5f28d23f9961d3837b7a885ea8"))

    def test_parent_workflow_references_use_reviewed_rollout_commit(self):
        self.assertEqual(set(fleet.parent_refs(ROOT)), {PARENT_WORKFLOW_PIN})

    def test_setup_precedes_first_plain_python_command(self):
        job = fleet.workflow_job(ROOT, "test")
        command = re.search(r"(?m)^\s+(?:-?\s*run:\s*)?python(?:3)?\s", job)
        self.assertIsNotNone(command)
        self.assertLess(job.index("uses: astral-sh/setup-uv@"), command.start(), "Configure matrix Python before stdlib checks")
        self.assertIn('python-version: "${{ matrix.python-version }}"', job)
        self.assertIn("activate-environment: true", job)

    def test_benchmark_mismatch_exits_under_optimization_without_reporting_success(self):
        patch = "g['reference'] = lambda *args: types.SimpleNamespace(find_matches=lambda *a, **kw: {})"
        result = fleet.run_main(ROOT, "scripts/benchmark_iupac_masks.py", patch, optimized=True, timeout=30)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("output mismatch", result.stderr)
        self.assertEqual(result.stdout, "")


if __name__ == "__main__":
    unittest.main()
