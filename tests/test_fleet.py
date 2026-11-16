"""Fleet contracts through the parent's shared checks (vendor/fleet_check.py): vendor lock and
xprobe provenance, parent pins, CI Python setup, benchmark guard; stdlib unittest, also under pytest.

Helper behavior tests belong to upstream xprobe; this repository checks provenance and consumer integration.
"""
import unittest

from vendor import fleet_check as fleet

ROOT = fleet.repo_root(__file__)
PARENT_WORKFLOW_PIN = "9e3da7ae1846e9604ebc24bb5178f7de094a4e77"
# Node sources and the shared checks stay on accepted upstream commits with exact bytes.
ACCEPTED = {
    "vendor/xprobe.py": "1ffbb7b984075c683690bbe4e9deb06199cdfe94",
    "vendor/xprobe_pytest.py": "326acd667e13b21bf53ccc1590af960edf8cbf6c",
    "vendor/fleet_check.py": "81f6aea6231d4490e581443fdeb43a09e6018ff0",
}


class FleetTests(unittest.TestCase):
    def test_vendor_lock_matches_checked_in_bytes_accepted_commits_and_sources(self):
        # check_lock verifies ref == commit, the placed bytes against blob_sha / sha256, and each row's origin.
        fleet.check_lock(ROOT, commits=ACCEPTED, sources={
            "vendor/xprobe.py": ("myon-bioinformatics/xprobe", "xprobe.py"),
            "vendor/xprobe_pytest.py": ("myon-bioinformatics/xprobe", "scripts/xprobe_pytest.py"),
            "vendor/xprobe-LICENSE": ("myon-bioinformatics/xprobe", "LICENSE"),
            "vendor/fleet_check.py": (fleet.PARENT, "fleet_check.py"),
            "vendor/myon-bioinformatics-LICENSE": (fleet.PARENT, "LICENSE"),
        })

    def test_pinned_runtime_helper(self):
        from vendor import xprobe
        self.assertTrue(callable(xprobe.pytest_receipt))

    def test_parent_workflow_references_use_reviewed_rollout_commit(self):
        self.assertEqual(set(fleet.parent_refs(ROOT)), {PARENT_WORKFLOW_PIN})

    def test_setup_precedes_first_plain_python_command(self):
        fleet.workflow_job(ROOT, "test", setup="astral-sh/setup-uv@",
                           inputs=('python-version: "${{ matrix.python-version }}"', "activate-environment: true"))

    def test_benchmark_mismatch_exits_under_optimization_without_reporting_success(self):
        patch = "g['reference'] = lambda *args: types.SimpleNamespace(find_matches=lambda *a, **kw: {})"
        result = fleet.run_main(ROOT, "scripts/benchmark_iupac_masks.py", patch, optimized=True, timeout=30)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("output mismatch", result.stderr)
        self.assertEqual(result.stdout, "")


if __name__ == "__main__":
    unittest.main()
