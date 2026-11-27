"""Fleet contracts through the parent's shared checks (vendor/fleet_check.py): vendor lock and
xprobe provenance, parent pins, action pins, CI Python setup, benchmark guard; stdlib unittest, also under pytest.

Helper behavior tests belong to upstream xprobe; this repository checks provenance and consumer integration.
"""
import unittest

from vendor import fleet_check as fleet

ROOT = fleet.repo_root(__file__)
PARENT_WORKFLOW_PIN = "27c46411c0a07d9b4d0f1ad5751ed1bb9e869b5c"
PARENT_PAGES_PIN = "27c46411c0a07d9b4d0f1ad5751ed1bb9e869b5c"
# Node sources and the shared checks stay on accepted upstream commits with exact bytes.
ACCEPTED = {
    "vendor/xprobe.py": "1ffbb7b984075c683690bbe4e9deb06199cdfe94",
    "vendor/xprobe_pytest.py": "326acd667e13b21bf53ccc1590af960edf8cbf6c",
    "vendor/fleet_check.py": "20de01921e50e9339aa9271b0547a5be7649b3b3",
}


class FleetTests(unittest.TestCase):
    def test_vendor_lock_matches_checked_in_bytes_accepted_commits_and_sources(self):
        # check_lock verifies ref == commit, the placed bytes against blob_sha / sha256, and each row's origin.
        rows = fleet.check_lock(ROOT, commits=ACCEPTED, sources={
            "vendor/xprobe.py": ("myon-bioinformatics/xprobe", "xprobe.py"),
            "vendor/xprobe_pytest.py": ("myon-bioinformatics/xprobe", "scripts/xprobe_pytest.py"),
            "vendor/xprobe-LICENSE": ("myon-bioinformatics/xprobe", "LICENSE"),
            "vendor/fleet_check.py": (fleet.PARENT, "fleet_check.py"),
            "vendor/myon-bioinformatics-LICENSE": (fleet.PARENT, "LICENSE"),
        })
        adapter = rows["vendor/xprobe_pytest.py"]
        self.assertEqual((adapter["sha256"], adapter["blob_sha"]),
                         ("c0ea71f9d971bf47a9f5f1794ef8a7f641dd17694ea8e7f29795d326c5229285",
                          "70fac53151e97f5f28d23f9961d3837b7a885ea8"))

    def test_pinned_runtime_helper(self):
        from vendor import xprobe
        self.assertTrue(callable(xprobe.pytest_receipt))

    def test_parent_workflow_references_use_reviewed_rollout_commit(self):
        pages = fleet.PARENT + "/.github/workflows/reusable-pages-static.yml"
        self.assertIn((pages, PARENT_PAGES_PIN, None), fleet.check_actions(ROOT, workflow="pages.yml"))
        pins = fleet.parent_refs(ROOT)
        pins.remove(PARENT_PAGES_PIN)
        self.assertEqual(set(pins), {PARENT_WORKFLOW_PIN})

    def test_every_action_is_pinned_by_full_sha_with_version_comment(self):
        fleet.check_actions(ROOT)

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


# 🧰 Scripts: 📝 each one describes itself and ⌨️ one that reads sys.argv declares its arguments (argparse /
# cli_args). Read with ast, never run; the portfolio's script catalog reads the same facts.

class ScriptsDescribeThemselves(unittest.TestCase):
    def test_every_script_has_a_description_and_declared_arguments(self):
        fleet.check_scripts(fleet.repo_root(__file__))

    def test_pull_request_ci_stays_quick(self):
        # ⏱ push limited to branches, superseded PR runs cancelled, every job timed (fleet_check.check_ci)
        fleet.check_ci(fleet.repo_root(__file__))

    def test_root_scripts_import_only_the_standard_library(self):
        root = fleet.repo_root(__file__)
        for path in sorted(root.glob("*.py")):
            if path.name != "conftest.py":
                fleet.check_imports(root, path.name)
