"""Fleet contracts through the parent's shared checks (vendor/fleet_check.py): vendor lock and
xprobe provenance, parent pins, action pins, CI Python setup, benchmark guard; stdlib unittest, also under pytest.

Helper behavior tests belong to upstream xprobe; this repository checks provenance and consumer integration.
"""
import importlib.util
import sys
from pathlib import Path
import unittest


def load_vendored(name, relpath):
    """Load a vendored file by path under a repo-unique name (cli_args README: Loading a vendored copy)."""
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parents[1] / relpath)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


fleet = load_vendored("search_seq_including_spaces_fleet_check", "vendor/fleet_check.py")

ROOT = fleet.repo_root(__file__)
PARENT_WORKFLOW_PIN = "main"  # 🔁 own sources follow main (2026-11-28)
PARENT_PAGES_PIN = "main"  # 🔁 own sources follow main (2026-11-28)
# Node sources and the shared checks stay on accepted upstream commits with exact bytes.
ACCEPTED = {
    "vendor/xprobe.py",
    "vendor/xprobe_pytest.py",
    "vendor/fleet_check.py",
    "vendor/cli_args.py",
    "vendor/cli_args-LICENSE",
}  # 🔁 locked destinations; they follow main


class FleetTests(unittest.TestCase):
    def test_vendor_lock_matches_checked_in_bytes_accepted_commits_and_sources(self):
        # check_lock verifies ref == commit, the placed bytes against blob_sha / sha256, and each row's origin.
        rows = fleet.check_lock(ROOT, sources={
            "vendor/xprobe.py": ("myon-bioinformatics/xprobe", "xprobe.py"),
            "vendor/xprobe_pytest.py": ("myon-bioinformatics/xprobe", "scripts/xprobe_pytest.py"),
            "vendor/xprobe-LICENSE": ("myon-bioinformatics/xprobe", "LICENSE"),
            "vendor/fleet_check.py": (fleet.PARENT, "fleet_check.py"),
            "vendor/myon-bioinformatics-LICENSE": (fleet.PARENT, "LICENSE"),
            "vendor/cli_args.py": ("myon-bioinformatics/cli_args", "cli_args.py"),
            "vendor/cli_args-LICENSE": ("myon-bioinformatics/cli_args", "LICENSE"),
        })
        self.assertEqual(rows["vendor/xprobe_pytest.py"]["ref"], "main")  # 🔁 bytes are checked by check_lock

    def test_pinned_runtime_helper(self):
        xprobe = load_vendored("search_seq_including_spaces_xprobe", "vendor/xprobe.py")
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

    def test_no_top_level_vendor_package_imports(self):
        assert fleet.check_no_vendor_package_imports(ROOT) == []
