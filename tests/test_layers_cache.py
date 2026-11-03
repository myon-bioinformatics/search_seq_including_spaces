"""Keep the layer-test source snapshot bounded to one module-suite invocation."""
import ast
import importlib.util
from pathlib import Path
import tempfile
import unittest


class TestLayerSnapshot(unittest.TestCase):
    def setUp(self):
        spec = importlib.util.spec_from_file_location(
            "isolated_layer_checks", Path(__file__).with_name("test_layers.py"))
        self.subject = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.subject)
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.subject.PACKAGE = Path(self.directory.name)
        self.subject.use_cache("suite")
        self.subject.setUpModule()
        self.addCleanup(self.subject.tearDownModule)

    def test_reuses_discovery_and_parsing_within_suite(self):
        path = self.subject.PACKAGE / "sample.py"
        path.write_text("value = 1\n", encoding="utf-8")
        self.assertEqual(self.subject.package_modules(), ("sample",))
        self.assertIs(self.subject.package_modules(), self.subject.package_modules())
        first = self.subject.parse_module("sample")
        self.assertIs(first, self.subject.parse_module("sample"))
        self.assertEqual(ast.dump(first), ast.dump(ast.parse("value = 1\n")))
        self.assertEqual(self.subject.SNAPSHOT.misses, 2)

    def test_next_suite_sees_modified_added_and_removed_sources(self):
        root = self.subject.PACKAGE
        (root / "old.py").write_text("old = 1\n", encoding="utf-8")
        (root / "sample.py").write_text("value = 1\n", encoding="utf-8")
        self.assertEqual(self.subject.package_modules(), ("old", "sample"))
        self.subject.parse_module("sample")
        self.subject.tearDownModule()
        self.assertEqual(self.subject.SNAPSHOT.cached_entries(), 0)
        (root / "old.py").unlink()
        (root / "new.py").write_text("new = 1\n", encoding="utf-8")
        (root / "sample.py").write_text("value = 2\n", encoding="utf-8")
        self.subject.setUpModule()
        self.assertEqual(self.subject.package_modules(), ("new", "sample"))
        self.assertEqual(ast.dump(self.subject.parse_module("sample")),
                         ast.dump(ast.parse("value = 2\n")))

    def test_setup_clears_cache_even_without_previous_teardown(self):
        root = self.subject.PACKAGE
        (root / "sample.py").write_text("value = 1\n", encoding="utf-8")
        self.subject.package_modules()
        self.subject.parse_module("sample")
        (root / "sample.py").write_text("value = 2\n", encoding="utf-8")
        self.subject.setUpModule()
        self.assertEqual(self.subject.SNAPSHOT.cached_entries(), 0)
        self.assertEqual(ast.dump(self.subject.parse_module("sample")),
                         ast.dump(ast.parse("value = 2\n")))


class TestCacheModes(unittest.TestCase):
    def setUp(self):
        spec = importlib.util.spec_from_file_location(
            "isolated_layer_modes", Path(__file__).with_name("test_layers.py"))
        self.subject = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.subject)
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.subject.PACKAGE = Path(self.directory.name)
        (self.subject.PACKAGE / "sample.py").write_text("value = 1\n", encoding="utf-8")

    def _edit_between_suites(self, mode):
        snap = self.subject.use_cache(mode)
        self.subject.setUpModule()
        first = self.subject.parse_module("sample")
        self.subject.tearDownModule()
        (self.subject.PACKAGE / "sample.py").write_text("value = 2\n", encoding="utf-8")
        self.subject.setUpModule()
        return snap, first, self.subject.parse_module("sample")

    def test_default_is_suite(self):
        self.assertEqual(self.subject.SNAPSHOT.cache, "suite")

    def test_off_reads_fresh_objects_every_call(self):
        snap, first, second = self._edit_between_suites("off")
        self.assertIsNot(self.subject.parse_module("sample"), second)
        self.assertEqual(snap.cached_entries(), 0)
        self.assertEqual(ast.dump(second), ast.dump(ast.parse("value = 2\n")))

    def test_reuse_keeps_stale_tree_across_suites(self):
        _, first, second = self._edit_between_suites("reuse")
        self.assertIs(first, second)  # documented risk: stale and shared

    def test_key_includes_package(self):
        self.subject.use_cache("reuse")
        self.assertEqual(self.subject.package_modules(), ("sample",))
        other = tempfile.TemporaryDirectory()
        self.addCleanup(other.cleanup)
        self.subject.PACKAGE = Path(other.name)
        self.assertEqual(self.subject.package_modules(), ())

    def test_rejects_unknown_mode(self):
        with self.assertRaises(ValueError):
            self.subject.use_cache("forever")
