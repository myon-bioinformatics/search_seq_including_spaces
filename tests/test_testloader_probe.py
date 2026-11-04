"""Contracts for the experiment's use of public unittest APIs."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments import testloader_probe as probe


class TestPublicLoader(unittest.TestCase):
    def tearDown(self):
        probe.subject.use_cache('suite')

    def test_discover_names_and_main_have_identical_ids_and_results(self):
        results = []
        for method in ('discover', 'names'):
            result = probe.run(probe.load(method))
            self.assertTrue(result.wasSuccessful(), result.errors)
            self.assertEqual(result.testsRun, 8)
            results.append(result.ids)
        captured = probe.RecordedResult()

        class Runner:
            def run(self, suite):
                suite.run(captured)
                return captured

        unittest.main(module=probe.subject, argv=['probe'], testRunner=Runner(), exit=False)
        self.assertEqual(results[0], results[1])
        self.assertEqual(results[0], captured.ids)
        self.assertTrue(captured.wasSuccessful())

    def test_selected_repeated_runs_call_real_module_hooks_and_clear_cache(self):
        for mode in ('off', 'suite', 'reuse'):
            snapshot = probe.subject.use_cache(mode)
            with patch.object(probe.subject, 'setUpModule', wraps=probe.subject.setUpModule) as setup, \
                 patch.object(probe.subject, 'tearDownModule', wraps=probe.subject.tearDownModule) as teardown:
                for _ in range(2):
                    suite = probe.load(patterns=[
                        'test_layers.TestLayers.test_only_standard_library_is_imported'])
                    self.assertEqual(suite.countTestCases(), 1)
                    result = probe.run(suite)
                    self.assertTrue(result.wasSuccessful())
                    self.assertEqual(result.testsRun, 1)
                    if mode != 'reuse':
                        self.assertEqual(snapshot.cached_entries(), 0)
                self.assertEqual(setup.call_count, 2)
                self.assertEqual(teardown.call_count, 2)
            if mode == 'reuse':
                self.assertGreater(snapshot.cached_entries(), 0)
            self.assertIsNot(probe.subject.use_cache(mode), snapshot)
            self.assertEqual(probe.subject.SNAPSHOT.cached_entries(), 0)

    def test_zero_selection_does_not_call_module_hooks_or_clear_snapshot(self):
        snapshot = probe.subject.use_cache('suite')
        snapshot.modules(probe.subject.PACKAGE)
        for method in ('discover', 'names'):
            with patch.object(probe.subject, 'setUpModule') as setup, \
                 patch.object(probe.subject, 'tearDownModule') as teardown:
                suite = probe.load(method, ['*.does_not_exist'])
                self.assertEqual(suite.countTestCases(), 0)
                result = probe.run(suite)
                self.assertEqual(result.testsRun, 0)
                self.assertTrue(result.wasSuccessful())  # callers must reject zero explicitly
                setup.assert_not_called()
                teardown.assert_not_called()
                self.assertGreater(snapshot.cached_entries(), 0)

    def test_explicit_name_uses_public_api_without_discovery(self):
        loader = unittest.TestLoader()
        name = 'TestLayers.test_only_standard_library_is_imported'
        suite = loader.loadTestsFromNames([name], module=probe.subject)
        self.assertFalse(loader.errors)
        result = probe.run(suite)
        self.assertTrue(result.wasSuccessful())
        self.assertEqual(result.ids, ['test_layers.' + name])

    def test_invalid_explicit_name_is_loader_error_not_zero_selection(self):
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromName('TestLayers.missing', module=probe.subject)
        self.assertTrue(loader.errors)
        result = probe.run(suite)
        self.assertEqual(result.testsRun, 1)
        self.assertEqual(len(result.errors), 1)

    def test_reuse_masks_real_source_regression_and_is_not_safe_to_adopt(self):
        for mode in ('off', 'suite', 'reuse'):
            with tempfile.TemporaryDirectory() as directory:
                package = Path(directory)
                source = package / 'residues.py'
                source.write_text('value = 1\n', encoding='utf-8')
                with patch.object(probe.subject, 'PACKAGE', package):
                    probe.subject.use_cache(mode)
                    pattern = ['*.test_only_standard_library_is_imported']
                    first = probe.run(probe.load(patterns=pattern))
                    self.assertTrue(first.wasSuccessful())
                    source.write_text('import imaginary_third_party\n', encoding='utf-8')
                    second = probe.run(probe.load(patterns=pattern))
                    self.assertEqual(second.wasSuccessful(), mode == 'reuse')
                    self.assertEqual(len(second.failures), 0 if mode == 'reuse' else 1)

    def test_measurement_reports_actual_scans_including_off(self):
        rows = []
        for mode in ('off', 'suite', 'reuse'):
            probe.subject.use_cache(mode)
            first = probe.measurement('discover', mode)
            second = probe.measurement('discover', mode)
            self.assertTrue(first['success'] and second['success'])
            rows.append((first, second))
        self.assertGreater(rows[0][0]['source_scans'], rows[1][0]['source_scans'])
        self.assertGreater(rows[0][0]['ast_parses'], rows[1][0]['ast_parses'])
        self.assertEqual(rows[1][1]['source_scans'], 1)
        self.assertEqual(rows[2][1]['source_scans'], 0)
        self.assertEqual(rows[2][1]['ast_parses'], 0)


if __name__ == '__main__':
    unittest.main()
