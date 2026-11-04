"""Contracts for the experiment's use of public unittest APIs."""
from pathlib import Path
import sys
import subprocess
import copy
import tempfile
import unittest
from unittest.mock import patch

with patch.object(sys, 'path', [str(Path(__file__).resolve().parents[1]), *sys.path]):
    from experiments import testloader_probe as probe


class TestPublicLoader(unittest.TestCase):
    def tearDown(self):
        probe.subject.use_cache('suite')

    def test_imports_and_loading_preserve_sys_path(self):
        script = """
import runpy, sys
from pathlib import Path
before = list(sys.path)
for name in ('experiments/testloader_probe.py', 'tests/test_testloader_probe.py'):
    runpy.run_path(name, run_name='import_probe')
    if sys.path != before:
        raise RuntimeError('import changed sys.path: ' + name)
from experiments import testloader_probe as probe
for method in ('discover', 'names'):
    probe.load(method)
    if sys.path != before:
        raise RuntimeError('loading changed sys.path: ' + method)
"""
        completed = subprocess.run([sys.executable, '-S', '-c', script], cwd=probe.ROOT,
                                   capture_output=True, text=True, timeout=30)
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_commands_have_timeout_and_propagate_timeout(self):
        with patch.object(probe.subprocess, 'run') as run_command:
            probe.run_command(['python', '--version'], check=True)
            self.assertEqual(run_command.call_args.kwargs['timeout'], 120)
            self.assertTrue(run_command.call_args.kwargs['check'])
            run_command.side_effect = subprocess.TimeoutExpired('python', 120)
            with self.assertRaises(subprocess.TimeoutExpired):
                probe.run_command(['python'])
        with patch.object(probe, 'run_command', side_effect=
                          subprocess.TimeoutExpired('python', 120)):
            with self.assertRaises(subprocess.TimeoutExpired):
                probe.entrypoints()

    def test_validation_rejects_corrupted_results(self):
        ids = ['case.' + str(i) for i in range(8)]
        external = {'main_ids': ids, 'commands': [{'selection': 'all', 'returncode': 0}]}
        recorded = {'ids': ids, 'exitcode': 0,
                    'reports': [{'outcome': 'passed'} for _ in ids]}
        rows = [{'selection': 'all', 'ids': ids, 'collected': 8,
                 'success': True, 'errors': 0, 'failures': 0}]
        probe.validate_results(rows, external, recorded)
        for target, key, value in [('external', 'main_ids', []),
                                   ('recorded', 'ids', []),
                                   ('recorded', 'exitcode', 1),
                                   ('recorded', 'reports', [{'outcome': 'failed'}] * 8),
                                   ('row', 'collected', 0), ('row', 'success', False)]:
            with self.subTest(target=target, key=key):
                copies = copy.deepcopy({'external': external, 'recorded': recorded,
                                        'row': rows[0]})
                copies[target][key] = value
                with self.assertRaises(RuntimeError):
                    probe.validate_results([copies['row']], copies['external'], copies['recorded'])

    def test_validation_still_fails_under_optimized_python(self):
        script = ("from experiments.testloader_probe import validate_results; "
                  "validate_results([], {'main_ids': []}, {})")
        completed = subprocess.run([sys.executable, '-O', '-S', '-c', script],
                                   cwd=probe.ROOT, capture_output=True, text=True, timeout=30)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn('Expected eight distinct baseline test IDs', completed.stderr)

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
