"""Reproducible experiment, not a replacement test runner. Run from repo root."""
import contextlib
import io
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
with patch.object(sys, 'path', [str(ROOT / 'tests'), *sys.path]):
    import test_layers as subject

COMMAND_TIMEOUT_SECONDS = 120


def load(method='discover', patterns=None):
    """Selection and traversal belong entirely to unittest's public loader."""
    loader = unittest.TestLoader()
    loader.testNamePatterns = patterns
    # discover may add its start directory; do not leak it to other tests.
    with patch.object(sys, 'path', list(sys.path)):
        if method == 'discover':
            suite = loader.discover(str(ROOT / 'tests'), pattern='test_layers.py')
        elif method == 'names':
            suite = loader.loadTestsFromNames(['test_layers'], module=None)
        else:
            raise ValueError(method)
    if loader.errors:
        raise RuntimeError('\n'.join(loader.errors))
    return suite


class RecordedResult(unittest.TestResult):
    def __init__(self):
        super().__init__()
        self.ids = []

    def startTest(self, test):
        self.ids.append(test.id())
        super().startTest(test)


def run(suite):
    result = RecordedResult()
    suite.run(result)
    return result


def measurement(method, mode, patterns=None):
    """Keep one snapshot across repeated calls; rebuild each consumed suite."""
    with patch.object(subject, '_scan_modules', wraps=subject._scan_modules) as scan, \
         patch.object(subject, '_parse_file', wraps=subject._parse_file) as parse:
        start = time.perf_counter_ns()
        suite = load(method, patterns)
        collected = suite.countTestCases()
        loaded = time.perf_counter_ns()
        result = run(suite)
        end = time.perf_counter_ns()
    return dict(method=method, mode=mode, collected=collected, ids=result.ids,
                success=result.wasSuccessful(), errors=len(result.errors),
                failures=len(result.failures), skipped=len(result.skipped),
                load_ms=(loaded-start)/1e6, run_ms=(end-loaded)/1e6,
                total_ms=(end-start)/1e6, source_scans=scan.call_count,
                ast_parses=parse.call_count,
                cache_entries=subject.SNAPSHOT.cached_entries())


def benchmark(repeats=21):
    rows = []
    # Alternate method and mode order to reduce systematic order effects.
    for label, patterns in [('all', None), ('selected', [
            'test_layers.TestLayers.test_only_standard_library_is_imported'])]:
        snapshots = {(method, mode): subject.SourceSnapshot(mode)
                     for method in ('discover', 'names') for mode in subject.CACHE_MODES}
        combinations = list(snapshots)
        for index in range(repeats):
            order = combinations if index % 2 == 0 else list(reversed(combinations))
            for method, mode in order:
                subject.SNAPSHOT = snapshots[method, mode]
                row = measurement(method, mode, patterns)
                row.update(selection=label, iteration=index)
                rows.append(row)
    subject.use_cache('suite')
    return rows


def run_command(command, *, check=False):
    return subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                          check=check, timeout=COMMAND_TIMEOUT_SECONDS)


def entrypoints():
    """Compare actual CLI entry points, plus their programmatic ID records."""
    commands = [
        [sys.executable, '-S', '-m', 'unittest', 'discover', '-s', 'tests',
         '-p', 'test_layers.py', '-v'],
        [sys.executable, '-S', 'tests/test_layers.py', '-v'],
        [sys.executable, '-m', 'pytest', 'tests/test_layers.py', '-v'],
    ]
    commands = [(selection, command + (['-k', pattern] if pattern else []))
                for selection, pattern in [('all', None),
                    ('selected', 'test_only_standard_library_is_imported'),
                    ('zero', 'does_not_exist')] for command in commands]
    outputs = []
    for selection, command in commands:
        start = time.perf_counter_ns()
        completed = run_command(command)
        outputs.append(dict(selection=selection, command=command, returncode=completed.returncode,
                            wall_ms=(time.perf_counter_ns()-start)/1e6,
                            stdout=completed.stdout, stderr=completed.stderr))
    subject.use_cache('suite')
    main_result = RecordedResult()

    class Runner:
        def run(self, suite):
            suite.run(main_result)
            return main_result

    program = unittest.main(module=subject, argv=['probe'], testRunner=Runner(), exit=False)
    return dict(commands=outputs, main_ids=program.result.ids)


def pytest_record():
    """Separate process: pytest's module fixtures must not overlap ours."""
    import pytest

    class Plugin:
        def __init__(self):
            self.ids = []
            self.reports = []

        def pytest_collection_modifyitems(self, items):
            self.ids = [f'{item.cls.__module__}.{item.cls.__qualname__}.{item.obj.__name__}'
                        for item in items]

        def pytest_runtest_logreport(self, report):
            if report.when == 'call' or report.failed:
                self.reports.append(dict(nodeid=report.nodeid, outcome=report.outcome,
                                         when=report.when))

    plugin = Plugin()
    with contextlib.redirect_stdout(io.StringIO()):
        code = pytest.main(['tests/test_layers.py', '-q'], plugins=[plugin])
    print(json.dumps(dict(ids=plugin.ids, reports=plugin.reports, exitcode=int(code))))


def validate_results(rows, external, recorded):
    """Fail closed even with Python optimization enabled."""
    expected = external['main_ids']
    if len(expected) != 8 or len(set(expected)) != 8:
        raise RuntimeError('Expected eight distinct baseline test IDs')
    if sorted(recorded['ids']) != sorted(expected):
        raise RuntimeError('pytest test IDs differ from unittest')
    if recorded['exitcode'] != 0 or len(recorded['reports']) != 8:
        raise RuntimeError('pytest did not complete eight tests successfully')
    if not all(r['outcome'] == 'passed' for r in recorded['reports']):
        raise RuntimeError('pytest reported a non-passing test')
    if not all(command['returncode'] == (5 if command['selection'] == 'zero' else 0)
               for command in external['commands']):
        raise RuntimeError('CLI outcome differs from the expected selection result')
    for row in rows:
        wanted = expected if row['selection'] == 'all' else [
            'test_layers.TestLayers.test_only_standard_library_is_imported']
        if row['ids'] != wanted or row['collected'] != len(wanted):
            raise RuntimeError('Measured test selection differs from expected IDs')
        if not row['success'] or row['errors'] or row['failures']:
            raise RuntimeError('Measured suite did not pass')


if __name__ == '__main__':
    if '--pytest-record' in sys.argv:
        pytest_record()
    else:
        rows = benchmark()
        external = entrypoints()
        completed = run_command([sys.executable, str(Path(__file__).resolve()),
                                 '--pytest-record'], check=True)
        recorded = json.loads(completed.stdout)
        validate_results(rows, external, recorded)
        output = dict(python=sys.version, platform=platform.platform(),
                      repeats=21, rows=rows, entrypoints=external,
                      pytest=recorded)
        print(json.dumps(output, indent=2))
