# In-process TestLoader experiment (2026-11-04)

Decision: use the standard library APIs directly; do not introduce a reusable
custom runner. This experiment adds no production runtime code, dependency,
CI pass or pytest conversion. The existing `suite` snapshot remains the safe
default. `reuse` is measured only as an unsafe control and is not recommended.

## Small reusable recipe

Run from the repository root. `discover` controls module loading; the public
`testNamePatterns` property controls method selection (full-ID glob patterns,
not CLI `-k` substring syntax). No directory walker or test-name parser is needed.

```python
import unittest

loader = unittest.TestLoader()
loader.testNamePatterns = ["*.test_only_standard_library_is_imported"]
suite = loader.discover("tests", pattern="test_layers.py")
if loader.errors:
    raise RuntimeError("\n".join(loader.errors))
if suite.countTestCases() == 0:
    raise RuntimeError("No tests selected")
result = unittest.TextTestRunner(verbosity=2).run(suite)
if not result.wasSuccessful():
    raise RuntimeError("Tests failed")
```

For a known module or exact case, put the test start directory on `sys.path`
(or use the project's existing import setup), then use
`loader.loadTestsFromNames(["test_layers.TestLayers.test_only_standard_library_is_imported"])`
instead of `discover`. This avoids directory discovery. No `tests/__init__.py`
is required. Keep `loader.errors` and zero-count checks. Use a fresh loader and
suite on every iteration; running a consumed TestSuite again is unsupported.
Do not add a generic cache around loaded suites or their TestCase instances.

For this repository's experiment, `test_layers.use_cache("suite")` initializes
the source cache once before repeated runs. Module setup/teardown then clears
it for every nonempty suite. This cache is for protein source inspection;
it is separate from TestLoader's search for test modules.

## Scope and correctness

Baseline: `5ab1885919c52c46a255c7944ef475ae8061c46d` (main).
The probe uses only public loader/suite/result APIs; mocks count the existing
source scanner and AST parser without changing their behavior.

- TestLoader discover, loadTestsFromNames, unittest.main(module=test_layers),
  and pytest collect the same eight logical IDs and all pass. Pytest runs the
  classes in source order; unittest sorts them. Direct `python tests/test_layers.py`
  reports `__main__` as its module name. The in-process main comparison passes
  the imported module explicitly; pytest IDs are recorded from item class and
  callable metadata, without parsing node-ID strings.
- All eight tests and the selected one are compared separately. Exact named
  loading and `testNamePatterns` selection both work. A nonexistent explicit
  name yields a loader error/failed test, not an empty valid suite.
- Repeated one-test runs call the real setUpModule and tearDownModule once per
  run. `suite` clears its cache on both hooks; `off` retains no entries.
  `use_cache` always creates a fresh snapshot, including for `reuse`.
- Zero selected tests call neither module hook. A preexisting snapshot is not
  cleared. A bare TestResult reports wasSuccessful() for zero tests; on this runtime,
  both unittest CLI entry points and pytest exit 5 (NO TESTS RAN). Do not
  generalize an older CLI exit-0 observation to this interpreter.
  The recipe explicitly rejects zero rather than treating it as verification.
- An isolated temporary package starts with a valid source, then gains an
  invalid third-party import between suites. `off` and `suite` fail the next
  real layer check. **`reuse` incorrectly passes** because it retains the old
  AST. No unsafe cache policy is adopted. Existing cache tests also cover
  added/removed files and mutation of shared trees.

Reproduce from the repo root (pytest is required only for comparison):

```sh
python -S -m unittest discover -s tests -p test_testloader_probe.py -v
python -m pytest tests/test_layers.py tests/test_layers_cache.py tests/test_testloader_probe.py -q
python experiments/testloader_probe.py > experiments/testloader-results.json
```

The JSON contains individual samples, IDs, results, actual source scan/parse
counts, cache entries and CLI logs. The probe asserts ID-set and result parity;
its benchmark never selects tests with a custom name parser. It remains an
experiment rather than an interface for CI or other repositories.

## Measurements and validation

Linux / CPython 3.12.14; pytest 9.2.0. Same source checkout, one process,
21 alternating-order samples per method/mode/selection (252 suite runs).
Each mode has a separate snapshot; reuse includes its first cold iteration.
Modules are already imported before timing. Load and execution times are
measured separately; totals below are medians in milliseconds. All runs
include scanner/parser counting instrumentation; these are local observations,
not CI speedup or cross-machine claims. No other test command was launched
during the final measurement run.

| Selection | API | Cache | Load ms | Total ms | Source scans | AST parses |
|---|---|---|---:|---:|---:|---:|
| all | discover | off | 0.570 | 121.341 | 7 | 35 |
| all | discover | suite | 0.639 | 86.557 | 1 | 7 |
| all | discover | reuse | 0.614 | 77.505 | 0 | 0 |
| all | names | off | 0.210 | 130.864 | 7 | 35 |
| all | names | suite | 0.189 | 92.931 | 1 | 7 |
| all | names | reuse | 0.207 | 78.905 | 0 | 0 |
| selected | discover | off | 0.519 | 16.503 | 1 | 7 |
| selected | discover | suite | 0.458 | 20.291 | 1 | 7 |
| selected | discover | reuse | 0.516 | 9.913 | 0 | 0 |
| selected | names | off | 0.189 | 17.736 | 1 | 7 |
| selected | names | suite | 0.188 | 18.256 | 1 | 7 |
| selected | names | reuse | 0.171 | 8.872 | 0 | 0 |

`all` is 8 tests; `selected` is the same one import-policy test for every mode.
Discovery overhead is below 1 ms here. Named loading avoids discovery, but
the large reduction for all eight tests comes from existing suite-scoped
source caching (7 scans/35 parses to 1/7), not a new loader optimization.
For the one selected test, off and suite both need 1 scan/7 parses; no
source-work reduction is expected. The unsafe reuse control reduces warm
source work to zero but fails freshness verification, so its timing is not
an acceptable optimization.

The JSON also records one fresh-process CLI sample for all/selected/zero.
These single observations include startup and are not a repeated speed comparison.
All/selected pass with 8/1 tests, while zero selects none and exits 5 in all
three CLIs. API loading, unittest.main and pytest have equal logical ID sets.

Validation: full pytest **346 passed, 225 subtests passed**; stdlib probe
**7 tests passed**. Existing runtime, workflow and dependency files are unchanged.
Added Python: **287 lines** (experiment and regression tests combined);
new reusable production runner: **0 lines**. The portable recipe above is
sufficient for other unittest repositories; no cross-repository rollout or
performance claim is made without repository-specific validation.
