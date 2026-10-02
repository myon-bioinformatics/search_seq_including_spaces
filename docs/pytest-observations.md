# Native pytest evidence

CI writes `reports/pytest-events.jsonl` through the shared xprobe pytest adapter
and preserves it with the JUnit report, including when tests fail. All Python matrix jobs use the adapter. The application runtime remains stdlib-only.

For local runs, choose a fresh destination (existing evidence is never overwritten):

```sh
python -m pytest -p vendor.xprobe_pytest --xprobe-jsonl=reports/local-001.jsonl --xprobe-repository=myon-bioinformatics/search_seq_including_spaces
```

The vendored adapter is byte-for-byte from xprobe commit
`326acd667e13b21bf53ccc1590af960edf8cbf6c`, `scripts/xprobe_pytest.py`.
Its license and SHA-256/Git blob provenance live alongside it in `vendor/`.
Update from a reviewed upstream commit and refresh the provenance together.

Records distinguish setup/call/teardown and collection, with native failure,
error, skip, xfail, xpass and strict XPASS outcomes. Records describe phases,
not one final result per test. A finish record's `complete` means no observations
were dropped; it does not mean tests passed. Missing finish means partial evidence.
The pytest exit status is preserved. Serial pytest only; xdist is not supported.

Repository identity is explicit. Commit SHA stays null until supplied from the
canonical metadata producer; GitHub's environment SHA is not substituted.
Tracebacks, captured output, marker reasons and parameter labels are omitted.
Test names remain visible, and node hashes are not encryption. Reports are
14-day Actions artifacts, not published Pages content.

## Explore downloaded reports

Use the pinned local helper; no separate xprobe checkout is required:

```sh
python -S -c "from pathlib import Path; from vendor.xprobe import corpus_from_json, merge_cases; rows = corpus_from_json(Path('reports/local-001.jsonl').read_text(encoding='utf-8'), jsonl=True); print(len(merge_cases(rows)))"
```

`vendor/xprobe.py` is byte-for-byte from upstream receipt-helper revision
`dbd56e22e91a1e7553053da4f6965758b9994e8e` ([paired xprobe PR #5](https://github.com/myon-bioinformatics/xprobe/pull/5)), with its own
SHA-256/Git blob provenance and the shared `vendor/xprobe-LICENSE`.
Use `vendor.xprobe.corpus_from_json(text, jsonl=True)` to validate observations,
`merge_cases(*corpora)` to deduplicate identical IDs (conflicts raise), and
`corpus_to_json(rows, jsonl=True)` to serialize them. Preserve run IDs when merging;
do not deduplicate distinct executions merely because test names match.
Upstream owns helper unit tests; this repository checks provenance and consumer
integration. Use `vendor.xprobe.pytest_receipt(raw_jsonl, expected_context=...)`
before corpus deduplication to validate one complete run. The shared helper owns
start/schema/finish, ID continuity, identity and phase/count validation. A valid
receipt does not make a failed run pass: compare its preserved exit status with
the producer process. Missing finish, dropped observations, mixed runs and
malformed/count-inconsistent evidence raise ValueError.

The vendor update and removal of local envelope/identity checks are part of this
same SearchSeq PR #24. The upstream receipt tests and adapter integration own
xfail/XPASS/setup/teardown and corrupt-evidence coverage; the local nested smoke
test runs a real reader assertion failure and checks the adapter/helper connection.
FASTA/FASTQ domain regressions remain here. This immutable upstream PR revision
is test-validated, not yet merged; any upstream revision must refresh this pin
and both provenance checks before either PR is considered ready.
Keep downloaded artifacts in separate job directories so identical filenames
do not overwrite each other. Run IDs distinguish observations across jobs.

Reproduction inputs and chat-text references must be explicit, reviewed records
linked by `origin_id`; this adapter does not infer causes or ingest chat history.

## CI gates and shared coverage

CodeQL is deliberately removed from this repository's CI. The behavior gates
are `python -S -m unittest discover -s tests` plus the full pytest suite on
Python 3.10–3.14; a failed test retains its nonzero status. Workflow changes are
checked separately by the existing shared actionlint workflow. Docs-only changes
do not trigger the Tests workflow.

| Coverage / evidence | Source and role |
| --- | --- |
| Reader behavior | `tests/test_sequence_readers.py`: FASTA/FASTQ, malformed input, compression, stream ownership, record-boundary blanks, provenance and 8192-character sniff boundaries |
| Existing domain regression | Anchor matcher, protein parser/features/rendering, CLI golden outputs and comparison tests in `tests/` |
| Cross-repository adapter integrity | `tests/test_native_observations.py::test_adapter_provenance`: checks the vendored xprobe adapter against its pinned upstream commit, blob and SHA-256 |
| Shared receipt validation | xprobe `tests/test_pytest_receipt.py` and `tests/test_pytest_adapter.py`: generic receipt corruption, all native outcomes, truncation/interruption and producer integration |
| Reader evidence integration | `test_native_failure_evidence`: runs a real reader assertion failure through the vendored adapter, calls the shared receipt validator, checks preserved exit status and omitted private diagnostics |
| Shared JUnit identity | `myon-bioinformatics/myon-bioinformatics/.github/workflows/reusable-junit-identity.yml@4dfda95d6573250477f991a0421fa6acb9bc0258`: collects the five exact matrix report paths, runs even after test failure and preserves collection evidence |

The shared collector checks report identity/completeness; pytest checks behavior.
Neither a successfully collected report nor a complete native finish receipt
turns a failed test into a pass. These checks do not claim CodeQL-equivalent
static security analysis. Reuse upstream helpers with provenance and meaningful
consumer regressions, rather than copying unrelated repositories' full suites.
