# Discover files and production entry points by convention

## Built-in sequence input discovery

```sh
python -S sequence_tool.py discover ./data --format fasta --format fastq --pattern 'sample*'
python -S sequence_tool.py discover tests --format fasta --format fastq --pattern '*.fa*'
```

This native command uses `sequence_readers.discover_sequences()`. It reuses the
reader's extension and compression rules (`.fa`, `.fasta`, `.fastq`, `.gbk`,
`.jsonl`, etc., with gzip/bz2/xz/lzma suffixes). Repeated `--format` and
`--pattern` select alternatives; the groups intersect. Patterns match basenames
or relative paths, case-sensitively. Results are JSON records with `path`,
`format` and `compression`, inside a `files`/`truncated` report.

Discovery lists filename candidates without reading/decompressing contents;
use the existing reader/analysis commands to validate actual sequence records.
Recursion is on by default (`--no-recursive` disables it); hidden paths are
opt-in via `--include-hidden`, and symlinks are skipped. `--max-files` defaults
to 10,000, with exact-limit detection. Exit 0 means complete (including empty),
2 means a bad input/walk error or a truncated report. The command runs with the
repository alone; the optional xprobe checkout below is for source inspection.

The second command was run against the repository fixture tree and selected
3 FASTA inputs with no truncation.

Validation: the full pytest suite passed **314 tests and 160 subtests** on
Python 3.12.14. New discovery cases cover reader extension/compression rules,
multiple naming rules, exact/global limits, hidden paths, symlinks, walk errors,
and standalone `python -S` execution with deliberately invalid file contents.

## Optional source-symbol inspection

Find sequence readers and validation entry points before adding a format or changing input handling.

Use a sibling xprobe checkout with the convention-discovery API from
[xprobe PR #28](https://github.com/myon-bioinformatics/xprobe/pull/28).
For a fresh checkout, select this tested revision:

```sh
git clone https://github.com/myon-bioinformatics/xprobe.git ../xprobe
git -C ../xprobe checkout --detach 045b54dd35baf43758a484e458712ed5ad47ba24
```

The commands below run from this repository root. If xprobe is elsewhere,
replace `../xprobe/xprobe.py` with that checkout's path. They use Python's
standard library and work with `python -S`.

```sh
python -S ../xprobe/xprobe.py --symbols . --include sequence_readers.py --include structured_readers.py --include 'protein/*.py' --symbol-pattern 'read_*' --symbol-pattern '*validate*' --brief
python -S ../xprobe/xprobe.py --discover --root . --include '*.yml' --include '*.toml' --include-hidden --brief
```

Repeat `--symbol-pattern` or `--include` to inject different naming conventions;
alternatives are ORed. `--symbol-kind class` (also `function`, `async_function`)
filters the selected definitions. For example, replace the name rules with
`--symbol-pattern '*Adapter' --symbol-pattern '*Reader' --symbol-kind class`
to look for class families. There is no required `tests/` path or `test_*` prefix.

Filename patterns match basenames or relative paths; symbol patterns match
short or lexical dotted names. Both are case-sensitive. Wildcard `*` spans
slashes/dots, so a parent-name match also selects nested definitions. Repeated
source definitions remain separate. This is a source inventory, not a count of
callable runtime members or an import/plugin registration mechanism.

JSON is the default for further processing; omit `--brief`. Exit status 0 means
complete matches, 1 means none, and 2 means errors, skipped inputs or truncated
output. Default bounds are 5,000 matched files, 1 MB per source and 1,000 symbols;
set `--max-files`, `--max-bytes`, `--max-matches` when needed. Hidden paths are
opt-in; dependency/build directories are excluded by default. Targets are parsed
without importing them. No `.gitignore` interpretation is implied.

## Observed trial

On 2026-10-21 at `725a3602327b083e60263a1eabc83a82436b84c3`, the first command parsed
9 files and found 6 definitions with no errors, skips or
truncation. Examples:

- `sequence_readers.py`: `read_sequences`
- `structured_readers.py`: `read_genbank`
- `structured_readers.py`: `read_json`
- `structured_readers.py`: `read_json.objects`

The second command selected 6 YAML/TOML files at the same revision.
These are observations for that snapshot, not future expected counts. Full
query options, source SHA-256 values and results are recorded in xprobe's
[`docs/evidence/convention-discovery.json`](https://github.com/myon-bioinformatics/xprobe/blob/045b54dd35baf43758a484e458712ed5ad47ba24/docs/evidence/convention-discovery.json).

This is an explicit developer checkout, separate from existing vendored copies.
Existing vendor locks and CI/deployment configuration are unchanged; promoting
these APIs to an installed consumer requires its normal source-pin review.
