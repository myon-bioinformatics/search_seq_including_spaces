# IUPAC compatibility with four-bit masks

The anchor and streaming matchers previously allocated Python sets for every
IUPAC-symbol comparison. They now share a 15-entry lookup of possible bases:
A=1, C=2, G=4, T=8, and ambiguity codes contain the union of those bits.
Compatibility is a nonzero bitwise AND. Input validation, literal matching,
reverse complements, coordinates, ordering, limits and output schemas are unchanged.
The existing `_BASES` strings remain available to candidate generation.

This is a four-bit **logical representation** used in Python integer objects.
It does not pack the sequence buffer, shrink serialized reports, or demonstrate
smaller builds. A fixed lookup is retained for the lifetime of the module; no
net process-memory reduction is claimed.

## Verification

- `python -S -m unittest discover -s tests -p test_iupac_bitmasks.py -v`
  checks every one of the 225 symbol pairs in each anchor position, including
  a chunk boundary, against independently declared base sets. Literal behavior
  is checked separately.
- `uv run --locked pytest -xq`: 340 passed, 675 subtests passed
  (CPython 3.14.9, Linux). Existing randomized tests cover both strands,
  source coordinates, variable gaps and chunk boundaries.
- `python -S -m scripts.benchmark_iupac_masks` loads the fixed pre-change
  commit, asserts equality of whole result dictionaries and streaming event
  lists, then alternates the two versions for nine samples each.

## Local timing experiment

Raw samples and environment are in [iupac-mask-measurement.json](iupac-mask-measurement.json).
CPython 3.12.14 on Linux, seed 20261104, 20,000 synthetic concrete DNA bases,
IUPAC anchors ACR/YGT, both strands, gap range 0–2; output construction and
normalization are included. Reference commit: `5ab1885919c52c46a255c7944ef475ae8061c46d`.

| Path | Previous median | Mask median |
| --- | ---: | ---: |
| Whole sequence | 70.946 ms | 55.308 ms |
| 1,024-character chunks | 69.416 ms | 52.450 ms |

These are about 22% and 24% reductions in this local workload. They are not
CI/build-time measurements or guarantees for other inputs. Literal-mode
performance, process peak memory and deployed artifact size were not benchmarked.
No dependencies or workflow changes are added.
