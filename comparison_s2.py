"""S2 local comparison helpers for already aligned sequence pairs.

Builds on comparison.py's S1 semantics without performing alignment.
Sliding windows use aligned 1-based inclusive coordinates and preserve the
same gap/ambiguity policies as compare_aligned().
"""

from __future__ import annotations

import hashlib

from comparison import (
    ComparisonError,
    _alphabet_tables,
    _candidate_set,
    _percent,
    _prepare_aligned,
)

SCHEMA = "sequence-comparison-windows/v1"


def sliding_identity(
    sequence_a: str,
    sequence_b: str,
    *,
    alphabet: str,
    window: int,
    step: int = 1,
    gap_policy: str = "exclude",
    ambiguity_policy: str = "exclude",
    start: int = 1,
    end: int | None = None,
) -> dict:
    """Return full sliding-window identity/coverage records in O(L).

    Only complete windows are emitted. A window larger than the selected
    region is rejected rather than silently producing a partial window.
    """
    a, b, end = _prepare_aligned(
        sequence_a,
        sequence_b,
        alphabet=alphabet,
        gap_policy=gap_policy,
        ambiguity_policy=ambiguity_policy,
        start=start,
        end=end,
    )
    if isinstance(window, bool) or not isinstance(window, int) or window < 1:
        raise ComparisonError(f"window must be an integer >= 1, got {window!r}")
    if isinstance(step, bool) or not isinstance(step, int) or step < 1:
        raise ComparisonError(f"step must be an integer >= 1, got {step!r}")
    selected = end - start + 1
    if window > selected:
        raise ComparisonError(
            f"window {window} exceeds selected aligned region length {selected}"
        )

    canonical, ambiguous = _alphabet_tables(alphabet)
    columns = [
        _column_metrics(left, right, canonical, ambiguous, gap_policy, ambiguity_policy)
        for left, right in zip(a, b)
    ]
    keys = (
        "matches",
        "mismatches",
        "compared_positions",
        "excluded_positions",
        "gap_columns",
        "gap_mismatches",
        "ambiguous_columns",
        "compatible_matches",
    )
    prefixes = {key: [0] for key in keys}
    for row in columns:
        for key in keys:
            prefixes[key].append(prefixes[key][-1] + row[key])

    windows = []
    first0 = start - 1
    last_start0 = end - window
    for left0 in range(first0, last_start0 + 1, step):
        right0 = left0 + window
        counts = {key: prefixes[key][right0] - prefixes[key][left0] for key in keys}
        compared = counts["compared_positions"]
        identity = _percent(counts["matches"], compared)
        compatible = (
            _percent(counts["compatible_matches"], compared)
            if ambiguity_policy == "compatible"
            else None
        )
        windows.append(
            {
                "start": left0 + 1,
                "end": right0,
                "columns": window,
                "counts": {
                    **counts,
                    "compatible_matches": (
                        counts["compatible_matches"]
                        if ambiguity_policy == "compatible"
                        else None
                    ),
                },
                "percent_identity": identity,
                "divergence_percent": None if identity is None else round(100.0 - identity, 6),
                "percent_compatible": compatible,
                "coverage_percent": _percent(compared, window),
            }
        )

    return {
        "schema": SCHEMA,
        "metric": "sliding_percent_identity",
        "alphabet": alphabet,
        "alignment": "pre-aligned-required",
        "coordinate_convention": "1-based-inclusive",
        "region": {"start": start, "end": end, "columns": selected},
        "window": window,
        "step": step,
        "windows": windows,
        "provenance": {
            "gap_policy": gap_policy,
            "ambiguity_policy": ambiguity_policy,
            "sequence_a_sha256": hashlib.sha256(a.encode("ascii")).hexdigest(),
            "sequence_b_sha256": hashlib.sha256(b.encode("ascii")).hexdigest(),
        },
    }


def render_ascii_diff(
    sequence_a: str,
    sequence_b: str,
    *,
    alphabet: str,
    gap_policy: str = "exclude",
    ambiguity_policy: str = "exclude",
    start: int = 1,
    end: int | None = None,
    width: int = 60,
) -> str:
    """Render a bounded aligned diff without changing comparison semantics.

    Marker legend:
    | exact compared match
    ~ compatible but non-identical (compatible policy only)
    . compared mismatch
    ^ one-sided gap counted as mismatch
    ? excluded ambiguity
      excluded gap-gap / gap column
    """
    a, b, end = _prepare_aligned(
        sequence_a,
        sequence_b,
        alphabet=alphabet,
        gap_policy=gap_policy,
        ambiguity_policy=ambiguity_policy,
        start=start,
        end=end,
    )
    if isinstance(width, bool) or not isinstance(width, int) or width < 10:
        raise ComparisonError(f"width must be an integer >= 10, got {width!r}")
    canonical, ambiguous = _alphabet_tables(alphabet)
    marks = "".join(
        _marker(left, right, canonical, ambiguous, gap_policy, ambiguity_policy)
        for left, right in zip(a[start - 1:end], b[start - 1:end])
    )
    left_seq = a[start - 1:end]
    right_seq = b[start - 1:end]

    blocks: list[str] = []
    for offset in range(0, len(left_seq), width):
        chunk_a = left_seq[offset:offset + width]
        chunk_m = marks[offset:offset + width]
        chunk_b = right_seq[offset:offset + width]
        pos = start + offset
        blocks.append(
            f"{pos:>8}  A {chunk_a}\n"
            f"{'':>8}    {chunk_m}\n"
            f"{pos:>8}  B {chunk_b}"
        )
    legend = (
        "legend: | exact  ~ compatible  . mismatch  ^ gap-mismatch  "
        "? excluded-ambiguity  [space] excluded-gap"
    )
    return legend + "\n\n" + "\n\n".join(blocks) + "\n"


def _column_metrics(
    left: str,
    right: str,
    canonical: frozenset[str],
    ambiguous: dict[str, frozenset[str]],
    gap_policy: str,
    ambiguity_policy: str,
) -> dict[str, int]:
    row = {
        "matches": 0,
        "mismatches": 0,
        "compared_positions": 0,
        "excluded_positions": 0,
        "gap_columns": 0,
        "gap_mismatches": 0,
        "ambiguous_columns": 0,
        "compatible_matches": 0,
    }
    if left == "-" or right == "-":
        row["gap_columns"] = 1
        if left == right or gap_policy == "exclude":
            row["excluded_positions"] = 1
            return row
        row["compared_positions"] = 1
        row["mismatches"] = 1
        row["gap_mismatches"] = 1
        return row

    left_amb = left in ambiguous
    right_amb = right in ambiguous
    if left_amb or right_amb:
        row["ambiguous_columns"] = 1
        if ambiguity_policy == "exclude":
            row["excluded_positions"] = 1
            return row

    row["compared_positions"] = 1
    if left == right:
        row["matches"] = 1
    else:
        row["mismatches"] = 1
    if ambiguity_policy == "compatible":
        if _candidate_set(left, canonical, ambiguous) & _candidate_set(
            right, canonical, ambiguous
        ):
            row["compatible_matches"] = 1
    return row


def _marker(
    left: str,
    right: str,
    canonical: frozenset[str],
    ambiguous: dict[str, frozenset[str]],
    gap_policy: str,
    ambiguity_policy: str,
) -> str:
    if left == "-" or right == "-":
        if left != right and gap_policy == "mismatch":
            return "^"
        return " "
    if (left in ambiguous or right in ambiguous) and ambiguity_policy == "exclude":
        return "?"
    if left == right:
        return "|"
    if ambiguity_policy == "compatible" and (
        _candidate_set(left, canonical, ambiguous)
        & _candidate_set(right, canonical, ambiguous)
    ):
        return "~"
    return "."
