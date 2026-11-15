"""Aligned pairwise sequence identity/coverage, sliding windows and ASCII diff (stdlib only).

This module does not align sequences. Inputs must already have equal aligned
lengths after whitespace removal. Coordinates are 1-based inclusive.

Identity is literal symbol equality among compared columns. Ambiguity handling
is explicit:
- exclude: columns containing an ambiguous symbol are excluded.
- strict: ambiguous symbols are compared literally (e.g. N == N).
- compatible: columns are compared literally for percent_identity and also
  evaluated by candidate-set overlap for percent_compatible.

Gap handling is explicit:
- exclude: any column containing '-' is excluded.
- mismatch: one-sided gaps are mismatches; gap-gap columns remain excluded.

No substitution matrix or evolutionary similarity score is used here.

S1 is compare_aligned(). S2 (sliding_identity(), render_ascii_diff()) builds on
the same aligned-input contract and gap/ambiguity policies; sliding windows use
aligned 1-based inclusive coordinates.
"""

from __future__ import annotations

import hashlib

SCHEMA = "sequence-comparison/v1"
ALPHABETS = ("dna", "rna", "protein")
GAP_POLICIES = ("exclude", "mismatch")
AMBIGUITY_POLICIES = ("exclude", "strict", "compatible")

_DNA_CANONICAL = frozenset("ACGT")
_RNA_CANONICAL = frozenset("ACGU")
_PROTEIN_CANONICAL = frozenset("ACDEFGHIKLMNPQRSTVWYOU")

_DNA_AMBIGUOUS = {
    "R": frozenset("AG"), "Y": frozenset("CT"), "S": frozenset("GC"),
    "W": frozenset("AT"), "K": frozenset("GT"), "M": frozenset("AC"),
    "B": frozenset("CGT"), "D": frozenset("AGT"), "H": frozenset("ACT"),
    "V": frozenset("ACG"), "N": _DNA_CANONICAL,
}
_RNA_AMBIGUOUS = {
    "R": frozenset("AG"), "Y": frozenset("CU"), "S": frozenset("GC"),
    "W": frozenset("AU"), "K": frozenset("GU"), "M": frozenset("AC"),
    "B": frozenset("CGU"), "D": frozenset("AGU"), "H": frozenset("ACU"),
    "V": frozenset("ACG"), "N": _RNA_CANONICAL,
}
_PROTEIN_AMBIGUOUS = {
    "B": frozenset("DN"), "Z": frozenset("EQ"), "J": frozenset("IL"),
    "X": frozenset("ACDEFGHIKLMNPQRSTVWY"),
}


class ComparisonError(ValueError):
    """Invalid or incompatible aligned-comparison input."""


def compare_aligned(
    sequence_a: str,
    sequence_b: str,
    *,
    alphabet: str,
    gap_policy: str = "exclude",
    ambiguity_policy: str = "exclude",
    start: int = 1,
    end: int | None = None,
) -> dict:
    """Compare two pre-aligned sequences and return a JSON-friendly record.

    Whitespace is removed and ASCII letters are upper-cased. The aligned
    lengths must then be equal. start/end use 1-based inclusive aligned
    coordinates. end=None means the final aligned column.
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

    canonical, ambiguous = _alphabet_tables(alphabet)
    matches = mismatches = compatible_matches = 0
    compared = excluded = gap_columns = ambiguous_columns = 0
    gap_mismatches = 0

    for left, right in zip(a[start - 1:end], b[start - 1:end]):
        if left == "-" or right == "-":
            gap_columns += 1
            if left == right:
                excluded += 1
                continue
            if gap_policy == "exclude":
                excluded += 1
                continue
            compared += 1
            mismatches += 1
            gap_mismatches += 1
            continue

        left_amb = left in ambiguous
        right_amb = right in ambiguous
        if left_amb or right_amb:
            ambiguous_columns += 1
            if ambiguity_policy == "exclude":
                excluded += 1
                continue

        compared += 1
        if left == right:
            matches += 1
        else:
            mismatches += 1

        if ambiguity_policy == "compatible":
            if _candidate_set(left, canonical, ambiguous) & _candidate_set(
                right, canonical, ambiguous
            ):
                compatible_matches += 1

    region_columns = end - start + 1
    identity = _percent(matches, compared)
    compatible = (
        _percent(compatible_matches, compared)
        if ambiguity_policy == "compatible"
        else None
    )
    coverage = _percent(compared, region_columns)

    return {
        "schema": SCHEMA,
        "metric": "percent_identity",
        "alphabet": alphabet,
        "alignment": "pre-aligned-required",
        "coordinate_convention": "1-based-inclusive",
        "region": {"start": start, "end": end, "columns": region_columns},
        "sequence_a_length": len(a),
        "sequence_b_length": len(b),
        "counts": {
            "matches": matches,
            "mismatches": mismatches,
            "compared_positions": compared,
            "excluded_positions": excluded,
            "gap_columns": gap_columns,
            "gap_mismatches": gap_mismatches,
            "ambiguous_columns": ambiguous_columns,
            "compatible_matches": (
                compatible_matches if ambiguity_policy == "compatible" else None
            ),
        },
        "percent_identity": identity,
        "percent_compatible": compatible,
        "coverage_percent": coverage,
        "provenance": {
            "gap_policy": gap_policy,
            "ambiguity_policy": ambiguity_policy,
            "sequence_a_sha256": hashlib.sha256(a.encode("ascii")).hexdigest(),
            "sequence_b_sha256": hashlib.sha256(b.encode("ascii")).hexdigest(),
        },
    }


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
        "schema": "sequence-comparison-windows/v1",
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


def _prepare_aligned(
    sequence_a: str,
    sequence_b: str,
    *,
    alphabet: str,
    gap_policy: str,
    ambiguity_policy: str,
    start: int,
    end: int | None,
) -> tuple[str, str, int]:
    """Normalize and validate the shared S1/S2 aligned-input contract."""
    if alphabet not in ALPHABETS:
        raise ComparisonError(f"alphabet must be one of {', '.join(ALPHABETS)}")
    if gap_policy not in GAP_POLICIES:
        raise ComparisonError(f"gap_policy must be one of {', '.join(GAP_POLICIES)}")
    if ambiguity_policy not in AMBIGUITY_POLICIES:
        raise ComparisonError(
            f"ambiguity_policy must be one of {', '.join(AMBIGUITY_POLICIES)}"
        )

    a = _normalize(sequence_a)
    b = _normalize(sequence_b)
    if len(a) != len(b):
        raise ComparisonError(
            f"aligned sequences must have equal length ({len(a)} != {len(b)})"
        )
    if not a:
        raise ComparisonError("aligned sequences must not be empty")

    _validate(a, alphabet, "sequence_a")
    _validate(b, alphabet, "sequence_b")

    if isinstance(start, bool) or not isinstance(start, int) or start < 1:
        raise ComparisonError(f"start must be an integer >= 1, got {start!r}")
    if end is None:
        end = len(a)
    if isinstance(end, bool) or not isinstance(end, int) or end < start or end > len(a):
        raise ComparisonError(
            f"end must be an integer in {start}..{len(a)}, got {end!r}"
        )
    return a, b, end

def _normalize(sequence: str) -> str:
    symbols: list[str] = []
    for char in sequence:
        if char.isspace():
            continue
        if not char.isascii():
            raise ComparisonError(f"non-ASCII symbol {char!r} is not allowed")
        symbols.append(char.upper())
    return "".join(symbols)


def _validate(sequence: str, alphabet: str, label: str) -> None:
    canonical, ambiguous = _alphabet_tables(alphabet)
    allowed = canonical | set(ambiguous) | {"-"}
    for position, symbol in enumerate(sequence, start=1):
        if symbol not in allowed:
            raise ComparisonError(
                f"{label}: invalid {alphabet} symbol {symbol!r} at aligned position {position}"
            )


def _alphabet_tables(alphabet: str):
    if alphabet == "dna":
        return _DNA_CANONICAL, _DNA_AMBIGUOUS
    if alphabet == "rna":
        return _RNA_CANONICAL, _RNA_AMBIGUOUS
    return _PROTEIN_CANONICAL, _PROTEIN_AMBIGUOUS


def _candidate_set(
    symbol: str,
    canonical: frozenset[str],
    ambiguous: dict[str, frozenset[str]],
) -> frozenset[str]:
    if symbol in ambiguous:
        return ambiguous[symbol]
    if symbol in canonical:
        return frozenset((symbol,))
    return frozenset()


def _percent(numerator: int, denominator: int) -> float | None:
    return None if denominator == 0 else round(numerator * 100.0 / denominator, 6)
