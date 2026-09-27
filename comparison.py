"""Aligned pairwise sequence identity/coverage (stdlib only).

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
    lengths must then be equal. \`start\`/\`end\` use 1-based inclusive aligned
    coordinates. \`end=None\` means the final aligned column.
    """
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
