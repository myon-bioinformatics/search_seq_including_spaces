"""Stdlib anchor/gap matching; coordinates are 1-based inclusive."""
from dataclasses import dataclass
import hashlib

SCHEMA = "anchor-gap/1"
__all__ = ["SCHEMA", "NormalizedSequence", "normalize_sequence", "reverse_complement", "find_matches"]

_BASES = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
          "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT",
          "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


@dataclass(frozen=True)
class NormalizedSequence:
    sequence: str
    # One original-text character position for every normalized residue.
    source_positions: tuple


def normalize_sequence(text):
    """Uppercase ASCII letters and remove whitespace, preserving source offsets.

    Input is sequence text only, not a FASTA/FASTQ file or alignment. Reject
    headers, digits, gap markers and non-ASCII letters rather than lose residues.
    Source coordinates count Python string characters, not encoded bytes.
    """
    if not isinstance(text, str):
        raise TypeError("sequence must be text")
    letters, positions = [], []
    for position, char in enumerate(text, 1):
        if char.isspace():
            continue
        if not ("a" <= char <= "z" or "A" <= char <= "Z"):
            raise ValueError("sequence must contain ASCII letters and whitespace only")
        letters.append(char.upper())
        positions.append(position)
    return NormalizedSequence("".join(letters), tuple(positions))


def _dna(sequence):
    if any(char not in _BASES for char in sequence):
        raise ValueError("DNA mode requires IUPAC DNA symbols (RNA U is not accepted)")


def reverse_complement(text):
    """Return normalized IUPAC DNA reverse complement; retains ambiguity symbols."""
    sequence = normalize_sequence(text).sequence
    _dna(sequence)
    return sequence.translate(_COMPLEMENT)[::-1]


def _integer(value, name, *, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(name + " must be an integer >= " + str(minimum))


def find_matches(text, left, right, *, min_gap=0, max_gap=None, iupac=False,
                 strand="+", max_matches=10000):
    """Return deterministic overlapping anchor pairs and oriented middle strings.

    max_gap=None means exactly min_gap. Literal mode compares normalized letters
    (including N) literally. IUPAC mode tests nonempty intersection of possible
    DNA base sets in both anchors and input: compatibility, not a proven base.
    Results use forward-reference coordinates even on '-'. '+' results precede
    '-' results; within a strand, scan oriented starts then ascending gap length.
    A pair is emitted once per strand; palindromes may occur on both strands.
    Exceeding max_matches raises instead of returning an apparently complete
    prefix. Runtime is proportional to sequence length, feasible gap range and
    anchor lengths; this baseline is not a genome-scale indexed search engine.
    """
    _integer(min_gap, "min_gap")
    max_gap = min_gap if max_gap is None else max_gap
    _integer(max_gap, "max_gap")
    _integer(max_matches, "max_matches", minimum=1)
    if max_gap < min_gap:
        raise ValueError("max_gap must be >= min_gap")
    if not isinstance(iupac, bool) or strand not in ("+", "-", "both"):
        raise ValueError("invalid IUPAC flag or strand")
    normalized = normalize_sequence(text)
    left = normalize_sequence(left).sequence
    right = normalize_sequence(right).sequence
    if not left or not right:
        raise ValueError("anchors must not be empty")
    if iupac or strand != "+":
        for value in (normalized.sequence, left, right):
            _dna(value)

    def matches(sequence, offset, anchor):
        observed = sequence[offset:offset + len(anchor)]
        if not iupac:
            return observed == anchor
        return all(set(_BASES[a]).intersection(_BASES[b])
                   for a, b in zip(observed, anchor))

    result = []
    length = len(normalized.sequence)
    for direction in (("+", "-") if strand == "both" else (strand,)):
        sequence = (normalized.sequence if direction == "+"
                    else reverse_complement(normalized.sequence))
        positions = (normalized.source_positions if direction == "+"
                     else normalized.source_positions[::-1])

        def span(start, end):
            if start == end:
                return None
            return [start + 1, end] if direction == "+" else [length - end + 1, length - start]

        def source_span(start, end):
            if start == end:
                return None
            return [min(positions[start], positions[end - 1]), max(positions[start], positions[end - 1])]

        for start in range(length - len(left) - min_gap - len(right) + 1):
            if not matches(sequence, start, left):
                continue
            middle_start = start + len(left)
            limit = min(max_gap, length - middle_start - len(right))
            for gap in range(min_gap, limit + 1):
                right_start = middle_start + gap
                end = right_start + len(right)
                if not matches(sequence, right_start, right):
                    continue
                if len(result) >= max_matches:
                    raise ValueError("match limit exceeded; narrow the search")
                result.append({"strand": direction, "span": span(start, end),
                               "source_span": source_span(start, end),
                               "left_span": span(start, middle_start),
                               "middle_span": span(middle_start, right_start),
                               "right_span": span(right_start, end),
                               "middle_source_span": source_span(middle_start, right_start),
                               "left": sequence[start:middle_start],
                               "middle": sequence[middle_start:right_start],
                               "right": sequence[right_start:end], "gap": gap})
    return {"schema": SCHEMA, "coordinate_convention": "1-based-inclusive",
            "source_coordinate_unit": "python-string-character", "sequence_length": length,
            "iupac": iupac, "matches": result, "complete": True,
            "provenance": {"source_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
                           "source_length": len(text), "left_anchor": left, "right_anchor": right,
                           "min_gap": min_gap, "max_gap": max_gap, "strand": strand}}
