"""Read one Stockholm alignment and locate anchor hits against a column annotation."""
from collections import defaultdict

from seqlib.anchor_matcher import find_matches


def read_stockholm(handle):
    """Read sequence rows and #=GC tracks from one blank-line-delimited alignment."""
    if handle.readline().strip() != "# STOCKHOLM 1.0":
        raise ValueError("expected # STOCKHOLM 1.0")
    parts = defaultdict(list)
    block = {}

    def finish_block():
        if not block:
            return
        if len({len(segment) for segment in block.values()}) != 1:
            raise ValueError("unequal widths within a Stockholm block")
        if parts and parts.keys() != block.keys():
            raise ValueError("each block must contain the same sequence IDs and #=GC tracks")
        for key, segment in block.items():
            parts[key].append(segment)
        block.clear()

    for line in handle:
        line = line.strip()
        if line == "//":
            finish_block()
            break
        if not line:
            finish_block()
            continue
        if line.startswith("#=GC ") or line.startswith("#=GC\t"):
            _, key, segment = line.split()
            kind = "track"
        elif line.startswith("#"):
            continue
        else:
            key, segment = line.split()
            kind = "sequence"
            if not all(c.isascii() and (c.isalpha() or c in "-.") for c in segment):
                raise ValueError(f"unsupported alignment symbol in {key}")
        if (kind, key) in block:
            raise ValueError(f"duplicate {kind} {key} within a Stockholm block")
        block[kind, key] = segment
    else:
        raise ValueError("missing Stockholm terminator //")
    if any(line.strip() for line in handle):
        raise ValueError("expected a single Stockholm alignment")
    sequences = {key: "".join(v).replace(".", "-") for (kind, key), v in parts.items() if kind == "sequence"}
    if not sequences:
        raise ValueError("alignment contains no sequence rows")
    tracks = {key: "".join(v) for (kind, key), v in parts.items() if kind == "track"}
    return sequences, tracks


def annotated_matches(handle, *, track, left, right, min_gap=0, max_gap=None):
    """Search ungapped rows and return the selected track across each hit's aligned span."""
    sequences, tracks = read_stockholm(handle)
    if track not in tracks:
        raise ValueError(f"unknown column annotation {track!r}; available: {', '.join(tracks)}")
    results = []
    for name, sequence in sequences.items():
        columns = [i for i, char in enumerate(sequence, 1) if char != "-"]
        report = find_matches(sequence.replace("-", ""), left, right, min_gap=min_gap, max_gap=max_gap)
        for hit in report["matches"]:
            start, end = hit["span"]
            first, last = columns[start - 1], columns[end - 1]
            results.append({"record_id": name, "span": [start, end], "columns": [first, last],
                            "aligned_sequence": sequence[first - 1:last],
                            "annotation": tracks[track][first - 1:last]})
    return {"schema": "annotated-matches/1", "coordinate_convention": "1-based-inclusive",
            "track": track, "records_searched": len(sequences), "results": results}
