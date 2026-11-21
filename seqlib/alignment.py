"""Global alignment of unaligned nucleotide records, retaining the supplied orientation."""
from seqlib.comparison import _normalize, _alphabet_tables


def _pair(anchor, sequence):
    """Needleman–Wunsch with linear gaps and deterministic diagonal/up/left tie order."""
    width = len(sequence) + 1
    trace = bytearray((len(anchor) + 1) * width)
    prev = [-2 * j for j in range(width)]
    trace[1:width] = bytes([3]) * len(sequence)
    for i, a in enumerate(anchor, 1):
        cur = [-2 * i] + [0] * len(sequence)
        offset = i * width
        trace[offset] = 2
        for j, b in enumerate(sequence, 1):
            match = (a == b and a in 'ACGTU') or (a in 'TU' and b in 'TU')
            diag = prev[j - 1] + (2 if match else -1)
            up, left = prev[j] - 2, cur[j - 1] - 2
            best = max(diag, up, left)
            cur[j] = best
            trace[offset + j] = 1 if best == diag else 2 if best == up else 3
        prev = cur

    # Bucket insertions by the number of anchor residues preceding them.
    aligned = ['-'] * len(anchor)
    insertions = [[] for _ in range(len(anchor) + 1)]
    i, j = len(anchor), len(sequence)
    while i or j:
        step = trace[i * width + j]
        if step == 3:
            j -= 1
            insertions[i].append(sequence[j])
        else:
            i -= 1
            if step == 1:
                j -= 1
                aligned[i] = sequence[j]
    return aligned, [''.join(reversed(chunk)) for chunk in insertions]


def align_records(records):
    """Align (id, ungapped DNA/RNA) pairs around a deterministic median-length anchor."""
    canonical, ambiguous = _alphabet_tables('dna')
    allowed = canonical | set(ambiguous) | {'U'}
    inputs = {}
    for name, raw in records:
        if not name or name.startswith('#') or name == '//' or any(c.isspace() for c in name):
            raise ValueError(f'invalid alignment record ID: {name!r}')
        if name in inputs:
            raise ValueError(f'duplicate alignment record ID: {name!r}')
        sequence = _normalize(raw)
        if not sequence or set(sequence) - allowed:
            raise ValueError(f'{name}: expected ungapped IUPAC DNA/RNA')
        inputs[name] = sequence
    if len(inputs) < 2:
        raise ValueError('alignment requires at least two sequences')

    by_length = sorted(inputs, key=lambda name: (len(inputs[name]), name))
    length = len(inputs[by_length[(len(inputs) - 1) // 2]])
    anchor_id = next(name for name in by_length if len(inputs[name]) == length)
    anchor = inputs[anchor_id]
    if sum((length + 1) * (len(seq) + 1) for name, seq in inputs.items() if name != anchor_id) > 8_000_000:
        raise ValueError('alignment exceeds 8 million dynamic-programming cells')
    widths = [0] * (length + 1)
    parts = {}
    for name, seq in inputs.items():
        if name == anchor_id:
            parts[name] = (anchor, [''] * (length + 1))
        else:
            bases, gaps = _pair(anchor, seq)
            parts[name] = bases, gaps
            widths = [max(width, len(gap)) for width, gap in zip(widths, gaps)]
    if (length + sum(widths)) * len(inputs) > 8_000_000:
        raise ValueError('alignment output exceeds 8 million symbols')
    rows = {}
    for name, (bases, gaps) in parts.items():
        pieces = []
        for i, base in enumerate(bases):
            pieces.extend((gaps[i].ljust(widths[i], '-'), base))
        pieces.append(gaps[-1].ljust(widths[-1], '-'))
        rows[name] = ''.join(pieces)
    return rows
