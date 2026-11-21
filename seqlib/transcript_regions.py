"""Local nucleotide similarity and antisense complementarity in transcript coordinates."""
from bisect import bisect_left
from collections import Counter

from seqlib.anchor_matcher import normalize_sequence, reverse_complement
from seqlib.result_summary import integer


def _dna(text):
    sequence = normalize_sequence(text).sequence.replace('U', 'T')
    if not sequence or any(base not in 'ACGTRYSWKMBDHVN' for base in sequence):
        raise ValueError('expected a nonempty RNA/DNA nucleotide sequence')
    return sequence


def _trace(query, target, trace, i, j, reverse, score):
    width = len(target) + 1
    end_i, end_j = i, j
    qrow, trow = [], []
    while (step := trace[i * width + j]):
        qrow.append(query[i - 1] if step != 3 else '-')
        trow.append(target[j - 1] if step != 2 else '-')
        i -= step != 3
        j -= step != 2
    qrow = ''.join(reversed(qrow)).replace('T', 'U')
    trow = ''.join(reversed(trow)).replace('T', 'U')
    matches = sum(a == b and a in 'ACGU' for a, b in zip(qrow, trow))
    gaps = qrow.count('-') + trow.count('-')
    span = [len(query) - end_i + 1, len(query) - i] if reverse else [i + 1, end_i]
    return {'score': score, 'query_span': span, 'target_span': [j + 1, end_j],
            'aligned_query': qrow, 'aligned_target': trow, 'matches': matches,
            'substitutions': len(qrow) - matches - gaps, 'gap_columns': gaps,
            'identity_fraction': matches / len(qrow), 'query_coverage': (end_i - i) / len(query)}


def _local(query, target, reverse, minimum, top):
    """Keep the best traceback ending at each target position; rank disjoint sites."""
    if minimum > 2 * min(len(query), len(target)):
        return {'hits': [], 'more_hits': False}
    width = len(target) + 1
    trace = bytearray(width * (len(query) + 1))
    previous = [0] * width
    best = [0] * width
    ends = [0] * width
    for i, a in enumerate(query, 1):
        current = [0] * width
        row_offset = i * width
        for j, b in enumerate(target, 1):
            diag = previous[j - 1] + (2 if a == b and a in 'ACGT' else -2)
            up = previous[j] - 3
            left = current[j - 1] - 3
            value = max(0, diag, up, left)
            current[j] = value
            if value:
                trace[row_offset + j] = 1 if diag == value else 2 if up == value else 3
            if value > best[j]:
                best[j], ends[j] = value, i
        previous = current

    taken, hits = [], []
    endpoints = (j for j in range(1, width) if best[j] >= minimum)
    for j in sorted(endpoints, key=lambda col: -best[col]):
        # Reject endpoints inside an accepted interval before tracing.
        k = bisect_left(taken, (j + 1,))
        if k and taken[k - 1][1] >= j:
            continue
        hit = _trace(query, target, trace, ends[j], j, reverse, best[j])
        start, end = hit['target_span']
        if k and taken[k - 1][1] >= start:
            continue
        if len(hits) == top:
            return {'hits': hits, 'more_hits': True}
        taken.insert(k, (start, end))
        hits.append(hit)
    return {'hits': hits, 'more_hits': False}


def _seed_sites(query, target, top):
    # Start with one site per occurrence of the reverse complement of bases 2..7.
    core = reverse_complement(query[1:7])
    extension = reverse_complement(query[7])
    sites = []
    counts = Counter()
    offset = target.find(core)
    while offset >= 0:
        end = offset + len(core)
        m8 = offset > 0 and target[offset - 1] == extension
        a1 = target[end:end + 1] == 'A'
        kind = ('8mer' if m8 and a1 else '7mer-m8' if m8 else '7mer-A1' if a1 else '6mer')
        counts[kind] += 1
        if counts[kind] <= top:
            sites.append({'class': kind, 'target_span': [offset + 1 - int(m8), end + int(a1)]})
        offset = target.find(core, offset + 1)
    return {'sites': sites, 'counts': dict(counts), 'omitted': sum(counts.values()) - len(sites)}


def transcript_regions(query, target, *, min_score=None, top=10):
    """Search both orientations; input and output spans are 1-based inclusive."""
    query, target = _dna(query), _dna(target)
    if (len(query) + 1) * (len(target) + 1) > 8_000_000:
        raise ValueError('local alignment exceeds 8 million cells; use a shorter query region or target')
    minimum = min(12, 2 * len(query)) if min_score is None else integer(min_score, 'min_score', 1)
    integer(top, 'top', 1)
    rc = reverse_complement(query)
    result = {'schema': 'transcript-regions/1', 'coordinate_convention': '1-based-inclusive',
              'query_length': len(query), 'target_length': len(target), 'min_score': minimum,
              'similar': _local(query, target, False, minimum, top),
              'reverse_complement': _local(rc, target, True, minimum, top)}
    # A short query can be any RNA fragment. The label describes a pattern, not a miRNA call.
    if 18 <= len(query) <= 26 and set(query[1:8]) <= set('ACGT'):
        result['seed_patterns'] = _seed_sites(query, target, top)
    else:
        result['seed_patterns'] = None
    return result
