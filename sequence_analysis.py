"""Deterministic nucleotide features, marker clues and a bounded RNA baseline."""
from collections import Counter
import hashlib
import math

from anchor_matcher import normalize_sequence, _BASES
from comparison import compare_aligned, _alphabet_tables
from result_summary import integer


def _nucleotide(text, alphabet):
    if alphabet not in ('dna', 'rna'):
        raise ValueError('alphabet must be dna or rna')
    sequence = normalize_sequence(text).sequence
    canonical, ambiguous = _alphabet_tables(alphabet)
    if not sequence or any(c not in canonical and c not in ambiguous for c in sequence):
        raise ValueError(f'expected nonempty {alphabet} sequence')
    return sequence, canonical, ambiguous


def nucleotide_features(text, *, alphabet='dna', k=2, codon_frame=None):
    integer(k, 'k', 1)
    if k > 6:
        raise ValueError('k must be <= 6')
    sequence, canonical, ambiguous = _nucleotide(text, alphabet)
    counts = Counter(sequence)
    codons = None
    if codon_frame is not None:
        integer(codon_frame, 'codon_frame')
        if codon_frame > 2: raise ValueError('codon_frame must be 0, 1 or 2')
        triplets = [sequence[i:i+3] for i in range(codon_frame, len(sequence)-2, 3)]
        valid = Counter(t for t in triplets if set(t) <= canonical)
        codons = {'frame': codon_frame, 'counts': dict(sorted(valid.items())),
                  'excluded_triplets': len(triplets)-sum(valid.values()),
                  'trailing_bases': max(0, len(sequence)-codon_frame) % 3}
    known = sum(counts[c] for c in canonical)
    kmers = Counter(sequence[i:i + k] for i in range(len(sequence) - k + 1)
                    if set(sequence[i:i + k]) <= canonical)
    runs, start = [], 0
    for i in range(1, len(sequence) + 1):
        if i == len(sequence) or sequence[i] != sequence[start]:
            if i - start > 1:
                runs.append({'symbol': sequence[start], 'start': start + 1, 'end': i, 'length': i - start})
            start = i
    entropy = -sum((counts[c] / known) * math.log2(counts[c] / known) for c in canonical if counts[c]) if known else None
    return {'schema': 'nucleotide-features/1', 'alphabet': alphabet, 'length': len(sequence),
            'counts': dict(sorted(counts.items())), 'canonical_count': known,
            'ambiguous_count': len(sequence) - known, 'ambiguous_fraction': (len(sequence) - known) / len(sequence),
            'gc_fraction': (counts['G'] + counts['C']) / known if known else None,
            'at_or_au_fraction': (counts['A'] + counts['T' if alphabet == 'dna' else 'U']) / known if known else None,
            'entropy_bits': entropy, 'entropy_denominator': 'canonical residues only',
            'k': k, 'kmers': dict(sorted(kmers.items())), 'codons': codons,
            'kmer_denominator': sum(kmers.values()), 'excluded_kmers': max(0, len(sequence) - k + 1) - sum(kmers.values()),
            'homopolymers': runs, 'coordinate_convention': '1-based-inclusive',
            'source_sha256': hashlib.sha256(text.encode('utf-8')).hexdigest()}


def nucleotide_pairwise(a, b, *, alphabet='dna', gap_policy='exclude', ambiguity_policy='exclude'):
    result = compare_aligned(a, b, alphabet=alphabet, gap_policy=gap_policy, ambiguity_policy=ambiguity_policy)
    if alphabet not in ('dna', 'rna'):
        raise ValueError('transition/transversion requires dna or rna')
    canonical, _ = _alphabet_tables(alphabet)
    a, b = ''.join(a.split()).upper(), ''.join(b.split()).upper()
    transitions = transversions = 0
    for left, right in zip(a, b):
        if left not in canonical or right not in canonical or left == right:
            continue
        if {left, right} in ({'A', 'G'}, {'C', 'T' if alphabet == 'dna' else 'U'}):
            transitions += 1
        else:
            transversions += 1
    result['substitutions'] = {'transitions': transitions, 'transversions': transversions,
                              'ratio': transitions / transversions if transversions else None,
                              'policy': 'canonical non-gap substitutions only'}
    return result


def marker_clues(text, *, kind, motif=None, min_length=6, position=None, allele=None,
                 repeat_length=4, min_gap=0, max_gap=100, max_results=10000):
    sequence, _, _ = _nucleotide(text, 'dna')
    integer(max_results, 'max_results', 1)
    results = []
    def add(row):
        if len(results) >= max_results:
            raise ValueError('preset result limit exceeded')
        results.append(row)
    if kind in ('poly-a', 'str'):
        integer(min_length, 'min_length', 1)
        motif = 'A' if kind == 'poly-a' else motif
        if not isinstance(motif, str):
            raise ValueError('str requires motif')
        motif, canonical, _ = _nucleotide(motif, 'dna')
        if not set(motif) <= canonical:
            raise ValueError('motif must be concrete DNA')
        width, i = len(motif), 0
        while i <= len(sequence) - width:
            if sequence[i:i + width] != motif:
                i += 1
                continue
            end = i + width
            while sequence[end:end + width] == motif:
                end += width
            copies = (end - i) // width
            if (end - i if kind == 'poly-a' else copies) >= min_length:
                add({'start': i + 1, 'end': end, 'motif': motif, 'copies': copies})
            i = end
    elif kind == 'snp':
        integer(position, 'position', 1)
        if position > len(sequence) or allele not in ('A', 'C', 'G', 'T'):
            raise ValueError('SNP requires an in-range 1-based position and concrete allele')
        observed = sequence[position - 1]
        add({'position': position, 'observed': observed, 'allele': allele,
             'exact_match': observed == allele, 'compatible': allele in _BASES[observed]})
    elif kind == 'tsd':
        integer(repeat_length, 'repeat_length', 1)
        integer(min_gap, 'min_gap')
        integer(max_gap, 'max_gap', min_gap)
        if max_gap > 10000 or len(sequence) * (max_gap - min_gap + 1) > 10000000:
            raise ValueError('TSD work limit exceeded; narrow input or gap range')
        for start in range(len(sequence) - 2 * repeat_length - min_gap + 1):
            left = sequence[start:start + repeat_length]
            if set(left) - set('ACGT'):
                continue
            for gap in range(min_gap, min(max_gap, len(sequence) - start - 2 * repeat_length) + 1):
                right = start + repeat_length + gap
                if sequence[right:right + repeat_length] == left:
                    add({'left_span': [start + 1, start + repeat_length],
                         'right_span': [right + 1, right + repeat_length], 'repeat': left, 'gap': gap})
    else:
        raise ValueError('unknown preset kind')
    return {'schema': 'sequence-marker-clues/1', 'kind': kind, 'coordinate_convention': '1-based-inclusive',
            'results': results, 'complete': True, 'interpretation': 'sequence clues only; no pseudogene, identity or ancestry inference',
            'source_sha256': hashlib.sha256(text.encode('utf-8')).hexdigest()}


def fold_rna(text, *, min_loop=3, max_length=300, wobble=True):
    """Maximum noncrossing pairs, O(n^3)/O(n^2), not a free-energy model.

    Ties prefer unpaired right endpoint, then the smallest partner. Input is
    canonical RNA only: ambiguity is not silently resolved.
    """
    integer(min_loop, 'min_loop')
    integer(max_length, 'max_length', 1)
    if max_length > 1000 or not isinstance(wobble, bool):
        raise ValueError('max_length <= 1000 and a boolean wobble flag are required')
    sequence, canonical, _ = _nucleotide(text, 'rna')
    n = len(sequence)
    if n > max_length or set(sequence) - canonical:
        raise ValueError('RNA exceeds length limit or contains ambiguous bases')
    allowed = {'AU', 'UA', 'GC', 'CG'} | ({'GU', 'UG'} if wobble else set())
    score = [[0] * n for _ in range(n)]
    trace = {}
    def value(i, j):
        return score[i][j] if i <= j else 0
    for width in range(1, n):
        for i in range(n - width):
            j = i + width
            best, partner = value(i, j - 1), None
            for k in range(i, j - min_loop):
                if sequence[k] + sequence[j] in allowed:
                    candidate = value(i, k - 1) + value(k + 1, j - 1) + 1
                    if candidate > best:
                        best, partner = candidate, k
            score[i][j] = best
            trace[i, j] = partner
    pairs, stack = [], [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if i >= j:
            continue
        k = trace[i, j]
        if k is None:
            stack.append((i, j - 1))
        else:
            pairs.append((k + 1, j + 1))
            stack.extend(((i, k - 1), (k + 1, j - 1)))
    pairs.sort()
    brackets = ['.'] * n
    for i, j in pairs:
        brackets[i - 1], brackets[j - 1] = '(', ')'
    pair_set = set(pairs)
    stems = []
    for i, j in pairs:
        if (i - 1, j + 1) in pair_set:
            continue
        width = 1
        while (i + width, j - width) in pair_set:
            width += 1
        stems.append({'outer_pair': [i, j], 'pair_count': width})
    loops = []
    for i, j in pairs:
        children = [(a, b) for a, b in pairs if i < a < b < j
                    and not any(i < x < a < b < y < j for x, y in pairs)]
        if not children:
            loops.append({'kind': 'hairpin', 'closing_pair': [i, j], 'unpaired': j - i - 1})
        elif len(children) == 1:
            a, b = children[0]
            left, right = a - i - 1, j - b - 1
            if left or right:
                loops.append({'kind': 'internal' if left and right else 'bulge',
                              'closing_pair': [i, j], 'unpaired_left': left, 'unpaired_right': right})
    return {'schema': 'rna-nussinov/1', 'sequence': sequence, 'pair_count': len(pairs),
            'pairs': [list(pair) for pair in pairs], 'dot_bracket': ''.join(brackets), 'stems': stems, 'loops': loops,
            'min_loop': min_loop, 'wobble': wobble, 'max_length': max_length,
            'coordinate_convention': '1-based-inclusive', 'model': 'maximum-pairs; pseudoknot-free; no energy parameters',
            'tie_break': 'unpaired-right then smallest partner'}
