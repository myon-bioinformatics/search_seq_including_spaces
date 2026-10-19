"""Deterministic nucleotide features, marker clues and a bounded RNA baseline."""
from collections import Counter
import hashlib
import math

from itertools import product
from anchor_matcher import normalize_sequence, _BASES
from comparison import compare_aligned, _alphabet_tables
from result_summary import integer

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

_CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


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


def translate_nucleotide(text, *, frame=0, table=1, strand='+', to_stop=False, stop_symbol='*'):
    if table != 1:
        raise ValueError('only standard genetic code table 1 is supported')
    if isinstance(frame, str) and frame.isdigit():
        frame = int(frame)
    if frame not in (0, 1, 2, 'all'):
        raise ValueError("frame must be 0, 1, 2 or 'all'")
    if strand not in ('+', '-', 'both'):
        raise ValueError("strand must be '+', '-' or 'both'")
    if not isinstance(to_stop, bool):
        raise ValueError('to_stop must be a boolean')
    if not isinstance(stop_symbol, str) or len(stop_symbol) != 1:
        raise ValueError('stop_symbol must be a single character')

    raw_seq = normalize_sequence(text).sequence.replace('U', 'T')
    if not raw_seq or any(c not in _BASES for c in raw_seq):
        raise ValueError('expected nonempty DNA/RNA nucleotide sequence')

    seq_len = len(raw_seq)
    strands = ('+', '-') if strand == 'both' else (strand,)
    frames = (0, 1, 2) if frame == 'all' else (frame,)
    translations = []

    for s in strands:
        oriented = raw_seq if s == '+' else raw_seq.translate(_COMPLEMENT)[::-1]
        for f in frames:
            amino_acids = []
            for i in range(f, seq_len - 2, 3):
                triplet = oriented[i:i + 3]
                bases = [_BASES[c] for c in triplet]
                targets = {_CODON_TABLE[b0 + b1 + b2] for b0, b1, b2 in product(*bases)}
                aa = list(targets)[0] if len(targets) == 1 else 'X'
                if to_stop and aa == '*':
                    break
                amino_acids.append(aa)

            # Stop identity comes from the code table, never the display symbol.
            protein_seq = ''.join(stop_symbol if aa == '*' else aa for aa in amino_acids)
            stop_count = amino_acids.count('*')
            if to_stop:
                has_internal = False
            else:
                has_internal = '*' in amino_acids[:-1]

            num_codons = len(amino_acids)
            if num_codons > 0:
                if s == '+':
                    span = [f + 1, f + num_codons * 3]
                else:
                    span = [seq_len - f - num_codons * 3 + 1, seq_len - f]
            else:
                span = None

            translations.append({
                'strand': s,
                'frame': f,
                'protein_sequence': protein_seq,
                'length': num_codons,
                'stop_count': stop_count,
                'has_internal_stop': has_internal,
                'nucleotide_span': span,
            })

    return {
        'schema': 'nucleotide-translation/1',
        'table': table,
        'to_stop': to_stop,
        'coordinate_convention': '1-based-inclusive',
        'translations': translations,
        'source_sha256': hashlib.sha256(text.encode('utf-8')).hexdigest(),
    }
