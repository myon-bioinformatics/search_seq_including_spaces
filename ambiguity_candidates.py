"""Bounded IUPAC candidate generation; pairwise coverage is separate from sampling."""
import hashlib
import math
import random
from itertools import combinations, product

from anchor_matcher import _BASES, normalize_sequence
from result_summary import integer

MODES = ('preserve', 'regex', 'exhaustive', 'sample', 'pairwise')


def generate_candidates(text, *, mode='preserve', limit=10000, seed=0):
    integer(limit, 'limit', 1)
    integer(seed, 'seed')
    sequence = normalize_sequence(text).sequence
    if not sequence or any(c not in _BASES for c in sequence):
        raise ValueError('expected a nonempty IUPAC DNA sequence')
    if mode not in MODES:
        raise ValueError('unsupported candidate mode')
    factors = [(i, _BASES[c]) for i, c in enumerate(sequence) if len(_BASES[c]) > 1]
    total = math.prod(len(values) for _, values in factors)
    report = {'schema': 'sequence-candidates/1', 'mode': mode, 'seed': seed, 'limit': limit,
              'source_sha256': hashlib.sha256(text.encode('utf-8')).hexdigest(),
              'ambiguous_positions': [i + 1 for i, _ in factors], 'possible_candidates': total,
              'coverage': None, 'candidates': []}
    if mode in ('preserve', 'regex'):
        report['sequence'] = sequence
        if mode == 'regex':
            report['pattern'] = ''.join(v if len(v) == 1 else '[' + v + ']' for v in map(_BASES.__getitem__, sequence))
        report.update(generated_candidates=0, omitted_candidates=None)
        return report
    rng = random.Random(seed)
    if mode == 'exhaustive':
        if total > limit:
            raise ValueError(f'exhaustive expansion needs {total} candidates, exceeds limit {limit}')
        assignments = product(*(values for _, values in factors))
    elif mode == 'sample':
        # Floyd sampling supports arbitrary-sized combination spaces without range/len overflow.
        numbers = set()
        for end in range(total - min(total, limit), total):
            candidate = rng.randrange(end + 1)
            numbers.add(end if candidate in numbers else candidate)
        assignments = []
        for number in sorted(numbers):
            values = []
            for _, levels in reversed(factors):
                number, offset = divmod(number, len(levels))
                values.append(levels[offset])
            assignments.append(tuple(reversed(values)))
    else:
        if len(factors) > 64:
            raise ValueError('pairwise supports at most 64 ambiguous positions')
        if len(factors) < 2:
            assignments = list(product(*(v for _, v in factors)))
            if len(assignments) > limit:
                raise ValueError('candidate limit cannot cover all levels')
            report['coverage'] = {'strength': min(2, len(factors)), 'required_pairs': 0, 'covered_pairs': 0, 'complete': True}
        else:
            uncovered = {(i, a, j, b) for i, j in combinations(range(len(factors)), 2)
                         for a in factors[i][1] for b in factors[j][1]}
            required = len(uncovered)
            assignments = []
            while uncovered:
                if len(assignments) >= limit:
                    raise ValueError('candidate limit reached before pairwise coverage; increase limit')
                i, a, j, b = min(uncovered)
                row = [rng.choice(values) for _, values in factors]
                row[i], row[j] = a, b
                # Fill other factors to cover as many outstanding obligations as possible.
                fixed = {i, j}
                for k in range(len(factors)):
                    if k in fixed:
                        continue
                    def score(value):
                        return sum((p, row[p], k, value) in uncovered if p < k
                                   else (k, value, p, row[p]) in uncovered for p in fixed)
                    row[k] = max(factors[k][1], key=score)
                    fixed.add(k)
                assignments.append(row)
                for p, q in combinations(range(len(factors)), 2):
                    uncovered.discard((p, row[p], q, row[q]))
            report['coverage'] = {'strength': 2, 'required_pairs': required, 'covered_pairs': required, 'complete': True}
    candidates = []
    for values in assignments:
        letters = list(sequence)
        for (position, _), value in zip(factors, values):
            letters[position] = value
        candidates.append(''.join(letters))
    report.update(candidates=candidates, generated_candidates=len(candidates), omitted_candidates=total - len(candidates))
    return report


def exact_hit_matrix(candidates, targets):
    """One int bitset per candidate. Bit j means exact literal hit in target j."""
    normalized = [normalize_sequence(t).sequence for t in targets]
    result = []
    for candidate in candidates:
        if not candidate or any(c not in 'ACGT' for c in candidate):
            raise ValueError('exact hit matrix requires concrete DNA candidates')
        result.append(sum(1 << j for j, target in enumerate(normalized) if candidate in target))
    return result


def select_representatives(candidates, targets, *, rank=8, seed=0):
    """Opt-in low-rank Boolean approximation of an EXACT hit matrix.

    Greedy biclique factors use integer bitsets. Factors never invent hits;
    uncovered true bits remain residual error. One representative per factor
    membership signature is selected. Full matrix is rechecked to report lost
    target coverage. This is not LittleBit's model quantization algorithm.
    """
    integer(rank, 'rank', 1)
    integer(seed, 'seed')
    if rank > 64 or len(candidates) > 10000 or len(targets) > 10000:
        raise ValueError('representative selection limit exceeded')
    matrix = exact_hit_matrix(candidates, targets)
    residual = matrix.copy()
    factors, memberships = [], [0] * len(matrix)
    rng = random.Random(seed)
    order = list(range(len(matrix)))
    rng.shuffle(order)
    for bit in range(rank):
        if not any(residual):
            break
        chosen = max(order, key=lambda i: residual[i].bit_count())
        basis = residual[chosen]
        factors.append(basis)
        for i, row in enumerate(matrix):
            if row & basis == basis:
                memberships[i] |= 1 << bit
                residual[i] &= ~basis
    groups = {}
    for i, signature in enumerate(memberships):
        groups.setdefault(signature, []).append(i)
    selected = sorted(max(indices, key=lambda i: (matrix[i].bit_count(), -i)) for indices in groups.values())
    covered, retained = 0, 0
    for row in matrix:
        covered |= row
    for i in selected:
        retained |= matrix[i]
    return {'schema': 'sequence-representatives/1', 'rank': rank, 'seed': seed,
            'method': 'greedy-boolean-bicliques', 'exact_evaluation': 'literal-substring',
            'selected_indices': selected, 'candidates': [candidates[i] for i in selected],
            'input_candidates': len(candidates), 'omitted_candidates': len(candidates) - len(selected),
            'residual_true_bits': sum(row.bit_count() for row in residual),
            'lost_target_hits': (covered & ~retained).bit_count(),
            'factor_count': len(factors), 'pairwise_coverage_preserved': False,
            'approximation': True}
