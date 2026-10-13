"""Versioned custom anchor/gap pattern packs; strict validation before I/O."""
from anchor_matcher import find_matches
from structured_readers import strict_json

SCHEMA = 'sequence-patterns/1'


def load_patterns(text):
    doc = strict_json(text)
    if not isinstance(doc, dict) or set(doc) != {'schema', 'patterns'} or doc['schema'] != SCHEMA:
        raise ValueError('expected sequence-patterns/1 with schema and patterns')
    patterns = doc['patterns']
    if not isinstance(patterns, list) or not 1 <= len(patterns) <= 100:
        raise ValueError('patterns must contain 1..100 definitions')
    seen = set()
    for pattern in patterns:
        allowed = {'id', 'left', 'right', 'min_gap', 'max_gap', 'iupac', 'strand', 'max_matches'}
        if not isinstance(pattern, dict) or set(pattern) - allowed or not {'id', 'left', 'right'} <= set(pattern):
            raise ValueError('invalid pattern keys; require id/left/right and matcher options only')
        identity = pattern['id']
        if not isinstance(identity, str) or not identity.strip() or identity in seen:
            raise ValueError('pattern ids must be unique nonempty strings')
        seen.add(identity)
        find_matches('', **{k: v for k, v in pattern.items() if k != 'id'})
    return patterns
