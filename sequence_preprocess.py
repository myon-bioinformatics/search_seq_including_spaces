"""Auditable sequence normalization and FASTA/JSONL/MEGA export, no aligner."""
from dataclasses import asdict, replace
import hashlib
import json
import re

from comparison import _alphabet_tables
from result_summary import integer


def preprocess(records, *, alphabet='dna', normalize_ut=False, deduplicate='none',
               min_length=1, max_length=None, max_ambiguity=1.0, aligned=False,
               include_ids=None, motif=None):
    canonical, ambiguous = _alphabet_tables(alphabet)
    integer(min_length, 'min_length')
    if max_length is not None:
        integer(max_length, 'max_length', min_length)
    if not isinstance(max_ambiguity, (int, float)) or isinstance(max_ambiguity, bool) or not 0 <= max_ambiguity <= 1:
        raise ValueError('max_ambiguity must be in [0, 1]')
    if deduplicate not in ('none', 'id', 'sequence'):
        raise ValueError('invalid deduplicate mode')
    if normalize_ut and alphabet == 'protein':
        raise ValueError('U/T normalization is only defined for nucleotides')
    if motif is not None and (not isinstance(motif, str) or not motif or set(motif) - canonical):
        raise ValueError('motif must use canonical uppercase symbols')
    output, receipts, seen, names, length = [], [], set(), set(), None
    for index, record in enumerate(records, 1):
        raw = record.sequence
        stripped = ''.join(raw.split())
        if any(not ('A' <= c <= 'Z' or 'a' <= c <= 'z' or c in '-?') for c in stripped):
            raise ValueError('sequence must contain ASCII letters, gaps or missing symbols')
        sequence = stripped.upper()
        if normalize_ut:
            sequence = sequence.replace('U', 'T') if alphabet == 'dna' else sequence.replace('T', 'U')
        if not sequence or set(sequence) - (canonical | set(ambiguous) | {'-', '?'}):
            raise ValueError(f'record {record.id!r} contains invalid {alphabet} symbols or is empty')
        ungapped = len(sequence) - sequence.count('-')
        ambiguous_count = sum(sequence.count(c) for c in ambiguous) + sequence.count('?')
        ratio = ambiguous_count / ungapped if ungapped else 1.0
        digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
        reason = None
        if len(sequence) < min_length or (max_length is not None and len(sequence) > max_length):
            reason = 'length'
        elif ratio > max_ambiguity:
            reason = 'ambiguity'
        elif include_ids is not None and record.id not in include_ids:
            reason = 'id-filter'
        elif motif is not None and motif not in sequence:
            reason = 'motif-filter'
        key = record.id if deduplicate == 'id' else digest
        if reason is None and deduplicate != 'none' and key in seen:
            reason = 'duplicate'
        if reason is None:
            if aligned and length is not None and len(sequence) != length:
                raise ValueError('retained aligned records must have equal lengths')
            length = len(sequence)
            seen.add(key)
        name = re.sub(r'[^A-Za-z0-9_.-]', '_', record.id) or 'record'
        if not name[0].isalpha():
            name = 'seq_' + name
        stem, suffix = name, 2
        while name in names:
            name = f'{stem}_{suffix}'
            suffix += 1
        if reason is None:
            names.add(name)
            output.append(replace(record, id=name, sequence=sequence, quality=None))
        receipts.append({'input_index': index, 'input_id': record.id, 'output_id': name if reason is None else None,
                         'included': reason is None, 'reason': reason,
                         'input_sha256': hashlib.sha256(raw.encode('utf-8')).hexdigest(),
                         'output_sha256': digest, 'source': record.source, 'read_provenance': asdict(record.provenance),
                         'transforms': {'remove_whitespace': raw != ''.join(raw.split()),
                                        'uppercase': ''.join(raw.split()) != ''.join(raw.split()).upper(),
                                        'normalize_ut': normalize_ut, 'quality_removed': record.quality is not None},
                         'length': len(sequence), 'gap_count': sequence.count('-'),
                         'missing_count': sequence.count('?'), 'ambiguous_count': ambiguous_count})
    return output, {'schema': 'sequence-preprocessing/1', 'alphabet': alphabet,
                    'input_records': len(receipts), 'output_records': len(output), 'aligned': aligned,
                    'deduplicate': deduplicate, 'receipts': receipts}


def export_records(records, *, format='fasta', alphabet='dna'):
    if format not in ('fasta', 'jsonl', 'mega'):
        raise ValueError('export supports fasta, jsonl or mega')
    records = list(records)
    if format == 'jsonl':
        return ''.join(json.dumps({'schema': 'sequence-record/1', **asdict(r)}, ensure_ascii=False, sort_keys=True) + '\n' for r in records)
    for record in records:
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_.-]*', record.id):
            raise ValueError('export requires sanitized ids; use preprocess first')
        canonical, ambiguous = _alphabet_tables(alphabet)
        if not record.sequence or set(record.sequence) - (canonical | set(ambiguous) | {'-', '?'}):
            raise ValueError('export requires validated normalized sequences')
    if len({r.id for r in records}) != len(records):
        raise ValueError('export requires unique record ids')
    if format == 'mega':
        if not records or len({len(r.sequence) for r in records}) != 1:
            raise ValueError('MEGA requires nonempty equally aligned records')
        datatype = {'dna': 'DNA', 'rna': 'RNA', 'protein': 'Protein'}[alphabet]
        prefix = f'#MEGA\n!Title Sequence export;\n!Format DataType={datatype} Indel=- Missing=?;\n\n'
        return prefix + ''.join(f'#{r.id}\n{r.sequence}\n' for r in records)
    return ''.join(f'>{r.id}\n{r.sequence}\n' for r in records)
