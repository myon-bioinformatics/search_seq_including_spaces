"""Incremental anchor matching for a single giant record with bounded overlap.

New event contract, independent of the existing record-envelope API. Matches
are emitted by chunk completion, not whole-record strand order. Success is
established only by exhausting the iterator and receiving the completion event.
"""
from contextlib import contextmanager
import hashlib
import os

from anchor_matcher import find_matches, normalize_sequence, reverse_complement, _BASES, _dna
from result_summary import integer
from sequence_readers import _extension, SequenceReadError


def match_chunks(chunks, left, right, *, min_gap=0, max_gap=None, iupac=False,
                 strand='+', max_matches=10000, max_span=1000000):
    """Yield match events then a completion event. Coordinates are record-local.

    Memory O(maximum match span + caller's chunk). Source positions refer to the
    concatenation of the caller's chunks, including all their whitespace.
    """
    integer(max_span, 'max_span', 1)
    template = find_matches('', left, right, min_gap=min_gap, max_gap=max_gap,
                            iupac=iupac, strand=strand, max_matches=max_matches)
    params = template['provenance']
    width = len(params['left_anchor']) + params['max_gap'] + len(params['right_anchor'])
    if width > max_span:
        raise ValueError('maximum pattern span exceeds streaming overlap limit')
    retained, positions, offset, source_offset, count = '', [], 0, 0, 0
    digest = hashlib.sha256()
    for chunk in chunks:
        normalized = normalize_sequence(chunk)
        digest.update(chunk.encode('utf-8'))
        previous_end = offset + len(retained)
        sequence = retained + normalized.sequence
        mapping = positions + [source_offset + p for p in normalized.source_positions]
        source_offset += len(chunk)
        if iupac or strand != '+':
            _dna(sequence)
        left_anchor, right_anchor = params['left_anchor'], params['right_anchor']
        for direction in (('+', '-') if strand == 'both' else (strand,)):
            oriented = sequence if direction == '+' else reverse_complement(sequence)
            length = len(sequence)
            def matches(start, anchor):
                observed = oriented[start:start + len(anchor)]
                return (observed == anchor if not iupac else
                        all(set(_BASES[a]) & set(_BASES[b]) for a, b in zip(observed, anchor)))
            def span(start, end):
                if start == end:
                    return None
                return [start + 1, end] if direction == '+' else [length - end + 1, length - start]
            def source_span(start, end):
                bounds = span(start, end)
                return [mapping[bounds[0] - 1], mapping[bounds[1] - 1]] if bounds else None
            for start in range(length - len(left_anchor) - min_gap - len(right_anchor) + 1):
                if not matches(start, left_anchor):
                    continue
                middle = start + len(left_anchor)
                for gap in range(min_gap, min(params['max_gap'], length - middle - len(right_anchor)) + 1):
                    right_start = middle + gap
                    end = right_start + len(right_anchor)
                    bounds = span(start, end)
                    if offset + bounds[1] <= previous_end or not matches(right_start, right_anchor):
                        continue
                    if count >= max_matches:
                        raise ValueError('match limit exceeded; no complete event will follow')
                    count += 1
                    match = {'strand': direction, 'span': bounds,
                             'left_span': span(start, middle), 'middle_span': span(middle, right_start),
                             'right_span': span(right_start, end), 'source_span': source_span(start, end),
                             'middle_source_span': source_span(middle, right_start),
                             'left': oriented[start:middle], 'middle': oriented[middle:right_start],
                             'right': oriented[right_start:end], 'gap': gap}
                    for name in ('span', 'left_span', 'middle_span', 'right_span'):
                        if match[name]:
                            match[name] = [p + offset for p in match[name]]
                    yield {'schema': 'sequence-match-events/1', 'event': 'match', 'match': match}
        keep = min(width - 1, len(sequence))
        retained = sequence[-keep:] if keep else ''
        positions = mapping[-keep:] if keep else []
        offset += len(sequence) - keep
    yield {'schema': 'sequence-match-events/1', 'event': 'complete', 'match_count': count,
           'sequence_length': offset + len(retained), 'source_length': source_offset,
           'source_sha256': digest.hexdigest(), 'coordinate_convention': '1-based-inclusive',
           'source_coordinate_unit': 'python-string-character', 'order': 'chunk then strand/start/gap'}


@contextmanager
def _open(input):
    if isinstance(input, (str, os.PathLike)):
        path = os.fspath(input)
        _, codec = _extension(path)
        opener = codec[1] if codec else open
        with opener(path, 'rt', encoding='utf-8-sig', newline=None) as handle:
            yield handle
    else:
        yield input


def match_large_file(input, left, right, *, format='fasta', chunk_size=65536, **options):
    """Stream FASTA/text paths (including compression) or caller-owned text streams.

    FASTA headers are limited to chunk_size characters. This explicit mode does
    not sniff formats or accept FASTQ. Line fragments, including unwrapped giant
    sequence lines, are bounded. FASTA source positions refer to concatenated
    sequence-line text (line endings included), excluding headers/comments.
    """
    integer(chunk_size, 'chunk_size', 1)
    if chunk_size > 1048576 or format not in ('fasta', 'fastq', 'text'):
        raise ValueError('large-file matching requires fasta/fastq/text and chunk_size <= 1 Mi')
    with _open(input) as handle:
        if format == "fastq":
            yield from _fastq_events(handle, left, right, chunk_size, options)
            return
        if format == 'text':
            chunks = iter(lambda: handle.read(chunk_size), '')
            for event in match_chunks(chunks, left, right, **options):
                event['record'] = {'id': 'record-1', 'index': 1}
                yield event
            return
        pending, index, line_number = None, 0, 0
        while True:
            header = pending
            pending = None
            if header is None:
                while True:
                    header = handle.readline(chunk_size + 1)
                    if not header:
                        return
                    line_number += 1
                    if line_number == 1:
                        header = header.removeprefix('\ufeff')
                    if header.startswith('>'):
                        break
                    if len(header) > chunk_size:
                        # Consume comment fragments, but reject unexpected sequence.
                        if not header.startswith(';') and header.strip():
                            raise SequenceReadError('sequence before FASTA header', line=line_number)
                        while header and not header.endswith(('\n', '\r')):
                            header = handle.readline(chunk_size)
                    elif header.strip() and not header.startswith(';'):
                        raise SequenceReadError('sequence before FASTA header', line=line_number)
            if len(header) > chunk_size or not header[1:].strip():
                raise SequenceReadError('empty or oversized FASTA header', line=line_number)
            identity = header[1:].strip().split(None, 1)[0]
            index += 1
            def chunks():
                nonlocal pending, line_number
                at_start, comment = True, False
                while True:
                    piece = handle.readline(chunk_size)
                    if not piece:
                        return
                    if at_start:
                        line_number += 1
                        if piece.startswith('>'):
                            if not piece.endswith(('\n', '\r')) and len(piece) == chunk_size:
                                raise SequenceReadError('oversized FASTA header', line=line_number)
                            pending = piece
                            return
                        comment = piece.startswith(';')
                    at_start = piece.endswith('\n')
                    if not comment:
                        yield piece
            for event in match_chunks(chunks(), left, right, **options):
                if event['event'] == 'complete' and event['sequence_length'] == 0:
                    raise SequenceReadError('empty FASTA sequence', line=line_number)
                event['record'] = {'id': identity, 'index': index}
                yield event
            if pending is None:
                return


def _fastq_events(handle, left, right, chunk_size, options):
    index = 0
    while True:
        header = handle.readline(chunk_size + 1)
        if index == 0:
            header = header.removeprefix('\ufeff')
        if not header:
            return
        if not header.strip():
            continue
        if len(header) > chunk_size or not header.startswith('@') or not header[1:].strip():
            raise SequenceReadError('expected bounded nonempty FASTQ @ header')
        identity = header[1:].strip().split(None, 1)[0]
        index += 1
        length = 0
        def chunks():
            nonlocal length
            at_start = True
            while True:
                piece = handle.readline(chunk_size)
                if not piece:
                    raise SequenceReadError('missing FASTQ + separator')
                if at_start and piece.startswith('+'):
                    if not piece.endswith(('\r', '\n')) and len(piece) == chunk_size:
                        raise SequenceReadError('oversized FASTQ + header')
                    repeat = piece[1:].strip()
                    if repeat and repeat != header[1:].strip():
                        raise SequenceReadError('FASTQ + header does not match @ header')
                    if not length:
                        raise SequenceReadError('empty FASTQ sequence')
                    return
                stripped = piece.rstrip('\r\n')
                # A newline-only fragment is valid after a full chunk of sequence.
                if (at_start and not stripped) or any(c.isspace() for c in stripped):
                    raise SequenceReadError('empty/whitespace FASTQ sequence line')
                length += len(stripped)
                at_start = piece.endswith('\n')
                yield piece
        for event in match_chunks(chunks(), left, right, **options):
            event['record'] = {'id': identity, 'index': index}
            if event['event'] != 'complete':
                yield event
                continue
            # Completion is held until the FASTQ quality is validated in full.
            size, at_start = 0, True
            while True:
                piece = handle.readline(chunk_size)
                if not piece:
                    raise SequenceReadError('FASTQ quality-length mismatch')
                stripped = piece.rstrip('\r\n')
                if (at_start and not stripped) or any(not 33 <= ord(c) <= 126 for c in stripped):
                    raise SequenceReadError('invalid FASTQ quality')
                size += len(stripped)
                if size > length:
                    raise SequenceReadError('FASTQ quality-length mismatch')
                at_start = piece.endswith('\n')
                if size == length:
                    if not at_start and len(piece) == chunk_size:
                        # Consume the rest of this line: a newline is permitted, more quality is not.
                        tail = handle.readline(chunk_size)
                        if tail.rstrip('\r\n'):
                            raise SequenceReadError('FASTQ quality-length mismatch')
                    yield event
                    break
