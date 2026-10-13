"""JSON/JSONL and GenBank ORIGIN subsets, shared SequenceRecord contract."""
import json
import re

from sequence_records import ReadProvenance, SequenceRecord

MAX_JSON_CHARS = 16 * 1024 * 1024


def strict_json(text):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f'duplicate JSON key {key!r}')
            result[key] = value
        return result
    def constant(value):
        raise ValueError(f'non-finite JSON number {value}')
    return json.loads(text, object_pairs_hook=pairs, parse_constant=constant)


def _path(value, path):
    if not isinstance(path, str) or not path:
        raise ValueError('JSON paths must be a nonempty key or RFC 6901 pointer')
    parts = path[1:].split('/') if path.startswith('/') else [path]
    for part in parts:
        if path.startswith('/'):
            if re.search(r'~(?![01])', part):
                raise ValueError('invalid JSON pointer escape')
            part = part.replace('~1', '/').replace('~0', '~')
        if isinstance(value, dict):
            value = value[part]
        elif isinstance(value, list) and re.fullmatch(r'0|[1-9][0-9]*', part):
            value = value[int(part)]
        else:
            raise ValueError(f'cannot resolve JSON path {path!r}')
    return value


def read_json(handle, source, compression, format, sequence_path, id_path, description_path):
    from sequence_readers import SequenceReadError
    sequence_path = 'sequence' if sequence_path is None else sequence_path
    for path in (sequence_path, id_path, description_path):
        if path is not None and (not isinstance(path, str) or not path):
            raise SequenceReadError('JSON paths must be nonempty strings', source)
    if format == 'jsonl':
        def objects():
            number = 0
            while True:
                line = handle.readline(MAX_JSON_CHARS + 1)
                if not line:
                    return
                number += 1
                if len(line) > MAX_JSON_CHARS:
                    raise ValueError('JSONL record exceeds 16 Mi characters')
                if number == 1:
                    line = line.removeprefix('\ufeff')
                if line.strip():
                    try:
                        value = strict_json(line)
                    except (ValueError, RecursionError) as error:
                        raise SequenceReadError(f"invalid jsonl: {error}", source, number) from None
                    yield value, number, number
    else:
        def objects():
            # JSON is deliberately bounded and document-based; JSONL is streaming.
            text = handle.read(MAX_JSON_CHARS + 1).removeprefix('\ufeff')
            if len(text) > MAX_JSON_CHARS:
                raise ValueError('JSON document exceeds 16 Mi characters; use JSONL')
            data = strict_json(text)
            records = data if isinstance(data, list) else [data]
            # All records share the document's bounding line span, not invented exact ranges.
            end = max(1, len(text.splitlines()))
            for value in records:
                yield value, 1, end
    index = 0
    try:
        for value, start, end in objects():
            index += 1
            if not isinstance(value, dict):
                raise ValueError('each JSON record must be an object')
            sequence = _path(value, sequence_path)
            if not isinstance(sequence, str) or not sequence.strip():
                raise ValueError('sequence path must select nonempty text')
            identity = _path(value, id_path) if id_path is not None else f'record-{index}'
            description = _path(value, description_path) if description_path is not None else ''
            if not isinstance(identity, str) or not identity.strip() or not isinstance(description, str):
                raise ValueError('id must be nonempty text and description must be text')
            yield SequenceRecord(identity, description, sequence, source,
                                 ReadProvenance(format, compression, index, start, end))
    except SequenceReadError:
        raise
    except (ValueError, KeyError, IndexError, RecursionError) as error:
        raise SequenceReadError(f'invalid {format} record {index or 1}: {error}', source) from None


def read_genbank(lines, source, compression):
    """GenBank/GenPept LOCUS..ORIGIN..// subset. No FEATURES/CONTIG semantics."""
    from sequence_readers import SequenceReadError
    index, record, origin, definition = 0, None, False, False
    for number, line in lines:
        if line.startswith('LOCUS '):
            if record is not None:
                raise SequenceReadError('unterminated GenBank record', source, number)
            parts = line.split()
            if len(parts) < 2:
                raise SequenceReadError('empty LOCUS', source, number)
            record = {'locus': parts[1], 'description': '', 'chunks': [], 'start': number}
            origin, definition = False, False
        elif record is None:
            if line.strip():
                raise SequenceReadError('expected GenBank LOCUS', source, number)
        elif line.strip() == '//':
            if not origin or not record['chunks']:
                raise SequenceReadError('GenBank record requires nonempty ORIGIN (CONTIG unsupported)', source, number)
            index += 1
            metadata = {k: record[k] for k in ('locus', 'accession', 'version') if k in record}
            yield SequenceRecord(record.get('version', record.get('accession', record['locus'])),
                                 record['description'], ''.join(record['chunks']), source,
                                 ReadProvenance('genbank', compression, index, record['start'], number),
                                 metadata=metadata)
            record, origin = None, False
        elif origin:
            if not re.fullmatch(r'\s*\d+\s+[A-Za-z\s]+', line):
                raise SequenceReadError('invalid ORIGIN sequence line', source, number)
            position, *parts = line.split()
            if int(position) != sum(map(len, record['chunks'])) + 1:
                raise SequenceReadError('nonconsecutive ORIGIN position', source, number)
            record['chunks'].append(''.join(parts))
        elif line.startswith('ORIGIN'):
            origin, definition = True, False
        elif line.startswith('DEFINITION '):
            record['description'], definition = line[12:].strip(), True
        elif definition and line.startswith(' ' * 12):
            record['description'] += ' ' + line.strip()
        else:
            definition = False
            if line.startswith(('ACCESSION ', 'VERSION ')):
                parts = line.split()
                if len(parts) < 2:
                    raise SequenceReadError('empty accession/version', source, number)
                record[parts[0].lower()] = parts[1]
    if record is not None:
        raise SequenceReadError('unterminated GenBank record (missing //)', source, record['start'])
