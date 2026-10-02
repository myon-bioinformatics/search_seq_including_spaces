"""Streaming sequence readers. See docs/SEQUENCE_READERS.md for the subset."""
import bz2
import csv
import gzip
import lzma
import os
from pathlib import Path

from sequence_records import ReadProvenance, SequenceRecord

__all__ = ["SequenceReadError", "read_sequences"]
_FORMATS = {".fa": "fasta", ".fas": "fasta", ".fasta": "fasta",
            ".fna": "fasta", ".ffn": "fasta", ".faa": "fasta",
            ".frn": "fasta", ".fq": "fastq", ".fastq": "fastq",
            ".txt": "text", ".csv": "csv", ".tsv": "tsv"}
_CODECS = {".gz": ("gzip", gzip.open), ".bz2": ("bz2", bz2.open),
           ".xz": ("xz", lzma.open), ".lzma": ("lzma", lzma.open)}
_SNIFF_LIMIT = 8192


class SequenceReadError(ValueError):
    """Malformed record or undetermined format, with source and line context."""
    def __init__(self, reason, source=None, line=None):
        self.reason, self.source, self.line = reason, source, line
        super().__init__(f"{source or '<stream>'}: "
                         + (f"line {line}: " if line is not None else "") + reason)


def _extension(name):
    suffixes = Path(name).suffixes if name else []
    codec = _CODECS.get(suffixes[-1].lower()) if suffixes else None
    if codec:
        suffixes.pop()
    return (_FORMATS.get(suffixes[-1].lower()) if suffixes else None), codec


def _detect(handle, format, source):
    if format is not None:
        if format not in ("fasta", "fastq", "text", "csv", "tsv"):
            raise SequenceReadError("format must be fasta, fastq, text, csv or tsv", source)
        return format
    extension, _ = _extension(source)
    if extension:
        return extension
    if not handle.seekable():
        raise SequenceReadError("non-seek stream needs explicit format or recognized filename", source)
    position = handle.tell()
    try:
        # Use the caller's own line splitting, just as the parser does. In
        # particular, Unicode separators are not implicit record boundaries.
        pieces, size = [], 0
        while size < _SNIFF_LIMIT:
            piece = handle.readline(_SNIFF_LIMIT - size)
            if not piece:
                break
            pieces.append(piece)
            size += len(piece)
    finally:
        handle.seek(position)
    sample = ''.join(pieces)
    lines = [piece.rstrip('\r\n') for piece in pieces]
    if lines:
        lines[0] = lines[0].removeprefix('\ufeff')
    fasta_start = 0
    while fasta_start < len(lines) and (not lines[fasta_start].strip()
                                       or lines[fasta_start].startswith(';')):
        fasta_start += 1
    if (fasta_start < len(lines) and lines[fasta_start].startswith('>')
            and lines[fasta_start][1:].strip()):
        return "fasta"
    # FASTQ permits blank record boundaries, but not FASTA-style comments.
    while lines and not lines[0].strip():
        lines.pop(0)
    # Only a complete simple four-line record is evidence for FASTQ sniffing.
    # Wrapped/long/ambiguous input needs an explicit format or extension.
    if (len(lines) >= 4 and lines[0].startswith('@') and lines[0][1:].strip()
            and lines[1] and not any(c.isspace() for c in lines[1])
            and lines[2].startswith('+') and len(lines[1]) == len(lines[3])
            and all(33 <= ord(c) <= 126 for c in lines[3])
            and (sample.endswith(('\n', '\r')) or len(sample) < _SNIFF_LIMIT
                 or len(lines) > 4)):
        return "fastq"
    raise SequenceReadError("cannot conservatively detect format; specify format", source)


def read_sequences(input, *, format=None, source=None, sequence_column=None,
                   id_column=None, name_column=None, description_column=None):
    """Yield SequenceRecord from a path or caller-owned text stream.

    Paths are UTF-8 (initial BOM accepted), compression selected by suffix.
    Streams must already be decoded/decompressed; they are never closed here.
    Errors are raised lazily while iterating. Close a partially consumed path
    iterator explicitly to release its file. Memory is proportional to one
    record, not the whole file. CSV/TSV use csv.field_size_limit(); other
    formats have no fixed per-record cap. Tabular column names are explicit;
    id_column and name_column are aliases and cannot be combined.
    """
    columns = (sequence_column, id_column, name_column, description_column)
    if isinstance(input, (str, os.PathLike)):
        path = os.fspath(input)
        _, codec = _extension(path)
        compression, opener = codec if codec else (None, open)
        with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
            yield from _read(handle, format, source if source is not None else path, compression, path, columns)
    else:
        name = source if source is not None else getattr(input, 'name', None)
        if not isinstance(name, str):
            name = None
        yield from _read(input, format, name, None, name, columns)


def _read(handle, format, source, compression, filename, columns):
    selected = _detect(handle, format, filename)
    if selected in ('csv', 'tsv'):
        yield from _tabular(handle, source, compression, selected, *columns)
        return
    if any(column is not None for column in columns):
        raise SequenceReadError("column options require csv or tsv", source)
    lines = ((number, line.rstrip('\r\n').removeprefix('\ufeff')
              if number == 1 else line.rstrip('\r\n'))
             for number, line in enumerate(handle, 1))
    parser = {'fasta': _fasta, 'fastq': _fastq, 'text': _text}[selected]
    yield from parser(iter(lines), source, compression)


def _record(header, chunks, quality, source, format, compression, index, start, end):
    parts = header.strip().split(None, 1)
    return SequenceRecord(parts[0], parts[1] if len(parts) > 1 else '',
                          ''.join(chunks), source,
                          ReadProvenance(format, compression, index, start, end), quality)


def _fasta(lines, source, compression):
    header, chunks, start, end, index = None, [], 0, 0, 0
    for number, line in lines:
        if line.startswith('>'):
            if header is not None:
                if not chunks:
                    raise SequenceReadError("empty FASTA sequence", source, start)
                index += 1
                yield _record(header, chunks, None, source, 'fasta', compression, index, start, end)
            header, chunks, start = line[1:], [], number
            if not header.strip():
                raise SequenceReadError("empty FASTA header", source, number)
        elif line.strip() and not line.startswith(';'):
            if header is None:
                raise SequenceReadError("sequence before FASTA header", source, number)
            chunks.append(''.join(line.split()))
        end = number
    if header is not None:
        if not chunks:
            raise SequenceReadError("empty FASTA sequence", source, start)
        yield _record(header, chunks, None, source, 'fasta', compression, index + 1, start, end)


def _fastq(lines, source, compression):
    index = 0
    for start, line in lines:
        if not line.strip():
            continue
        if not line.startswith('@') or not line[1:].strip():
            raise SequenceReadError("expected nonempty FASTQ @ header", source, start)
        header, chunks, length = line[1:], [], 0
        for number, line in lines:
            if line.startswith('+'):
                if line[1:] and line[1:].strip() != header.strip():
                    raise SequenceReadError("FASTQ + header does not match @ header", source, number)
                break
            if not line or any(c.isspace() for c in line):
                raise SequenceReadError("empty/whitespace FASTQ sequence line", source, number)
            chunks.append(line)
            length += len(line)
        else:
            raise SequenceReadError("missing FASTQ + separator", source, start)
        if not length:
            raise SequenceReadError("empty FASTQ sequence", source, start)
        quality, size = [], 0
        for number, line in lines:
            if not line or any(not 33 <= ord(c) <= 126 for c in line):
                raise SequenceReadError("invalid FASTQ quality (expected ASCII 33..126)", source, number)
            quality.append(line)
            size += len(line)
            if size >= length:
                break
        if size != length:
            raise SequenceReadError("FASTQ quality-length mismatch", source, start)
        index += 1
        yield _record(header, chunks, ''.join(quality), source, 'fastq', compression, index, start, number)


def _text(lines, source, compression):
    chunks, start, end = [], None, 0
    for number, line in lines:
        sequence = ''.join(line.split())
        if sequence:
            if start is None:
                start = number
            chunks.append(sequence)
        end = number
    if start is not None:
        yield SequenceRecord('record-1', '', ''.join(chunks), source,
                             ReadProvenance('text', compression, 1, start, end))


def _tabular(handle, source, compression, format, sequence_column, id_column,
             name_column, description_column):
    if not isinstance(sequence_column, str) or not sequence_column:
        raise SequenceReadError("csv/tsv requires a nonempty sequence_column name", source)
    if id_column is not None and name_column is not None:
        raise SequenceReadError("select id_column or name_column, not both", source)
    identity = id_column if id_column is not None else name_column
    for column in (identity, description_column):
        if column is not None and (not isinstance(column, str) or not column):
            raise SequenceReadError("column names must be nonempty strings", source)
    # Keep terminators for csv.reader: quoted fields may contain newlines.
    lines = (line.removeprefix('\ufeff') if number == 1 else line
             for number, line in enumerate(handle, 1))
    reader = csv.reader(lines, delimiter=',' if format == 'csv' else '\t', strict=True)
    try:
        header = next(reader, None)
        if header is None:
            return
        if len(set(header)) != len(header) or any(not name for name in header):
            raise SequenceReadError("empty or duplicate column name", source, 1)
        requested = (sequence_column, identity, description_column)
        for column in requested:
            if column is not None and column not in header:
                raise SequenceReadError(f"missing column {column!r}", source, 1)
        indices = [header.index(column) if column is not None else None for column in requested]
        index = 0
        while True:
            start = reader.line_num + 1
            row = next(reader, None)
            if row is None:
                break
            if not row:  # blank physical rows are boundaries, not empty sequences
                continue
            if len(row) != len(header):
                raise SequenceReadError("row column count does not match header", source, start)
            sequence = ''.join(row[indices[0]].split())
            if not sequence:
                raise SequenceReadError("empty tabular sequence", source, start)
            index += 1
            identifier = row[indices[1]] if indices[1] is not None else f'record-{index}'
            if not identifier.strip():
                raise SequenceReadError("empty tabular id/name", source, start)
            description = row[indices[2]] if indices[2] is not None else ''
            yield SequenceRecord(identifier, description, sequence, source,
                                 ReadProvenance(format, compression, index, start, reader.line_num))
    except csv.Error as error:
        raise SequenceReadError(f"malformed {format}: {error}", source, reader.line_num) from error
