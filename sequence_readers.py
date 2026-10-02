"""Streaming FASTA/FASTQ readers. See docs/SEQUENCE_READERS.md for the subset."""
import bz2
import gzip
import lzma
import os
from pathlib import Path

from sequence_records import ReadProvenance, SequenceRecord

__all__ = ["SequenceReadError", "read_sequences"]
_FORMATS = {".fa": "fasta", ".fas": "fasta", ".fasta": "fasta",
            ".fna": "fasta", ".ffn": "fasta", ".faa": "fasta",
            ".frn": "fasta", ".fq": "fastq", ".fastq": "fastq"}
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
        if format not in ("fasta", "fastq"):
            raise SequenceReadError("format must be 'fasta' or 'fastq'", source)
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


def read_sequences(input, *, format=None, source=None):
    """Yield SequenceRecord from a path or caller-owned text stream.

    Paths are UTF-8 (initial BOM accepted), compression selected by suffix.
    Streams must already be decoded/decompressed; they are never closed here.
    Errors are raised lazily while iterating. Close a partially consumed path
    iterator explicitly to release its file. Memory is proportional to one
    record, not the whole file; no fixed per-record size cap is imposed.
    """
    if isinstance(input, (str, os.PathLike)):
        path = os.fspath(input)
        _, codec = _extension(path)
        compression, opener = codec if codec else (None, open)
        with opener(path, 'rt', encoding='utf-8-sig', newline=None) as handle:
            yield from _read(handle, format, source if source is not None else path, compression, path)
    else:
        name = source if source is not None else getattr(input, 'name', None)
        if not isinstance(name, str):
            name = None
        yield from _read(input, format, name, None, name)


def _read(handle, format, source, compression, filename):
    selected = _detect(handle, format, filename)
    lines = ((number, line.rstrip('\r\n').removeprefix('\ufeff')
              if number == 1 else line.rstrip('\r\n'))
             for number, line in enumerate(handle, 1))
    parser = _fasta if selected == 'fasta' else _fastq
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
