# Record-by-record matching (roadmap #2, Phase 4 adapter)

`sequence_matcher.find_sequence_matches()` connects `read_sequences()` to
`anchor_matcher.find_matches()` without changing either core contract. It is a
stdlib-only Python library API, not a CLI or a new pattern language.

```python
from contextlib import closing
from sequence_matcher import find_sequence_matches

with closing(find_sequence_matches(
        'input.fasta.gz', 'AC', 'GT', min_gap=1, max_matches=1000)) as results:
    for envelope in results:
        print(envelope['record']['id'], envelope['record']['provenance'])
        for match in envelope['result']['matches']:
            print(match['source_span'], match['middle'])
```

For text/tabular inputs the existing reader options apply:

```python
with closing(find_sequence_matches(
        'table.csv.bz2', 'AC', 'GT', min_gap=1,
        sequence_column='residues', id_column='accession',
        description_column='notes')) as results:
    for envelope in results:
        print(envelope['record']['id'], envelope['result']['matches'])
```

`input`, `format`, `source`, `sequence_column`, `id_column`, `name_column`, and
`description_column` are passed to the existing reader. Format selection,
compression, column validation and supported subsets remain as documented in
[SEQUENCE_READERS.md](SEQUENCE_READERS.md). `left`, `right`, `min_gap`, `max_gap`,
`iupac`, `strand`, and `max_matches` are passed unchanged to the existing matcher.
Its defaults and validation remain unchanged; `max_matches` defaults to 10000.

## Envelope and coordinates

Each input record yields one independent JSON-friendly dictionary, in input
order, **including records with zero matches**. Records are never concatenated,
so anchors cannot span records. Duplicate IDs are retained, not merged; use the
source and reader `record_index` to distinguish them.

```text
schema: "sequence-matches/1"
record:
  id: original record id
  description: original description
  source: input path / caller label / null
  provenance:
    format, compression, record_index, start_line, end_line
result:
  the complete, unmodified return value of find_matches(record.sequence, ...)
  (schema "anchor-gap/1")
```

`record` is a metadata projection, **not** a complete `sequence-record/1`
serialization: raw `sequence` and FASTQ `quality` are not exported. Quality is
validated by the reader but is not searched. Reader and matcher provenance are
separate dictionaries; no fields are injected into the `anchor-gap/1` result.

The nested `result` retains its 1-based inclusive coordinate semantics:
normalized spans refer to the normalized sequence; `source_span` and
`middle_source_span` count Python characters in **`record.sequence`**, not the
original file, its header, or encoded/compressed bytes. Reverse-strand matches
still use forward-reference coordinates. Reader whitespace removal has already
happened before matching. Reader `start_line` / `end_line` are decompressed-text
bounding line spans, not per-residue file-coordinate mappings. The matcher's
`source_sha256` hashes `record.sequence`, not the file or record header.

## Laziness, errors and ownership

Calling the function constructs a generator; reading and matching begin when
it is advanced. The adapter requests one record and matches it before requesting
the next. It does not listify, prefetch or retain earlier envelopes. The reader
may still perform its documented bounded sniffing, FASTA next-header lookahead,
and ordinary I/O buffering; this is not a byte-exact no-read-ahead promise.

`max_matches` is **per record**, not a whole-input budget. Exactly the limit is
allowed. Overflow raises the existing `ValueError` instead of silently
truncating, and no envelope for the offending record is yielded. Earlier
complete records may already have been yielded. Matcher option validation is
also per record; an empty reader yields nothing without calling the matcher.

Reader errors (including source/line context), matcher errors, I/O errors,
codec errors and decoding errors propagate unchanged and terminate iteration.
There is no skip/resync/recovery policy. Nested `result.complete` only describes
that one record; there is no outer whole-input success flag. Exhaust the
iterator successfully to establish that all input was processed.

The adapter closes its reader on exhaustion, error or explicit generator
`close()`. Use `contextlib.closing` when stopping early or when consumer code
can raise: `break` alone does not close a generator still referenced by the
caller. Reader-owned path handles are released; caller-owned decoded,
decompressed text streams are never closed. A stream's compression remains
`None` as in the reader, even when its source label ends in `.gz`.

## Memory guarantee and remaining scope

Memory is proportional to **the current one record plus its matches**, including
reader/codec buffers, FASTQ quality, normalized sequence/position maps and
matcher temporary workspace. Neither records nor results from previous records
are accumulated by this adapter. Consumer code that saves envelopes (for
example, `list(results)`) adds its own storage and loses that streaming property.

This is **not constant memory for a giant single FASTA record**. Plain text is
also one whole record. `find_matches()` materializes the matches for each record;
`max_matches` bounds their count, not the size of captured middle strings or the
size of the sequence. There is no new per-record byte cap, genome-scale indexing,
or incremental matching within a record. Phase 4 therefore remains partial.

`anchor_matcher.py`, `sequence_readers.py`, `sequence_records.py`,
`definition_main.py`, and `sequence_tool.py` behavior is unchanged. Phase 2
custom-pattern CLI/schema, JSON/JSONL readers, GenBank, exporters and biological
presets are separate work and are not implemented by this adapter.
