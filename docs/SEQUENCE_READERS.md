# Streaming sequence records (roadmap #2, FASTA/FASTQ foundation and text/tabular increment)

`sequence_records.py` defines frozen stdlib dataclasses. `dataclasses.asdict()`
produces JSON-friendly fields; serialized envelopes may use schema
`sequence-record/1`. There is no network, analysis or matcher dependency.

| Field | Type / semantics |
| --- | --- |
| id | str, nonempty; derived per format (FASTA/FASTQ: first whitespace-delimited header token) |
| description | str; FASTA/FASTQ remaining header text, or empty string |
| sequence | str; sequence lines concatenated, case/symbols preserved |
| quality | str for FASTQ, None for other formats; FASTQ requires equal length to sequence; no Phred decoding |
| source | str input path/caller label, or None if unknown |
| provenance | ReadProvenance: format, compression (gzip/bz2/xz/lzma or None), record_index, start_line, end_line |

Record index and line ranges are 1-based inclusive in **decompressed text**,
relative to the stream's current position. They are not file byte offsets or
genomic coordinates. A range is a **bounding span**, not a claim that every
line in it contributes residues or belongs exclusively to that record. FASTA
ranges include trailing/intervening blank/comment lines up to the next header;
FASTQ ranges end at the last quality line and exclude boundary blanks.
`anchor_matcher` spans are also 1-based inclusive: normalized spans refer to
the normalized sequence, while `source_span` / `middle_source_span` refer to
Python character positions in `record.sequence`, never the original file.
No whole-file hash or per-residue file-coordinate map is promised. FASTA removes
sequence whitespace; FASTQ removes only line endings and rejects whitespace in
sequence lines. Neither reader uppercases, removes N/gaps/stops, infers alphabet,
or validates DNA/RNA/protein symbols. Those are explicit downstream policies.

### Contract evolution

The current `ReadProvenance.format` vocabulary is lowercase `fasta` / `fastq` / `text` / `csv` / `tsv`.
Each future reader must document its new canonical format name and how it derives
`id` and `description`; the FASTA header rule is not universal. New trailing
fields with defaults and documented new format names may be added within
`sequence-record/1`, preserving existing positional construction and existing
field meanings. Consumers should tolerate additional serialized keys. Removing,
renaming, reordering or changing existing field meanings requires a new schema
version. There is no generic metadata field in this increment; future formats
must define a typed/defaulted extension before emitting extra metadata, rather
than hiding structured data in `description`.

Current readers always provide meaningful positive line ranges; there is no
unknown-coordinate sentinel. A future non-contiguous record may use a documented
bounding span. A format without meaningful line coordinates needs an explicit
contract decision before implementation.

## Existing implementation map

| Existing surface | Contract / relationship |
| --- | --- |
| anchor_matcher.NormalizedSequence | Uppercase letters plus original positions within the supplied sequence text; remains unchanged |
| anchor_matcher.find_matches | `anchor-gap/1` dict; pass `record.sequence` explicitly; source spans refer to that concatenated string, not the file |
| protein.parse.read_fasta | Existing list of `(header, raw)` tuples, including plain-text fallback; unchanged for CLI compatibility |
| protein.parse.Prepared / protein features | Existing protein validation, transform positions and outputs; later adapter may consume SequenceRecord |
| sequence_tool.py | Existing CLI and input behavior unchanged; this increment adds a library API only |
| definition_main.py | Legacy script and CSV behavior unchanged |

## Supported input subset

```python
from sequence_readers import read_sequences

for record in read_sequences('input.fasta.gz'):
    print(record.id, len(record.sequence))

# Caller owns stream lifetime; decoded/decompressed text streams only.
for record in read_sequences(text_stream, format='fastq', source='stdin'):
    print(record.id, record.quality)
```

- FASTA: multiline sequence, nonempty headers/sequences, LF/CRLF, initial BOM,
  blank lines and column-one `;` comments. Headers must start in column one.
- FASTQ: multiline sequence and quality, nonempty `@` header/sequence, `+`
  separator (optional repeated complete header must match), ASCII quality
  33..126 and exact sequence/quality length. Quality lines starting `@` or `+`
  are quality, not delimiters. Short/long quality, missing separator, whitespace
  sequence and invalid quality raise `SequenceReadError` with source/line.
  Blank/whitespace-only lines before, between and after complete records are
  ignored. Blank lines inside sequence or quality remain errors. Column-one `;`
  comments are FASTA-only and are rejected at FASTQ record boundaries.
  Wrapped quality consumes lines until sequence length; missing data can consume
  a subsequent header before reporting an error, so there is no recovery/resync.
- Paths: UTF-8; transparent `.gz`, `.bz2`, `.xz`, `.lzma` selected by suffix.
  One compression suffix is stripped before extension detection. No magic-byte
  detection or nested compression. Stream compression is caller-managed and
  reported as None (unknown/unmanaged), even if its label ends in `.gz`.
- Empty explicitly selected input yields no records. OSError, codec errors and
  UnicodeDecodeError propagate; record/detection errors use SequenceReadError.

Format selection is **explicit `format` → extension → conservative sniffing**.
FASTA extensions: fa/fas/fasta/fna/ffn/faa/frn; FASTQ: fq/fastq;
plain text: txt; tabular: csv/tsv (case insensitive).
For paths the physical filename selects the format even when `source` is a
logical provenance label. For streams the source/name supplies the extension.
An explicit unsupported format fails rather than falling through.

Sniffing reads at most 8192 decoded characters on a seekable stream, restores
its current position, and recognizes a nonempty FASTA header or a complete
simple four-line FASTQ record with matching quality length. It does not guess
plain sequences, wrapped/long FASTQ, or ambiguous/truncated samples. Use explicit
format for those cases. **Non-seek streams require explicit format or a recognized
filename; no sniff read or seek occurs.** No stdin/binary-stream CLI adapter is
added here.

Sniffing uses bounded `readline()` calls and the parser uses iteration over the
same handle, preserving its line-splitting policy; Unicode separators such as
U+0085/U+2028 are not separately split by the sniffer. Paths recognize universal
newlines (LF, CRLF and lone CR), retaining terminators for CSV quoted fields. Caller-owned text streams retain their configured
newline policy, so use universal-newline decoding for path-equivalent line
numbers/results. Leading blank lines are allowed for either format; leading `;`
comments are skipped only for FASTA. A complete quality line ending exactly at
character 8192 is sniffable; a quality line truncated at that limit, or ending
there without a line terminator (EOF cannot be established within the limit),
requires explicit format/extension.

Parsing yields records lazily with memory proportional to the largest current
record/line, not total file size. A single chromosome-sized record still occupies
memory; there is no FASTA/FASTQ/text per-record cap or constant-memory chromosome claim.
CSV/TSV retain the stdlib `csv.field_size_limit()` (usually 131072 characters);
oversized fields raise contextual SequenceReadError. This library does not mutate
that process-global limit. Callers may explicitly configure it for large fields. Iteration
can yield valid records before a later error; exhaust it to validate the whole
input. Paths are closed on exhaustion/error or generator.close(); close partially
consumed path iterators explicitly. Caller-provided streams remain open.

## Follow-ups (not complete in this PR)

JSON/JSONL (key/path selection), GenBank best-effort and broader exports require
separate normal/error fixtures and documented contracts. No best-effort format is currently claimed. The protein/CLI migration,
custom pattern CLI and feature result schema are separate increments. Roadmap
#2's multi-format reader checkbox remains open until those requirements and
main/test evidence are satisfied.


## Plain text and CSV/TSV increment

```python
# A whole plain-text file is one sequence; record id is generated.
record, = read_sequences('sequence.txt.gz')

# Header names are exact/case-sensitive, never guessed.
records = read_sequences('table.csv.bz2', sequence_column='residues',
                         id_column='accession')
# name_column is an alias for selecting the identity, not a description.
records = read_sequences(text_stream, format='tsv', sequence_column='sequence',
                         name_column='name', description_column='notes')
```

- `text`: concatenate all nonempty sequence lines after removing whitespace,
  preserving case, Unicode and symbols. No header/comment syntax, alphabet
  detection or record splitting. ID is `record-1`, description is empty, quality
  is None. Empty/whitespace-only input yields no records. Span starts at the first
  non-whitespace line and ends at the last physical line, including trailing
  blank lines; internal blank lines are part of this bounding span.
- `csv` / `tsv`: first logical row is the required header. `sequence_column` is
  mandatory even with a recognized extension. Optional `id_column` or
  `name_column` (mutually exclusive) selects a nonempty identity; otherwise IDs
  are `record-1`, `record-2`, etc. Optional `description_column` copies a field
  verbatim, otherwise description is empty. Identity values are copied verbatim
  but must contain a non-whitespace character. Quality is None.
- Delimiter is comma for CSV and tab for TSV, with stdlib double-quote escaping
  and strict CSV parsing. Duplicate/empty header names, missing selected columns,
  short/extra-field rows, whitespace-only sequences and empty selected identity
  fields are errors. Sequence whitespace, including quoted newlines, is removed;
  other characters and case are preserved. Blank physical rows outside quotes
  are skipped; a quoted empty sequence is an error. Empty input yields no records
  after validating the required column options.
- A tabular span covers the physical lines consumed for that logical row,
  including quoted multiline fields. It excludes the header and skipped boundary
  blanks. Errors include source and physical line context; CSV syntax errors
  report the last physical line consumed. Initial BOM is accepted.
- Use `newline=''` when opening caller-owned CSV/TSV streams, so quoted CR/LF
  terminators survive and line counts match physical input. Path readers already
  do this. All existing compression suffixes also work for txt/csv/tsv.
- Text/CSV/TSV are **never sniffed**: arbitrary words and tables are ambiguous.
  Choose an explicit format or recognized extension. Existing conservative
  FASTA/FASTQ sniffing is unchanged. Column options on other formats are errors.

## Deferred preprocessing / export follow-up

This increment adds readers only. `sequence_conversion.py`, `fasta_to_text()`
and `text_to_fasta()` are deferred to roadmap #2's consolidated former #8
preprocessing / multi-format export track. They are output consumers of
`read_sequences()`, not dependencies of the reader API or `SequenceRecord`
contract, so FASTA ↔ TSV conversion is not needed to deliver text/CSV/TSV input.

A separate small export PR should define its supported representation, validation,
semantic round-trip guarantees (not original bytes/layout), streaming failure and
output ownership policy, and machine-readable transformation provenance/receipts
required by former #8. Its tests should cover compressed inputs, empty inputs,
Unicode/quoted fields, parsed field/order preservation and unrepresentable FASTA
fields. The removed prototype and its tests remain available in PR #25's original
commit `cf2750ec3e054ec1b6c956df0ab97b1f593ab5f8` as a reference, not an API promised
by this reader increment. Former #8 stays closed as a tracking consolidation;
implementation progress is tracked in #2.

The foundation was merged in #24 (`c9d5878a9c4601d38d654f836f9d474c92370709`).
This follow-up adds text/tabular readers only; JSON/JSONL, GenBank and
preprocessing/export remain roadmap work.
