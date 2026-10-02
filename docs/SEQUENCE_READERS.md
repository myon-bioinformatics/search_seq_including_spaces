# Streaming sequence records (roadmap #2, first readers increment)

`sequence_records.py` defines frozen stdlib dataclasses. `dataclasses.asdict()`
produces JSON-friendly fields; serialized envelopes may use schema
`sequence-record/1`. There is no network, analysis or matcher dependency.

| Field | Type / semantics |
| --- | --- |
| id | str; first whitespace-delimited header token, nonempty |
| description | str; remaining header text, or empty string |
| sequence | str; sequence lines concatenated, case/symbols preserved |
| quality | str for FASTQ, None for FASTA; equal length to sequence; no Phred decoding |
| source | str input path/caller label, or None if unknown |
| provenance | ReadProvenance: format, compression (gzip/bz2/xz/lzma or None), record_index, start_line, end_line |

Record index and line ranges are 1-based inclusive in **decompressed text**,
relative to the stream's current position. They are not file byte offsets or
genomic coordinates. FASTA ranges include intervening blank/comment lines.
No whole-file hash or per-residue file-coordinate map is promised. FASTA removes
sequence whitespace; FASTQ removes only line endings and rejects whitespace in
sequence lines. Neither reader uppercases, removes N/gaps/stops, infers alphabet,
or validates DNA/RNA/protein symbols. Those are explicit downstream policies.

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
  Wrapped quality consumes lines until sequence length; missing data can consume
  a subsequent header before reporting an error, so there is no recovery/resync.
- Paths: UTF-8; transparent `.gz`, `.bz2`, `.xz`, `.lzma` selected by suffix.
  One compression suffix is stripped before extension detection. No magic-byte
  detection or nested compression. Stream compression is caller-managed and
  reported as None (unknown/unmanaged), even if its label ends in `.gz`.
- Empty explicitly selected input yields no records. OSError, codec errors and
  UnicodeDecodeError propagate; record/detection errors use SequenceReadError.

Format selection is **explicit `format` → extension → conservative sniffing**.
FASTA extensions: fa/fas/fasta/fna/ffn/faa/frn; FASTQ: fq/fastq (case insensitive).
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

Parsing yields records lazily with memory proportional to the largest current
record/line, not total file size. A single chromosome-sized record still occupies
memory; there is no per-record cap or constant-memory chromosome claim. Iteration
can yield valid records before a later error; exhaust it to validate the whole
input. Paths are closed on exhaustion/error or generator.close(); close partially
consumed path iterators explicitly. Caller-provided streams remain open.

## Follow-ups (not complete in this PR)

Plain text, CSV/TSV (column selection), JSON/JSONL (key/path selection), GenBank
best-effort and exports require separate normal/error fixtures and documented
contracts. No best-effort format is currently claimed. The protein/CLI migration,
custom pattern CLI and feature result schema are separate increments. Roadmap
#2's multi-format reader checkbox remains open until those requirements and
main/test evidence are satisfied.
