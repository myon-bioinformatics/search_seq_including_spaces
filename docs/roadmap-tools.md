# Sequence roadmap tools

This increment addresses open issues #29 and #30 and implements the executable
baseline across roadmap #2. Runtime remains Python 3.10+ standard library only;
all existing protein golden files and the legacy `definition_main.py` are retained.

## Delivery map

| Track | Interface | Verification / boundary |
| --- | --- | --- |
| #2 Phase 2 | `match`, `sequence_patterns.load_patterns` | strict versioned pattern packs, CLI JSON/JSONL, reader provenance |
| #2 Phase 3 | `preset`, `sequence_analysis.marker_clues` | explicit poly-A, tandem repeat, SNP allele and repeated-anchor TSD clues |
| #2 Phase 4 / former #4 | extended `read_sequences`, `match --chunk-size` | JSON/JSONL, GenBank ORIGIN subset; bounded overlap FASTA/FASTQ/text search |
| Former #5 | `fetch`, `acquisition.Client` | NCBI search/page/summary/sequence/PubMed; injected offline transports |
| Former #6 | `nucleotide-features`, `nucleotide_pairwise`, `protein-ionization` | composition/entropy/k-mers/codons, transition/transversion, pH charge/pI |
| Former #7 | `rna-fold`, `fold_rna` | bounded Nussinov maximum-pair baseline, deterministic dot-bracket and loop descriptions |
| Former #8 | `preprocess`, `export_records` | filtering, names, deduplication and transform receipts; FASTA/JSONL/MEGA |
| Former #9 | `structure`, `Client.structure` | bounded coordinate subsets, spatial-bin contacts; RCSB acquisition / MMDB lookup link |
| #2 Phases 5–6 | `sequence_service.py`, `examples/service` | HTTP/webhook call endpoint, MCP stdio, optional Docker definition |
| #2 Phase 7 / former #10 | `sequence_router`, `jev_adapter`, offline orchestration example | deterministic-first routing, official Jev Choice wire contract, explicit thresholds, failure fallback |
| #29 | `summarize-hits`, `compare --window ... --summary` | bounded S3 output with counts, provenance and source-row references |
| #30 | `candidates`, `select_representatives` | pairwise is a **separate fifth mode**; representative selection is an independent opt-in |

Implementation and fixture validation do not constitute deployment or scientific
validation. Live NCBI/RCSB/Jev, model quality/calibration, Docker execution, MEGA GUI
import and a running llama.cpp/Pro backend are separate operational checks. The
parent roadmap should remain open for these gates and advanced optional research
rather than treating a merged prototype as proof that every future goal is done.

## CLI examples

```sh
python -S sequence_tool.py match input.fa --left AC --right GT --min-gap 1 --jsonl
python -S sequence_tool.py match input.jsonl --sequence-path /payload/sequence --id-path id --patterns examples/patterns.json --jsonl
python -S sequence_tool.py match genome.fa.gz --input-format fasta --chunk-size 65536 --left AC --right GT --min-gap 1 --jsonl
python -S sequence_tool.py nucleotide-features input.fa --alphabet dna --k 3 --codon-frame 0 --jsonl
python -S sequence_tool.py translate input.fa --frame 0 --strand + --jsonl
python -S sequence_tool.py preset input.fa --kind poly-a --min-length 6 --jsonl
python -S sequence_tool.py preset input.fa --kind tsd --repeat-length 4 --min-gap 5 --max-gap 50 --jsonl
python -S sequence_tool.py rna-fold input.fa --min-loop 3 --jsonl
python -S sequence_tool.py protein-ionization input.fa --ph 7 --jsonl
python -S sequence_tool.py preprocess input.fa --alphabet dna --deduplicate sequence --aligned --export-format mega --receipt new-receipt.json > output.meg
python -S sequence_tool.py structure input.cif --input-format cif --cutoff 4 --top 20
python -S sequence_tool.py fetch --operation search --database pubmed --query 'RNA folding' --limit 20
python -S sequence_tool.py fetch --operation records --database structure --ids 1CRN --format cif
python -S sequence_tool.py fetch --operation records --database structure --ids 1CRN --format pdb
```

New record commands use `--input-format` plus the common reader's column/path
options. `--json` (default) requires exactly one result, including one pattern;
`--jsonl` emits results incrementally. A later input error exits 2 and may leave
previously emitted JSONL lines; a prefix does **not** prove whole-input success.
Malformed input, bad option combinations, limits and I/O errors exit 2 without
a Python traceback. The legacy script and old protein output schemas are intact.

A pattern pack is `{"schema":"sequence-patterns/1","patterns":[...]}`. Each
pattern requires unique nonempty `id`, `left`, `right`; optional fields are
`min_gap`, `max_gap`, `iupac`, `strand`, `max_matches`. Unknown keys, duplicate JSON
keys, nonfinite numbers, booleans as integers and invalid matcher options fail.
At most 100 patterns / 1 Mi characters of pattern JSON are accepted. Pattern
files and inline anchors are mutually exclusive. CLI `--max-matches` caps each
record/pattern (also caps a larger limit supplied in a pack).

## Input contracts and giant records

`SequenceRecord` gains a trailing defaulted `metadata: dict[str, str]`. Existing
positional fields and schema remain compatible. GenBank uses it for `locus`,
`accession`, `version`; it does not encode metadata inside a description.

JSON accepts one object or an array of objects, with a 16 Mi-character document
cap. All records share the document's enclosing line span. JSONL streams nonblank
lines with a 16 Mi-character per-line cap and exact line provenance. `sequence_path`
defaults to the literal key `sequence`; `id_path` and `description_path` are opt-in.
A path is a literal key or an RFC 6901 pointer (`/payload/sequence`, `/items/0/seq`;
`~0` and `~1` escape `~` and `/`). No dotted-path guessing or JSON sniffing. Sequence
text is preserved, including whitespace, for downstream explicit normalization.
Default ids are `record-N`. Duplicate keys, non-string fields and malformed
records fail. Compression uses the same existing extension-based codecs.

GenBank/GenPept is a **best-effort ORIGIN subset**: LOCUS, wrapped DEFINITION,
ACCESSION/VERSION, numbered sequence lines and `//` termination. ORIGIN numbering
must be consecutive. Features, joins/complements, CONTIG-only records, declared
length validation and full GenBank semantics are not implemented. Extensions are
`.gb`, `.gbk`, `.genbank`, `.gp`; explicit format name is `genbank`.

The original `find_sequence_matches` still streams **whole records**. The new
`streaming_matcher.match_chunks` / `match_large_file` APIs keep only a pattern-span
lookbehind plus one chunk. The maximum pattern span defaults to 1,000,000;
file chunks are at most 1 Mi characters (default 65,536). Long unwrapped lines
are read in bounded pieces. FASTA/FASTQ headers and FASTQ `+` headers must fit a
chunk. This explicit mode supports FASTA/FASTQ/text, compression, both strands
and IUPAC, but no pattern packs or tabular inputs.

The new `sequence-match-events/1` output has `match` and per-record `complete`
events. Coordinates remain 1-based inclusive in the forward reference. Source
positions count supplied sequence chunks, excluding FASTA/FASTQ headers,
comments and quality; path input uses universal-newline text. The source hash
covers that supplied sequence-line text, not the whole file. Quality length and
ASCII validity must pass before a FASTQ completion event is emitted. Match
order is chunk completion then strand/start/gap, **not** the original full-record
strand order. Every exact hit is retained, including overlaps and boundary hits.
A limit/error terminates without a completion event for that record. The caller
must exhaust successfully to establish success for the entire file.

## S3 summaries (#29)

```sh
python -S sequence_tool.py summarize-hits hits.csv --gap-merge 10000 --top 3
python -S sequence_tool.py compare a.fa b.fa --alphabet dna --window 21 --step 5 --summary --identity-threshold 90 --min-coverage 50 --top 10
```

CSV clustering uses **distance between consecutive sorted start coordinates**,
with inclusive `<= gap_merge`; `gap_merge=0` disables merging even for duplicates.
Intervals cover the minimum start and maximum input end. The checked-in 5,966-row
example yields 1,246 regions at 10,000 and top counts 49/38/35. Its original file
is never written. The legacy coordinate default is `0-based-half-open`; callers
can select `1-based-inclusive`. Saved coordinates are preserved as-is; this does
not correct the legacy script's historical N-removal/offset behavior or certify
them as genomic coordinates.

Regions rank by count descending, then start/end/source row. Output includes
omitted region/item counts, top back-element frequencies with omitted counts,
and density per Mb. Supply `--sequence-length` for whole-sequence density;
otherwise the denominator is explicitly the observed hit bounding span (not a
whole-chromosome estimate). Empty input has null density. Source row references
are 1-based **data** row enclosing ranges; coordinate sorting may make membership
noncontiguous, so ranges are not exact lists. CSV SHA-256 identifies the source.
Clustering memory is proportional to input hits; only output size is top-bounded.

S2 summaries retain counts for `conserved`, `divergent`, `unknown`. Null identity
or coverage below the threshold is unknown. Consecutive same-class windows join
only when their covered columns touch/overlap. Opposite-class overlapping
windows may yield overlapping regions: these are window-evidence summaries,
not a disjoint biological partition or a homology inference. Source rows refer
to the exact contiguous S2 window range; original sequence hashes are retained.

## Ambiguity candidates and approximation (#30)

```sh
python -S sequence_tool.py candidates input.fa --mode pairwise --limit 1000 --seed 42
python -S sequence_tool.py candidates input.fa --mode sample --limit 100 --seed 42
python -S sequence_tool.py candidates input.fa --mode pairwise --targets targets.fa --rank 4 --seed 42
```

`preserve` keeps IUPAC text; `regex` emits a safe concrete-DNA pattern (no claim
that literal regex matching handles ambiguous target symbols); `exhaustive`
expands all combinations or errors at the cap; `sample` picks unique mixed-radix
indices reproducibly without allocating the full combination space. `pairwise`
is a separate fifth mode, not a seeded-sampling setting. Every ambiguous position
is a factor and each permitted concrete base is a level. All two-position level
pairs are covered; zero/one factors produce the sole sequence/all levels.
Pairwise supports at most 64 factors and fails if the candidate cap prevents
complete coverage. It does not guarantee three-way/full combinatorial coverage
or an optimal minimum candidate set. Seed and generated/omitted counts are saved.

Optional representative selection first computes exact literal candidate hits
against targets as integer bitsets. Greedy Boolean biclique factors yield a
bounded-rank approximation. One candidate per factor-membership signature is
retained. Output records rank/seed, omitted candidates, residual true matrix
bits and exact lost target coverage. This can lose pairwise coverage and genuine
candidate hits, explicitly flagged; it is a screening aid, never a substitute
for exact search. The implementation borrows binary low-rank ideas, not LittleBit
model quantization, scale correction, official kernels or model-quality claims.
FASTA/MEGA syntax is unaffected: generated concrete candidates can be wrapped
in SequenceRecords and passed through the same preprocess/export API.

## Analysis and preprocessing boundaries

Nucleotide GC/AT(U)/entropy denominators exclude ambiguous symbols. K-mers
(`1 <= k <= 6`) skip any window containing ambiguity. Optional codon-frame counts
nonoverlapping triplets with explicit frame 0/1/2 and leftover/excluded counts.
`translate_nucleotide` translates DNA/RNA to protein amino acids using standard genetic code
Table 1 across reading frames (0, 1, 2 or `all`), strands (`+`, `-` or `both`), and optional
`to_stop` truncations, resolving unambiguous IUPAC degenerate codons (such as ACN -> T).
The library `stop_symbol` option accepts a single character (default `*`) for rendering;
stop detection and counts always follow the genetic code, independent of that symbol.
`nucleotide_pairwise` reuses aligned comparison and counts transitions/transversions only for canonical substitutions.
RNA fold accepts canonical RNA, AU/GC and optional GU pairing, default minimum
loop 3, default length cap 300 (library ceiling 1,000). Ties prefer an unpaired
right endpoint, then the smallest partner. It maximizes noncrossing pairs;
it is not a free-energy model, pseudoknot predictor, or validated structure.
Small exhaustive tests compare scores against an independent recurrence.

Protein ionization is a separate output, leaving old simplified-charge outputs
unchanged. It uses Bjellqvist pKa constants, free termini, Henderson–Hasselbalch
fractions and 60 bisection steps over pH 0–14; no canonical root gives null pI.
Only canonical protein symbols are accepted; no ambiguous-residue charge guesses.
Constants and reference examples are cross-checked against
[Bio.SeqUtils.IsoelectricPoint](https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/IsoelectricPoint.py).
No BioPython runtime dependency or copied implementation is used.

Preprocessing normalizes ASCII symbols, optionally U/T, sanitizes unique ids,
filters lengths/ambiguity/id/motif, deduplicates id or sequence hash, and checks
retained aligned lengths when requested. Length filters count aligned columns;
ambiguity fractions exclude gap positions and include missing `?`. Each input
has a transform/filter receipt with hashes and original identity. FASTQ quality
is explicitly removed. CLI requires a **new** receipt path, opened exclusively.
MEGA requires equal nonempty lengths; it never aligns. FASTA↔TSV conversion
remains excluded by the owner's earlier decision.

## Acquisition, structure, and adapters

NCBI client operations have timeout, transfer/decompressed byte caps, <=3/s
pacing, transient-only bounded retry, optional cache, and injectable clock,
sleep and transport. Search pages retain `complete=false` and `next_start` at a
caller cap; duplicate pages/error payloads fail. PubMed ESearch is bounded to
10,000 ids; partition queries beyond that. XML supports article titles, abstracts
and cross-reference ids, not book records or internal entity declarations.
`Client.structure` accepts four-character legacy RCSB ids; extended PDB ids and
MMDB ASN.1 are not implied. Live access is optional, never a CI requirement.

For `fetch --operation records --database structure`, supply exactly one legacy
PDB id; `--format` is `cif` (default) or `pdb`. Output is a JSON
`structure-download/1` envelope containing the downloaded text and provenance,
not a raw coordinate file or geometry report. Sequence records retain default
`fasta` or explicit `genbank`. Other database/format combinations fail before
network I/O. `--format` is not accepted for PubMed, search or summary. Structure
search/summary still use NCBI's structure database and its numeric UIDs, not the
RCSB record-download endpoint or PDB ids.

Structure parsers support PDB fixed-column ATOM/HETATM with explicit elements,
and one mmCIF atom-site loop/data block with quoted/multiline tokens. Missing or
nonfinite coordinates, malformed rows and atom limits fail. No biological
assembly reconstruction or full dictionary validation. Geometry separates
models and incompatible alternate locations; spatial bins find bounded contacts,
residue/contact counts, chain interfaces, model centroids/bounds/extents, water/
metal/hetero contacts, CYS SG disulfide candidates (1.8–2.3 Å), charged side-chain
salt-bridge candidates (<=4 Å). A narrow H-bond candidate subset requires explicit
H attached to N/O (0.5–1.3 Å), oxygen acceptor, D–A <=3.5 Å, H–A <=2.5 Å, angle
>=120°. These are documented heuristic geometry thresholds, not proven bonds or
complete donor/acceptor chemistry. Missing atoms cannot be inferred from absent
records. Source hashes and reported atom references preserve evidence identity.

HTTP: `python -S sequence_service.py http` exposes `POST /tools/call` on loopback
port 8765, body `{"name":"match","arguments":{...}}`. It is a simple webhook/API,
not the MCP Streamable HTTP transport. Browser origins and transfer-encoding are
rejected. MCP: `python -S sequence_service.py stdio` supports protocol 2025-11-25,
initialize/initialized, ping, tools/list and tools/call, no unsolicited messages.
Both share strict schemas and computation/input/output caps, accept inline data,
and cannot read files or fetch URLs. Tool execution errors are structured;
notifications produce no replies. Docker is an optional example, not deployment.

`route_item` uses an explicit known alphabet hint when deterministic validation
suffices. Otherwise an optional classifier gets one atomic question. Important
items escalate, low confidence/unknown/malformed/failed decisions fall back. A
model cannot validate a tool invocation or execute a side effect. `JevClassifier`
implements the official Choice request/answer shape with bounded urllib I/O and
redirect refusal; a model name and API key must be supplied explicitly for live
calls. No key is read or used by tests. Choice is used for route selection;
Noul probabilities and Score rubrics are **not** converted into confidence.
Thresholds must be chosen from downstream evaluation; no default quality claim.
The example shows injected local/Pro handlers; actual llama.cpp/model deployment
and a domain benchmark remain operational/evaluation work under #2.

## References

- [NCBI GenBank sample record](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
- [NLM E-utilities parameter guide](https://www.nlm.nih.gov/dataguide/eutilities/utilities.html)
- [Nussinov RNA folding recurrence](https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Nussinov)
- [PDBx/mmCIF atom-site dictionary](https://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Categories/atom_site.html)
- [MEGA format manual](https://megasoftware.net/mega4/mega4.pdf)
- [MCP lifecycle](https://modelcontextprotocol.io/specification/2025-11-25/basic/lifecycle)
- [MCP tools](https://modelcontextprotocol.io/specification/2025-11-25/server/tools)
- [Official TypeSafe API](https://docs.typesafe.ai/api)
