# Anchor/gap matcher — Phase 1

`anchor_matcher.py` is one stdlib-only importable module. It does not read files,
configuration or the network at import time. The legacy `definition_main.py`
remains intact. CLI/file readers and pseudogene/TSD presets are later phases.

```python
import json
from anchor_matcher import find_matches

report = find_matches("a c\nN t g", "AC", "TG", min_gap=1)
assert report["matches"][0]["middle"] == "N"
print(json.dumps(report))
```

Use `min_gap=2` for exactly two residues, or `min_gap=2, max_gap=5` for the
inclusive range 2–5. Anchors are nonempty literal letter strings, not regexes.
Whitespace is removed and ASCII lowercase is uppercased in both input and
anchors. N is retained as a residue. Headers, digits, alignment gap markers,
non-ASCII letters and controls other than whitespace are rejected. Supply raw
sequence text, not an entire FASTA/FASTQ file.

## Coordinates and multiplicity

All spans are **1-based inclusive**, matching the existing comparison outputs:

- `span`, `left_span`, `middle_span`, `right_span` refer to the normalized
  forward input sequence, even for a reverse-strand match.
- `source_span` and `middle_source_span` refer to original Python string
  character positions, including removed whitespace in the enclosing interval.
  They are not byte offsets or base positions in a file.
- Zero-length middle strings have null middle spans, not an invalid end/start
  interval. `normalize_sequence()` separately exposes an exact residue-to-source
  mapping through its `source_positions` tuple.

Scanning starts from the left anchor only. Every valid (start, gap, strand)
tuple occurs once; overlapping hits and different valid gap lengths are kept.
`strand="both"` emits forward hits before reverse hits; palindromic hits on
opposite strands remain separate observations. Within each strand, oriented
start then gap length determines order. Reverse-strand left/right anchors are
defined in the oriented reverse-complement string, so their forward-reference
spans appear in reverse order. `middle` always follows the searched orientation.

## DNA ambiguity and reverse complement

Literal forward mode supports generic ASCII letter sequences, including
proteins. `iupac=True` enables DNA-code compatibility: possible base sets must
intersect in both input and anchor. An ambiguous compatible match does not
establish which base is present. Reverse search (`strand="-"` or `"both"`) and
`reverse_complement()` accept IUPAC DNA symbols only; RNA U is rejected.

Base sets follow the [NCBI nucleotide-code table](https://www.ncbi.nlm.nih.gov/books/NBK44863/table/sequencesquickstart.Td/).
Complement pairs preserve ambiguity sets (R/Y, K/M, B/V, D/H; S/W/N are
self-complementary). Tests check the complement operation against the base
sets and its double-application identity.

## Output and limits

The `anchor-gap/1` dict is JSON serializable. Provenance contains a hash and
length of the original text plus normalized anchors/gap/strand configuration;
it does not infer a repository, organism or scientific interpretation. A hash
is an identity aid, not encryption.

`max_matches` defaults to 10000. Exceeding it raises ValueError instead of
returning a prefix marked complete. `complete` indicates completion of this
requested search, not biological validity. The feasible gap range is clamped
to the available sequence. This direct baseline scales with sequence length,
gap range and anchor lengths; it is not an indexed genome-scale search engine.

No CLI integration, alignment, edit-distance matching, TSD/pseudogene verdict
or network acquisition is added in this phase.
