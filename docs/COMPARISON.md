# Sequence comparison S1

S1 provides deterministic, stdlib-only comparison of **already aligned**
DNA, RNA, or protein sequences.

It does **not** align sequences and does not compute evolutionary similarity.
Needleman-Wunsch / Smith-Waterman, BLOSUM/PAM, MSA, sliding-window identity,
ASCII diff, and conservation/divergence ROI are deferred.

## CLI

```bash
python -S sequence_tool.py compare a.fasta b.fasta --alphabet dna
```

Optional controls:

```text
--gap-policy exclude|mismatch
--ambiguity-policy exclude|strict|compatible
--start N
--end N
```

Coordinates are **1-based inclusive aligned columns**.

Both files must contain exactly one FASTA record (headerless plain sequence
text is also accepted by the shared FASTA splitter), and the normalized
aligned lengths must be equal.

## Metrics

`percent_identity` is literal symbol equality among compared columns:

```text
matches / compared_positions * 100
```

`coverage_percent` is:

```text
compared_positions / selected_aligned_columns * 100
```

If there are zero comparable positions, `percent_identity` is JSON
`null`, not 0.

### Gap policy

- `exclude` (default): any column containing `-` is excluded.
- `mismatch`: a one-sided gap is a mismatch; a gap-gap column is still
  excluded because it carries no residue/base comparison.

### Ambiguity policy

- `exclude` (default): columns containing an ambiguity symbol are excluded.
- `strict`: ambiguity symbols are compared literally, so `N` vs `N`
  is a literal match and `N` vs `A` is a mismatch.
- `compatible`: `percent_identity` remains literal identity, while the
  separate `percent_compatible` reports candidate-set overlap.

For example, under DNA `compatible`:

```text
R vs A  compatible
R vs G  compatible
R vs C  not compatible
N vs A  compatible
```

This compatible fraction is **not** exact identity.

## Alphabets

DNA supports canonical `A C G T` and IUPAC ambiguity
`R Y S W K M B D H V N`.

RNA supports canonical `A C G U` with the corresponding U-based IUPAC
candidate sets.

Protein supports the 20 canonical residues plus defined `U/O`, with
`B/Z/J/X` as ambiguous symbols. No substitution matrix is used in S1.

## Output contract

The JSON document includes:

- schema: `sequence-comparison/v1`
- alphabet
- aligned selected region
- match / mismatch / gap / ambiguity counts
- percent identity
- compatibility percentage when requested
- comparison coverage
- gap and ambiguity policies
- SHA-256 of both normalized aligned sequences
- FASTA record IDs when invoked through the CLI

The coordinate convention and the fact that the input must already be
aligned are recorded explicitly.

## Scientific boundary

S1 answers: “For these aligned columns, how often are the symbols identical,
and how much of the requested region was actually compared?”

It does not infer homology, evolutionary distance, functional equivalence,
or structural similarity. A high identity over low coverage must not be
presented as high whole-sequence similarity.

## Next

S2 may add sliding-window identity, local divergence, and ASCII diff.

S3 may add conservation/divergence regions and bounded multi-record summaries.
