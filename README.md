# search_seq_including_spaces
## Summary
![DNA_IMAGE](line_dna_short.webp)

- Narrowly, Find sequences including nucleotide(or amino acid) spaces and Output the results to every csv file.
- Broadly, you can find the speific sequences of front and back while ignoring the middle sequences)
- Language: Python (CI-tested on 3.10–3.14)
>__Note__ You can confirm for Bioinformatics Tips at Wiki.

![GitHub license](https://img.shields.io/github/license/myon-bioinformatics/search_seq_including_spaces)
![GitHub last commit](https://img.shields.io/github/last-commit/myon-bioinformatics/search_seq_including_spaces)
[![CodeQL](https://github.com/myon-bioinformatics/search_seq_including_spaces/actions/workflows/codeql.yml/badge.svg)](https://github.com/myon-bioinformatics/search_seq_including_spaces/actions/workflows/codeql.yml)

[![GitHub followers](https://img.shields.io/github/followers/myon-bioinformatics?style=social)](https://github.com/myon-bioinformatics)
[![Reddit User Karma](https://img.shields.io/reddit/user-karma/combined/myon_reddit?style=social)](https://www.reddit.com/user/myon_reddit/)
[![Twitter Follow](https://img.shields.io/twitter/follow/myonitbusiness?style=social)](https://twitter.com/myonitbusiness)



## Protein sequence tools (in development)
Standard library only. CI-tested on Python 3.10–3.14. Values are sequence-derived hints, not structure.

```sh
python sequence_tool.py fasta validate in.fasta
python sequence_tool.py protein features in.fasta            # JSON (protein-hints/1)
python sequence_tool.py protein profile in.fasta --window 9  # per-residue TSV (or --json / --jsonl)
python sequence_tool.py protein roi in.fasta --ascii         # regions of interest (heuristic)
```

Ambiguity codes B/Z/J/X are excluded by default (`--ambiguity mean-of-candidates` to opt in). Exit code 2 means an input error. Details, thresholds and limitations: [docs/PROTEIN.md](docs/PROTEIN.md).

## Aligned sequence comparison

The stdlib-only comparison tools accept already-aligned DNA, RNA, or protein
pairs. They do not perform alignment or infer homology.

```sh
python -S sequence_tool.py compare a.fasta b.fasta --alphabet dna
python -S sequence_tool.py compare a.fasta b.fasta --alphabet dna --window 21 --step 5
python -S sequence_tool.py compare a.fasta b.fasta --alphabet protein --ascii --width 60
```

S1 reports exact identity and comparison coverage. S2 adds sliding-window
identity/local divergence and an ASCII aligned diff while preserving the same
gap and ambiguity policies. Details and scientific boundaries:
[docs/COMPARISON.md](docs/COMPARISON.md).

## SEARCH_SEQUENCE_IMAGE
![SEARCH_IMAGE](SEARCH_SEQUENCE_INCLUDING_SPACES_IMAGE.webp)

## References
- About Nucleotides:https://www.ncbi.nlm.nih.gov/nuccore/
- About amino-acid: https://www.ncbi.nlm.nih.gov/protein/


## Test evidence

Install test-only tools with `python -m pip install -r tests/requirements.txt`.
CI preserves native pytest JSONL and JUnit while retaining the `python -S`
runtime isolation lane. See [recording and exploration](docs/pytest-observations.md).
