# Test fixtures

Public protein sequences from UniProtKB, retrieved 2026-09-26. The sha256
values are of the residue string only (upper case, no header, no line breaks).

| File | Source | Residues | Length | sha256 |
| --- | --- | --- | --- | --- |
| `P00698_mature.fasta` | https://rest.uniprot.org/uniprotkb/P00698.fasta (hen egg-white lysozyme C) | 19-147 of the 147-aa precursor (mature chain; signal peptide 1-18 removed) | 129 | `eaf1f8d7f5bf78bbefe1100d912a09893361c57dda057f634adb5c5af6dfe79c` |
| `P37840.fasta` | https://rest.uniprot.org/uniprotkb/P37840.fasta (human alpha-synuclein) | 1-140 (full entry) | 140 | `9a06387610edd099466d614941c42a433f444feb6a8f25375bb0b6ad8ec9c6f4` |
| `P37840_bom.fasta` | byte-for-byte copy of `P37840.fasta` with a UTF-8 byte order mark (EF BB BF) prepended | same as `P37840.fasta` | 140 | same as `P37840.fasta` |

Molecular-weight cross-check against UniProt sequence masses (sum of the
average residue masses in `protein/residues.py` plus one water, rounded
to 1 Da):

- Automated: P37840 gives 14,460 Da, equal to UniProt
  (`tests/test_residues.py`).
- Checked once by hand at retrieval, not automated: the P00698 precursor
  (147 aa, not a fixture) gives 16,239 Da, equal to UniProt.

The first two sequences are also the fixtures of the ROI threshold-sensitivity
measurement in the design document, so keep them unchanged; later PRs use
them as regression goldens.
