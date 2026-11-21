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

## Reference alignment and structural records

- `tRNA5.sto` is the unmodified two-block, five-row Stockholm file from
  [Infernal](https://github.com/EddyRivasLab/infernal/blob/f80cd7618636e44e3963f65625cea78c38a8d823/tutorial/tRNA5.sto).
  The reference test verifies the file SHA-256. Its column annotation is a
  consensus across the supplied alignment. The upstream notice accompanies
  the fixture in [INFERNAL-LICENSE](INFERNAL-LICENSE).
- `1crn-sulfur.pdb` retains the three SSBOND rows and six cysteine SG coordinates
  from [PDB entry 1CRN](https://www.rcsb.org/structure/1CRN), using this
  [PDB file](https://github.com/rdkit/rdkit-tutorials/blob/master/data/1CRN.pdb)
  (Git blob `e05c7c363fe90c17c19eb654c77d15defcdff478`). Trailing padding
  was removed. SSBOND records are labels for evaluating geometry.

## RNA neighbourhood coordinates

[1EHZ](https://www.rcsb.org/structure/1EHZ) is the yeast phenylalanine tRNA
crystal structure. `1ehz-g19-c56.pdb` contains the named base nitrogen and oxygen
ATOM records for chain A residues G19 and C56 from the
[deposited PDB](https://files.rcsb.org/download/1EHZ.pdb), with trailing padding
removed. After replacing modified residues in the SEQRES sequence with their
MODRES parents, its first 73 residues match the ungapped Infernal tRNA1 row;
the deposited chain also has a terminal CCA. The companion JSON explicitly
maps just the two fixture residues. Distances are coordinate observations.
