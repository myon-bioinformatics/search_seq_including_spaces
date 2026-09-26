# Golden outputs

Expected CLI output for the two fixtures in `tests/fixtures/`, compared byte
for byte by `tests/test_cli.py`. The Python version in the provenance block
is stored as `<python>` and substituted before the comparison.

| File | Command |
| --- | --- |
| `<name>.features.json` | `python -S sequence_tool.py protein features tests/fixtures/<name>.fasta` |
| `<name>.profile.tsv` | `python -S sequence_tool.py protein profile tests/fixtures/<name>.fasta` |

Checked against independent values when the files were created (2026-09-27):

- Simplified net charge: lysozyme (mature) +8, alpha-synuclein -9, as
  recorded in the design document's measurements.
- Molecular weight: alpha-synuclein 14,460 Da and the lysozyme precursor
  16,239 Da equal the UniProt sequence masses (automated in
  `tests/test_features.py`).
- Cys, Pro, Gly and aromatic positions agree with the UniProt sequences
  (lysozyme Cys 6, 30, 64, 76, 80, 94, 115, 127 in mature numbering).
- `hydropathy_window`, `charge_window` and `entropy_window` agree at every
  position with the separate prototype used for the ROI threshold
  measurements (window 9, entropy window 12).

The outputs are identical on Python 3.11 to 3.14. Regenerate a file only
when a change is meant to alter the output, and say why in the pull request.
