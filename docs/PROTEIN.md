# Protein sequence tools

Standard-library-only tools (Python 3.11+) that compute descriptive values
from a protein sequence and point to regions worth a closer look. Every
value is labelled with where it comes from. Nothing here predicts
structure, binding or function.

```sh
python sequence_tool.py fasta validate in.fasta
python sequence_tool.py protein features in.fasta              # JSON
python sequence_tool.py protein profile in.fasta --window 9    # TSV (--json, --jsonl)
python sequence_tool.py protein roi in.fasta --ascii           # JSON by default
```

`FILE` may be `-` for standard input. Exit code 0 means success and 2 means
an input error (unreadable file, invalid sequence, bad option). The same
input gives the same output bytes; JSON keys are sorted.

## Input

| Class | Symbols | Default handling |
| --- | --- | --- |
| canonical | `ACDEFGHIKLMNPQRSTVWY` | table values |
| rare-but-defined | `U` (Sec), `O` (Pyl) | counted and weighed; no Kyte-Doolittle or simplified charge value, so left out of those windows |
| ambiguous | `B` (D/N), `Z` (E/Q), `J` (I/L), `X` (any) | `--ambiguity` policy, below |
| stop | `*` | a trailing one is removed; internal ones are errors unless `--allow-internal-stop` |
| gap | `-` | error unless `--strip-gaps` |
| invalid | anything else, including non-ASCII letters | error naming the position and symbol |

Positions are 1-based. `input_position` is the position in the input as
written (whitespace not counted), so it stays correct after gaps are
stripped.

`--ambiguity`:

- `exclude` (default): ambiguous residues get no value, are left out of
  windows, and molecular weight is null. Nothing is filled with 0 or a mean.
- `mean-of-candidates`: hydropathy, molecular weight and composition
  fractions use the mean over the candidates (`X`: the plain mean of the 20
  canonical residues).
- `error`: stop at the first ambiguous residue.

Charge never uses a candidate mean. The net charge gives `range: [min, max]`
over the candidates and a `value` only when the range is one number. Every
document records `ambiguity_policy`, `ambiguous_count` and
`excluded_positions`.

## Output documents (`protein-hints/1`)

Each document has `record` (id, length, sha256 of the residues, symbol
classes, ambiguity), `provenance` (tool and Python version, table version,
parameters), `field_classes`, the data and `limitations`.

Every data field has one of five evidence classes, set only in
`protein/provenance.py` (`FIELD_CLASSES`):

| Class | Meaning | Used for |
| --- | --- | --- |
| sequence-derived | computed directly from the sequence | features, profiles |
| heuristic | rule-based pointers, no physical claim | ROI regions |
| predicted | statistical or propensity model | (later: secondary structure) |
| structure-derived-experimental | measured on an experimental structure | (later) |
| structure-derived-predicted | measured on a predicted structure | (later) |

## Features and profiles

| Field | Definition | Source |
| --- | --- | --- |
| `length`, `composition` | residue counts (stop markers excluded) | - |
| `fraction_aromatic`, `fraction_sulfur`, `fraction_pro_gly` | (F+W+Y), (C+M), (P+G) over residues | - |
| `molecular_weight_avg` | sum of average residue masses + one water; unmodified, no disulfides | Expasy average residue masses |
| `net_charge_simplified` | D,E = -1; K,R = +1; H and termini = 0 | tool convention (about pH 7) |
| `cys_positions` etc. | input positions of C, P, G, F/W/Y | - |
| `hydropathy`, `hydropathy_window` | Kyte-Doolittle value, mean over a centred odd window (default 9) | Kyte & Doolittle 1982 |
| `charge`, `charge_window` | simplified charge; window net charge = window x mean over usable residues (the plain sum when all are usable) | - |
| `entropy_window` | Shannon entropy (bits) of residues of known identity over 12 residues starting 6 before the position | - |

Windows are null within half a window of either end and when fewer than
half of their residues are usable. Window sums use exact integer
arithmetic, so results are the same on every Python version.

## Regions of interest (ROI)

A region of interest is a run of residues where at least one value crosses
an absolute threshold. It says where to look first and always lists why.
It is a heuristic, not a functional, binding, interface or drug site, and no
score is given.

| Trigger | Rule (default) | Threshold source |
| --- | --- | --- |
| `charge_cluster` | \|window net charge, w=9\| >= 3 | tool default |
| `charge_transition` | net charge of 5 residues left and 5 right: opposite signs, difference >= 4 | tool default |
| `hydropathy_peak` | mean KD, w=9 >= 1.5 | tool default (Kyte-Doolittle's 1.6 is for window 19, not reused) |
| `hydropathy_gradient` | \|mean KD 5 right - 5 left\| >= 3.0 | tool default |
| `cys_cluster` | literal known-residue Cys count in w=9 >= 2 | tool default |
| `aromatic_cluster` | literal known-residue F/W/Y count in w=9 >= 3 | tool default |
| `pro_gly_cluster` | literal known-residue P+G count in w=9 >= 4 | tool default |
| `low_complexity` | entropy, w=12 <= 2.2 bits | SEG-like defaults (window 12, 2.2); not SEG itself |

For the three count triggers, ambiguous residues are excluded and the known
matching residues are counted literally; missing/ambiguous positions are not
extrapolated to a full window.

Triggered positions are joined across gaps of at most `--merge-gap` (2)
residues, runs shorter than `--min-len` (5) are dropped, and a stop marker
ends a run. Each region reports its triggers with residue counts, the rule
and threshold source, the sequence context (6 residues each side), the
simplified net charge inside it, `boundary_uncertainty` (half a window) and
`edge_truncated` when a trigger that actually fired in that region touches
that trigger's own evaluable sequence boundary. Coverage above 50% adds a warning. Override thresholds with
`--threshold NAME=VALUE`; the provenance then marks them as `user`.

`--ascii` prints tracks per block (`--width`, default 60, minimum 20; every
line fits in width + 8 characters, no colour):

```text
        61        71        81        91        101       111
        |         |         |         |         |         |
SEQ     EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP
CHG9    ..............................................--..----------
HYD9    .....+#++#+###++++.++..+++++++#++....................+++....
ROI      [===========]                                [=============
         ROI-03                                       ROI-04 ->
```

`CHG9`: `+` net >= +3, `-` <= -3, `.` otherwise. `HYD9`: `#` >= 1.5, `+` 0
to 1.5, `.` < 0. `?` marks positions without a value.

### Threshold sensitivity [measured]

Measured on 2026-09-26 on the two test fixtures (hen lysozyme mature chain,
P00698 residues 19-147; human alpha-synuclein, P37840), and reproduced by
`tests/test_roi.py`:

| Setting | Lysozyme: ROIs / coverage / longest | Alpha-synuclein |
| --- | --- | --- |
| defaults | 2 / 10% / 7 (75-80, 111-117) | 4 / 40% / 30 (14-20, 47-52, 62-74, 107-136) |
| charge >= 2 | 2 / 14% / 10 | 5 / 49% / 31 |
| charge >= 4 | 1 / 5% / 6 | 4 / 24% / 13 |
| h9 >= 1.2 | 2 / 12% / 8 | 5 / 48% / 30 |
| window 7 (thresholds not rescaled) | 2 / 11% / 7 | 5 / 37% / 13 |
| merge gap 0 | 1 / 5% / 6 | 2 / 23% / 26 |
| merge gap 4 | 2 / 10% / 7 | 3 / 73% / 66 |
| min len 3 | 3 / 12% / 7 | 5 / 42% / 30 |

Every setting changes the result, the merge gap most of all. These are
regression values for two sequences, not tuned or validated thresholds; no
comparison with known sites has been made.

## Scientific limitations

- The simplified charge ignores His, the termini and pH. It is not a
  pKa-based charge.
- Hydropathy is a residue-scale average. It does not show burial or surface
  exposure, and folding is driven mainly by the hydrophobic effect (Dill
  1990), not by charge alone.
- Window smoothing blurs region boundaries by about half a window.
- "Hotspot" means residues whose alanine mutation strongly weakens binding
  at an interface (Bogan & Thorn 1998). This tool never uses the word for
  its regions.
- Conservation from a multiple alignment predicts functional residues
  better than physicochemical values (Capra & Singh 2007). An ROI is not a
  substitute.
- Aggregation-prone regions need dedicated, validated methods (TANGO 2004,
  AGGRESCAN 2007). A hydropathy peak is not an aggregation region.
- Pockets need a structure (fpocket 2009, P2Rank 2018).

## References

- Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132. [PMID 7108955](https://pubmed.ncbi.nlm.nih.gov/7108955/)
- Pommie C et al. (2004) J Mol Recognit 17:17-32 (IMGT hydrogen-bond classes). [PMID 14872534](https://pubmed.ncbi.nlm.nih.gov/14872534/)
- Expasy average residue masses: https://web.expasy.org/findmod/findmod_masses.html
- Dill KA (1990) Biochemistry 29:7133-7155. [PMID 2207096](https://pubmed.ncbi.nlm.nih.gov/2207096/)
- Bogan AA, Thorn KS (1998) J Mol Biol 280:1-9. [PMID 9653027](https://pubmed.ncbi.nlm.nih.gov/9653027/)
- Capra JA, Singh M (2007) Bioinformatics 23:1875-1882. [PMID 17519246](https://pubmed.ncbi.nlm.nih.gov/17519246/)
- Fernandez-Escamilla AM et al. (2004) Nat Biotechnol 22:1302-1306 (TANGO). [PMID 15361882](https://pubmed.ncbi.nlm.nih.gov/15361882/)
- Conchillo-Sole O et al. (2007) BMC Bioinformatics 8:65 (AGGRESCAN). [PMID 17324296](https://pubmed.ncbi.nlm.nih.gov/17324296/)
- Le Guilloux V et al. (2009) BMC Bioinformatics 10:168 (fpocket). [PMID 19486540](https://pubmed.ncbi.nlm.nih.gov/19486540/)
- Krivak R, Hoksza D (2018) J Cheminform 10:39 (P2Rank). [PMID 30109435](https://pubmed.ncbi.nlm.nih.gov/30109435/)
- Wootton JC, Federhen S (1993) Comput Chem 17:149-163 (SEG). [doi:10.1016/0097-8485(93)85006-X](https://doi.org/10.1016/0097-8485(93)85006-X)
