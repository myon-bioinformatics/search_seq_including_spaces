"""Canonical amino-acid residue table (layer L0, stdlib only).

This module is the single source of per-residue properties. Other modules
must read values from here instead of defining their own tables.

Sources (see SOURCES for the machine-readable form):

* kd              Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132,
                  PMID 7108955. Defined for the 20 canonical residues only.
* mw_residue_avg  Average residue masses (Da) as listed by Expasy FindMod,
                  https://web.expasy.org/findmod/findmod_masses.html
                  (residue masses, i.e. free amino acid minus H2O).
* hbond_donor /   Side-chain hydrogen-bond donor/acceptor classes from
  hbond_acceptor  Pommie C et al. (2004) J Mol Recognit 17:17-32,
                  PMID 14872534 (IMGT physicochemical classes). Defined for
                  the 20 canonical residues only.
* charge_simplified
                  Tool convention for the "simplified" charge model:
                  D, E = -1; K, R = +1; every other canonical residue = 0
                  (His is treated as neutral). Not defined for U/O, whose
                  charge state is outside the simplified model.
* WATER_AVG       Average mass of H2O derived from the atomic weights
                  H = 1.00794 and O = 15.9994 (2 x 1.00794 + 15.9994).

Values that a cited source does not define are None, never a guess.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

TABLE_VERSION = "aa-props/1"

SOURCES: Mapping[str, str] = MappingProxyType({
    "kd": "Kyte & Doolittle 1982, J Mol Biol 157:105-132 (PMID 7108955)",
    "mw_residue_avg": "Expasy FindMod average residue masses "
                      "(https://web.expasy.org/findmod/findmod_masses.html)",
    "hbond": "Pommie et al. 2004, J Mol Recognit 17:17-32 (PMID 14872534), "
             "IMGT hydrogen donor/acceptor classes "
             "(https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html)",
    "charge_simplified": "tool convention: D,E=-1; K,R=+1; others 0; U/O undefined",
    "water_avg": "derived from atomic weights H = 1.00794, O = 15.9994 "
                 "(2 x 1.00794 + 15.9994 = 18.01528)",
})

# Average mass of water, added once per chain when computing a molecular
# weight from residue masses (used from PR B on). Derived value, see SOURCES.
WATER_AVG = 18.01528

STOP = "*"
GAP = "-"

SYMBOL_CLASSES = (
    "canonical",
    "rare-but-defined",
    "ambiguous",
    "stop",
    "gap",
    "invalid",
)


@dataclass(frozen=True)
class Residue:
    code: str
    name3: str
    name: str
    mw_residue_avg: float
    kd: float | None
    charge_simplified: int | None
    aromatic: bool
    sulfur: bool
    hbond_donor: bool | None
    hbond_acceptor: bool | None


def _r(code, name3, name, mw, kd, charge, aromatic, sulfur, donor, acceptor):
    return Residue(code, name3, name, mw, kd, charge, aromatic, sulfur, donor, acceptor)


# aromatic: F, W, Y (His is not counted as aromatic in this tool).
# sulfur: C, M. Selenocysteine (U) carries selenium, not sulfur.
_ROWS = (
    #  code  name3  name              mw        kd    chg   arom   S      donor  acceptor
    _r("A", "Ala", "Alanine",        71.0788,  1.8,  0,    False, False, False, False),
    _r("R", "Arg", "Arginine",      156.1875, -4.5,  1,    False, False, True,  False),
    _r("N", "Asn", "Asparagine",    114.1038, -3.5,  0,    False, False, True,  True),
    _r("D", "Asp", "Aspartic acid", 115.0886, -3.5, -1,    False, False, False, True),
    _r("C", "Cys", "Cysteine",      103.1388,  2.5,  0,    False, True,  False, False),
    _r("E", "Glu", "Glutamic acid", 129.1155, -3.5, -1,    False, False, False, True),
    _r("Q", "Gln", "Glutamine",     128.1307, -3.5,  0,    False, False, True,  True),
    _r("G", "Gly", "Glycine",        57.0519, -0.4,  0,    False, False, False, False),
    _r("H", "His", "Histidine",     137.1411, -3.2,  0,    False, False, True,  True),
    _r("I", "Ile", "Isoleucine",    113.1594,  4.5,  0,    False, False, False, False),
    _r("L", "Leu", "Leucine",       113.1594,  3.8,  0,    False, False, False, False),
    _r("K", "Lys", "Lysine",        128.1741, -3.9,  1,    False, False, True,  False),
    _r("M", "Met", "Methionine",    131.1926,  1.9,  0,    False, True,  False, False),
    _r("F", "Phe", "Phenylalanine", 147.1766,  2.8,  0,    True,  False, False, False),
    _r("P", "Pro", "Proline",        97.1167, -1.6,  0,    False, False, False, False),
    _r("S", "Ser", "Serine",         87.0782, -0.8,  0,    False, False, True,  True),
    _r("T", "Thr", "Threonine",     101.1051, -0.7,  0,    False, False, True,  True),
    _r("W", "Trp", "Tryptophan",    186.2132, -0.9,  0,    True,  False, True,  False),
    _r("Y", "Tyr", "Tyrosine",      163.1760, -1.3,  0,    True,  False, True,  True),
    _r("V", "Val", "Valine",         99.1326,  4.2,  0,    False, False, False, False),
    # Rare but defined residues: real amino acids, not ambiguity codes.
    _r("U", "Sec", "Selenocysteine", 150.0388, None, None, False, False, None,  None),
    _r("O", "Pyl", "Pyrrolysine",    237.3018, None, None, False, False, None,  None),
)

RESIDUES: Mapping[str, Residue] = MappingProxyType({row.code: row for row in _ROWS})

CANONICAL: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")
RARE_BUT_DEFINED: frozenset[str] = frozenset("UO")

# Ambiguity codes map to the set of canonical residues they may stand for.
# They are candidate sets, not values: how to turn them into numbers is a
# policy decided by the caller (default "exclude", from PR B on).
AMBIGUOUS: Mapping[str, frozenset[str]] = MappingProxyType({
    "B": frozenset("DN"),   # Asx
    "Z": frozenset("EQ"),   # Glx
    "J": frozenset("IL"),   # Xle
    "X": CANONICAL,         # unknown / unspecified
})


def symbol_class(ch: str) -> str:
    """Classify one sequence symbol (case-insensitive, ASCII only).

    Returns one of SYMBOL_CLASSES. Non-ASCII characters are always
    'invalid': Unicode case mapping would otherwise turn e.g. dotless i
    (U+0131) into 'I' or long s (U+017F) into 'S'. Raises ValueError unless
    ``ch`` is a single character.
    """
    if not isinstance(ch, str) or len(ch) != 1:
        raise ValueError(f"expected a single character, got {ch!r}")
    if not ch.isascii():
        return "invalid"
    up = ch.upper()
    if up in CANONICAL:
        return "canonical"
    if up in RARE_BUT_DEFINED:
        return "rare-but-defined"
    if up in AMBIGUOUS:
        return "ambiguous"
    if up == STOP:
        return "stop"
    if up == GAP:
        return "gap"
    return "invalid"
