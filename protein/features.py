"""Sequence-derived features and window profiles (layer L2, stdlib only).

protein_features() and protein_profiles() each return one protein-hints/1
document (a dict ready for json.dumps) for one sequence.

Ambiguity codes B, Z, J, X (see residues.AMBIGUOUS):

- ambiguity="exclude" (default): ambiguous residues get no hydropathy value,
  are left out of every window, and molecular weight is null when any are
  present. Values are never filled with 0 or with a mean.
- ambiguity="mean-of-candidates" (opt-in): hydropathy, molecular weight and
  the composition fractions use the mean over the candidate residues (for X,
  the plain mean of the 20 canonical residues, independent of composition).
- ambiguity="error": the first ambiguous residue raises SequenceError.

Two rules do not depend on the policy. Charge never uses a candidate mean:
ambiguous residues get no per-residue charge, and the whole-sequence net
charge reports the [min, max] over their candidates. Window entropy counts
only residues of known identity (canonical, U, O).

U and O are real residues: they count in composition and molecular weight,
but Kyte-Doolittle and the simplified charge do not define them, so those
values are null and they are left out of those windows.

Hydropathy and charge windows are centred and odd; within (window - 1) / 2
residues of either end the value is null (edge="none"). Any window with fewer
than half of its residues usable is null. Window sums are exact (integer arithmetic on the
table values), so results do not depend on summation order or on the
Python version.

Floats are rounded to 4 decimal places.
"""

import math
from collections import Counter
from collections.abc import Mapping, Sequence
from fractions import Fraction
from functools import lru_cache
from itertools import accumulate

from .parse import Prepared, SequenceError, prepare
from .provenance import FIELD_CLASSES, SCHEMA, field_classes, provenance, sequence_sha256
from .residues import AMBIGUOUS, CANONICAL, RESIDUES, TABLE_VERSION, WATER_AVG

AMBIGUITY_POLICIES = ("exclude", "mean-of-candidates", "error")
SCALES = ("kyte-doolittle",)
CHARGE_MODELS = ("simplified",)
EDGES = ("none",)
DECIMALS = 4

_DEFINED = ("canonical", "rare-but-defined")
_COMPOSITION_KEYS = tuple(sorted(CANONICAL)) + ("O", "U") + tuple(sorted(AMBIGUOUS))
_CHARGE_ASSUMPTION = "D,E = -1; K,R = +1; H and termini = 0 (about pH 7)"
_LIMITATIONS = (
    "simplified charge ignores His, the termini and pH; it is not a pKa-based charge",
    "hydropathy is a residue-scale average; it does not show burial or surface exposure",
    "values describe the sequence only; they are not structure, interaction or function",
)


def protein_features(seq: str | Prepared, *, ambiguity: str = "exclude") -> dict:
    """Whole-sequence features: composition, fractions, MW, net charge, positions."""
    p = _as_prepared(seq)
    _check_ambiguity(p, ambiguity)
    features = [
        _feature("length", _length(p), "count of residues (stop markers excluded)"),
        _feature("composition", _composition(p), "count of each symbol",
                 "every canonical, rare-but-defined and ambiguous symbol is listed"),
        _feature("fraction_aromatic", _fraction(p, ambiguity, lambda c: RESIDUES[c].aromatic),
                 "(F + W + Y) / residues", _fraction_assumption(ambiguity)),
        _feature("fraction_sulfur", _fraction(p, ambiguity, lambda c: RESIDUES[c].sulfur),
                 "(C + M) / residues", _fraction_assumption(ambiguity)),
        _feature("fraction_pro_gly", _fraction(p, ambiguity, lambda c: c in "PG"),
                 "(P + G) / residues", _fraction_assumption(ambiguity)),
        _molecular_weight(p, ambiguity),
        _net_charge(p),
        _feature("cys_positions", _positions(p, "C"), "input positions of C",
                 "ambiguous residues are not included"),
        _feature("pro_positions", _positions(p, "P"), "input positions of P",
                 "ambiguous residues are not included"),
        _feature("gly_positions", _positions(p, "G"), "input positions of G",
                 "ambiguous residues are not included"),
        _feature("aromatic_positions", _positions(p, "FWY"), "input positions of F, W, Y",
                 "ambiguous residues are not included"),
    ]
    params = {"ambiguity": ambiguity, "charge_model": "simplified", "pH": None}
    return _document(p, ambiguity, params, [f["name"] for f in features],
                     features=features)


def protein_profiles(seq: str | Prepared, *, window: int = 9,
                     scale: str = "kyte-doolittle", charge_model: str = "simplified",
                     edge: str = "none", ambiguity: str = "exclude",
                     entropy_window: int = 12) -> dict:
    """Per-residue values and window profiles.

    hydropathy_window  mean Kyte-Doolittle value over the centred window
    charge_window      local net charge: window x mean simplified charge
                       (equals the plain sum when every residue is usable)
    entropy_window     Shannon entropy in bits over `entropy_window` residues
                       starting entropy_window // 2 before the residue
    """
    _check_choice("scale", scale, SCALES)
    _check_choice("charge_model", charge_model, CHARGE_MODELS)
    _check_choice("edge", edge, EDGES)
    if not _is_int(window) or window < 1 or window % 2 == 0:
        raise ValueError(f"window must be a positive odd integer, got {window!r}")
    if not _is_int(entropy_window) or entropy_window < 1:
        raise ValueError(f"entropy_window must be a positive integer, got {entropy_window!r}")
    p = _as_prepared(seq)
    _check_ambiguity(p, ambiguity)

    # The class of a symbol is fixed, so per-residue values are looked up per symbol.
    symbol_classes = dict(zip(p.residues, p.classes))
    hydropathy = {c: _hydropathy(c, cls, ambiguity) for c, cls in symbol_classes.items()}
    charge = {c: _charge(c, cls) for c, cls in symbol_classes.items()}
    hydropathy_out = {c: _round(v) for c, v in hydropathy.items()}
    hydropathy_window = window_values(p.residues, hydropathy, window)
    charge_window = window_values(p.residues, charge, window, scale_to_window=True)
    entropy = entropy_values(p.residues, p.classes, entropy_window)

    rows = [
        {
            "position": i + 1,
            "input_position": p.input_positions[i],
            "residue": code,
            "class": p.classes[i],
            "hydropathy": hydropathy_out[code],
            "charge": charge[code],
            "hydropathy_window": hydropathy_window[i],
            "charge_window": charge_window[i],
            "entropy_window": entropy[i],
        }
        for i, code in enumerate(p.residues)
    ]
    params = {"window": window, "entropy_window": entropy_window, "scale": scale,
              "edge": edge, "charge_model": charge_model, "pH": None,
              "ambiguity": ambiguity}
    names = ("hydropathy", "charge", "hydropathy_window", "charge_window", "entropy_window")
    return _document(p, ambiguity, params, names, residues=rows)


# ---------------------------------------------------------------- windows

def window_values(symbols: Sequence[str], value_of: Mapping[str, Fraction | int | None],
                  window: int, *, scale_to_window: bool = False) -> list[float | None]:
    """Centred window mean of value_of[symbol] over the symbols with a value.

    O(n) with prefix sums. The values are exact rationals (table values are
    exact decimals), scaled to integers by their common denominator, so the
    sums are exact. None where the window runs past either end or fewer
    than half of its symbols have a value. With scale_to_window the mean is
    multiplied by the window length.
    """
    need = (window + 1) // 2
    exact = {s: None if v is None else Fraction(v) for s, v in value_of.items()}
    denominator = math.lcm(*(v.denominator for v in exact.values() if v is not None))
    scaled = {s: None if v is None else int(v * denominator) for s, v in exact.items()}
    per_symbol = [scaled[symbol] for symbol in symbols]
    sums = list(accumulate((0 if v is None else v for v in per_symbol), initial=0))
    counts = list(accumulate((v is not None for v in per_symbol), initial=0))
    half = window // 2
    factor = window if scale_to_window else 1
    out: list[float | None] = [None] * len(symbols)
    for i in range(half, len(symbols) - half):
        used = counts[i + half + 1] - counts[i - half]
        if used >= need:
            total = sums[i + half + 1] - sums[i - half]
            out[i] = _round(total * factor / (used * denominator))
    return out


def entropy_values(residues: str, classes: tuple[str, ...],
                   window: int) -> list[float | None]:
    """Shannon entropy (bits) of the residues of known identity in each window.

    The window for position i covers i - window // 2 .. i - window // 2 +
    window - 1 (0-based). None past either end or when fewer than half of
    the window's residues are of known identity.
    """
    n = len(residues)
    need = (window + 1) // 2
    out: list[float | None] = [None] * n
    counts: Counter[str] = Counter()
    known = [cls in _DEFINED for cls in classes]
    # counts holds the known residues of residues[lo:hi]; total is their number.
    lo = hi = total = 0
    for i in range(n):
        start = i - window // 2
        if start < 0 or start + window > n:
            continue
        while hi < start + window:
            if known[hi]:
                counts[residues[hi]] += 1
                total += 1
            hi += 1
        while lo < start:
            if known[lo]:
                counts[residues[lo]] -= 1
                total -= 1
                if not counts[residues[lo]]:
                    del counts[residues[lo]]
            lo += 1
        if total >= need:
            out[i] = _entropy(tuple(sorted(counts.values())))
    return out


@lru_cache(maxsize=4096)
def _entropy(counts: tuple[int, ...]) -> float:
    total = sum(counts)
    return _round(math.fsum(-(c / total) * math.log2(c / total) for c in counts))


# ---------------------------------------------------------------- per residue

def _hydropathy(code: str, cls: str, ambiguity: str) -> Fraction | None:
    if cls in _DEFINED:
        return _kd(code)
    if cls == "ambiguous" and ambiguity == "mean-of-candidates":
        return _candidate_mean_kd(code)
    return None


def _charge(code: str, cls: str) -> int | None:
    return RESIDUES[code].charge_simplified if cls in _DEFINED else None


@lru_cache(maxsize=None)
def _kd(code: str) -> Fraction | None:
    kd = RESIDUES[code].kd
    return None if kd is None else Fraction(str(kd))  # exact decimal table value


@lru_cache(maxsize=None)
def _candidate_mean_kd(code: str) -> Fraction:
    candidates = sorted(AMBIGUOUS[code])
    return sum((_kd(c) for c in candidates), Fraction(0)) / len(candidates)


# ---------------------------------------------------------------- features

def _length(p: Prepared) -> int:
    return sum(cls != "stop" for cls in p.classes)


def _composition(p: Prepared) -> dict[str, int]:
    counts = Counter(c for c, cls in zip(p.residues, p.classes) if cls != "stop")
    return {code: counts.get(code, 0) for code in _COMPOSITION_KEYS}


def _fraction(p: Prepared, ambiguity: str, has_property) -> float | None:
    numerator, denominator = Fraction(0), 0
    for code, count in Counter(p.residues).items():
        if code in RESIDUES:
            numerator += count * bool(has_property(code))
            denominator += count
        elif code in AMBIGUOUS and ambiguity == "mean-of-candidates":
            candidates = AMBIGUOUS[code]
            numerator += count * Fraction(sum(bool(has_property(c)) for c in candidates),
                                          len(candidates))
            denominator += count
    return None if denominator == 0 else _round(numerator / denominator)


def _fraction_assumption(ambiguity: str) -> str:
    if ambiguity == "mean-of-candidates":
        return "ambiguous residues count as the candidate mean (X: mean of 20 residues)"
    return "ambiguous residues are left out of numerator and denominator"


def _molecular_weight(p: Prepared, ambiguity: str) -> dict:
    method = "sum of average residue masses + one water (Expasy average masses)"
    if "stop" in p.classes:
        return _feature("molecular_weight_avg", None, method,
                        "not computed: internal stop markers split the sequence",
                        unit="Da")
    masses = []
    for code, cls in zip(p.residues, p.classes):
        if cls in _DEFINED:
            masses.append(RESIDUES[code].mw_residue_avg)
        elif ambiguity == "mean-of-candidates":
            candidates = sorted(AMBIGUOUS[code])
            masses.append(math.fsum(RESIDUES[c].mw_residue_avg for c in candidates)
                          / len(candidates))
        else:
            return _feature("molecular_weight_avg", None, method,
                            "not computed: ambiguous residues are excluded", unit="Da")
    assumption = "unmodified residues; no disulfides"
    if "ambiguous" in p.classes:
        assumption += "; ambiguous residues use the candidate mean mass (X: mean of 20)"
    value = _round(math.fsum(masses + [WATER_AVG]))
    return _feature("molecular_weight_avg", value, method, assumption, unit="Da")


def _net_charge(p: Prepared) -> dict:
    low = high = 0
    unmodeled = []
    for code, cls, position in zip(p.residues, p.classes, p.input_positions):
        if cls == "canonical":
            low += RESIDUES[code].charge_simplified
            high += RESIDUES[code].charge_simplified
        elif cls == "ambiguous":
            charges = [RESIDUES[c].charge_simplified for c in AMBIGUOUS[code]]
            low += min(charges)
            high += max(charges)
        elif cls == "rare-but-defined":
            unmodeled.append(position)
    return _feature(
        "net_charge_simplified", low if low == high else None, "simplified",
        _CHARGE_ASSUMPTION + "; ambiguous residues add the [min, max] of their "
        "candidates (value is null when the range is wider than one number); "
        "U and O are not modeled and not counted",
        range=[low, high], unmodeled_positions=unmodeled)


def _positions(p: Prepared, codes: str) -> list[int]:
    return [pos for c, cls, pos in zip(p.residues, p.classes, p.input_positions)
            if cls == "canonical" and c in codes]


# ---------------------------------------------------------------- helpers

def _document(p: Prepared, ambiguity: str, params: dict, names, **body) -> dict:
    return {
        "schema": SCHEMA,
        "record": _record(p, ambiguity),
        "provenance": provenance(table_version=TABLE_VERSION, params=params),
        "field_classes": field_classes(names),
        **body,
        "limitations": list(_LIMITATIONS),
    }


def _record(p: Prepared, ambiguity: str) -> dict:
    ambiguous = [pos for pos, cls in zip(p.input_positions, p.classes) if cls == "ambiguous"]
    tokens = p.header.split()
    return {
        "id": tokens[0] if tokens else "",
        "header": p.header,
        "length": _length(p),
        "sequence_sha256": sequence_sha256(p.residues),
        "residue_classes": {cls: p.classes.count(cls)
                            for cls in ("canonical", "rare-but-defined", "ambiguous", "stop")},
        "ambiguity": {
            "ambiguity_policy": ambiguity,
            "ambiguous_count": len(ambiguous),
            "excluded_positions": ambiguous if ambiguity == "exclude" else [],
        },
        "modifications": {"stripped_terminal_stop": p.stripped_terminal_stop,
                          "stripped_gaps": p.stripped_gaps},
    }


def _feature(name: str, value, method: str, assumption: str = "", **extra) -> dict:
    return {"name": name, "value": value, "evidence_class": FIELD_CLASSES[name],
            "method": method, "assumption": assumption, **extra}


def _as_prepared(seq: str | Prepared) -> Prepared:
    if isinstance(seq, Prepared):
        return seq
    if isinstance(seq, str):
        return prepare("", seq)
    raise TypeError(f"expected str or Prepared, got {type(seq).__name__}")


def _check_ambiguity(p: Prepared, ambiguity: str) -> None:
    _check_choice("ambiguity", ambiguity, AMBIGUITY_POLICIES)
    if ambiguity == "error":
        for code, cls, position in zip(p.residues, p.classes, p.input_positions):
            if cls == "ambiguous":
                raise SequenceError("ambiguous symbol (ambiguity policy 'error')",
                                    header=p.header, position=position, symbol=code)


def _check_choice(name: str, value: str, allowed: tuple[str, ...]) -> None:
    if value not in allowed:
        raise ValueError(f"{name} must be one of {', '.join(allowed)}; got {value!r}")


def _is_int(value) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _round(value) -> float | None:
    if value is None:
        return None
    rounded = round(float(value), DECIMALS)
    return 0.0 if rounded == 0 else rounded  # no "-0.0" in output
