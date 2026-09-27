"""Regions of interest by absolute thresholds (layer L3, stdlib only).

A region of interest (ROI) is a run of residues where at least one
sequence-derived value crosses a fixed threshold (mode "simple"). It points
to where to look first and always carries its reasons. It is a heuristic:
not a functional, binding, interface ("hotspot" in the Bogan & Thorn sense)
or drug site. No score is given in this mode.

Input is a protein_profiles() document, so the values used here are the
values shown there. Triggers (defaults in DEFAULT_THRESHOLDS):

    charge_cluster       |charge_window| >= 3
    charge_transition    net charge of the `half_window` residues on the
                         left and on the right have opposite signs and
                         differ by >= 4
    hydropathy_peak      hydropathy_window >= 1.5
    hydropathy_gradient  |mean KD right - mean KD left| >= 3.0 (half windows)
    cys_cluster          Cys in the window >= 2
    aromatic_cluster     F/W/Y in the window >= 3
    pro_gly_cluster      P + G in the window >= 4
    low_complexity       entropy_window <= 2.2 bits (SEG-like, not SEG)

A position whose value is null (ends of the sequence, too few usable
residues) never triggers. Triggered positions are joined when at most
`merge_gap` untriggered residues separate them, runs shorter than
`min_len` are dropped, and a stop marker always ends a run.
"""

import math
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType

from .provenance import SCHEMA, field_classes
from .residues import RESIDUES

ROI_METHOD = "simple-threshold/1"
LABEL = "physicochemically notable region (not a functional-site prediction)"
CONTEXT = 6
COVERAGE_WARNING = 0.5


@dataclass(frozen=True)
class Threshold:
    value: float
    op: str      # ">=" or "<="
    source: str
    rule: str    # formatted with value, window, half_window, entropy_window


DEFAULT_THRESHOLDS = MappingProxyType({
    "charge_cluster": Threshold(3, ">=", "tool-default",
                                "|net charge, w={window}| >= {value}"),
    "charge_transition": Threshold(4, ">=", "tool-default",
                                   "net charge of {half_window} left vs {half_window} right: "
                                   "opposite signs, difference >= {value}"),
    "hydropathy_peak": Threshold(1.5, ">=", "tool-default",
                                 "mean KD, w={window} >= {value}"),
    "hydropathy_gradient": Threshold(3.0, ">=", "tool-default",
                                     "|mean KD {half_window} right - {half_window} left| >= {value}"),
    "cys_cluster": Threshold(2, ">=", "tool-default", "Cys in w={window} >= {value}"),
    "aromatic_cluster": Threshold(3, ">=", "tool-default", "F/W/Y in w={window} >= {value}"),
    "pro_gly_cluster": Threshold(4, ">=", "tool-default", "P+G in w={window} >= {value}"),
    "low_complexity": Threshold(2.2, "<=", "seg-like-default",
                                "entropy, w={entropy_window} <= {value} bits"),
})

_LIMITATIONS = (
    "ROI is a heuristic from the sequence alone; it is not a functional, binding, "
    "interface or drug site",
    "thresholds are tool defaults; results change strongly with thresholds, window "
    "and merge gap",
    "window smoothing blurs region boundaries by about half a window",
)


def regions_of_interest(profiles: dict, *, thresholds: dict | None = None,
                        merge_gap: int = 2, min_len: int = 5,
                        half_window: int = 5) -> dict:
    """Detect ROIs in a protein_profiles() document; return a new document.

    `thresholds` overrides default values by trigger name, for example
    {"charge_cluster": 4}. The result is the profile document plus
    "regions", "summary" and provenance["roi"].
    """
    _check_profiles(profiles)
    for name, value, lowest in (("merge_gap", merge_gap, 0), ("min_len", min_len, 1),
                                ("half_window", half_window, 1)):
        if not _is_int(value) or value < lowest:
            raise ValueError(f"{name} must be an integer >= {lowest}, got {value!r}")
    active = _thresholds(thresholds or {})
    params = profiles["provenance"]["params"]
    rows = profiles["residues"]
    sizes = {"window": params["window"], "half_window": half_window,
             "entropy_window": params["entropy_window"]}

    hits = _trigger_masks(rows, active, sizes)
    stops = [row["class"] == "stop" for row in rows]
    mask = [any(values[i] for values in hits.values()) for i in range(len(rows))]
    runs = _runs(mask, stops, merge_gap, min_len)

    first, last = _evaluable_range(len(rows), sizes)
    digits = max(2, len(str(len(runs))))
    regions = [_region(f"ROI-{k:0{digits}d}", start, end, rows, hits, active, sizes, first, last)
               for k, (start, end) in enumerate(runs, start=1)]

    residue_count = len(rows) - sum(stops)
    covered = sum(r["length"] for r in regions)
    coverage = round(covered / residue_count, 4) if residue_count else 0.0
    warnings = []
    if coverage > COVERAGE_WARNING:
        warnings.append(f"ROI coverage {coverage:.0%} is above 50%: regions are not "
                        "selective; raise the thresholds")

    # New containers for what changes; the residue rows are shared, not modified.
    result = dict(profiles)
    result["provenance"] = dict(profiles["provenance"])
    result["field_classes"] = dict(profiles["field_classes"])
    result["provenance"]["roi"] = {
        "roi_method": ROI_METHOD,
        "threshold_method": "absolute",
        "thresholds": {name: t.value for name, t in active.items()},
        "threshold_ops": {name: t.op for name, t in active.items()},
        "threshold_sources": {name: t.source for name, t in active.items()},
        "half_window": half_window,
        "merge_gap": merge_gap,
        "min_len": min_len,
        "weights": None,
        "top_regions": None,
    }
    result["field_classes"].update(field_classes(["regions", "summary"]))
    result["regions"] = regions
    result["summary"] = {"roi_count": len(regions), "roi_coverage": coverage,
                         "warnings": warnings}
    result["limitations"] = list(profiles["limitations"]) + list(_LIMITATIONS)
    return result


# ---------------------------------------------------------------- triggers

def _trigger_masks(rows: list[dict], active: dict, sizes: dict) -> dict[str, list[bool]]:
    n = len(rows)
    window, half = sizes["window"], sizes["half_window"]
    r = window // 2
    codes = [row["residue"] for row in rows]
    known = [row["class"] in ("canonical", "rare-but-defined") for row in rows]
    charge = _Prefix([row["charge"] for row in rows])
    kd = _Prefix([row["hydropathy"] for row in rows])
    counts = {
        "cys_cluster": _Prefix([int(c == "C") if k else None for c, k in zip(codes, known)]),
        "aromatic_cluster": _Prefix([int(RESIDUES[c].aromatic) if k else None
                                     for c, k in zip(codes, known)]),
        "pro_gly_cluster": _Prefix([int(c in "PG") if k else None
                                    for c, k in zip(codes, known)]),
    }
    columns = {"charge_cluster": "charge_window", "hydropathy_peak": "hydropathy_window",
               "low_complexity": "entropy_window"}

    def exact(name: str, i: int) -> tuple[int, int] | None:
        """Value as an exact ratio (numerator, positive denominator), or None."""
        if name in ("charge_transition", "hydropathy_gradient"):
            prefix = charge if name == "charge_transition" else kd
            left, right = prefix.ratio(i - half, i), prefix.ratio(i + 1, i + 1 + half)
            if left is None or right is None:
                return None
            (a, b), (c, d) = left, right
            if name == "charge_transition":
                if a * c >= 0:
                    return None
                return half * abs(a * d - c * b), b * d
            return abs(c * b - a * d), b * d
        ratio = counts[name].ratio(i - r, i + r + 1)
        return None if ratio is None else (ratio[0] * window, ratio[1])

    masks = {}
    for name, t in active.items():
        at_least = t.op == ">="
        if name in columns:
            # Profile columns are compared as shown in the profile.
            key, bound = columns[name], float(t.value)
            values = (row[key] for row in rows)
            if name == "charge_cluster":
                values = (None if v is None else abs(v) for v in values)
            masks[name] = [v is not None and (v >= bound if at_least else v <= bound)
                           for v in values]
            continue
        limit = Fraction(str(t.value))
        mask = []
        for i in range(n):
            v = exact(name, i)
            if v is None:
                mask.append(False)
            else:
                lhs, rhs = v[0] * limit.denominator, limit.numerator * v[1]
                mask.append(lhs >= rhs if at_least else lhs <= rhs)
        masks[name] = mask
    return masks


class _Prefix:
    """Exact range means of numbers (int, float taken as its decimal, None skipped)."""

    def __init__(self, values: list):
        exact = {v: Fraction(str(v)) if isinstance(v, float) else Fraction(v)
                 for v in set(values) if v is not None}
        self.denominator = math.lcm(*(v.denominator for v in exact.values()))
        scaled = {v: int(f * self.denominator) for v, f in exact.items()}
        self.sums, self.counts = [0], [0]
        for v in values:
            self.sums.append(self.sums[-1] + (0 if v is None else scaled[v]))
            self.counts.append(self.counts[-1] + (v is not None))

    def ratio(self, lo: int, hi: int) -> tuple[int, int] | None:
        """Mean of values[lo:hi] as (numerator, denominator); None past either
        end or when fewer than half of the range is usable."""
        if lo < 0 or hi > len(self.sums) - 1:
            return None
        used = self.counts[hi] - self.counts[lo]
        if used < (hi - lo + 1) // 2 or used == 0:
            return None
        return self.sums[hi] - self.sums[lo], used * self.denominator

    def mean(self, lo: int, hi: int) -> Fraction | None:
        ratio = self.ratio(lo, hi)
        return None if ratio is None else Fraction(*ratio)


def _runs(mask: list[bool], breaks: list[bool], merge_gap: int,
          min_len: int) -> list[tuple[int, int]]:
    """0-based inclusive (start, end) runs of True, gaps <= merge_gap joined."""
    runs: list[tuple[int, int]] = []
    start = last = None
    for i, (hit, stop) in enumerate(zip(mask, breaks)):
        if start is not None and (stop or (hit and i - last - 1 > merge_gap)):
            runs.append((start, last))
            start = None
        if hit and not stop:
            if start is None:
                start = i
            last = i
    if start is not None:
        runs.append((start, last))
    return [(s, e) for s, e in runs if e - s + 1 >= min_len]


# ---------------------------------------------------------------- regions

def _region(region_id: str, start: int, end: int, rows: list[dict], hits: dict,
            active: dict, sizes: dict, first: int, last: int) -> dict:
    triggers = []
    for order, (name, mask) in enumerate(hits.items()):
        count = sum(mask[start:end + 1])
        if count:
            t = active[name]
            triggers.append((-count, order, {
                "feature": name,
                "residues_triggered": count,
                "threshold": f"{t.op}{t.value}",
                "threshold_source": t.source,
                "rule": t.rule.format(value=t.value, **sizes),
            }))
    codes = "".join(row["residue"] for row in rows[max(0, start - CONTEXT):end + 1 + CONTEXT])
    offset = max(0, start - CONTEXT)
    inside = rows[start:end + 1]
    negative = sum(row["charge"] == -1 for row in inside)
    positive = sum(row["charge"] == 1 for row in inside)
    at_start, at_end = start <= first, end >= last
    return {
        "region_id": region_id,
        "start": start + 1,
        "end": end + 1,
        "length": end - start + 1,
        "input_start": rows[start]["input_position"],
        "input_end": rows[end]["input_position"],
        "evidence_class": "heuristic",
        "detection": ROI_METHOD,
        "triggers": [entry for _, _, entry in sorted(triggers, key=lambda x: x[:2])],
        "score": None,
        "contributions": None,
        "boundary_uncertainty": sizes["window"] // 2,
        "edge_truncated": ("both" if at_start and at_end else "start" if at_start
                           else "end" if at_end else None),
        "context": {"before": codes[:start - offset],
                    "region": codes[start - offset:end + 1 - offset],
                    "after": codes[end + 1 - offset:]},
        "net_charge_simplified": {"value": positive - negative, "negative": negative,
                                  "positive": positive,
                                  "not_counted": sum(row["charge"] is None for row in inside)},
        "label": LABEL,
    }


def _evaluable_range(n: int, sizes: dict) -> tuple[int, int]:
    """0-based first and last positions where any trigger can be evaluated."""
    window, half, ew = sizes["window"], sizes["half_window"], sizes["entropy_window"]
    first = min(window // 2, half, ew // 2)
    last = max(n - 1 - window // 2, n - 1 - half, n - ew + ew // 2)
    return first, last


# ---------------------------------------------------------------- checks

def _check_profiles(profiles: dict) -> None:
    if not isinstance(profiles, dict) or profiles.get("schema") != SCHEMA \
            or "residues" not in profiles:
        raise ValueError("expected a protein_profiles() document")
    if profiles["provenance"]["params"].get("edge") != "none":
        raise ValueError("ROI detection supports profiles with edge='none' only")


def _thresholds(overrides: dict) -> dict[str, Threshold]:
    unknown = set(overrides) - set(DEFAULT_THRESHOLDS)
    if unknown:
        raise ValueError(f"unknown trigger(s): {', '.join(sorted(unknown))}; "
                         f"known: {', '.join(DEFAULT_THRESHOLDS)}")
    active = {}
    for name, default in DEFAULT_THRESHOLDS.items():
        if name not in overrides:
            active[name] = default
            continue
        value = overrides[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)) \
                or not math.isfinite(value):
            raise ValueError(f"threshold for {name} must be a finite number, got {value!r}")
        active[name] = Threshold(value, default.op, "user", default.rule)
    return active


def _is_int(value) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)
