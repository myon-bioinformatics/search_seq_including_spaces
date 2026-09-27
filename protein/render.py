"""Plain-text rendering of a regions_of_interest() document (output sink).

render_ascii() takes the result dict and returns a string; it imports
nothing from the package and does no science. Tracks, per block of
`width` residues:

    ruler   position every 10 residues, printed only if it fits the block
    SEQ     residues
    CHG<w>  window net charge: + >= +T, - <= -T, . otherwise, ? no value
    HYD<w>  window mean KD:    # >= T, + 0 to T, . < 0, ? no value
    ROI     [====] per region, with the region id on the line(s) below

Every line fits in width + 8 characters. No colour and no escape codes;
with charset "ascii" every character is ASCII.
"""

import bisect
import math
import textwrap

CHARSETS = ("ascii",)
MIN_WIDTH = 20
MARGIN = 8  # track-name column


def render_ascii(result: dict, *, width: int = 60, charset: str = "ascii") -> str:
    if charset not in CHARSETS:
        raise ValueError(f"charset must be one of {', '.join(CHARSETS)}; got {charset!r}")
    if isinstance(width, bool) or not isinstance(width, int) or width < MIN_WIDTH:
        raise ValueError(f"width must be an integer >= {MIN_WIDTH}, got {width!r}")
    limit = width + MARGIN
    lines = _header(result, limit)
    rows = result["residues"]
    ends = [region["end"] for region in result["regions"]]
    for block in range(0, len(rows), width):
        lines.append("")
        lines.extend(_block(result, block, min(block + width, len(rows)), ends))
    for region in result["regions"]:
        lines.append("")
        lines.extend(_frame(result, region, limit))
    return "\n".join(_ascii(line.rstrip()) for line in lines) + "\n"


def _header(result: dict, limit: int) -> list[str]:
    params = result["provenance"]["params"]
    roi = result["provenance"]["roi"]
    summary = result["summary"]
    changed = sorted(k for k, s in roi["threshold_sources"].items() if s == "user")
    thresholds = ("thresholds=defaults" if not changed else
                  "thresholds=user-set " + ", ".join(f"{k}={roi['thresholds'][k]}" for k in changed))
    charge_t = _number(roi["thresholds"]["charge_cluster"])
    kd_t = _number(roi["thresholds"]["hydropathy_peak"])
    w = params["window"]
    texts = [
        "ROI: heuristic, sequence-only; not a functional/binding site",
        f"method={roi['roi_method']} window={w} half={roi['half_window']} "
        f"entropy_w={params['entropy_window']} merge_gap={roi['merge_gap']} "
        f"min_len={roi['min_len']} ambiguity={params['ambiguity']}",
        f"{thresholds} (all values and sources in JSON provenance)",
        f"{_ascii(result['record']['id']) or '(no id)'} len={result['record']['length']} "
        f"ROIs={summary['roi_count']} coverage={summary['roi_coverage']:.0%}",
        f"CHG{w}  + net >= +{charge_t}   - net <= -{charge_t}   . otherwise   ? n/a",
        f"HYD{w}  # >= {kd_t}   + 0 to {kd_t}   . < 0   ? n/a",
        *(f"WARNING: {text}" for text in summary["warnings"]),
    ]
    lines = []
    for text in texts:
        lines.extend(textwrap.wrap(text, width=limit, initial_indent="# ",
                                   subsequent_indent="#   ", break_on_hyphens=False))
    return lines


def _block(result: dict, lo: int, hi: int, ends: list[int]) -> list[str]:
    rows = result["residues"][lo:hi]
    size = hi - lo
    w = result["provenance"]["params"]["window"]
    charge_t = result["provenance"]["roi"]["thresholds"]["charge_cluster"]
    kd_t = result["provenance"]["roi"]["thresholds"]["hydropathy_peak"]

    ruler, ticks = [" "] * size, [" "] * size
    for col in range(0, size, 10):
        ticks[col] = "|"
        label = str(lo + col + 1)
        if col + len(label) <= size:
            ruler[col:col + len(label)] = label

    def chg(v):
        return "?" if v is None else "+" if v >= charge_t else "-" if v <= -charge_t else "."

    def hyd(v):
        return "?" if v is None else "#" if v >= kd_t else "+" if v >= 0 else "."

    roi = [" "] * size
    labels = []
    for region in _overlapping(result["regions"], ends, lo, hi):
        start, end = region["start"] - 1, region["end"] - 1
        a, b = max(start, lo) - lo, min(end, hi - 1) - lo
        roi[a:b + 1] = "=" * (b - a + 1)
        if start >= lo:
            roi[a] = "["
        if end < hi:
            roi[b] = "]"
        text = ("<- " if start < lo else "") + region["region_id"] + (" ->" if end >= hi else "")
        labels.append((a, text))

    lines = [
        "".ljust(MARGIN) + "".join(ruler),
        "".ljust(MARGIN) + "".join(ticks),
        "SEQ".ljust(MARGIN) + "".join(r["residue"] for r in rows),
        f"CHG{w}".ljust(MARGIN) + "".join(chg(r["charge_window"]) for r in rows),
        f"HYD{w}".ljust(MARGIN) + "".join(hyd(r["hydropathy_window"]) for r in rows),
        "ROI".ljust(MARGIN) + "".join(roi),
    ]
    lines.extend("".ljust(MARGIN) + row for row in _label_rows(labels, size))
    return lines


def _overlapping(regions: list[dict], ends: list[int], lo: int, hi: int) -> list[dict]:
    """Regions overlapping 0-based [lo, hi); regions are sorted and disjoint."""
    found = []
    for region in regions[bisect.bisect_left(ends, lo + 1):]:
        if region["start"] - 1 >= hi:
            break
        found.append(region)
    return found


def _label_rows(labels: list[tuple[int, str]], size: int) -> list[str]:
    """Place each label at its column (shifted left to fit), one space apart."""
    rows: list[list[str]] = []
    ends: list[int] = []
    for col, text in labels:
        if len(text) > size:
            text = text[:max(0, size - 1)] + "~"
        col = max(0, min(col, size - len(text)))
        k = next((k for k, end in enumerate(ends) if col > end), None)
        if k is None:
            rows.append([" "] * size)
            ends.append(-1)
            k = len(rows) - 1
        rows[k][col:col + len(text)] = text
        ends[k] = col + len(text)
    return ["".join(row) for row in rows]


def _frame(result: dict, region: dict, limit: int) -> list[str]:
    w = result["provenance"]["params"]["window"]
    ctx = region["context"]
    more_before = region["start"] - 1 > len(ctx["before"])
    more_after = len(result["residues"]) - region["end"] > len(ctx["after"])
    context = ("..." if more_before else "") + f"{ctx['before']}[{ctx['region']}]{ctx['after']}" \
        + ("..." if more_after else "")
    net = region["net_charge_simplified"]
    detail = (f"net charge {_signed(net['value'])} (simplified): "
              f"D/E={net['negative']}, K/R={net['positive']}")
    if net["not_counted"]:
        detail += f"; {net['not_counted']} residue(s) without a charge value not counted"
    edges = f"+/-{region['boundary_uncertainty']} (window {w})"
    if region["edge_truncated"] in ("start", "both"):
        edges += f"; left edge truncated at {region['start']}"
    if region["edge_truncated"] in ("end", "both"):
        edges += f"; right edge truncated at {region['end']}"

    lines = textwrap.wrap(
        f"{region['region_id']}  {region['start']}-{region['end']}  len={region['length']}  "
        f"class={region['evidence_class']}  mode=simple", width=limit, subsequent_indent="    ")
    lines += _field("context", context, limit)
    for k, trigger in enumerate(region["triggers"]):
        filled = math.ceil(10 * trigger["residues_triggered"] / region["length"])
        bar = ("#" * filled).ljust(10)
        lines += _field("reason" if k == 0 else "",
                        f"{trigger['feature']}  {bar}  "
                        f"{trigger['residues_triggered']}/{region['length']} residues", limit)
        lines += _field("", f"rule {trigger['rule']} ({trigger['threshold_source']})", limit)
    lines += _field("detail", detail, limit)
    lines += _field("edges", edges, limit)
    return lines


def _field(name: str, text: str, limit: int) -> list[str]:
    head = f"  {name:<7}  "
    return textwrap.wrap(text, width=limit, initial_indent=head,
                         subsequent_indent=" " * len(head), break_on_hyphens=False)


def _signed(value: int) -> str:
    return f"{value:+d}" if value else "0"


def _number(value) -> str:
    return f"{value:g}"


def _ascii(line: str) -> str:
    return line.encode("ascii", "backslashreplace").decode("ascii")
