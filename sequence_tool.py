#!/usr/bin/env python3
"""Command-line entry point for the sequence tools (stdlib only).

    python sequence_tool.py compare FILE_A FILE_B --alphabet {dna,rna,protein} [--json]
    python sequence_tool.py fasta validate FILE [--strip-gaps] [--allow-internal-stop] [--json]
    python sequence_tool.py protein features FILE [--ambiguity POLICY] [--json | --jsonl]
    python sequence_tool.py protein profile FILE [--window 9] [--entropy-window 12]
                                            [--ambiguity POLICY] [--tsv | --json | --jsonl]
    python sequence_tool.py protein roi FILE [--window 9] [--merge-gap 2] [--min-len 5]
                                        [--threshold NAME=VALUE ...] [--json | --jsonl | --ascii]
                                        [--width 60]

FILE may be "-" for standard input. Each subcommand only parses arguments,
calls one function of the protein package and writes the result; the
analysis itself lives in the package. Exit codes: 0 = ok, 2 = input error
(unreadable file, invalid sequence, bad option).

JSON output is UTF-8 with sorted keys, so the same input gives the same
bytes. --json writes one document and needs exactly one FASTA record;
--jsonl writes one compact document per record.
"""

import argparse
import json
import math
import sys

from comparison import (ALPHABETS, AMBIGUITY_POLICIES as COMPARE_AMBIGUITY_POLICIES,
                        GAP_POLICIES, ComparisonError, compare_aligned)
from comparison_s2 import render_ascii_diff, sliding_identity
from protein.features import AMBIGUITY_POLICIES, protein_features, protein_profiles
from protein.parse import SequenceError, prepare, read_fasta
from protein.render import MIN_WIDTH, render_ascii
from protein.roi import DEFAULT_THRESHOLDS, regions_of_interest

EXIT_OK = 0
EXIT_INPUT_ERROR = 2
PROFILE_COLUMNS = ("position", "input_position", "residue", "class", "hydropathy", "charge",
                   "hydropathy_window", "charge_window", "entropy_window")


class InputError(Exception):
    """A problem with the user's input, reported on stderr with exit code 2."""


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        output, status = args.handler(args)
    except (InputError, SequenceError, ComparisonError) as error:
        print(f"{parser.prog}: error: {error}", file=sys.stderr)
        return EXIT_INPUT_ERROR
    sys.stdout.buffer.write(output.encode("utf-8"))
    sys.stdout.flush()
    return status


# ---------------------------------------------------------------- subcommands

def _compare(args) -> tuple[str, int]:
    if args.file_a == "-" and args.file_b == "-":
        raise InputError("compare cannot read both inputs from stdin")
    if args.ascii and args.window is not None:
        raise InputError("--ascii and --window cannot be used together")
    if args.ascii and args.step != 1:
        raise InputError("--step is only meaningful with --window")
    if not args.ascii and args.width != 60:
        raise InputError("--width is only meaningful with --ascii")
    if args.window is None and args.step != 1:
        raise InputError("--step is only meaningful with --window")
    header_a, raw_a = _one_raw_record(args.file_a)
    header_b, raw_b = _one_raw_record(args.file_b)
    if args.ascii:
        body = render_ascii_diff(
            raw_a, raw_b, alphabet=args.alphabet, gap_policy=args.gap_policy,
            ambiguity_policy=args.ambiguity_policy, start=args.start, end=args.end,
            width=args.width
        )
        record_a = header_a.split()[0] if header_a.split() else ""
        record_b = header_b.split()[0] if header_b.split() else ""
        context = (
            f"records: a={record_a} b={record_b}  alphabet={args.alphabet}  "
            f"gap={args.gap_policy}  ambiguity={args.ambiguity_policy}\n"
        )
        return context + body, EXIT_OK
    if args.window is not None:
        doc = sliding_identity(
            raw_a, raw_b, alphabet=args.alphabet, gap_policy=args.gap_policy,
            ambiguity_policy=args.ambiguity_policy, start=args.start, end=args.end,
            window=args.window, step=args.step
        )
    else:
        doc = compare_aligned(
            raw_a, raw_b, alphabet=args.alphabet, gap_policy=args.gap_policy,
            ambiguity_policy=args.ambiguity_policy, start=args.start, end=args.end
        )
    doc["records"] = {
        "a": {"id": header_a.split()[0] if header_a.split() else ""},
        "b": {"id": header_b.split()[0] if header_b.split() else ""},
    }
    return _dump(doc, indent=2) + "\n", EXIT_OK


def _fasta_validate(args) -> tuple[str, int]:
    """Report every record; exit 2 if any is invalid (errors go to stderr)."""
    rows, reports, failed = [], [], False
    for header, raw in _read_records(args.file):
        record_id = header.split()[0] if header.split() else ""
        try:
            p = prepare(header, raw, strip_gaps=args.strip_gaps,
                        allow_internal_stop=args.allow_internal_stop)
        except SequenceError as error:
            failed = True
            print(error, file=sys.stderr)
            reports.append({"id": record_id, "valid": False,
                            "error": {"reason": error.reason, "position": error.position,
                                      "symbol": error.symbol}})
            continue
        counts = {cls: p.classes.count(cls)
                  for cls in ("canonical", "rare-but-defined", "ambiguous", "stop")}
        length = len(p.residues) - counts["stop"]
        rows.append([record_id, length, *counts.values(),
                     str(p.stripped_terminal_stop).lower(), p.stripped_gaps])
        reports.append({"id": record_id, "valid": True, "length": length,
                        "residue_classes": counts,
                        "modifications": {"stripped_terminal_stop": p.stripped_terminal_stop,
                                          "stripped_gaps": p.stripped_gaps}})
    if args.json:
        output = _dump({"valid": not failed, "records": reports}, indent=2)
    else:
        header = ["id", "length", "canonical", "rare_but_defined", "ambiguous", "stop",
                  "stripped_terminal_stop", "stripped_gaps"]
        output = "".join("\t".join(map(str, row)) + "\n" for row in [header, *rows])
    return output, EXIT_INPUT_ERROR if failed else EXIT_OK


def _protein_features(args) -> tuple[str, int]:
    docs = [protein_features(p, ambiguity=args.ambiguity) for p in _prepared(args)]
    return _write_documents(docs, args.format), EXIT_OK


def _protein_profile(args) -> tuple[str, int]:
    docs = [protein_profiles(p, window=args.window, entropy_window=args.entropy_window,
                             ambiguity=args.ambiguity) for p in _prepared(args)]
    if args.format != "tsv":
        return _write_documents(docs, args.format), EXIT_OK
    lines = ["\t".join(("id",) + PROFILE_COLUMNS)]
    for doc in docs:
        for row in doc["residues"]:
            lines.append("\t".join([doc["record"]["id"]]
                                   + [_tsv_cell(row[c]) for c in PROFILE_COLUMNS]))
    return "\n".join(lines) + "\n", EXIT_OK


def _protein_roi(args) -> tuple[str, int]:
    results = [
        regions_of_interest(
            protein_profiles(p, window=args.window, entropy_window=args.entropy_window,
                             ambiguity=args.ambiguity),
            thresholds=dict(args.threshold), merge_gap=args.merge_gap, min_len=args.min_len)
        for p in _prepared(args)
    ]
    if args.format == "ascii":
        return "\n".join(render_ascii(r, width=args.width) for r in results), EXIT_OK
    return _write_documents(results, args.format), EXIT_OK


# ---------------------------------------------------------------- helpers

def _read_records(path: str) -> list[tuple[str, str]]:
    try:
        if path == "-":
            text = sys.stdin.buffer.read().decode("utf-8")
        else:
            with open(path, encoding="utf-8") as handle:
                text = handle.read()
    except OSError as error:
        raise InputError(f"cannot read {path}: {error.strerror or error}") from None
    except UnicodeDecodeError as error:
        raise InputError(f"{path} is not UTF-8 text ({error.reason} at byte {error.start})") from None
    records = read_fasta(text)
    if not records:
        raise InputError(f"no FASTA records in {path}")
    return records


def _one_raw_record(path: str) -> tuple[str, str]:
    records = _read_records(path)
    if len(records) != 1:
        raise InputError(f"compare expects exactly one record in {path}; found {len(records)}")
    return records[0]


def _prepared(args) -> list:
    return [prepare(header, raw, strip_gaps=args.strip_gaps,
                    allow_internal_stop=args.allow_internal_stop)
            for header, raw in _read_records(args.file)]


def _write_documents(docs: list[dict], fmt: str) -> str:
    if fmt == "jsonl":
        return "".join(_dump(doc) + "\n" for doc in docs)
    if len(docs) != 1:
        raise InputError(f"--json writes one document but the input has {len(docs)} "
                         "records; use --jsonl")
    return _dump(docs[0], indent=2) + "\n"


def _dump(obj, indent: int | None = None) -> str:
    separators = None if indent else (",", ":")
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, indent=indent,
                      separators=separators, allow_nan=False)


def _tsv_cell(value) -> str:
    return "NA" if value is None else str(value)


def _odd_positive_int(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {text!r}") from None
    if value < 1 or value % 2 == 0:
        raise argparse.ArgumentTypeError(f"must be a positive odd integer: {value}")
    return value


def _int_at_least(lowest: int):
    def parse(text: str) -> int:
        try:
            value = int(text)
        except ValueError:
            raise argparse.ArgumentTypeError(f"not an integer: {text!r}") from None
        if value < lowest:
            raise argparse.ArgumentTypeError(f"must be an integer >= {lowest}: {value}")
        return value
    return parse


def _threshold(text: str) -> tuple[str, int | float]:
    name, sep, raw = text.partition("=")
    if not sep or name not in DEFAULT_THRESHOLDS:
        raise argparse.ArgumentTypeError(
            f"expected NAME=VALUE with NAME one of {', '.join(DEFAULT_THRESHOLDS)}: {text!r}")
    try:
        value = int(raw)
    except ValueError:
        try:
            value = float(raw)
        except ValueError:
            raise argparse.ArgumentTypeError(f"not a number: {raw!r}") from None
    if not math.isfinite(value):
        raise argparse.ArgumentTypeError(f"not a finite number: {raw!r}")
    return name, value


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="sequence_tool.py",
                                     description="Sequence tools (stdlib only).")
    groups = parser.add_subparsers(dest="group", required=True,
                                   metavar="{compare,fasta,protein}")

    compare = groups.add_parser(
        "compare", help="aligned pairwise identity/coverage (does not align sequences)")
    compare.add_argument("file_a", help='first aligned FASTA/plain sequence, or "-"')
    compare.add_argument("file_b", help="second aligned FASTA/plain sequence")
    compare.add_argument("--alphabet", choices=ALPHABETS, required=True,
                         help="sequence alphabet; no auto-detection")
    compare.add_argument("--gap-policy", choices=GAP_POLICIES, default="exclude",
                         help="exclude any gap column or count one-sided gaps as mismatch")
    compare.add_argument("--ambiguity-policy", choices=COMPARE_AMBIGUITY_POLICIES,
                         default="exclude",
                         help="exclude, literal strict, or compatible candidate-set handling")
    compare.add_argument("--start", type=_int_at_least(1), default=1,
                         help="1-based inclusive aligned start (default: 1)")
    compare.add_argument("--end", type=_int_at_least(1),
                         help="1-based inclusive aligned end (default: final column)")
    compare.add_argument("--window", type=_int_at_least(1),
                         help="S2 sliding-window size; emits only complete windows")
    compare.add_argument("--step", type=_int_at_least(1), default=1,
                         help="S2 sliding-window step (default: 1)")
    compare.add_argument("--ascii", action="store_true",
                         help="render an aligned ASCII diff instead of JSON")
    compare.add_argument("--width", type=_int_at_least(10), default=60,
                         help="aligned columns per block for --ascii (default: 60)")
    compare.add_argument("--json", action="store_true",
                         help="JSON output (default unless --ascii)")
    compare.set_defaults(handler=_compare)

    input_options = argparse.ArgumentParser(add_help=False)
    input_options.add_argument("file", help='FASTA file, or "-" for standard input')
    input_options.add_argument("--strip-gaps", action="store_true",
                               help="remove '-' gap symbols instead of rejecting them")
    input_options.add_argument("--allow-internal-stop", action="store_true",
                               help="keep internal '*' as a stop marker instead of rejecting it")

    fasta = groups.add_parser("fasta", help="FASTA checks")
    fasta_commands = fasta.add_subparsers(dest="command", required=True)
    validate = fasta_commands.add_parser("validate", parents=[input_options],
                                         help="validate protein FASTA records")
    validate.add_argument("--json", action="store_true", help="JSON report instead of TSV")
    validate.set_defaults(handler=_fasta_validate)

    protein = groups.add_parser("protein", help="protein sequence analysis")
    protein_commands = protein.add_subparsers(dest="command", required=True)
    analysis_options = argparse.ArgumentParser(add_help=False, parents=[input_options])
    analysis_options.add_argument("--ambiguity", choices=AMBIGUITY_POLICIES, default="exclude",
                                  help="how B, Z, J, X are treated (default: exclude)")

    features = protein_commands.add_parser("features", parents=[analysis_options],
                                           help="whole-sequence features (JSON)")
    _format_options(features, ("json", "jsonl"), default="json")
    features.set_defaults(handler=_protein_features)

    window_options = argparse.ArgumentParser(add_help=False, parents=[analysis_options])
    window_options.add_argument("--window", type=_odd_positive_int, default=9,
                                help="hydropathy and charge window, odd (default: 9)")
    window_options.add_argument("--entropy-window", type=_int_at_least(1), default=12,
                                help="entropy window (default: 12)")

    profile = protein_commands.add_parser("profile", parents=[window_options],
                                          help="per-residue window profiles")
    _format_options(profile, ("tsv", "json", "jsonl"), default="tsv")
    profile.set_defaults(handler=_protein_profile)

    roi = protein_commands.add_parser(
        "roi", parents=[window_options],
        help="regions of interest (heuristic; not functional or binding sites)")
    roi.add_argument("--merge-gap", type=_int_at_least(0), default=2,
                     help="join regions separated by at most this many residues (default: 2)")
    roi.add_argument("--min-len", type=_int_at_least(1), default=5,
                     help="drop regions shorter than this (default: 5)")
    roi.add_argument("--threshold", type=_threshold, action="append", default=[],
                     metavar="NAME=VALUE",
                     help="override a trigger threshold, e.g. charge_cluster=4 (repeatable)")
    roi.add_argument("--width", type=_int_at_least(MIN_WIDTH), default=60,
                     help=f"residues per line for --ascii (default: 60, minimum {MIN_WIDTH})")
    _format_options(roi, ("json", "jsonl", "ascii"), default="json")
    roi.set_defaults(handler=_protein_roi)
    return parser


def _format_options(parser: argparse.ArgumentParser, formats: tuple[str, ...],
                    default: str) -> None:
    group = parser.add_mutually_exclusive_group()
    for fmt in formats:
        group.add_argument(f"--{fmt}", dest="format", action="store_const", const=fmt,
                           help=f"{fmt.upper()} output" + (" (default)" if fmt == default else ""))
    parser.set_defaults(format=default)


if __name__ == "__main__":
    sys.exit(main())
