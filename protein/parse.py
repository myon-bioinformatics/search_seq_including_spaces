"""FASTA reading and symbol validation (layer L1, stdlib only).

Positions are 1-based and refer to the input sequence as written, counting
every non-whitespace symbol (residues, gaps and stops alike). They stay
attached to each residue after gaps or a terminal stop are removed, so later
layers can always report positions in the user's own coordinates.
"""

from dataclasses import dataclass

from .residues import STOP, symbol_class


class SequenceError(ValueError):
    """Invalid or unusable sequence input.

    Attributes: reason (str), header (str), position (int | None, 1-based
    input position), symbol (str | None, the offending symbol as written).
    """

    def __init__(self, reason: str, *, header: str = "",
                 position: int | None = None, symbol: str | None = None):
        self.reason = reason
        self.header = header
        self.position = position
        self.symbol = symbol
        parts = [reason]
        if symbol is not None:
            parts.append(f"symbol {symbol!r}")
        if position is not None:
            parts.append(f"at position {position}")
        if header:
            parts.append(f"in record {header!r}")
        super().__init__(", ".join(parts))


@dataclass(frozen=True)
class Prepared:
    """A validated protein sequence.

    residues         upper-case symbols kept for analysis (canonical,
                     rare-but-defined, ambiguous; plus '*' only when
                     internal stops were explicitly allowed)
    input_positions  1-based input position of each kept symbol
    classes          symbol class of each kept symbol (see
                     residues.SYMBOL_CLASSES)
    """

    header: str
    residues: str
    input_positions: tuple[int, ...]
    classes: tuple[str, ...]
    stripped_terminal_stop: bool
    stripped_gaps: int


def read_fasta(text: str) -> list[tuple[str, str]]:
    """Split FASTA text into (header, raw_sequence) pairs.

    Header lines start with '>'; the header is the rest of the line,
    stripped. Blank lines and ';' comment lines are ignored. Text without
    any header is read as a single record with an empty header. Sequence
    text before the first header, when headers exist, is an error.
    Raw sequences are returned unvalidated; pass them to prepare().
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    preamble: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if stripped.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            elif preamble:
                raise SequenceError("sequence data before the first '>' header")
            header = stripped[1:].strip()
            chunks = []
        elif header is None:
            preamble.append(stripped)
        else:
            chunks.append(stripped)
    if header is not None:
        records.append((header, "".join(chunks)))
    elif preamble:
        records.append(("", "".join(preamble)))
    return records


def prepare(header: str, raw: str, *, strip_gaps: bool = False,
            allow_internal_stop: bool = False) -> Prepared:
    """Validate one raw sequence and return it as Prepared.

    - Whitespace is removed; symbols are upper-cased.
    - A single trailing '*' is removed and recorded.
    - An internal '*' is an error unless allow_internal_stop=True, in which
      case it is kept as a 'stop' marker for later layers to handle.
    - '-' is an error unless strip_gaps=True, in which case gaps are removed
      and counted.
    - Any other unknown symbol is an error reporting its input position.
    - Ambiguity codes (B, Z, J, X) are kept with class 'ambiguous'; turning
      them into numbers is a later-layer policy (default: exclude).
    """
    symbols = [c for c in raw if not c.isspace()]
    terminal_stop = bool(symbols) and symbols[-1] == STOP
    end = len(symbols) - 1 if terminal_stop else len(symbols)

    residues: list[str] = []
    positions: list[int] = []
    classes: list[str] = []
    gaps = 0
    for index in range(end):
        written = symbols[index]
        position = index + 1
        cls = symbol_class(written)
        if cls == "invalid":
            raise SequenceError("invalid symbol", header=header,
                                position=position, symbol=written)
        if cls == "gap":
            if strip_gaps:
                gaps += 1
                continue
            raise SequenceError("gap symbol (use strip_gaps=True to remove gaps)",
                                header=header, position=position, symbol=written)
        if cls == "stop" and not allow_internal_stop:
            raise SequenceError("internal stop (use allow_internal_stop=True to keep it)",
                                header=header, position=position, symbol=written)
        residues.append(written.upper())
        positions.append(position)
        classes.append(cls)

    if all(cls == "stop" for cls in classes):  # also true when nothing is left
        raise SequenceError("empty sequence", header=header)

    return Prepared(
        header=header,
        residues="".join(residues),
        input_positions=tuple(positions),
        classes=tuple(classes),
        stripped_terminal_stop=terminal_stop,
        stripped_gaps=gaps,
    )
