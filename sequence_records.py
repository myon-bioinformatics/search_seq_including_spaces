"""Minimal stdlib reader contract; independent of matching and analysis."""
from dataclasses import dataclass

SCHEMA = "sequence-record/1"
__all__ = ["SCHEMA", "SequenceRecord", "ReadProvenance"]


@dataclass(frozen=True)
class ReadProvenance:
    format: str
    compression: str | None
    record_index: int
    start_line: int
    end_line: int


@dataclass(frozen=True)
class SequenceRecord:
    id: str
    description: str
    sequence: str
    source: str | None
    provenance: ReadProvenance
    quality: str | None = None
