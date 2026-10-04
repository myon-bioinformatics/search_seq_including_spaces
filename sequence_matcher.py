"""Lazy reader-to-matcher adapter; see docs/SEQUENCE_MATCHER.md."""
from contextlib import closing
from dataclasses import asdict

from anchor_matcher import find_matches
from sequence_readers import read_sequences

SCHEMA = "sequence-matches/1"
__all__ = ["SCHEMA", "find_sequence_matches"]


def find_sequence_matches(input, left, right, *, format=None, source=None,
                          sequence_column=None, id_column=None, name_column=None,
                          description_column=None, min_gap=0, max_gap=None,
                          iupac=False, strand="+", max_matches=10000):
    """Yield one metadata envelope and unchanged anchor-gap/1 result per record.

    Reader options go to read_sequences; matcher options go to find_matches.
    Records are independent, including records with no matches. max_matches is
    a per-record limit, not a total budget; overflow raises without yielding
    a partial result for that record. Reader/matcher/I/O errors propagate, so
    previously yielded envelopes do not establish whole-input completeness.
    Work and validation are lazy; matcher validation occurs on each record.

    source_span and middle_source_span remain character positions within
    record.sequence, NOT file coordinates. Reader provenance stays in the
    outer record metadata; raw sequence and FASTQ quality are not copied there.

    Memory is proportional to the current record plus its matches (including
    reader buffers and matcher normalization/workspace), not constant for a
    huge single record. No earlier records/results are accumulated here.
    Exhaust or close this generator to release reader-owned files; use
    contextlib.closing when breaking early. Caller-owned streams stay open.
    """
    with closing(read_sequences(
            input, format=format, source=source, sequence_column=sequence_column,
            id_column=id_column, name_column=name_column,
            description_column=description_column)) as records:
        for record in records:
            yield {
                "schema": SCHEMA,
                "record": {
                    "id": record.id,
                    "description": record.description,
                    "source": record.source,
                    "provenance": asdict(record.provenance),
                },
                "result": find_matches(
                    record.sequence, left, right, min_gap=min_gap, max_gap=max_gap,
                    iupac=iupac, strand=strand, max_matches=max_matches),
            }
            # Release the previous sequence before asking for the next record.
            del record
