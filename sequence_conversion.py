"""Semantic FASTA <-> TSV text conversion; no byte/layout round-trip claim."""
import csv

from sequence_readers import read_sequences, SequenceReadError

__all__ = ['fasta_to_text', 'text_to_fasta']


def fasta_to_text(input, output, *, source=None):
    """Write id/description/sequence TSV to a caller-owned text stream.

    Input accepts the reader's paths (including compression) or text streams.
    Open output with newline='' to preserve csv quoting on every platform.
    Headers, sequence case/symbols and record order survive the reverse helper;
    comments, wrapping, source/provenance and original newline bytes do not.
    """
    writer = csv.writer(output, delimiter='\t', lineterminator='\n')
    writer.writerow(['id', 'description', 'sequence'])
    for record in read_sequences(input, format='fasta', source=source):
        writer.writerow([record.id, record.description, record.sequence])


def text_to_fasta(input, output, *, source=None):
    """Convert the TSV schema from fasta_to_text to canonical one-line FASTA.

    Reject fields that cannot survive FASTA parsing instead of silently losing
    data. Like readers, conversion is streaming and may write a valid prefix
    before a later error. Output streams remain open; no CLI behavior changes.
    """
    for record in read_sequences(input, format='tsv', source=source,
                                 sequence_column='sequence', id_column='id',
                                 description_column='description'):
        if (any(c.isspace() for c in record.id)
                or record.description != record.description.strip()
                or '\r' in record.description or '\n' in record.description
                or record.sequence.startswith(('>', ';'))):
            raise SequenceReadError('text fields cannot round-trip through FASTA',
                                    record.source, record.provenance.start_line)
        header = record.id + (' ' + record.description if record.description else '')
        output.write(f'>{header}\n{record.sequence}\n')
