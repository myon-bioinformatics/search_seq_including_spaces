"""CLI adapters for roadmap tools. Algorithms and I/O contracts stay importable."""
from contextlib import closing
from dataclasses import asdict
import json
import sys

from ambiguity_candidates import MODES, generate_candidates, select_representatives
from anchor_matcher import find_matches
from result_summary import summarize_csv
from sequence_patterns import load_patterns
from sequence_readers import read_sequences
from streaming_matcher import match_large_file
from result_summary import integer
from sequence_analysis import nucleotide_features, marker_clues, fold_rna
from sequence_preprocess import preprocess, export_records
from protein_ionization import ionization
from structure_tools import parse_structure, geometry
from acquisition import Client


def dump(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(',', ':')) + '\n'


def records(args):
    options = {key: getattr(args, key) for key in ('sequence_column', 'id_column', 'description_column',
                                                   'sequence_path', 'id_path', 'description_path')}
    return read_sequences(sys.stdin if args.file == '-' else args.file,
                          format=args.input_format, **options)


def match_documents(args):
    integer(args.max_matches, "max_matches", 1)
    if args.patterns:
        if args.left is not None or args.right is not None or args.min_gap is not None or args.max_gap is not None or args.iupac or args.strand != '+':
            raise ValueError('--patterns cannot be combined with inline pattern options')
        with open(args.patterns, encoding='utf-8') as handle:
            text = handle.read(1024 * 1024 + 1)
        if len(text) > 1024 * 1024:
            raise ValueError('pattern file exceeds 1 Mi characters')
        patterns = load_patterns(text)
    else:
        if args.left is None or args.right is None:
            raise ValueError('match requires --patterns or both --left and --right')
        patterns = [{'id': 'inline', 'left': args.left, 'right': args.right,
                     'min_gap': args.min_gap or 0, 'max_gap': args.max_gap, 'iupac': args.iupac,
                     'strand': args.strand, 'max_matches': args.max_matches}]
        # Validate even on an empty input, before opening the reader.
        find_matches('', **{k: v for k, v in patterns[0].items() if k != 'id'})
    with closing(records(args)) as inputs:
        for record in inputs:
            for pattern in patterns:
                options = {k: v for k, v in pattern.items() if k != 'id'}
                # CLI cap also applies to packs, never silently overrides a smaller pack limit.
                options['max_matches'] = min(options.get('max_matches', args.max_matches), args.max_matches)
                result = find_matches(record.sequence, **options)
                yield {'schema': 'sequence-matches/1', 'pattern_id': pattern['id'],
                       'record': {'id': record.id, 'description': record.description,
                                  'source': record.source, 'provenance': asdict(record.provenance),
                                  'metadata': record.metadata}, 'result': result}


def _output(documents, jsonl):
    if jsonl:
        def lines():
            with closing(documents):
                for doc in documents:
                    yield dump(doc)
        return lines()
    with closing(documents) as docs:
        first = next(docs, None)
        if first is None or next(docs, None) is not None:
            raise ValueError('--json requires exactly one result; use --jsonl')
        return dump(first)


def run_match(args):
    if args.chunk_size is not None:
        if args.patterns or not args.jsonl or args.input_format not in ("fasta", "fastq", "text"):
            raise ValueError("--chunk-size requires inline anchors, --jsonl and explicit fasta/fastq/text --input-format")
        if any(getattr(args, key) is not None for key in ("sequence_column", "id_column", "description_column", "sequence_path", "id_path", "description_path")):
            raise ValueError("column/path options do not apply to chunk mode")
        if args.left is None or args.right is None:
            raise ValueError("chunk mode requires --left and --right")
        events = match_large_file(sys.stdin if args.file == "-" else args.file, args.left, args.right,
                                  format=args.input_format, chunk_size=args.chunk_size,
                                  min_gap=args.min_gap or 0, max_gap=args.max_gap, iupac=args.iupac,
                                  strand=args.strand, max_matches=args.max_matches)
        def output():
            with closing(events):
                for event in events:
                    yield dump(event)
        return output(), 0
    return _output(match_documents(args), args.jsonl), 0


def run_candidates(args):
    if args.targets and args.mode in ('preserve', 'regex'):
        raise ValueError('representative selection requires concrete candidate generation')
    if args.rank is not None and not args.targets:
        raise ValueError('--rank requires --targets')
    targets = None
    if args.targets:
        with closing(read_sequences(args.targets)) as target_records:
            targets=[]
            for target in target_records:
                if len(targets)>=10000: raise ValueError('target limit exceeded')
                targets.append(target.sequence)
    def documents():
        with closing(records(args)) as inputs:
            for record in inputs:
                report = generate_candidates(record.sequence, mode=args.mode, limit=args.limit, seed=args.seed)
                if targets is not None:
                    report['representatives'] = select_representatives(report['candidates'], targets, rank=8 if args.rank is None else args.rank, seed=args.seed)
                report['record'] = {'id': record.id, 'source': record.source, 'provenance': asdict(record.provenance)}
                yield report
    return _output(documents(), args.jsonl), 0


def run_summary(args):
    return dump(summarize_csv(args.file, gap_merge=args.gap_merge, top=args.top,
                              coordinate_convention=args.coordinates, sequence_length=args.sequence_length)), 0


def reader_arguments(parser):
    parser.add_argument('file', help='sequence file, or - for decoded stdin')
    parser.add_argument('--input-format', choices=('fasta', 'fastq', 'text', 'csv', 'tsv', 'json', 'jsonl', 'genbank'))
    for name in ('sequence-column', 'id-column', 'description-column', 'sequence-path', 'id-path', 'description-path'):
        parser.add_argument('--' + name)


def output_arguments(parser):
    output = parser.add_mutually_exclusive_group()
    output.add_argument('--jsonl', action='store_true', help='stream one JSON document per result')
    output.add_argument('--json', action='store_true', help='require one result (default)')


def add_commands(groups):
    for command, operation in (('nucleotide-features', 'nucleotide'), ('rna-fold', 'rna'), ('preset', 'preset'), ('protein-ionization', 'ionization')):
        analysis = groups.add_parser(command)
        reader_arguments(analysis)
        output_arguments(analysis)
        analysis.set_defaults(handler=run_analysis, operation=operation)
        if operation == 'nucleotide':
            analysis.add_argument('--alphabet', choices=('dna', 'rna'), default='dna')
            analysis.add_argument('--k', type=int, default=2)
            analysis.add_argument('--codon-frame', type=int)
        elif operation == 'ionization':
            analysis.add_argument('--ph', type=float, default=7.0)
        elif operation == 'rna':
            analysis.add_argument('--min-loop', type=int, default=3)
            analysis.add_argument('--max-length', type=int, default=300)
            analysis.add_argument('--no-wobble', action='store_true')
        else:
            analysis.add_argument('--kind', choices=('poly-a', 'tsd', 'str', 'snp'), required=True)
            analysis.add_argument('--motif')
            analysis.add_argument('--min-length', type=int, default=6)
            analysis.add_argument('--position', type=int)
            analysis.add_argument('--allele')
            analysis.add_argument('--repeat-length', type=int, default=4)
            analysis.add_argument('--min-gap', type=int, default=0)
            analysis.add_argument('--max-gap', type=int, default=100)
            analysis.add_argument('--max-results', type=int, default=10000)
    prep = groups.add_parser('preprocess', help='normalize/filter/deduplicate with a required transform receipt')
    reader_arguments(prep)
    prep.add_argument('--alphabet', choices=('dna', 'rna', 'protein'), default='dna')
    prep.add_argument('--normalize-ut', action='store_true')
    prep.add_argument('--deduplicate', choices=('none', 'id', 'sequence'), default='none')
    prep.add_argument('--min-length', type=int, default=1)
    prep.add_argument('--max-length', type=int)
    prep.add_argument('--max-ambiguity', type=float, default=1.0)
    prep.add_argument('--aligned', action='store_true')
    prep.add_argument('--motif')
    prep.add_argument('--export-format', choices=('fasta', 'jsonl', 'mega'), default='fasta')
    prep.add_argument('--receipt', required=True, help='new JSON file for transformation evidence')
    prep.set_defaults(handler=run_preprocess)
    match = groups.add_parser('match', help='search sequence records using anchors or a JSON pattern pack')
    reader_arguments(match)
    output_arguments(match)
    match.add_argument('--chunk-size', type=int, help='bounded single-record mode; explicit fasta/fastq/text and --jsonl')
    match.add_argument('--patterns', help='sequence-patterns/1 JSON file')
    match.add_argument('--left')
    match.add_argument('--right')
    match.add_argument('--min-gap', type=int)
    match.add_argument('--max-gap', type=int)
    match.add_argument('--iupac', action='store_true')
    match.add_argument('--strand', choices=('+', '-', 'both'), default='+')
    match.add_argument('--max-matches', type=int, default=10000)
    match.set_defaults(handler=run_match)
    candidates = groups.add_parser('candidates', help='preserve, regex, exhaustive, sample or pairwise IUPAC candidates')
    reader_arguments(candidates)
    output_arguments(candidates)
    candidates.add_argument('--mode', choices=MODES, default='preserve')
    candidates.add_argument('--limit', type=int, default=10000)
    candidates.add_argument('--seed', type=int, default=0)
    candidates.add_argument('--targets', help='optional sequence file for approximate representative selection')
    candidates.add_argument('--rank', type=int)
    candidates.set_defaults(handler=run_candidates)
    structure = groups.add_parser('structure', help='bounded PDB/mmCIF coordinate geometry')
    structure.add_argument('file')
    structure.add_argument('--input-format', choices=('cif', 'pdb'), default='cif')
    structure.add_argument('--max-atoms', type=int, default=20000)
    structure.add_argument('--cutoff', type=float, default=4.0)
    structure.add_argument('--top', type=int, default=100)
    structure.set_defaults(handler=run_structure)
    fetch = groups.add_parser('fetch', help='explicit optional NCBI/RCSB acquisition (uses network)')
    fetch.add_argument('--operation', choices=('search', 'summary', 'records'), required=True)
    fetch.add_argument('--database', choices=('pubmed','nuccore','protein','structure'), default='pubmed')
    fetch.add_argument('--ids', nargs='+')
    fetch.add_argument('--query')
    fetch.add_argument('--format', choices=('fasta','genbank','cif','pdb'),
                       help='records only: fasta/genbank for sequences (default fasta), '
                            'cif/pdb for one RCSB structure (default cif)')
    fetch.add_argument('--timeout', type=float, default=15)
    fetch.add_argument('--max-bytes', type=int, default=8*1024*1024)
    fetch.add_argument('--limit', type=int, default=100)
    fetch.add_argument('--email')
    fetch.set_defaults(handler=run_fetch)
    summary = groups.add_parser('summarize-hits', help='add a bounded summary of a legacy hit CSV')
    summary.add_argument('file')
    summary.add_argument('--gap-merge', type=int, default=10000)
    summary.add_argument('--top', type=int, default=10)
    summary.add_argument('--sequence-length', type=int)
    summary.add_argument('--coordinates', choices=('0-based-half-open', '1-based-inclusive'), default='0-based-half-open')
    summary.set_defaults(handler=run_summary)


def run_analysis(args):
    def documents():
        with closing(records(args)) as inputs:
            for record in inputs:
                if args.operation == 'nucleotide':
                    report = nucleotide_features(record.sequence, alphabet=args.alphabet, k=args.k, codon_frame=args.codon_frame)
                elif args.operation == 'ionization':
                    report = ionization(record.sequence, ph=args.ph)
                elif args.operation == 'rna':
                    report = fold_rna(record.sequence, min_loop=args.min_loop, max_length=args.max_length, wobble=not args.no_wobble)
                else:
                    report = marker_clues(record.sequence, kind=args.kind, motif=args.motif, min_length=args.min_length,
                                          position=args.position, allele=args.allele, repeat_length=args.repeat_length,
                                          min_gap=args.min_gap, max_gap=args.max_gap, max_results=args.max_results)
                report['record'] = {'id': record.id, 'source': record.source, 'provenance': asdict(record.provenance)}
                yield report
    return _output(documents(), args.jsonl), 0


def run_preprocess(args):
    with closing(records(args)) as inputs:
        output, receipt = preprocess(inputs, alphabet=args.alphabet, normalize_ut=args.normalize_ut,
                                     deduplicate=args.deduplicate, min_length=args.min_length,
                                     max_length=args.max_length, max_ambiguity=args.max_ambiguity,
                                     aligned=args.aligned or args.export_format == 'mega', motif=args.motif)
    body = export_records(output, format=args.export_format, alphabet=args.alphabet)
    # Explicit receipt path uses exclusive creation, protecting existing evidence.
    with open(args.receipt, 'x', encoding='utf-8') as handle:
        handle.write(dump(receipt))
    return body, 0


def run_structure(args):
    with open(args.file, encoding='utf-8') as handle:
        text=handle.read(16*1024*1024+1)
    if len(text)>16*1024*1024: raise ValueError('structure text exceeds 16 Mi characters')
    return dump(geometry(parse_structure(text,format=args.input_format,max_atoms=args.max_atoms),cutoff=args.cutoff,top=args.top)),0


def run_fetch(args):
    if args.format is not None and (args.operation != 'records' or args.database == 'pubmed'):
        raise ValueError('--format is only supported for sequence/structure records')
    client=Client(timeout=args.timeout,max_bytes=args.max_bytes,email=args.email)
    if args.operation=='search':
        if args.ids: raise ValueError('search uses --query, not --ids')
        return dump(client.search(args.query,db=args.database,limit=args.limit)),0
    if not args.ids or args.query: raise ValueError('summary/records requires --ids, without --query')
    if args.operation=='summary': return dump(client.summary(args.ids,db=args.database)),0
    if args.database=='pubmed': return dump(client.pubmed(args.ids)),0
    if args.database=='structure':
        if len(args.ids) != 1:
            raise ValueError('structure records requires exactly one legacy PDB id')
        return dump(client.structure(args.ids[0],format=args.format or 'cif')),0
    return dump({'schema':'ncbi-sequences/1','records':[asdict(r) for r in client.sequences(args.ids,db=args.database,format=args.format or 'fasta')]}),0
