"""S3 bounded, count-preserving summaries; never rewrites the source results."""
from collections import Counter
import csv
import hashlib
import math


def integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return value


def _bounded(regions, top, count_key):
    selected = sorted(regions, key=lambda r: (-r[count_key], r['start'], r['end'], r['source_rows'][0]))[:top]
    return {'total_regions': len(regions), 'regions': selected,
            'omitted_regions': len(regions) - len(selected),
            'omitted_items': sum(r[count_key] for r in regions) - sum(r[count_key] for r in selected)}


def summarize_hits(rows, *, gap_merge=10000, top=10,
                   coordinate_convention='0-based-half-open', sequence_length=None):
    """Summarize legacy CSV rows or start/end/back_elements mappings.

    gap_merge is distance between consecutive STARTS (single-linkage). Zero
    disables merging, even for duplicate hits. Sort by start, then source row.
    source_rows is an enclosing 1-based DATA-row range; rows need not be
    contiguous after coordinate sorting. No genomic coordinate repair is done.
    """
    integer(gap_merge, 'gap_merge')
    integer(top, 'top', 1)
    if coordinate_convention not in ('0-based-half-open', '1-based-inclusive'):
        raise ValueError('unsupported coordinate convention')
    if sequence_length is not None:
        integer(sequence_length, 'sequence_length', 1)
    hits, frequencies = [], Counter()
    minimum = 0 if coordinate_convention == '0-based-half-open' else 1
    for index, row in enumerate(rows, 1):
        try:
            start, end = row['start_num'], row['end_num']
            # CSV values are integers written as decimal strings, not floats.
            start = int(start) if isinstance(start, str) else start
            end = int(end) if isinstance(end, str) else end
            integer(start, 'start_num', minimum)
            integer(end, 'end_num', start + (minimum == 0))
            label = row['back_elements']
            if not isinstance(label, str):
                raise ValueError('back_elements must be text')
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f'invalid hit data row {index}: {error}') from None
        if sequence_length is not None and end > sequence_length:
            raise ValueError(f'hit data row {index} exceeds sequence_length')
        hits.append((start, end, index))
        frequencies[label] += 1
    hits.sort(key=lambda h: (h[0], h[2]))
    regions, previous_start = [], None
    for start, end, index in hits:
        if not regions or not gap_merge or start - previous_start > gap_merge:
            regions.append({'start': start, 'end': end, 'hit_count': 1,
                            'source_rows': [index, index]})
        else:
            current = regions[-1]
            current['end'] = max(current['end'], end)
            current['hit_count'] += 1
            current['source_rows'][0] = min(current['source_rows'][0], index)
            current['source_rows'][1] = max(current['source_rows'][1], index)
        previous_start = start
    covered_span = (max(h[1] for h in hits) - hits[0][0] + minimum) if hits else 0
    denominator = sequence_length if sequence_length is not None else covered_span
    common = sorted(frequencies.items(), key=lambda p: (-p[1], p[0]))[:top]
    return {'schema': 'sequence-hit-summary/1', 'coordinate_convention': coordinate_convention,
            'merge_rule': 'consecutive-start-distance; zero disables merging',
            'gap_merge': gap_merge, 'top': top, 'total_hits': len(hits),
            **_bounded(regions, top, 'hit_count'),
            'back_elements': [{'value': k, 'count': v} for k, v in common],
            'omitted_back_elements': len(frequencies) - len(common),
            'omitted_back_element_hits': len(hits) - sum(v for _, v in common),
            'density_per_mb': len(hits) * 1000000 / denominator if denominator else None,
            'density_denominator_bp': denominator,
            'density_scope': 'sequence' if sequence_length is not None else 'observed-bounding-span',
            'source_row_convention': '1-based data rows; enclosing range, not exact membership'}


def summarize_csv(path, **options):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, strict=True)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not {'start_num', 'end_num', 'back_elements'} <= set(fields):
            raise ValueError('CSV requires unique start_num/end_num/back_elements columns')
        report = summarize_hits(reader, **options)
    # A separate bounded streaming pass records identity of the untouched source.
    with open(path, 'rb') as handle:
        digest = hashlib.sha256()
        for block in iter(lambda: handle.read(65536), b''):
            digest.update(block)
    report['source'] = {'path': str(path), 'sha256': digest.hexdigest()}
    return report


def summarize_windows(report, *, threshold=90.0, min_coverage=0.0, top=10):
    """Coalesce consecutive S2 windows of the same classification.

    Overlapping windows in different classes may produce overlapping regions;
    these are window evidence, not a partition or inference of homology.
    Null identity or insufficient coverage forms an explicit unknown region.
    """
    integer(top, 'top', 1)
    for name, value in (('threshold', threshold), ('min_coverage', min_coverage)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 100:
            raise ValueError(f'{name} must be a finite percentage in [0, 100]')
    if report.get('schema') != 'sequence-comparison-windows/v1':
        raise ValueError('expected S2 window report')
    regions, previous_start, counts = [], None, Counter()
    step = integer(report['step'], 'step', 1)
    for index, window in enumerate(report['windows'], 1):
        identity, coverage = window['percent_identity'], window['coverage_percent']
        kind = ('unknown' if identity is None or coverage is None or coverage < min_coverage
                else 'conserved' if identity >= threshold else 'divergent')
        counts[kind] += 1
        start, end = window['start'], window['end']
        if (regions and kind == regions[-1]['classification'] and start == previous_start + step
                and start <= regions[-1]['end'] + 1):
            regions[-1]['end'] = max(regions[-1]['end'], end)
            regions[-1]['window_count'] += 1
            regions[-1]['source_rows'][1] = index
        else:
            regions.append({'start': start, 'end': end, 'classification': kind,
                            'window_count': 1, 'source_rows': [index, index]})
        previous_start = start
    return {'schema': 'sequence-window-summary/1', 'coordinate_convention': report['coordinate_convention'],
            'threshold': threshold, 'min_coverage': min_coverage, 'top': top,
            'total_windows': sum(counts.values()), 'classification_counts': dict(sorted(counts.items())),
            **_bounded(regions, top, 'window_count'), 'provenance': dict(report['provenance']),
            'source_row_convention': '1-based S2 windows; exact contiguous range',
            'interpretation': 'window identity evidence; regions may overlap; not homology'}
