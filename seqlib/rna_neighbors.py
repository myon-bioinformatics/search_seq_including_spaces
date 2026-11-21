"""Local RNA neighbourhoods in a supplied alignment and explicitly mapped coordinates."""
from collections import defaultdict
import math

from seqlib.result_summary import integer
from seqlib.stockholm import read_stockholm


def _partners(track):
    """Pair WUSS brackets by family, including upper/lowercase pseudoknot labels."""
    closes = dict(zip('>)]}', '<([{'))
    stacks = defaultdict(list)
    pairs = {}
    for col, symbol in enumerate(track, 1):
        if symbol in '<([{' or 'A' <= symbol <= 'Z':
            stacks[symbol].append(col)
        elif symbol in closes or 'a' <= symbol <= 'z':
            stack = stacks[closes.get(symbol, symbol.upper())]
            if not stack:
                raise ValueError(f'unmatched {symbol!r} at column {col}')
            left = stack.pop()
            pairs[left], pairs[col] = col, left
        elif symbol not in '._-,:~':
            raise ValueError(f'unsupported structure symbol {symbol!r} at column {col}')
    if any(stacks.values()):
        raise ValueError('unclosed structure annotation')
    return pairs


def _contacts(sequence, position, structure, mapping, record, cutoff):
    if not math.isfinite(cutoff) or not 0 < cutoff <= 20:
        raise ValueError('cutoff must be in (0,20] angstrom')
    if not isinstance(mapping, dict) or mapping.get('record') != record:
        raise ValueError('mapping must name the selected Stockholm record')
    model, chain = mapping.get('model'), mapping.get('chain')
    if not isinstance(model, str) or not isinstance(chain, str) or not isinstance(mapping.get('residues'), list):
        raise ValueError('mapping requires string model and chain and a residues list')
    # Restrict geometry to named base N/O atoms. Sugar/phosphate atoms are excluded.
    names = {'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'O2', 'O4', 'O6'}
    atoms = defaultdict(list)
    for atom in structure['atoms']:
        if (atom.model == model and atom.chain == chain and atom.name in names
                and atom.element == atom.name[0]):
            atoms[atom.residue_id, atom.residue].append(atom)
    mapped, ids = {}, set()
    for row in mapping['residues']:
        if not isinstance(row, dict):
            raise ValueError('each mapped residue must be an object')
        pos = integer(row.get('position'), 'mapped position', 1)
        rid, component = row.get('residue_id'), row.get('component')
        if not isinstance(rid, str) or not isinstance(component, str):
            raise ValueError('residue_id and component must be strings')
        if pos > len(sequence) or pos in mapped or rid in ids:
            raise ValueError('mapping must be one-to-one and in sequence bounds')
        # A modified residue needs an explicit parent base from the user's source.
        base = component if component in ('A', 'C', 'G', 'U') else row.get('parent_base')
        if base != sequence[pos - 1] or base not in ('A', 'C', 'G', 'U'):
            raise ValueError(f'base mismatch or missing parent_base at position {pos}')
        if (rid, component) not in atoms:
            raise ValueError(f'no supported base atoms at {model}/{chain}/{rid}/{component}')
        mapped[pos] = (rid, component)
        ids.add(rid)
    if position not in mapped:
        raise ValueError('focus position missing from the supplied mapping')
    hits = []
    for pos, key in sorted(mapped.items()):
        if pos == position:
            continue
        nearest = min(((math.dist((a.x, a.y, a.z), (b.x, b.y, b.z)), a.name, b.name, a.id, b.id)
                       for a in atoms[mapped[position]] for b in atoms[key]
                       if not a.altloc or not b.altloc or a.altloc == b.altloc), default=None)
        if nearest and nearest[0] <= cutoff:
            distance, atom_a, atom_b, id_a, id_b = nearest
            hits.append((pos, {'distance_angstrom': round(distance, 4),
                               'atom_names': [atom_a, atom_b], 'atom_ids': [id_a, id_b],
                               'residue_ids': [mapped[position][0], key[0]],
                               'components': [mapped[position][1], key[1]]}))
    return hits, {'model': model, 'chain': chain, 'mapped_positions': len(mapped),
                  'focus_residue_id': mapped[position][0], 'focus_component': mapped[position][1],
                  'cutoff_angstrom': cutoff, 'source_sha256': structure['source_sha256']}


def neighbors(handle, record, position, *, structure=None, mapping=None, cutoff=4.0):
    """Return ONTO-compatible nodes and edges incident to one ungapped RNA position."""
    rows, tracks = read_stockholm(handle)
    if record not in rows:
        raise ValueError(f'unknown Stockholm record {record!r}')
    aligned = rows[record].upper()
    sequence = aligned.replace('-', '')
    if any(c not in 'ACGURYSWKMBDHVN' for c in sequence):
        raise ValueError('expected an RNA alignment row')
    columns = [col for col, base in enumerate(aligned, 1) if base != '-']
    if integer(position, 'position', 1) > len(columns):
        raise ValueError('position exceeds ungapped record length')
    included, edges = {position}, []

    def link(a, b, kind, basis, **extra):
        included.update((a, b))
        edges.append({'source': str(a), 'target': str(b), 'kind': kind, 'basis': basis, **extra})

    if position > 1:
        link(position - 1, position, 'next', 'sequence')
    if position < len(sequence):
        link(position, position + 1, 'next', 'sequence')
    column = columns[position - 1]
    track = tracks.get('SS_cons')
    partner = _partners(track).get(column) if track is not None else None
    if partner is None:
        status = 'track-unavailable' if track is None else 'no-annotated-pair'
    elif aligned[partner - 1] == '-':
        status = 'partner-gap'
    else:
        other = columns.index(partner) + 1
        link(*sorted((position, other)), 'paired', 'SS_cons')
        status = 'consensus-pair'
    coordinates = None
    if (structure is None) != (mapping is None):
        raise ValueError('structure and mapping are required together')
    if structure is not None:
        contacts, coordinates = _contacts(sequence, position, structure, mapping, record, cutoff)
        for other, evidence in contacts:
            link(position, other, 'spatial-contact', 'coordinate-distance', **evidence)
    nodes = [{'id': str(p), 'name': f'{sequence[p-1]}{p}', 'position': p,
              'column': columns[p-1], 'base': sequence[p-1]} for p in sorted(included)]
    return {'schema': 'rna-neighbors/1', 'record': record, 'focus': str(position), 'nodes': nodes,
            'edges': edges, 'pair_status': status, 'partner_column': partner, 'coordinates': coordinates}
