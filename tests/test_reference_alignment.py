"""Pinned external references for annotation search, projection and structural candidates."""
import hashlib
from io import StringIO
import json
from pathlib import Path
import subprocess
import sys
import unittest

from seqlib.comparison import project_interval
from seqlib.sequence_analysis import fold_rna
from seqlib.stockholm import annotated_matches, read_stockholm
from seqlib.structure_tools import geometry, parse_structure

FIXTURES = Path(__file__).parent / 'fixtures'
STOCKHOLM = FIXTURES / 'tRNA5.sto'


def rna_pair_counts():
    with STOCKHOLM.open() as handle:
        rows, tracks = read_stockholm(handle)
    # The pinned reference uses <> for nested pairs and . for unpaired columns.
    stack, pairs = [], []
    for i, c in enumerate(tracks['SS_cons']):
        if c == '<':
            stack.append(i)
        elif c == '>':
            pairs.append((stack.pop(), i))
        elif c != '.':
            raise ValueError(f'unexpected bracket: {c}')
    if stack:
        raise ValueError('incomplete reference structure')
    counts = []
    for seq in rows.values():
        positions = {col: pos for pos, col in enumerate((i for i, c in enumerate(seq) if c != '-'), 1)}
        reference = {(positions[i], positions[j]) for i, j in pairs if i in positions and j in positions}
        predicted = set(map(tuple, fold_rna(seq.replace('-', ''))['pairs']))
        counts.append((len(predicted & reference), len(predicted - reference), len(reference - predicted)))
    return counts


class TestReferenceAlignment(unittest.TestCase):
    def test_annotated_search_cli(self):
        self.assertEqual(hashlib.sha256(STOCKHOLM.read_bytes()).hexdigest(),
                         'c495d7d4c1b6c0851b9d0f51b43899158b2d52969214e26bc7757742e32f15d8')
        output = subprocess.check_output([sys.executable, '-S', str(FIXTURES.parents[1] / 'sequence_tool.py'),
                                          'annotated-match', str(STOCKHOLM), '--left', 'GUUC', '--right', 'GA'])
        report = json.loads(output)
        self.assertEqual(report['records_searched'], 5)
        # Independently transcribed from source alignment columns 54–59.
        expected = [('tRNA1', [53, 58]), ('tRNA2', [52, 57]), ('tRNA4', [52, 57]), ('tRNA5', [53, 58])]
        self.assertEqual([(r['record_id'], r['span']) for r in report['results']], expected)
        for hit in report['results']:
            self.assertEqual((hit['columns'], hit['annotation']), ([54, 59], '<.....'))

    def test_gap_crossing_hit_and_projection(self):
        with STOCKHOLM.open() as handle:
            rows, _ = read_stockholm(handle)
        with STOCKHOLM.open() as handle:
            report = annotated_matches(handle, track='SS_cons', left='GGG', right='AGA')
        hit, = report['results']
        self.assertEqual((hit['record_id'], hit['span'], hit['columns'], hit['annotation'],
                          hit['aligned_sequence']), ('tRNA1', [18, 23], [18, 24], '.....>>', 'GGG-AGA'))
        result = project_interval(rows['tRNA1'], rows['tRNA2'], alphabet='rna',
                                  start=hit['span'][0], end=hit['span'][1])
        self.assertEqual(result['blocks'], [
            {'a': [18, 20], 'b': [17, 19], 'columns': [18, 20]},
            {'a': [21, 23], 'b': [21, 23], 'columns': [22, 24]}])

    def test_malformed_alignment_and_missing_track(self):
        for data in ('a AA\n', '# STOCKHOLM 1.0\na A\n', '# STOCKHOLM 1.0\na AC\n#=GC X .\n//\n',
                     '# STOCKHOLM 1.0\na A\n\nb G\n//\n', '# STOCKHOLM 1.0\na A\na U\n//\n',
                     '# STOCKHOLM 1.0\na A\n//\n# STOCKHOLM 1.0\n//\n'):
            with self.subTest(data=data), self.assertRaises(ValueError):
                read_stockholm(StringIO(data))
        with STOCKHOLM.open() as handle, self.assertRaisesRegex(ValueError, 'unknown column'):
            annotated_matches(handle, track='absent', left='G', right='U')

    def test_fold_comparison_to_consensus(self):
        # Each tuple is TP, FP, FN for paired ungapped positions on one row.
        # A changed predictor should trigger re-evaluation of the documented baseline.
        self.assertEqual(rna_pair_counts(), [(10, 18, 11), (7, 20, 14), (4, 23, 17),
                                            (10, 18, 11), (5, 22, 16)])

    def test_crambin_disulfide_candidates_against_deposited_pairs(self):
        text = (FIXTURES / '1crn-sulfur.pdb').read_text()
        expected = {tuple(sorted((int(line[17:21]), int(line[31:35]))))
                    for line in text.splitlines() if line.startswith('SSBOND')}
        self.assertEqual(expected, {(3, 40), (4, 32), (16, 26)})
        report = geometry(parse_structure(text, format='pdb'), cutoff=3)
        refs = report['atom_references']
        observed = {tuple(sorted(int(refs[str(i)]['residue_id']) for i in contact['atoms']))
                    for contact in report['contacts'] if 'disulfide-candidate' in contact['labels']}
        self.assertEqual(observed, expected)

    def test_coordinate_graph_and_model_selection(self):
        from seqlib.structure_tools import structure_graph
        for fixture, expected in [('1crn-sulfur.pdb', 2.003), ('1ehz-g19-c56.pdb', 2.817)]:
            raw = (FIXTURES / fixture).read_text()
            parsed = parse_structure(raw, format='pdb')
            graph = structure_graph(parsed)
            atoms = [n for n in graph['nodes'] if n['type'] == 'Atom']
            self.assertEqual([(n['x'], n['y'], n['z']) for n in atoms],
                             [(a.x, a.y, a.z) for a in parsed['atoms']])
            contacts = [e for e in graph['edges'] if e['relation'] == 'nearAtom']
            self.assertTrue(all(e['symmetric'] for e in contacts))
            self.assertAlmostEqual(contacts[0]['distance_angstrom'], expected, delta=.01)
            small = structure_graph(parsed, top=1)
            self.assertEqual(small['nodes'], graph['nodes'])
            self.assertEqual(small['omitted_contacts'], len(contacts) - 1)
        two = parse_structure('MODEL        1\n' + raw + 'ENDMDL\nMODEL        2\n' + raw, format='pdb')
        with self.assertRaisesRegex(ValueError, 'select --model'):
            structure_graph(two)
        selected = structure_graph(two, model='2')
        self.assertTrue(all(n['model'] == '2' for n in selected['nodes'] if n['type'] == 'Atom'))
