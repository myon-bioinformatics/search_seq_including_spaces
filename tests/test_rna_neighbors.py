import json
from io import StringIO
from pathlib import Path
import subprocess
import sys
import unittest

from seqlib.rna_neighbors import neighbors
from seqlib.structure_tools import parse_structure

FIX = Path(__file__).parent / 'fixtures'


class TestNeighbors(unittest.TestCase):
    def test_focus_partners_and_column_offsets(self):
        with (FIX / 'tRNA5.sto').open() as f:
            graph = neighbors(f, 'tRNA1', 31)
        self.assertEqual([n['name'] for n in graph['nodes']], ['G30', 'A31', 'C32', 'U39'])
        self.assertEqual((graph['pair_status'], graph['partner_column']), ('consensus-pair', 40))
        self.assertIsNone(graph['coordinates'])

    def test_crossing_pairs_gap_and_absent_track(self):
        def inspect(seq, track, pos):
            annotation = f'#=GC SS_cons {track}\n' if track else ''
            return neighbors(StringIO(f'# STOCKHOLM 1.0\nr {seq}\n{annotation}//\n'), 'r', pos)
        self.assertEqual(inspect('ACGU', '<A>a', 4)['partner_column'], 2)
        self.assertEqual(inspect('AC-U', '<A>a', 1)['pair_status'], 'partner-gap')
        self.assertEqual(inspect('A', None, 1)['pair_status'], 'track-unavailable')
        for track in ('<<<', '.>.', '..?'):
            with self.subTest(track=track), self.assertRaises(ValueError):
                inspect('ACG', track, 1)

    def test_real_contact_via_cli(self):
        command = [sys.executable, '-S', str(FIX.parents[1] / 'sequence_tool.py'), 'base-neighbors',
                   str(FIX / 'tRNA5.sto'), '--record', 'tRNA1', '--position', '19',
                   '--structure', str(FIX / '1ehz-g19-c56.pdb'), '--mapping', str(FIX / '1ehz-g19-c56-map.json')]
        result = json.loads(subprocess.check_output(command))
        edge, = [e for e in result['edges'] if e['kind'] == 'spatial-contact']
        self.assertEqual((edge['target'], edge['atom_names'], edge['distance_angstrom']),
                         ('56', ['N2', 'O2'], 2.8174))
        result = json.loads(subprocess.check_output(command + ['--cutoff', '2']))
        self.assertEqual(result['coordinates']['mapped_positions'], 2)
        self.assertFalse(any(e['kind'] == 'spatial-contact' for e in result['edges']))

    def test_invalid_coordinate_mappings(self):
        structure = parse_structure((FIX / '1ehz-g19-c56.pdb').read_text(), format='pdb')
        for change in ('wrong-chain', 'wrong-base', 'duplicate', 'unmapped-focus'):
            mapping = json.loads((FIX / '1ehz-g19-c56-map.json').read_text())
            if change == 'wrong-chain':
                mapping['chain'] = 'B'
            elif change == 'wrong-base':
                mapping['residues'][0]['component'] = 'C'
            elif change == 'duplicate':
                mapping['residues'].append(mapping['residues'][0])
            else:
                mapping['residues'].pop(0)
            with self.subTest(change=change), (FIX / 'tRNA5.sto').open() as f, self.assertRaises(ValueError):
                neighbors(f, 'tRNA1', 19, structure=structure, mapping=mapping)
