from dataclasses import replace
import math
import unittest
from structure_tools import Atom, parse_structure, geometry

CIF='''data_demo
_title
;multiline metadata
not atom coordinates
;
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.label_alt_id
_atom_site.pdbx_PDB_model_num
ATOM 1 S SG CYS A 1 0 0 0 . 1
ATOM 2 S SG CYS B 2 2.05 0 0 . 1
ATOM 3 S SG CYS B 2 0 0 0 . 2
#
'''


class StructureTests(unittest.TestCase):
    def test_cif_models_distance_candidates_and_evidence(self):
        structure=parse_structure(CIF)
        doc=geometry(structure,cutoff=3)
        self.assertEqual(doc['atom_count'],3)
        self.assertEqual(doc['contact_count'],1)
        self.assertEqual(doc['contacts'][0]['distance_angstrom'],2.05)
        self.assertIn('disulfide-candidate',doc['contacts'][0]['labels'])
        self.assertEqual(len(doc['models']),2)
        self.assertEqual(doc['residue_contact_pairs'],1)

    def test_pdb_fixed_columns_and_invalid_coordinates(self):
        line=f'ATOM  {1:5d} {"SG":>4s} {"CYS":3s} A{1:4d}    {0:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{10:6.2f}          {"S":>2s}'
        structure=parse_structure(line,format='pdb')
        self.assertEqual(structure['atoms'][0].name,'SG')
        self.assertEqual(structure['atoms'][0].element,'S')
        with self.assertRaises(ValueError): parse_structure('ATOM  1',format='pdb')
        with self.assertRaises(ValueError): parse_structure(CIF.replace('2.05','nan'))
        with self.assertRaises(ValueError): parse_structure(CIF,max_atoms=1)
        with self.assertRaises(ValueError): parse_structure(CIF.replace(' . 2',' .'))

    def test_geometry_matches_bruteforce_and_bounds(self):
        atoms=[Atom(str(i),'CA','C','ALA',str(i%2),str(i),'1','',i*1.1,(i%3)*1.2,0) for i in range(15)]
        structure={'atoms':atoms,'source_sha256':'synthetic'}
        expected=sum(math.dist((a.x,a.y,a.z),(b.x,b.y,b.z))<=3 for i,a in enumerate(atoms) for b in atoms[:i])
        doc=geometry(structure,cutoff=3,top=2)
        self.assertEqual(doc['contact_count'],expected)
        self.assertEqual(doc['omitted_contacts'],expected-2)
        with self.assertRaises(ValueError): geometry(structure,max_contacts=1)
        atoms[0]=replace(atoms[0],altloc='A');atoms[1]=replace(atoms[1],altloc='B')
        filtered=geometry(structure,cutoff=3)
        self.assertEqual(filtered['contact_count'],expected-1)

    def test_explicit_hydrogen_angle_candidate(self):
        atoms=[Atom('1','N','N','ALA','A','1','1','',0,0,0),
               Atom('2','H','H','ALA','A','1','1','',1,0,0),
               Atom('3','O','O','ALA','A','2','1','',2.8,0,0)]
        report=geometry({'atoms':atoms,'source_sha256':'fixture'})
        heavy=next(c for c in report['contacts'] if c['atoms']==[0,2])
        self.assertEqual(heavy['hydrogen_angle_degrees'],180)
        self.assertIn('hydrogen-bond-candidate',heavy['labels'])
        self.assertEqual(report['atom_references']['0']['id'],'1')
