"""Bounded PDB/mmCIF atom-site subsets and geometry-based contact candidates."""
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
import hashlib
import itertools
import math
import re

from result_summary import integer


@dataclass(frozen=True)
class Atom:
    id: str
    name: str
    element: str
    residue: str
    chain: str
    residue_id: str
    model: str
    altloc: str
    x: float
    y: float
    z: float
    hetero: bool = False

    @property
    def residue_key(self):
        return self.model, self.chain, self.residue_id, self.residue


def _tokens(text):
    """CIF1 data tokens: quotes terminate at whitespace, semicolon blocks at col 1."""
    lines=iter(text.splitlines())
    for line in lines:
        if line.startswith(';'):
            block=[line[1:]]
            for continuation in lines:
                if continuation.startswith(';'):
                    if continuation[1:].strip(): raise ValueError('text after CIF multiline terminator')
                    break
                block.append(continuation)
            else:
                raise ValueError('unterminated CIF multiline value')
            yield '\n'.join(block), True
            continue
        index=0
        while index<len(line):
            if line[index].isspace(): index+=1;continue
            if line[index]=='#': break
            start=index
            if line[index] in ('"', "'"):
                quote=line[index];index+=1;start=index
                while index<len(line) and not (line[index]==quote and (index+1==len(line) or line[index+1].isspace())):
                    index+=1
                if index==len(line): raise ValueError('unterminated CIF quoted value')
                yield line[start:index],True
                index+=1
            else:
                while index<len(line) and not line[index].isspace(): index+=1
                yield line[start:index],False


def _cif_rows(text):
    tokens=iter(_tokens(text));pending=None;blocks=0
    while True:
        token=pending if pending is not None else next(tokens,None)
        pending=None
        if token is None: return
        value,quoted=token
        if not quoted and value.lower().startswith('data_'):
            blocks+=1
            if blocks>1: raise ValueError('mmCIF subset supports a single data block')
        elif not quoted and value=='loop_':
            names=[]
            token=next(tokens,None)
            while token is not None and not token[1] and token[0].startswith('_'):
                names.append(token[0]);token=next(tokens,None)
            if not names or len(set(names))!=len(names): raise ValueError('invalid CIF loop fields')
            atom_loop=any(name.startswith('_atom_site.') for name in names)
            if atom_loop and not all(name.startswith('_atom_site.') for name in names): raise ValueError('mixed atom-site loop')
            row=[]
            while token is not None:
                value,quoted=token
                boundary=not quoted and (value.startswith('_') or value in ('loop_','stop_','global_') or value.startswith(('data_','save_')))
                if boundary:
                    if row: raise ValueError('incomplete CIF loop row')
                    pending=token;break
                row.append(value)
                if len(row)==len(names):
                    if atom_loop: yield dict(zip(names,row))
                    row=[]
                token=next(tokens,None)
            if row: raise ValueError('incomplete CIF loop row')
        elif not quoted and value.startswith('_'):
            token=next(tokens,None)
            if token is None or (not token[1] and token[0].startswith('_')):
                raise ValueError('missing CIF scalar value')
        elif not quoted and value.startswith(('save_','global_')):
            raise ValueError('CIF save/global frames unsupported')


def parse_structure(text, *, format='cif', max_atoms=20000):
    integer(max_atoms,'max_atoms',1)
    if format not in ('cif','pdb'): raise ValueError('format must be cif or pdb')
    atoms=[]
    def add(atom):
        if len(atoms)>=max_atoms: raise ValueError('atom limit exceeded')
        if not all(math.isfinite(v) for v in (atom.x,atom.y,atom.z)):
            raise ValueError('nonfinite atom coordinates')
        atoms.append(atom)
    if format=='pdb':
        model='1'
        for number,line in enumerate(text.splitlines(),1):
            if line.startswith('MODEL '):
                model=line[10:14].strip()
                if not model: raise ValueError('empty MODEL number')
            if line[:6].strip() not in ('ATOM','HETATM'): continue
            try:
                if len(line)<54: raise ValueError('short atom line')
                element=line[76:78].strip().upper()
                if not element: raise ValueError('PDB subset requires explicit element column')
                add(Atom(line[6:11].strip(),line[12:16].strip(),element,line[17:20].strip(),
                         line[21:22].strip(),line[22:27].strip(),model,line[16:17].strip(),
                         float(line[30:38]),float(line[38:46]),float(line[46:54]),line.startswith('HETATM')))
            except ValueError as error: raise ValueError(f'PDB line {number}: {error}') from None
    else:
        for row in _cif_rows(text):
            def get(*names,default=None):
                for name in names:
                    value=row.get('_atom_site.'+name)
                    if value is not None and value not in ('.','?'): return value
                if default is not None: return default
                raise ValueError('missing atom_site '+names[0])
            add(Atom(get('id'),get('auth_atom_id','label_atom_id'),get('type_symbol').upper(),
                     get('auth_comp_id','label_comp_id'),get('auth_asym_id','label_asym_id'),
                     get('auth_seq_id','label_seq_id')+get('pdbx_PDB_ins_code',default=''),
                     get('pdbx_PDB_model_num',default='1'),get('label_alt_id',default=''),
                     float(get('Cartn_x')),float(get('Cartn_y')),float(get('Cartn_z')),get('group_PDB',default='ATOM')=='HETATM'))
    if not atoms: raise ValueError('structure has no supported atomic coordinates')
    return {'schema':'structure-atoms/1','format':format,'source_sha256':hashlib.sha256(text.encode()).hexdigest(),
            'atoms':atoms,'subset':'coordinate atom records only; no full dictionary/biological assembly reconstruction'}


def geometry(structure, *, cutoff=4.0, top=100, max_contacts=100000):
    if isinstance(cutoff,bool) or not isinstance(cutoff,(int,float)) or not math.isfinite(cutoff) or not 0<cutoff<=20:
        raise ValueError('cutoff must be finite in (0,20] angstrom')
    integer(top,'top',1);integer(max_contacts,'max_contacts',1)
    atoms=structure['atoms'];cells=defaultdict(list)
    contacts=[];types=Counter();interfaces=Counter();residue_pairs=set();nearest={}
    by_residue=defaultdict(list)
    for atom in atoms:
        by_residue[atom.residue_key].append(atom)
    def distance(a,b): return math.dist((a.x,a.y,a.z),(b.x,b.y,b.z))
    def charged(a):
        if (a.residue=='LYS' and a.name=='NZ') or (a.residue=='ARG' and a.name in ('NE','NH1','NH2')): return 1
        if (a.residue=='ASP' and a.name in ('OD1','OD2')) or (a.residue=='GLU' and a.name in ('OE1','OE2')): return -1
        return 0
    def hydrogen_angle(donor,acceptor):
        # Deliberately narrow typing: explicit H on N/O, oxygen acceptors only.
        if donor.element not in ('N','O') or acceptor.element!='O': return None
        angles=[]
        for hydrogen in by_residue[donor.residue_key]:
            if hydrogen.element!='H' or (hydrogen.altloc and donor.altloc and hydrogen.altloc!=donor.altloc): continue
            if hydrogen.altloc and acceptor.altloc and hydrogen.altloc!=acceptor.altloc: continue
            dh=distance(donor,hydrogen);ha=distance(hydrogen,acceptor)
            if not 0.5<=dh<=1.3 or not 0<ha<=2.5: continue
            u=(donor.x-hydrogen.x,donor.y-hydrogen.y,donor.z-hydrogen.z)
            v=(acceptor.x-hydrogen.x,acceptor.y-hydrogen.y,acceptor.z-hydrogen.z)
            angle=math.degrees(math.acos(max(-1,min(1,sum(x*y for x,y in zip(u,v))/(dh*ha)))))
            if angle>=120: angles.append(angle)
        return max(angles) if angles else None
    for i,a in enumerate(atoms):
        cell=tuple(math.floor(v/cutoff) for v in (a.x,a.y,a.z))
        for delta in itertools.product((-1,0,1),repeat=3):
            neighbor=tuple(v+d for v,d in zip(cell,delta))
            for j in cells[(a.model,*neighbor)]:
                b=atoms[j]
                if a.residue_key==b.residue_key or (a.altloc and b.altloc and a.altloc!=b.altloc): continue
                d=distance(a,b)
                if d>cutoff: continue
                if len(contacts)>=max_contacts: raise ValueError('contact limit exceeded; narrow cutoff/input')
                labels=['geometry-contact']
                angles=[angle for angle in (hydrogen_angle(a,b),hydrogen_angle(b,a)) if angle is not None] if d<=3.5 else []
                if angles: labels.append('hydrogen-bond-candidate')
                if a.name==b.name=='SG' and a.residue==b.residue=='CYS' and 1.8<=d<=2.3: labels.append('disulfide-candidate')
                if charged(a)*charged(b)==-1 and d<=4.0: labels.append('salt-bridge-candidate')
                if a.residue in ('HOH','WAT') or b.residue in ('HOH','WAT'): labels.append('water-contact')
                if a.element in ('ZN','FE','MG','MN','CA','CU','CO','NI') or b.element in ('ZN','FE','MG','MN','CA','CU','CO','NI'): labels.append('metal-contact')
                if a.hetero or b.hetero: labels.append('hetero-contact')
                if a.chain!=b.chain:
                    labels.append('interface-contact');interfaces[tuple(sorted((a.chain,b.chain)))]+=1
                contacts.append({'atoms':[j,i],'distance_angstrom':round(d,6),'labels':labels, 'hydrogen_angle_degrees':round(max(angles),6) if angles else None})
                types.update(labels)
                residue_pairs.add(tuple(sorted((a.residue_key,b.residue_key))))
                for source,target in ((i,j),(j,i)):
                    if source not in nearest or (d,target)<nearest[source]: nearest[source]=(d,target)
        cells[(a.model,*cell)].append(i)
    # Centroid/extent must not mix independent models.
    models={}
    for model in sorted({a.model for a in atoms}):
        selected=[a for a in atoms if a.model==model]
        points=[(a.x,a.y,a.z) for a in selected]
        centroid=[sum(p[k] for p in points)/len(points) for k in range(3)]
        models[model]={'atoms':len(points),'centroid':centroid,
                       'bounds':[[min(p[k] for p in points),max(p[k] for p in points)] for k in range(3)],
                       'extent':max(math.dist(p,centroid) for p in points)}
    contacts.sort(key=lambda c:(c['distance_angstrom'],c['atoms']))
    selected=contacts[:top]
    referenced=sorted({i for contact in selected for i in contact['atoms']})
    return {'schema':'structure-geometry/1','source_sha256':structure['source_sha256'],
            'atom_count':len(atoms),'residue_count':len({a.residue_key for a in atoms}),
            'chain_count':len({(a.model,a.chain) for a in atoms}),'models':models,
            'cutoff_angstrom':cutoff,'contact_count':len(contacts),'contacts':selected,
            'atom_references':{str(i):asdict(atoms[i]) for i in referenced},
            'omitted_contacts':max(0,len(contacts)-top),'contact_types':dict(types),
            'residue_contact_pairs':len(residue_pairs),
            'interfaces':[{'chains':list(pair),'contacts':count} for pair,count in sorted(interfaces.items())],
            'alternate_location_atoms':sum(bool(a.altloc) for a in atoms),
            'missing_coordinates':'not inferred from absent records',
            'interpretation':'geometry-based candidates only; hydrogen subset uses explicit H on N/O and oxygen acceptors; no proven bonds'}
