"""Optional pKa-based sequence proxy; leaves existing protein golden outputs intact."""
from collections import Counter
import math
from anchor_matcher import normalize_sequence

# Bjellqvist model constants; provenance and limits in docs/ROADMAP_TOOLS.md.
POSITIVE = {'Nterm':7.5,'K':10.0,'R':12.0,'H':5.98}
NEGATIVE = {'Cterm':3.55,'D':4.05,'E':4.45,'C':9.0,'Y':10.0}
N_TERMINAL = {'A':7.59,'M':7.0,'S':6.93,'P':8.36,'T':6.82,'V':7.44,'E':7.7}
C_TERMINAL = {'D':4.55,'E':4.75}


def ionization(text, *, ph=7.0):
    if isinstance(ph,bool) or not isinstance(ph,(int,float)) or not math.isfinite(ph) or not 0<=ph<=14:
        raise ValueError('pH must be finite in [0,14]')
    sequence=normalize_sequence(text).sequence
    if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWY'):
        raise ValueError('ionization requires a nonempty canonical protein sequence')
    counts=Counter(sequence);counts['Nterm']=counts['Cterm']=1
    positive={**POSITIVE,'Nterm':N_TERMINAL.get(sequence[0],POSITIVE['Nterm'])}
    negative={**NEGATIVE,'Cterm':C_TERMINAL.get(sequence[-1],NEGATIVE['Cterm'])}
    def charge(value):
        return (sum(counts[name]/(1+10**(value-pka)) for name,pka in positive.items())-
                sum(counts[name]/(1+10**(pka-value)) for name,pka in negative.items()))
    lo,hi=0.0,14.0
    bracketed=charge(lo)>=0>=charge(hi)
    for _ in range(60):
        mid=(lo+hi)/2
        if charge(mid)>0: lo=mid
        else: hi=mid
    return {'schema':'protein-ionization/1','ph':ph,'net_charge':charge(ph),
            'estimated_pi':(lo+hi)/2 if bracketed else None,
            'model':'Bjellqvist pKa sequence proxy; free termini',
            'pka':{'positive':positive,'negative':negative},
            'limitations':'no solvent, folding, modifications or residue-interaction correction'}
