"""Protein sequence core (stdlib only).

Layering, lowest first (a module may import only from lower layers):

    L0  residues, provenance   data and utilities, no package imports
    L1  parse                  FASTA reading and symbol classification
    L2  features               sequence-derived profiles           (PR B)
    L3  hints, roi             heuristics                           (PR C/E)
    L4  ss                     predicted secondary structure        (PR D)
        render                 output sink: result dict -> text; imports
                               nothing from the package

Nothing in this package performs I/O or prints at import time.
"""
