#!/usr/bin/env python3

import os
import re
import sys
import numpy as np
import argparse as arg

def skiplines(openfile, nlines=0):
    '''Skips nlines + 1 lines in openfile. In other words, if nlines=0 it will
    go to the next line.'''

    for i in range(nlines):
        next(openfile)

    return next(openfile)


def extract(filename):

    with open(filename) as f:

        atoms = []
        structure = []

        next(f)
        for line in f:
            data = line.strip().split()
            atom = data[1][0]
            coords = list(map(float, data[2:5]))
            atoms.append(atom)
            structure.append(coords)

    structure = np.array(structure)

    return atoms, structure


def key_to_idxs(keyfile):

    idxs = []
    with open(keyfile) as f:
        for line in f:
            if line.strip().startswith("MM") or line.strip().startswith("QM"):
                data = line.split()
                try:
                    start_idx = np.abs(int(data[1]))
                    end_idx = int(data[2])
                    tmpidxs = list(range(start_idx, end_idx + 1))
                    idxs.extend(tmpidxs)
                except:
                    pass
            elif line.strip().startswith("LA"):
                data = line.split()
                tmpidx = int(data[-1])
                idxs.extend([ tmpidx ])

    return idxs

def tk2xyz(tk, key):
    '''
    '''
    idxs = key_to_idxs(key)
    idxs = np.array(idxs) - 1
    atoms, structure = extract(tk)
    atoms = [ atoms[i] for i in idxs ]
    structure = structure[idxs]

    data = []
    for i, atom in enumerate(atoms):
        if atom == "L":
            atom = "H"

        data.append(tuple( [ atom ] + structure[i].tolist() ))

    return data[9:]

if __name__ == '__main__':

        pass
