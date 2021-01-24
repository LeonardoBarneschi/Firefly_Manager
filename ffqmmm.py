#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import MDAnalysis as mda
import os
import shutil
import sys
import warnings

from subprocess import check_output
from tkff2chg import tkff2chg
from tk2xyz import tk2xyz 


if not sys.warnoptions:
            warnings.simplefilter("ignore")

basis_set = ['6-31gd', 'cc-pvtz']
leveloftheory = ['casscf', 'xmcqdpt2']

def options():

    parser = argparse.ArgumentParser()

    parser.add_argument('-lot', 
                        choices=leveloftheory,
                        required=True,
                        nargs='?',
                        help='Level of theory')

    parser.add_argument('-bas', 
                        choices=basis_set, 
                        required=True,
                        nargs='?',
                        help='Basis set')

    args = vars(parser.parse_args())

    return args

def flatten(lst):

    flattened = sum( ([x] if not isinstance(x, list)
                     else flatten(x) for x in lst), [] )

    return flattened

def checkey(key):
    fstln = 'parameters melacu63.prm'

    if not os.path.exists(key):
        print(f'File {key} was not found. Aborting.')
        return False 

    with open(key,'r') as k:
        lns = k.read().splitlines()
        if lns[0] != fstln:
            print(f'File {key} not in the expected format.')
            return False
    
        else:
            return

if __name__ == "__main__":

    args = options()
    lot = args['lot']
    bas = args['bas']
    
    fftemps = os.environ.get('FFTEMPS')
    method = f'{lot}-{bas}.inp'

    currpath = os.getcwd()
    tklst = [f for f in os.listdir(currpath) if f.endswith('.xyz')]
    keylst = [f for f in os.listdir(currpath) if f.endswith('.key')]
  
    if len(tklst) == 1:
        try:
            # check for format (MDA-style)
            u = mda.Universe(tklst[0], topology_format='TXYZ')
        except:
            print(f'File {tklst[0]} not in the proper format')
            sys.exit()
        
        tk = tklst[0]
    else:
        sys.exit()

    if len(keylst) == 1 :
        try:
            checkey(keylst[0])
        except:
           print(f'File {keylst[0]} not in the expected format.')
           sys.exit()
        key = keylst[0]
    
    # Quick input name generation.
    inpname = tk.split('.')[0]+'.inp'

    # Generate xyz qm coords for Firefly.
    qm = tk2xyz(tk, key)

    # Generate fragment potential
    frg = tkff2chg(tk, key)

    prodict = { 'H' : 1.0, 'C' : 6.0, 'N' : 7.0 }
    newgeom = [ '$DATA\n qmmm chgs cloud\n C1' ]
    elem = [ x[0] for x in qm ]
    elenum = [ x[0]+str(i+1) for i, x in enumerate(qm) ]
    coords = [ list(map(float, x[1:])) for x in qm ]
    atomic = list(map(prodict.get, elem))
    row = list(zip(elenum, atomic, coords))

    oldinpath = os.path.join(fftemps, method)
    newinpath = os.path.join(currpath, inpname)
    
    if not os.path.exists(newinpath):

        shutil.copy(oldinpath, newinpath) 

        with open(newinpath, 'a') as n:

            n.write(' $DATA\n')
            n.write(' point charges cloud\n')
            n.write(' C1\n')
            for a in row:
                n.write(f' {a[0]:<5}{a[1]:<7}{a[2][0]:^14.6f}{a[2][1]:^14.6f}{a[2][2]:^14.6f}\n')

            n.write(' $END\n')
            orb = os.path.join(currpath, 'PUNCH')
            
            if os.path.exists(orb):
                with open(orb, 'r') as o:
                    for l in o:
                        n.write(l)

            if os.path.exists('MM_charges'):
                with open('MM_charges', 'r') as m:
                    for r in m:
                        n.write(r)

