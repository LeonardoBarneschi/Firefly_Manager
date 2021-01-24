#!/usr/bin/env python

import numpy as np
import pandas as pd
import MDAnalysis as mda

import os
import shutil
import sys
import warnings
from subprocess import check_output

if not sys.warnoptions:
        warnings.simplefilter("ignore")

def linkidx(key):
    '''
    Parameters
    ----------
    key:
        file with QMMM specifics for Molcas.

    Returns
    -------
    idxs: list
        indices of atoms with special params.
    chgs: list
        charges of atoms with special params.

    '''

    idxs = []
    chgs = []
    if os.path.exists(key):
        with open(key ,'r') as k:
            for line in k:
                if line.lower().startswith('charge'):
                    idxs.append(int(line.split()[1][1:]))
                    chgs.append(line.split()[2])

    return idxs, chgs

def parse_tk(tk, ff, key=None):
    '''
    Parameters
    ----------
    tk: file obj.
        TXYZ file.
    ff: file obj.
        Force field Tinker format.

    Returns
    -------
    '''

    if os.path.exists(tk):
        u = mda.Universe(tk, topology_format='TXYZ')
    else:
        sys.exit()

    types = u.atoms.types.astype(int)
    coord = u.atoms.positions

    tkchg = []
    tkcls = []
    if os.path.exists(ff):
        with open(ff, 'r') as f:
            for line in f:
                if line.startswith('charge'):
                    tkcls.append(int(line.split()[1]))
                    tkchg.append(float(line.split()[2]))

    tkchgs = dict(zip(tkcls, tkchg))
    df = pd.DataFrame(types)
    df.index = np.arange(1, len(df) + 1)
    chgs = df.iloc[:,0].map(tkchgs).values

    if key:
        idxs, keychgs = linkidx(key)
        keychgs = list(map(float, keychgs))
        keychgs = dict(zip(idxs, keychgs))
        dfchgs = pd.DataFrame()
        dfchgs['idx'] = df.index
        dfchgs['chg'] = chgs
        dfchgs['check'] = dfchgs['idx'].map(keychgs)
        dfchgs.loc[dfchgs['check'].isnull(),'check'] = dfchgs['chg']
        charges = dfchgs['check'].values
        idxs = np.array(idxs).astype(int)
        coord = np.c_[coord, charges]
        coord = coord[~np.isnan(coord[:,-1])]
        coord = coord[coord[:,-1] != 0]

    return coord

def tkff2chg(tk, key):
    '''
    tk: path.
        Tinker-xyz file.
    key: path.
        key (Tinker) file.
    '''

    binpath = os.environ.get('PACKAGES')
    exe = os.path.join(binpath, 'prep.exe')
    ff  = os.path.join(binpath, 'melacu63.prm')
    coord = parse_tk(tk, ff, key)
    np.savetxt('MMpart_x_y_z_charge', coord, fmt='%14.6f')
    cmd = '{0}'.format(exe)
    out = check_output([cmd]).decode("utf-8") 
    os.remove('MMpart_x_y_z_charge')
    if not os.path.exists('MM_charges'):
        with open('MM_charges', 'w') as m:
            m.write(out)

    else:
        print('MM_charges already exists.')
    
    return out

if __name__ == "__main__":
    
    tk  = sys.argv[1]
    key = sys.argv[2]
    tkff2chg(tk, key)
