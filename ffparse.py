#!/usr/bin/env python3

import csv
import os
import re
import sys
import numpy as np
import pandas as pd

def skiplines(openfile, nlines=0):
    '''Skips nlines + 1 lines in openfile. In other words, if nlines=0 it will
    go to the next line.'''

    for i in range(nlines):
        next(openfile)

    return next(openfile)

def parse_firefly(ffout):
    '''
    ffout: str or path obj.
       Firefly (PC GAMESS) output file.
    '''
    casscf = []
    xmcqdpt2 = []
    charges = []

    casscfen = "E(MCSCF)"
    qdpt2en = "XMC-QDPT2 ENERGIES"
    pop = "TOTAL MULLIKEN"

    with open(ffout, 'r') as f:
        for line in f:

            # Energies
            if qdpt2en.lower() in line.lower():
                line = skiplines(f, 2)

                while len(line.split()) == 5:

                    # 2nd Order
                    xmcqdpt2.append(float(line.split()[-1]))

                    # 1st Order
                    casscf.append(float(line.split()[2]))
                    line = skiplines(f)

            # Charges / Population
            if pop.lower() in line.lower():

                line = skiplines(f, 1)

                while len(line.split()) == 6:

                    # Mulliken
                    charges.append(float(line.split()[3]))
                    line = skiplines(f)

    casscf = np.array(casscf)
    xmcqdpt2 = np.array(xmcqdpt2)
    charges= np.array(charges)

    charges = charges.reshape(-1, casscf.shape[0], order='F')

    return casscf, xmcqdpt2, charges

def parsedconcat(casscf, xmcqdpt2, charges):

    roots = charges.shape[1]
    idxs = [0, 1, 2, 3, 4, 18, 19, 36, 37, 53]

    k1 = [f'CASSCF S{i} (a.u.)' for i in range(roots)]
    k2 = [f'XMCQDPT2 S{i} (a.u.)' for i in range(roots)]
    #k4 = [f'XMCQDPT2 S{i} CHARGE' for i in range(roots)]
    k3 = [f'XMCQDPT2 S{i} PSB CHARGE' for i in range(roots)]

    df1 = pd.DataFrame( [ dict(zip(k1, casscf)) ] )
    df2 = pd.DataFrame( [ dict(zip(k2, xmcqdpt2)) ] )
    #df4 = pd.DataFrame( [ dict(zip(k4, charges.T)) ] )
    df3 = pd.DataFrame( [ dict(zip(k3, charges.T[:, idxs].sum(axis=1) )) ] )

    df = pd.concat( [df1, df2, df3 ], axis=1 )

    return df

if __name__ == "__main__":

    path = os.getcwd()
    paths = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.endswith(".out"):
                 paths.append(os.path.join(root, f))

    dfs = []
    for p in paths: 
        dfs.append(parsedconcat(*parse_firefly(p)))
   
    dfparsed = pd.concat(dfs).reset_index(drop=True)
    dfparsed = dfparsed.sort_values( 'CASSCF S0 (a.u.)' )
  
    pattern = re.compile(r'CASSCF S\d{1}')
    entries = list(dfparsed.columns.values)
    roots = len([entry for entry in entries if pattern.search(entry)])

    k1 = [f'CASSCF S{i} (kcal)' for i in range(roots)]
    k2 = [f'XMCQDPT2 S{i} (kcal)' for i in range(roots)]
    df1 = (dfparsed.iloc[:, :roots] - dfparsed.iloc[:, 0].min()) * 627.509
    df2 = (dfparsed.iloc[:, :roots] - dfparsed.iloc[:, 0].min()) * 627.509
    df1.columns = k1
    df2.columns = k2
  
    df = pd.concat( [dfparsed, df1, df2], axis=1 )
    df.to_csv('parsed.csv', quoting=csv.QUOTE_NONNUMERIC)

    print(df)
