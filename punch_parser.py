#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import re
import sys
import pandas as pd
import numpy as np

ang2au = 1/0.52917724924

def txt_parse(txtfile, pattern1, pattern2):
    '''
    Reads a text file and extract the text between
    two patterns as a list of lists (each row goes
    into a list).

    pattern1: upper bound line.
    pattern2: lower bound line.

    :returns txtlines: a list of lists.
    '''

    rf = 0
    gf = 0

    txt_lines = []
    if os.path.exists(txtfile):
        with open(txtfile, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if re.search(r'{0}'.format(pattern1), line):
                    gf = 1
                elif re.search(r'{0}'.format(pattern2), line) and gf == 1:
                    rf = 1
                if gf == 1 and rf == 0:
                    txt_lines.append([line])
    return txt_lines

def pandas_csv_to_txt(file_name, text):
    '''
    Read list of lists to pandas dataframe, then writes
    the content to a new text file.

    param file_name: type str, name of the new file.
    param text: type lst, text to be written fo file.

    returns pd_dataframe: content of txt file as Pandas Df.
    '''
    if os.path.exists(file_name):
        df = pd.DataFrame(text)
        with open(file_name,'a') as fn:
            to_write = df.to_csv(sep=' ',  index=False, header=False)
            to_write = to_write.replace('"','')
            fn.write(to_write)
    else:
        raise('The file already exists and I am not allowed to overwrite it. Bye Bye.')
        sys.exit(0)

def firefly_check_input(firefly_input):
    '''
    :param firefly_input: type str, name of the file input for FireFly.
    '''
    if os.path.exists(firefly_input):
        with open(firefly_input,'r') as fi:
            for line in fi:
                print(line)
    else:
        raise('The file {0} does not exists in the current path. Bye Bye.'.format(firefly_input))

def firefly_check_geometry(file_xyz):
    '''
    Check that the input xyz file is existing
    and that it is indeed a "pure" xyz file.
    '''
    qm_atoms = []

    if os.path.exists(file_xyz):
        # Standard RET+LYS systems have 63 atoms.
        qmmm_atoms = 63

        # The number of rows should thus be 65.
        rows = 0
        with open(file_xyz) as fx:
            for i, l in enumerate(fx):
                rows = i + 1
            if rows == 65:
                print('The xyz file is in the correct format: [✓]')
                return True
            else:
                raise('The xyz file has {0} atoms. It should be 63: [X]'.format(rows-2))
                return False
    else:
        raise('The file {0} does not exists in the current path.'.format(file_xyz))
        return False

def firefly_check_orbitals(file_orb):
    '''
    Check that the input PUNCH file is existing
    and that it has the correct format.
    '''
    flag = 0

    if os.path.exists(file_orb):
        with open(file_orb) as fo:
            for line in fo:
                if 'OPTIMIZED MCSCF MO-S' in line:
                    flag = 1
                    print('The orb file is in the correct format: [✓]')
                    return True
            if flag == 0:
                raise('The orbital file is corrupted or has an incorrect format: [X]')
                return False
    else:
        raise('The file {0} does not exists in the current path. Bye Bye.'.format(file_orb))
        return False

def firefly_check_mmcharge(file_chg):
    '''
    Check that the file including MM charges exists
    and that it has the correct format.
    '''
    flag = 0

    if os.path.exists(file_chg):
        with open(file_chg) as fc:
            for line in fc:
                if '$FRG001' in line:
                    flag = 1
                    print('The chg file is in the correct format: [✓]')
                    return True
            if flag == 0:
                raise('The MM charges file is corrupted or has an incorrect format: [X]')
                return False
    else:
        raise('The file {0} does not exists in the current path. Bye Bye.'.format(file_chg))
        return False

def prepare_firefly_qmmm_input(template, file_xyz, file_orb=None, file_chg=None):
    '''
    Prepare an input file for a "QMMM" calculation with
    FireFly.
    '''

    if firefly_check_geometry(file_xyz):
        geometry = txt_parse(file_xyz,'63','NOEND')[11:]
        newgeom = [' $DATA\n qmmm chgs cloud\n C1']
        anames = [' '+x[0].split()[0] for x in geometry]
        coords = [list(map(float, x[0].split()[1:])) for x in geometry]
        coords = np.array(coords)

        for i,n in enumerate(range(len(coords))):
            if anames[i].strip(' ') == 'C':
                m = ' 6.0 '
            elif anames[i].strip(' ') == 'H':
                m = ' 1.0 '
            elif anames[i].strip(' ') == 'N':
                m = ' 7.0 '

            n += 1
            c = '%14.6f' * len(coords[i].tolist())
            c = c % tuple(coords[i].tolist())
            newatom = anames[i] + '%i ' % n + m + c
            newgeom.append(newatom)

        geometry = newgeom + [' $END']
        pandas_csv_to_txt(template, geometry)

    if file_orb:
        if firefly_check_orbitals(file_orb):
            orbitals = txt_parse(file_orb,'OPTIMIZED','END')
            pandas_csv_to_txt(template, orbitals)

    if file_chg:
        if firefly_check_mmcharge(file_chg):
            mmcharge = txt_parse(file_chg, 'FRG001', 'NOEND')
            pandas_csv_to_txt(template, mmcharge)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--a',help='template', required=True)
    parser.add_argument('-b', '--b',help='geometry', required=True)
    parser.add_argument('-c', '--c',help='orbitals', required=False)
    parser.add_argument('-d', '--d',help='mmcharge', required=False)
    args = vars(parser.parse_args())

    template = args['a']
    file_xyz = args['b']
    file_orb = args['c']
    file_chg = args['d']
    prepare_firefly_qmmm_input(template, file_xyz, file_orb, file_chg)
