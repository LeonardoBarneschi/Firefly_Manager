#!/usr/bin/env python3

from pprint import pprint

import csv
import collections
import numpy as np
import sys
import time
import ff_get_xyz_labels as gxl

def bonds(c):
    """
    Calculates vector(bond) of 2 atoms.
    """
    v0 = c[0] - c[1]

    return np.linalg.norm(v0)

def angle(c):
    """
    Calculates angle of 3 atoms with the cosine formula.
    """
    # Find vector 1 and 2
    v0 = c[0] - c[1]
    v1 = c[2] - c[1]

    # Calculate dot and magnitudes product.
    w = np.dot(v0,v1)
    x = np.linalg.norm(v0)*np.linalg.norm(v1)
    cosθ = w / x

    return np.degrees(np.arccos(cosθ))

def dihed(c):
    """
    The formula has been taken from Wikipedia article on "Dihedral angle"
    """

    # If not would yield the supplementary.
    v0 = c[0] - c[1]
    v1 = c[2] - c[1]
    v2 = c[3] - c[2]

    # Two cross products to find the orthogonal vectors.
    v0_x_v1 = np.cross(v0, v1)
    v1_x_v2 = np.cross(v2, v1)

    # Cross product of former cross products.
    w = np.cross(v0_x_v1, v1_x_v2)
    y = np.dot(w, v1)*(1.0/np.linalg.norm(v1))
    x = np.dot(v0_x_v1, v1_x_v2)

    return np.degrees(np.arctan2(y, x))

def atoms_array_list(atomdict):
    """
    Given a dictionary with the structure:
    dic = {key:'[[xa,ya,za],[xb,yb,zb], .. [xn,yn,zn]]'}
    i.e. dic[key][0][0] would access "xa".
    where for all the keys I have the same number of nested
    entries "([[],[]])" creates a numpy array from the set
    of 3 coordinates. Hence 1 array == 1 atom.
    """
    list_of_array = []
    length = len(list(atomdict.values())[0])
    for i in range(0,length):
        for key in atomdict:
            xyz_coord = [float(j) for j in atomdict[key][i]]
            list_of_array.append(np.array(xyz_coord))

    return list_of_array

def bla_setup(lst):
    """
    Turns a list of coordinates into ordinate pairs
    corresponding to bonds.
    """
    c,j = 0,0
    bonds_list = []

    # upper bound for the loop.
    length = len(lst)

    for i in range(len(lst) - 2):
        bonds_list.append(np.array([lst[j], lst[j+1]]))
        j += 1

    return bonds_list

def calculate_bla(bonds_list):
    """
    1) Calculates full-BLA
    2) Calculates bion-BLA
    3) Calculates rpsb-BLA
    """
    bla_list = []

    # full BLA = (avg double bonds length - avg single bonds length)
    avg_double_bond = (bonds_list[0]+bonds_list[2]+bonds_list[4]\
                      +bonds_list[6]+bonds_list[8]+bonds_list[10])/6
    avg_single_bond = (bonds_list[1]+bonds_list[3]+bonds_list[5]\
                      +bonds_list[7]+bonds_list[9])/5

    full_bla = (avg_single_bond - avg_double_bond)

    # B-ionone BLA = (avg double bonds length  - avg single bonds length) until C13
    avg_double_bond_bion = (bonds_list[0]+bonds_list[2]+bonds_list[4]+bonds_list[6])/4
    avg_single_bond_bion = (bonds_list[1]+bonds_list[3]+bonds_list[5]+bonds_list[7])/4

    bion_bla = (avg_single_bond_bion - avg_double_bond_bion)

    # PSB BLA = (C14-C15) - (C15-N)
    psb_bla = (bonds_list[9] - bonds_list[10])

    bla_list.extend([full_bla, bion_bla, psb_bla])

    return bla_list

def calc(prm_lst, mode=None):
    """
    Does the calculation required depending on the
    keyword argument.
    """
    measure = []
    if mode == 'bonds':
        for i in range(len(prm_lst)):
            measure.append(bonds(prm_lst[i]))
    elif mode == 'angle':
        for i in range(len(prm_lst)):
            measure.append(angle(prm_lst[i]))
    elif mode == 'dihed':
        for i in range(len(prm_lst)):
            measure.append(dihed(prm_lst[i]))
    return measure

def return_ret_bla(file_xyz):
    """
    Calculates the bond-length-alternation.
    """
    mld_dic = gxl.get_xyz_dictionary(file_xyz)
    atoms = atoms_array_list(mld_dic)
    bonds_setup = bla_setup(atoms)
    bonds_values = calc(bonds_setup, 'bonds')
    bla_values = calculate_bla(bonds_values)

    return bla_values 

def return_ret_hoop_debug(file_xyz, config):

    tmp_dic = gxl.get_xyz_coordinates(file_xyz)
    dihed_list = []
    dihed_list = gxl.select_hoop_atoms(tmp_dic, config)

    tors_list = []
    for atom in dihed_list:
        np_atom = np.array(atom[2:5],dtype='float')
        tors_list.append(np_atom)

    tors_value = calc([tors_list], 'dihed')

    return at.formatter(tors_value, 3)

def return_ret_tors(file_xyz, config):

    tmp_dic = gxl.get_xyz_coordinates(file_xyz)
    dihed_list = []
    dihed_list = gxl.select_reactive_torsion(tmp_dic, config)
    tors_list = []
    for atom in dihed_list:
        np_atom = np.array(atom[2:5],dtype='float')
        tors_list.append(np_atom)

    tors_value = calc([tors_list], 'dihed')

    return tors_value

def main():

    file_xyz = sys.argv[1]
    config = '13C-AT'

    # ---> TESTING
    bla = return_ret_bla(file_xyz)

    print('Bond Length Alternation: ')
    print(bla)
    print('\n')

if __name__ == '__main__':

    main()
