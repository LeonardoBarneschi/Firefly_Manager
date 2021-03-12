#!/usr/bin/python3

import collections
import os
import sys
from operator import itemgetter

def get_xyz_coordinates(file_xyz):
    """
    Initializes a dictionary where several entries are created. Those entries
    corresponds to selection that the user may wish to store and/or use. They
    are built retrieveng data from an xyz file.
    """
    tk  = 'file'
    c2r = 'xyz_atom_c2r_ret'
    c3r = 'xyz_atom_c3r_ret'
    hr  = 'xyz_atom_hr_ret'
    lah = 'xyz_atom_lah'

    xyz_tk_dict = {lah:[],hr:[],c3r:[],c2r:[],tk:[]}

    with open (file_xyz, "r") as xyz_tinker_template_file:
        for line in xyz_tinker_template_file:
            ls = line.split()

            xyz_tk_dict[tk].append(ls[0:])

            if "C3R" in line:
                xyz_tk_dict[c3r].append(list(ls[0:]))

            if "C2R" in line:
                xyz_tk_dict[c2r].append(list(ls[0:]))

            if "HR" in line:
                xyz_tk_dict[hr].append(list(ls[0:]))

            if "LAH" in line:
                xyz_tk_dict[lah].append(list(ls[0:]))

    return xyz_tk_dict

def get_xyz_coordinates_charge_vacuo(file_xyz, frag_cut):

    xyz_list = []
    chg_index = ()

    # >>>> Totally empiric, I hate it
    if frag_cut == 13:

        ### 51 atom file
        #chg_index = (0, 1, 2, 16, 17, 34, 35)
        ### 54 atom file
        chg_index = (0, 1, 2, 3, 4, 18, 19, 36, 37, 53)

    elif frag_cut == 11:
        ### 51 atom file
        chg_index = (0, 1, 2, 14, 15, 16, 17, 33, 34, 35, 22, 48, 49, 50)

    with open (file_xyz, "r") as xyz_vacuo_file:
        for line in xyz_vacuo_file:
            ls = line.split()
            xyz_list.append(ls[1:])

    xyz_list = xyz_list[2:]

    return [xyz_list[i] for i in chg_index]

def get_xyz_coordinates_vacuo(file_xyz):

    xyz_list = []

    # >>>> Totally empiric, I hate it

    ### Works with 51 atom file
    #tors_index = (16, 20)

    ### Works with 54 atom file
    tors_index = (18, 22)

    with open (file_xyz, "r") as xyz_vacuo_file:
        for line in xyz_vacuo_file:
            ls = line.split()
            xyz_list.append(ls[1:])


    return xyz_list[tors_index[0]:tors_index[1]]

def get_xyz_coordinates_bla_vacuo(file_xyz):

    xyz_list = []

    # >>>> Totally empiric, I hate it
    bla_index = (11, 22, 3, 2)

    with open (file_xyz, "r") as xyz_vacuo_file:
        for line in xyz_vacuo_file:
            ls = line.split()
            xyz_list.append(ls[1:])

    #print( *xyz_list[bla_index[0]:bla_index[1]]+[xyz_list[bla_index[2]]]+[xyz_list[bla_index[3]]] , sep='\n')

    return xyz_list[bla_index[0]:bla_index[1]]+[xyz_list[bla_index[2]]]+[xyz_list[bla_index[3]]]

def get_xyz_dictionary(file_xyz):
    '''
    Store in a dictionary the dihedral angles of
    the retinal pi-chain.
    '''

    xyz_dict = collections.OrderedDict()
    tmp_dict = get_xyz_coordinates(file_xyz)

    xyz_dict.update({
                    'c4':  [],
                    'c5':  [],
                    'c6':  [],
                    'c7':  [],
                    'c8':  [],
                    'c9':  [],
                    'c10': [],
                    'c11': [],
                    'c12': [],
                    'c13': [],
                    'c14': [],
                    'c15': [],
                    'n':   [],
                    'ce':  []
                    })

    dihed_list = select_dihedrals(tmp_dict)

    for i, key in enumerate(xyz_dict.keys()):
        xyz_dict[key].append(dihed_list[i][2:5])

    return xyz_dict

def select_qm_charge(xyz_tk_dict, cut=None):
    """
    The method "get_xyz_coordinates" creates a dictionary and entries taking data
    directly from an xyz file. Here values are taken from those entries to create
    new entries in the same dictionary, or eventually a new dictionary.
    """

    # new selection name
    mk = 'mk_espf_charges'

    # add new empty entry to the dictionary
    xyz_tk_dict[mk] = []

    #print(*xyz_tk_dict)
    # the specific task is in this case to create a susbet corresponding to the QM atoms
    # used for calculating Mulliken or Espf charges. The chosen reference atoms are here 
    # C14 and C15.
    incr = len(xyz_tk_dict['xyz_atom_c2r_ret'])-1
    xyz_tk_dict[mk].append(xyz_tk_dict['xyz_atom_c2r_ret'][incr-1])
    xyz_tk_dict[mk].append(xyz_tk_dict['xyz_atom_c2r_ret'][incr])

    # Tricky. In the tinker format, after the coordinates, connectivity is specified.
    # Here, basically, I start from C15 and I may select the following atom/s because 
    # those are specified at C15 entry. And so on for the following.

    if cut=='11':
        args = [(xyz_tk_dict,'file',mk,1,6),(xyz_tk_dict,'file',mk,2,6),(xyz_tk_dict,'file',mk,2,7),\
                (xyz_tk_dict,'file',mk,3,8),(xyz_tk_dict,'file',mk,3,9),(xyz_tk_dict,'file',mk,0,8),\
                (xyz_tk_dict,'file',mk,1,8),(xyz_tk_dict,'file',mk,0,6),(xyz_tk_dict,'file',mk,9,8),\
                (xyz_tk_dict,'file',mk,10,7),(xyz_tk_dict,'file',mk,10,8),(xyz_tk_dict,'file',mk,10,9),\
                (xyz_tk_dict,'file',mk,9,6),(xyz_tk_dict,'file',mk,14,8)]

    elif cut=='13':
        args = [(xyz_tk_dict,'file',mk,1,6),(xyz_tk_dict,'file',mk,2,6),(xyz_tk_dict,'file',mk,2,7),\
                (xyz_tk_dict,'file',mk,3,8),(xyz_tk_dict,'file',mk,3,9),(xyz_tk_dict,'file',mk,0,8),\
                (xyz_tk_dict,'file',mk,1,8)]

    for arg in args:
        iterate(*arg)

    if xyz_tk_dict['xyz_atom_lah']:
        xyz_tk_dict[mk].append(xyz_tk_dict['xyz_atom_lah'][0])

    qm_charges_list = xyz_tk_dict[mk]

    return qm_charges_list

def select_dihedrals(xyz_tk_dict):
    """
    Updates the dictionary with a new entry for atoms included
    in the dihedral calculation.
    """

    # >>>> The N of the proton Schiff base (psb) is specified (in Tinker format) 
    # >>>> at position "6" of the entry of the last carbon of the C2R chain.
    psb = 6

    # >>>> Initialize to empty "[]" slot for this new key.
    xyz_tk_dict['dihedrals'] = xyz_tk_dict['xyz_atom_c2r_ret']

    # >>>> the first C2R carbon
    ref_c2r = xyz_tk_dict['dihedrals'][0][6]

    # >>>> to include C4=C5 is necessary to grab 1 C3R atom.
    for c3r in xyz_tk_dict['xyz_atom_c3r_ret']:
        if c3r[0] == ref_c2r:
            xyz_tk_dict['dihedrals'] = [c3r] + xyz_tk_dict['xyz_atom_c2r_ret']
            break

    # >>>> To find the entry representing the last carbon of the C2R chain.
    last_elem = len(xyz_tk_dict['dihedrals'])-1

    # >>>> Args to be passed to "iterate" func. with tuple unpacking.
    args_dihed = [(xyz_tk_dict,'file','dihedrals',last_elem,psb), \
                  (xyz_tk_dict,'file','dihedrals',last_elem+1,psb)]

    # >>>> Tuple unpacking => 1
    for arg_dihed in args_dihed:
        iterate(*arg_dihed)

    dihedrals_chain_list = xyz_tk_dict['dihedrals']

    return dihedrals_chain_list

def select_reactive_torsion(xyz_tk_dict, config=None):
    """
    Given the Retinal configuration (9-11-13C or AT), find the
    corresponding atom labels in the Tinker-xyz file.
    """
    reactive_torsion = []

    # Initialize to empty "[]" slot for this new key.
    xyz_tk_dict['dihedrals'] = xyz_tk_dict['xyz_atom_c2r_ret']

    if config == '9C':
        # Torsional Angle relative to C9 = C10 bond.
        reactive_torsion.extend([xyz_tk_dict['dihedrals'][3]] + \
                                [xyz_tk_dict['dihedrals'][4]] + \
                                [xyz_tk_dict['dihedrals'][5]] + \
                                [xyz_tk_dict['dihedrals'][6]])
    if config == '11C':
        # Torsional Angle relative to C11 = C12 bond.
        reactive_torsion.extend([xyz_tk_dict['dihedrals'][5]] + \
                                [xyz_tk_dict['dihedrals'][6]] + \
                                [xyz_tk_dict['dihedrals'][7]] + \
                                [xyz_tk_dict['dihedrals'][8]])

    if config == 'AT' or config == '13C':
        # Torsional Angle relative to C13 = C14 bond.
        reactive_torsion.extend([xyz_tk_dict['dihedrals'][7]] + \
                                [xyz_tk_dict['dihedrals'][8]] + \
                                [xyz_tk_dict['dihedrals'][9]] + \
                                [xyz_tk_dict['dihedrals'][10]])

    return reactive_torsion

def select_ret_constr(xyz_tk_dict, config=None):
    """
    Given the Retinal configuration (9-11-13C or AT), find the
    atoms to be constrained in a relaxed scan calculation.

    :param xyz_tk_dict:

    # Import my package:return constraints: dictionary where keys are the atom forming
                         the dihedral angle and values are the
                         corresponding values.
    """

    xyz_tk_dict['config'] = xyz_tk_dict['xyz_atom_c2r_ret']

    n_constr = 4

    constraints = {}

    d1 = {}
    d2 = {}
    d3 = {}
    d4 = {}

    for i in range(1, n_constr+1):
        constraints['dihedral_{0}'.format(i)] = []

    if config == '13C-AT':

        # ---> Rc (C12-C13=C14-C15)
        d1['C12'] = [xyz_tk_dict['config'][7][2:5]]
        d1['C13'] = [xyz_tk_dict['config'][8][2:5]]
        d1['C14'] = [xyz_tk_dict['config'][9][2:5]]
        d1['C15'] = [xyz_tk_dict['config'][10][2:5]]

        # ---> Cc1 (C12-C13=C14-H)
        d2['C12'] = [xyz_tk_dict['config'][7][2:5]]
        d2['C13'] = [xyz_tk_dict['config'][8][2:5]]
        d2['C14'] = [xyz_tk_dict['config'][9][2:5]]
        d2['H14'] = [find_coordinates(xyz_tk_dict, 'file', 'config',9,8)[0][2:5]]

        # ---> Cc2 (C13=C14-C15-H)
        d3['CT13'] = [find_coordinates(xyz_tk_dict, 'file', 'config',8,8)[0][2:5]]
        d3['C13']  = [xyz_tk_dict['config'][8][2:5]]
        d3['C14']  = [xyz_tk_dict['config'][9][2:5]]
        d3['H14']  = [find_coordinates(xyz_tk_dict, 'file', 'config',9,8)[0][2:5]]

        # ---> Cc3 (H-C14-C15-H)
        d4['CT13'] = [find_coordinates(xyz_tk_dict, 'file', 'config',8,8)[0][2:5]]
        d4['C13']  = [xyz_tk_dict['config'][8][2:5]]
        d4['C14'] = [xyz_tk_dict['config'][9][2:5]]
        d4['C15'] = [xyz_tk_dict['config'][10][2:5]]

        constraints['dihedral_1'].append(d1)
        constraints['dihedral_2'].append(d2)
        constraints['dihedral_3'].append(d3)
        constraints['dihedral_4'].append(d4)

    #elif config == '11C-AT':
    #    # ---> dihedral 1 (C10-C11=C12-C13)
    #    d1['C10'] = [xyz_tk_dict['config'][5][2:5]]
    #    d1['C11'] = [xyz_tk_dict['config'][6][2:5]]
    #    d1['C12'] = [xyz_tk_dict['config'][7][2:5]]
    #    d1['C13'] = [xyz_tk_dict['config'][8][2:5]]

    #    # ---> dihedral 2 (C10-C11=C12-H12)
    #    d2['C10'] = [xyz_tk_dict['config'][5][2:5]]
    #    d2['C11'] = [xyz_tk_dict['config'][6][2:5]]
    #    d2['C12'] = [xyz_tk_dict['config'][7][2:5]]
    #    d2['H12'] = [find_coordinates(xyz_tk_dict, 'file', 'config',7,8)[0][2:5]]

    #    # ---> dihedral 3 (C11=C12-C13-CT13)
    #    d3['C11']  = [xyz_tk_dict['config'][6][2:5]]
    #    d3['C12']  = [xyz_tk_dict['config'][7][2:5]]
    #    d3['C13']  = [xyz_tk_dict['config'][8][2:5]]
    #    d3['CT13'] = [find_coordinates(xyz_tk_dict, 'file', 'config',8,8)[0][2:5]]

    #    # ---> dihedral 4 (H12-C12-C13-CT13)
    #    d4['H12']  = [find_coordinates(xyz_tk_dict, 'file', 'config',7,8)[0][2:5]]
    #    d4['C12']  = [xyz_tk_dict['config'][7][2:5]]
    #    d4['C13']  = [xyz_tk_dict['config'][8][2:5]]
    #    d4['CT13'] = [find_coordinates(xyz_tk_dict, 'file', 'config',8,8)[0][2:5]]

    #    constraints['dihedral_1'].append(d1)
    #    constraints['dihedral_2'].append(d2)
    #    constraints['dihedral_3'].append(d3)
    #    constraints['dihedral_4'].append(d4)

    #elif config == '9C-AT':
    #    # ---> dihedral 1 (C10-C11=C12-C13)
    #    d1['C8']  = [xyz_tk_dict['config'][3][2:5]]
    #    d1['C9']  = [xyz_tk_dict['config'][4][2:5]]
    #    d1['C10'] = [xyz_tk_dict['config'][5][2:5]]
    #    d1['C11'] = [xyz_tk_dict['config'][6][2:5]]

    #    # ---> dihedral 2 (C10-C11=C12-H12)
    #    d2['C8']  = [xyz_tk_dict['config'][3][2:5]]
    #    d2['C9']  = [xyz_tk_dict['config'][4][2:5]]
    #    d2['C10'] = [xyz_tk_dict['config'][5][2:5]]
    #    d2['H10'] = [find_coordinates(xyz_tk_dict,'file','config',5,8)[0][2:5]]

    #    # ---> dihedral 3 (C11=C12-C13-CT13)
    #    d3['C9']  = [xyz_tk_dict['config'][4][2:5]]
    #    d3['C10'] = [xyz_tk_dict['config'][5][2:5]]
    #    d3['C11'] = [xyz_tk_dict['config'][6][2:5]]
    #    d3['H11'] = [find_coordinates(xyz_tk_dict,'file','config',6,8)[0][2:5]]

    #    # ---> dihedral 4 (H12-C12-C13-CT13)
    #    d4['H10'] = [find_coordinates(xyz_tk_dict,'file','config',5,8)[0][2:5]]
    #    d4['C10'] = [xyz_tk_dict['config'][5][2:5]]
    #    d4['C11'] = [xyz_tk_dict['config'][6][2:5]]
    #    d4['H11'] = [find_coordinates(xyz_tk_dict,'file','config',6,8)[0][2:5]]

    #    constraints['dihedral_1'].append(d1)
    #    constraints['dihedral_2'].append(d2)
    #    constraints['dihedral_3'].append(d3)
    #    constraints['dihedral_4'].append(d4)

    return constraints

def select_hoop_atoms(xyz_tk_dict, config=None):
    """
    Given the Retinal configuration (9-11-13C or AT), find the
    corresponding - HOOP - atom labels in the Tinker-xyz file.
    """

    hoop_atoms = []

    # ---> Initialize to empty "[]" slot for this new key.
    xyz_tk_dict['dihedrals'] = xyz_tk_dict['xyz_atom_c2r_ret']

    if config == '9C':
        atom_1 = find_coordinates(xyz_tk_dict,'file','dihedrals',4,8)
        atom_2 = [xyz_tk_dict['dihedrals'][4]]
        atom_3 = [xyz_tk_dict['dihedrals'][5]]
        atom_4 = find_coordinates(xyz_tk_dict,'file','dihedrals',5,8)

    if config == '11C':
        atom_1 = find_coordinates(xyz_tk_dict,'file','dihedrals',6,8)
        atom_2 = [xyz_tk_dict['dihedrals'][6]]
        atom_3 = [xyz_tk_dict['dihedrals'][7]]
        atom_4 = find_coordinates(xyz_tk_dict,'file','dihedrals',7,8)

    if config == 'AT' or config == '13C':
        atom_1 = find_coordinates(xyz_tk_dict,'file','dihedrals',8,8)
        atom_2 = [xyz_tk_dict['dihedrals'][8]]
        atom_3 = [xyz_tk_dict['dihedrals'][9]]
        atom_4 = find_coordinates(xyz_tk_dict,'file','dihedrals',9,8)

    hoop_atoms = atom_1 + atom_2 + atom_3 + atom_4

    return hoop_atoms

def iterate(dic, list_name1, list_name2, ind1, ind2):
    """
    Requires the knowledge of tinker format. Grab a new
    atom/s based upon the connectivity of the input one.
    """
    # this entry of the dictionary (dic['file']) is the whole 
    # xyz_tinker file.

    for item in dic[list_name1]:
        if dic[list_name2][ind1][ind2] in item[0]:
            dic[list_name2].append(item)

def find_coordinates(dic, list_name1, list_name2, ind1, ind2):
    """
    Requires the knowledge of tinker format. Grab a new
    atom/s based upon the connectivity of the input one.
    """
    # this entry of the dictionary (dic['file']) is the whole
    # xyz_tinker file.

    to_find = []

    for item in dic[list_name1]:
        if dic[list_name2][ind1][ind2] in item[0]:
            to_find.append(item)
    return to_find

if __name__ == "__main__":

    file_xyz=sys.argv[1]
    cut = '13'

    #a = get_xyz_coordinates_charge_vacuo(file_xyz, cut)
    #print(a)
    # >>>> DEBUG (f: get_xyz_coordinates_charge_vacuo)
    #print(*get_xyz_coordinates_bla_vacuo(file_xyz), sep='\n')
    #print(*get_xyz_coordinates_charge_vacuo(file_xyz, 11), sep='\n')

    # >>>> DEBUG (f: get_xyz_coordinates)
    #print(get_xyz_coordinates(file_xyz))

    # >>>> DEBUG (f: get_xyz_dictionary)
    #xyz_dict = get_xyz_dictionary(file_xyz)

    # >>>> DEBUG (f: get_xyz_coordinates_vacuo)
    #get_xyz_coordinates_vacuo(file_xyz)

    # >>>> DEBUG (f: get_xyz_coordinates)
    tmp_dic = get_xyz_coordinates(file_xyz)
    mk_list = select_qm_charge(tmp_dic, cut)
    print(*mk_list, sep='\n')
