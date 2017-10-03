#!/usr/bin/env python

### Converts a pdb file (requires that the correct columns is used)
### to V3 format with correct naming based on stereo-chemistry and estimated bond-lengths
### assuming that N, C, CA, CB naming, as well as a few others, is correct.

import string
import sys
import numpy as np
import math

# maximum covalent bond lengths involving H, or not. Variables contains maximum bond length squared in AA^2
h_cov_bond_cutoff = 1.6
other_cov_bond_cutoff = 2.7
# Enable debug mode by setting this to True
debug = True


# Easily cut a series of characters from a string
def cut_string(string, begin, end):

        # Initialize the fragment
        fragment = ""

        # Make sure that last character is inside 
        # the string
        if end > len(string):
                end = len(string)

        # Assemble string       
        for i in range(begin, end+1):
                fragment = fragment + string[i-1]

        return fragment

# get atom coordinates
def get_coords(line):
    tokens = cut_string(line, 31, 54).split()
    return np.array([float(x) for x in tokens])

# get dihedral from coordinates
def get_dihedral(a,b,c,d):
    v1 = (b-a)
    v1/np.sqrt(v1.dot(v1))
    v2 = (c-b)
    v2/np.sqrt(v2.dot(v2))
    v3 = (d-c)
    v3/np.sqrt(v3.dot(v3))

    n1 = np.cross(v1,v2)
    n2 = np.cross(v2,v3)
    m1 = np.cross(n1,v2)
    return math.atan2(np.dot(m1,n2), np.dot(n1,n2))


# convert atom names to pdb3 format based on their coordinates
def convert_atom_names(line_data,index_dict, first_residue_number, second_residue_number):
    # for all residues
    for index_dict_key in index_dict.keys():
        assigned_coords = {}
        unassigned = {}
        # for all atom names in residue
        for name in index_dict[index_dict_key].keys():
            index = index_dict[index_dict_key][name]
            residue_name = line_data[index][2]
            coords = line_data[index][4]
            # do we believe the current naming?
            if name in trusted_name_dict[residue_name]:
                assigned_coords[name.strip()] = coords[:]
            else: # if not put in placeholder dictionary
                atom_type = line_data[index][5]
                if atom_type not in unassigned: unassigned[atom_type] = {}
                unassigned[atom_type][name] = coords[:]

        # Try to name any unassigned atoms
        if debug:
            print "REMARK Residue %s %d" % (residue_name, index_dict_key)
        assign_residue_atoms(assigned_coords, unassigned, connect[residue_name], line_data, index_dict, index_dict_key)

        # assign term atoms
        if index_dict_key == first_residue_number:
            assign_term(' H  ', 'H1', index_dict, index_dict_key, atom_type, assigned_coords)
            # X is just a dummy since it's not needed here
            new_assignments = False
            assign_name('H2', ['N', 'C', 'H1','H3'], assigned_coords, unassigned, connect[residue_name], line_data, index_dict, index_dict_key, new_assignments)
            assign_name('H3', ['N', 'C', 'H2','H1'], assigned_coords, unassigned, connect[residue_name], line_data, index_dict, index_dict_key, new_assignments)


def assign_residue_atoms(assigned_coords, unassigned,connections, line_data, index_dict, index_dict_key):
    new_assignments = False
    for name in connections.keys():
        if name in assigned_coords.keys():
            continue
        # only continue if there's a covalently bound atom already assigned
        for cov in connections[name]:
            new_assignments = assign_name(name, cov, assigned_coords, unassigned, connections, line_data, index_dict, index_dict_key, new_assignments)

    if new_assignments:
        assign_residue_atoms(assigned_coords, unassigned, connections, line_data, index_dict, index_dict_key)
    else:
        for name in connections.keys():
            if name not in assigned_coords.keys():
                # assume that HE2 and HD2 are missing due to protonation state, and only print warning for these
                # if debugging is enabled
                if debug or name not in ['HE2','HD2']:
                    print "REMARK Warning: %s does not exist in pdb. (Could just be due to protonation state)" % name

def assign_name(name, cov, assigned_coords, unassigned, connections, line_data, index_dict, index_dict_key, new_assignments):
    cov_name = cov[0]
    if cov_name not in assigned_coords.keys():
        return False
    atom_type = name[0]
    cov_atom_type = cov_name[0]
    #check if there's only 1 atom within a generous range for covalent bond length
    possible_assignments = []
    for uname in unassigned[atom_type]:
        distance_sq = sum((assigned_coords[cov_name] - unassigned[atom_type][uname])**2)
        if "H" in [atom_type,cov_atom_type]:
            cutoff = h_cov_bond_cutoff
        else:
            cutoff = other_cov_bond_cutoff
        
        #print distance_sq, name, cov_name, uname
        if distance_sq < cutoff:
            possible_assignments.append(uname)
    
    # assume that this is the correct covalently bound neighbour
    #print possible_assignments. If there's 3 possible atoms, then the specific assignment most likely won't matter (eg. HE1/2/3)
    if len(possible_assignments) in [1,3]:
        uname = possible_assignments[0]
        assign_single_name(uname, assigned_coords, name, unassigned, atom_type, index_dict, index_dict_key)
        return True
    # If there's 2 possible atoms to assign, then choose the correct one from stereochemistry.
    elif len(possible_assignments) == 2:
        #print possible_assignments, cov
        for uname_index, uname in enumerate(possible_assignments):
            # check for the correct stereo assignment. We only need two atoms connected to the neighbour to calculate this:
            # check if the neighbour is sp3
            if len(cov) == 4:
                if cov[1] in assigned_coords.keys() and cov[2] in assigned_coords.keys():
                    c1 = assigned_coords[cov[1]] - assigned_coords[cov[0]]
                    c2 = assigned_coords[cov[2]] - assigned_coords[cov[0]]
                elif cov[1] in assigned_coords.keys() and cov[3] in assigned_coords.keys():
                    c1 = assigned_coords[cov[1]] - assigned_coords[cov[0]]
                    c2 = assigned_coords[cov[3]] - assigned_coords[cov[0]]
                elif cov[2] in assigned_coords.keys() and cov[3] in assigned_coords.keys():
                    c1 = assigned_coords[cov[2]] - assigned_coords[cov[0]]
                    c2 = assigned_coords[cov[3]] - assigned_coords[cov[0]]
                else:
                    break
                if np.dot(np.cross(c1, c2), unassigned[atom_type][uname]-assigned_coords[cov[0]]) > 0:
                    assign_single_name(uname, assigned_coords, name, unassigned, atom_type, index_dict, index_dict_key)
                    new_assignments = True
                    break
            else: #sp2
                # If we are trying to assign HH11, then this will be the hydrogen atom where the
                # (HH11 - NH1) and (CZ - NE) dot product is maximized
                # The E-Z definitions of these atoms doesn't seem to be used in pdb3
                sec_cov_name = connections[cov_name][0][0]
    
                if sec_cov_name in assigned_coords.keys():
                    # We need the heavy atom linked by 3 covalent bonds
                    third_cov_name = connections[sec_cov_name][0][0]
                    if third_cov_name in assigned_coords.keys():
                        #print sec_cov_name, third_cov_name, uname, cov_name
                        dihedral = get_dihedral(assigned_coords[third_cov_name],
                                                assigned_coords[sec_cov_name],
                                                assigned_coords[cov_name],
                                                unassigned[atom_type][uname])
                        #v1 = (assigned_coords[sec_cov_name] - assigned_coords[third_cov_name])
                        #v1 = v1/np.sqrt(v1.dot(v1))
                        #v2 = (unassigned[atom_type][uname] - assigned_coords[cov_name])
                        #v2 = v2/np.sqrt(v2.dot(v2))
                        #v3 = (unassigned[atom_type][possible_assignments[uname_index+1]] - assigned_coords[cov_name])
                        #v3 = v3/np.sqrt(v3.dot(v3))
                        #first_dot = np.dot(v1,v2)
                        #second_dot = np.dot(v1,v3)
                        # 'Z' and 'E' is used loosely here
                        if abs(dihedral) < math.pi/2:
                            Z_conformer_name = uname
                            E_conformer_name = possible_assignments[uname_index+1]
                        else:
                            E_conformer_name = uname
                            Z_conformer_name = possible_assignments[uname_index+1]
    
                        # The specific aminoacids is handles individually here, since 1) the pdb format is different from iupac
                        # and 2) We avoid iterating down the chain to check the priority since.
                        # For iupac we can remove "NH2".
                        if cov_name in ["NH2", "NE2"]: 
                            Z_conformer_name, E_conformer_name = E_conformer_name, Z_conformer_name
    
                        assign_single_name(Z_conformer_name, assigned_coords, name, unassigned, atom_type, index_dict, index_dict_key)
                        new_assignments = True
                        break
                else:
                    break
    
    
    return new_assignments

def assign_single_name(uname, assigned_coords, name, unassigned, atom_type, index_dict, index_dict_key):
    assigned_coords[name] = unassigned[atom_type][uname][:]
    unassigned[atom_type].pop(uname, None)
    index = index_dict[index_dict_key][uname]
    if len(name) < 4:
        name = " " + name
    if len(name) < 4:
        name += " "
    if len(name) < 4:
        name += " "
    line_data[index][1] = name
    index_dict[index_dict_key][name] = index_dict[index_dict_key][uname]
    if debug:
        if uname != name:
           print "REMARK",uname, "assigned as", name

def assign_term(old_name, new_name, index_dict, index_dict_key, atom_type, assigned_coords):
    assigned_coords[new_name] = assigned_coords[old_name.strip()][:]
    index = index_dict[index_dict_key][old_name]
    if len(new_name) < 4:
        new_name = " " + new_name
    if len(new_name) < 4:
        new_name += " "
    if len(new_name) < 4:
        new_name += " "
    line_data[index][1] = new_name
    if debug:
        if old_name != new_name:
            print "REMARK",old_name, "assigned as", new_name




# Get filename from command line
if (len(sys.argv) != 2):
        if  (len(sys.argv) > 2):
                print "Too many arguments!"
        print "Usage: $ ./convert_pdb pdbfile.pdb [> converted_pdbfile.pdb]"
        exit()

filename = sys.argv[1]


# Test if file exists
try:
        file_open = open(filename)
except IOError as e:
        print "Could not open file: ", filename

file_text = file_open.readlines()


# Create dictionary with atomnames that are assumed correct
trusted_name_dict = {
            "ALA": [" HA "," C  "," CA "," CB "," N  "," O  "],
            "ARG": [" HA "," HE "," C  "," CA "," CB "," CD "," CZ "," CG "," N  "," NE "," NH1"," NH2"," O  "],
            "ASP": [" HA "," C  "," CA "," CB "," CG "," N  "," O  "," OD1"," OD2"],
            "ASN": [" HA "," C  "," CA "," CB "," CG "," N  "," ND2"," O  "," OD1"],
            "CYS": [" HA "," C  "," CA "," CB "," N  "," O  "," SG "],
            "GLU": [" HA "," C  "," CA "," CB "," CG "," CD "," N  "," O  "," OE1"," OE2"],
            "GLN": [" HA "," C  "," CA "," CB "," CG "," CD "," N  "," NE2"," O  "," OE1"],
            "GLY": [" C  "," CA "," N  "," O  "],
            "HIS": [" HA "," HD1"," HD2"," HE1"," C  "," CA "," CB "," CG "," CD2"," CE1"," N  "," ND1"," NE2"," O  "],
            "ILE": [" HA "," HB "," C  "," CA "," CB "," CG1"," CG2"," CD1"," N  "," O  "],
            "LEU": [" HA "," HG "," C  "," CA "," CB "," CG "," CD1"," CD2"," N  "," O  "],
            "LYS": [" HA "," C  "," CA "," CB "," CG "," CD "," CE "," N  "," NZ "," O  "],
            "MET": [" HA "," C  "," CA "," CB "," CG "," CE "," N  "," O  "," SD "],
            "PHE": [" HA "," HD1"," HD2"," HE1"," HE2"," HZ "," C  "," CA "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," N  "," O  "],
            "PRO": [" HA "," C  "," CA "," CB "," CG "," CD "," N  "," O  "],
            "SER": [" HA "," C  "," CA "," CB "," N  "," O  "," OG "],
            "THR": [" HA "," HB "," HG1"," C  "," CA "," CB "," CG2"," N  "," O  "," OG1"],
            "TRP": [" HA "," HD1"," HE3"," HZ2"," HZ3"," HH2"," C  "," CA "," CB "," CG "," CD1"," CD2"," CE2"," CE3"," CZ2"," CZ3"," CG2"," N  "," NE1"," O  "],
            "TYR": [" HA "," HB "," C  "," CA "," CB "," N  "," O  "],
            "VAL": [" HA "," HB "," C  "," CA "," CB "," N  "," O  "]
            }

# Connection definitions in pdb_v3. (Added newlines to avoid issues with vim syntax highlighting)
connect = {'CYS': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'CB': [['CA', 'HA', 'C', 'N'], ['SG', 'HG']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'SG']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'H': [['N', 'CA', 'pC']], 'HA': [['CA', 'N', 'C', 'CB']], 'SG': [['CB', 'CA', 'HB2', 'HB3']], 'HG': [['SG', 'CB']], 'HB3': [['CB', 'CA', 'SG', 'HB2']], 'HB2': [['CB', 'HB3', 'SG', 'CA']]},
 'ASP': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'HD2': [['OD2', 'CG']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'OD2', 'OD1']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['OD2', 'HD2']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'OD1': [['CG', 'CB', 'OD2']], 'H': [['N', 'CA', 'pC']], 'HA': [['CA', 'N', 'C', 'CB']], 'OD2': [['CG', 'OD1', 'CB']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]},
 'SER': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'OG': [['CB', 'CA', 'HB2', 'HB3']], 'CB': [['CA', 'HA', 'C', 'N'], ['OG', 'HG']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'OG']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'H': [['N', 'CA', 'pC']], 'HA': [['CA', 'N', 'C', 'CB']], 'HG': [['OG', 'CB']], 'HB3': [['CB', 'CA', 'OG', 'HB2']], 'HB2': [['CB', 'HB3', 'OG', 'CA']]},
 'GLN': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'HG3', 'HG2', 'CD']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['CD', 'NE2', 'OE1']], 'HG3': [['CG', 'CB', 'CD', 'HG2']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CD': [['CG', 'CB', 'HG2', 'HG3'], ['NE2', 'HE22', 'HE21']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'H': [['N', 'CA', 'pC']], 'HE22': [['NE2', 'HE21', 'CD']], 'HE21': [['NE2', 'CD', 'HE22']], 'HA': [['CA', 'N', 'C', 'CB']], 'NE2': [['CD', 'OE1', 'CG']], 'OE1': [['CD', 'CG', 'NE2']], 'HG2': [['CG', 'HG3', 'CD', 'CB']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]},
 'LYS': {'HA': [['CA', 'N', 'C', 'CB']], 'HE2': [['CE', 'HE3', 'NZ', 'CD']], 'HE3': [['CE', 'CD', 'NZ', 'HE2']], 'HG2': [['CG', 'HG3', 'CD', 'CB']], 'HG3': [['CG', 'CB', 'CD', 'HG2']], 'NZ': [['CE', 'CD', 'HE2', 'HE3']], 'HZ1': [['NZ', 'CE', 'HZ2', 'HZ3']], 'HZ3': [['NZ', 'CE', 'HZ1', 'HZ2']], 'HZ2': [['NZ', 'HZ3', 'HZ1', 'CE']], 'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'HG3', 'HG2', 'CD']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['CD', 'HD3', 'HD2', 'CE']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CE': [['CD', 'CG', 'HD2', 'HD3'], ['NZ', 'HZ3', 'HZ2', 'HZ1']], 'CD': [['CG', 'CB', 'HG2', 'HG3'], ['CE', 'HE3', 'HE2', 'NZ']], 'HD3': [['CD', 'CG', 'CE', 'HD2']], 'HD2': [['CD', 'HD3', 'CE', 'CG']], 'H': [['N', 'CA', 'pC']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]},
 'ASN': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'ND2': [['CG', 'OD1', 'CB']], 'HD22': [['ND2', 'HD21', 'CG']], 'HD21': [['ND2', 'CG', 'HD22']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'ND2', 'OD1']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['ND2', 'HD22', 'HD21']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'OD1': [['CG', 'CB', 'ND2']], 'H': [['N', 'CA', 'pC']], 'HA': [['CA', 'N', 'C', 'CB']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]},
 'PRO': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCD']], 'HD3': [['CD', 'HD2', 'CG', 'N']], 'HD2': [['CD', 'N', 'CG', 'HD3']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'HG3', 'HG2', 'CD']], 'CA': [['N', 'H', 'CD'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'HG2': [['CG', 'HG3', 'CD', 'CB']], 'HG3': [['CG', 'CB', 'CD', 'HG2']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['CD', 'HD2', 'HD3', 'CG']], 'CD': [['N', 'CA', 'H'], ['CG', 'CB', 'HG2', 'HG3']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['CD', 'N', 'HD3', 'HD2']], 'HA': [['CA', 'N', 'C', 'CB']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]},
 'THR': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'HG22': [['CG2', 'HG23', 'HG21', 'CB']], 'CB': [['CA', 'HA', 'C', 'N'], ['OG1', 'HG1'], ['CG2', 'HG23', 'HG22', 'HG21']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB', 'CG2', 'OG1']], 'OG1': [['CB', 'CA', 'CG2', 'HB']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'HB': [['CB', 'CA', 'OG1', 'CG2']], 'CG2': [['CB', 'HB', 'OG1', 'CA']], 'HG1': [['OG1', 'CB']], 'H': [['N', 'CA', 'pC']], 'HG21': [['CG2', 'CB', 'HG22', 'HG23']], 'HG23': [['CG2', 'CB', 'HG21', 'HG22']], 'HA': [['CA', 'N', 'C', 'CB']]},
 'PHE': {'HZ': [['CZ', 'CE2', 'CE1']], 'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'HE2': [['CE2', 'CZ', 'CD2']], 'HD1': [['CD1', 'CE1', 'CG']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'CD2', 'CD1']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['CD1', 'HD1', 'CE1'], ['CD2', 'HD2', 'CE2']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CZ': [['CE1', 'CD1', 'HE1'], ['CE2', 'CD2', 'HE2']], 'CE2': [['CD2', 'CG', 'HD2'], ['CZ', 'CE1', 'HZ']], 'HD2': [['CD2', 'CE2', 'CG']], 'CD1': [['CG', 'CB', 'CD2'], ['CE1', 'HE1', 'CZ']], 'CE1': [['CD1', 'CG', 'HD1'], ['CZ', 'HZ', 'CE2']], 'H': [['N', 'CA', 'pC']], 'HE1': [['CE1', 'CZ', 'CD1']], 'HA': [['CA', 'N', 'C', 'CB']], 'CD2': [['CG', 'CD1', 'CB'], ['CE2', 'HE2', 'CZ']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]},
 'ALA': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'CB': [['CA', 'HA', 'C', 'N']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB2', 'HB3', 'HB1']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'H': [['N', 'CA', 'pC']], 'HA': [['CA', 'N', 'C', 'CB']], 'HB1': [['CB', 'CA', 'HB3', 'HB2']], 'HB3': [['CB', 'HB2', 'HB1', 'CA']], 'HB2': [['CB', 'CA', 'HB1', 'HB3']]},
 'HIS': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'HE2': [['NE2', 'CE1', 'CD2']], 'HD1': [['ND1', 'CE1', 'CG']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'CD2', 'ND1']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'HB2': [['CB', 'HB3', 'CG', 'CA']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['ND1', 'HD1', 'CE1'], ['CD2', 'HD2', 'NE2']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CE1': [['ND1', 'CG', 'HD1'], ['NE2', 'CD2', 'HE2']], 'H': [['N', 'CA', 'pC']], 'ND1': [['CG', 'CB', 'CD2'], ['CE1', 'HE1', 'NE2']], 'HE1': [['CE1', 'NE2', 'ND1']], 'HA': [['CA', 'N', 'C', 'CB']], 'NE2': [['CD2', 'CG', 'HD2'], ['CE1', 'ND1', 'HE1']], 'CD2': [['CG', 'ND1', 'CB'], ['NE2', 'HE2', 'CE1']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HD2': [['CD2', 'NE2', 'CG']]},
 'GLY': {'C': [['CA', 'N', 'HA2', 'HA3'], ['nN', 'nH', 'nCA']], 'H': [['N', 'CA', 'pC']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA3', 'HA2', 'C'], ['pC', 'pO', 'pCA']], 'HA2': [['CA', 'HA3', 'C', 'N']], 'HA3': [['CA', 'N', 'C', 'HA2']]},
 'ILE': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'HG22': [['CG2', 'HG23', 'HG21', 'CB']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG2', 'HG23', 'HG22', 'HG21'], ['CG1', 'HG13', 'HG12', 'CD1']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB', 'CG1', 'CG2']], 'CG2': [['CB', 'CA', 'CG1', 'HB']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CG1': [['CB', 'HB', 'CG2', 'CA'], ['CD1', 'HD11', 'HD13', 'HD12']], 'HB': [['CB', 'CA', 'CG2', 'CG1']], 'HG12': [['CG1', 'HG13', 'CD1', 'CB']], 'HG13': [['CG1', 'CB', 'CD1', 'HG12']], 'HD13': [['CD1', 'HD11', 'HD12', 'CG1']], 'CD1': [['CG1', 'CB', 'HG12', 'HG13']], 'H': [['N', 'CA', 'pC']], 'HG21': [['CG2', 'CB', 'HG22', 'HG23']], 'HD12': [['CD1', 'CG1', 'HD13', 'HD11']], 'HG23': [['CG2', 'CB', 'HG21', 'HG22']], 'HA': [['CA', 'N', 'C', 'CB']], 'HD11': [['CD1', 'CG1', 'HD12', 'HD13']]},
 'LEU': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'HD22': [['CD2', 'CG', 'HD21', 'HD23']], 'HD23': [['CD2', 'HD22', 'HD21', 'CG']], 'HD21': [['CD2', 'CG', 'HD23', 'HD22']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'HG', 'CD2', 'CD1']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['CD1', 'HD12', 'HD13', 'HD11'], ['CD2', 'HD22', 'HD23', 'HD21']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CD1': [['CG', 'CB', 'CD2', 'HG']], 'CD2': [['CG', 'HG', 'CD1', 'CB']], 'H': [['N', 'CA', 'pC']], 'HD13': [['CD1', 'HD12', 'HD11', 'CG']], 'HD12': [['CD1', 'CG', 'HD11', 'HD13']], 'HD11': [['CD1', 'CG', 'HD13', 'HD12']], 'HA': [['CA', 'N', 'C', 'CB']], 'HG': [['CG', 'CB', 'CD1', 'CD2']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]},
 'ARG': {'NE': [['CD', 'CG', 'HD2', 'HD3'], ['CZ', 'NH2', 'NH1']], 'HA': [['CA', 'N', 'C', 'CB']], 'HE': [['NE', 'CZ', 'CD']], 'HG2': [['CG', 'HG3', 'CD', 'CB']], 'HG3': [['CG', 'CB', 'CD', 'HG2']], 'HH22': [['NH2', 'CZ', 'HH21']], 'HH21': [['NH2', 'HH22', 'CZ']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'HG3', 'HG2', 'CD']], 'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'H': [['N', 'CA', 'pC']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['CD', 'HD3', 'HD2', 'NE']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CZ': [['NE', 'CD', 'HE'], ['NH1', 'HH11', 'HH12'], ['NH2', 'HH21', 'HH22']], 'NH1': [['CZ', 'NE', 'NH2']], 'NH2': [['CZ', 'NH1', 'NE']], 'CD': [['CG', 'CB', 'HG2', 'HG3'], ['NE', 'HE', 'CZ']], 'HD3': [['CD', 'CG', 'NE', 'HD2']], 'HD2': [['CD', 'HD3', 'NE', 'CG']], 'HH12': [['NH1', 'CZ', 'HH11']], 'HH11': [['NH1', 'CZ', 'HH12']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]},
 'TRP': {'HH2': [['CH2', 'CZ3', 'CZ2']], 'CZ2': [['CE2', 'NE1', 'CD2'], ['CH2', 'HH2', 'CZ3']], 'CZ3': [['CE3', 'CD2', 'HE3'], ['CH2', 'CZ2', 'HH2']], 'CD1': [['CG', 'CB', 'CD2'], ['NE1', 'HE1', 'CE2']], 'CD2': [['CG', 'CD1', 'CB'], ['CE2', 'CZ2', 'NE1'], ['CE3', 'HE3', 'CZ3']], 'HA': [['CA', 'N', 'C', 'CB']], 'HE1': [['NE1', 'CE2', 'CD1']], 'HE3': [['CE3', 'CZ3', 'CD2']], 'CH2': [['CZ2', 'CE2', 'HZ2'], ['CZ3', 'CE3', 'HZ3']], 'HZ3': [['CZ3', 'CH2', 'CE3']], 'HZ2': [['CZ2', 'CH2', 'CE2']], 'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'CD2', 'CD1']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['CD1', 'HD1', 'NE1'], ['CD2', 'CE3', 'CE2']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CE3': [['CD2', 'CE2', 'CG'], ['CZ3', 'HZ3', 'CH2']], 'CE2': [['CD2', 'CG', 'CE3'], ['NE1', 'CD1', 'HE1'], ['CZ2', 'HZ2', 'CH2']], 'HD1': [['CD1', 'NE1', 'CG']], 'H': [['N', 'CA', 'pC']], 'HB2': [['CB', 'HB3', 'CG', 'CA']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'NE1': [['CD1', 'CG', 'HD1'], ['CE2', 'CD2', 'CZ2']]},
 'VAL': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'CG2': [['CB', 'CA', 'CG1', 'HB']], 'HG22': [['CG2', 'HG23', 'HG21', 'CB']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG2', 'HG23', 'HG22', 'HG21'], ['CG1', 'HG12', 'HG13', 'HG11']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB', 'CG1', 'CG2']], 'HB': [['CB', 'CA', 'CG2', 'CG1']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CG1': [['CB', 'HB', 'CG2', 'CA']], 'HG11': [['CG1', 'CB', 'HG13', 'HG12']], 'HG12': [['CG1', 'CB', 'HG11', 'HG13']], 'HG13': [['CG1', 'HG12', 'HG11', 'CB']], 'H': [['N', 'CA', 'pC']], 'HG21': [['CG2', 'CB', 'HG22', 'HG23']], 'HG23': [['CG2', 'CB', 'HG21', 'HG22']], 'HA': [['CA', 'N', 'C', 'CB']]},
 'GLU': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'HE2': [['OE2', 'CD']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'HG3', 'HG2', 'CD']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'HG2': [['CG', 'HG3', 'CD', 'CB']], 'HG3': [['CG', 'CB', 'CD', 'HG2']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CD': [['CG', 'CB', 'HG2', 'HG3'], ['OE2', 'HE2']], 'H': [['N', 'CA', 'pC']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['CD', 'OE2', 'OE1']], 'HA': [['CA', 'N', 'C', 'CB']], 'OE2': [['CD', 'OE1', 'CG']], 'OE1': [['CD', 'CG', 'OE2']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]},
 'TYR': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'HE2': [['CE2', 'CZ', 'CD2']], 'HD1': [['CD1', 'CE1', 'CG']], 'OH': [['CZ', 'CE2', 'CE1']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'CD2', 'CD1']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['CD1', 'HD1', 'CE1'], ['CD2', 'HD2', 'CE2']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'CZ': [['CE1', 'CD1', 'HE1'], ['CE2', 'CD2', 'HE2'], ['OH', 'HH']], 'HH': [['OH', 'CZ']], 'HD2': [['CD2', 'CE2', 'CG']], 'CD1': [['CG', 'CB', 'CD2'], ['CE1', 'HE1', 'CZ']], 'CE1': [['CD1', 'CG', 'HD1'], ['CZ', 'OH', 'CE2']], 'H': [['N', 'CA', 'pC']], 'HE1': [['CE1', 'CZ', 'CD1']], 'CE2': [['CD2', 'CG', 'HD2'], ['CZ', 'CE1', 'OH']], 'HA': [['CA', 'N', 'C', 'CB']], 'CD2': [['CG', 'CD1', 'CB'], ['CE2', 'HE2', 'CZ']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]},
 'MET': {'C': [['CA', 'N', 'CB', 'HA'], ['nN', 'nH', 'nCA']], 'HE1': [['CE', 'SD', 'HE3', 'HE2']], 'HE2': [['CE', 'SD', 'HE1', 'HE3']], 'HE3': [['CE', 'HE2', 'HE1', 'SD']], 'CB': [['CA', 'HA', 'C', 'N'], ['CG', 'HG3', 'HG2', 'SD']], 'CA': [['N', 'pC', 'H'], ['C', 'nN', 'O'], ['CB', 'HB3', 'HB2', 'CG']], 'HG2': [['CG', 'HG3', 'SD', 'CB']], 'HG3': [['CG', 'CB', 'SD', 'HG2']], 'O': [['C', 'CA', 'nN']], 'N': [['CA', 'HA', 'CB', 'C'], ['pC', 'pO', 'pCA']], 'SD': [['CG', 'CB', 'HG2', 'HG3'], ['CE', 'HE2', 'HE3', 'HE1']], 'CE': [['SD', 'CG']], 'H': [['N', 'CA', 'pC']], 'CG': [['CB', 'CA', 'HB2', 'HB3'], ['SD', 'CE']], 'HA': [['CA', 'N', 'C', 'CB']], 'HB3': [['CB', 'CA', 'CG', 'HB2']], 'HB2': [['CB', 'HB3', 'CG', 'CA']]}}



#Determine first and last residue name (for detecting O -> OT1 conversion)
last_residue_number = 0
first_residue_number = np.inf

for line in file_text:
    if "ATOM" not in line:
        continue
    residue_number = int(cut_string(line, 24, 26))
    if residue_number < first_residue_number:
        first_residue_number = residue_number


    if residue_number > last_residue_number:
        last_residue_number = residue_number

# Read lines from input file and store for later use
line_data = []
index_dict = {}
i = 0
for line in file_text:
    if "ATOM" not in line:
        continue

    # Read information from lines
    #atom_number = cut_string(line, 1, 12)
    atom_number = line[0:12]
    #atom_name = cut_string(line, 13, 16)
    atom_name = line[12:16]
    #residue_name = cut_string(line, 18, 20)
    residue_name = line[17:20]
    #residue_number = int(cut_string(line, 24, 26))
    residue_number = int(line[23:26].split()[0])
    if len(line) > 77 and line[77] in ["O","S","C","N","H"]:
        atom_type = line[77]
    elif "O" in atom_name:
        atom_type = "O"
    elif "S" in atom_name:
        atom_type = "S"
    elif "C" in atom_name:
        atom_type = "C"
    elif "N" in atom_name and "HN" not in atom_name:
        atom_type = "N"
    else:
        atom_type = "H"
    #rest = cut_string(line, 21, 1000)
    rest = line[20:]
    # if no chain is given, set to A
    if rest[1] == " ":
        rest = rest[0] + "A" + rest[2:]
    coords = get_coords(line)
    # store for later
    line_data.append([atom_number, atom_name, residue_name, residue_number,
                     coords, atom_type, rest])

    # create dictionary for mapping residue number and atom name to line_data.
    if residue_number not in index_dict: index_dict[residue_number] = {}

    index_dict[residue_number][atom_name] = i
    i += 1

convert_atom_names(line_data,index_dict, first_residue_number, last_residue_number)


for item in line_data:
    # Print the shebang
    atom_number = item[0]
    atom_name = item[1]
    residue_name = item[2]
    rest = item[-1]
    print atom_number + atom_name + " "+ residue_name + rest,

