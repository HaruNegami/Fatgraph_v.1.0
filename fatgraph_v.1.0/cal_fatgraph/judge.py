# coding: utf-8

import os
import numpy
import operate_file.operate_file as operate_file
import configuration.classFilePath as classFilePath
import error_check
import cal_fatgraph.normalize as normalize

def dot(vec1, vec2):
    return numpy.dot(vec1, vec2)

def cal_rotation(vec_i,vec_j, base_i):
    twist_info = "u"
    change_cyclic_order = 0
    y_dot_z = dot(base_i[1],base_i[2])
    len_y = numpy.linalg.norm(base_i[1])
    len_z = numpy.linalg.norm(base_i[2])
    rotation_range = dot(vec_i[1],vec_j[1]) + dot(vec_i[2],vec_j[2])
    if y_dot_z >= 0:
        if rotation_range <= 0:
            twist_info = "t"
    else:
        change_cyclic_order = 1
        if rotation_range >= 0:
            twist_info = "t"
    return twist_info, change_cyclic_order

def get_h_bond_coordinate(h_bond, orth_b, base):
    index = [0, 0]
    flag = 2
    for i, pair in enumerate(orth_b):
        print(pair)
        if pair[4]+1 == h_bond[0] and pair[5] == h_bond[1]:
            vec_i = pair
            base_i = base[i]
            index[0] = pair[7]
            flag -= 1
        if pair[4] == h_bond[2] and pair[5] == h_bond[3]:
            vec_j = pair
            index[1] = pair[7]
            flag -= 1
    if flag != 0:
        return 0,0,'no edge',0
    return vec_i, vec_j, base_i, index

def twisted_or_untwisted_for_h_bond(h_bonds, orth_b, base):
    twisted_info = []
    for h_bond in h_bonds:
        bond_pair = []
        vec_i, vec_j, base_i, index = get_h_bond_coordinate(h_bond, orth_b, base)
        if vec_i == 0:
            continue
        else:
            bond_pair = []
            bond_pair.extend(vec_i[3:])
            bond_pair.append(index[0])
            bond_pair.extend(vec_j[3:])
            bond_pair.append(index[1])
            bond_pair.extend(cal_rotation(vec_i, vec_j, base_i))
            twisted_info.append(bond_pair)
    return twisted_info

def twisted_or_untwisted_for_pep_bond(orth_b, basis):
    twisted_info = []
    count_chain_change = 0
    chain = orth_b[0][5]
    for i in range(len(orth_b)-1):
        bond_pair = []
        bond_pair.extend(orth_b[i][3:])
        bond_pair.extend(cal_rotation(orth_b[i], orth_b[i+1], basis[i]))
        twisted_info.append(bond_pair)
    return twisted_info


def twisted_or_untwisted_for_pep_bond_test(orth_b, basis):
    twisted_info = []
    count_chain_change = 0
    chain = orth_b[0][5]
    orth_i1 = [i[1] for i in orth_b[:-1]]
    orth_i2 = [i[1] for i in orth_b[:-1]]
    orth_j1 = [i[2] for i in orth_b[1:]]
    orth_j2 = [i[2] for i in orth_b[1:]]
    base_1 = [i[1] for i in basis[:-1]]
    base_2 = [i[2] for i in basis[:-1]]
    bond_pair = [i[3:] for i in orth_b[:-1]]
    bond_dot = numpy.einsum('ij, ij->i', base_1, base_2)
    orth_i1_dot = numpy.einsum('ij, ij->i', orth_i1, orth_j1)
    orth_i2_dot = numpy.einsum('ij, ij->i', orth_i2, orth_j2)
    orth_dot = [x+y for x, y in zip(orth_i1_dot, orth_i2_dot)]
    change_cyclic_order = 0
    for i in range(len(bond_pair)):
        if bond_dot[i] >= 0:
            if orth_dot[0] <= 0:
                bond_pair[i].extend(["t", 0])
        else:
            if rotation_range >= 0:
                bond_pair[i].extend(["t", 1])
            else:
                bond_pair[i].extend(["u", 1])
    return bond_pair


def rotations(extracted_data, hydrogen_bonds, file, threshold):
    twisted_data_dict = dict()
    if error_check.lack_of_data(extracted_data["pdb"], file, 1, "no atom data") == 0:
        return 0, 'no atom data'
    else:
        orth_b, basis = normalize.cal_base(extracted_data["pdb"], file)
        if orth_b == 0:
            return 0, 'cal base error'
        else:
            twisted_data_dict["peptide"] = twisted_or_untwisted_for_pep_bond(orth_b, basis)
            twisted_data_dict["hydrogen"] = twisted_or_untwisted_for_h_bond(hydrogen_bonds, orth_b, basis)
            if twisted_data_dict['hydrogen'] == 'no edge':
                return 0, 'no edge'
            twisted_data = twisted_data_dict["peptide"] + twisted_data_dict["hydrogen"]
            operate_file.print_output(twisted_data, file.path[file.twisted][threshold], file)
            return twisted_data_dict, 'OK'
