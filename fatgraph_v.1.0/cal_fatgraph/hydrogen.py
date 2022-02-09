# coding: utf-8

import operate_file.operate_file as operate_file
import configuration.configuration as configuration
import shutil
import  error_check


def get_distance_of_two_CAs(axes_0, axes_1):
    distance = 0.0
    for i in range(3):
        distance += (axes_0[i] - axes_1[i]) * (axes_0[i] - axes_1[i])
    return distance


def get_bond_candidate(data_dssp, threshold):
    # thresholdごとにまとめる
    # 閾値までエッジを追加する

    atomNH = [] 
    atomO = [] 
    bond_list_pre = []
    dssp_id_with_residue_and_chain =[]
    for j in range(1,5):
        for residue in data_dssp:
            if residue == "no data":
                continue
            else:
                if j == 1:
                    dssp_id_with_residue_and_chain.append([residue[0],residue[5],residue[6]])
                if residue[j][2] == 0 and residue[j][0] != 0 and residue[j][1]<= threshold:
                    dssp_id_0 = residue[0]
                    dssp_id_1 = dssp_id_0 + residue[j][0]
                    if (dssp_id_1 <= 0 or len(data_dssp) <= dssp_id_1) or data_dssp[dssp_id_1 - 1] == "no data":
                        continue
                    else:
                        distance_of_two_CAs = get_distance_of_two_CAs(residue[-1], data_dssp[dssp_id_1 - 1][-1])
                        if j%2 ==1:
                            bond_list_pre.append([dssp_id_0,dssp_id_1,residue[j][1], distance_of_two_CAs])
                            if dssp_id_0 not in atomNH:
                                atomNH.append(residue[0])
                        else:
                            bond_list_pre.append([dssp_id_1,dssp_id_0,residue[j][1], distance_of_two_CAs])
                            if dssp_id_0 not in atomO:
                                atomO.append(residue[0])
    bond_list = get_bond_list(bond_list_pre)
    atom_dict = {"NH": atomNH, "O": atomO}
    return bond_list, atom_dict , dssp_id_with_residue_and_chain

def get_bond_list(bond_list_pre):
    bond_list = []
    for items in bond_list_pre:
        if items not in bond_list:
            bond_list.append(items)
    return bond_list

def select_stronger_bond(atom, bond_list, atom_type):
    bond_list_selected = []
    for item in atom:
        min = 100.0
        tempList = []
        for items in bond_list:
            if item == items[atom_type]:
                tempList.append(items)
        for itemTemplist in tempList:
            if itemTemplist[2] < min:
                min = itemTemplist[2]
        for item in tempList:
            if item[2] == min:
                bond_list_selected.append(item)
    return bond_list_selected


def unique(array):
    unique_array = []
    for item in array:
        if item not in unique_array:
            unique_array.append(item)
    return unique_array

def get_max_bond_id(bond_list_selected):
    max_bond_id = 0
    for bond_pair in bond_list_selected:
        if bond_pair[0] > max_bond_id or bond_pair[1] > max_bond_id:
            max_bond_id = max(bond_pair[0],bond_pair[1])
    return max_bond_id

def get_duplicated_bonds(max_bond_id, bond_list_selected):
    duplicated_bonds = []
    for i in range(max_bond_id+1):
        duplicated_bonds.append([0,0,[],100.0,100.0])
    ## duplicated_bonds はリスト [Nに結合しているbondの数, Oに結合しているbondの数, 結合しているbondのリスト, それぞれのbond間の距離]
    ## 100.0は十分遠いことを示す
    for bond_pair in bond_list_selected:
        if bond_pair not in duplicated_bonds[bond_pair[0]][2]:
            duplicated_bonds[bond_pair[0]][0] += 1
            duplicated_bonds[bond_pair[0]][2].append(bond_pair)
        if bond_pair not in duplicated_bonds[bond_pair[1]][2]:
            duplicated_bonds[bond_pair[1]][2].append(bond_pair)
            duplicated_bonds[bond_pair[1]][1] += 1
    return duplicated_bonds


def check_uniqueness_of_bonds_with_lowest_potential(duplicated_bonds, unique_bond_list):
    for i in range(len(duplicated_bonds)):
        for j in range(2):
            if duplicated_bonds[i][j] >= 2:
                for k in range(len(duplicated_bonds[i][2])):
                    if duplicated_bonds[i][j] >= 2 and duplicated_bonds[i][2][k][j] == i:
                        if duplicated_bonds[i][2][k][2] > duplicated_bonds[i][3+j]:
                            if duplicated_bonds[i][2][k] in unique_bond_list:
                                unique_bond_list.remove(duplicated_bonds[i][2][k])
                            duplicated_bonds[i][2][k].append("delete")
                            duplicated_bonds[i][j] -= 1
    return(duplicated_bonds, unique_bond_list)

def get_hydrogen_bond_list_filtered_by_potential(bond_list, atom_dict):
    bond_list_selected_O = select_stronger_bond(atom_dict["O"], bond_list, 1)
    bond_list_selected_NH = select_stronger_bond(atom_dict["NH"], bond_list, 0)
    bond_list_selected = bond_list_selected_O + bond_list_selected_NH
    unique_bond_list = unique(bond_list_selected)
    max_bond_id = get_max_bond_id(bond_list_selected)
    duplicated_bonds = get_duplicated_bonds(max_bond_id, unique_bond_list)
    duplicated_bonds, unique_bond_list = check_uniqueness_of_bonds_with_lowest_potential(duplicated_bonds, unique_bond_list)
    return duplicated_bonds, unique_bond_list, max_bond_id

def get_hydrogen_bond_list_filtered_by_potential_and_distance(data,duplicated_bonds,unique_bond_list,max_bond_id, file):
    bond_list = []
    file_dssp = file.path[file.dssp]
    file_with_flaw = file.path[file.flawed]
    for i in range(len(duplicated_bonds)):
        for j in range(2):
            if duplicated_bonds[i][j] == 1:
                if duplicated_bonds[i][2][0][:3] not in bond_list:
                    bond_list.append(duplicated_bonds[i][2][0][:3])
            if duplicated_bonds[i][j] > 1:
                bond_set_candidate = []
                minimum_distance = duplicated_bonds[i][2][0]
                for k in range(1,len(duplicated_bonds[i][2])):
                    if minimum_distance > duplicated_bonds[i][2][k]:
                        minimum_distance = duplicated_bonds[i][2][k]
                for k in range(len(duplicated_bonds[i][2])):
                    if minimum_distance == duplicated_bonds[i][2][k]:
                        bond_set_candidate.append(duplicated_bonds[i][2][k][:3])
                if len(bond_set_candidate) > 1:
                    shutil.move(file_dssp, file_with_flaw)
                    raise Exception ("we need further valiable to filter bonds")
                else:
                    if bond_set_candidate[0] not in bond_list:
                        bond_list.append(bond_set_candidate[0])
    return bond_list

def change_dssp_id_to_pdb_id(bond_list_dssp_id, data_pdb, dssp_id_with_residue_and_chain):
    bond_list_pdb_id = []
    for items in bond_list_dssp_id:
        append_list = [[],[]]
        for candidate in dssp_id_with_residue_and_chain:
            if items[0] == candidate[0]:
                append_list[0] = [candidate[1], "N ", candidate[2]]
            if items[1] == candidate[0]:
                append_list[1] = [candidate[1], "O ", candidate[2]]
        bond_list_pdb_id.append(append_list)
    return bond_list_pdb_id


def count_max_edge_length(connected_set):
    count = 0
    for item in connected_set:
        if item.count('_') > count:
            count = item.count('_')
    return count


def select_bond(model, target, flag):
    new_bond = []
    for bond in target:
        if model == 'chain' or 'whole':
            if bond[1] not in flag and bond[2] not in flag:
                new_bond.append(bond)
                flag.add(bond[1])
                flag.add(bond[2])
        elif model == 'n_v':
            new_bond.append(bond)
            flag.add(bond[1])
            flag.add(bond[2])
    return new_bond, flag


def dssp_id_to_pdb_id(data_d_p, bonds):
    bond_pdb = [[] for i in range(len(bonds))]
    for i, bond in enumerate(bonds):
        bond_pdb[i].extend(data_d_p[bond[1]])
        bond_pdb[i].extend(data_d_p[bond[2]])
    return bond_pdb


def get_h_bond(model, data, file):
    hydrogen_bonds = dict()
    threshold_list = configuration.threshold()
    flag = set()
    data_dssp = data['dssp']
    data_d_p = data['dssp_pdb']
    if len(data_dssp) > 0:
        if len(data_dssp[0]) > 0:
            start_potential = data_dssp[0][0]
            count = 0
            target_bond = [[]]
            for bond in data_dssp:
                if start_potential < bond[0]:
                    start_potential = bond[0]
                    target_bond.append([])
                    count += 1
                if start_potential == bond[0]:
                    if len(data_d_p[bond[1]]) > 1 and len(data_d_p[bond[1]]) > 1:
                        if data_d_p[bond[1]][1] == data_d_p[bond[2]][1]:
                            target_bond[count].append(bond)
            for target in target_bond:
                if target != []:
                    potential = target[0][0]
                    bonds, flag = select_bond(model, target, flag)
                    for threshold in threshold_list:
                        if float(threshold) >= potential:
                            bond_pdb = dssp_id_to_pdb_id(data_d_p, bonds)
                            hydrogen_bonds.setdefault(threshold, []).extend(bond_pdb)
        else:
            shutil.move(file.path[file.dssp], '/Users/harunegami/Desktop/'+file.pdbID+'_error.dssp')
    else:
        shutil.move(file.path[file.dssp], '/Users/harunegami/Desktop/'+file.pdbID+'_error.dssp')
    for threshold in threshold_list:
        if threshold in hydrogen_bonds.keys():
            operate_file.print_output(hydrogen_bonds[threshold], file.path[file.tau][str(threshold)], file)
    return hydrogen_bonds


def get_hydrogen_bond_list_filtered_by_potential_t(bond_list, atom_dict):
    bond_list_selected_O = select_stronger_bond(atom_dict["O"], bond_list, 1)
    bond_list_selected_NH = select_stronger_bond(atom_dict["NH"], bond_list, 0)
    bond_list_selected = bond_list_selected_O + bond_list_selected_NH
    unique_bond_list = unique(bond_list_selected)
    max_bond_id = get_max_bond_id(bond_list_selected)
    duplicated_bonds = get_duplicated_bonds(max_bond_id, unique_bond_list)
    duplicated_bonds, unique_bond_list = check_uniqueness_of_bonds_with_lowest_potential(duplicated_bonds, unique_bond_list)
    return duplicated_bonds, unique_bond_list, max_bond_id
