# coding: utf-8

import os
import shutil
import cal_fatgraph.judge as judge
import cal_fatgraph.hydrogen as hydrogen
import cal_fatgraph.generate_fatgraph as fatgraph
import cal_fatgraph.track_boundary as boundary
import configuration.classFilePath as FilePath
import configuration
import error_check
import operate_file.operate_file as operate_file


def get_coordinate_data(model, file):
    data = operate_file.get_data(file, file.path[file.dssp])
    h_bonds = []
    if error_check.lack_of_data(data, file, [1, 1], 'atom error') == 1:
        h_bonds = hydrogen.get_h_bond(model, data, file)
        return data, h_bonds
    else:
        #shutil.move(file.path[file.dssp], file.path[file.flawed])
        return 0, 0


def get_activesite(file, input_file):
    data = operate_file.get_data(file, input_file)
    return data


def make_activesite(input_file):
    file, input_file = extract_file_paths(input_file)
    data = get_activesite(file, input_file)
    if data != 0:
        check_active_site(data['activesite'], file)


def pdb_loop(model, input_file):
    file = FilePath.FileClass(model, input_file)
    data, h_bonds = get_coordinate_data(model, file)
    if error_check.lack_of_data(data, file, 1, 'error file') == 1:
        result = calculate_invariants(data, h_bonds, file, model)
        invariant, length_of_invariants, error_type = result
        if error_check.lack_of_data(invariant, file, 1, error_type) == 1:
            if model != 'whole':
                with open(file.path[file.invhead], 'w') as f:
                    for item in invariant[1]:
                        for lengths in item[1:]:
                            for i, length in enumerate(lengths):
                                str_chain = data['chain'][i]
                                str_inv = [str(n) for n in length]
                                inv = (' ').join(str_inv)
                                str_list = str_chain + ':' + inv
                                f.write(str_list)
                                f.write('\n')
            else:
                with open(file.path[file.invhead], 'w') as f:
                    for item in invariant[1]:
                        for lengths in item[1:]:
                            str_inv = [str(n) for n in lengths[0]]
                            inv = (' ').join(str_inv)
                            f.write(inv)
                            f.write('\n')
            with open(file.path[file.invariants], 'w') as f:
                for item in length_of_invariants:
                    str_list = [str(n) for n in item[1:]]
                    f.write(' '.join(str_list))
                    f.write('\n')
            write_result(invariant[0], length_of_invariants, file, data)


def write_result(invariant, length_of_invariants, file, data):
    thresholds = configuration.configuration.threshold()
    boundary = list(zip(thresholds, invariant))
    for i, threshold in enumerate(thresholds):
        path = file.path[file.topology][threshold]
        operate_file.print_output(boundary[i], path, file)


def calculate_invariants(data, h_bonds, file, model):
    thresholds = configuration.configuration.threshold()
    boundaries = [[data["header"]] for i in range(len(thresholds))]
    num_boundary = [[data["header"]] for i in range(len(thresholds))]
    invariants = []
    #check_active_site(data["activesite"], file)
    for i, threshold in enumerate(thresholds):
        if threshold in h_bonds.keys():
            result = judge.rotations(data, h_bonds[threshold], file, threshold)
        else:
            result = judge.rotations(data, [], file, threshold)
        [twisted_data, error_type] = result
        if error_check.lack_of_data(twisted_data, file, 1, error_type) == 1:
            fat = fatgraph.get_fatgraph(data, twisted_data, file, model)
            if fat[0] == 0:
                return 0, 0, fat[1]
            bc = boundary.count(fat, data, file, threshold, boundaries, model)
            invariant, length_of_invariants, error_type = bc
            if error_check.lack_of_data(invariant, file, 1, 'track') == 1:
                invariants.append(invariant[0])
                num_boundary[i].append(length_of_invariants)
                boundaries[i].append(invariant[1])
                error_type = 'none'
                invariants_list = [invariants, boundaries]
            else:
                return 0, 0, 'track boundary'
    return invariants_list, num_boundary, error_type


def check_active_site(active_sites, file):
    site = dict()
    print(file.pdbID)
    for key, values in active_sites.items():
        residue = []
        chain = []
        residue_num = []
        query = []
        for val in values:
            print(val)
            for i in range(4):
                if val[0+11*i:3+11*i] != '   ':
                    residue.append(val[0+11*i:3+11*i])
                    chain.append(val[4+11*i])
                    residue_num.append(str(int(val[5+11*i:9+11*i])))
        for i in range(len(residue)):
            query.extend([residue[i], chain[i], residue_num[i]])
        site.setdefault(key,[]).append(('_').join(query))
    path = '/Users/harunegami/Desktop/activesite/'+file.pdbID+'.txt'
    with open(path, 'w') as f_site:
        for key, item in site.items():
            str_value = ('=').join(item)
            f_site.write(key+':'+ str_value+'\n')
    print(site)





