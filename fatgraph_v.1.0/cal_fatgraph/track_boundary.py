# coding: utf-8

import math
import error_check
import shutil
import configuration.configuration as configuration


def count(conversion_function_list, extracted_data, file, threshold, length_of_invariants_for_each_threshold, model):
    N = configuration.get_number_of_node(model)
    fatgraph_dict, M, number_of_residue = conversion_function_list
    boundary_component = [[]]
    half_edge = fatgraph_dict['half_edge']
    concat_edge = fatgraph_dict['concat_edge']
    cyclic_order = fatgraph_dict['cyclic_order']
    flag = set([i for i in range(len(half_edge))])
    count = 0
    start = half_edge[0]
    length_of_boundary_component = [[] for i in range(len(extracted_data['chain']))]
    number_of_boundary_component = 0
    start_index_for_each_chain = [0 for i in range(len(number_of_residue))]
    for i in range(1, len(number_of_residue)):
        start_index_for_each_chain[i] = start_index_for_each_chain[i-1] + number_of_residue[i-1]
    chain = set()
    for i in range(len(half_edge)):
        halfedge_id = str(half_edge[count])
        concat_id = str(concat_edge[count])
        cycle_id = str(cyclic_order[count])
        chain.add(halfedge_id[0])
        chain.add(concat_id[0])
        chain.add(cycle_id[0])
        boundary_component[number_of_boundary_component].append(half_edge[count])
        boundary_component[number_of_boundary_component].append(concat_edge[count])
        boundary_component[number_of_boundary_component].append(cyclic_order[count])
        flag.remove(count)
        if start == cyclic_order[count] and i < len(half_edge)-1:
            length = len(set(boundary_component[number_of_boundary_component]))
            if type(length) == list:
                length.sort()
            if len(chain) > 1:
                print('**************************************************error')
                shutil.move(file.path[file.dssp],'/Users/harunegami/Desktop/ligand/'+file.pdbID+'.dssp')
                return 0,0, 'track boundary'
            else:
                chain_index = int(chain.pop())
            length_of_boundary_component[chain_index-1].append(length)
            number_of_boundary_component += 1
            boundary_component.append([])
            count = flag.pop()
            flag.add(count)
            start = half_edge[count]
            chain = set()
        else:
            count = find_index(cyclic_order[count], M, int(len(half_edge)/2),start_index_for_each_chain, N)
    length = len(set(boundary_component[number_of_boundary_component]))
    if type(length) == list:
        length.sort()
    if len(chain) > 1:
        print('**************************************************error')
        shutil.move(file.path[file.dssp],'/Users/harunegami/Desktop/ligand/'+file.pdbID+'.dssp')
        return 0,0, 'track boundary'
    else:
        chain_index = int(chain.pop())
    length_of_boundary_component[chain_index-1].append(length)
    number_of_boundary_component += 1
    print(length_of_boundary_component)
    return [boundary_component,length_of_boundary_component], number_of_boundary_component,'OK'


def find_index(value, M, num,start_index_for_each_chain, N):
    chain, updown_halfedge = divmod(value, 10*M)
    updown, residue_place = divmod(updown_halfedge, M)
    residue, place = divmod(residue_place, 2 * N)
    residue_with_chain = start_index_for_each_chain[chain - 1] + residue
    next_index = residue_with_chain * 2 * N + place + updown * num
    return next_index


def print_invariants(file, threshold, conversion_function, tau_with_chain, boundary_components, M, extracted_data):
    invariant = []
    pdbID = file.pdbID
    edge = conversion_function["edge"]
    invariant.append('******* ' + pdbID + ' *******')
    invariant.append('***** ' + threshold + ' *****')
    invariant.append("# of halfedge: " + str(len(edge)))
    invariant.append("# of tau_hydrogen_bond_with_chain_info: "+ str(len(tau_with_chain["hydrogen"])))
    invariant.append("# of edge: " + str(len(tau_with_chain["peptide"])))
    for i, chain in enumerate(extracted_data["chain"]):
        invariant.append("# of boundary component "+ chain + " : " + str(len(boundary_components[i])))
    for items in invariant:
        print(items)
    return invariant
