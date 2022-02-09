# coding: utf-8

import math
import error_check
import configuration.configuration as configuration


def get_fatgraph(extracted_data, twisted_data_dict, file, model):
    fatgraph_dict = dict()
    N = configuration.get_number_of_node(model)
    M, sigma_total, number_of_residue = generate_half_edge(extracted_data, file, N)
    fatgraph_list = generate_fatgraph(extracted_data,number_of_residue, twisted_data_dict, sigma_total,N)
    if fatgraph_list == 'generate fatgraph error':
        return 0, fatgraph_list
    fatgraph_dict.setdefault('half_edge', fatgraph_list[0])
    fatgraph_dict.setdefault('concat_edge', fatgraph_list[1])
    fatgraph_dict.setdefault('cyclic_order', fatgraph_list[2])
    return [fatgraph_dict, M, number_of_residue]


def generate_half_edge(extracted_data, file, N):
    chain_whole = extracted_data["chain"]
    number_of_residue = count_residue(extracted_data)
    half_edge_total, M = make_sigma_list(chain_whole, number_of_residue, N)
    return M, half_edge_total, number_of_residue


def generate_fatgraph(extracted_data,number_of_residue, twisted_data, sigma_total, N):
    peptide_bond = twisted_data["peptide"]
    hydrogen_bond = twisted_data["hydrogen"]
    fatgraph = extend_sigma_total(sigma_total, N)
    up = int(len(fatgraph[0])/2)
    count = 0
    start_residue = 0
    length = len(number_of_residue)
    for chain in range(length):
        start_residue += number_of_residue[chain]
        for i in range(number_of_residue[chain]):
            # amino
            if i == 0 or i == number_of_residue[chain]-1:
                residue_0 = start_residue - number_of_residue[chain]
                residue_end = start_residue - 1
                if i == 0:
                    end = [residue_0, 'u', 0]
                else:
                    end = [residue_end, 'u', 0]
                    fatgraph = concatenate_half_edge(-(N-1),0,fatgraph,end, up,count, N)
                    count += 1
                fatgraph = concatenate_half_edge(N-1,N,fatgraph,end, up,count, N)
            else:
                if count >= len(peptide_bond):
                    return 'generate fatgraph error'
                peptide_type = peptide_bond[count][-1]
                if peptide_type == 1:
                    fatgraph = change_res_direction(peptide_bond[count][-3], up, fatgraph, N)
                fatgraph = concatenate_half_edge(N+1,2*N,fatgraph,peptide_bond[count], up,count, N)
                fatgraph = concatenate_half_edge(3*N-1,3*N,fatgraph,peptide_bond[count], up,count, N)
                count += 1
    for h_bond in hydrogen_bond:
        fatgraph = concatenate_half_edge_h(1, 5+(N-3), fatgraph, h_bond, up,count, N)
    return fatgraph


def count_residue(extracted_data):
    extracted_data_pdb = extracted_data["pdb"]
    chain_whole = extracted_data["chain"]
    number_of_residue_for_each_residue = [0 for i in range(len(chain_whole))]
    count = 0
    for items in extracted_data_pdb:
        count = chain_whole.index(items[2])
        number_of_residue_for_each_residue[count] += 1
    return number_of_residue_for_each_residue

def make_sigma_list(chain_whole, number_of_residue, N):
    half_edge_up = []
    half_edge_down = []
    sigma_inverse_list = []
    order = int(math.log10(sum(number_of_residue) * 2*N) + 1) + 1
    M = 10 ** order
    for i in range(len(chain_whole)):
        chainINDEX = (i+1)*10*M
        for j in range(number_of_residue[i]):
            for k in range(2*N):
                half_edge_up.append(chainINDEX+ 2*N*j+k)
                half_edge_down.append(chainINDEX+ 2*N*j+k+M)
    half_edge_total = half_edge_up + half_edge_down
    return half_edge_total, M


def extend_sigma_total(sigma_total, N):
    edge = []
    tau_edge = []
    sigma_edge = []
    for i in range(len(sigma_total)):
        if i < int(len(sigma_total)/2):
            res, place = divmod(i, N)
            edge.append(sigma_total[i])
            tau_edge.append(sigma_total[i])
            place2 = (place+1)%N
            sigma_edge.append(sigma_total[N*res+place2])
        else:
            res, place = divmod(i, N)
            edge.append(sigma_total[i])
            tau_edge.append(sigma_total[i])
            place2 = (place -1 +N) % N
            sigma_edge.append(sigma_total[N*res+place2])
    return [edge, tau_edge, sigma_edge]


def switch_half_edge(index_1, index_2, array):
    temp = array[index_2]
    array[index_2] = array[index_1]
    array[index_1] = temp
    return array


def get_start_index_for_each_chain(chain, number_of_residue):
    start_index = [0 for i in range(len(chain))]
    for j in range(1, len(chain)):
        start_index[j] = start_index[j-1] + number_of_residue[j-1]
    return start_index


def switch_two_sides(index1, index2, residue1, residue2, up_or_down, twist_type, fatgraph, N):
    place_1 = index1 + 2*N*residue1+(0+twist_type)*up_or_down
    place_2 = index2 + 2*N*residue2+0*up_or_down
    fatgraph[1] = switch_half_edge(place_1, place_2, fatgraph[1])
    fatgraph[2] = switch_half_edge(place_1, place_2, fatgraph[2])
    place_1 = index1 + 2*N*residue1+(1-twist_type)*up_or_down
    place_2 = index2 + 2*N*residue2+1*up_or_down
    fatgraph[1] = switch_half_edge(place_1, place_2, fatgraph[1])
    fatgraph[2] = switch_half_edge(place_1, place_2, fatgraph[2])
    return fatgraph


def concatenate_half_edge(index1, index2,fatgraph,twist_data, up_or_down,count, N):
    if index1 == (N-1) or index1 == -(N-1):
        twist_type = 0
        residue1 = twist_data[0]
    elif twist_data[-2] == 't':
        twist_type = 1
        residue1= twist_data[-3] - 1
    else:
        twist_type = 0
        residue1 = twist_data[-3] - 1
    fatgraph = switch_two_sides(index1, index2, residue1, residue1, up_or_down, twist_type, fatgraph, N)
    return fatgraph

def concatenate_half_edge_h(index1, index2,fatgraph,twist_data, up_or_down,count, N):
    if twist_data[-2] == 't':
        twist_type = 1
    else:
        twist_type = 0
    residue1 = twist_data[4]
    residue2 = twist_data[-3]
    fatgraph = switch_two_sides(index1, index2, residue1, residue2, up_or_down, twist_type, fatgraph, N)
    return fatgraph


def change_res_direction(res, up, fatgraph, N):
    start_res = res * 2 * N
    temp = fatgraph[2][start_res + 0]
    temp_down = fatgraph[2][start_res + N - 1 + up]
    for i in range(N-1):
        place_next = i + 1
        fatgraph[2][start_res + i] = fatgraph[2][start_res + place_next]
    for i in reversed(range(1, N)):
        next2 = i - 1
        fatgraph[2][start_res + i + up] = fatgraph[2][start_res + next2 + up]
    fatgraph[2][start_res + N - 1] = temp
    fatgraph[2][start_res + 0 + up] = temp_down
    return fatgraph

