# coding: utf-8

import numpy
import error_check


def cal_base(residue_list, file):
    orthonormal_basis = []
    basis = []
    chain = residue_list[0][2]
    count = 0
    for i in range(1, len(residue_list)):
        residue_1 = residue_list[i-1]
        residue_2 = residue_list[i]
        chain_1 = residue_1[2]
        chain_2 = residue_2[2]
        if chain_1 == chain_2:
            if chain_1!= chain:
                count += 1
                chain = chain_1
            place_vec = extract_vector(residue_1, residue_2, i, file)
            vector_on_a_plane = calculate_axes(place_vec)
            basis.append(vector_on_a_plane)
            flame = normalize_basis(vector_on_a_plane)
            flame.extend(residue_1[:4])
            flame.append(count)
            flame.extend(residue_2[:4])
            flame.append(count+1)
            orthonormal_basis.append(flame)
            count += 1
    return orthonormal_basis, basis

def extract_vector(residue_1, residue_2, i, file):
    check = ["N ", "CA", "C ","O "]
    vector = [[] for i in range(5)]
    for i in range(len(residue_2)):
        for j in range(2):
            if residue_2[i] == check[j]:
                vector[j] = residue_2[i+1:i+4]
    for i in range(len(residue_1)):
        for j in range(2,5):
            if residue_1[i] == check[j-1]:
                vector[j] = residue_1[i+1:i+4]
    return vector


def calculate_axes(vector):
    [N_2, CA_2, CA_1, C_1, O_1] = vector
    xi = [x - y for (x, y) in zip(N_2, C_1)]
    yi = [x - y for (x, y) in zip(C_1, CA_1)]
    zi = [x - y for (x, y) in zip(CA_2, N_2)]
    return [xi,yi,zi]

def normalize_basis(basis):
    lengU = numpy.linalg.norm(basis[0])
    ui = [x/lengU for x in basis[0]]
    ui_dot_vi = numpy.dot(ui,basis[1])
    vi_to_be_normalized = [y - (ui_dot_vi)*u for (y, u) in zip(basis[1], ui)]
    lengV = numpy.linalg.norm(vi_to_be_normalized)
    vi = [v/lengV for v in vi_to_be_normalized]
    wi = list(numpy.cross(ui, vi))
    return [ui,vi,wi]

