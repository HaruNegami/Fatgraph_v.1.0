# coding: utf-8


def threshold():
    threshold_list = ["-3.5", "-3.0", "-2.5", "-2.0", "-1.5", "-1.0", "-0.5"]
    return threshold_list


def get_number_of_node(model):
    if model == 'n_v':
        N = 4
    else:
        N = 3
    return N
