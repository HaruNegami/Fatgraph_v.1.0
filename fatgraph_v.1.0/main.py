# coding: utf-8

import argparse
from importlib import import_module
from multiprocessing import Pool
import multiprocessing as multi
import os
import configuration.classFilePath as FilePath
import pdb_loop

"""
Usage:
    python main.py model_option test_option update_option

Usage example:
    python main.py chain test difference

model_option:
    chain : calculate fatgraph invariants for each chain
    whole : calculate fatgraph invariants of whole structure
    n_v   : calculate fatgraph with n half-edges for each node

test_option:
    test       : calculate fatgraph invariants of file at ./test_dssp/
    production : calculate fatgraph invariants of data at ./dssp/

update_option:
    overwrite  : calculate all files
    difference : calculate iff not yet calculated

"""


def wrapper(args):
    return pdb_loop.pdb_loop(*args)


def multi_process(model, file_paths):
    files = [[] for i in range(28)]
    for file in file_paths:
        size = os.path.getsize(file)
        for i in range(28):
            if 100000 * i <= size and size <= 100000 * (i + 1):
                files[i].append((model, file))
    for i in range(28):
        if i < 15:
            multiProcess = Pool(multi.cpu_count())
            multiProcess.map(wrapper, files[i])
            multiProcess.close()
        else:
            multiProcess = Pool(multi.cpu_count()-4)
            multiProcess.map(wrapper, files[i])
            multiProcess.close()




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate invarinats')
    parser.add_argument('model_option')
    parser.add_argument('test_option')
    parser.add_argument('update_option')

    args = parser.parse_args()
    model_option = args.model_option
    test_option = args.test_option
    update_option = args.update_option
    file_paths = FilePath.get_paths(model_option, test_option, update_option)
    if test_option == 'test':
        for file in file_paths:
            pdb_loop.pdb_loop(model_option, file)
    else:
        multi_process(model_option, file_paths)
