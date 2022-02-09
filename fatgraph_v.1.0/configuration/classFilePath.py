# coding: utf-8

import sys
import os
import glob
import copy
import shutil
import configuration.configuration as configuration


####################################################
### settings #######################################
### please edit the file_path or extensions here ###
####################################################

## -*- absolute path of the dir of the files -*- ##

ABS_PATH = '/Users/harunegami/Desktop/'
TEST_ABS_PATH = ABS_PATH + "/test/"
## -*- absolute path of the error log files -*- ##
ABS_PATH_ERROR = ABS_PATH + "/error_"
TEST_ABS_PATH_ERROR = ABS_PATH + "/test/error_"

####################################################
## -**- extensions of the files              -**- ##

# -*-  dir for flawed data                     -*- #

flawed = 'flawed'

# -*-  extension for twisted/ untwisted data  -*- #

twisted = 'twisted'

# -*-  extension for twisted/ untwisted data  -*- #

tau = 'tau'

# -*-  extension for boundary components      -*- #

topology = 'topology'

# -*-  extension for the summary of the invariants#

invariants = 'invariants'

# -*-  extension for the summary of the length_of boundary components with header#

invhead = 'invhead'

# -*-  extension for finished files           -*- #

done_dssp = 'done_dssp'

# -*-  extension for error files           -*- #

error = 'error'

####################################################
### end of the settings ############################
### please don't edit below  #######################
####################################################

# extension for dssp file
dssp = 'dssp'

# extensin for pdb file
pdb = 'pdb'


extensions_integrated = [dssp, pdb, invariants, flawed, invhead, done_dssp, error]
extensions_for_each_threshold = [twisted, tau, topology]
#extensions_for_each_seq = [invhead]


def make_path(arr):
    path = ('/').join(arr)
    return path


def get_paths(model, test_option, update_option):
    if test_option == 'test':
        path = ABS_PATH + 'test/*.dssp'
    else:
        path = ABS_PATH + dssp + '/*/*.dssp'
    filenames = glob.glob(path)
    filenames_new = []
    for filename in filenames:
        file = filename.split("/")
        if len(file[-1]) == 9:
            if update_option == 'difference':
                pdb = file[-1].replace('.dssp', '')
                dirs = [ABS_PATH, model, invhead, pdb[1:3], pdb, '.invhead']
                if os.path.exists(make_path(dirs)) is False:
                    filenames_new.append(filename)
            else:
                return filenames
    return filenames_new


def get_error_proteins():
    path = ABS_PATH_ERROR + '*.txt'
    filenames = glob.glob(path)
    error_pdbID = []
    for filename in filenames:
        with open(filename, "r") as error_f:
            for line in error_f:
                if len(line) > 4:
                    error_pdbID.append(line[:4])
    return error_pdbID


def pdbID_from_input_file(abs_path):
    dir_list = abs_path.split('/')
    pdbID_with_extension = dir_list[-1]
    pdbID_extension = pdbID_with_extension.split('.')
    pdbID = pdbID_extension[0]
    return pdbID


def subfoldername_from_pdbID(pdbID):
    foldername = pdbID[1:3]
    return foldername


def make_file_paths(model, input_file):

    file_dict = dict()
    pdbID = pdbID_from_input_file(input_file)
    file_dict.setdefault("pdbID", pdbID)
    folder = subfoldername_from_pdbID(pdbID).lower()
    threshold_list = configuration.threshold()

    for extension in extensions_integrated:
        if extension == 'dssp':
            file_dict.setdefault(extension, input_file)
        else:
            if extension == 'pdb':
                path = '/Users/harunegami/Documents/pdb/' + folder
            else:
                path = ABS_PATH + model + '/' + extension + '/' + folder
            os.makedirs(path, exist_ok=True)
            file = path + '/' + pdbID + '.' + extension
            file_dict.setdefault(extension, file)
    for extension in extensions_for_each_threshold:
        for TH_str in threshold_list:
            path_list = [ABS_PATH, model, extension, folder, TH_str]
            path = make_path(path_list)
            os.makedirs(path, exist_ok=True)
            file = path + '/' + pdbID + '.' + extension
            TH_PATH = {TH_str: file}
            file_dict.setdefault(extension, TH_PATH).update(TH_PATH)
    return file_dict


class FileClass:
    def __init__(self, model, input_file):
        self.pdbID = pdbID_from_input_file(input_file)
        self.path = make_file_paths(model, input_file)
        self.dssp = dssp
        self.pdb = pdb
        self.invariants = invariants
        self.flawed = flawed
        self.tau = tau
        self.twisted = twisted
        self.topology = topology
        self.invhead = invhead
        self.done_dssp = done_dssp
        self.error = error

