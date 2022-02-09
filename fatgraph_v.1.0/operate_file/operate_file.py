# coding: utf-8

import sys
import os
import glob
from operate_file import ReadFile
from operate_file import PrintFile
import copy
import shutil

####################################################
### read/ write data ###############################
####################################################


def get_data_from_dssp(file_dssp):
    protein_from_dssp = ReadFile.dssp_FileRead(file_dssp)
    data_dssp = protein_from_dssp.read_text(file_dssp)
    extracted_data_dssp, d_p_dict = protein_from_dssp.read_hydrogenBond(data_dssp)
    extracted_data_dssp.sort(key=lambda x: x[0]) #sort by potential
    return extracted_data_dssp, d_p_dict

def get_data_from_pdb(file):
    protein_from_pdb = ReadFile.pdb_FileRead(file.path[file.pdb])
    data_pdb, chain_whole, header, active_site = protein_from_pdb.read_atom_text(file.path[file.pdb])
    if data_pdb == 0:
        #shutil.move(file.path[file.dssp], '/Users/harunegami/Desktop/ligand/'+file.pdbID+'.dssp')
        return 0, 0, 0, 0
    extracted_data_pdb = protein_from_pdb.read_coordinates(data_pdb, file.path[file.pdb],file)
    return extracted_data_pdb, chain_whole, header, active_site

def get_data_from_tauset(file):
    tau = file.tau
    file_dict = file.path
    protein_from_tauset = ReadFile.tauset_FileRead(file_dict[tau])
    data_HB = protein_from_tauset.read_text(file_dict[tau])
    extracted_data_tauset = protein_from_tauset.read_tauEdge(data_HB)
    return extracted_data_tauset


def get_data(file, input_file):
    pdb = file.pdb
    file_dict = file.path
    extracted_data_dict = dict()
    extracted_data_dssp, d_p_dict =  get_data_from_dssp(input_file)
    extracted_data_pdb, chain_whole, header, active_site =  get_data_from_pdb(file)
    extracted_data_dict.setdefault("dssp", extracted_data_dssp)
    extracted_data_dict.setdefault("dssp_pdb", d_p_dict)
    extracted_data_dict.setdefault("pdb", extracted_data_pdb)
    extracted_data_dict.setdefault("chain", chain_whole)
    extracted_data_dict.setdefault("header", header)
    extracted_data_dict.setdefault("activesite", active_site)
    return extracted_data_dict

def get_data_chain_only(file):
    pdb = file.pdb
    file_dict = file.path
    chain_whole =  get_data_from_pdb(file_dict[pdb])
    return chain_whole

def print_output(output_data, output_file, file):
    data_with_format = []
    printData = PrintFile.Write_to_file(output_data,output_file)
    printData.write_text(output_data,data_with_format,output_file, file)
