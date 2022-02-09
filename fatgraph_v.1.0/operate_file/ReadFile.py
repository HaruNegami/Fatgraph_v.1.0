# coding: utf-8

import sys
import os
from operate_file import abstractReadFile
import re
import error_check

def make_HydrogenBond(line, extracted_data, d_p_dict):
    dssp_ID = int(line[:5])
    pdb_ID = int(line[5:10])
    chain = line[11]
    d_p_dict.setdefault(dssp_ID, [pdb_ID, chain])
    axes_CA = [float(line[115:122]),float(line[122:129]),float(line[129:136])]
    N_1 = [float(line[46:50]), dssp_ID, dssp_ID+int(line[39:45]), chain, axes_CA]
    N_2 = [float(line[68:72]), dssp_ID, dssp_ID+int(line[61:67]), chain, axes_CA]
    O_1 = [float(line[57:61]), dssp_ID+int(line[50:56]), dssp_ID, chain, axes_CA]
    O_2 = [float(line[79:83]), dssp_ID+int(line[72:78]), dssp_ID, chain, axes_CA]
    if N_1[1] != N_1[2]:
        extracted_data.append(N_1)
    if N_2[1] != N_2[2]:
        extracted_data.append(N_2)
    if O_1[1] != O_1[2]:
        extracted_data.append(O_1)
    if O_2[1] != O_2[2]:
        extracted_data.append(O_2)
    return extracted_data, d_p_dict

def make_TauSet(line):
    atom1 = (line[19])
    residue1 = int(line[9:13])
    chain1 = (line[29])
    atom2 = (line[54])
    residue2 =int(line[43:47])
    chain2 = (line[65])
    tempLine=["0",residue1,chain1,atom1,"0",residue2,chain2,atom2]
    return(tempLine)

def input_coordinates(line, extracted_data, seen_keys,countResidue, next):
    if next == 0:
        pass
    else :
        atomNum = int(line[6:11])
        tempLine = []
        atom = line[13:15]
        residue = line[17:20]
        numResidue = int(line[22:26])
        chain = line[21]
        key = numResidue, chain
        placeX = float(line[30:38])
        placeY = float(line[38:46])
        placeZ = float(line[46:54])
        if key not in seen_keys:
            seen_keys.add(key)
            extracted_data.append(
                [residue,numResidue,chain,atom,placeX,placeY,placeZ])
            countResidue += 1
        else:
            extracted_data[-1].extend(
                [residue,numResidue,chain,atom,placeX,placeY,placeZ])
    return extracted_data,seen_keys

def make_twistedEdge(line, extracted_data_HB, chainNuM):
    residue1 = line[9:12]
    numResidue1 = int(line[26:32])
    chain1 = line[39]
    atom1 = line[50:52]
    residue2 = line[9+55:12+55]
    numResidue2 = int(line[26+55:32+55])
    chain2 = line[39+55]
    atom2 = line[50+55:52+55]
    twist = line[52+55+10]
    extracted_data_HB.append([residue1,numResidue1,chain1,atom1,residue2,numResidue2,chain2,atom2,twist])
    if chain1 not in chainNuM:
        chainNuM.append(chain1)
    elif chain2 not in chainNuM:
        chainNuM.append(chain2)
    return(extracted_data_HB, chainNuM)



class dssp_FileRead(abstractReadFile.Read_text):
    def __init__(self,input_file):
        super().__init__(input_file)

    def read_text(self, input_file):
        print(input_file)
        raw_data = super().read_text(input_file)
        return raw_data

    def read_hydrogenBond(self,raw_data):
        extracted_data = []
        d_p_dict = dict()
        signal = 0
        for i in range(len(raw_data)):
            line = raw_data[i]
            if signal > 0:
                if line[13]!= "!":  # in the chain
                    extracted_data, d_p_dict = make_HydrogenBond(line, extracted_data, d_p_dict)
                    signal += 1
            if "#  RESIDUE AA" in raw_data[i]: # read data below this row
                signal = 1
        return extracted_data, d_p_dict


class pdb_FileRead(abstractReadFile.Read_pdb_text):
    def __init__(self, input_file):
        super().__init__(input_file)

    def read_atom_text(self, input_file):
        raw_data, chain_whole, header, active_site = super().read_atom_text(input_file)
        return raw_data, chain_whole, header, active_site
    def read_coordinates(self, raw_data, input_file, file):
        extracted_data = []
        countResidue = 0
        next_atom = 1
        seen_keys = set()
        if raw_data == 0:
            return 0
        for i in range(len(raw_data)):
            line = raw_data[i]
            if line[16] == " ":
                next_atom = 1
            elif line[16]!= " " and i+1 < len(raw_data):
                nextAtom = raw_data[i+1]
                if ord(nextAtom[16])> ord(line[16]) and line[13:15] == nextAtom[13:15]:
                    next_atom = 1
                else:
                    next_atom = 0
            else:
                next_atom = 0
            extracted_data, seen_keys = input_coordinates(line, extracted_data,seen_keys,countResidue,next_atom)
        for items in extracted_data:
            if (error_check.lack_of_data(list(set(items)&set(["N ", "C ", "CA", "O "])), file, 4, "lack of atom 1"))==0:
                print('lack of atom here: ', file.pdbID, items)
                return 0
        return extracted_data


class FileRead(abstractReadFile.Read_text):
    def __init__(self,input_file):
        super().__init__(input_file)
    def read_text(self, input_file):
        raw_data = super().read_text(input_file)
        return raw_data
    def read_tauEdge(self, raw_data):
        extracted_data_HB = []
        self.chainNuM = []
        for i in range(len(raw_data)):
            line = raw_data[i]
            extracted_data_HB.append(make_TauSet(line))
        return extracted_data_HB
