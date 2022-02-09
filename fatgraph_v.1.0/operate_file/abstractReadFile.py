# coding: utf-8

import sys
import os
import datetime
import time
import re
import shutil

def make_Chain(line, chainWhole):
    pre = line.split(":")
    pre2 = pre[1].split(";")
    pre3 = pre2[0].split(",")
    for i in range(len(pre3)):
        for j in range(len(pre3[i])):
            if re.match(r"[A-Z]+",pre3[i][j]):
                chainWhole.append(pre3[i][j])
    return chainWhole


class Read_text(object) :
    def __init__(self, input_file) :
        pass
    # テキストを読み込む
    def read_text(self, input_filename) :
        raw_data = []
        with open (input_filename, 'r') as input:
            for line in input:
                line = line.strip("\n")
                line = re.sub('\t', '    ', line)
                raw_data.append(line)
        return(raw_data)

class Read_pdb_text(object) :
    def __init__(self, input_file) :
        pass
    # テキストを読み込む
    def read_text(self, input_filename) :
        raw_data = []
        with open (input_filename, 'r') as input:
            for line in input:
                line = line.strip("\n")
                line = re.sub('\t', '    ', line)
                raw_data.append(line)
        return(raw_data)

    def read_atom_text(self, input_filename) :
        raw_data = []
        chain = []
        header = ""
        active_site = dict()
        with open (input_filename, 'r') as input:
            for line in input:
                line = line.strip("\n")
                line = re.sub('\t', '    ', line)
                count = len('ANISOU   29  N  A')
                if line[:4] == "ATOM" or ('HETATM' == line[:6] and 'HOH' not in line[count:count+4]):
                    if line[13:15] in ["N ", "CA", "C ", "O "]:
                        raw_data.append(line)
                        if line[21] not in chain:
                            chain.append(line[21])
                elif line[:6] == "HEADER":
                    line = line.strip("\n")
                    line = re.sub('\t', '    ', line)
                    line = line[10:50]
                    header = line
                elif line[:4] == 'SITE':
                    residue = line[18:61]
                    active_site.setdefault(line[11:14],[]).append(residue)
                #if 'HETNAM' == line[:6]:
                #    return 0, 0, 0, 'ligand'
        chain.sort()
        return raw_data, chain, header, active_site


    def read_chain_text(self, input_filename) :
        data = []
        with open (input_filename, 'r') as input:
            for line in input:
                line = line.strip("\n")
                line = re.sub('\t', '    ', line)
                if (line[11:16]=="CHAIN") and (line[:6] == "COMPND"):
                    data.append(line)
        chainWhole = []
        for line in data:
            pre = line.split(":")
            pre2 = pre[1].split(";")
            pre3 = pre2[0].split(",")
            for i in range(len(pre3)):
                for j in range(len(pre3[i])):
                    if re.match(r"[A-Z]+",pre3[i][j]):
                        chainWhole.append(pre3[i][j])
        chainWhole.sort()
        return chainWhole
