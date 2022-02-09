# coding: utf-8

import sys
import os


class Write_text_previous(object) :
    def __init__(self, data_with_format,output_file) :
        pass
    # テキストを読み込む
    def write_text(self,output_data,data_with_format,output_file) :
        with open (output_file, 'a') as output:
            for list in data_with_format:
                output.writelines(list)
                output.write("\n")

class Write_text(object) :
    def __init__(self, data_with_format,output_file) :
        pass
    # テキストを読み込む
    def write_text(self,output_data,data_with_format,output_file) :
        with open(output_file, 'w') as f:
            f.write( "\n".join(data_with_format) )

def format_of_the_file(output_data, process_type, file):
    #if process_type == file.tau:
    #    format ='{0[0]:>4} {0[1]:>6} {0[2]:>7} {1[0]:>5} {1[1]:>6} {1[2]:>7}'.format(output_data[0],output_data[1])
    #    return format
    #if process_type == file.twisted:
    #    format ='{0[0]:>4} {0[1]:>8} {0[2]:>7} {0[3]:>6} {0[4]:>5} {0[5]:>8} {0[6]:>7} {0[7]:>6} {0[8]:>2}'.format(output_data)
    #    return format
    if process_type == file.topology:
        format = str(output_data)
        return format
    elif process_type == file.invariants:
        format = str(output_data)
    else:
        format = str(output_data)
        return format

class Write_to_file(Write_text):
    def __init__(self, output_data, output_file):
        super().__init__(output_data, output_file)
    def write_text(self,output_data,data_with_format, output_file, file):
        print(output_file+": done")
        process_type = output_file.split(".")
        data_with_format = []
        if process_type[-1] == file.tau:
            data_with_format.append('RES1  ATOM1  CHAIN1  RES2  ATOM2  CHAIN2')
        elif process_type[-1] == file.twisted:
            data_with_format.append('RES1  NUMRES1  CHAIN1  ATOM1  RES2  NUMRES2  CHAIN2  ATOM2  t/u')
        for data in output_data:
            data_with_format.append(format_of_the_file(data, process_type[-1], file))
        super().write_text(self,data_with_format,output_file)
