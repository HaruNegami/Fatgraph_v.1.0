# coding: utf-8

import sys
import os


class Write_text(object) :
    def __init__(self, data_with_format,output_file) :
        pass
    # テキストを読み込む
    def write_text(self,output_data,data_with_format,output_file) :
        with open (output_file, 'w') as output:
            for list in data_with_format:
                output.writelines(list)
                output.write("\n")
