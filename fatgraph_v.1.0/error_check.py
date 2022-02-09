# coding: utf-8

import shutil
import os
import configuration.classFilePath as classFilePath


def error_categorize(input_file):
    error = []
    file, input_file = extract_file_paths(input_file)

    if os.path.exists(file.path[file.error]):
        with open(file.path[file.error], 'r') as f:
            for line in f:
                error.append(line.rstrip())
    return error


def extract_file_paths(input_file):
    if "test" in input_file:
        len_path = len(input_file)
        input_file = input_file[: len_path - 4]
        file = test_classFilePath.FileClass(input_file)
    else:
        file = classFilePath.FileClass(input_file)
        print(file.path[file.pdb])
    return file, input_file


def categorize(conversion_function):
    if "peptide" not in conversion_function.keys():
        return 0, 0, "no keys of peptide"
    elif "sigma" not in conversion_function.keys():
        return 0, 0, "no keys of sigma"
    elif conversion_function["peptide"] == "error":
        return 0, 0, "peptide error"
    elif conversion_function["sigma"] == 0:
        return 0, 0, "sigma error"
    else:
        return 1, 0, "no error"


def write_error(file, detail):
    with open(file.path[file.error], "a") as f_error:
        print(detail)
        f_error.write(str(file.pdbID) + ': ERROR ' + str(detail[0]) + "\n")
        print(str(file.pdbID) + ': ERROR ' + str(detail[0]) + "\n")


def check_error(requirement, arr, file, detail):
    if type(arr) == 0:
        write_error(file, detail)
        return 0
    elif type(arr) == list:
        if len(arr) < requirement:
            write_error(file, detail)
            return 0
    return 1


def lack_of_data(data, file, requirement_length, *detail):
    check_next = 1
    if type(requirement_length) == int:
        check_next = check_error(requirement_length, data, file, detail)
    elif type(requirement_length) == list:
        for requirement in requirement_length:
            check_next = check_error(requirement, data, file, detail)
            if check_next == 1:
                if type(data) == list:
                    data = data[0]
                elif type(data) == dict:
                    if file.pdb in data.keys():
                        data = data['pdb']
                    else:
                        print('uncategorized error')
            else:
                print(data, file.pdbID)
    return check_next


def list_difference(list1, list2):
    for item1, item2 in zip(list1, list2):
        if item1 != item2:
            print(item1, item2)


def test_same_result(previous, test, msg):
    if previous == test:
        print("\n-------*------- test OK -------*------- ", msg, "\n\n")
        return 1
    else:
        print("\n-------*------- test NG -------*------- ", msg, "\n\n")
        for items in list_difference(previous, test):
            print(items)
        return 0
