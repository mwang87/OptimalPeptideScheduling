#!/usr/bin/python

import os

def get_header_mappings(header_str):
    header_mapping = {}
    header_splits = header_str.rstrip().split("\t")
    index_count = 0
    for header in header_splits:
        header_mapping[header] = index_count
        index_count += 1
    return header_mapping

#Parses a filename and returns 2 things
#first is the number of lines, and then a map to lists with the key being the column.
def parse_table_with_headers(filename):
    input_file = open(filename, "r")

    line_count = 0
    headers = []
    index_to_header_map = {}
    column_values = {}
    for line in input_file:
        line_count += 1
        if line_count == 1:
            headers = line.rstrip().split("\t")
            header_idx = 0
            for header in headers:
                index_to_header_map[header_idx] = header
                header_idx += 1
                if len(header) > 0:
                    column_values[header] = []
            continue

        line_splits = line.rstrip().split("\t")
        column_count = 0
        for line_split in line_splits:
            header_name = index_to_header_map[column_count]
            if len(header_name) < 1:
                continue
            column_values[header_name].append(line_split)
            column_count += 1

    return (line_count-1, column_values)

def parse_table_without_headers(filename):
    input_file = open(filename, "r")

    line_count = 0
    column_values = {}
    for line in input_file:
        line_splits = line.rstrip().split("\t")

        line_count += 1
        if line_count == 1:
            for i in range(len(line_splits)):
                column_values[i] = []


        column_count = 0
        for line_split in line_splits:
            column_values[column_count].append(line_split)
            column_count += 1

    return (line_count-1, column_values)


#Lists only files in directory, with prefix
def list_files_in_dir(directory):
    onlyfiles = [ os.path.join(directory,f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory,f)) ]
    return onlyfiles

#Lists only folders in directory, with prefix
def list_folders_in_dir(directory):
    onlyfolders = [ os.path.join(directory,f) for f in os.listdir(directory) if os.path.isdir(os.path.join(directory,f)) ]
    return onlyfolders

#List all files in dir
def list_all_in_dir(directory):
    everything = [ os.path.join(directory,f) for f in os.listdir(directory) ]
    return everything

#Returns the leaf filename
def get_only_leaf_filename(path):
    return os.path.basename(path)

def get_filename_without_extension(path):
    return os.path.splitext(path)[0]

def is_path_present(path):
    return os.path.exists(path)

#Making sure directory exists, if not make it.
def make_sure_path_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
