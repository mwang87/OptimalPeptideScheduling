#!/usr/bin/python

import sys
import getopt
import os
import ming_fileio_library

RT_BINS_IN_SECONDS = 5.0

def usage():
    print "<input results file> <input all whole peptides>"

def parse_identification_file(input_results_filename):
    input_results_file = open(input_results_filename, "r")
    peptide_to_rt_map = {}

    grouped_columns_information = []

    lines_to_skip = 11
    header_line = 11

    line_count = 0
    for line in input_results_file:
        line_count += 1
        if line_count < lines_to_skip:
            continue
        splits = line.rstrip().split("\t")
        #Getting header columns
        if line_count == header_line:
            mz_column = -1
            peptide_column = -1
            rt_column = -1
            #Determine column information
            for i in range(len(splits)):
                if splits[i] == "m/z":
                    mz_column = i
                if splits[i] == "DB Peptide":
                    peptide_column = i
                if splits[i] == "RT":
                    rt_column = i
                    grouped_columns_information.append([mz_column, peptide_column, rt_column])
            continue

        #Doing normal lines
        for column_indices in grouped_columns_information:
            mz = splits[column_indices[0]]
            peptide = splits[column_indices[1]]
            rt = splits[column_indices[2]]
            if len(peptide) < 5:
                continue
            rt = float(rt) * 60
            if not peptide in peptide_to_rt_map:
                peptide_to_rt_map[peptide] = []
            peptide_to_rt_map[peptide].append(rt)
    return peptide_to_rt_map

def map_products_to_peptide_rt(peptide_to_rt_map, all_peptides):
    full_peptide_rt_map = {}
    full_peptide_to_products_map = {}
    for product_peptide in peptide_to_rt_map:
        #Finding which full peptide it can be a substring for
        for full_peptide in all_peptides:
            if full_peptide.find(product_peptide) != -1:
                #print product_peptide + " in " + full_peptide
                if not full_peptide in full_peptide_rt_map:
                    full_peptide_rt_map[full_peptide] = {}
                    full_peptide_to_products_map[full_peptide] = []
                full_peptide_rt_map[full_peptide][product_peptide] = peptide_to_rt_map[product_peptide]
                full_peptide_to_products_map[full_peptide].append(product_peptide)
        #print product_peptide + "\t" + str(peptide_to_rt_map[product_peptide])
    #print full_peptide_to_products_map
    #print full_peptide_rt_map
    return full_peptide_rt_map

#This is passed a list of retention times, and in each cell, it is a list of product peptides
def count_number_of_unique_products(rt_timeline_list):
    unique_cell_products = []
    for rt_cell in rt_timeline_list:
        if len(rt_cell) == 1:
            unique_cell_products.append(list(rt_cell)[0])

    return len(set(unique_cell_products))

def count_number_of_empty_rt_slots(rt_timeline_list):
    empty_cell_count = 0
    cell_count = 0
    for rt_cell in rt_timeline_list:
        cell_count += 1
        if len(rt_cell) == 0:
            empty_cell_count += 1
        print str(cell_count) + "\t" + str(len(rt_cell))
    return empty_cell_count


#Add a full peptide to a rt_timeline_list
def add_peptide_to_rt_timeline_list(rt_timeline_list, peptide_products):
    for product in peptide_products:
        for rt in peptide_products[product]:
            rt_index = int(rt/RT_BINS_IN_SECONDS)
            rt_timeline_list[rt_index].add(product)

def get_number_of_unique_rt_windows(list_of_rts):
    list_of_rt_indexes = set()
    for rt in list_of_rts:
        rt_index = int(rt/RT_BINS_IN_SECONDS)
        list_of_rt_indexes.add(rt_index)
    return len(list_of_rt_indexes)



def partition_peptides(full_peptides_to_rt, number_partitions):
    number_of_peptides = len(full_peptides_to_rt)
    print "Number of Peptides: " + str(number_of_peptides)
    partition_peptides_unsorted(full_peptides_to_rt, number_partitions)


#Returns the rt_timeline_list
def partition_peptides_unsorted(full_peptides_to_rt, number_partitions):
    rt_timeline_list = []

    for i in range(int(60 * 120 / RT_BINS_IN_SECONDS)):
        rt_timeline_list.append(set())

    #Lets do a greedy thing
    for peptide in full_peptides_to_rt:
        #print peptide + "\t" + str(full_peptides_to_rt[peptide])
        add_peptide_to_rt_timeline_list(rt_timeline_list, full_peptides_to_rt[peptide])
        #print peptide + "\t" + str(get_number_of_unique_rt_windows(full_peptides_to_rt[peptide]))

    print count_number_of_unique_products(rt_timeline_list)

    return rt_timeline_list

def main():
    input_results_filename = sys.argv[1]
    input_peptide_list_filename = sys.argv[2]

    peptide_to_rt_map = parse_identification_file(input_results_filename)
    line_counts, table_data = ming_fileio_library.parse_table_with_headers(input_peptide_list_filename)
    all_peptides = table_data["Peptides"]

    full_peptides_to_rt = map_products_to_peptide_rt(peptide_to_rt_map, all_peptides)
    partition_peptides(full_peptides_to_rt, 2)








if __name__ == "__main__":
    main()
