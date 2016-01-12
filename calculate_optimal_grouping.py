#!/usr/bin/python

import sys
import getopt
import os
import ming_fileio_library
import maxflow

RT_BINS_IN_SECONDS = 10.0

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


#Returns partition of peptides as a list, randomly assigning them
def partition_peptides_random(full_peptides_to_rt, number_partitions):
    peptide_lists = []
    for i in range(number_partitions):
        peptide_lists.append([])

    #Randomly partition
    peptide_number = 0
    for peptide in full_peptides_to_rt:
        partition_index = peptide_number % number_partitions
        peptide_number += 1
        peptide_lists[partition_index].append(peptide)

    return peptide_lists

#Returns partition of peptides as a list, alternating most abundant products
def partition_peptides_number_products(full_peptides_to_rt, number_partitions):
    peptide_lists = []
    for i in range(number_partitions):
        peptide_lists.append([])

    peptide_sorting_list = []
    #Randomly partition
    for peptide in full_peptides_to_rt:
        peptide_sorting_list.append([peptide, len(full_peptides_to_rt[peptide])])

    peptide_sorting_list = sorted(peptide_sorting_list, key=lambda peptide: peptide[1], reverse=True)
    print peptide_sorting_list

    peptide_number = 0
    for peptide in peptide_sorting_list:
        partition_index = peptide_number % number_partitions
        peptide_number += 1
        peptide_lists[partition_index].append(peptide[0])

    return peptide_lists


def count_number_of_acquireable_products(peptide_list, peptide_to_product_rt):
    graph = maxflow.GraphInt()

    #Doing bipartite matching
    product_to_node_idx_mapping = {}
    rt_to_node_idx_mapping = {}
    for peptide in peptide_list:
        if not peptide in peptide_to_product_rt:
            continue

        product_list = peptide_to_product_rt[peptide]
        for product in product_list:
            if not product in product_to_node_idx_mapping:
                #insert it into the graph
                node_added = graph.add_nodes(1)
                #print node_added[0]
                product_to_node_idx_mapping[product] = node_added[0]

            for rt in product_list[product]:
                #print rt
                rt_index = int(rt/RT_BINS_IN_SECONDS)
                #print rt_index
                if not rt_index in rt_to_node_idx_mapping:
                    node_added = graph.add_nodes(1)
                    rt_to_node_idx_mapping[rt_index] = node_added[0]

                #Create edge between product and the rt
                graph.add_edge(rt_to_node_idx_mapping[rt_index], product_to_node_idx_mapping[product], 1, 1)

    #Adding connections to sources
    for product in product_to_node_idx_mapping:
        graph.add_tedge(product_to_node_idx_mapping[product], 1, 0)

    #Adding connections to sink
    for rt_index in rt_to_node_idx_mapping:
        graph.add_tedge(rt_to_node_idx_mapping[rt_index], 0, 1)

    return graph.maxflow()

def main():
    input_results_filename = sys.argv[1]
    input_peptide_list_filename = sys.argv[2]

    products_to_rt_map = parse_identification_file(input_results_filename)
    line_counts, table_data = ming_fileio_library.parse_table_with_headers(input_peptide_list_filename)
    all_peptides = table_data["Peptides"]

    full_peptides_to_rt = map_products_to_peptide_rt(products_to_rt_map, all_peptides)
    partitioned_peptide_list = partition_peptides_random(full_peptides_to_rt, 3)
    #partitioned_peptide_list = partition_peptides_number_products(full_peptides_to_rt, 3)

    print "Total Products: " + str(len(products_to_rt_map))
    total_detectable_products = 0
    for peptide_list in partitioned_peptide_list:
        number_products_detectable = count_number_of_acquireable_products(peptide_list, full_peptides_to_rt)
        #print number_products_detectable
        total_detectable_products += number_products_detectable
    print "Total Products Detectable: " + str(total_detectable_products)

    for peptide_list in partitioned_peptide_list:
        print "Partition================="
        for peptide in peptide_list:
            print peptide









if __name__ == "__main__":
    main()
