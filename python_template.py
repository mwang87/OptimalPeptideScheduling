#!/usr/bin/python


import sys
import getopt
import os
import ming_fileio_library
from pyteomics import mass

def usage():
    print "<input file>"


def find_all_substring_of_length(input_string, length):
    return_list = []
    for i in range(len(input_string) - length + 1):
        return_list.append(input_string[i:i+length])

    return return_list


def main():
    input_filename = sys.argv[1]
    line_counts, table_data = ming_fileio_library.parse_table_with_headers(input_filename)

    all_sub_peptides = []

    for i in range(line_counts):
        #print table_data["Peptides"][i]
        for length in range(10):
            peptide = table_data["Peptides"][i]
            substrings = find_all_substring_of_length(peptide, length + 4)
            #print peptide + "\t" + str(substrings)
            all_sub_peptides += substrings

    #print len(all_sub_peptides)
    all_sub_peptides = list(set(all_sub_peptides))
    #print len(all_sub_peptides)
    for peptide in all_sub_peptides:
        print peptide + "\t" + "2" + "\t" + str(mass.calculate_mass(sequence=peptide, ion_type='M', charge=2))
        print peptide + "\t" + "3" + "\t" + str(mass.calculate_mass(sequence=peptide, ion_type='M', charge=3))
        print peptide + "\t" + "4" + "\t" + str(mass.calculate_mass(sequence=peptide, ion_type='M', charge=4))






if __name__ == "__main__":
    main()
