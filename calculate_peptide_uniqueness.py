#!/usr/bin/python


import sys
import getopt
import os
import ming_fileio_library
from pyteomics import mass

def usage():
    print "<input file> <ppm tolerance>"


def find_all_substring_of_length(input_string, length):
    return_list = []
    for i in range(len(input_string) - length + 1):
        return_list.append(input_string[i:i+length])

    return return_list

def find_resolveable_peptides(peptide_mass_map, ppm_tolerance):
    print "ppm_tolerance " + str(ppm_tolerance)

    peptide_mass_list = []
    for peptide in peptide_mass_map:
        peptide_mass_list.append([peptide, peptide_mass_map[peptide]])

    sorted_peptide_mass_list = sorted(peptide_mass_list, key=lambda pep_obj: pep_obj[1])
    ambiguous_peptides = []
    resolveable_peptides = []
    #Iterate and check if adjacent elements are more than PPM apart
    for i in range(len(sorted_peptide_mass_list)):
        current_object = sorted_peptide_mass_list[i]

        if i > 0:
            #Check Left
            other_object = sorted_peptide_mass_list[i-1]
            mass_difference = abs(other_object[1] - current_object[1])
            ppm_difference = mass_difference / other_object[1] * 1000000
            if ppm_difference < ppm_tolerance:
                ambiguous_peptides.append(current_object)
                continue
            #print "LEFT: " + str(ppm_difference)

        if i < len(sorted_peptide_mass_list) - 1:
            #Check Right
            other_object = sorted_peptide_mass_list[i+1]
            mass_difference = abs(other_object[1] - current_object[1])
            ppm_difference = mass_difference / other_object[1] * 1000000
            if ppm_difference < ppm_tolerance:
                ambiguous_peptides.append(current_object)
                continue
            #print "RIGHT: " + str(ppm_difference)

        #Its valid
        resolveable_peptides.append(current_object)

    for peptide in resolveable_peptides:
        print "RESOLVEABLE\t" + peptide[0] + "\t" + str(peptide[1])

    for peptide in ambiguous_peptides:
        print "AMBIGUOUS\t" + peptide[0] + "\t" + str(peptide[1])


def main():
    input_filename = sys.argv[1]
    ppm_tolerance = float(sys.argv[2])
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
    peptide_mass_map = {}
    for peptide in all_sub_peptides:
        peptide_key = peptide + ".2"
        peptide_mass = mass.calculate_mass(sequence=peptide, ion_type='M', charge=2)
        peptide_mass_map[peptide_key] = peptide_mass

        peptide_key = peptide + ".3"
        peptide_mass = mass.calculate_mass(sequence=peptide, ion_type='M', charge=3)
        peptide_mass_map[peptide_key] = peptide_mass

        peptide_key = peptide + ".4"
        peptide_mass = mass.calculate_mass(sequence=peptide, ion_type='M', charge=4)
        peptide_mass_map[peptide_key] = peptide_mass


        #print peptide + "\t" + "2" + "\t" + str(mass.calculate_mass(sequence=peptide, ion_type='M', charge=2))
        #print peptide + "\t" + "3" + "\t" + str(mass.calculate_mass(sequence=peptide, ion_type='M', charge=3))
        #print peptide + "\t" + "4" + "\t" + str(mass.calculate_mass(sequence=peptide, ion_type='M', charge=4))


    #Determine uniqueness
    find_resolveable_peptides(peptide_mass_map, ppm_tolerance)






if __name__ == "__main__":
    main()
