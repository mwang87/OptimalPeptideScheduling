#!/usr/bin/python


import sys
import getopt
import os
import ming_fileio_library
from pyteomics import mass

class Feature:
    def __init__(self, rt, mz, intensity):
        self.rt = rt
        self.mz = mz
        self.intensity = intensity

def usage():
    print "<input masses> <input features> <ppm tolerance>"

def load_masses(input_filename):
    masses_list = []
    line_counts, table_data = ming_fileio_library.parse_table_with_headers(input_filename)

    for i in range(line_counts):
        masses_list.append([table_data["Resolaveability"][i], table_data["Peptide"][i], float(table_data["m/z"][i])])

    #Sort this mofo
    sorted_peptide_mass_list = sorted(masses_list, key=lambda pep_obj: pep_obj[2])

    return sorted_peptide_mass_list

def load_features_table(input_filename):
    feature_list = []
    line_counts, table_data = ming_fileio_library.parse_table_with_headers(input_filename)
    for i in range(line_counts):
        feature = Feature(float(table_data["#rt"][i]), float(table_data["mz"][i]), float(table_data["intensity"][i]))
        feature_list.append(feature)
    return feature_list

def identify_features(masses_list, features_list, ppm_tolerance):
    unambiguously_identified = []

    for feature in features_list:
        print "Looking for: " + str(feature.mz)
        mz = feature.mz

        identifications = identify_mass(mz, masses_list, ppm_tolerance)

        if len(identifications) > 0:
            print identifications

        if len(identifications) == 1:
            unambiguously_identified += identifications





        #for mass_obj in masses_list:
            #print mass_obj[2]
            #print str(mz) + '\t' + str(mass_obj[2])
        #    if mz < mass_obj[2]:
        #        print str(mz) + '\t' + str(mass_obj[2])
        #        print "FOUND"

    for identification in unambiguously_identified:
        print "IDENTIFICATION\t" + identification


def identify_mass(mass_input, masses_list, ppm_tolerance):
    number_of_masses = len(masses_list)
    identifications = []
    for i in range(number_of_masses):
        mass_delta = abs(masses_list[i][2] - mass_input)
        ppm_error = mass_delta / mass_input * 1000000
        if ppm_error < ppm_tolerance:
            identifications.append(masses_list[i][1])
    return identifications

def main():
    input_masses_filename = sys.argv[1]
    input_features_filename = sys.argv[2]
    ppm_tolerance = float(sys.argv[3])

    masses_list = load_masses(input_masses_filename)
    features_list = load_features_table(input_features_filename)

    identify_features(masses_list, features_list, ppm_tolerance)



if __name__ == "__main__":
    main()
