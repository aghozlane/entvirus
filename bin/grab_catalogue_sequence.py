#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Extract a list of sequence from the catalogue."""


from __future__ import print_function
import argparse
import sys
import os
import csv
import bisect
import textwrap

__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__credits__ = ["Amine Ghozlane"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    # from Jonathan Barnoud
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.set_defaults(results=".{0}".format(os.sep))
    parser.add_argument('-i', dest='list_sequences_file', type=isfile,
                        required=True,
                        help='List of sequence to extract.')
    parser.add_argument('-d', dest='catalogue_file', type=isfile,
                        required=True, help='Database query.')
    parser.add_argument('-n', dest='not_in_database', action='store_true',
                        help='Select instead elements which are not in the'
                        ' list.')
    parser.add_argument('-id', dest='identity', type=str, default=None,
                        help="Add sample identity")
    parser.add_argument('-e', dest='identity_threshold', type=float,
                        default=0.0, help='Identity threshold (default 0.0)')
    parser.add_argument('-v', dest='coverage_threshold', type=float,
                        default=0.0, help='Coverage threshold (default 0.0)')
    parser.add_argument('-a', dest='vp1_id_file', type=isfile, default=None,
                        help='VP1 database annotation', required=True)
    parser.add_argument('-o', dest='output_file', type=str,
                        help='Output file.')
    parser.add_argument('-r', dest='results', type=isdir,
                        help='Path to result directory.')
    return parser.parse_args()


def load_vp1_id(vp1_id_file):
    """Load vp1 id file
    """
    vp1_id_dict = {}
    try:
        with open(vp1_id_file, "rt") as vp1_id:
            vp1_id_reader = csv.reader(vp1_id, delimiter="\t")
            for line in vp1_id_reader:
                #print(line)
                vp1_id_dict[line[0]] = float(line[2])
    except IOError:
        sys.exit("Error cannot open {0}".format(vp1_id_file))
    return vp1_id_dict


def extract_interest_elements(list_sequences_file, identity_threshold,
                              coverage_threshold, vp1_id_dict):
    """Get a list of the element of interest
    """
    list_sequences = []
    list_reverse_comp = []
    try:
        with open(list_sequences_file, "rt") as list_seq:
            list_sequences_reader = csv.reader(list_seq, delimiter="\t")
            for line in list_sequences_reader:
                if (float(line[8]) >= identity_threshold and
                    round(100.0 * float(line[3])/vp1_id_dict[line[1]]) >= coverage_threshold):
                    list_sequences.append(line[0])
                    # check if the sequence need to be reversed_complemented
                    if int(line[6]) > int(line[7]):
                        list_reverse_comp.append(line[0])
            #assert(len(list_sequences) > 0)
            list_sequences.sort()
            list_reverse_comp.sort()
    except IOError:
        sys.exit("Error cannot the file : {0}".format(list_sequences_file))
    #except AssertionError:
    #    sys.exit("Error no element detected in the file : {0}"
    #             .format(list_sequences_file))
    return list_sequences, list_reverse_comp


# def is_selected(header, list_sequences):
#     """
#     """
#     for element in list_sequences:
#         if element in header:
#             return element
#     return None


def get_element(name, input_list):
    """Search name in input list
      Arguments:
        input_list: List
        name: Search criteria
    """
    # Searching the node with its name
    i = bisect.bisect_left(input_list, name)
    # Object has been found
    if(i != len(input_list) and input_list[i] == name):
        return True #input_list[i]
    return False #None


def extract_catalogue_sequence(list_sequences, catalogue_file, not_in_database):
    """Extract the sequence from the multifasta
    """
    grab_sequence = False
    interest_sequence = {}
    title = ""
    try:
        with open(catalogue_file, "rt") as catalogue:
            for line in catalogue:
                if line[0] == ">":
                    grab_sequence = False
                    title = line[1:].replace("\n", "").replace("\r", "")
                    if " " in title:
                        title = title.split(" ")[0]
                    #print(title)
                    selection = get_element(title, list_sequences)
                    if selection and not not_in_database:
                        interest_sequence[title] = ""
                        grab_sequence = True
                    elif not selection and not_in_database:
                        interest_sequence[title] = ""
                        grab_sequence = True
                elif grab_sequence and len(line) > 0:
                    interest_sequence[title] += line.replace("\n", "").replace("\r", "")
            #assert(len(interest_sequence) > 0)
    except IOError:
        sys.exit("Error cannot the file : {0}".format(catalogue_file))
    #except AssertionError:
    #    sys.exit("Error no element detected in the file : {0}"
    #             .format(catalogue_file))
    return interest_sequence


def fill(text, width=80):
    """Split text"""
    return os.linesep.join(text[i:i+width] for i in xrange(0, len(text), width))

def reverse_complement(dna, complement):
    """Reverse complement a DNA sequence
    """
    return ''.join([complement[base] for base in dna[::-1]])

def rev_comp(interest_sequence, list_reverse_comp):
    """
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    for id_comp in list_reverse_comp:
        interest_sequence[id_comp] = reverse_complement(
                interest_sequence[id_comp], complement)
    return interest_sequence

def write_interest_sequence(interest_sequence, identity, output_file):
    """Write extracted sequence
    """
    idname = ""
    if identity:
        idname = identity + "_"
    try:
        with open(output_file, "wt") as output:
            for key in interest_sequence:
                output.write(">{2}{1}{0}{3}{0}".format(
                                os.linesep, key, idname,
                                fill(interest_sequence[key])))
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get the arguments
    args = getArguments()
    # Load Vp1
    print("Load vp1")
    vp1_id_dict = load_vp1_id(args.vp1_id_file)
    # Get List of sequence of interest
    print("Load the list of sequence of interest ...")
    list_sequences, list_reverse_comp = extract_interest_elements(
        args.list_sequences_file, args.identity_threshold,
        args.coverage_threshold, vp1_id_dict)
    print("{0} sequences to search".format(len(list_sequences)))
    # Extract catalogue sequence
    print("Extract sequences from the catalogue...")
    interest_sequence = extract_catalogue_sequence(list_sequences,
                                                   args.catalogue_file,
                                                   args.not_in_database)
    print("{0} extracted sequences".format(len(interest_sequence)))
    #
    if not args.not_in_database:
        rev_comp(interest_sequence, list_reverse_comp)
    print()
    # Write sequences
    if not args.output_file:
       args.output_file = "extracted_sequence.fasta"
    print("Write sequences to {0}".format(args.output_file))
    write_interest_sequence(interest_sequence, args.identity, args.output_file)
    print("Done.")


if __name__ == '__main__':
    main()
