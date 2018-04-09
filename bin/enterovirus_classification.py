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

"""Extract all process result specific to vp1."""

from __future__ import print_function
import os
import sys
import argparse
import csv


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2017, Institut Pasteur"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@pasteur.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
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
            msg = "{0} is a file.".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

#===================
# parameters
#===================
def get_arguments():
    """Extract program options
    """
    parser = argparse.ArgumentParser(description=__doc__,
            usage="{0} -h [options] [arg]"
            .format(sys.argv[0]))
    parser.add_argument('-i', dest='resume_file', type=isfile, required=True,
                        help='Path to the resume file.')
    parser.add_argument('-f', dest='fasta_file', type=isfile, required=True,
                        help='Path to the database file.')
    parser.add_argument('-t', dest='tag', type=str, required=True,
                        choices=["Contigs_with_VP1", "P1_sequences", "VP1_sequences"],
                        help='Tag possible value.')
    parser.add_argument('-r', dest='root_seq_file', type=isfile,
                        help='Path to the sequence to root the tree.')
    parser.add_argument('-o', dest='results', type=isdir,
                        default=os.curdir + os.sep,
                        help='Path to result directory.')
    return parser.parse_args()

def get_unique(seq):
    # Not order preserving
    return {}.fromkeys(seq).keys()

def get_query(query_file, tag):
    """Load resume file
    """
    classify_list = []
    query_dict = {}
    try:
        with open(query_file, "rt") as query:
            query_reader = csv.reader(query, delimiter="\t")
            header = query_reader.next()
            interest_seq_posit = header.index(tag)
            interest_serotype_posit = header.index("Serotype_VP1")
            interest_specie_posit = header.index("Specie_VP1")
            for line in query_reader:
                line_len = len(line)
                #print(line)
                if line_len > interest_seq_posit and line_len > interest_serotype_posit:
                    if line[interest_seq_posit] != "":
                        #print(line)
                        assert(line[interest_seq_posit] not in query_dict)
                        query_dict[line[interest_seq_posit]] = [
                            line[interest_serotype_posit].upper(),
                            line[interest_specie_posit].upper()]
                        # Serotype and specie
                        classify_list += [line[interest_serotype_posit].upper(),
                                          line[interest_specie_posit].upper()]
            classify_list.sort()
            classify_list = get_unique(classify_list)
    except IOError:
        sys.exit("Error cannot open {0}".format(query_file))
    except AssertionError:
        print(line[interest_seq_posit])
        print(query_dict)
        sys.exit("Strange value in {0}".format(query_file))
    return query_dict, classify_list


def get_sequence(query_dict, fasta_file):
    """
    """
    result = {}
    title = None
    try:
        with open(fasta_file, "rt") as target:
            for line in target:
                if title and line[0] != ">":
                    result[title] += line[0:].strip().replace("\n", "").replace("\r", "")
                if line.startswith(">"):
                    if len(result) == len(query_dict):
                        break
                    title = None
                    title = line[1:].replace("\n", "").replace("\r", "")
                    if " " in title:
                        title = title.split(" ")[0]
                    if title in query_dict:
                        result[title] = ""
            assert(len(result) == len(query_dict))
    except IOError:
        sys.exit("Error : cannot open {0}".format(fasta_file))
    except AssertionError:
        print("The program has not find every query sequence, "
            "check the following elements :", file=sys.stderr)
    return result


def fill(text, width=80):
    """Split text"""
    return os.linesep.join(text[i:i+width] for i in xrange(0, len(text), width))


def write_sequence(results, sequence_data, query_dict, classify_list, tag,
                   root_seq_file):
    """
    """
    root_lines = ""
    output_list = []
    try:
        for clas in classify_list:
            output_file = results + os.sep + clas + "_" + tag + ".fasta"
            if root_seq_file:
                output_file = results + os.sep + clas + "_" + tag + "_rooted.fasta"
                with open(root_seq_file, "rt") as root_seq:
                    root_lines = root_seq.readlines()
            output_list += [open(output_file, "wt")]
        for seq in sequence_data:
            serotype_posit = classify_list.index(query_dict[seq][0])
            output_list[serotype_posit].write(">{1}{0}{2}{0}".format(
                os.linesep, seq, fill(sequence_data[seq])))
            specie_posit = classify_list.index(query_dict[seq][1])
            output_list[specie_posit].write(">{1}{0}{2}{0}".format(
                os.linesep, seq, fill(sequence_data[seq])))
        
        for output in output_list:
            output.writelines(root_lines)
            output.close()
    except IOError:
        sys.exit("Error : cannot open {0}".format(output_file))


def main():
    """
    """
    tag_dict = {"Contigs_with_VP1":"contigs", "P1_sequences":"p1",
                "VP1_sequences":"vp1"}
    args = get_arguments()
    # Load query elements
    print("Load resume file")
    query_dict, classify_list = get_query(args.resume_file, args.tag)
    # Grab query sequence in the database
    print("Load database sequence")
    sequence_data = get_sequence(query_dict, args.fasta_file)
    # Write the new fasta file
    print("Write the new fasta")
    write_sequence(args.results, sequence_data, query_dict, 
                   classify_list, tag_dict[args.tag], args.root_seq_file)
    print("Done")

if __name__ == '__main__':
    main()
