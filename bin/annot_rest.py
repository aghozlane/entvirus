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

"""Extract the annotation in case of blast on imomi database."""

from __future__ import print_function
import os
import sys
import argparse
import csv
import glob

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
      :Parameters:
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
    parser.add_argument('-i', dest='blast_dir',
                        type=isdir, required=True, default="./",
                        help='Input blast result file, in m8 mode.')
    parser.add_argument('-db', dest='db',
                        type=isfile, required=True,
                        help='Input blast result file, in m8 mode.')
    parser.add_argument('-id', dest='identity', type=str, default="",
                        help="Sample identity")
    parser.add_argument('-o', '--output_file', dest='output_file', type=str,
                        help='Output file')
    return parser.parse_args()

def load_database(db_file):
    """
    """
    db_dict = {}
    try:
        with open(db_file, "rt")  as db:
            for line in db:
                if line.startswith(">"):
                    data=line.split(" ")
                    db_dict[data[0][1:].replace("\n", "")] = " ".join(data[1:]).replace("\n", "")
    except IOError:
        sys.exit("Err")
    return db_dict

def load_annotation(blast_file, name, annot_dict):
    try:
        with open(blast_file, "rt") as blast:
            blast_reader = csv.reader(blast, delimiter="\t")
            for line in blast_reader:
                contigname = name + "_" +line[0]
                if  contigname in annot_dict:
                    annot_dict[contigname] += [[line[1], float(line[10]), float(line[11])]]
                else:
                    annot_dict[contigname] = [[line[1], float(line[10]), float(line[11])]]
    except IOError:
        sys.exit("Error")
    return annot_dict

def write_annot(annot_dict, db_dict, output_file):
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter="\t")
            output_writer.writerow(["ContigName", "match", "identity", "coverage", "annotation"])
            for contigs in annot_dict:
                annot_dict[contigs].sort(key=lambda x: x[1] + x[2], reverse=True)
                #print(annot_dict[contigs][0])
                if annot_dict[contigs][0][0] in db_dict:
                    output_writer.writerow([contigs] + annot_dict[contigs][0] + [db_dict[annot_dict[contigs][0][0]]])
                else:
                    #pass
                    print(annot_dict[contigs][0][0])
    except IOError:
        sys.exit("Err")

def check_directory(path):
    """Get the list of files in the directory
      Arguments:
          pdb_dir: Path to the directory that contains PDB files
      Returns: List of PDB files.
    """
    return glob.glob('{0}{1}*_polston.tsv'.format(path, os.sep))


def main():
    """
    """
    args = get_arguments()
    blast_result_file = check_directory(args.blast_dir)
    db_dict = load_database(args.db)
    annot_dict = {}
    for blast_file in blast_result_file:
        name = os.path.basename(blast_file).replace("_polston.tsv", "")
        annot_dict = load_annotation(blast_file, name, annot_dict)
    write_annot(annot_dict, db_dict, args.output_file)
if __name__ == "__main__":
    main()

