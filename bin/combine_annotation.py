#!/usr/bin/env python2
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


"""Find best annotation."""

from __future__ import print_function
import argparse
import os
import sys
import csv


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2017, Institut Pasteur"
__license__ = "GPL"
__version__ = "0.0.1"
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


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument('-i', dest='annotation_list', type=isfile,
                        required=True, nargs='+', default=[],
                        help="NCBI annotation file")
    parser.add_argument('-o', dest='output_file', type=str, required=True,
                        help='Output file')
    return parser.parse_args()


def load_blast(annotation_file, annotation_dict):
    """Load ruffly the annotation
    """
    try:
        with open(annotation_file, "rt") as annotation:
            annotation_reader = csv.reader(annotation, delimiter="\t")
            annotation_reader.next()
            for line in annotation_reader:
                if line[0] in annotation_dict:
                    annotation_dict[line[0]] += [line[1:-4] +
                                                [float(line[-4]), float(line[-3]),
                                                float(line[-2]), float(line[-1])]]
                else:
                    annotation_dict[line[0]] = [line[1:-4] +
                                                [float(line[-4]), float(line[-3]),
                                                float(line[-2]), float(line[-1])]]
    except IOError:
        sys.exit("Error cannot open {0}".format(annotation_file))
    return annotation_dict


def write_result(annotation_dict, output_file):
    """Write annotation results
    """
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter="\t")
            #output_writer.writerow(["ContigName", "GI", "superkingdom", "kingdom",
            #                        "phylum", "class", "order", "family",
            #                        "genus","species", "PourcID",
            #                        "Coverage", "evalue", "bitscore"])
            for contigs in annotation_dict:
                annotation_dict[contigs].sort(key=lambda x: x[-4] + x[-3], reverse=True)
                output_writer.writerow([contigs] +
                                       annotation_dict[contigs][0])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    annotation_dict = {}
    # Get arguments
    args = get_arguments()
    # Get the best hit blast info
    for annotation_file in args.annotation_list:
        annotation_dict = load_blast(annotation_file, annotation_dict)
    write_result(annotation_dict, args.output_file)

if __name__ == '__main__':
    main()
