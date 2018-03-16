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
from __future__ import print_function
import os
import sys
import argparse


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2015, Institut Pasteur"
__credits__ = ["Amine Ghozlane"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@pasteur.fr"
__status__ = "Developpement"


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument('-i', dest='fasta_file', type=str, required=False,
                        help='Path to the fasta file.')
    parser.add_argument('-s', dest='sample_name', type=str, default="",
                        help='New name.', required=True)
    parser.add_argument('-o', dest='output_file', type=str, default=None,
                        help='Output file.')
    args = parser.parse_args()
    return args


def fill(text, width=80):
    """Split text"""
    return os.linesep.join(text[i:i+width] for i in xrange(0, len(text), width))



def rename_fasta(fasta_file, sample_name, output_file):
    """
    """
    if not output_file:
        output = sys.stdout
    else:
        output = open(output_file, "wt")
    seq = ""
    header = ""
    line = ""
    try:
        with open(fasta_file, "rt") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    if header != "":
                        print(">{0}|{1}\n{2}".format(sample_name, header,
                                fill(seq)), file=output)
                        seq = ""
                    header=line[1:].replace("\n","").replace("\r", "")
                elif len(line) > 0:
                    seq += line.replace("\n", "").replace("\r", "")
            if len(line) > 0:
                print(">{0}|{1}\n{2}".format(sample_name, header,
                      fill(seq)), file=output)
    except IOError:
        sys.exit("Error cannot open {0}".format(fasta_file))
    if output_file:
        output.close()


def main():
    """Main program
    """
    args = getArguments()
    if os.path.isfile(args.fasta_file):
        rename_fasta(args.fasta_file, args.sample_name, args.output_file)
    else:
        # Open an empty file
        open(args.output_file, "wt").close()
        open(args.fasta_file, "wt").close()

if __name__ == '__main__':
    main()
