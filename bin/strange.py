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

"""Extract a list of sequence from the catalogue."""


from __future__ import print_function
import argparse
import sys
import os
import csv
import bisect

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


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.set_defaults(results=".{0}".format(os.sep))
    parser.add_argument('-i', dest='blast_file', type=isfile,
                        required=True,
                        help='List of sequence to extract.')
    parser.add_argument('-c', dest='interest_file', type=str,
                        help="Add sample identity")
    parser.add_argument('-o', dest='output_file', type=str, required=True,
                        help='Output file.')
    return parser.parse_args()

def get_element(input_list, name):                                               
    """Search name in input list                                                 
      Arguments:                                                                 
        input_list: List                                                         
        name: Search criteria                                                    
    """                                                                          
    # Searching the node with its name                                           
    i = bisect.bisect_left(input_list, name)                                     
    # Object has been found                                                      
    if(i != len(input_list) and input_list[i] == name):                          
        return True                                                              
    return False 

def get_keys(interest_file):
    """
    """
    interest_tab = []
    try:                                                                         
        with open(interest_file, "rt") as interest:                             
            interest_reader = csv.reader(interest, delimiter="\t")                   
            for line in interest_reader:                                             
                interest_tab.append(line[0])                                         
            interest_tab.sort()                                                      
    except IOError:                                                              
        sys.exit("Error cannot open {0}".format(interest_file))
    return interest_tab


def write_result(interest_tab, blast_file, output_file):
    """
    """
    try:
        with open(blast_file, "rt") as blast:                                   
            with open(output_file, "wt") as output:
                blast_reader = csv.reader(blast, delimiter="\t")
                output_writer = csv.writer(output, delimiter="\t")
                for line in blast_reader:                                                
                    if get_element(interest_tab, line[0]):                               
                        output_writer.writerow(line)
    except IOError:
        sys.exit("Error cannot open {0}".format(blast_file))

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get the arguments
    args = getArguments()
    interest_tab = get_keys(args.interest_file)
    write_result(interest_tab, args.blast_file, args.output_file)


if __name__ == '__main__':
    main()
