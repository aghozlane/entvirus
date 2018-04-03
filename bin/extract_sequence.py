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


"""Extract EggNogg annotation."""



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
    parser.add_argument('-q', dest='contigs_file', type=isfile, required=True,
                        help="Contigs fasta file")
    parser.add_argument('-b', dest='blast_result_file', type=isfile,
                        required=True, help="Blast result file")
    parser.add_argument('-a', dest='vp1_info_file', type=isfile,
                        required=True, help="Blast result file")
    parser.add_argument('-id', dest='identity', type=str, default=None,
                        help="Sample identity")
    parser.add_argument('-i', dest='identity_threshold', type=float, default=0.0,
                        help="Identity threshold (default >=0.0)")
    parser.add_argument('-c', dest='coverage_threshold', type=float, default=0.0,
                        help="Coverage threshold (default >=0.0)")
    parser.add_argument('-o', dest='output_file', type=str, required=True,
                        help='Output file')
    return parser.parse_args()


def load_blast(blast_result_file, identity_threshold, coverage_threshold,
               vp1_id_dict):
    """Load blast position identification and filter
    """
    position_dict = {}
    reverse_comp = False
    diff_length = 0
    try:
        with open(blast_result_file, "rt") as blast_result:
            blast_reader = csv.reader(blast_result, delimiter='\t')
            for line in blast_reader:
                #print(line)
                #print("id:" + str(round(float(line[8]),1)))
                #print("cov:" +str(round(100.0 * float(line[3])/vp1_id_dict[line[1]][1])))
                if (round(float(line[8]),1) >= identity_threshold and
                    round(100.0 * float(line[3])/vp1_id_dict[line[1]][1]) >= coverage_threshold):
                    sstart = int(line[6])
                    send = int(line[7])
                    qlen = int(line[2])
                    length = int(line[3])
                    qstart = int(line[4])
                    qend = int(line[5])
                    diff_length = 0
                    start_position = qstart
                    end_position = qend
                    # Contigs is in reverse position
                    if sstart > send:
                        reverse_comp = True
                        # check if vp1 need to be extended
                        if vp1_id_dict[line[1]][1] > sstart:
                            diff_length = vp1_id_dict[line[1]][1] - sstart
                            # need to extend query start
                            extended = qend + diff_length
                            if extended > qlen:
                                diff_length = qlen - length
                            start_position = qstart - diff_length
                    else:
                        reverse_comp = False
                        if vp1_id_dict[line[1]][1] > send:
                            diff_length = vp1_id_dict[line[1]][1] - send
                            # need to extend query end
                            extended = qend + diff_length
                            if extended > qlen:
                                diff_length = qlen - length
                            end_position = qend + diff_length
                    # report feature
                    position_dict[line[0]] = {
                        "position":[start_position-1, end_position],
                        "target":line[1],
                        "evalue":line[10],
                        "identity":round(float(line[8]), 1),
                        "coverage":round(100.0 * float(length)/float(vp1_id_dict[line[1]][1]), 1),
                        "reverse_comp": reverse_comp}
                        #"coverage":round(100.0 * float(line[3]) / float(line[2]),1)}
            #print(position_dict)
    except IOError:
        sys.exit("Error cannot open {0}".format(blast_result_file))
    return position_dict


def load_sequence(contigs_file):
    """Extract vp1 sequence from contigs
    """
    head = ""
    contig = ""
    try:
        with open(contigs_file, "rt") as contigs:
            for line in contigs:
                if line.startswith(">"):
                    if head != "":
                        yield [head, contig]
                    head = line[1:].replace("\n", "")
                    contig = ""
                elif len(line) > 0:
                    contig += line.replace("\n", "")
            if len(contig) > 0:
                yield [head, contig]
    except IOError:
        sys.exit("Error cannot open {0}".format(contigs_file))


def fill(text, width=80):
    """Split text"""
    return os.linesep.join(text[i:i+width] for i in xrange(0, len(text), width))


def load_vp1_id(vp1_id_file):
    """Load vp1 id file
    """
    vp1_id_dict = {}
    try:
        with open(vp1_id_file, "rt") as vp1_id:
            vp1_id_reader = csv.reader(vp1_id, delimiter="\t")
            for line in vp1_id_reader:
                #print(line)
                vp1_id_dict[line[0]] = [line[1], int(line[2])]
    except IOError:
        sys.exit("Error cannot open {0}".format(vp1_id_file))
    return vp1_id_dict

def rev_comp(dna, complement):
    """Reverse complement a DNA sequence
    """
    return ''.join([complement[base] for base in dna[::-1]])

def extract_sequence(position_dict, contigs_file, output_file, identity):
    """
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    contig_seq = ""
    idname = ""
    if identity:
        idname = identity + "_"
    try:
        with open(output_file, "wt") as output:
            for contigid, contig in load_sequence(contigs_file):
                if contigid in position_dict:
                    if position_dict[contigid]["reverse_comp"]:
                        contig_seq = rev_comp(
                                contig[position_dict[contigid]["position"][0]:position_dict[contigid]["position"][1]],
                                complement)
                    else:
                        contig_seq = contig[position_dict[contigid]["position"][0]:position_dict[contigid]["position"][1]]
                    output.write(">{1}_vp1 match:{2} identity:{3}% "
                                 "coverage:{4}% evalue:{5}{0}{6}{0}".format(
                        os.linesep, idname + contigid,
                        position_dict[contigid]["target"],
                        position_dict[contigid]["identity"],
                        position_dict[contigid]["coverage"],
                        position_dict[contigid]["evalue"],
                        fill(contig_seq)))
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Load vp1
    vp1_id_dict = load_vp1_id(args.vp1_info_file)
    #print(vp1_id_dict)
    # Load blast info
    position_dict = load_blast(args.blast_result_file, args.identity_threshold,
                               args.coverage_threshold, vp1_id_dict)
    #print(position_dict)
    # Extract sequence
    if len(position_dict) > 0:
        extract_sequence(position_dict, args.contigs_file, args.output_file,
                        args.identity)
    else:
        open(args.output_file, 'w').close()

if __name__ == '__main__':
    main()
