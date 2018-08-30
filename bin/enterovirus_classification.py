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

"""Extract all process result specific to vp1."""

from __future__ import print_function
import os
import sys
import argparse
import csv
import random
from colour import Color


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
                        choices=["Contigs_with_VP1", "P1_sequences",
                                 "VP1_sequences",  "5UTR_sequences",
                                 "3D_sequences"],
                        help='Tag possible value.')
    parser.add_argument('-e', dest='ent_serotype_file', type=isfile,
                        help='Path to enterovirus specie serotype association')
    parser.add_argument('-r', dest='template_seq_file', type=isfile,
                        help='Path to the sequence to root the tree.')
    parser.add_argument('-o', dest='results', type=isdir,
                        default=os.curdir + os.sep,
                        help='Path to result directory.')
    parser.add_argument('-itol', dest='itol_dir',  type=isdir,
                        default=os.curdir + os.sep,
                        help='Path to itol_label file.')
    parser.add_argument('-c', dest='incomplete', action='store_false',
                        default=True, help='Consider only complete sequence')
    return parser.parse_args()

def get_unique(seq):
    # Not order preserving
    return {}.fromkeys(seq).keys()

def get_query(query_file, tag, incomplete):
    """Load resume file
    """
    serotype_list = []
    classify_list = []
    classify_specie_list = []
    query_dict = {}
    try:
        with open(query_file, "rt") as query:
            query_reader = csv.reader(query, delimiter="\t")
            header = query_reader.next()
            interest_seq_posit = header.index(tag)
            #print(header)
            interest_serotype_posit = header.index("Serotype_VP1")
            interest_specie_posit = header.index("Specie_VP1")
            matched_5utr = header.index('Matched_5UTR')
            matched_3d = header.index('Matched_3D')
            for line in query_reader:
                line_len = len(line)
                #print(line)
                if line_len > interest_seq_posit and line_len > interest_serotype_posit:
                    if line[interest_seq_posit] != "":
                        assert(line[interest_seq_posit] not in query_dict)
                        if line[matched_5utr] != "" and line[matched_3d] != "" or incomplete:
                            #print(line)
                            query_dict[line[interest_seq_posit]] = [
                             line[interest_serotype_posit].upper(),
                             line[interest_specie_posit].upper()]
                            # Serotype and specie
                            classify_list += [line[interest_serotype_posit].upper(),
                                          line[interest_specie_posit].upper()]
                            classify_specie_list += [line[interest_specie_posit].upper()]
                            serotype_list += [line[interest_serotype_posit].upper()]
            classify_list.sort()
            classify_list = get_unique(classify_list)
            classify_specie_list.sort()
            classify_specie_list = get_unique(classify_specie_list)
    except IOError:
        sys.exit("Error cannot open {0}".format(query_file))
    except AssertionError:
        print(line[interest_seq_posit])
        print(query_dict)
        sys.exit("Strange value in {0}".format(query_file))
    # print("here")
    # print(query_dict)
    # print("classify_list")
    # print(classify_list)
    # print("classify specie list")
    # print(classify_specie_list)
    # print("serotype")
    # print(serotype_list)
    return query_dict, classify_list, classify_specie_list, serotype_list


def get_sequence(query_dict, fasta_file):
    """
    """
    result = {}
    title = None
    try:
        with open(fasta_file, "rt") as target:
            for line in target:
                if title and line[0] != ">":
                    if title in result:
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
                   ref_seq, ent_spe_sero):
    """
    """
    output_list = []
    #save_association = {}
    try:
        for clas in classify_list:
            output_file = results + os.sep + clas + "_" + tag + ".fasta"
            root_lines = ""
            if ref_seq:
                output_file = results + os.sep + clas + "_" + tag + "_rooted.fasta"
                output_association_file = results + os.sep + clas + "_" + tag + "_association.txt"
                # case specie
                if clas in ref_seq:
                    specie_list = ref_seq.keys()
                    specie_list.remove(clas)
                    choosen_specie = random.choice(specie_list)
                    list_serotype = ref_seq[choosen_specie].keys()[0:5]
                    for key in list_serotype:
                        root_lines += ">{1}{0}{2}{0}".format(os.linesep, key[1], fill(ref_seq[choosen_specie][key]))
                    #save_association[clas] = [i[1] for i in list_serotype]
                #print(clas)
                else:
                    specie_list = ref_seq.keys()
                    specie_list.remove(ent_spe_sero[clas])
                    choosen_specie = random.choice(specie_list)
                    list_serotype = ref_seq[choosen_specie].keys()[0:5]
                    for key in list_serotype:
                        root_lines += ">{1}{0}{2}{0}".format(os.linesep, key[1], fill(ref_seq[choosen_specie][key]))
                    #save_association[clas] = [i[1] for i in list_serotype]
                with open(output_association_file, "wt") as output_association:
                    for line in list_serotype:
                        output_association.write("{0}{1}".format(line[1], os.linesep))
            output=open(output_file, "wt")
            output.writelines(root_lines)
            output_list += [output]
        for seq in sequence_data:
            serotype_posit = classify_list.index(query_dict[seq][0])
            output_list[serotype_posit].write(">{1}{0}{2}{0}".format(
                os.linesep, seq, fill(sequence_data[seq])))
            specie_posit = classify_list.index(query_dict[seq][1])
            output_list[specie_posit].write(">{1}{0}{2}{0}".format(
                os.linesep, seq, fill(sequence_data[seq])))
        for output in output_list:
            #output.writelines(root_lines)
            output.close()
    except IOError:
        sys.exit("Error : cannot open {0}".format(output_file))
    #return save_association


def write_itol_label(results, sequence_data, query_dict, classify_list, tag):#, save_association):
    """Write label with sample name in itol
    """
    root_lines = None
    output_list = []
    try:
        for clas in classify_list:
            output_file = results + os.sep + clas + "_" + tag + "_LABELS.txt"
            #if ref_seq:
            #    output_file = results + os.sep + clas + "_" + tag + "_rooted_LABELS.txt"

                #with open(root_seq_file, "rt") as root_seq:
                #    root_lines = [line[1:].split(" ")[0].replace("\n", "").replace("\r", "")
                #                 for line in root_seq if line.startswith(">")]
            output_list += [open(output_file, "wt")]
        for output in output_list:
            output.write("LABELS{0}SEPARATOR TAB{0}DATA{0}".format(os.linesep))
        for seq in sequence_data:
            serotype_posit = classify_list.index(query_dict[seq][0])
            output_list[serotype_posit].write("{1}\t{2}{0}".format(
                os.linesep, seq.replace("|", "_"), seq.split("|")[0]))
            specie_posit = classify_list.index(query_dict[seq][1])
            output_list[specie_posit].write("{1}\t{2}{0}".format(
                os.linesep, seq.replace("|", "_"), seq.split("|")[0]))
        for output in output_list:
            #if root_lines:
            #    for root in root_lines:
            #        output.write("{1}\t{1}{0}".format(os.linesep, root))
            output.close()
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


def write_itol_tree_color(results, sequence_data, query_dict, classify_list, serotype_list,
                          tag):#, ref_seq):
    """
    """
    root_lines = None
    red = Color("red")
    blue = Color("lime")
    output_list = []
    #if ref_seq:
    #    color_list = list(red.range_to(blue, len(serotype_list) + 1))
    #else:
    try:
        color_list = list(red.range_to(blue, len(serotype_list)))
        for clas in classify_list:
            output_file = results + os.sep + clas + "_" + tag + "_TREE_COLORS.txt"
            #if ref_seq:
            #    output_file = results + os.sep + clas + "_" + tag + "_rooted_TREE_COLORS.txt"
                #with open(root_seq_file, "rt") as root_seq:
                #    root_lines =[line[1:].split(" ")[0].replace("\n", "").replace("\r", "")
                #                 for line in root_seq if line.startswith(">")]
            output_list += [open(output_file, "wt")]
        for output in output_list:
            output.write("TREE_COLORS{0}SEPARATOR TAB{0}DATA{0}".format(os.linesep))
        for seq in sequence_data:
            specie_posit = classify_list.index(query_dict[seq][1])
            output_list[specie_posit].write("{1}\trange\t{2}\t{3}{0}".format(
                os.linesep, seq.replace("|", "_"), color_list[serotype_list.index(query_dict[seq][0])].hex_l, query_dict[seq][0]))
        for output in output_list:
            #if root_lines:
            #    for root in root_lines:
            #        output.write("{1}\tlabel\t{2}{0}".format(os.linesep, root, "#000000"))
            output.close()
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))
    except ValueError:
        pass


def load_spe_sero(ent_serotype_file):
    """Specie serotype association
    """
    ent_spe_sero = {}
    try:
        with open(ent_serotype_file, "rt") as ent_serotype:
            ent_serotype_reader = csv.reader(ent_serotype, delimiter="\t")
            # pass header
            ent_serotype_reader.next()
            for line in ent_serotype_reader:
                ent_spe_sero[line[1]] = line[0]
            assert(len(ent_spe_sero) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(ent_serotype_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(ent_serotype_file))
    return ent_spe_sero


def get_template_sequence(template_seq_file, ent_spe_sero):
    """Load template sequence
    """
    ref_seq = {}
    seq = ""
    try:
        with open(template_seq_file, "rt") as template_seq:
            for line in template_seq:
                if line.startswith(">"):
                    data_seq = line.split(" ")
                    serotype = data_seq[1].replace("\n", "").replace("\r", "").strip()
                    if len(seq) > 0:
                        if ent_spe_sero[serotype] in ref_seq:
                            ref_seq[ent_spe_sero[serotype]].update({(serotype, data_seq[0][1:]):seq})
                        else:
                            ref_seq[ent_spe_sero[serotype]] = {(serotype, data_seq[0][1:]):seq}
                        seq = ""
                elif len(line) > 0:
                    seq += line.replace("\n", "").replace("\r", "")
    except IOError:
        sys.exit("Error cannot open {0}".format(template_seq_file))
    return ref_seq


def main():
    """Reclassify sequence and prepare phylogeny
    """
    ref_seq = {}
    ent_spe_sero = {}
    tag_dict = {"Contigs_with_VP1":"contigs", "P1_sequences":"p1",
            "VP1_sequences":"vp1", "5UTR_sequences":"5utr", "3D_sequences":"3d"}
    args = get_arguments()
    # Load query elements
    print("Load resume file")
    (query_dict, classify_list,
     classify_specie_list, serotype_list) = get_query(args.resume_file,
                                                      args.tag,
                                                      args.incomplete)
    print("{} descriptions loaded".format(len(query_dict)))
    # Load specie association
    if args.ent_serotype_file and args.template_seq_file:
        # Load enterovirus serotype
        print("Load enterovirus serotype association")
        ent_spe_sero = load_spe_sero(args.ent_serotype_file)
        # Load template sequence
        print("Load template sequence")
        ref_seq = get_template_sequence(args.template_seq_file, ent_spe_sero)
    # Grab query sequence in the database
    print("Load database sequence")
    sequence_data = get_sequence(query_dict, args.fasta_file)
    print("{} sequences loaded".format(len(sequence_data)))
    # Write the new fasta file
    print("Write the new fasta")
    write_sequence(args.results, sequence_data, query_dict, classify_list,
                   tag_dict[args.tag], ref_seq, ent_spe_sero)
    #print(save_association)
    print("Write the itol label")
    write_itol_label(args.itol_dir, sequence_data, query_dict, classify_list,
                     tag_dict[args.tag])
    print("Write the itol tree color")
    write_itol_tree_color(args.itol_dir, sequence_data, query_dict, classify_specie_list, serotype_list,
                          tag_dict[args.tag])
    print("Done")

if __name__ == '__main__':
    main()
