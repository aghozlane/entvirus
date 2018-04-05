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
import glob
import gzip
import tempfile
import subprocess

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
            msg = "{0} is a file".format(path)
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
    parser.add_argument('-i', dest='data_dir', type=isdir, required=True,
                        help='Annot virus directory')
    parser.add_argument('-r', dest='raw_reads_dir', type=isdir,
                        help='Raw reads directory.')
    parser.add_argument('-a', dest='vp1_id_file', type=isfile, default=None,
                        help='VP1 database annotation')
    parser.add_argument('-p', dest='p1_id_file', type=isfile, default=None,
                        help='P1 database annotation')
    parser.add_argument('-vp1', dest='vp1_annotation_file', type=isfile,
                        default=None, help='VP1 annotation')
    parser.add_argument('-c', dest='count_matrix_file', type=isfile,
                        default=None, help='Count matrix')
    parser.add_argument('-l', dest='annotated', type=str,
                        default="no", help='Annotated read files (default no)')
    parser.add_argument('-d', dest='identity_threshold', type=float,
                        default=0.0, help='Identity threshold (default 0.0)')
    parser.add_argument('-v', dest='coverage_threshold', type=float,
                        default=0.0, help='Coverage threshold (default 0.0)')
    parser.add_argument('-o', dest='output_file', type=str, required=True,
                        help='Output file')
    return parser.parse_args()


def check_file(screen):
    """Get the list of files in the directory
    """
    list_file = []
    try:
        glob_file = glob.glob(screen)
        for file in glob_file:
            if os.path.getsize(file) > 0:
                list_file.append(file)
        assert(len(list_file) > 0)
    except AssertionError:
        print("No file found in {0}".format(screen), file=sys.stderr)
        list_file = [""]
    return list_file


# Sample contigs contigs_vp1 length_vp1 abundance_in_sample annotation
def parse_fastq(fastq_file):
    """Get length of each read
    """
    seq_len_tab = []
    ext = os.path.splitext(fastq_file)[1]
    temp_out = None
    try:
        if ext == ".gz":
            fastq = gzip.open(fastq_file, "rt")
        elif ext == ".dsrc2":
            temp_out = (tempfile.gettempdir() + os.sep +
                        os.path.basename(fastq_file).replace(".dsrc2", ""))
            cmd="dsrc d {0} {1}".format(fastq_file, temp_out)
            retcode = subprocess.call(cmd, shell=True)
            # Case of no return
            if retcode is None:
                sys.exit("Child was terminated")
            fastq = open(temp_out, "rt")
        else:
            fastq = open(fastq_file, "rt")
        for line in fastq:
            # Get the sequence
            seq_len_tab.append(len(fastq.next()))
            # Pass separator
            fastq.next()
            # Pass quality
            fastq.next()
        fastq.close()
        if temp_out:
            os.remove(temp_out)
    except IOError:
        sys.exit("Error cannot open fastq file : {0}".format(fastq_file))
    except OSError as error:
        sys.exit("Execution failed: {0}".format(error))
    return seq_len_tab


def parse_fasta(fasta_file, tag=None):
    """Parse fasta sequence
    """
    header = ""
    sequence = ""
    seq_len_tab = []
    seq_info = {}
    try:
        with open(fasta_file, "rt") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    if len(header) > 0:
                        len_seq = len(sequence)
                        seq_len_tab.append(len_seq)
                        if tag == 'vp1_contigs':
                            seq_info[header] = [len_seq, sequence]
                        else:
                            seq_info[header] = [len_seq]
                        sequence = ""
                    header = line[1:].replace("\n", "").replace("\r", "")
                else:
                    sequence += line.replace("\n", "").replace("\r", "")
            len_seq = len(sequence)
            if len_seq > 0:
                #if tag == 'vp1_contigs':
                #    seq_info[header] = [len_seq, sequence]
                #else:
                seq_info[header] = [len_seq]
                seq_len_tab.append(len_seq)
            #else:
            #     print("WTF: {0}".format(fasta_file))
    except IOError:
        sys.exit("Error cannot open {0}".format(fasta_file))
    return seq_info, seq_len_tab


def get_size_info(seq_len_tab):
    """Get number, mean length and median_length
    """
    return [len(seq_len_tab), sum(seq_len_tab)/len(seq_len_tab),
            sorted(seq_len_tab)[len(seq_len_tab)//2]]


def get_reads_data(sample_data, list_reads, tag):
    """Get reads info
    """
    for i in xrange(len(list_reads[0])):
        # Get sample name
        name,ext = os.path.splitext(os.path.basename(list_reads[0][i]))
        # second round with gz files
        if ext == ".gz" or ext == ".dsrc2":
            name = os.path.splitext(name)[0]
        name = name.replace("_R1","")
        seq_len_tab_fwd = parse_fastq(list_reads[0][i])
        seq_len_tab_rev = parse_fastq(list_reads[1][i])
        if name in sample_data:
            sample_data[name].update({tag+"_fwd":get_size_info(seq_len_tab_fwd)})
            sample_data[name].update({tag+"_rev":get_size_info(seq_len_tab_rev)})
        else:
            sample_data[name] = {tag+"_fwd":get_size_info(seq_len_tab_fwd)}
            sample_data[name].update({tag+"_rev":get_size_info(seq_len_tab_rev)})
    return sample_data


def get_fasta_data(sample_data, list_file, tag):
    """Get contigs data
    """
    #print("get fasta data: " + tag)
    for sample in list_file:
        #print(sample)
        # Get sample name
        name,ext = os.path.splitext(os.path.basename(sample))
        name = name.replace("_spades", "").replace("_clc", "").replace("_minia", "")
        name = name.replace("_megahit", "").replace("_metacompass", "").replace("_ray", "")
        name = name.replace("_vp1_contigs", "")
        if name in sample_data:
            seq_info, seq_len_tab = parse_fasta(sample, tag)
            sample_data[name].update({tag:seq_info})
        else:
            print("Sample {0} not seen before".format(name), file=sys.stderr)
    return sample_data


def load_id(id_file):
    """Load id file
    """
    id_dict = {}
    try:
        with open(id_file, "rt") as id:
            id_reader = csv.reader(id, delimiter="\t")
            for line in id_reader:
                id_dict[line[0]] = line[1:3]
    except IOError:
        sys.exit("Error cannot open {0}".format(id_file))
    return id_dict


def get_blast(blast_file):
    """Get blast annotation
    """
    blast_dict = {}
    try:
        with open(blast_file, "rt") as blast:
            blast_reader = csv.reader(blast, delimiter="\t")
            for line in blast_reader:
                blast_dict[line[0]]= line[1:2] + [round(float(line[8]), 1),
                                             float(line[3])]
    except IOError:
        sys.exit("Error cannot open {0}".format(blast_file))
    return blast_dict


def associate_vp1(sample_data, blast_vp1_file, vp1_id_dict, tag,
                  identity_threshold, coverage_threshold, type_seq):
    """Associate vp1 contigs and their annotation
    """
    for sample in blast_vp1_file:
        # Get sample name
        name,ext = os.path.splitext(os.path.basename(sample))
        if ext == ".gz":
            name = os.path.splitext(name)[0]
        name = name.replace("_vp1","")
        name = name.replace("_p1","")
        print(name)
        # Get target length
        vp1_dict = get_blast(sample)
        for vp1 in vp1_dict:
            if (round(float(vp1_dict[vp1][1]),1) >= float(identity_threshold) and
                round(100.0 * float(vp1_dict[vp1][2])/float(vp1_id_dict[vp1_dict[vp1][0]][1])) >= float(coverage_threshold)):
                # Gives id and coverage against vp1 db
                #print(vp1_dict)
                #print(vp1_id_dict[vp1_dict[vp1][0]])
                sample_data[name][tag][vp1] += (
                    [vp1 + "_" + type_seq, vp1_dict[vp1][0]] + [vp1_id_dict[vp1_dict[vp1][0]][0]] + 
                    [vp1_dict[vp1][1]] +
                    [round(vp1_dict[vp1][2]/float(vp1_id_dict[vp1_dict[vp1][0]][1])*100.0,1)])
            elif tag in sample_data[name]:
                if vp1 in sample_data[name][tag]:
                    sample_data[name][tag][vp1] += [""] * 5

    return sample_data

def load_vp1_annotation(vp1_annotation_file, sample_data):
    """Load the taxonomical annotation
    """
    try:
        with open(vp1_annotation_file, "rt") as vp1_annotation:
            vp1_annotation_reader = csv.reader(vp1_annotation, delimiter="\t")
            for line in vp1_annotation_reader:
                if line[0] in sample_data[line[0].split("|")[0]]["vp1_contigs"]:
                    sample_data[line[0].split("|")[0]]["vp1_contigs"][line[0]] += [line[1], ",".join(line[2:-4]), line[-4], line[-3]]
    except IOError:
        sys.exit("Error cannot open {0}".format(vp1_annotation_file))
    return sample_data


def get_abundance(sample_data, count_matrix_file):
    """
    """
    try:
        with open(count_matrix_file, "rt") as count_matrix:
            count_matrix_reader = csv.reader(count_matrix, delimiter="\t")
            # Pass header
            sample_list = count_matrix_reader.next()[2:]
            total_abundance = [0.0] * len(sample_list)
            # sum column
            for line in count_matrix_reader:
                total_abundance = [total_abundance[i] + float(line[2 + i])
                                   for i in xrange(len(line[2:]))]
        with open(count_matrix_file, "rt") as count_matrix:
            count_matrix_reader = csv.reader(count_matrix, delimiter="\t")
            # Pass header
            sample_list = count_matrix_reader.next()[2:]
            for line in count_matrix_reader:
                abundances = [float(val) for val in line[2:]]
                #print(sample_data[sample_list[0]]['vp1_contigs'][line[0]])
                for i in xrange(len(sample_list)):
                    if 'vp1_contigs' in  sample_data[sample_list[i]]:
                        if line[0] in sample_data[sample_list[i]]['vp1_contigs']:
                            if total_abundance[i] > 0.0:
                                # Raw abundance
                                sample_data[sample_list[i]]['vp1_contigs'][line[0]] += [abundances[i], abundances[i] / total_abundance[i]]
                            else:
                                sample_data[sample_list[i]]['vp1_contigs'][line[0]] += [0.0]
                            #print(sample_data[sample_list[i]]['vp1_contigs'][line[0]])
    except IOError:
        sys.exit("Error cannot open {0}".format(count_matrix_file))
    return sample_data


def write_result(sample_data, output_file, annotated, count_matrix_file,
                 raw_reads_dir):
    """
    """
    #id_tag = {"SEW":"Sewage", "SUP":"Supernatant", "CON":"Concentrate",
    #          "TOL":"Toliara", "MAH":"Mahajanga", "ANT":"Antananarivo"}
    info = ["Sample"]
    main_info = ["Processed_read_fwd", "Mean_length_proc_fwd",
                 "Processed_read_rev", "Mean_length_proc_rev", "Number_contigs",
                 "Number_of_contigs_with_VP1", "Contigs_with_VP1", "Length_contigs_with_VP1",
                 "VP1_sequences", "Matched_VP1", "Annotation_VP1", "Identity_VP1",
                 "Coverage_VP1", "P1_sequences", "Matched_P1", "Annotation_P1", "Identity_P1",
                 "Coverage_P1", "Map_NCBI", "Annotation_NCBI", "Identity_NCBI",
                 "Coverage_NCBI"]
    if raw_reads_dir:
        main_info = ["Raw_read_fwd", "Mean_length_raw_fwd", "Raw_read_rev",
                     "Mean_length_raw_rev"] + main_info
    if annotated == "yes":
        info += ["ID1", "ID2", "ID3", "ID4"]
    if count_matrix_file:
        main_info += ["Raw abundance", "Relative abundance"]
    
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter="\t")            
            output_writer.writerow(info + main_info)
            #print(sample_data.keys())
            for sample in sample_data:
                if annotated == "yes":
                    tag = [sample] + sample.split("_")[0].split("-")[0:4]
                else:
                    tag = [sample]
                # print(sample)
                # Get distric, group
                if "assembly" in sample_data[sample]:
                    if "vp1_contigs" in sample_data[sample]:
                        if len(sample_data[sample]["vp1_contigs"]) == 0:
                            output_writer.writerow(
                                    tag + sample_data[sample]["raw_fwd"][0:2] +
                                    sample_data[sample]["raw_rev"][0:2] +
                                    sample_data[sample]["proc_fwd"][0:2] +
                                    sample_data[sample]["proc_rev"][0:2] +
                                    [len(sample_data[sample]["assembly"]),
                                     0])
                        # print("vp1_contigs")
                        for vp1 in sample_data[sample]["vp1_contigs"]:
                            #print(sample_data[sample]["vp1_contigs"])
                            #print(vp1)
                            #print(len(vp1))
                            #print(sample_data[sample]["vp1_contigs"][vp1])
                            if "raw_fwd" in sample_data[sample]:
                                output_writer.writerow(
                                    tag + sample_data[sample]["raw_fwd"][0:2] +
                                    sample_data[sample]["raw_rev"][0:2] +
                                    sample_data[sample]["proc_fwd"][0:2] +
                                    sample_data[sample]["proc_rev"][0:2] +
                                    [len(sample_data[sample]["assembly"]),
                                     len(sample_data[sample]["vp1_contigs"]),
                                     vp1] +
                                    sample_data[sample]["vp1_contigs"][vp1])
                            else:
                                #print("nono: " + sample)
                                output_writer.writerow(
                                    tag + sample_data[sample]["proc_fwd"][0:2] +
                                    sample_data[sample]["proc_rev"][0:2] +
                                    [len(sample_data[sample]["assembly"]),
                                     len(sample_data[sample]["vp1_contigs"]),
                                     vp1] +
                                    sample_data[sample]["vp1_contigs"][vp1])
                    else:
                        if "raw_fwd" in sample_data[sample]:
                            output_writer.writerow(
                                tag + sample_data[sample]["raw_fwd"][0:2] +
                                sample_data[sample]["raw_rev"][0:2] +
                                sample_data[sample]["proc_fwd"][0:2] +
                                sample_data[sample]["proc_rev"][0:2] +
                                [len(sample_data[sample]["assembly"]), 0])
                        else:
                            output_writer.writerow(
                                tag + sample_data[sample]["proc_fwd"][0:2] +
                                sample_data[sample]["proc_rev"][0:2] +
                                [len(sample_data[sample]["assembly"]), 0])
                else:
                    if "raw_fwd" in sample_data[sample]:
                        output_writer.writerow(
                                tag + sample_data[sample]["raw_fwd"][0:2] +
                                sample_data[sample]["raw_rev"][0:2] +
                                sample_data[sample]["proc_fwd"][0:2] +
                                sample_data[sample]["proc_rev"][0:2] + [0])
                    else:
                        print("no assembly for sample {0}".format(sample))
                        output_writer.writerow(
                                tag + sample_data[sample]["proc_fwd"][0:2] +
                                sample_data[sample]["proc_rev"][0:2] + [0])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


def main():
    """Main program
    """
    args = get_arguments()
    sample_data = {}
    # Get raw data
    if args.raw_reads_dir:
        print("Raw reads")
        list_r1 = check_file(args.raw_reads_dir + os.sep + "*R1*.f*q*")
        reads_file = ([list_r1] + [[r.replace("R1", "R2") for r in list_r1]])
        sample_data = get_reads_data(sample_data, reads_file, "raw")
    # Get cleaned data
    processed_reads_dir = args.data_dir + os.sep + "cleaned_reads" + os.sep
    if os.path.isdir(processed_reads_dir):
        print("Clean reads")
        list_r1 = check_file(processed_reads_dir + "*R1*.f*q*")
        cleaned_reads_file = ([list_r1] + [[r.replace("R1", "R2") for r in list_r1]])
        sample_data = get_reads_data(sample_data, cleaned_reads_file, "proc")

    #print(sample_data["SEW-SUP-MAD-MAH-CMA-17-004_S74"])
    # Check contigs
    contigs_dir = args.data_dir + os.sep + "assembly" + os.sep
    if os.path.isdir(contigs_dir):
        print("Contigs")
        contigs_file = check_file(contigs_dir + "*.fasta")
        sample_data = get_fasta_data(sample_data, contigs_file, "assembly")
    #print(sample_data["SEW-SUP-MAD-MAH-CMA-17-004_S74"])
    # Check vp1 contigs
    vp1_contigs_dir = args.data_dir + os.sep + "vp1_contigs" + os.sep
    if os.path.isdir(vp1_contigs_dir):
        print("vp1_contigs")
        print(sample_data["EVB-2_S265"])
        vp1_contigs_file = check_file(vp1_contigs_dir + "*.fasta")
        #print(vp1_contigs_file)
        sample_data = get_fasta_data(sample_data, vp1_contigs_file,
                                     "vp1_contigs")
        print("after")
        print(sample_data["EVB-2_S265"])
    
    #print(sample_data["SEW-SUP-MAD-MAH-CMA-17-004_S74"])
    # Check annotation
    blast_dir = args.data_dir + os.sep + "blast" + os.sep
    if os.path.isdir(blast_dir) and args.vp1_id_file:
        print("VP1 Annotation")
        vp1_id_dict = load_id(args.vp1_id_file)
        blast_vp1_file = check_file(blast_dir + "*_vp1.tsv")
        #print(blast_vp1_file)
        print(sample_data["EVB-2_S265"])
        sample_data = associate_vp1(sample_data, blast_vp1_file,
                                    vp1_id_dict, "vp1_contigs",
                                    args.identity_threshold,
                                    args.coverage_threshold,
                                    "vp1")
        print("after")
        print(sample_data["EVB-2_S265"])
    if os.path.isdir(blast_dir) and args.p1_id_file:
        print("P1 Annotation")
        p1_id_dict = load_id(args.p1_id_file)
        blast_p1_file = check_file(blast_dir + "*_p1.tsv")
        #print(blast_vp1_file)
        print(sample_data["EVB-2_S265"])
        sample_data = associate_vp1(sample_data, blast_p1_file,
                                   p1_id_dict, "vp1_contigs",
                                   args.identity_threshold,
                                   args.coverage_threshold,
                                   "p1")
        print(sample_data["EVB-2_S265"])
    #print(sample_data["SEW-SUP-MAD-MAH-CMA-17-004_S74"])
    # Load vp1 annotation
    if args.vp1_annotation_file:
        print("avant")
        print(sample_data["EVB-2_S265"])
        sample_data = load_vp1_annotation(args.vp1_annotation_file, sample_data)
        print("apres")
        print(sample_data["EVB-2_S265"])
    print("vp1_contigs abundance")
    #print(sample_data["SEW-SUP-MAD-MAH-CMA-17-004_S74"])
    # Load vp1 abundance
    if args.count_matrix_file:
        sample_data = get_abundance(sample_data, args.count_matrix_file)
    print("Final")
    #print(sample_data["SEW-SUP-MAD-MAH-CMA-17-004_S74"])
    # Write result
    write_result(sample_data, args.output_file, args.annotated,
                 args.count_matrix_file, args.raw_reads_dir)


if __name__ == "__main__":
    main()
