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

"""Compute compatibility matrix."""

from __future__ import print_function
import os
import sys
import argparse
import subprocess
import multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

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
    parser.add_argument('-i', dest='aligned_multifasta_file', type=isfile, required=True,
                        help='Aligned multifasta file')
    parser.add_argument('-d', dest='tree_file', type=isfile,
                        help='Tree file corresponding to the aligned multifasta')
    parser.add_argument('-p', dest='phylogeny_software', type=str,
                        choices=["FastTree", "IQ-TREE"],
                        default="IQ-TREE", help='Choose which software to '
                        'perform phylogeny (default FastTree)')
    parser.add_argument('-s', dest='segment_size', type=int, default=300,
                        help='Segment size (default 300)')
    parser.add_argument('-si', dest='segment_increment_size', type=int,
                        default=100, help='Segment increment size (default 100)')
    parser.add_argument('-sc', dest='scale_window', type=int, default=1000,
                        help='Scale for heatmap presentation (default 1000)')
    parser.add_argument('-st', dest='step', type=int, default=10,
                        help='Scale for heatmap presentation (default 10)')
    parser.add_argument('-t', dest='temp_dir', type=isdir, required=True,
                        help='Temporary output dir')
    parser.add_argument('-n', dest='num_cpu', type=int, default=mp.cpu_count(),
                        help='Number of cpu (default {0}).'.format(mp.cpu_count()))
    parser.add_argument('-r', dest='support_threshold', type=float, default=0.7,
                        help='Support threshold (default 0.7)')
    parser.add_argument('-o', dest='res_png', type=str, default='result.png',
                         help='Output png file')
    parser.add_argument('-oz', dest='res_zoomed_png', type=str,
                        default='result_zoomed.png',
                        help='Output zoomed png file')
    parser.add_argument('-om', dest='output_matrix', type=str,
                        default='matrix.tsv', help='Output matrix for real data')
    parser.add_argument('-obm', dest='output_binmatrix', type=str,
                        default='binmatrix.tsv', help='Output binary matrix for real data')
    parser.add_argument('-os', dest='sim_png', type=str,
                        default='simulation.png', help='Output png file')
    parser.add_argument('-osm', dest='output_sim_matrix', type=str,
                        default='matrix_simuation.tsv', help='Output matrix for simulation')
    return parser.parse_args()


def alignment_parameters(aligned_multifasta_file, segment_size, segment_increment_size):
    """Check alignment parameters
    """
    seq = ""
    #check_size = False
    count_seq = 0
    try:
        with open(aligned_multifasta_file, "rt") as aligned_multifasta:
            for line in aligned_multifasta:
                if line.startswith(">"):
                    count_seq += 1
                #elif line.startswith(">") and check_size:
                #    break
                elif len(line) > 0 and count_seq == 1:
                    seq += line.replace("\n", "").replace("\r", "")
            size_alignment = len(seq)
            assert(segment_size < size_alignment)
    except IOError:
        sys.exit("Ali_param: Error cannot open {0}".format(aligned_multifasta_file))
    except AssertionError:
        sys.exit("Error segment size is superior to aligment size: {} {}"
                 .format(segment_size, size_alignment))
    if(count_seq < 4):
        sys.exit("The program need at least 4 sequences in the alignment")
    return size_alignment, count_seq


def divide_alignment(size_alignment, segment_size, segment_increment_size):
    """Define window in the alignment
    """
    i = 0
    newi = segment_size
    segment = [[i, newi]]
    while newi < size_alignment:
        newi = newi + segment_increment_size
        i = i + segment_increment_size
        if newi < size_alignment:
            segment.append([i, newi])
        # We keep segment of the right length for the moment
        elif newi > size_alignment:
            segment.append([size_alignment - segment_size, size_alignment])
    return segment


def fill(text, width=80):
    """Split text"""
    return os.linesep.join(text[i:i+width] for i in xrange(0, len(text), width))


def fasta_seq(multifasta_file):
    header = ""
    seq = ""
    try:
        with open(multifasta_file, "rt") as multifasta:
                for line in multifasta:
                    if line.startswith(">"):
                        if len(header) > 0 and len(seq) > 0:
                            yield [header, seq]
                            header = ""
                            seq = ""
                        header = line[1:]
                        if " " in header:
                            header = header.split(" ")[0]
                        header = line[1:].replace("\n", "")
                    elif len(line) > 0:
                        seq += line.replace("\n", "")
                if len(header) > 0 and len(seq) > 0:
                    yield [header, seq]
    except IOError:
       sys.exit("Fasta_seq:Error cannot open {}".format(multifasta_file))


def extract_segment(aligned_multifasta_file, segment, temp_dir, segment_size,
                    support_threshold, phylogeny_software):
    """Extrat segment from the aligned multifasta file
    """
    seg_section = []
    list_run = []
    list_tree = []
    remove_segment = []
    try:
        for i in segment:
            name = "{}{}segment_{}_{}".format(temp_dir, os.sep, i[0], i[1])
            seg_section.append(open(name + ".fasta", "wt"))
            if phylogeny_software == "FastTree":
                list_run.append("FastTree -gtr -gamma -nt {0}.fasta 2> /dev/null "
                    "|gotree collapse support -s {1} -o {0}.tree".format(
                        name, support_threshold))
            else:
                list_run.append("iqtree -redo -m GTR+I+G4 -s {0}.fasta > /dev/null "
                    "&& gotree collapse support -i {0}.fasta.treefile -s {1} -o {0}.tree".format(
                        name, support_threshold))
            list_tree.append("{0}".format(name))
        #list_run.sort()
        for fasta in fasta_seq(aligned_multifasta_file):
            for i,seg in enumerate(segment):
                assert(len(fasta[1][seg[0]: seg[1]]) == segment_size)
                seq = fasta[1][seg[0]: seg[1]]
                ratio_gaps = seq.count('-') / segment_size
                if ratio_gaps > 0.5:
                    remove_segment.append([seg[0], seg[1]])
                seg_section[i].write(">{}\n{}\n".format(fasta[0],
                    fill(seq)))
        for seg in seg_section:
            seg.close()
        if len(remove_segment) > 0:
            # get unique couple
            remove_segment = [list(x) for x in set(tuple(x) for x in remove_segment)]
            for rm_seg in remove_segment:
                name = "{}{}segment_{}_{}".format(temp_dir, os.sep, rm_seg[0], rm_seg[1])
                if phylogeny_software == "FastTree":
                    cmd = ("FastTree -gtr -gamma -nt {0}.fasta 2> /dev/null "
                           "|gotree collapse support -s {1} -o {0}.tree".format(
                            name, support_threshold))
                else:
                    cmd = ("iqtree -redo -m GTR+I+G4 -s {0}.fasta > /dev/null "
                          "&& gotree collapse support -i {0}.fasta.treefile -s {1} -o {0}.tree".format(
                            name, support_threshold))
                list_run.remove(cmd)
                list_tree.remove("{0}".format(name))
                #segment.remove(rm_seg)
    except IOError:
       sys.exit("extract_segment: Error cannot open {0}".format(aligned_multifasta_file))
    except AssertionError:
        print(fasta)
        print(seg)
        sys.exit("One section does not respect the segment size criteria")
    # except ValueError:
    #     print(rm_seg)
    #     print(name)
    #     print("iqtree -redo -m GTR+I+G4 -s {0}.fasta > /dev/null "
    #           "&& gotree collapse support -i {0}.fasta.treefile -s {1} -o {0}.tree".format(
    #             name, support_threshold))
        # print(list_run)
    #     sys.exit("Stuck")
    return list_run, list_tree, segment


def run_command(cmd):
    """Run command
      Arguments:
          cmd: Command to run
    """
    try:
        # Execution de la ligne de commande
        retcode = subprocess.call(cmd, shell=True)
        # Cas aucun retour du soft
        if retcode == None:
            sys.exit("Child was terminated")
    except OSError, e:
        sys.exit("Execution failed: {0}".format(e))
    except:
        sys.exit("There is something wrong with the command: {0}".format(cmd))


def run_compair(list_run, num_cpus):
    """Run all calculation in parallel
    """
    pool = mp.Pool(processes=num_cpus)
    asyncResult = pool.map_async(run_command, list_run)
    asyncResult.get()
    # for  i in list_run:
    #     run_command(i)


def compute_distance(list_tree):
    list_run = []
    list_seg_dist =[]
    for i in range(len(list_tree)):
        for j in range(i+1, len(list_tree)):
            list_run.append("{3}/quartet_dist -v {0}.tree {1}.tree |cut -f 6 > {0}_{2}.dist".format(list_tree[i], list_tree[j], os.path.basename(list_tree[j]), os.path.dirname(os.path.abspath(__file__))))
            list_seg_dist += ["{0}_{1}.dist".format(list_tree[i], os.path.basename(list_tree[j]))]
    return list_run, list_seg_dist


def get_interest_segment(posit, segment):
    interest_segment = []
    for i in range(len(segment)):
        #print("posit {} segment_{}_{}".format(posit, segment[i][0], segment[i][1]))
        if posit >= segment[i][0] and posit <= segment[i][1]:
            interest_segment.append(segment[i])
        elif not posit > segment[i][1]:
            break
    return interest_segment


def parse_segment_distance(list_seg_dist, segment, segment_increment_size, size_alignment, temp_dir):
    matrix_dist = np.zeros((len(segment), len(segment)))
    posit = 0
    i = 0
    # identify segment for each step
    while i < len(segment):
        list_interest_segment = get_interest_segment(posit, segment)
        # print("posit {} num_segment {}".format(posit, len(list_interest_segment)))
        # print(list_interest_segment)
        val = np.zeros(len(list_interest_segment))
        j = 0
        for seg in segment:
            val = []
            for seg_int in list_interest_segment:
                # print("Segment d'interet")
                # print(seg_int)
                # print(seg)
                name = "{}segment_{}_{}_segment_{}_{}.dist".format(temp_dir, seg_int[0], seg_int[1], seg[0], seg[1])
                name_oposit = "{}segment_{}_{}_segment_{}_{}.dist".format(temp_dir, seg[0], seg[1], seg_int[0], seg_int[1])
                # print(name)
                if os.path.isfile(name):
                    # dat = open(name, "rt").read().splitlines()[0]
                    # print(dat)
                    val.append(1.0 - float(open(name, "rt").read().splitlines()[0]))
                    # print("posit {} segment_{}_{} file = {}".format(posit, seg_int[0], seg_int[1], name))
                elif os.path.isfile(name_oposit):
                    # dat = open(name_oposit, "rt").read().splitlines()[0]
                    # print(dat)
                    val.append(1.0 - float(open(name_oposit, "rt").read().splitlines()[0]))
                # else:
                #     print("0.0")
                #     val.append(0.0)
            
            if len(val) > 0:
                # pass
                matrix_dist[i][j] = np.mean(val)
            # print("Ecriture dans la matrice")
            # print(i, j)
            # print(val)
            # print(np.mean(val))
            j+=1
        # print(posit)
        posit = posit + segment_increment_size
        i += 1
        # print("next")
    return matrix_dist


def plot_heatmap(matrix_dist, output_file, segment_increment_size,
                 size_alignment, count_seq, scale_window, step, min_value,
                 max_value, vmin=0.0, vmax=1.0):
    """Plot half heatmap
    """
    fig, ax = plt.subplots()
    # im = ax.imshow(matrix_dist)
    mask = np.zeros_like(matrix_dist, dtype=np.bool)
    mask[np.tril_indices_from(mask)] = True
    # print(mask)
    sns_plot = sns.heatmap(matrix_dist, xticklabels=step, yticklabels=step,
                           mask=mask, vmin=vmin, vmax=vmax)
    #square=True,
    sns_plot.invert_yaxis()
    scale_labels = np.arange(0, size_alignment + segment_increment_size, scale_window)
    print("scale labels")
    print(scale_labels)
    sns_plot.set_xticklabels(scale_labels)
    sns_plot.set_yticklabels(scale_labels)
    #sns.axes_style("white")
    #sns_plot.set_yticks(all_posit)
    sns_plot.set_title('Phylogenetic compatibility matrix\n'
                       '(n={}, min value={:.2g}, max value={:.2g})'.format(
                        count_seq, min_value, max_value))
    sns_plot.set_ylabel('Sequence position')
    sns_plot.set_xlabel('Sequence position')
    #fig.add_axes(sns_plot, label='axes1')
    #fig.add_axes(sns_plot, label='axes2')
    figure = sns_plot.get_figure()
    figure.savefig(output_file, dpi=500)
    # We want to show all ticks...
    # all_posit = np.arange(0, 7000, 1000)
    # #ax.set_xticks(all_posit)
    # # ax.set_yticks(all_posit)
    # # ... and label them with the respective list entries
    # ax.set_xticklabels(np.arange(0, 8000, 1000))
    # #ax.set_ylabel(all_posit)

    # ax.set_title("Phylogenetic similarity")
    #fig.tight_layout()
    #plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    #cax = plt.axes(all_posit)
    # plt.colorbar(cax=matrix_dist)
    #plt.show()
    #plt.savefig(output_file, dpi=500)


def run_computation(aligned_multifasta_file, segment_size, segment_increment_size,
    size_alignment, count_seq, scale_window, step, temp_dir, phylogeny_software, support_threshold,
    num_cpu, output_file, output_zoomed_file, matrix_file, binmatrix_file, vmin=0.0):
    # Extract alignment
    segment = divide_alignment(size_alignment, segment_size,
                               segment_increment_size)
    print("Segment:")
    print(segment)
    list_run, list_tree, segment = extract_segment(
        aligned_multifasta_file, segment, temp_dir, segment_size,
        support_threshold, phylogeny_software)
    # Generate and collapse tree depending on support value
    print("Compute trees")
    run_compair(list_run, num_cpu)
    # Prepare distance between every topology
    list_run, list_seg_dist = compute_distance(list_tree)
    # Compute distance between topologies
    print("Compute distance")
    run_compair(list_run, num_cpu)
    # Generate matrix
    print("Generate matrix")
    matrix_dist = parse_segment_distance(list_seg_dist, segment,
        segment_increment_size, size_alignment, temp_dir)
    np.savetxt(matrix_file, matrix_dist, delimiter="\t")
    min_value = matrix_dist.min()
    max_value = matrix_dist.max()
    # Plot matrix
    print("Plot matrix")
    plot_heatmap(matrix_dist, output_file, segment_increment_size,
                 size_alignment, count_seq, scale_window, step, min_value, max_value)
    print("vmin {0}".format(vmin))
    if vmin > 0.0 and vmin < 1.0:
        np.savetxt(binmatrix_file, (matrix_dist >= vmin).astype(int), delimiter="\t")
        plot_heatmap(matrix_dist, output_zoomed_file, segment_increment_size,
                 size_alignment, count_seq, scale_window, step, vmin, max_value, vmin, max_value)
    return max_value


def main():
    """Main program
    """
    args = get_arguments()
    # Get alignment size first
    size_alignment, count_seq = alignment_parameters(
        args.aligned_multifasta_file, args.segment_size,
        args.segment_increment_size)
    print("Size alignment {0} nt".format(size_alignment))
    # Run simulation
    print("Let's go for simulation first")
    output_name = args.temp_dir + os.path.splitext(os.path.basename(args.aligned_multifasta_file))[0]
    if args.tree_file:
        # pass
        run_command("{4}/seq-gen -l {2} -of -m GTR {0} | mafft  --thread {3} "
                    "--maxiterate 1000 --localpair - > {1}.ali".format(
                        args.tree_file, output_name, size_alignment,
                        args.num_cpu,
                        os.path.dirname(os.path.abspath(__file__))))
    else:
        if phylogeny_software == "FastTree":
            run_command("FastTree -gtr -gamma -nt {0}.fasta 2> /dev/null "
                        "> {0}.treefile  && {4}/seq-gen -l {2} -of -m GTR {0}.treefile "
                        "| mafft  --thread {3} --maxiterate 1000 --localpair - > {1}.ali".format(
                            args.aligned_multifasta_file, output_name,
                            size_alignment, args.num_cpu,
                            os.path.dirname(os.path.abspath(__file__))))
        else:
            run_command("iqtree -nt AUTO -m GTR+I+G4 -s {0} -redo > /dev/null "
                        "&& {4}/seq-gen -l {2} -of -m GTR {0}.treefile "
                        "| mafft  --thread {3} --maxiterate 1000 --localpair - > {1}.ali".format(
                            args.aligned_multifasta_file, output_name,
                            size_alignment, args.num_cpu,
                            os.path.dirname(os.path.abspath(__file__))))
    # run computation for simulation data
    temp_dir = args.temp_dir + os.sep + "simulation/"
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    vmin = run_computation(output_name + ".ali", args.segment_size,
        args.segment_increment_size, size_alignment, count_seq, args.scale_window,
        args.step, temp_dir, args.phylogeny_software, args.support_threshold,
        args.num_cpu, args.sim_png, args.res_zoomed_png, args.output_sim_matrix,
        "")

    # run computation for real data
    print("Now real data")
    temp_dir = args.temp_dir + os.sep + "result/"
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    run_computation(args.aligned_multifasta_file, args.segment_size,
        args.segment_increment_size, size_alignment, count_seq, args.scale_window,
        args.step, temp_dir, args.phylogeny_software, args.support_threshold, 
        args.num_cpu, args.res_png, args.res_zoomed_png, args.output_matrix,
        args.output_binmatrix, vmin)

if __name__ == "__main__":
    main()
