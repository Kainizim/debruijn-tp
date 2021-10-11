#!/bin/env python3
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

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Hippolyte MIZINIAK"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"  
__maintainer__ = "Hippolyte MIZINIAK"
__email__ = "hippolyte.miziniak@gmail.com"
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
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    i = 0
    with open(fastq_file) as filin:
        for line in filin:
            if i%4 == 1:
                yield line[:-1]
            i = i+1


def cut_kmer(read, kmer_size):
    for i in range(len(read) - kmer_size + 1):
        yield read[i: i + kmer_size]



def build_kmer_dict(fastq_file, kmer_size):
    k_mer_dict = {}
    generator_reads = read_fastq(fastq_file)
    for read in generator_reads:
        generator_kmers = cut_kmer(read, kmer_size)
        for k_mer in generator_kmers:
            if k_mer not in k_mer_dict.keys():
                k_mer_dict[k_mer] = 1
            else:
                k_mer_dict[k_mer] += 1
    print(len(k_mer_dict))
    return k_mer_dict


def build_graph(kmer_dict):
    G = nx.DiGraph()
    for key, value in kmer_dict.items():
        G.add_edge(key[:-1], key[1:], weight = value)
    print(G) 
    return G


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_node_from(path)
        elif not delete_entry_node and delete_sink_node:
            graph.remove_node_from(path[1:])
        elif delete_entry_node and not delete_sink_node:
            graph.remove_node_from(path[:-1])
        else:
            graph.remove_node_from(path[1:-1])
    return graph

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    weight_total = 0
    for  _,_, weight in graph.subgraph(path).edges(data=True):
        print(weight)
        weight_total += weight["weight"]

    return weight_total / (len(path) - 1)


def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    list_start = []
    for node in graph:
        list_it = list(graph.predecessors(node))
        if len(list_it) == 0:
            list_start.append(node)
    return list_start

def get_sink_nodes(graph):
    list_end = []
    for node in graph:
        list_it = list(graph.successors(node))
        if len(list_it) == 0:
            list_end.append(node)
    return list_end

def get_contigs(graph, starting_nodes, ending_nodes):
    list_contigs = []
    seq = ""
    for start in starting_nodes:
        for end in ending_nodes:
            if nx.has_path(graph, start, end):
                for path in nx.all_simple_paths(graph, start, end):
                    for elem in path:
                        seq += elem[0]
                    seq += elem[1:]
                    list_contigs.append((seq, len(seq)))
                    seq = ""
    print(list_contigs)
    return list_contigs

def save_contigs(contigs_list, output_file):
    with open("contigs_file.fasta", "w") as filout:
        for contings in contigs_list:
            a, b = contings
            filout.write(a)
            filout.write(str(b))
            fill("contigs_file.fasta") 
    return contigs_list


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    #fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    #plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    kmer = build_kmer_dict(args.fastq_file, args.kmer_size)

    graph = build_graph(kmer)

    iterator = get_starting_nodes(graph)
    iterator2 = get_sink_nodes(graph)

    starting_nodes = iterator
    ending_nodes = iterator2
    
    list_contigs = get_contigs(graph, starting_nodes, ending_nodes)

    path = save_contigs(list_contigs, "assemblage.txt")

    path_average_weight(graph, path)
    
    graph_1 = nx.DiGraph()
    graph_1.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    graph_1 = select_best_path(graph_1, [[1,2], [3,2]], [1, 1], [5, 10], delete_entry_node=True)
    
    remove_paths(graph, path, delete_entry_node, delete_sink_node)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    
    main()
