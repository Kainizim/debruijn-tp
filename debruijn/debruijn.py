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
__copyright__ = "Universite de Paris"
__credits__ = ["Hippolyte MIZINIAK"]
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


def read_fastq(fastq_f):
    """ Identification des k-mer uniques
    :Parameters:
         fastq_file : Path to the file
    """
    with open(fastq_f, "rt") as fillin:
        for line in fillin:
            yield next(fillin).replace("\n", "")
            next(fillin)
            next(fillin)


def cut_kmer(read, kmer_size):
    """Generator cuting kmer contained in the sequence.
      :Parameters:
         read : sequence
         kmer_size : size of the kmer
    """
    i = 0
    count = len(read)
    size = kmer_size 
    while(size <= count):
        kmer = read[i:i+kmer_size]
        i += 1
        size = i + kmer_size
        yield(kmer)


def build_kmer_dict(fastq_file, kmer_size):
    """Create a kmer dictionnary based on the sequences of a file with the specified size.
      :Parameters:
         fastq_file : Path of the file
         kmer_size : size of the kmer
    """
    kmers_list = []
    for line in read_fastq(fastq_file):
        for kmers in cut_kmer(line, kmer_size):
            kmers_list.append(kmers)
    dico = {kmer_key:kmers_list.count(kmer_key) for kmer_key in kmers_list}
    return dico


def build_graph(kmer_dict):
    """Return the corresponding oriented graph of a kmer dictionnary.
      :Parameters:
         kmer_dict : kmer dictionnary
    """
    graph = nx.DiGraph()
    for kmer_key in kmer_dict.keys():
        start = kmer_key[:-1]
        end = kmer_key[1:]
        graph.add_edge(start, end, weight= kmer_dict[kmer_key])
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """qui prend un graphe et une liste de chemin,
    la variable booléenne delete_entry_node pour indiquer si les noeuds d’entrée
    seront supprimés et
    la variable booléenne delete_sink_node pour indiquer si les noeuds de sortie
    seront supprimés et retourne un graphe nettoyé des chemins indésirables."""
    for path in path_list:
        for node in path:
            if path.index(node) == len(path)-1:
                if delete_sink_node:
                    graph.remove_node(node)
                else:
                    continue
            elif path.index(node) == 0:
                if delete_entry_node:
                    graph.remove_node(node)
                else:
                    continue
            else:
                graph.remove_node(node)
    return graph


def std(data):
    """Function which takes a list of values, which returns the standard deviation
    """
    ecart_type = statistics.stdev(data)
    
    return ecart_type


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """Function which takes a graph, a list of paths, a list giving the length
    for each path, a list giving the average weight of each path,
    delete_entry_node to indicate if input nodes will be deleted
    and delete_sink_node to indicate if the output nodes will be deleted
     and returns a graph cleaned of unwanted paths.
    By default delete_entry_node and delete_sink_node here will be False
    """
    id_weight = [path for path, weight in enumerate(weight_avg_list) if weight == max(weight_avg_list)]
    paths_fin = []
    if len(weight_avg_list) > 1:
        paths_fin = [path_list[i] for i in id_weight]
    if len(path_length) > 1:
        path_length = [path_length[i] for i in id_weight]
        idx_length = [path for path, length in enumerate(path_length) if length == max(path_length)]
        paths_fin = [paths_fin[i] for i in idx_length]
    if len(paths_fin) > 1:
        id_random = random.randint(0, len(paths_fin))
        paths_fin = paths_fin[id_random]
    path_list.remove(paths_fin[0])
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

    return graph


def path_average_weight(graph, path):
    """Return the avergae weight of a path.
      :Parameters:
         graph : graph
         path : path within the graph
    """
    somme = sum([graph[path[i]][path[i+1]]['weight'] for i in range(len(path)-1)])

    return somme/(len(path)-1)


def solve_bubble(graph, ancestor_node, descendant_node):
    """Function which takes a graph, an ancestor node, a descendant node and
    returns a clean graph of the bubble located between these two nodes in
    using previously developed functions
    """
    paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    length, weights_average = ([], [])
    for path in paths:
        length.append(len(path))
        weights_average.append(path_average_weight(graph, path))
    graph = select_best_path(graph, paths, length, weights_average)

    return graph


def simplify_bubbles(graph):
    """Function that takes a graph and returns a graph without bubble
    """
    ancestor, descendant = ([], [])
    for node in graph.nodes:
        ans = [succ for succ in graph.successors(node)]
        des = [desc for desc in graph.predecessors(node)]
        if len(ans) >= 2:
            ancestor.append(node)
        if len(des) >= 2:
            descendant.append(node)
    for a, d in zip(ancestor, descendant):
        graph = solve_bubble(graph, a, d)

    return graph


def solve_entry_tips(graph, starting_nodes):
    """Function that takes a graph and a list of input nodes and returns a graph 
    with no unwanted input path
    """
    graph = simplify_bubbles(graph)
    path_list, path_length, path_avg_weight = ([], [], [])
    for node in starting_nodes:
         descs = list(nx.descendants(graph, node))
         for d in descs:
             preds = list(graph.predecessors(d))
             if len(preds) >= 2:
                 paths = list(nx.all_simple_paths(graph, node, d))
                 for p in paths:
                     path_list.append(p)
                     path_length.append(len(p))
                     path_avg_weight.append(path_average_weight(graph, p))
    graph = select_best_path(graph, path_list, path_length, path_avg_weight,
                            delete_entry_node=True, delete_sink_node=False)

    return graph


def solve_out_tips(graph, ending_nodes):
    """Function which takes a graph and a list of output nodes and returns graph without
    unwanted exit path
    """
    graph = simplify_bubbles(graph)
    path_list, path_avg_weight, path_length = ([], [], [])
    final_b = []
    for node in ending_nodes:
        all_nodes = list(graph.nodes)
        for i in range(len(all_nodes)):
            succs = list(graph.successors(all_nodes[i]))
            if len(succs) > 1:
                s = all_nodes[i]
                final_b.append([s, node])
    for b in final_b:
        for path in nx.all_simple_paths(graph, source=b[0], target=b[1]):
            path_list.append(path)
            path_avg_weight.append(path_average_weight(graph, path))
            path_length.append(len(path))
    graph = select_best_path(graph, path_list, path_length, path_avg_weight,
                            delete_entry_node=False, delete_sink_node=True)

    return graph


def get_starting_nodes(graph):
    """Return nodes of the graph with no ancestor.
      :Parameters:
         graph : the graph
    """
    list_starting_nodes = [n for n in graph.nodes() if graph.in_degree(n)==0]

    return list_starting_nodes


def get_sink_nodes(graph):
    """Return nodes of the graph with no successors.
      :Parameters:
         graph : the graph
    """
    list_sink_nodes = [n for n in graph.nodes() if not graph.out_degree(n)]

    return list_sink_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    """Return the contig of the graph for the specified starting node and ending node
       or 'FALSE' if theree is no such contig.
      :Parameters:
         graph : the graph
         starting_node : the node that will be used as start position in the graph to
         find contig
         ending : the node that will be used as end position in the graph to find contig
    """
    contigs = []
    all_length = dict(nx.all_pairs_shortest_path_length(graph))
    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            path = ''
            for paths in nx.all_simple_paths(graph, source=starting_node, target=ending_node):
                contig  = [v for k, v in enumerate(paths,1) if k%2!=0]
                for p in contig:
                    path  += p
            contigs.append( (path, all_length[starting_node][ending_node]+2) )

    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    """
    :Parameters:
         contigs_list : the contig list
         output_file : Path of the file
    """
    with open(output_file, "w+") as f:
        count = 0
        for i in range(0, len(contigs_list)):
            fill_ = fill(contigs_list[i][0])
            count += 1
            print(">contig_" + str(count) + "len=" + str(contigs_list[i][1]) + "\n" )
            print(fill_)
            print("\n")


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
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = simplify_bubbles(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    list_contigs = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(list_contigs, args.output_file)
    draw_graph(graph, 'mygraph')
    
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