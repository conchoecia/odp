#!/usr/bin/env python
"""
This reads in the edge stats and node stats and plots the phylogenetic tree.
Does an easy algorithm and sorts by the lineage.
"""

import argparse
import matplotlib.pyplot as plt
import ete3
import os
import pandas as pd
import sys

# import from other program
from Newick_to_common_ancestors import TaxIDtree
from odp_plotting_functions import format_matplotlib


def parse_args():
    """
    The args that we need are:
      - Node stats file
      - Edge stats file
      - Per_species presence_fusion file - this is used to make a pixel plot

    We just check that these exist. Currently there are no rules that we need to define for plotting.

    Returns:
        args: The parsed arguments
    """
    parser = argparse.ArgumentParser(description='Plot the phylogenetic tree')
    parser.add_argument('-n', '--node_stats', type=str, required=True, help='The node stats file')
    parser.add_argument('-e', '--edge_stats', type=str, required=True, help='The edge stats file')
    parser.add_argument('-p', '--presence_fusion', type=str, required=True, help='The presence fusion file')
    args = parser.parse_args()

    if not os.path.exists(args.node_stats):
        sys.exit('Node stats file does not exist: {}'.format(args.node_stats))

    if not os.path.exists(args.edge_stats):
        sys.exit('Edge stats file does not exist: {}'.format(args.edge_stats))

    return args

def main():
    args = parse_args()

    perspchromdf = pd.read_csv(args.presence_fusion, sep='\t')

    edgedf = pd.read_csv(args.edge_stats, sep='\t')
    nodedf = pd.read_csv(args.node_stats, sep='\t')

    #print("fusions: max={}, min={}".format(edgedf["num_fusions_per_my_this_branch"].max(), edgedf["num_fusions_per_my_this_branch"].min()))
    #print("losses: max={}, min={}".format(edgedf["num_losses_per_my_this_branch"].max(), edgedf["num_losses_per_my_this_branch"].min()))

    # make a tree
    tree = TaxIDtree()
    # update the tree with the nodes
    tree.ingest_node_edge(args.node_stats, args.edge_stats)
    print("node0 is ")
    print(tree.nodes[1])
    sys.exit()

    # just for QC, print the first 10 noes of the sorted nodes and their lineages
    for i in range(len(tree.leaf_order)):
        node = tree.leaf_order[i]
        print(node, tree.nodes[node].lineage)
        if i == 5:
            break

    format_matplotlib()
    # make a simple matplotlib figure with one axis onto which we can plot
    # make 4 panels. Top is a standard tree
    # Second is plotted with the number of fusions
    # third is plotted with the number of losses
    # 4th is a composite
    # set the dimensions for 20in per plot
    fig, ax = plt.subplots(4, 1, figsize=(20, 80))
    # plot the tree
    ax[0] = tree.plot_tree(ax[0], sort = "ascending", text_older_than = 50)
    # set max y to 800
    ax[0].set_ylim(0, 800)
    # fusions
    ax[1] = tree.plot_tree(ax[1], variable = "num_fusions_per_my_this_branch", sort = "ascending")
    ax[1].set_ylim(0, 800)
    # losses
    ax[2] = tree.plot_tree(ax[2], variable = "num_losses_per_my_this_branch", sort = "ascending")
    ax[2].set_ylim(0, 800)
    # save as a pdf
    plt.savefig('tree.pdf')

    # get the leaf order of the tree
    node_order = tree.leaf_order
    from plot_ALG_fusions_v3 import standard_plot_out
    standard_plot_out(perspchromdf, "matches_tree", taxid_order = node_order)

if __name__ == '__main__':
    main()