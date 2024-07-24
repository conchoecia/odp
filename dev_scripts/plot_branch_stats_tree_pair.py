#!/usr/bin/env python

"""
This plots a phylogenetic tree, plots the distances for a single pair across the phylogenetic tree
"""

import argparse
import matplotlib.pyplot as plt
import ete3
import os
import pandas as pd
# import rbh tools
import rbh_tools
import sys

# import from other program
from Newick_to_common_ancestors import TaxIDtree
from odp_plotting_functions import format_matplotlib

import matplotlib.patches as patches
from   matplotlib.patches import Rectangle


def parse_args():
    """
    The args that we need are:
      - df file - a dataframe that has the lineage information for each node.
      - Node stats file
      - Edge stats file
      - Ortholog pair
      - a directory of RBH files

    We just check that these exist. Currently there are no rules that we need to define for plotting.

    Returns:
        args: The parsed arguments
    """
    parser = argparse.ArgumentParser(description='Plot the phylogenetic tree')
    parser.add_argument('-d', '--df_file',     type=str, required=True, help='The dataframe file')
    parser.add_argument('-n', '--node_stats',  type=str, required=True, help='The node stats file')
    parser.add_argument('-e', '--edge_stats',  type=str, required=True, help='The edge stats file')
    parser.add_argument('-p', '--pairs_file',  type=str, required=True, help='Comma-separated pair of orthologs.')
    parser.add_argument('-m', '--min_species', type=int, default=1000,  help='The minimum number of species to include in the pair')
    parser.add_argument('-a', '--alg_name',    type=str, default="BCnSSimakov2022", help="The name of the ALG set we are using.")
    parser.add_argument('-r', '--rbh_dir',     type=str, required=True, help='The directory of RBH files')
    args = parser.parse_args()

    if not os.path.exists(args.node_stats):
        sys.exit('Node stats file does not exist: {}'.format(args.node_stats))

    if not os.path.exists(args.edge_stats):
        sys.exit('Edge stats file does not exist: {}'.format(args.edge_stats))

    # check that the pairs file exists
    if not os.path.exists(args.pairs_file):
        sys.exit('Pairs file does not exist: {}'.format(args.pairs_file))

    # check that the rbh directory exists
    if not os.path.exists(args.rbh_dir):
        sys.exit('RBH directory does not exist: {}'.format(args.rbh_dir))

    return args

def main():
    args = parse_args()

    # df of genomes:
    genomedf = pd.read_csv(args.df_file, sep='\t')
    genomedf["taxid_list"] = genomedf["taxid_list"].apply(eval)
    # sort by taxid_list
    genomedf.sort_values(by="taxid_list", inplace=True)
    print(genomedf)
    print(genomedf.columns)

    # load up the pairsdf
    pairsdf = pd.read_csv(args.pairs_file, sep="\t")
    # print the columns
    print(pairsdf.columns)
    pairsdf = pairsdf[pairsdf["sd_in_out_ratio_log_sigma"] < -2]
    pairsdf.sort_values(by="sd_in_out_ratio_log_sigma", inplace=True, ascending=True)
    pairsdf = pairsdf[(pairsdf["num_species_in"] >= args.min_species) | ((pairsdf["nodename"] == "Deuterostomia") & ((pairsdf["rbh1"] == "A1a") | (pairsdf["rbh2"] == "A1a") | (pairsdf["rbh1"] == "A1b") | (pairsdf["rbh2"] == "A1b")))]
    #pairsdf = pairsdf[(pairsdf["occupancy_out"] >= 0.3) | (pairsdf["rbh1"] == "B2") | (pairsdf["rbh2"] == "B2") | (pairsdf["rbh1"] == "A1a") | (pairsdf["rbh2"] == "A1a")]
    # drop the rows with NaN
    pairsdf.dropna(inplace=True)
    print(pairsdf)
    # print the first row of pairsdf

    # go through the RBH files and load them into a dictionary
    rbh_dict = {}
    # first we go through the directory that contains the RBH files, and check them a
    rbh_files = [f for f in os.listdir(args.rbh_dir) if f.endswith(".rbh")]
    counter = 0
    rbhlen = len(rbh_files)
    for rbh_file in rbh_files:
        print(f"    Loading the {counter}/{rbhlen} file: {rbh_file}                    ", end="\r")
        thisrbh = rbh_tools.parse_rbh(os.path.join(args.rbh_dir, rbh_file))
        # get the non-ALG name
        samplename = [x for x in thisrbh.columns if ("_scaf" in x) and (args.alg_name not in x)][0].split("_")[0]
        # get the sample string
        rbh_dict[samplename] = thisrbh
        counter += 1

    format_matplotlib()
    # For each pair, go through the genomes that have those pairs on the same chromosome.
    #  Save those pairs to keep, and add those nodes and edges to the tree
    for index, row in pairsdf.iterrows():
        print()
        print(row)
        ortholog1 = row["ortholog1"]
        alg1      = row["rbh1"]
        ortholog2 = row["ortholog2"]
        alg2      = row["rbh2"]
        taxid     = row["taxid"]
        clade     = row["nodename"]
        sigma = row["sd_in_out_ratio_log_sigma"]
        alg1, alg2 = sorted([alg1, alg2])
        outfile = f"{clade}_S{sigma}_{alg1}_{alg2}_{ortholog1}_{ortholog2}.pdf"
        # dont' bother if the file already exists
        if os.path.exists(outfile):
            continue

        in_values  = {}
        out_values = {}
        for sample_i, sample_row in genomedf.iterrows():
            sample = sample_row["sample"]
            if taxid in sample_row["taxid_list"]:
                inclade = True
            else:
                inclade = False
            # get the rbh file
            if sample in rbh_dict:
                tempdf = rbh_dict[sample][(rbh_dict[sample]["rbh"] == ortholog1) | (rbh_dict[sample]["rbh"] == ortholog2)]
                if len(tempdf) == 2:
                    # check if the scaffold is the same
                    if tempdf[f"{sample}_scaf"].values[0] == tempdf[f"{sample}_scaf"].values[1]:
                        # get the absolute distance between the two genes
                        dist = abs(tempdf[f"{sample}_pos"].values[0] - tempdf[f"{sample}_pos"].values[1])
                        if inclade:
                            in_values[sample] = dist
                        else:
                            out_values[sample] = dist


        # make a pdf with the distances
        fig, ax = plt.subplots(1, 1, figsize=(20, 10))
        # first we plot all the in_values.
        # split them so that they are centerered around 0
        y      = []
        x_min  = []
        x_max  = []
        labels = []
        for i, (sample, value) in enumerate(in_values.items()):
            y.append(i)
            x_min.append(-value/2)
            x_max.append(value/2)
            labels.append(sample)
        ax.axhline(y=i, color="red", linewidth=2)
        # now continue with the out_values
        y_last = i + 1
        for i, (sample, value) in enumerate(out_values.items()):
            y.append(y_last+i)
            x_min.append(-value/2)
            x_max.append(value/2)
            labels.append(sample)
        # plot as rectangles, height 1
        for i in range(len(y)):
            bottom_left = (x_min[i], y[i])
            width = x_max[i] - x_min[i]
            height = 1
            rect = Rectangle(bottom_left, width, height, color="black", alpha=1, linewidth = 0)
            ax.add_patch(rect)
        # set the labels
        ax.set_yticks(range(len(y)))
        ax.set_yticklabels(labels, fontsize = 1)
        ax.set_xlabel("Distance between orthologs")
        ax.set_ylabel("Sample")
        # set a title for the plot with the clade name, the sigma value, the orthologs, and the ALGs
        ax.set_title(f"{clade} S{sigma} {ortholog1} {alg1} {ortholog2} {alg2}")

        # add text to the bottom left corner of the nodename
        in_n_genomes  = len(in_values)
        out_n_genomes = len(out_values)
        ax.text(1.1*min(x_min), 0,     f"Genomes in {clade}, n={in_n_genomes}", fontsize=10, ha="left", va="bottom")
        # add test to the top left corner "not in clade"
        ax.text(1.1*min(x_min), y[-1], f"Genomes not in clade, n={out_n_genomes}", fontsize=10, ha="left", va="top")
        # flip the 0 on y axis

        # set the xlim to bet 1.1 * min and 1.1 * max
        ax.set_xlim(1.1 * min(x_min), 1.1 * max(x_max))
        plt.tight_layout()
        # save this as a pdf
        plt.savefig(outfile)
        plt.close()



    #edgedf = pd.read_csv(args.edge_stats, sep='\t')
    #nodedf = pd.read_csv(args.node_stats, sep='\t')

    #print("This is the edgedf")
    #print(edgedf)
    #print("This is the nodedf")
    #print(nodedf)

    #    ortholog1 = row["ortholog1"]
    #    ortholog2 = row["ortholog2"]
    #    rbh_keys_with_pairs = []


    #    # just for QC, print the first 10 noes of the sorted nodes and their lineages
    #    for i in range(len(tree.leaf_order)):
    #        node = tree.leaf_order[i]
    #        print(node, tree.nodes[node].lineage)
    #        if i == 5:
    #            break

    #    # make a simple matplotlib figure with one axis onto which we can plot
    #    # make 4 panels. Top is a standard tree
    #    # Second is plotted with the number of fusions
    #    # third is plotted with the number of losses
    #    # 4th is a composite
    #    # set the dimensions for 20in per plot
    #    fig, ax = plt.subplots(2, 1, figsize=(20, 40))
    #    # plot the tree
    #    ax[0] = tree.plot_tree(ax[0], sort = "ascending", text_older_than = 50)
    #    # set max y to 800
    #    ax[0].set_ylim(0, 800)
    #    # save as a pdf
    #    plt.savefig('tree.pdf')

    #    # get the leaf order of the tree
    #    node_order = tree.leaf_order
    #    print(node_order)

if __name__ == '__main__':
    main()