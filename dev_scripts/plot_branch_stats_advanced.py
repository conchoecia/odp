#!/usr/bin/env python

"""
The point of this script is to plot the output of the script plot_branch_stats_vs_time.py.
The above script outputs stats about number of fusions and fissions over time, on different edges.
This script will explore those outputs.
"""

import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import numpy as np
import os
import pandas as pd

def parse_args():
    """
    For now, we just need the modified node and modified edge files.
    Parameters:
      #-e --edgefile: The edge file that contains the fission and fusion stats.
      -n --nodefile: The node file that contains the annotations for the nodes
    """
    parser = argparse.ArgumentParser(description="Plot the fission and fusion stats over time for different edges.")
    parser.add_argument("-e", "--edgefile", help="The edge file that contains the fission and fusion stats.")
    parser.add_argument("-n", "--nodefile", help="The node file that contains the annotations for the nodes.")
    args = parser.parse_args()
    # check that the edgefile and nodefile exist
    for arg in [args.edgefile, args.nodefile]:
        if not os.path.exists(arg):
            raise ValueError(f"File {arg} does not exist.")
    return args

def plot_median_chrom_vs_correlations(df, outfilepath):
    """
    This makes four plots. In each plot, it plots the median chromsize vs the correlation of that value.
    The plots are:
      - median chromsize vs Spearman's R for fusion for extinction
      - median chromsize vs Spearman's R for loss for extinction
      - median chromsize vs Spearman's R for fusion for origination
      - median chromsize vs Spearman's R for loss for origination
    """
    if not outfilepath.endswith(".pdf"):
        raise ValueError("The outfilepath must end in .pdf")

    red = "#D22C16"
    blue = "#3054A3"

    df["dotsize"] = np.sqrt(df["dist_crown"])/np.sqrt(max(df["dist_crown"])) * 20
    # create a normalization object for the nodeage property
    norm = Normalize(vmin=0, vmax=800)
    colormap = cm.plasma
    #df["agecolor"] = [colormap(1 - norm(x)) for x in df["nodeage"]]
    df["agecolor"] = [colormap(norm(x)) for x in df["nodeage"]]


    fig, ax = plt.subplots(4, 1, figsize=(8, 12))
    # get the columns we care about. Get the taxid, median chromsize
    target_x = "chromsize_median"
    target_y = "extinction_fusion_spearman_r"
    title    = "Chromosome size vs Spearman's r for Fusion-ExtinctInten"
    color    = blue
    # get the subdf where there aren't missing values in either of these columns
    plotdf = df[[target_x, target_y, "dotsize", "agecolor"]].dropna()
    ax[0].set_title(title)
    ax[0].scatter(plotdf[target_x], plotdf[target_y], color=plotdf["agecolor"],
                  s = plotdf["dotsize"], edgecolors="none", linewidth = 0)

    target_y = "extinction_losses_spearman_r"
    title    = "Chromosome size vs Spearman's r for Loss-ExtinctInten"
    color    = red
    plotdf = df[[target_x, target_y, "dotsize", "agecolor"]].dropna()
    ax[1].set_title(title)
    ax[1].scatter(plotdf[target_x], plotdf[target_y], color=plotdf["agecolor"],
                  s = plotdf["dotsize"], edgecolors="none", linewidth = 0)

    target_y = "origination_fusion_spearman_r"
    title    = "Chromosome size vs Spearman's r for Fusion-OrigInten"
    color    = blue
    plotdf = df[[target_x, target_y, "dotsize", "agecolor"]].dropna()
    ax[2].set_title(title)
    ax[2].scatter(plotdf[target_x], plotdf[target_y], color=plotdf["agecolor"],
                  s = plotdf["dotsize"], edgecolors = "none", linewidth = 0)

    target_y = "origination_losses_spearman_r"
    title    = "Chromosome size vs Spearman's r for Loss-OrigInten"
    color    = red
    plotdf = df[[target_x, target_y, "dotsize", "agecolor"]].dropna()
    ax[3].set_title(title)
    ax[3].scatter(plotdf[target_x], plotdf[target_y], color=plotdf["agecolor"],
                  s = plotdf["dotsize"], edgecolors = "none", linewidth = 0)

    # set the axis values for all the subplots
    for i in range(4):
        ax[i].set_xlabel("Median Chromosome Size")
        ax[i].set_ylabel("Spearman's R")

    # change the alpha
    for i in range(4):
        for j in range(len(ax[i].collections)):
            ax[i].collections[j].set_alpha(0.5)

    # add more space between each plot
    plt.subplots_adjust(hspace = 0.5)

    # add a colorbar to the side of the plot
    sm = cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical',
                        shrink = 0.5, aspect = 10, pad = 0.1)
    cbar.set_label("Node Age")


    #plt.tight_layout()
    plt.savefig(outfilepath)


def plot_median_chrom_vs_clade_evo_rate(df, outfilepath):
    """
    In this plot, we plot the median chromsize vs the clade evolutionary rate.
    The top plot shows the rate of fusions.
    The bottom plot shows the rate of losses.

    To calculate the rate, we take the number of events in that clade divided by the amount of evolutionary time in that clade.
      More explicitely, this is the "fusions_in_this_clade" column divided by the "dist_crown" column for the fusions panel.
      For the losses panel, it is the "losses_in_this_clade" column divided by the "dist_crown" column.
    """
    if not outfilepath.endswith(".pdf"):
        raise ValueError("The outfilepath must end in .pdf")

    df["dotsize"] = np.sqrt(df["dist_crown"])/np.sqrt(max(df["dist_crown"])) * 20
    df["fusion_rate"] = df["fusions_in_this_clade"]/df["dist_crown"]
    df["loss_rate"] = df["losses_in_this_clade"]/df["dist_crown"]
    # sort by largest fusion rate, descending
    df = df.sort_values(by="fusion_rate", ascending=False)
    # remove the rows where dist_crown is 0
    df = df[df["dist_crown"] != 0]
    print(df[["taxid", "fusions_in_this_clade", "losses_in_this_clade", "fusion_rate", "loss_rate", "dist_crown"]])


    norm = Normalize(vmin=0, vmax=800)
    colormap = cm.plasma
    #df["agecolor"] = [colormap(1 - norm(x)) for x in df["nodeage"]]
    df["agecolor"] = [colormap(norm(x)) for x in df["nodeage"]]

    fig, ax = plt.subplots(2, 1, figsize=(5, 8))
    # get the columns we care about. Get the taxid, median chromsize
    target_x = "chromsize_median"
    target_y = "fusion_rate"
    title    = "Chromosome size vs Fusion Rate"
    # get the subdf where there aren't missing values in either of these columns
    plotdf = df[[target_x, target_y, "dotsize", "agecolor"]].dropna()
    # remove the rows with inf
    plotdf = plotdf[plotdf[target_y] != np.inf]
    print(plotdf)
    ax[0].set_title(title)
    ax[0].scatter(plotdf[target_x], np.log2(plotdf[target_y]), color=plotdf["agecolor"],
                  s = plotdf["dotsize"], edgecolors="none", linewidth = 0)
    print(max(plotdf[target_y]))
    #ax[0].set_ylim(0, max(plotdf[target_y]) * 1.1)

    target_x = "chromsize_median"
    target_y = "loss_rate"
    title    = "Chromosome size vs Loss Rate"
    # get the subdf where there aren't missing values in either of these columns
    plotdf = df[[target_x, target_y, "dotsize", "agecolor"]].dropna()
    plotdf = plotdf[plotdf[target_y] != np.inf]
    print(plotdf)
    ax[1].set_title(title)
    ax[1].scatter(plotdf[target_x], np.log2(plotdf[target_y]), color=plotdf["agecolor"],
                  s = plotdf["dotsize"], edgecolors="none", linewidth = 0)
    #ax[1].set_ylim(0, max(plotdf[target_y]) * 1.1)

    # set the axis values for all the subplots
    for i in range(2):
        ax[i].set_xlabel("Median Chromosome Size")
        ax[i].set_ylabel("Number of Events/Million Years")

    # change the alpha
    for i in range(2):
        for j in range(len(ax[i].collections)):
            ax[i].collections[j].set_alpha(0.5)

    # add more space between each plot
    plt.subplots_adjust(hspace = 0.5)

    # add a colorbar to the side of the plot
    sm = cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical',
                        shrink = 0.5, aspect = 10, pad = 0.1)
    cbar.set_label("Node Age")

    #plt.tight_layout()
    plt.savefig(outfilepath)

def plot_all_variables_node(df, outfilepath):
    """
    This plots the following variables against each other.
    The columns we are interested in are:
      - fusion_rate
      - loss_rate
      - dist_crown
      - dist_crown_plus_root
      - chromsize_median
      - chromsize_mean
      - nodeage
      - fusions_in_this_clade
      - losses_in_this_clade
      - extinction_fusion_spearman_r
      - extinction_losses_spearman_r
      - origination_fusion_spearman_r
      - origination_losses_spearman_r
    """
    from pandas.plotting import scatter_matrix
    if not outfilepath.endswith(".pdf"):
        raise ValueError("The outfilepath must end in .pdf")

    norm = Normalize(vmin=0, vmax=800)
    colormap = cm.plasma
    #df["agecolor"] = [colormap(1 - norm(x)) for x in df["nodeage"]]
    df["agecolor"] = [colormap(norm(x)) for x in df["nodeage"]]

    df["dotsize"] = np.sqrt(df["dist_crown"])/np.sqrt(max(df["dist_crown"])) * 40
    df["fusion_rate"] = np.log2(df["fusions_in_this_clade"]/df["dist_crown"])
    df["loss_rate"]   = np.log2(df["losses_in_this_clade"]/df["dist_crown"])
    # sort by largest fusion rate, descending
    df = df.sort_values(by="fusion_rate", ascending=False)
    # remove the rows where dist_crown is 0
    df = df[df["dist_crown"] != 0]
    print(df[["taxid", "fusions_in_this_clade", "losses_in_this_clade", "fusion_rate", "loss_rate", "dist_crown"]])

    scatter_columns = ["fusion_rate", "loss_rate", "dist_crown",
                          "chromsize_median", "chromsize_mean", "nodeage",
                          "fusions_in_this_clade", "losses_in_this_clade",
                          "extinction_fusion_spearman_r", "extinction_losses_spearman_r",
                          "origination_fusion_spearman_r", "origination_losses_spearman_r", "dotsize", "agecolor"]
    plotdf = df[scatter_columns].dropna()
    # plotdf2 is just those that aren't dotsize or agecolor
    plotdf2 = plotdf.drop(columns=["dotsize", "agecolor"])
    axes = scatter_matrix(plotdf2, figsize=(20, 20), diagonal='kde', alpha=0.2,
                          color = plotdf["agecolor"],  s=plotdf["dotsize"], edgecolors="none", linewidth=0)

    # Rotate x-axis labels by 45 degrees
    for ax in axes.flatten():
        ax.xaxis.label.set_rotation(45)
        ax.yaxis.label.set_rotation(0)
        ax.yaxis.label.set_ha('right')

    # rotate the x-axis labels 45 degrees
    plt.xticks(rotation=45)
    # save as pdf
    plt.savefig(outfilepath)

def plot_all_variables_edge(df, outfilepath):
    """
    This plots the following variables against each other.
    The columns we are interested in are:
      - fusion_rate
      - loss_rate
      - dist_crown
      - dist_crown_plus_root
      - chromsize_median
      - chromsize_mean
      - nodeage
      - fusions_in_this_clade
      - losses_in_this_clade
      - extinction_fusion_spearman_r
      - extinction_losses_spearman_r
      - origination_fusion_spearman_r
      - origination_losses_spearman_r
    """
    pass


def main():
    args = parse_args()
    nodedf = pd.read_csv(args.nodefile, sep='\t', header=0, index_col=None)
    print(nodedf)
    for col in nodedf.columns:
        print(col)

    # plot the median chromsize vs the correlations
    outfilepath = "chromsize_vs_correlations.pdf"
    plot_median_chrom_vs_correlations(nodedf, outfilepath)

    # plot the median chromsize vs the clade evolutionary rate
    outfilepath = "chromsize_vs_clade_evo_rate.pdf"
    plot_median_chrom_vs_clade_evo_rate(nodedf, outfilepath)

    # plot all variables against each other
    outfilepath = "all_variables.pdf"
    #plot_all_variables_node(nodedf, outfilepath)

    # now plot the edge information
    edgedf = pd.read_csv(args.edgefile, sep='\t', header=0, index_col=None)
    print("This is the edgedf")
    print(edgedf)
    for col in edgedf.columns:
        print(col)

if __name__ == '__main__':
    main()