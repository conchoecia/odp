#!/usr/bin/env python

"""
This script plots the number of chromosomes in a sample vs the number of changes.
The user must provide to
sys.arvg[1] the path to the 'per_species_ALG_presence_fusions.tsv' file
sys.arvg[2] the path to the 'species_chrom_counts.tsv' file
"""

from perspchrom_df_to_tree import parse_gain_loss_from_perspchrom_df
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
import sys

def panel_chromsize_vs_changes(ax, labels, chromnum, colocs, losses):
    """
    generates a single panel for the plot of the number of chromosomes vs the number of changes
    """

    # On the panel plot the number of fusions (y1) on the positive y axis
    #  and the number of losses (y2) on the negative y axis. They both use x.
    ax.scatter(chromnum, colocs, color="blue", lw = 0, alpha = 0.02)
    ax.scatter(chromnum, losses, color="red",  lw = 0, alpha = 0.02)

    # make a black vertical line between the two points
    for i in range(len(chromnum)):
        ax.plot([chromnum[i], chromnum[i]],
                [colocs[i], losses[i]], color="black", lw=2, alpha=0.01)

    # add a legend of the two types of dots
    # make the dots' alpha 1 so they are visible
    ax.scatter([], [], color="blue", label="fusions", lw = 0, alpha = 1)
    ax.scatter([], [], color="red", label="losses", lw = 0, alpha = 1)
    ax.legend(loc="upper right")

    # make the x-axis go between 0 and 100
    ax.set_xlim(0, 100)
    # make the x-axis label
    ax.set_xlabel("Number of chromosomes")
    ax.set_ylabel("number of changes. colocalizations are blue, losses are negative")
    ax.set_title("Number of chromosomes vs number of changes")
    return ax

def panel_colocalization_vs_losses(ax, labels, chromnum, colocs, losses):
    """
    This panel is an x-y of the number of colocalizations vs the number of losses
    """
    # make a pandas df of colocs and losses
    df = pd.DataFrame({"colocs": colocs, "losses": losses})
    # remove the rows where either colocs or losses is 0
    keepdf = df[(df["colocs"] != 0) & (df["losses"] != 0)].reset_index(drop=True)
    # show the rows where either colocs or losses is 0
    remdf = df[(df["colocs"] == 0) | (df["losses"] == 0)].reset_index(drop=True)

    # On the panel plot the number of fusions (y1) on the positive y axis
    #  and the number of losses (y2) on the negative y axis. They both use x.
    ax.scatter(keepdf["colocs"], keepdf["losses"], color="black", lw = 0, alpha = 0.05)
    ax.scatter(remdf["colocs"],  remdf["losses"],  color="red", lw = 0, alpha = 0.05)

    colocs = keepdf["colocs"]
    losses = keepdf["losses"]

    # Calculate Spearman's rank correlation coefficient
    spearman_corr, p_value = scipy.stats.spearmanr(colocs, losses)

    # Calculate Kendall's tau
    kendall_tau, kp_value = scipy.stats.kendalltau(colocs, losses)

    print(f"Kendall's tau: {kendall_tau}")
    print(f"P-value: {p_value}")

    print(f"Spearman's rank correlation coefficient: {spearman_corr}")
    print(f"P-value: {kp_value}")
    # Put the spearman's rank correlation coefficient on the plot. round to 3 decimal places
    # Can you put it in the lower right.
    ax.text(0.5, 0.2, "Spearman's coeff.: {:.3f}\nP-value: {:e}\n\nKendall's tau coeff: {:.3f}\nP-value: {:e}".format(
        spearman_corr, p_value, kendall_tau, kp_value), transform=ax.transAxes)

    # add axis labels
    ax.set_xlabel("Number of colocalizations")
    ax.set_ylabel("Number of losses")
    # add a legend of the two types of dots
    # make the dots' alpha 1 so they are visible
    ax.scatter([],[], label="cor. counted", color="black", lw = 0, alpha = 1)
    ax.scatter([],[], label="cor. ignored", color="red", lw = 0, alpha = 1)
    ax.legend(loc="lower right")
    ax.set_title("Correlation between number of colocalizations and losses")
    return ax

def panel_rank_coloc_losses(ax, labels, chromnum, colocs, losses):
    """
    Plots the ranks of the number of colocalizations vs the number of losses
    """
    # make a pandas df of colocs and losses
    df = pd.DataFrame({"colocs": colocs, "losses": losses})
    # remove the rows where either colocs or losses is 0
    keepdf = df[(df["colocs"] != 0) & (df["losses"] != 0)].reset_index(drop=True)

    # Plotting ranks
    rank_x = scipy.stats.rankdata(keepdf["colocs"])
    rank_y = scipy.stats.rankdata(keepdf["losses"])

    ax.scatter(rank_x, rank_y)
    ax.set_xlabel("Rank of X")
    ax.set_ylabel("Rank of Y")
    ax.set_title("Rank of number of colocalizations vs rank of number of losses")
    return ax

def plot_chrom_number_vs_changes(changesfilename, chromsizefilename, outfilename):
    """
    saves a pdf of the plot of the number of chromosomes vs the number of changes
    """
    # make a dictionary of the sample id to the number of chromosomes
    df = pd.read_csv(chromsizefilename, sep='\t')
    chromsize_to_chromnumber = dict(zip(df["sample"], df["chromosomes"]))
    del df

    changedf = parse_gain_loss_from_perspchrom_df(pd.read_csv(changesfilename, sep='\t'))

    labels   = []
    chromnum = [] # number of chromosomes
    colocs   = [] # number of fusions
    losses   = [] # number of losses
    for sample in chromsize_to_chromnumber:
        labels.append(sample)
        chromnum.append(chromsize_to_chromnumber[sample])
        subdf = changedf[changedf["samplename"] == sample]
        # y1 is the total size of the lists in the colocalizations column
        colocs.append(sum([len(x) for x in subdf["colocalizations"]]))
        # y2 is the total size of the lists in the losses column
        y2val = sum([len(x) for x in subdf["losses"]])
        losses.append(0 if y2val == 0 else -1*y2val)

    fw = 20
    fh = 32
    fig = plt.figure(figsize=(fw, fh))
    axes = []

    #for aligning all the panels
    left1   = 0.6
    left2   = 6.5
    left3   = 12.5
    left4   = 18.5

    axes = []

    bottom1 = 0.6
    paneldim = 5

    # This panel is the number of chromosomes vs the number of changes
    # we start with a single panel
    plot_params = [left1   /fw, # left offset
                   bottom1 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_chromsize_vs_changes(axes[-1], labels, chromnum, colocs, losses)

    # this panel is x-y of number of colocalizations vs number of losses
    plot_params = [left2   /fw, # left offset
                   bottom1 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_colocalization_vs_losses(axes[-1], labels, chromnum, colocs, losses)

    # add a panel of the rank plot of the numver of colocalizations and the number of losses
    plot_params = [left3   /fw, # left offset
                   bottom1 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = panel_rank_coloc_losses(axes[-1], labels, chromnum, colocs, losses)

    ## add a panel of the rank plot of the numver of colocalizations and the number of losses
    #plot_params = [left4   /fw, # left offset
    #               bottom1 /fh, # bottom offset
    #               paneldim/fw, # width
    #               paneldim/fh] # height
    #axes.append(fig.add_axes(plot_params))
    #axes[-1] = panel_rankmatrix_coloc_losses(axes[-1], labels, chromnum, colocs, losses)



    fig.savefig(outfilename, bbox_inches="tight")

def main():
    plot_chrom_number_vs_changes(sys.argv[1], sys.argv[2], "chrom_number_vs_changes.pdf")

if __name__ == "__main__":
    main()