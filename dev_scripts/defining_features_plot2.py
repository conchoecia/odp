#!/usr/bin/env python

"""
This script takes in a directory that contains .tsv.gz files with information about in-clade vs out-clade measurements.

For now, it outputs a table of:
  - pairs that are absolutely unique to specific clades
  - pairs that are exceptionally stable in specific clades
  - pairs that are well-represented and exceptionally stable in specific clades
"""

import argparse
from PhyloTreeUMAP import algcomboix_file_to_dict
import numpy as np
import os
import pandas as pd
import rbh_tools
import scipy.stats as stats
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
sns.set_theme(style="ticks")

def parse_args():
    """
    The things we need for this are:
      - a directory that contains the .tsv.gz files of stats for clades. The format of the files is:
        {nodename}_{clade}_unique_pair_df.tsv.gz
      - path to alg_rbh file that has the info for the pairs
      - a path to a file that has the combinations of indices to the pairs, to the contributing features for the pairs
    """
    parser = argparse.ArgumentParser(description='Define features for clades')
    parser.add_argument('--clade_stats_dir', type=str, help='Directory containing .tsv.gz files with stats for clades')
    parser.add_argument("--rbh_file", type=str, help="Path to the rbh file", required=True)
    parser.add_argument("--pair_combination_path", type=str, help="Path to the pair combination file", required=True)
    parser.add_argument("--sigma",                type=int, help="For the left-side of the distribution, what sigma cutoff to use", default=2)

    args = parser.parse_args()
    # make sure that the rbh file exists
    if not os.path.exists(args.rbh_file):
        raise ValueError(f"{args.rbh_file} does not exist")
    # check that the pair combination file exists
    if not os.path.exists(args.pair_combination_path):
        raise ValueError(f"{args.pair_combination_path} does not exist")
    return parser.parse_args()

def make_marginal_plot(x, y, xlabel, ylabel, max_points = 1000, fontsize = 6,
                       xrange = [], yrange = [], vertical_lines = [], vertical_line_labels = []):
    """
    Just a lil plot to make a marginal plot of an x and y
    """
    # Start with a square Figure.
    fig = plt.figure(figsize=(10, 10))

    # make a title with the xlabel and ylabel on the top of the title
    figtitle = f"{xlabel} (x) vs {ylabel} (y)"
    fig.suptitle(figtitle, fontsize=fontsize+2)
    # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
    # the size of the marginal Axes and the main Axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(3, 3,  width_ratios=(4, 1, 4), height_ratios=(4, 1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.15, hspace=0.15)
    # Create the Axes.
    ax = fig.add_subplot(gs[2, 0])
    ax_histx = fig.add_subplot(gs[1, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[2, 1], sharey=ax)
    # Draw the scatter plot and marginals.
    scatter_hist(x, y, xlabel, ylabel, ax, ax_histx, ax_histy,
                 xrange = xrange, yrange = yrange,
                 vertical_lines = vertical_lines,
                 vertical_line_labels = vertical_line_labels,
                 fontsize = fontsize)

    # we also want to plot the theoretical vs the emperical for the x-axis, because it should be Gaussian
    # theoretical ax
    tax = fig.add_subplot(gs[0, 0])
    qq_plot(tax, x, xlabel, max_points = max_points, fontsize = fontsize)

    # make a qq plot for the other variable too, even though it may not be normally distributed
    tax2 = fig.add_subplot(gs[2, 2])
    qq_plot(tax2, y, ylabel, max_points = max_points, fontsize = fontsize)

    return fig

def qq_plot(ax, x, xlabel, max_points = 1000, fontsize = 6):
    """
    This makes a qqplot in an axis
    """
    ax.set_xlabel("Theoretical Quantiles", fontsize=fontsize)
    ax.set_ylabel("Observed Quantiles", fontsize=fontsize)
    ax.set_title(f"Q-Q Plot of {xlabel}", fontsize=fontsize)
    # set the tick label sizes
    ax.tick_params(axis='both', which='major', labelsize=fontsize)

    if np.std(x) == 0:
        pass
    else:
        # Calculate the z-scores for x
        z_scores = (x - np.mean(x)) / np.std(x)

        # Subset the z-scores if the number of points is too large
        if len(z_scores) > max_points:
            z_scores = np.random.choice(z_scores, size=max_points, replace=False)

        # Get theoretical quantiles
        quantiles = np.linspace(0, 1, len(z_scores))
        theoretical_quantiles = stats.norm.ppf(quantiles)

        # Sort z-scores to get observed quantiles
        observed_quantiles = np.sort(z_scores)

        # Plot the theoretical vs observed quantiles
        ax.scatter(theoretical_quantiles, observed_quantiles, alpha=0.5)
        ax.plot(theoretical_quantiles, theoretical_quantiles, color='red', linestyle='--')  # y=x line

def scatter_hist(x, y, xlabel, ylabel, ax, ax_histx, ax_histy,
                 fontsize = 6, xrange =[], yrange=[],
                 vertical_lines = [], vertical_line_labels = []):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # make a 2d histogram. There are too many points
    # make the alpha low
    #ax.scatter(x, y, alpha=0.1)
    # use a white(low) to black(high) colormap
    ax.hexbin(x, y, gridsize=50, cmap='Greys', bins='log')

    # set the y-limit to between 0 and 1
    if len(yrange) == 2:
        ax.set_ylim(yrange[0], yrange[1])
    # add the labels
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)

    # make sure that the vertical lines and the vertical line labels are the same length
    if len(vertical_lines) != len(vertical_line_labels):
        raise ValueError("The vertical lines and the vertical line labels are not the same length")
    # plot vertical lines if they are given
    for i in range(len(vertical_lines)):
        line      = vertical_lines[i]
        linelabel = vertical_line_labels[i]
        ax.axvline(x=line, color='r', linestyle='--')
        # add the label text
        # get the top of the axis for y
        ylim = ax.get_ylim()
        ypos = ylim[1] - 0.05
        ax.text(line, ypos, linelabel, rotation=90, fontsize=fontsize)

    # use ax limits to determine bins for histograms
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Adjust bin width as needed
    binwidthx = (xlim[1] - xlim[0]) / 50  # Example bin width calculation
    binwidthy = (ylim[1] - ylim[0]) / 50  # Adjust number of bins as needed

    binsx = np.arange(xlim[0], xlim[1] + binwidthx, binwidthx)
    binsy = np.arange(ylim[0], ylim[1] + binwidthy, binwidthy)

    ax_histx.hist(x, bins=binsx)
    ax_histy.hist(y, bins=binsy, orientation='horizontal')

    # change the tick label sizes
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax_histx.tick_params(axis='both', which='major', labelsize=fontsize)
    ax_histy.tick_params(axis='both', which='major', labelsize=fontsize)

def main():
    args = parse_args()

    # load in the rbh file as a df
    rbh = rbh_tools.parse_rbh(args.rbh_file)
    rbh_to_ALG = dict(zip(rbh["rbh"], rbh["gene_group"]))

    # ALGcomboix
    # the type of the key is a tuple with two strings, the type of the value is an int
    ALGcomboix = algcomboix_file_to_dict(args.pair_combination_path)
    # print the first 10 items
    print(list(ALGcomboix.items())[:10])
    # the type of the key is an int, the type of the value is a tuple with two strings
    ALGcomboix_reverse = {v:k for k,v in ALGcomboix.items()}
    print(list(ALGcomboix_reverse.items())[:10])

    # make sure that the directory exists
    if not os.path.isdir(args.clade_stats_dir):
        raise ValueError(f'{args.clade_stats_dir} is not a directory')
    # find all of the files in the director that end in unique_pair_df.tsv.gz
    filelist = [x for x in os.listdir(args.clade_stats_dir) if x.endswith('unique_pair_df.tsv.gz')]

    # Let's set up a list to store the dfs
    summary_stats = []
    clade_to_unique_paircount_any = {}
    sd_number = args.sigma
    for file in filelist:
        nodename = file.split('_')[0]
        taxid    = file.split('_')[1]
        print(f"Looking at the node {nodename} and the file {file}")
        filepath = os.path.join(args.clade_stats_dir, file)
        df       = pd.read_csv(filepath, sep='\t')
        df["nodename"]  = nodename
        df["taxid"]     = taxid
        # make sure that the type of "pair", "notna_in", and "notna_out" is an int
        df["pair"] = df["pair"].astype(str)
        df["notna_in"] = df["notna_in"].astype(int)
        df["notna_out"] = df["notna_out"].astype(int)
        # Infer the number of species in the clade using the notna_in and occupancy_in columns.
        # Get the index of the max value of the notna_in column, divide that value by the occupancy_in column
        maxinix = df["notna_in"].idxmax()
        num_genomes_in = round(df["notna_in"][maxinix] / df["occupancy_in"][maxinix])
        # do it for num_genomes_out
        maxoutix = df["notna_out"].idxmax()
        num_genomes_out = round(df["notna_out"][maxoutix] / df["occupancy_out"][maxoutix])
        df["num_genomes_in"]  = num_genomes_in
        df["num_genomes_out"] = num_genomes_out
        print(f"  - The number of species in the clade is {num_genomes_in} and outside is {num_genomes_out}")
        df["unique_to_clade"] = np.where(df["notna_out"] == 0, 1, 0)
        clade_to_unique_paircount_any[nodename] = df["unique_to_clade"].sum()
        print("The unique pairs for this clade are:\n", df[df["unique_to_clade"] == 1])
        # add a 1 to all the values that could cause a return ratio of 0
        df["mean_in"]  = df["mean_in"]  + 1
        df["mean_out"] = df["mean_out"] + 1
        df["sd_in"]    = df["sd_in"]    + 1
        df["sd_out"]   = df["sd_out"]   + 1
        # now that we have some pseudocounts, the logs will not fail
        df["sd_in_out_ratio"]             = df["sd_in"]   / df["sd_out"]
        df["sd_in_out_ratio_log"]         = np.log10(df["sd_in_out_ratio"])
        df["mean_in_out_ratio"]           = df["mean_in"] / df["mean_out"]
        df["mean_in_out_ratio_log"]       = np.log10(df["mean_in_out_ratio"])

        print("  - Plotting the sd_in_out_ratio_log vs occupancy in")
        outjointplot = f"{nodename}_{taxid}_sdratiolog_occupancy.pdf"
        with PdfPages(outjointplot) as pdf:
            xcol = "sd_in_out_ratio_log"
            ycol = "occupancy_in"
            plotdf = df[[xcol, ycol]]
            # drop the nans
            plotdf = plotdf.dropna()
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol, yrange=[0, 1])
            pdf.savefig(fig)
            plt.close(fig)
            # now plot the same thing, but subset the dataframe s.t. occupancy is g.t.e.q. 0.5
            plotdf = plotdf[plotdf["occupancy_in"] >= 0.5]
            sd_cutoff = plotdf["sd_in_out_ratio_log"].mean() - (sd_number * plotdf["sd_in_out_ratio_log"].std())
            print("the sd cutoff is")
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol, yrange=[0.5, 1], vertical_lines=[sd_cutoff], vertical_line_labels=[f"-{sd_number} sigma"])
            pdf.savefig(fig)
            plt.close(fig)

        print("  - Plotting the sd_in_out_ratio_log vs occupancy in")
        outjointplot = f"{nodename}_{taxid}_sdratiolog_occupancy.pdf"
        with PdfPages(outjointplot) as pdf:
            xcol = "sd_in_out_ratio_log"
            ycol = "occupancy_in"
            plotdf = df[[xcol, ycol]]
            # drop the nans
            plotdf = plotdf.dropna()
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol,  yrange=[0, 1])
            pdf.savefig(fig)
            plt.close(fig)
            # now plot the same thing, but subset the dataframe s.t. occupancy is g.t.e.q. 0.5
            plotdf = plotdf[plotdf["occupancy_in"] >= 0.5]
            sd_cutoff = plotdf[xcol].mean() - (sd_number * plotdf[xcol].std())
            sd_high   = plotdf[xcol].mean() + (sd_number * plotdf[xcol].std())
            print(f"the sd cutoff is: {sd_cutoff}")
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol,
                                     yrange=[0.5, 1], vertical_lines=[sd_cutoff, sd_high], vertical_line_labels=[f"-{sd_number} sigma", f"{sd_number} sigma"])
            pdf.savefig(fig)
            plt.close(fig)

        print("  - Plotting the mean_in_out_ratio_log vs the occupancy in")
        outjointplot = f"{nodename}_{taxid}_meanratiolog_occupancy.pdf"
        with PdfPages(outjointplot) as pdf:
            xcol = "mean_in_out_ratio_log"
            ycol = "occupancy_in"
            plotdf = df[[xcol, ycol]]
            # drop the nans
            plotdf = plotdf.dropna()
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol, yrange = [0, 1])
            pdf.savefig(fig)
            plt.close(fig)
            # now plot the same thing, but subset the dataframe s.t. occupancy is g.t.e.q. 0.5
            plotdf = plotdf[plotdf["occupancy_in"] >= 0.5]
            sd_cutoff = plotdf[xcol].mean() - (sd_number * plotdf[xcol].std())
            print(f"the sd cutoff is: {sd_cutoff}")
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol, yrange=[0.5, 1], vertical_lines=[sd_cutoff], vertical_line_labels=[f"-{sd_number} sigma"])
            pdf.savefig(fig)
            plt.close(fig)

        print("  - Plotting the mean_in_out_ratio_log vs sd_in_out_ratio_log")
        outjointplot = f"{nodename}_{taxid}_meanratiolog_sdratiolog.pdf"
        with PdfPages(outjointplot) as pdf:
            xcol = "mean_in_out_ratio_log"
            ycol = "sd_in_out_ratio_log"
            plotdf = df[[xcol, ycol]]
            # drop the nans
            plotdf = plotdf.dropna()
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol)
            pdf.savefig(fig)
            plt.close(fig)

        print("  - Plotting the mean_in vs occupancy_in")
        outjointplot = f"{nodename}_{taxid}_meandist_occupancy.pdf"
        with PdfPages(outjointplot) as pdf:
            xcol = "mean_in"
            ycol = "occupancy_in"
            plotdf = df[[xcol, ycol]]
            # drop the nans
            plotdf = plotdf.dropna()
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol, yrange = [0, 1])
            pdf.savefig(fig)
            plt.close(fig)
            # now plot the same thing, but subset the dataframe s.t. occupancy is g.t.e.q. 0.5
            plotdf = plotdf[plotdf["occupancy_in"] >= 0.5]
            sd_cutoff = plotdf[xcol].mean() - (sd_number * plotdf[xcol].std())
            print(f"the sd cutoff is: {sd_cutoff}")
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol,
                                     yrange=[0.5, 1], vertical_lines=[sd_cutoff], vertical_line_labels=[f"-{sd_number} sigma"])
            pdf.savefig(fig)
            plt.close(fig)

        # get the df down to where occupancy_in is at least 0.5
        df = df[df["occupancy_in"] >= 0.5]
        print("The unique pairs for this clade are:\n", df[df["unique_to_clade"] == 1])
        df["sd_in_out_ratio_log_sigma"]   = (df["sd_in_out_ratio_log"] - df["sd_in_out_ratio_log"].mean()) / df["sd_in_out_ratio_log"].std()
        df["mean_in_out_ratio_log_sigma"] = (df["mean_in_out_ratio_log"] - df["mean_in_out_ratio_log"].mean()) / df["mean_in_out_ratio_log"].std()
        # if the sd_in_out_ratio_log_sigma is less than -2, and the occupancy_in is at least 0.5, then it is stable, append "stable_in_clade" to the pair_type
        df["stable_in_clade"]   = np.where((df["sd_in_out_ratio_log_sigma"]   < -sd_number) & (df["occupancy_in"] >= 0.5), 1, 0)
        df["unstable_in_clade"] = np.where((df["sd_in_out_ratio_log_sigma"]   >  sd_number) & (df["occupancy_in"] >= 0.5), 1, 0)
        df["close_in_clade"]    = np.where((df["mean_in_out_ratio_log_sigma"] < -sd_number) & (df["occupancy_in"] >= 0.5), 1, 0)
        df["distant_in_clade"]  = np.where((df["mean_in_out_ratio_log_sigma"] >  sd_number) & (df["occupancy_in"] >= 0.5), 1, 0)
        columns_of_interest = ["close_in_clade", "distant_in_clade", "stable_in_clade", "unstable_in_clade", "unique_to_clade"]
        # filter the df to only keep the rows where at least one of the clades of interest is 1
        df = df[df[columns_of_interest].sum(axis=1) > 0]
        summary_stats.append(df)

    # make a composite df
    df = pd.concat(summary_stats)
    df["ortholog1"] = df["pair"].apply(lambda x: ALGcomboix_reverse[int(x)][0])
    df["rbh1"]      = df["ortholog1"].apply(lambda x: rbh_to_ALG[x])
    df["ortholog2"] = df["pair"].apply(lambda x: ALGcomboix_reverse[int(x)][1])
    df["rbh2"]      = df["ortholog2"].apply(lambda x: rbh_to_ALG[x])
    # move these columns to the front:
    priority_cols = ["nodename", "taxid", "ortholog1", "rbh1", "ortholog2", "rbh2"] + columns_of_interest
    df = df[priority_cols + [x for x in df.columns if x not in priority_cols]]
    print(df)
    outpdf = "unique_pairs.tsv"
    df.to_csv(outpdf, sep='\t', index=False)

    # for each nodename, taxid, num_genomes_in, num_genomes_out, get the stats of the pairs
    summary = []
    groupcols = ["nodename", "taxid", "num_genomes_in", "num_genomes_out"]
    gb = df.groupby(groupcols)
    # iterate through the groupby object
    for name, group in gb:
        nodename, taxid, num_genomes_in, num_genomes_out = name
        # get the number of unique pairs
        num_unique_pairs_any     = clade_to_unique_paircount_any[nodename]
        num_unique_pairs_highocc = group["unique_to_clade"].sum()
        num_close_pairs  = group["close_in_clade"].sum()
        num_distant_pairs = group["distant_in_clade"].sum()
        num_stable_pairs = group["stable_in_clade"].sum()
        num_unstable_pairs = group["unstable_in_clade"].sum()
        # get the number of pairs
        num_pairs = len(group)
        # get the number of pairs that are close
        entry = {"nodename": nodename,
                 "taxid":    taxid,
                 "num_genomes_in":  num_genomes_in,
                 "num_genomes_out": num_genomes_out,
                 "num_pairs_unique_to_clade":   num_unique_pairs_any,
                 "num_pairs_highocc_all_types":         num_pairs,
                 "num_pairs_unique_to_clade_highocc":   num_unique_pairs_highocc,
                 "num_pairs_close_in_clade_highocc":    num_close_pairs,
                 "num_pairs_distant_in_clade_highocc":  num_distant_pairs,
                 "num_pairs_stable_in_clade_highocc":   num_stable_pairs,
                 "num_pairs_unstable_in_clade_highocc": num_unstable_pairs}
        summary.append(entry)
    # make a df
    summary_stats_df = pd.DataFrame(summary)
    outpath = "summary_stats.tsv"
    summary_stats_df.to_csv(outpath, sep='\t', index=False)

if __name__ == '__main__':
    main()
