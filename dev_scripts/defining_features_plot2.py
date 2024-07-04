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

    args = parser.parse_args()
    # make sure that the rbh file exists
    if not os.path.exists(args.rbh_file):
        raise ValueError(f"{args.rbh_file} does not exist")
    # check that the pair combination file exists
    if not os.path.exists(args.pair_combination_path):
        raise ValueError(f"{args.pair_combination_path} does not exist")
    return parser.parse_args()

def make_marginal_plot(x, y, xlabel, ylabel, outpdf,
                       xrange = [], yrange = [], vertical_lines = []):
    """
    Just a lil plot to make a marginal plot of an x and y
    """
    # Start with a square Figure.
    fig = plt.figure(figsize=(6, 6))
    # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
    # the size of the marginal Axes and the main Axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    # Create the Axes.
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    # Draw the scatter plot and marginals.
    scatter_hist(x, y, xlabel, ylabel, ax, ax_histx, ax_histy,
                 xrange = xrange, yrange = yrange, vertical_lines = vertical_lines)
    # return this as a fig
    return fig

def scatter_hist(x, y, xlabel, ylabel, ax, ax_histx, ax_histy,
                 xrange =[], yrange=[], vertical_lines = []):
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
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # plot vertical lines if they are given
    for line in vertical_lines:
        ax.axvline(x=line, color='r', linestyle='--')

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

def main():
    args = parse_args()

    # load in the rbh file as a df
    rbh = rbh_tools.parse_rbh(args.rbh_file)
    rbh_to_ALG = dict(zip(rbh["rbh"], rbh["gene_group"]))

    # ALGcomboix
    # the type of the key is a tuple with two strings, the type of the value is an int
    ALGcomboix = algcomboix_file_to_dict(args.pair_combination_path)
    # the type of the key is an int, the type of the value is a tuple with two strings
    ALGcomboix_reverse = {v:k for k,v in ALGcomboix.items()}

    # make sure that the directory exists
    if not os.path.isdir(args.clade_stats_dir):
        raise ValueError(f'{args.clade_stats_dir} is not a directory')
    # find all of the files in the director that end in unique_pair_df.tsv.gz
    filelist = [x for x in os.listdir(args.clade_stats_dir) if x.endswith('unique_pair_df.tsv.gz')]

    # First, let's make an out table that will contain the information about the pairs
    summary_stats = []
    entries       = []
    sd_number = 3
    for file in filelist:
        nodename = file.split('_')[0]
        taxid    = file.split('_')[1]
        print(f"Looking at the node {nodename} and the file {file}")
        filepath = os.path.join(args.clade_stats_dir, file)
        df       = pd.read_csv(filepath, sep='\t')
        # make sure that the type of "pair", "notna_in", and "notna_out" is an int
        df["pair"] = df["pair"].astype(str)
        df["notna_in"] = df["notna_in"].astype(int)
        df["notna_out"] = df["notna_out"].astype(int)
        # Infer the number of species in the clade using the notna_in and occupancy_in columns.
        # Get the index of the max value of the notna_in column, divide that value by the occupancy_in column
        maxinix = df["notna_in"].idxmax()
        num_species_in = round(df["notna_in"][maxinix] / df["occupancy_in"][maxinix])
        # do it for num_species_out
        maxoutix = df["notna_out"].idxmax()
        num_species_out = round(df["notna_out"][maxoutix] / df["occupancy_out"][maxoutix])
        print(f"  - The number of species in the clade is {num_species_in} and outside is {num_species_out}")
        # first, get the things that are unique to the clade
        tempdf = df[df["notna_out"] == 0]
        summary_entry = {
            "nodename": nodename,
            "taxid": taxid,
            "num_species_in": num_species_in,
            "num_species_out": num_species_out,
            "num_pairs_unique_to_clade": len(tempdf)}
        # The point of this is adding all of the pairs that are unique to this clade.
        for i, row in tempdf.iterrows():
            thispairix = eval(row["pair"])
            thispair   = ALGcomboix_reverse[thispairix]
            ortholog1  = thispair[0]
            rbh1       = rbh_to_ALG[ortholog1]
            ortholog2  = thispair[1]
            rbh2       = rbh_to_ALG[ortholog2]
            entry = {
                "nodename": nodename,
                "taxid": taxid,
                "pair_type": "unique_to_clade",
                "pair": row["pair"],
                # put the new info here
                "ortholog1": ortholog1,
                "rbh1":      rbh1,
                "ortholog2": ortholog2,
                "rbh2":      rbh2,
                "mean_in": row["mean_in"],
                "mean_out": row["mean_out"],
                "sd_in": row["sd_in"],
                "sd_out": row["sd_out"],
                "num_species_in": num_species_in,
                "num_species_out": num_species_out}
            entries.append(entry)
        # add a 1 to all the values that could cause a return ratio of 0
        df["mean_in"] = df["mean_in"] + 1
        df["mean_out"] = df["mean_out"] + 1
        df["sd_in"] = df["sd_in"] + 1
        df["sd_out"] = df["sd_out"] + 1
        # now that we have some pseudocounts, the logs will not fail
        df["sd_in_out_ratio"]       = df["sd_in"]   / df["sd_out"]
        df["sd_in_out_ratio_log"]   = np.log10(df["sd_in_out_ratio"])
        df["mean_in_out_ratio"]     = df["mean_in"] / df["mean_out"]
        df["mean_in_out_ratio_log"] = np.log10(df["mean_in_out_ratio"])

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
            fig = make_marginal_plot(x, y, xcol, ycol, outjointplot, yrange=[0, 1])
            pdf.savefig(fig)
            plt.close(fig)
            # now plot the same thing, but subset the dataframe s.t. occupancy is g.t.e.q. 0.5
            plotdf = plotdf[plotdf["occupancy_in"] >= 0.5]
            sd_cutoff = plotdf["sd_in_out_ratio_log"].mean() - (sd_number * plotdf["sd_in_out_ratio_log"].std())
            print("the sd cutoff is")
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol, outjointplot, yrange=[0.5, 1], vertical_lines=[sd_cutoff])
            pdf.savefig(fig)
            plt.close(fig)

        # figure out why the sd is 0
        df = df.sort_values("sd_in_out_ratio", ascending=True)

        # only get the things that are greater than or equal to  0.5
        tempdf = df[df["occupancy_in"] >= 0.5]
        sd_cutoff = tempdf["sd_in_out_ratio_log"].mean() - (sd_number * tempdf["sd_in_out_ratio_log"].std())
        # get the cutoff for the occupancy, too
        occupancy_cutoff = 0.5
        print("The cutoff for the sd_in_out_ratio_log is", sd_cutoff)
        print("The cutoff for the occupancy is", occupancy_cutoff)
        print("filtering")
        # get the rows where sd_in_out_ratio_log is less than the cutoff
        #   and the occupancy_in is greater than the cutoff
        tempdf = tempdf[tempdf["sd_in_out_ratio_log"] < sd_cutoff]
        print("Clade filtered for sd_in_out_ratio_log and occupancy_in (stable)")
        print(tempdf)
        summary_entry["num_pairs_stable_in_clade"] = len(tempdf)
        for i, row in tempdf.iterrows():
            thispairix = eval(row["pair"])
            thispair   = ALGcomboix_reverse[thispairix]
            ortholog1  = thispair[0]
            rbh1       = rbh_to_ALG[ortholog1]
            ortholog2  = thispair[1]
            rbh2       = rbh_to_ALG[ortholog2]
            entry = {
                "nodename": nodename,
                "taxid": taxid,
                "pair_type": "stable_in_clade",
                "pair": row["pair"],
                # put the new info here
                "ortholog1": ortholog1,
                "rbh1":      rbh1,
                "ortholog2": ortholog2,
                "rbh2":      rbh2,
                "mean_in": row["mean_in"],
                "mean_out": row["mean_out"],
                "sd_in": row["sd_in"],
                "sd_out": row["sd_out"],
                "num_species_in": num_species_in,
                "num_species_out": num_species_out}
            entries.append(entry)

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
            fig = make_marginal_plot(x, y, xcol, ycol, outjointplot, yrange=[0, 1])
            pdf.savefig(fig)
            plt.close(fig)
            # now plot the same thing, but subset the dataframe s.t. occupancy is g.t.e.q. 0.5
            plotdf = plotdf[plotdf["occupancy_in"] >= 0.5]
            sd_cutoff = plotdf[xcol].mean() - (sd_number * plotdf[xcol].std())
            sd_high   = plotdf[xcol].mean() + (sd_number * plotdf[xcol].std())
            print(f"the sd cutoff is: {sd_cutoff}")
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol, outjointplot,
                                     yrange=[0.5, 1], vertical_lines=[sd_cutoff, sd_high])
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
            fig = make_marginal_plot(x, y, xcol, ycol, outjointplot, yrange = [0, 1])
            pdf.savefig(fig)
            plt.close(fig)
            # now plot the same thing, but subset the dataframe s.t. occupancy is g.t.e.q. 0.5
            plotdf = plotdf[plotdf["occupancy_in"] >= 0.5]
            sd_cutoff = plotdf[xcol].mean() - (sd_number * plotdf[xcol].std())
            print(f"the sd cutoff is: {sd_cutoff}")
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol, outjointplot, yrange=[0.5, 1], vertical_lines=[sd_cutoff])
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
            fig = make_marginal_plot(x, y, xcol, ycol, outjointplot)
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
            fig = make_marginal_plot(x, y, xcol, ycol, outjointplot, yrange = [0, 1])
            pdf.savefig(fig)
            plt.close(fig)
            # now plot the same thing, but subset the dataframe s.t. occupancy is g.t.e.q. 0.5
            plotdf = plotdf[plotdf["occupancy_in"] >= 0.5]
            sd_cutoff = plotdf[xcol].mean() - (sd_number * plotdf[xcol].std())
            print(f"the sd cutoff is: {sd_cutoff}")
            x = plotdf[xcol].to_list()
            y = plotdf[ycol].to_list()
            fig = make_marginal_plot(x, y, xcol, ycol, outjointplot,
                                     yrange=[0.5, 1], vertical_lines=[sd_cutoff])
            pdf.savefig(fig)
            plt.close(fig)

        # get the pairs that are especially unstable in the clade
        occupancy_cutoff = 0.5
        tempdf = df[df["occupancy_in"] >= occupancy_cutoff]
        sd_cutoff = tempdf["sd_in_out_ratio_log"].mean() + (sd_number * tempdf["sd_in_out_ratio_log"].std())
        print("filtering")
        # get the rows where sd_in_out_ratio_log is less than the cutoff
        #   and the occupancy_in is greater than the cutoff
        tempdf = tempdf[tempdf["sd_in_out_ratio_log"] > sd_cutoff]
        print("Clade filtered for sd_in_out_ratio_log and occupancy_in (unstable)")
        print(tempdf)
        summary_entry["num_pairs_unstable_in_clade"] = len(tempdf)
        for i, row in tempdf.iterrows():
            thispairix = eval(row["pair"])
            thispair   = ALGcomboix_reverse[thispairix]
            ortholog1  = thispair[0]
            rbh1       = rbh_to_ALG[ortholog1]
            ortholog2  = thispair[1]
            rbh2       = rbh_to_ALG[ortholog2]
            entry = {
                "nodename": nodename,
                "taxid": taxid,
                "pair_type": "unstable_in_clade",
                "pair": row["pair"],
                # put the new info here
                "ortholog1": ortholog1,
                "rbh1":      rbh1,
                "ortholog2": ortholog2,
                "rbh2":      rbh2,
                "mean_in": row["mean_in"],
                "mean_out": row["mean_out"],
                "sd_in": row["sd_in"],
                "sd_out": row["sd_out"],
                "num_species_in": num_species_in,
                "num_species_out": num_species_out}
            entries.append(entry)

        # now find the ones that are close in a lot of the clade
        occupancy_cutoff = 0.5
        tempdf = df[df["occupancy_in"] >= occupancy_cutoff]
        mean_cutoff = tempdf["mean_in"].mean() - (sd_number * tempdf["mean_in"].std())
        print("filtering")
        tempdf = tempdf[tempdf["mean_in"] < mean_cutoff]
        print("Clade filtered for mean_in and occupancy_in (close means)")
        print(tempdf)
        summary_entry["especially close pairs"] = len(tempdf)
        for i, row in tempdf.iterrows():
            thispairix = eval(row["pair"])
            thispair   = ALGcomboix_reverse[thispairix]
            ortholog1  = thispair[0]
            rbh1       = rbh_to_ALG[ortholog1]
            ortholog2  = thispair[1]
            rbh2       = rbh_to_ALG[ortholog2]
            entry = {
                "nodename": nodename,
                "taxid": taxid,
                "pair_type": "close_in_clade",
                "pair": row["pair"],
                # put the new info here
                "ortholog1": ortholog1,
                "rbh1":      rbh1,
                "ortholog2": ortholog2,
                "rbh2":      rbh2,
                "mean_in":  row["mean_in"],
                "mean_out": row["mean_out"],
                "sd_in":    row["sd_in"],
                "sd_out":   row["sd_out"],
                "num_species_in": num_species_in,
                "num_species_out": num_species_out}
            entries.append(entry)
        summary_stats.append(summary_entry)

    # convert the entries to a df
    df = pd.DataFrame(entries)
    # save as a tsv
    outpath = "unique_pairs.tsv"
    df.to_csv(outpath, sep='\t', index=False)
    print(df)
    # now save the summary stats
    summary_stats_df = pd.DataFrame(summary_stats)
    outpath = "summary_stats.tsv"
    summary_stats_df.to_csv(outpath, sep='\t', index=False)

if __name__ == '__main__':
    main()
