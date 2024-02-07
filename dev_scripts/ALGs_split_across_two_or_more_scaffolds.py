#!/usr/bin/env python

"""
Program  : ALGs_split_across_two_or_more_scaffolds.py
Language : python 3
Date     : 2024-02-05
Author   : Darrin T. Schultz
Email    : darrin.schultz@univie.ac.at
Github   : https://github.com/conchoecia/odp
Support  : For issues or questions, please search if the topic has been discussed already
           on github and open a new issue if not: https://github.com/conchoecia/odp/issues
License  : GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007. See the LICENSE file.
Citation : If you use this software for your scientific publication, please cite:
           Schultz, DT; Haddock, SHD; Bredeson, JV; Green, RE; Simakov, O & Rokhsar, DS
           Ancient gene linkages support ctenophores as sister to other animals. Nature (2023).
           https://doi.org/10.1038/s41586-023-05936-6

Description:
  - This program takes in a directory of .rbh files and makes a datastructure of all of the significantly-occurring gene groups on different chromosomes.

Usage instructions:
  - See https://github.com/conchoecia/odp#getting-started
"""

# odp stuff to format the plot
import os
import sys
# ODP-specific imports
thisfile_path = os.path.dirname(os.path.realpath(__file__))
scripts_path = os.path.join(thisfile_path, "../scripts")
sys.path.insert(1, scripts_path)
import odp_plotting_functions as odp_plot
import rbh_tools
import argparse
import pandas as pd

# matplotlib stuff
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

def parse_args():
    """
    The things we need to know are:
      - -d --rbh_directory - the directory of the rbh files that we will look at
      - -m --minsig        - the number of the minimum significance value for the whole_FET column in the rbh files. This is used to filter the rbh files.
      - -a --alg           - the name of the ALG set that we are looking at. This should be something that is the header of the rbh files, for example "BCnSSimakov2022"
    """
    parser = argparse.ArgumentParser(description="This program takes in a directory of .rbh files and makes a datastructure of all of the significantly-occurring gene groups on different chromosomes.")
    parser.add_argument("-d", "--rbh_directory", help="The directory of the rbh files that we will look at")
    parser.add_argument("-m", "--minsig", type = float, default = 0.005, help="The minimum significance value for the whole_FET column in the rbh files. This is used to filter the rbh files.")
    parser.add_argument("-a", "--alg", help="The name of the ALG set that we are looking at. This should be something that is the header of the rbh files, for example 'BCnSSimakov2022'")
    args = parser.parse_args()

    # Check that the directory exists
    if not os.path.exists(args.rbh_directory):
        print("The directory you provided does not exist.")
    return args

def plot_chrom_number_vs_number_ALGs_split(ax, splitsdf, min_splits, inferredchromsize):
    """
    Saves a pdf of the plot of the number of ALGs split across two or more scaffolds vs the number of chromosomes.
    returns an axis
    """
    # we need to go through and find all the samples that have at least min_splits
    # first we groupby the sample
    gb = splitsdf.groupby("sample")
    entries = []
    for name, group in gb:
        # samplename
        samplename = group["sample"].unique()[0]
        # get the number of ALGs that are significant at all
        present_ALG_num = len(group["gene_group"].unique())
        # num_ALGs on at least n scaffolds
        num_ALGs = len([x for x in group["gene_group"].value_counts() if x >= min_splits])
        entries.append({"sample": samplename,
                        "num_chroms": inferredchromsize[samplename],
                        "num_ALGs": present_ALG_num,
                        "num_ALGs_min_splits": num_ALGs})
    # make a dataframe of all of the entries
    df = pd.DataFrame(entries)
    print(df)
    # for now make a simple scatter plot
    ax.scatter(df["num_chroms"], df["num_ALGs_min_splits"], alpha = 0.1, lw = 0)
    ax.set_xlabel("Number of chromosomes")
    ax.set_ylabel(f"Number of ALGs split across {min_splits} or more scaffolds")
    return ax

def plot_chrom_number_vs_number_ALGs_perchrom(ax, splitsdf, inferredchromsize):
    """
    For every genome, plots the number of chromosomes (x) vs the number of ALGs on each chromosome (y)
    For every genome, we will plot every genome.
    """
    # we need to go through and find all the samples that have at least min_splits
    # first we groupby the sample
    gb = splitsdf.groupby("sample")
    entries = []
    for name, group in gb:
        # samplename
        samplename = group["sample"].unique()[0]
        # scafcounts
        scafcounts = group["scaffold"].value_counts()
        for chrom in scafcounts.index:
            entries.append({"sample": samplename,
                            "chrom": chrom,
                            "num_chroms": inferredchromsize[samplename],
                            "ALGs_on_chrom": scafcounts[chrom]})
    # make a dataframe of all of the entries
    df = pd.DataFrame(entries)
    #ax.scatter(df["num_chroms"], df["ALGs_on_chrom"], alpha = 0.01, lw = 0)
    #ax.set_xlabel("Number of chromosomes")
    #ax.set_ylabel(f"Number of ALGs on each chromosome")
    # I don't like this plot.

    # We do a little more processing, instead try plotting the mean number of ALGs on each chromosome
    entries = []
    gb = df.groupby(["sample"])
    for name, group in gb:
        entries.append({"sample": name,
                        "num_chroms": group["num_chroms"].unique()[0],
                        "mean_ALGs_on_chrom": group["ALGs_on_chrom"].mean()})
    df = pd.DataFrame(entries)
    print(df)

    # for now make a simple scatter plot
    ax.scatter(df["num_chroms"], df["mean_ALGs_on_chrom"], alpha = 0.1, lw = 0)
    ax.set_xlabel("Number of chromosomes")
    ax.set_ylabel(f"Mean number of ALGs on each chromosome")
    return ax

def main():
    args = parse_args()

    rbh_files = [os.path.join(args.rbh_directory, f) for f in os.listdir(args.rbh_directory) if f.endswith(".rbh")]
    # for testing purposes, just get the top 100 files
    #rbh_files = rbh_files[100]
    #rbh_files = [x for x in rbh_files if "Lepisosteus" in x]

    # we need something of sample to chromnum
    sample_to_chromnum = {}
    # now we go through the files
    entries = []
    for i in range(len(rbh_files)):
        # make a counter that goes back to the beginning of the line
        print(f"\r  analyzing {i+1}/{len(rbh_files)}", end = "")
        rbhfile = rbh_files[i]
        rbhdf = rbh_tools.parse_rbh(rbhfile)
        splitdf, samplename = rbh_tools.rbhdf_to_alglocdf(rbhdf, args.minsig, args.alg)
        chromnum = rbh_tools.rbh_to_scafnum(rbhdf, samplename)
        sample_to_chromnum[samplename] = chromnum
        entries.append(splitdf)
    print()
    print(f"\r  Done analyzing {len(rbh_files)}/{len(rbh_files)}", end = "")

    # make a dataframe of all of the entries
    splitsdf = pd.concat(entries)

    # make a plot
    # CALL THIS TO GET THE VISUAL STYLE WE NEED
    odp_plot.format_matplotlib()
    fw = 10
    fh = 20

    fig = plt.figure(figsize=(fw, fh))
    axes = []

    #for aligning all the panels
    left1   = 0.6
    left2   = 6.5
    left3   = 12.5
    left4   = 18.5

    axes = []

    # This panel is the number of chromosomes vs the number of changes
    bottom1 = 0.6
    bottom2 = 7
    paneldim  = 5
    # we start with a single panel
    plot_params = [left1    /fw, # left offset
                   bottom1  /fh, # bottom offset
                   paneldim /fw, # width
                   paneldim /fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = plot_chrom_number_vs_number_ALGs_split(axes[-1], splitsdf,
                                                      2, sample_to_chromnum)

    # This panel is the number of chromosomes vs the number of ALGs on each chromosome
    paneldim  = 5
    # we start with a single panel
    plot_params = [left1   /fw, # left offset
                   bottom2 /fh, # bottom offset
                   paneldim/fw, # width
                   paneldim/fh] # height
    axes.append(fig.add_axes(plot_params))
    axes[-1] = plot_chrom_number_vs_number_ALGs_perchrom(axes[-1], splitsdf, sample_to_chromnum)

    outfilename = "ALGs_split_across_two_or_more_scaffolds.pdf"
    fig.savefig(outfilename, bbox_inches="tight")

if __name__ == "__main__":
    main()