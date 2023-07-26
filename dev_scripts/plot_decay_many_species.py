#!/usr/bin/env python

"""
This script takes the ALG decay values of many species and plots the correlation between
ALG size (recovered specifically for that species) and the decay value.

The example data for each species looks like this:

ALG   conserved scattered total
A1a     110       89       199
G       103       80       183
D       100       81       181
H       34        119      153
F       79        56       135
M       40        84       124
C1      72        50       122
I       30        85       115
K       65        49       114
L       69        44       113
N       33        71       104
Ea      67        33       100
J1      24        67       91
B1      54        31       85
P       52        28       80
J2      26        31       57

The observed number of genes from the ALG will not be the same in all of the species we look at, so we must scale each sample individually.
In this dataset we should probably also record the chromosome on which this ALG was found, the size of that chromosome, and the number of genes on that chromosome.
"""

import matplotlib.pyplot as plt
import argparse
import os
import pandas as pd
import sys

# set up argparse method to get the directory of the .tsv files we want to plot
def parse_args():
    """
    Argparser to get a directory of .tsv files to plot the decay of ALGs for many species.
    """
    parser = argparse.ArgumentParser(description="Plot the decay of ALGs for many species.")
    parser.add_argument("-d", "--directory", type=str, required=True, help="The directory containing the .tsv files for each species.")
    args = parser.parse_args()

    # check that the directory actually exists
    if not os.path.isdir(args.directory):
        raise Exception("The directory you provided does not exist.")
    return args

def plot_decay_log(ax, df_dict):
    """
    makes a log plot of the decay values.
    The log is only on the y-axis.
    """

    ## print the df dict
    #print("DFDICT")
    #print(df_dict)
    #print()
    #sys.exit()

    # Extract two columns from the dataframe
    x_column = "total"
    y_column = "fraction_conserved"

    # Loop through each line and plot
    for thisfile in df_dict:
        df = df_dict[thisfile]
        # make a line plot of each line
        ax.plot(df[x_column], df[y_column], label = thisfile,
                 alpha=0.25, color ="black", lw = 0.5)
    
        # Add labels and legend
        ax.set_xlabel("ALG size (total genes)")
        ax.set_ylabel("fraction conserved (log2 scale)")
    
        # Using log scale for x-axis
        # plt.xscale('log')
        ax.set_yscale('log', base=2)
    return ax    

def plot_decay_median(ax, df_dict):
    """
    makes a plot of the decay values, with the values based on difference from median.

    The y-axis is linear
    """
    # Extract two columns from the dataframe
    x_column = "total"
    y_column = "diff_from_median"

    # Loop through each line and plot
    for thisfile in df_dict:
        df = df_dict[thisfile]
        # make a line plot of each line
        ax.plot(df[x_column], df[y_column], label = thisfile,
                 alpha=0.25, color ="black", lw = 1.0)
    
        # Add labels and legend
        ax.set_xlabel('ALG size (total genes)')
        ax.set_ylabel('difference from median of fraction conserved')
    
    ax.set_ylim(-1, 1)

    return ax    

def filelist_to_plot_data_structure(filelist, bins):
    """
    This function takes the list of fiels and returns a data structure that can be used for plotting

    - Filelist is the list of tsv files with the decay results
    - bins is the number of bins to split the data into. The bins are based on the ALGs
    - mindif is the minimum percent conservation that the ALGs must have in order to be retained in this dataset.
      0.9, for example means that more than 90 of the genes must be retained in the ALG on the same chromosome.
    - maxdif is the maximum percent conservation that the ALGs must have in order to be retained in this dataset.
      1.0, for example means that all of the genes must be retained in the ALG on the same chromosome.

    Mindif is exclusive, while maxdif is inclusive. (i.e. mindif <= x < maxdif) or [mindif, maxdif)

    We make an exception if mindif is 0.0, in which case we do not filter the ALGs at all.
    We also make an exception if maxdif is 1.0, in which case we allow including everything up to 1.0.
    """

    # convert the number of bins into a set of min cutoff values scaled between 0 and 1
    min_cutoffs = [x / bins for x in range(bins)]
    
    # The data structure is like so:
    # plot_these = {  (0,    0.33):   {"file1": df, "file2": df},
    #                 (0.33, 0.66):   {},
    #                 (0.66, 1.0):    {}
    #               }
    #  The keys of the dictionary are the min and max cutoffs for the ALG conservation values
    #   in the diff_from_median_abs column. The keys of those dictionaries are the filenames.
    #   The values of those dictionaries are the dataframes that will be plotted.

    plot_these = {}
    for i in range(len(min_cutoffs)-1):
        plot_these[(min_cutoffs[i], min_cutoffs[i+1])] = {}
    plot_these[(min_cutoffs[-1], 1.0)] = {}


    # go through the files and add the plottable info to a data structure
    for thisfile in filelist:
        # get the basename of the file
        basename = os.path.basename(thisfile)
        # load the results
        thisdf = pd.read_csv(thisfile, sep="\t")
        thisdf["fraction_conserved"] = (thisdf["conserved"] / thisdf["total"]) + 0.000000001
        #get the median conservation value for this species 
        median_conservation = thisdf["fraction_conserved"].median()

        thisdf["diff_from_median"]       = thisdf["fraction_conserved"] - median_conservation
        thisdf["diff_from_median_abs"]   = abs(thisdf["diff_from_median"])

        # we make a filtered DF that only keeps the 50% biggest ALGs
        # this is to perform screens on only the bigger ALGs, which behave most consistently
        filtdf = thisdf[thisdf["total"] > thisdf["total"].quantile(0.3333333)]

        #get the median conservation value for the biggest ALGs
        large_ALG_median_conservation = filtdf["fraction_conserved"].median()

        # classify to which group this genome belongs in terms of conservation
        for thismin, thismax in sorted(plot_these.keys(), reverse= True):
            if large_ALG_median_conservation >= thismin and large_ALG_median_conservation <= thismax:
                plot_these[(thismin, thismax)][basename] = thisdf.copy()
                break
    
    # now that we have classified everything, print out which genomes end up in which group
    for thismin, thismax in sorted(plot_these.keys(), reverse= True):
        print (thismin, thismax, plot_these[(thismin, thismax)].keys())
        print()
    return plot_these


def main():
    # parse the arguments
    args = parse_args()
    # First we must get all of the files in the directory that has the .tsv files that we want to plot the decay for.
    # only get the files that end in .tsv
    filelist = [os.path.abspath(os.path.join(args.directory, f))
                for f in os.listdir(args.directory)
                if f.endswith(".tsv")]

    NUMBER_OF_BINS = 5
    plot_these = filelist_to_plot_data_structure(filelist, NUMBER_OF_BINS)
    # use i = 0 to get the highest-value cutoff, -1 is no cutoff
    index_to_bin = [x for x in sorted(plot_these.keys(), reverse=True)]

    # now we plot the results as many lines
    # use matplotlib plt

    # Each column of plotting is a single bin with a minimum and maximum cutoff for the ALG conservation values
    # We must label each column with the bins.
    # COORDINATE SYSTEM:
    #  x-axis: Bins in decreasing value, can access
    fig, axes = plt.subplots(2, NUMBER_OF_BINS, figsize = (5 * NUMBER_OF_BINS, 15))
    fig.suptitle('Horizontally stacked subplots')


    # now we make the individual plots for each cutoff
    for i in range(NUMBER_OF_BINS-1, -1, -1):
        print("We are in bin {}. Bin values are {}".format(i, index_to_bin[i]))
        plotdfs = plot_these[index_to_bin[i]]
        if len(plotdfs) > 0:
            axes[0, i] = plot_decay_log(    axes[0, i], plotdfs)
            axes[1, i] = plot_decay_median( axes[1, i], plotdfs)

    # Save the plot as a jpeg
    plt.savefig('output_plot.jpg', format='jpeg')
    
    # Show the plot (optional, you can comment this line if you don't want to display the plot)
    plt.show()

if __name__ == "__main__":
    main()