#!/usr/bin/env python

"""
This script makes plots related to plotting the changes on phylogenetic branches over time.

Takes in a file called "Jun240604.edge_information.tsv" that has the following columns:
    - source: the parent node
    - target: the child node
    - source_age: the age of the parent node
    - target_age: the age of the child node
    - branch_length: the branch length between the parent and child node.
    - source_ages: the age counter of the parent node
    - target_ages: the age counter of the child node

Also takes in the df_stats file that has fusions or losses.
"""

import argparse
import numpy as np
from collections import Counter
import ete3
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.stats import spearmanr,kendalltau # I use scipy for the correlation functions

from odp_plotting_functions import format_matplotlib

def parse_args():
    """
    We need two files:
      - the edge_information file.
      - the statsdf file.
      - the intensity of extinction file.
      - an optional list of taxids to omit, anything in the clade that has one of these taxids as a parent (or as the actual taxid) will be omitted.

    We will check that both files exist.
    """
    parser = argparse.ArgumentParser(description='Plot branch stats vs time')
    parser.add_argument("-e", "--edge_information", type=str, help='The edge information file')
    parser.add_argument("-s", "--statsdf", type=str, help='The df_stats file')
    parser.add_argument("-i", "--intensity_of_extinction", type=str, help='The intensity of extinction file. Must have columns: Time (Ma)       All Genera      Short-Lived     Long-Lived      Well-Resolved   Intensity (%)   Intensity(%)')
    parser.add_argument("-o", "--omit_taxids", type=str, help='A string with a list of taxids to omit, separated by commas.')
    args = parser.parse_args()

    if not os.path.exists(args.edge_information):
        raise ValueError('The edge information file does not exist')

    if not os.path.exists(args.statsdf):
        raise ValueError('The statsdf file does not exist')

    # we can optionally plot the intensity of extinction
    if args.intensity_of_extinction:
        if not os.path.exists(args.intensity_of_extinction):
            raise ValueError('The intensity of extinction file does not exist')

    # parse the omit_taxids into a list
    if args.omit_taxids:
        if type(args.omit_taxids) == str:
            args.omit_taxids = [int(x) for x in args.omit_taxids.split(',')]
    return args

def plot_intensity_of_extinction(count_df, intensity_of_extinction_filepath):
    """
    This makes another pdf of the intensity of extinction plus a df of the counts.
    This dataframe should be extracted from the supplementary information of this paper: https://www.nature.com/articles/nature03339

    Rohde, Robert A., and Richard A. Muller. "Cycles in fossil diversity." Nature 434.7030 (2005): 208-210.
    """
    # read in the intensity of extinction file
    intensity_of_extinction_df = pd.read_csv(intensity_of_extinction_filepath, sep='\t')

    # here we will make 6 plots.
    # Top row is fusions and losses plotted together
    # Second row is fusions
    # Third row is losses
    # For the first column, the y-axis (ratio of changes) is in linear scale
    # For the second column, the y-axis (ratio of changes) is in log2 scale
    # The x-axis will be the intensity of extinction

    # first, we need to get the intensity of extinction for each age.
    mya_to_intensity = {row['Time (Ma)']: row['Intensity (%)'] for i, row in intensity_of_extinction_df.iterrows()}
    # map the ages, but if the ages aren't present, give a value of -1
    count_df["age_positive"] = count_df["age"].apply(lambda x: -1 * x)
    count_df["intensity_of_extinction"] = count_df["age_positive"].map(mya_to_intensity).fillna(-1)

    # only plot the values where intensity of extinction is not -1
    plotdf = count_df[count_df["intensity_of_extinction"] != -1]
    plotdf["intensity_normalized"] = plotdf["intensity_of_extinction"] / 100
    # We want to measure the spearman correlation between the intensity of extinction and the fusions and losses. Use log for the fusions and losses.
    #  - It is not appropriate to use Pearson correlation, because the data are not linear (the percent, ranging between 0-100 or 0-1, cannot be linearly related to the number of fusions or losses).
    fusion_spearman_corr, fusion_spearman_p = spearmanr(plotdf["intensity_normalized"], np.log2(plotdf["fusions_ratio"]))
    loss_spearman_corr, loss_spearman_p     = spearmanr(plotdf["intensity_normalized"], np.log2(plotdf["losses_ratio"]))
    fusion_kendalltau_corr, fusion_kendalltau_p = kendalltau(plotdf["intensity_normalized"], np.log2(plotdf["fusions_ratio"]))
    loss_kendalltau_corr,   loss_kendalltau_p   = kendalltau(plotdf["intensity_normalized"], np.log2(plotdf["losses_ratio"]))
    print(f"The Spearman correlation between the intensity of extinction and the log2 of the fusions ratio is {fusion_spearman_corr}")
    print(f"  - the p-value is {fusion_spearman_p}")
    print(f"The Kendall tau correlation between the intensity of extinction and the log2 of the fusions ratio is {fusion_kendalltau_corr}")
    print(f"  - the p-value is {fusion_kendalltau_p}")
    print(f"The Spearman correlation between the intensity of extinction and the log2 of the losses ratio is {loss_spearman_corr}")
    print(f"  - the p-value is {loss_spearman_p}")
    print(f"The Kendall tau correlation between the intensity of extinction and the log2 of the losses ratio is {loss_kendalltau_corr}")
    print(f"  - the p-value is {loss_kendalltau_p}")

    # now make the plots
    fig, ax = plt.subplots(3, 2, figsize=(10, 10))
    fontsize = 8
    # first, combined
    ax[0][0].set_title("Fusions and Losses/branch vs Intensity of Extinction")
    ax[0][0].scatter(plotdf["intensity_of_extinction"], plotdf["fusions_ratio"], color='blue')
    ax[0][0].scatter(plotdf["intensity_of_extinction"], plotdf["losses_ratio"], color='red')
    ax[0][0].set_xlabel("Intensity of Extinction (%)")
    ax[0][0].set_ylabel("Changes/branch")
    # second, just fusions
    ax[1][0].scatter(plotdf["intensity_of_extinction"], plotdf["fusions_ratio"], color='blue')
    ax[1][0].set_title("Fusions/branch vs Intensity of Extinction")
    ax[1][0].set_xlabel("Intensity of Extinction (%)")
    ax[1][0].set_ylabel("Fusions/branch")
    # third, just losses
    ax[2][0].scatter(plotdf["intensity_of_extinction"], plotdf["losses_ratio"], color='red')
    ax[2][0].set_title("Losses/branch vs Intensity of Extinction")
    ax[2][0].set_xlabel("Intensity of Extinction (%)")
    ax[2][0].set_ylabel("Losses/branch")

    # now make the same thing, but log2 of the ratio
    # first plot, combined on one axis
    ax[0][1].set_title("Fusions and Losses/branch vs Intensity of Extinction")
    ax[0][1].scatter(plotdf["intensity_of_extinction"], np.log2(plotdf["fusions_ratio"]), color='blue')
    ax[0][1].scatter(plotdf["intensity_of_extinction"], np.log2(plotdf["losses_ratio"]), color='red')
    ax[0][1].set_xlabel("Intensity of Extinction (%)")
    ax[0][1].set_ylabel("Changes/branch (log2)")
    # second, just fusions
    ax[1][1].scatter(plotdf["intensity_of_extinction"], np.log2(plotdf["fusions_ratio"]), color='blue')
    ax[1][1].set_title("Fusions/branch (log2) vs Intensity of Extinction")
    ax[1][1].set_xlabel("Intensity of Extinction (%)")
    ax[1][1].set_ylabel("Fusions/branch (log2)")
    # third, just losses
    ax[2][1].scatter(plotdf["intensity_of_extinction"], np.log2(plotdf["losses_ratio"]), color='red')
    ax[2][1].set_title("Losses Ratio/branch (log2) vs Intensity of Extinction")
    ax[2][1].set_xlabel("Intensity of Extinction (%)")
    ax[2][1].set_ylabel("Losses/branch (log2)")

    # change the fontsize of the axes and the titles
    fontsize = 8
    for axi in [0,1,2]:
        for axj in [0,1]:
            ax[axi][axj].tick_params(axis='both', which='major', labelsize=fontsize)
            ax[axi][axj].set_title(ax[axi][axj].get_title(), fontsize=fontsize)

    # force change the dot size
    dotsize = 1
    for axi in [0,1,2]:
        for axj in [0,1]:
            for coll in ax[axi][axj].collections:
                coll.set_sizes([dotsize])

    # change the opacity to 0.25
    opacity = 0.25
    for axi in [0,1,2]:
        for axj in [0,1]:
            for coll in ax[axi][axj].collections:
                coll.set_alpha(opacity)

    # increase the horizontal and vertical space between the panels
    plt.subplots_adjust(hspace=0.5, wspace=0.5)

    # save to a pdf
    outfile = "intensity_of_extinction_vs_branch_stats.pdf"
    plt.savefig(outfile)

def main():
    args = parse_args()

    # read in the dataframes
    edgedf  = pd.read_csv(args.edge_information, sep='\t')
    eventdf = pd.read_csv(args.statsdf, sep='\t')

    # verify that there every value in the target column is unique
    if len(edgedf['target'].unique()) != len(edgedf):
        raise ValueError('The target column is not unique')
    else:
        print("Every value in the target column is unique. This is a legal DAG.")
    # we need a taxid_to_parent dictionary to get the mapping to id the branch length
    # Actually, we can just get the index of this row, because we now know that the target column is unique,
    #   and we will want to work with multiple values from this row.
    taxid_to_parent = {row['target']: row['source'] for i, row in edgedf.iterrows()}
    taxid_to_row_index = {row['target']: i for i, row in edgedf.iterrows()}

    # build a structured vector to hold the counts at each age.
    # get the unique source ages
    # make a datastructure with the counts of the source ages
    columns = {"total_branches_at_this_age": 0,
               "total_losses_at_this_age":   0,
               "total_fusions_at_this_age":  0}
    # get all of the unique categories without the parentheses from the eventdf
    ALGs         = sorted(set([x for x in eventdf['change'].unique() if x[0] != '(']))
    for thisALG in ALGs:
        columns[thisALG] = 0

    # Get all of the possible combinations
    Combinations = sorted(set([x for x in eventdf['change'].unique() if x[0] == '(']))
    for thisCombo in Combinations:
        columns[thisCombo] = 0

    # get the oldest age in source_age in the edgedf
    oldest_age = int(edgedf['source_age'].max())
    # populate the data structure with one row per integer
    count_struct = {i: columns.copy() for i in range(oldest_age + 1)}

    if args.omit_taxids is not None:
        omit_count_structs = {taxid: count_struct.copy() for taxid in args.omit_taxids}

    # Now go through and record the total number of branches at each age. This is done just using the edgedf.
    # This value is what we
    for i, row in edgedf.iterrows():
        for j in range(int(row['target_age']), int(row['source_age'])+1):
            count_struct[j]['total_branches_at_this_age'] += 1

    # set up an NCBI object to get the lineage of each taxid
    NCBI = ete3.NCBITaxa()

    # Now go through the events and record the number of event
    missed_edges = 0
    for i, row in eventdf.iterrows():
        # This is a single event. We need to know the branch on which it occurred.
        event_branch_id = row["target_taxid"]

        if args.omit_taxids is not None:
            # get the NCBI lineage of this branch
            lineage = NCBI.get_lineage(event_branch_id)
            # delete everything up until the event_branch_id
            lineage = lineage[lineage.index(event_branch_id):]

        # get the row info for this branch. We need this info to plot later
        if event_branch_id not in taxid_to_row_index:
            missed_edges += 1
            continue
        else:
            edge_row_index = taxid_to_row_index[event_branch_id]
            branchrow = edgedf.iloc[edge_row_index]
            # determine whether the event was a fusion or a loss
            changetype = None
            if row["change"][0] == "(":
                changetype = "fusions"
            else:
                changetype = "losses"
            # We now must span over the branch distance and record the event at each age.
            for j in range(int(branchrow['target_age']), int(branchrow['source_age']) + 1):
                count_struct[j]['total_' + changetype + '_at_this_age'] += 1
                # increment the count of the specific event
                count_struct[j][row['change']] += 1

    print("There were {} missed edges".format(missed_edges))
    # all of the data should be in the struct now. Turn it into a pandas df
    resultsdf = pd.DataFrame(count_struct).T

    for event in ["fusions", "losses"]:
        resultsdf[f"{event}_ratio"] = resultsdf[f"total_{event}_at_this_age"] / resultsdf["total_branches_at_this_age"]
    resultsdf["age"] = -1 * resultsdf.index

    print(resultsdf)

    # call the function to properly format the text
    format_matplotlib()

    # Now make three plots, each top of one another.
    # The 0th row shows both fusions and fissions ratios together in the same axis.
    # The 1st row shows just the fusion ratio on the axis.
    # The 2nd row shows just the fission ratio on the axis.
    # Use two different panels, plot the fusions in the top panel, and the losses in the bottom panel. Color the fusions blue and the losses red.
    # We also make a second axis which is the log2 of the ratio of fusions to losses.

    fig, ax = plt.subplots(3, 2, figsize=(10, 10))
    # first, combined
    ax[0][0].set_title("Fusions and Losses/branch vs Age")
    ax[0][0].plot(resultsdf["age"], resultsdf["fusions_ratio"], color='blue')
    ax[0][0].plot(resultsdf["age"], resultsdf["losses_ratio"], color='red')
    ax[0][0].set_xlabel("Age")
    ax[0][0].set_ylabel("Changes/branch")
    # second, just fusions
    ax[1][0].plot(resultsdf["age"], resultsdf["fusions_ratio"], color='blue')
    ax[1][0].set_title("Fusions/branch vs Age")
    ax[1][0].set_xlabel("Age")
    ax[1][0].set_ylabel("Fusions/branch")
    # third, just losses
    ax[2][0].plot(resultsdf["age"], resultsdf["losses_ratio"], color='red')
    ax[2][0].set_title("Losses/branch vs Age")
    ax[2][0].set_xlabel("Age")
    ax[2][0].set_ylabel("Losses/branch")

    # now make the same thing, but log2 of the ratio
    # first plot, combined on one axis
    ax[0][1].set_title("Fusions and Losses/branch vs Age")
    ax[0][1].plot(resultsdf["age"], np.log2(resultsdf["fusions_ratio"]), color='blue')
    ax[0][1].plot(resultsdf["age"], np.log2(resultsdf["losses_ratio"]), color='red')
    ax[0][1].set_xlabel("Age")
    ax[0][1].set_ylabel("Changes/branch (log2)")
    # second, just fusions
    ax[1][1].plot(resultsdf["age"], np.log2(resultsdf["fusions_ratio"]), color='blue')
    ax[1][1].set_title("Fusions/branch (log2) vs Age")
    ax[1][1].set_xlabel("Age")
    ax[1][1].set_ylabel("Fusions/branch (log2)")
    # third, just losses
    ax[2][1].plot(resultsdf["age"], np.log2(resultsdf["losses_ratio"]), color='red')
    ax[2][1].set_title("Losses Ratio/branch (log2) vs Age")
    ax[2][1].set_xlabel("Age")
    ax[2][1].set_ylabel("Losses/branch (log2)")

    # make the axis limits go from -1200 to 0
    for axi in [0,1,2]:
        for axj in [0,1]:
            ax[axi][axj].set_xlim(-1200, 0)
    # save the plot as a pdf
    plt.savefig("branch_stats_vs_time.pdf")

    # now plot the intensity of extinction with the change types depending on whether they are present.
    if args.intensity_of_extinction is not None:
        plot_intensity_of_extinction(resultsdf, args.intensity_of_extinction)

if __name__== '__main__':
    main()