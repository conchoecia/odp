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
import matplotlib.patches as patches
import os
import pandas as pd
from scipy.stats import spearmanr,kendalltau # I use scipy for the correlation functions
import sys

from odp_plotting_functions import format_matplotlib

def parse_args():
    """
    We need two files:
      - the edge_information file.
      - the node information file
      - the chromosome information file
        - From this we want the median number of chromosomes for each taxid
      - the statsdf file.
      - the intensity of extinction file.
      - an optional list of taxids to omit, anything in the clade that has one of these taxids as a parent (or as the actual taxid) will be omitted.
      - a flag R to just rerun the incomplete parts of the analysis. Doesn't bother to output a node file.

    We will check that both files exist.
    """
    parser = argparse.ArgumentParser(description='Plot branch stats vs time')
    parser.add_argument("-e", "--edge_information", type=str, help='The edge information file')
    parser.add_argument("-n", "--node_information", type=str, help='The node information file')
    parser.add_argument("-s", "--statsdf", type=str, help='The df_stats file')
    parser.add_argument("-S", "--suppress_plotting", action='store_true', help='Suppress plotting')
    parser.add_argument("-i", "--intensity_of_extinction", type=str, help="The intensity of extinction file. Must have columns: \'Time (Ma)\' \'Diversity All Genera\' \'Diversity Short-Lived\' \'Diversity Long-Lived\' \'Diversity Well-Resolved\' \'Extinction Intensity (%)\' \'Origination Intensity(%)\'")
    parser.add_argument("-o", "--omit_taxids", type=str, help='A string with a list of taxids to omit, separated by commas.')
    args = parser.parse_args()

    if not os.path.exists(args.edge_information):
        raise ValueError('The edge information file does not exist')

    if not os.path.exists(args.node_information):
        raise ValueError('The node information file does not exist')

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

def plot_fusions_per_branch_vs_time(outprefix, resultsdf, intensity_of_extinction_filepath = None):
    """
    makes a plot of the number of fusions per branch vs time.
    """
    # call the function to properly format the text
    format_matplotlib()

    # Now make three plots, each top of one another.
    # The 0th row shows both fusions and fissions ratios together in the same axis.
    # The 1st row shows just the fusion ratio on the axis.
    # The 2nd row shows just the fission ratio on the axis.
    # Use two different panels, plot the fusions in the top panel, and the losses in the bottom panel. Color the fusions blue and the losses red.
    # We also make a second axis which is the log2 of the ratio of fusions to losses.

    num_rows = 3
    if intensity_of_extinction_filepath is not None:
        num_rows = 5
        # read in the intensity of extinction file
        intensity_of_extinction_df = pd.read_csv(intensity_of_extinction_filepath, sep='\t')
        print(intensity_of_extinction_df)
        # first, we need to get the intensity of extinction for each age.

    red = "#D22C16"
    blue = "#3054A3"

    fig, ax = plt.subplots(num_rows, 2, figsize=(10, 10))
    # first, combined
    ax[0][0].set_title("Fusions and Losses/branch vs Age")
    ax[0][0].plot(resultsdf["age"], resultsdf["fusions_ratio"], color=blue)
    ax[0][0].plot(resultsdf["age"], resultsdf["losses_ratio"], color=red)
    ax[0][0].set_xlabel("Age")
    ax[0][0].set_ylabel("Changes/branch")
    # second, just fusions
    ax[1][0].plot(resultsdf["age"], resultsdf["fusions_ratio"], color=blue)
    ax[1][0].set_title("Fusions/branch vs Age")
    ax[1][0].set_xlabel("Age")
    ax[1][0].set_ylabel("Fusions/branch")
    # third, just losses
    ax[2][0].plot(resultsdf["age"], resultsdf["losses_ratio"], color=red)
    ax[2][0].set_title("Losses/branch vs Age")
    ax[2][0].set_xlabel("Age")
    ax[2][0].set_ylabel("Losses/branch")

    if intensity_of_extinction_filepath is not None:
        # This is intensity of extinction
        for i, row in intensity_of_extinction_df.iterrows():
            left_x = -1 * row['Time (Ma)']
            left_y = 0
            height = row['Extinction Intensity (%)']
            width = 1
            rectangle = patches.Rectangle((left_x, left_y), width, height, fill="#555555", edgecolor=None)
            ax[3][0].add_patch(rectangle)
            rectangle = patches.Rectangle((left_x, left_y), width, height, fill="#555555", edgecolor=None)
            ax[3][1].add_patch(rectangle)
        # Add some text to the plots "Rhode & Miller (2005) Extinction Intensity"
        # Add it to the top left.
        ax[3][0].text(0.05, 0.92, "Rhode & Miller (2005) Extinction Intensity", fontsize=8, transform=ax[3][0].transAxes)
        ax[3][1].text(0.05, 0.92, "Rhode & Miller (2005) Extinction Intensity", fontsize=8, transform=ax[3][1].transAxes)
        # This is intensity of origination
        for i, row in intensity_of_extinction_df.iterrows():
            left_x = -1 * row['Time (Ma)']
            left_y = 0
            height = row['Origination Intensity (%)']
            width = 1
            rectangle = patches.Rectangle((left_x, left_y), width, height, fill="#555555", edgecolor=None)
            ax[4][0].add_patch(rectangle)
            rectangle = patches.Rectangle((left_x, left_y), width, height, fill="#555555", edgecolor=None)
            ax[4][1].add_patch(rectangle)
        # Add some text to the plots "Rhode & Miller (2005) Extinction Intensity"
        # Add it to the top left.
        ax[4][0].text(0.05, 0.92, "Rhode & Miller (2005) Origination Intensity", fontsize=8, transform=ax[4][0].transAxes)
        ax[4][1].text(0.05, 0.92, "Rhode & Miller (2005) Origination Intensity", fontsize=8, transform=ax[4][1].transAxes)

    # now make the same thing, but log2 of the ratio
    # first plot, combined on one axis
    ax[0][1].set_title("Fusions and Losses/branch vs Age")
    ax[0][1].plot(resultsdf["age"], np.log2(resultsdf["fusions_ratio"]), color=blue)
    ax[0][1].plot(resultsdf["age"], np.log2(resultsdf["losses_ratio"]), color=red)
    ax[0][1].set_xlabel("Age")
    ax[0][1].set_ylabel("Changes/branch (log2)")
    # second, just fusions
    ax[1][1].plot(resultsdf["age"], np.log2(resultsdf["fusions_ratio"]), color=red)
    ax[1][1].set_title("Fusions/branch (log2) vs Age")
    ax[1][1].set_xlabel("Age")
    ax[1][1].set_ylabel("Fusions/branch (log2)")
    # third, just losses
    ax[2][1].plot(resultsdf["age"], np.log2(resultsdf["losses_ratio"]), color=red)
    ax[2][1].set_title("Losses Ratio/branch (log2) vs Age")
    ax[2][1].set_xlabel("Age")
    ax[2][1].set_ylabel("Losses/branch (log2)")

    # make the axis limits go from -1200 to 0
    xmin = -900
    xmax = 0
    for axi in range(num_rows):
        for axj in [0,1]:
            ax[axi][axj].set_xlim(xmin, xmax)

    # if we have the intensity of extinction, we need to scale the y-axis
    if intensity_of_extinction_filepath is not None:
        # get the values where 'Time (Ma)' is between xmin and xmax
        subdf = intensity_of_extinction_df[(intensity_of_extinction_df['Time (Ma)'] >= xmax) & (intensity_of_extinction_df['Time (Ma)'] <= -1 * xmin)]
        # the ylim is 1.1 times the maximum value
        ymax = 1.1 * subdf['Extinction Intensity (%)'].max()
        ax[3][0].set_ylim(0, ymax)
        ax[3][1].set_ylim(0, ymax)

        # Now do the same for origination intensity
        ymax = 1.1 * subdf['Origination Intensity (%)'].max()
        ax[4][0].set_ylim(0, ymax)
        ax[4][1].set_ylim(0, ymax)

    # change the fontsize of the axes and the titles
    fontsize = 8
    for axi in range(num_rows):
        for axj in [0,1]:
            ax[axi][axj].tick_params(axis='both', which='major', labelsize=fontsize)
            ax[axi][axj].set_title(ax[axi][axj].get_title(), fontsize=fontsize)
            ax[axi][axj].set_xlabel(ax[axi][axj].get_xlabel(), fontsize=fontsize)
            ax[axi][axj].set_ylabel(ax[axi][axj].get_ylabel(), fontsize=fontsize)

    # increase the horizontal and vertical space between the panels
    plt.subplots_adjust(hspace=0.5, wspace=0.5)

    # save the plot as a pdf
    pdfout = outprefix.rstrip(".pdf") + ".pdf"
    plt.savefig(pdfout)
    # close the figure to free up memory
    plt.close(fig)

def _subdf_no_missing(df, col1name, col2name) -> pd.DataFrame:
    """
    Given two column names and a df, extract those two columns, then remove rows with missing values, or values that are -inf or inf.
    A view of the original df is fine since we do not modify the subdf.
    Returns a view of the original df.
    """
    statsdf = df[[col1name, col2name]]
    return statsdf[~statsdf.isin([np.nan, np.inf, -np.inf]).any(axis=1)]

def _full_sk_stats(df, col1name, col2name, spearman_or_kendall) -> tuple:
    """
    Given a df and two column names,
    Returns a 3-element tuple of correlation, the p-value, and n (sample size)
    """
    statsdf = _subdf_no_missing(df, col1name, col2name)
    if spearman_or_kendall == "spearman":
        corr, p = spearmanr(statsdf[col1name], statsdf[col2name])
    elif spearman_or_kendall == "kendall":
        corr, p = kendalltau(statsdf[col1name], statsdf[col2name])
    else:
        raise ValueError(f"Unknown correlation type {spearman_or_kendall}")
    return (corr, p, len(statsdf))

def plot_intensity_of_extinction(outprefix, count_df, intensity_of_extinction_filepath, suppress_plotting = False):
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
    mya_to_extinction_intensity  = {row['Time (Ma)']: row['Extinction Intensity (%)'] for i, row in intensity_of_extinction_df.iterrows()}
    mya_to_origination_intensity = {row['Time (Ma)']: row['Origination Intensity (%)'] for i, row in intensity_of_extinction_df.iterrows()}

    # map the ages, but if the ages aren't present, give a value of -1
    count_df["age_positive"] = count_df["age"].apply(lambda x: -1 * x)
    count_df["intensity_of_extinction"]  = count_df["age_positive"].map(mya_to_extinction_intensity).fillna(-1)
    count_df["intensity_of_origination"] = count_df["age_positive"].map(mya_to_origination_intensity).fillna(-1)

    # only plot the values where intensity of extinction is not -1
    plotdf = count_df[count_df["intensity_of_extinction"] != -1].copy()
    plotdf["extinction_intensity_normalized"]  = plotdf["intensity_of_extinction"] / 100
    plotdf["origination_intensity_normalized"] = plotdf["intensity_of_origination"] / 100
    # We want to measure the spearman correlation between the intensity of extinction and the fusions and losses. Use log for the fusions and losses.
    #  - It is not appropriate to use Pearson correlation, because the data are not linear (the percent, ranging between 0-100 or 0-1, cannot be linearly related to the number of fusions or losses).
    # get the values that have no nans or inf values
    fusion_spearman_corr,   fusion_spearman_p,   fusion_spearman_n   = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusions_ratio", "spearman")
    loss_spearman_corr,     loss_spearman_p,     loss_spearman_n     = _full_sk_stats(plotdf, "extinction_intensity_normalized", "losses_ratio", "spearman")
    fusion_kendalltau_corr, fusion_kendalltau_p, fusion_kendalltau_n = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusions_ratio", "kendall")
    loss_kendalltau_corr,   loss_kendalltau_p,   loss_kendalltau_n   = _full_sk_stats(plotdf, "extinction_intensity_normalized", "losses_ratio", "kendall")
    statsresults = {}
    statsresults["extinction_fusion_spearman_r"]   = fusion_spearman_corr
    statsresults["extinction_fusion_spearman_p"]   = fusion_spearman_p
    statsresults["extinction_fusion_spearman_n"]   = fusion_spearman_n
    statsresults["extinction_fusion_kendalltau_r"] = fusion_kendalltau_corr
    statsresults["extinction_fusion_kendalltau_p"] = fusion_kendalltau_p
    statsresults["extinction_fusion_kendalltau_n"] = fusion_kendalltau_n
    statsresults["extinction_losses_spearman_r"]   = loss_spearman_corr
    statsresults["extinction_losses_spearman_p"]   = loss_spearman_p
    statsresults["extinction_losses_spearman_n"]   = loss_spearman_n
    statsresults["extinction_losses_kendalltau_r"] = loss_kendalltau_corr
    statsresults["extinction_losses_kendalltau_p"] = loss_kendalltau_p
    statsresults["extinction_losses_kendalltau_n"] = loss_kendalltau_n


    ofusion_spearman_corr,   ofusion_spearman_p,   ofusion_spearman_n   = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusions_ratio", "spearman")
    oloss_spearman_corr,     oloss_spearman_p,     oloss_spearman_n     = _full_sk_stats(plotdf, "origination_intensity_normalized", "losses_ratio" , "spearman")
    ofusion_kendalltau_corr, ofusion_kendalltau_p, ofusion_kendalltau_n = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusions_ratio", "kendall")
    oloss_kendalltau_corr,   oloss_kendalltau_p,   oloss_kendalltau_n   = _full_sk_stats(plotdf, "origination_intensity_normalized", "losses_ratio" , "kendall")
    statsresults["origination_fusion_spearman_r"]   = ofusion_spearman_corr
    statsresults["origination_fusion_spearman_p"]   = ofusion_spearman_p
    statsresults["origination_fusion_spearman_n"]   = ofusion_spearman_n
    statsresults["origination_fusion_kendalltau_r"] = ofusion_kendalltau_corr
    statsresults["origination_fusion_kendalltau_p"] = ofusion_kendalltau_p
    statsresults["origination_fusion_kendalltau_n"] = ofusion_kendalltau_n
    statsresults["origination_losses_spearman_r"]   = oloss_spearman_corr
    statsresults["origination_losses_spearman_p"]   = oloss_spearman_p
    statsresults["origination_losses_spearman_n"]   = oloss_spearman_n
    statsresults["origination_losses_kendalltau_r"] = oloss_kendalltau_corr
    statsresults["origination_losses_kendalltau_p"] = oloss_kendalltau_p
    statsresults["origination_losses_kendalltau_n"] = oloss_kendalltau_n

    if not suppress_plotting:
        # EXTINCTION INTENSITY
        # now make the plots
        red = "#D22C16"
        blue = "#3054A3"
        fig, ax = plt.subplots(3, 4, figsize=(20, 10))
        fontsize = 8
        # first, combined
        ax[0][0].set_title("Fusions and Losses/branch vs Extinction Intensity")
        ax[0][0].scatter(plotdf["intensity_of_extinction"], plotdf["fusions_ratio"], color=blue)
        ax[0][0].scatter(plotdf["intensity_of_extinction"], plotdf["losses_ratio"], color=red)
        ax[0][0].set_xlabel("Extinction Intensity (%)")
        ax[0][0].set_ylabel("Changes/branch")
        # second, just fusions
        ax[1][0].scatter(plotdf["intensity_of_extinction"], plotdf["fusions_ratio"], color=blue)
        ax[1][0].set_title("Fusions/branch vs Extinction Intensity")
        ax[1][0].set_xlabel("Extinction Intensity (%)")
        ax[1][0].set_ylabel("Fusions/branch")
        # add the spearman and kendall tau information
        ax[1][0].text(0.5, 0.5, f"Spearman correlation: {fusion_spearman_corr:.2f}\nSpearman p={fusion_spearman_p:.4e}\nSpearman n={fusion_spearman_n}\n\nKendall tau: {fusion_kendalltau_corr:.2f}\nKendall p={fusion_kendalltau_p:.4e}\nKendall n={fusion_kendalltau_n}", fontsize=fontsize, transform=ax[1][0].transAxes)
        # third, just losses
        ax[2][0].scatter(plotdf["intensity_of_extinction"], plotdf["losses_ratio"], color=red)
        ax[2][0].set_title("Losses/branch vs Extinction Intensity")
        ax[2][0].set_xlabel("Extinction Intensity (%)")
        ax[2][0].set_ylabel("Losses/branch")
        ax[2][0].text(0.5, 0.5, f"Spearman correlation: {loss_spearman_corr:.2f}\nSpearman p={loss_spearman_p:.4e}\nSpearman n={loss_spearman_n}\nKendall tau: {loss_kendalltau_corr:.2f}\nKendall p={loss_kendalltau_p:.4e}\nKendall n={loss_kendalltau_n}", fontsize=fontsize, transform=ax[2][0].transAxes)

        # now make the same thing, but log2 of the ratio
        # first plot, combined on one axis
        ax[0][1].set_title("Fusions and Losses/branch vs Extinction Intensity")
        ax[0][1].scatter(plotdf["intensity_of_extinction"], np.log2(plotdf["fusions_ratio"]), color=blue)
        ax[0][1].scatter(plotdf["intensity_of_extinction"], np.log2(plotdf["losses_ratio"]), color=red)
        ax[0][1].set_xlabel("Extinction Intensity (%)")
        ax[0][1].set_ylabel("Changes/branch (log2)")
        # second, just fusions
        ax[1][1].scatter(plotdf["intensity_of_extinction"], np.log2(plotdf["fusions_ratio"]), color=blue)
        ax[1][1].set_title("Fusions/branch (log2) vs Extinction Intensity")
        ax[1][1].set_xlabel("Extinction Intensity (%)")
        ax[1][1].set_ylabel("Fusions/branch (log2)")
        # third, just losses
        ax[2][1].scatter(plotdf["intensity_of_extinction"], np.log2(plotdf["losses_ratio"]), color=red)
        ax[2][1].set_title("Losses Ratio/branch (log2) vs Extinction Intensity")
        ax[2][1].set_xlabel("Extinction Intensity (%)")
        ax[2][1].set_ylabel("Losses/branch (log2)")

        # ORIGINATION INTENSITY
        # first, combined
        ax[0][2].set_title("Fusions and Losses/branch vs Origination Intensity")
        ax[0][2].scatter(plotdf["intensity_of_origination"], plotdf["fusions_ratio"], color=blue)
        ax[0][2].scatter(plotdf["intensity_of_origination"], plotdf["losses_ratio"], color=red)
        ax[0][2].set_xlabel("Origination Intensity (%)")
        ax[0][2].set_ylabel("Changes/branch")
        # second, just fusions
        ax[1][2].scatter(plotdf["intensity_of_origination"], plotdf["fusions_ratio"], color=blue)
        ax[1][2].set_title("Fusions/branch vs Origination Intensity")
        ax[1][2].set_xlabel("Origination Intensity (%)")
        ax[1][2].set_ylabel("Fusions/branch")
        # add the spearman and kendall tau information
        ax[1][2].text(0.5, 0.5, f"Spearman correlation: {ofusion_spearman_corr:.2f}\nSpearman p={ofusion_spearman_p:.4e}\nSpearman n={ofusion_spearman_n}\nKendall tau: {ofusion_kendalltau_corr:.2f}\nKendall p={ofusion_kendalltau_p:.4e}\nKendall n={ofusion_kendalltau_n}", fontsize=fontsize, transform=ax[1][2].transAxes)
        # third, just losses
        ax[2][2].scatter(plotdf["intensity_of_origination"], plotdf["losses_ratio"], color=red)
        ax[2][2].set_title("Losses/branch vs Origination Intensity")
        ax[2][2].set_xlabel("Origination Intensity (%)")
        ax[2][2].set_ylabel("Losses/branch")
        ax[2][2].text(0.5, 0.5, f"Spearman correlation: {oloss_spearman_corr:.2f}\nSpearman p={oloss_spearman_p:.4e}\nSpearman n={oloss_spearman_n}\nKendall tau: {oloss_kendalltau_corr:.2f}\nKendall p={oloss_kendalltau_p:.4e}\nKendall n={oloss_kendalltau_n}", fontsize=fontsize, transform=ax[2][2].transAxes)

        # now make the same thing, but log2 of the ratio
        # first plot, combined on one axis
        ax[0][3].set_title("Fusions and Losses/branch vs Origination Intensity")
        ax[0][3].scatter(plotdf["intensity_of_origination"], np.log2(plotdf["fusions_ratio"]), color=blue)
        ax[0][3].scatter(plotdf["intensity_of_origination"], np.log2(plotdf["losses_ratio"]), color=red)
        ax[0][3].set_xlabel("Origination Intensity (%)")
        ax[0][3].set_ylabel("Changes/branch (log2)")
        # second, just fusions
        ax[1][3].scatter(plotdf["intensity_of_origination"], np.log2(plotdf["fusions_ratio"]), color=blue)
        ax[1][3].set_title("Fusions/branch (log2) vs Origination Intensity")
        ax[1][3].set_xlabel("Origination Intensity (%)")
        ax[1][3].set_ylabel("Fusions/branch (log2)")
        # third, just losses
        ax[2][3].scatter(plotdf["intensity_of_origination"], np.log2(plotdf["losses_ratio"]), color=red)
        ax[2][3].set_title("Losses Ratio/branch (log2) vs Origination Intensity")
        ax[2][3].set_xlabel("Origination Intensity (%)")
        ax[2][3].set_ylabel("Losses/branch (log2)")

        # change the fontsize of the axes and the titles
        fontsize = 8
        for axi in [0,1,2]:
            for axj in [0,1,2,3]:
                ax[axi][axj].tick_params(axis='both', which='major', labelsize=fontsize)
                ax[axi][axj].set_title(ax[axi][axj].get_title(), fontsize=fontsize)

        # force change the dot size
        dotsize = 1
        for axi in [0,1,2]:
            for axj in [0,1,2,3]:
                for coll in ax[axi][axj].collections:
                    coll.set_sizes([dotsize])

        # change the opacity to 0.25
        opacity = 0.25
        for axi in [0,1,2]:
            for axj in [0,1,2,3]:
                for coll in ax[axi][axj].collections:
                    coll.set_alpha(opacity)

        # increase the horizontal and vertical space between the panels
        plt.subplots_adjust(hspace=0.5, wspace=0.5)

        # save to a pdf
        outpdf =  outprefix.rstrip(".pdf").rstrip(".tsv") + ".pdf"
        plt.savefig(outpdf)
        # close the figure to free up memory
        plt.close(fig)
    return statsresults

def get_edge_stats_single_taxid(taxid, edgedf, nodedf, eventdf) -> pd.DataFrame:
    """
    This function takes the edgedf and builds up statistics about the edges.
    The things we will calculate:
      - The rate of fusions on this specific branch
      - The rate of fusions on this and child branches

    Returns:
      - A dataframe with event counts at each age.
    """
    # we need to filter the edgedf where the lineage of the child_taxid is the taxid
    keep_indices = []
    for i, row in edgedf.iterrows():
        if taxid in eval(row["child_lineage"]):
            keep_indices.append(i)
    edgedf = edgedf.iloc[keep_indices]
    edgedf.reset_index(inplace=True, drop=True)

    # build a structured vector to hold the counts at each age.
    # get the unique source ages
    # make a datastructure with the counts of the source ages
    columns = {"total_branches_at_this_age": 0,
               "total_losses_at_this_age":   0,
               "total_fusions_at_this_age":  0}

    # get the oldest age in source_age in the edgedf
    oldest_age = int(edgedf['parent_age'].max())
    # populate the data structure with one row per integer
    count_struct = {i: columns.copy() for i in range(oldest_age + 1)}

    total_fusions_in_this_branch = 0
    total_losses_in_this_branch = 0
    # Now go through and record the total number of branches at each age. This is done just using the edgedf.
    # This value is what we
    for i, row in edgedf.iterrows():
        for j in range(int(row['child_age']), int(row['parent_age'])+1):
            count_struct[j]['total_branches_at_this_age'] += 1
            # increment the total number of fusions and losses
            if row["num_fusions_this_branch"] > 0:
                count_struct[j]['total_fusions_at_this_age'] += row["num_fusions_this_branch"]
                total_fusions_in_this_branch += row["num_fusions_this_branch"]
            if row["num_losses_this_branch"] > 0:
                count_struct[j]['total_losses_at_this_age'] += row["num_losses_this_branch"]
                total_losses_in_this_branch += row["num_losses_this_branch"]

    # all of the data should be in the struct now. Turn it into a pandas df
    resultsdf = pd.DataFrame(count_struct).T

    for event in ["fusions", "losses"]:
        resultsdf[f"{event}_ratio"] = resultsdf[f"total_{event}_at_this_age"] / resultsdf["total_branches_at_this_age"]
    resultsdf["age"] = -1 * resultsdf.index
    return resultsdf, (total_fusions_in_this_branch, total_losses_in_this_branch)

def add_events_to_edge_df(edgedf, eventdf) -> pd.DataFrame:
    """
    All that this does is tracks the types and quantity of changes on each edge.
    This information will later be used in the other functions for plotting.
    """
    # we need a taxid_to_parent dictionary to get the mapping to id the branch length
    # Actually, we can just get the index of this row, because we now know that the target column is unique,
    #   and we will want to work with multiple values from this row.
    taxid_to_parent = {row['child_taxid']: row['parent_taxid'] for i, row in edgedf.iterrows()}

    # Initialize a dictionary to hold the new columns
    new_columns = {"num_fusions_this_branch": 0,
                   "num_losses_this_branch": 0,
                   "num_fusions_per_my_this_branch": 0,
                   "num_losses_per_my_this_branch": 0}
    # get all of the unique categories without the parentheses from the eventdf
    ALGs         = sorted(set([x for x in eventdf['change'].unique() if x[0] != '(']))
    for thisALG in ALGs:
        new_columns["num_" + thisALG + "_this_branch"] = 0

    # Get all of the possible combinations
    Combinations = sorted(set([x for x in eventdf['change'].unique() if x[0] == '(']))
    for thisCombo in Combinations:
        combostr = "+".join(eval(thisCombo))
        new_columns["num_" + combostr + "_this_branch"] = 0
    new_columns_df = pd.DataFrame(new_columns, index=edgedf.index)
    # Concatenate the new columns to the original DataFrame
    edgedf = pd.concat([edgedf, new_columns_df], axis=1)

    missing_events = 0
    for i, row in eventdf.iterrows():
        # This is a single event. We need to know the branch on which it occurred.
        event_branch_id = row["target_taxid"]
        if event_branch_id not in taxid_to_parent:
            missing_events += 1
            continue
        else:
            parent_taxid = taxid_to_parent[event_branch_id]
            edge = (parent_taxid, event_branch_id)

            # get the indices of the edgedf that this matches,
            edge_row = edgedf[(edgedf['parent_taxid'] == parent_taxid) & (edgedf['child_taxid'] == event_branch_id)]
            if len(edge_row) > 1:
                raise ValueError('For some reason we found this edge more than once. There is probably an error with the output software. {} {}'.format(
                    parent_taxid, event_branch_id))
            elif len(edge_row) == 0:
                missing_events += 1
                continue
            else:
                # determine whether the event was a fusion or a loss
                colnames = []
                changetype = None
                if row["change"][0] == "(":
                    changetype = "fusions"
                    colnames.append("num_fusions_this_branch")
                    combostr = "+".join(eval(row["change"]))
                    colnames.append("num_" + combostr + "_this_branch")
                else:
                    changetype = "losses"
                    colnames.append("num_losses_this_branch")
                    colnames.append("num_" + row["change"] + "_this_branch")
                for thiscol in colnames:
                    edgedf.loc[edge_row.index, thiscol] += 1
    print("There were {} missing events".format(missing_events))

    edgedf["num_fusions_per_my_this_branch"] = edgedf["num_fusions_this_branch"] / edgedf["branch_length"]
    edgedf["num_losses_per_my_this_branch"] = edgedf["num_losses_this_branch"] / edgedf["branch_length"]
    return edgedf

def main():
    args = parse_args()

    # read in the dataframes
    edgedf  = pd.read_csv(args.edge_information, sep='\t')
    nodedf  = pd.read_csv(args.node_information, sep='\t')
    # enforce that node parent is int. First, convert NaNs to type None
    nodedf['parent'] = nodedf['parent'].apply(lambda x: -1 if pd.isna(x) else int(x)).astype(int)
    eventdf = pd.read_csv(args.statsdf, sep='\t')

    print("This is nodedf")
    print(nodedf)
    print("This is edgedf")
    print(edgedf)

    edgedf = add_events_to_edge_df(edgedf, eventdf)
    print("Edgedf has been marked up")
    # save this to a new file, as it has the number of changes per branch
    edgedf.to_csv("modified_edge_list.tsv", sep='\t', index=False)

    # print the columns that have _this_branch

    # verify that there every value in the target column is unique
    if len(edgedf['child_taxid'].unique()) != len(edgedf):
        raise ValueError('The target column (child_taxid) is not unique')
    else:
        print("Every value in the target (child_taxid) column is unique. This is a legal DAG.")

    #for node in nodedf['taxid']:
    nodedf["fusions_in_this_clade"] = -1
    nodedf["losses_in_this_clade"] = -1
    # change the types of the previous two to int
    nodedf["fusions_in_this_clade"] = nodedf["fusions_in_this_clade"].astype(int)
    nodedf["losses_in_this_clade"] = nodedf["losses_in_this_clade"].astype(int)
    nodedf["fusions_in_this_clade_div_dist_crown"] = -1.00000000001
    nodedf["losses_in_this_clade_div_dist_crown"]  = -1.00000000001
    nodedf["fusions_in_this_clade_div_dist_crown_plus_root"] = -1.00000000001
    nodedf["losses_in_this_clade_div_dist_crown_plus_root"]  = -1.00000000001
    if args.intensity_of_extinction is not None:
        nodedf["extinction_fusion_spearman_r"]    = -1.00000000001
        nodedf["extinction_fusion_spearman_p"]    = -1.00000000001
        nodedf["extinction_fusion_spearman_n"]    = -1.00000000001
        nodedf["extinction_fusion_kendalltau_r"]  = -1.00000000001
        nodedf["extinction_fusion_kendalltau_p"]  = -1.00000000001
        nodedf["extinction_fusion_kendalltau_n"]  = -1.00000000001
        nodedf["extinction_losses_spearman_r"]    = -1.00000000001
        nodedf["extinction_losses_spearman_p"]    = -1.00000000001
        nodedf["extinction_losses_spearman_n"]    = -1.00000000001
        nodedf["extinction_losses_kendalltau_r"]  = -1.00000000001
        nodedf["extinction_losses_kendalltau_p"]  = -1.00000000001
        nodedf["extinction_losses_kendalltau_n"]  = -1.00000000001
        nodedf["origination_fusion_spearman_r"]   = -1.00000000001
        nodedf["origination_fusion_spearman_p"]   = -1.00000000001
        nodedf["origination_fusion_spearman_n"]   = -1.00000000001
        nodedf["origination_fusion_kendalltau_r"] = -1.00000000001
        nodedf["origination_fusion_kendalltau_p"] = -1.00000000001
        nodedf["origination_fusion_kendalltau_n"] = -1.00000000001
        nodedf["origination_losses_spearman_r"]   = -1.00000000001
        nodedf["origination_losses_spearman_p"]   = -1.00000000001
        nodedf["origination_losses_spearman_n"]   = -1.00000000001
        nodedf["origination_losses_kendalltau_r"] = -1.00000000001
        nodedf["origination_losses_kendalltau_p"] = -1.00000000001
        nodedf["origination_losses_kendalltau_n"] = -1.00000000001

    # We want some clade names
    NCBI = ete3.NCBITaxa()
    for i, row in nodedf.iterrows():
        node = row["taxid"]
        clade_name = NCBI.get_taxid_translator([node])[node]
        clade_name = clade_name[0].upper() + clade_name[1:]

        if " " in clade_name:
            clade_name = "".join([x.capitalize() for x in clade_name.split(" ")])
        print("  - We are processing node {} / {}".format(i, len(nodedf)))
        resultsdf, (fusions, losses) = get_edge_stats_single_taxid(node, edgedf, nodedf, eventdf)
        # update the nodedf with the number of fusions and losses
        # we can use iloc because we have i
        nodedf.loc[i, 'fusions_in_this_clade'] = fusions
        nodedf.loc[i, 'losses_in_this_clade'] = losses
        # if the name field is empty or Nan, add it
        if nodedf.loc[i, "name"] == "" or pd.isna(nodedf.loc[i, "name"]):
            nodedf.loc[i, "name"] = clade_name

        outprefix1 = f"{clade_name}_{node}_changes_vs_time"
        if not args.suppress_plotting:
            # plot the changes per time
            if args.intensity_of_extinction is not None:
                plot_fusions_per_branch_vs_time(outprefix1, resultsdf, intensity_of_extinction_filepath = args.intensity_of_extinction)
            else:
                plot_fusions_per_branch_vs_time(outprefix1, resultsdf)

        outprefix2 = f"{clade_name}_{node}_changes_vs_intensity"
        # now plot the intensity of extinction with the change types depending on whether they are present.
        if args.intensity_of_extinction is not None:
            statsdf = plot_intensity_of_extinction(outprefix2, resultsdf, args.intensity_of_extinction, suppress_plotting = args.suppress_plotting)
            for key in statsdf:
                nodedf.loc[nodedf['taxid'] == node, key] = statsdf[key]

    # now calculate the rates
    nodedf["fusions_in_this_clade_div_dist_crown"] = nodedf["fusions_in_this_clade"] / nodedf["dist_crown"]
    nodedf["losses_in_this_clade_div_dist_crown"]  = nodedf["losses_in_this_clade"] / nodedf["dist_crown"]
    nodedf["fusions_in_this_clade_div_dist_crown_plus_root"] = nodedf["fusions_in_this_clade"] / nodedf["dist_crown_plus_root"]
    nodedf["losses_in_this_clade_div_dist_crown_plus_root"]  = nodedf["losses_in_this_clade"] / nodedf["dist_crown_plus_root"]
    # now save the modified node df
    nodedf.to_csv("modified_node_list.tsv", sep='\t', index=False)

if __name__== '__main__':
    main()