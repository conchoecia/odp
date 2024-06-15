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
import copy
import ete3
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import pandas as pd
from scipy.stats import pearsonr,spearmanr,kendalltau # I use scipy for the correlation functions
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

    # The third column uses the mean of the fusion rate at each age.
    # The fourth column uses the log2 of the fusion rate at each age.

    num_rows = 3
    if intensity_of_extinction_filepath is not None:
        num_rows = 5
        # read in the intensity of extinction file
        intensity_of_extinction_df = pd.read_csv(intensity_of_extinction_filepath, sep='\t')
        print(intensity_of_extinction_df)
        # first, we need to get the intensity of extinction for each age.

    red = "#D22C16"
    blue = "#3054A3"

    fig, ax = plt.subplots(num_rows, 4, figsize=(10, 10))
    # first, combined
    print(list(resultsdf.columns))
    print(resultsdf[["fusions_ratio", "losses_ratio", "fusion_rate_at_this_age_mean", "loss_rate_at_this_age_mean"]])
    for coli in [0,1,2,3]:
        # The age axis is the same in all of the plots
        for rowi in range(num_rows):
            ax[rowi][coli].set_xlabel("Age")

        if coli < 2:
            yfusion = "fusions_ratio"
            yloss   = "losses_ratio"
            ylabel = "Changes/branch"
        else:
            yfusion = "fusion_rate_at_this_age_mean"
            yloss   = "loss_rate_at_this_age_mean"
            ylabel  = "Changes/million years"
        # Determine if we're doing the log2 of the values
        label_append = ""
        if coli % 2 == 0:
            # we use the native values
            ax[0][coli].plot(resultsdf["age"], resultsdf[yfusion], color=blue)
            ax[0][coli].plot(resultsdf["age"], resultsdf[yloss]  ,   color=red)
            ax[1][coli].plot(resultsdf["age"], resultsdf[yfusion], color=blue)
            ax[2][coli].plot(resultsdf["age"], resultsdf[yloss]  , color=red)
        else:
            # we use the log2 of the values
            label_append = "(log2)"
            ax[0][coli].plot(resultsdf["age"], np.log2(resultsdf[yfusion].astype('float64')), color=blue )
            ax[0][coli].plot(resultsdf["age"], np.log2(resultsdf[yloss].astype('float64')  ), color=red  )
            ax[1][coli].plot(resultsdf["age"], np.log2(resultsdf[yfusion].astype('float64')), color=blue )
            ax[2][coli].plot(resultsdf["age"], np.log2(resultsdf[yloss].astype('float64')  ), color=red  )

        ax[0][coli].set_title("Fusions and Losses/branch vs Age")
        ax[0][coli].set_ylabel("Changes/branch" + label_append)
        # second, just fusions
        ax[1][coli].set_title("Fusions/branch vs Age")
        ax[1][coli].set_ylabel("Fusions/branch" + label_append)
        # third, just losses
        ax[2][coli].set_title("Losses/branch vs Age")
        ax[2][coli].set_ylabel("Losses/branch" + label_append)

        if intensity_of_extinction_filepath is not None:
            # This is intensity of extinction
            for i, row in intensity_of_extinction_df.iterrows():
                left_x = -1 * row['Time (Ma)']
                left_y = 0
                height = row['Extinction Intensity (%)']
                width = 1
                rectangle = patches.Rectangle((left_x, left_y), width, height, fill="#555555", edgecolor=None)
                ax[3][coli].add_patch(rectangle)
            # Add some text to the plots "Rhode & Miller (2005) Extinction Intensity"
            # Add it to the top left.
            ax[3][coli].text(0.05, 0.92, "Rhode & Miller (2005) Extinction Intensity", fontsize=8, transform=ax[3][0].transAxes)
            # This is intensity of origination
            for i, row in intensity_of_extinction_df.iterrows():
                left_x = -1 * row['Time (Ma)']
                left_y = 0
                height = row['Origination Intensity (%)']
                width = 1
                rectangle = patches.Rectangle((left_x, left_y), width, height, fill="#555555", edgecolor=None)
                ax[4][coli].add_patch(rectangle)
            # Add some text to the plots "Rhode & Miller (2005) Extinction Intensity"
            # Add it to the top left.
            ax[4][coli].text(0.05, 0.92, "Rhode & Miller (2005) Origination Intensity", fontsize=8, transform=ax[4][0].transAxes)

        # make the axis limits go from -1200 to 0
        xmin = -900
        xmax = 0
        for axi in range(num_rows):
            ax[axi][coli].set_xlim(xmin, xmax)

        # if we have the intensity of extinction, we need to scale the y-axis
        if intensity_of_extinction_filepath is not None:
            # get the values where 'Time (Ma)' is between xmin and xmax
            subdf = intensity_of_extinction_df[(intensity_of_extinction_df['Time (Ma)'] >= xmax) & (intensity_of_extinction_df['Time (Ma)'] <= -1 * xmin)]
            # the ylim is 1.1 times the maximum value
            ymax = 1.1 * subdf['Extinction Intensity (%)'].max()
            ax[3][coli].set_ylim(0, ymax)

            # Now do the same for origination intensity
            ymax = 1.1 * subdf['Origination Intensity (%)'].max()
            ax[4][coli].set_ylim(0, ymax)

        # change the fontsize of the axes and the titles
        fontsize = 8
        for axi in range(num_rows):
            ax[axi][coli].tick_params(axis='both', which='major', labelsize=fontsize)
            ax[axi][coli].set_title( ax[axi][coli].get_title(),   fontsize=fontsize)
            ax[axi][coli].set_xlabel(ax[axi][coli].get_xlabel(),  fontsize=fontsize)
            ax[axi][coli].set_ylabel(ax[axi][coli].get_ylabel(),  fontsize=fontsize)

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
    elif spearman_or_kendall == "pearson":
        corr, p = pearsonr(np.log2(statsdf[col1name]), np.log2(statsdf[col2name]))
    else:
        raise ValueError(f"Unknown correlation type {spearman_or_kendall}")
    return (corr, p, len(statsdf))

def plot_intensity_of_extinction(outprefix, count_df, intensity_of_extinction_filepath, suppress_plotting = False, fontsize = 8):
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
    # set the correct columns to type float64
    for thiscol in ["fusions_ratio", "losses_ratio", "fusion_rate_at_this_age_mean", "loss_rate_at_this_age_mean"]:
        plotdf[thiscol] = plotdf[thiscol].astype('float64')
    # We want to measure the spearman correlation between the intensity of extinction and the fusions and losses. Use log for the fusions and losses.
    #  - It is not appropriate to use Pearson correlation, because the data are not linear (the percent, ranging between 0-100 or 0-1, cannot be linearly related to the number of fusions or losses).
    # get the values that have no nans or inf values
    fusion_spearman_corr,   fusion_spearman_p,   fusion_spearman_n   = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusions_ratio", "spearman")
    fusion_kendalltau_corr, fusion_kendalltau_p, fusion_kendalltau_n = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusions_ratio", "kendall" )
    fusion_pearson_corr,    fusion_pearson_p,    fusion_pearson_n    = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusions_ratio", "pearson" )
    loss_spearman_corr,     loss_spearman_p,     loss_spearman_n     = _full_sk_stats(plotdf, "extinction_intensity_normalized", "losses_ratio",  "spearman" )
    loss_kendalltau_corr,   loss_kendalltau_p,   loss_kendalltau_n   = _full_sk_stats(plotdf, "extinction_intensity_normalized", "losses_ratio",  "kendall"  )
    loss_pearson_corr,      loss_pearson_p,      loss_pearson_n      = _full_sk_stats(plotdf, "extinction_intensity_normalized", "losses_ratio",  "pearson"  )
    statsresults = {}
    statsresults["extinction_fusion_numpbranch_spearman_r"]   = fusion_spearman_corr
    statsresults["extinction_fusion_numpbranch_spearman_p"]   = fusion_spearman_p
    statsresults["extinction_fusion_numpbranch_spearman_n"]   = fusion_spearman_n
    statsresults["extinction_fusion_numpbranch_kendalltau_r"]  = fusion_kendalltau_corr
    statsresults["extinction_fusion_numpbranch_kendalltau_p"]  = fusion_kendalltau_p
    statsresults["extinction_fusion_numpbranch_kendalltau_n"]  = fusion_kendalltau_n
    statsresults["extinction_fusion_numpbranch_pearson_r"]      = fusion_pearson_corr
    statsresults["extinction_fusion_numpbranch_pearson_p"]      = fusion_pearson_p
    statsresults["extinction_fusion_numpbranch_pearson_n"]      = fusion_pearson_n
    statsresults["extinction_losses_numpbranch_spearman_r"]   = loss_spearman_corr
    statsresults["extinction_losses_numpbranch_spearman_p"]   = loss_spearman_p
    statsresults["extinction_losses_numpbranch_spearman_n"]   = loss_spearman_n
    statsresults["extinction_losses_numpbranch_kendalltau_r"]  = loss_kendalltau_corr
    statsresults["extinction_losses_numpbranch_kendalltau_p"]  = loss_kendalltau_p
    statsresults["extinction_losses_numpbranch_kendalltau_n"]  = loss_kendalltau_n
    statsresults["extinction_losses_numpbranch_pearson_r"]      = loss_pearson_corr
    statsresults["extinction_losses_numpbranch_pearson_p"]      = loss_pearson_p
    statsresults["extinction_losses_numpbranch_pearson_n"]      = loss_pearson_n


    ofusion_spearman_corr,   ofusion_spearman_p,   ofusion_spearman_n   = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusions_ratio", "spearman")
    ofusion_kendalltau_corr, ofusion_kendalltau_p, ofusion_kendalltau_n = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusions_ratio", "kendall")
    ofusion_pearson_corr,    ofusion_pearson_p,    ofusion_pearson_n    = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusions_ratio", "pearson")
    oloss_spearman_corr,     oloss_spearman_p,     oloss_spearman_n     = _full_sk_stats(plotdf, "origination_intensity_normalized", "losses_ratio" , "spearman")
    oloss_kendalltau_corr,   oloss_kendalltau_p,   oloss_kendalltau_n   = _full_sk_stats(plotdf, "origination_intensity_normalized", "losses_ratio" , "kendall")
    oloss_pearson_corr,      oloss_pearson_p,      oloss_pearson_n      = _full_sk_stats(plotdf, "origination_intensity_normalized", "losses_ratio" , "pearson")

    statsresults["origination_fusion_numpbranch_spearman_r"]   = ofusion_spearman_corr
    statsresults["origination_fusion_numpbranch_spearman_p"]   = ofusion_spearman_p
    statsresults["origination_fusion_numpbranch_spearman_n"]   = ofusion_spearman_n
    statsresults["origination_fusion_numpbranch_kendalltau_r"]  = ofusion_kendalltau_corr
    statsresults["origination_fusion_numpbranch_kendalltau_p"]  = ofusion_kendalltau_p
    statsresults["origination_fusion_numpbranch_kendalltau_n"]  = ofusion_kendalltau_n
    statsresults["origination_fusion_numpbranch_pearson_r"]      = ofusion_pearson_corr
    statsresults["origination_fusion_numpbranch_pearson_p"]      = ofusion_pearson_p
    statsresults["origination_fusion_numpbranch_pearson_n"]      = ofusion_pearson_n
    statsresults["origination_losses_numpbranch_spearman_r"]   = oloss_spearman_corr
    statsresults["origination_losses_numpbranch_spearman_p"]   = oloss_spearman_p
    statsresults["origination_losses_numpbranch_spearman_n"]   = oloss_spearman_n
    statsresults["origination_losses_numpbranch_kendalltau_r"]  = oloss_kendalltau_corr
    statsresults["origination_losses_numpbranch_kendalltau_p"]  = oloss_kendalltau_p
    statsresults["origination_losses_numpbranch_kendalltau_n"]  = oloss_kendalltau_n
    statsresults["origination_losses_numpbranch_pearson_r"]      = oloss_pearson_corr
    statsresults["origination_losses_numpbranch_pearson_p"]      = oloss_pearson_p
    statsresults["origination_losses_numpbranch_pearson_n"]      = oloss_pearson_n

    # Now these are on the rates at this age mean
    rate_fusion_spearman_corr,   rate_fusion_spearman_p,   rate_fusion_spearman_n   = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusion_rate_at_this_age_mean", "spearman")
    rate_fusion_kendalltau_corr, rate_fusion_kendalltau_p, rate_fusion_kendalltau_n  = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusion_rate_at_this_age_mean", "kendall" )
    rate_fusion_pearson_corr,    rate_fusion_pearson_p,    rate_fusion_pearson_n      = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusion_rate_at_this_age_mean", "pearson" )
    rate_loss_spearman_corr,     rate_loss_spearman_p,     rate_loss_spearman_n     = _full_sk_stats(plotdf, "extinction_intensity_normalized", "loss_rate_at_this_age_mean",   "spearman" )
    rate_loss_kendalltau_corr,   rate_loss_kendalltau_p,   rate_loss_kendalltau_n    = _full_sk_stats(plotdf, "extinction_intensity_normalized", "loss_rate_at_this_age_mean",   "kendall"  )
    rate_loss_pearson_corr,      rate_loss_pearson_p,      rate_loss_pearson_n        = _full_sk_stats(plotdf, "extinction_intensity_normalized", "loss_rate_at_this_age_mean",   "pearson"  )
    statsresults["extinction_fusion_ratepmya_spearman_r"]   = rate_fusion_spearman_corr
    statsresults["extinction_fusion_ratepmya_spearman_p"]   = rate_fusion_spearman_p
    statsresults["extinction_fusion_ratepmya_spearman_n"]   = rate_fusion_spearman_n
    statsresults["extinction_fusion_ratepmya_kendalltau_r"]  = rate_fusion_kendalltau_corr
    statsresults["extinction_fusion_ratepmya_kendalltau_p"]  = rate_fusion_kendalltau_p
    statsresults["extinction_fusion_ratepmya_kendalltau_n"]  = rate_fusion_kendalltau_n
    statsresults["extinction_fusion_ratepmya_pearson_r"]      = rate_fusion_pearson_corr
    statsresults["extinction_fusion_ratepmya_pearson_p"]      = rate_fusion_pearson_p
    statsresults["extinction_fusion_ratepmya_pearson_n"]      = rate_fusion_pearson_n
    statsresults["extinction_losses_ratepmya_spearman_r"]   = rate_loss_spearman_corr
    statsresults["extinction_losses_ratepmya_spearman_p"]   = rate_loss_spearman_p
    statsresults["extinction_losses_ratepmya_spearman_n"]   = rate_loss_spearman_n
    statsresults["extinction_losses_ratepmya_kendalltau_r"]  = rate_loss_kendalltau_corr
    statsresults["extinction_losses_ratepmya_kendalltau_p"]  = rate_loss_kendalltau_p
    statsresults["extinction_losses_ratepmya_kendalltau_n"]  = rate_loss_kendalltau_n
    statsresults["extinction_losses_ratepmya_pearson_r"]      = rate_loss_pearson_corr
    statsresults["extinction_losses_ratepmya_pearson_p"]      = rate_loss_pearson_p
    statsresults["extinction_losses_ratepmya_pearson_n"]      = rate_loss_pearson_n


    rate_ofusion_spearman_corr,   rate_ofusion_spearman_p,   rate_ofusion_spearman_n   = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusion_rate_at_this_age_mean", "spearman")
    rate_ofusion_kendalltau_corr, rate_ofusion_kendalltau_p, rate_ofusion_kendalltau_n  = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusion_rate_at_this_age_mean", "kendall")
    rate_ofusion_pearson_corr,    rate_ofusion_pearson_p,    rate_ofusion_pearson_n      = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusion_rate_at_this_age_mean", "pearson")
    rate_oloss_spearman_corr,     rate_oloss_spearman_p,     rate_oloss_spearman_n     = _full_sk_stats(plotdf, "origination_intensity_normalized", "loss_rate_at_this_age_mean",   "spearman")
    rate_oloss_kendalltau_corr,   rate_oloss_kendalltau_p,   rate_oloss_kendalltau_n    = _full_sk_stats(plotdf, "origination_intensity_normalized", "loss_rate_at_this_age_mean",   "kendall")
    rate_oloss_pearson_corr,      rate_oloss_pearson_p,      rate_oloss_pearson_n        = _full_sk_stats(plotdf, "origination_intensity_normalized", "loss_rate_at_this_age_mean",   "pearson")

    statsresults["origination_fusion_ratepmya_spearman_r"]   = rate_ofusion_spearman_corr
    statsresults["origination_fusion_ratepmya_spearman_p"]   = rate_ofusion_spearman_p
    statsresults["origination_fusion_ratepmya_spearman_n"]   = rate_ofusion_spearman_n
    statsresults["origination_fusion_ratepmya_kendalltau_r"]  = rate_ofusion_kendalltau_corr
    statsresults["origination_fusion_ratepmya_kendalltau_p"]  = rate_ofusion_kendalltau_p
    statsresults["origination_fusion_ratepmya_kendalltau_n"]  = rate_ofusion_kendalltau_n
    statsresults["origination_fusion_ratepmya_pearson_r"]      = rate_ofusion_pearson_corr
    statsresults["origination_fusion_ratepmya_pearson_p"]      = rate_ofusion_pearson_p
    statsresults["origination_fusion_ratepmya_pearson_n"]      = rate_ofusion_pearson_n
    statsresults["origination_losses_ratepmya_spearman_r"]   = rate_oloss_spearman_corr
    statsresults["origination_losses_ratepmya_spearman_p"]   = rate_oloss_spearman_p
    statsresults["origination_losses_ratepmya_spearman_n"]   = rate_oloss_spearman_n
    statsresults["origination_losses_ratepmya_kendalltau_r"]  = rate_oloss_kendalltau_corr
    statsresults["origination_losses_ratepmya_kendalltau_p"]  = rate_oloss_kendalltau_p
    statsresults["origination_losses_ratepmya_kendalltau_n"]  = rate_oloss_kendalltau_n
    statsresults["origination_losses_ratepmya_pearson_r"]      = rate_oloss_pearson_corr
    statsresults["origination_losses_ratepmya_pearson_p"]      = rate_oloss_pearson_p
    statsresults["origination_losses_ratepmya_pearson_n"]      = rate_oloss_pearson_n

    def gen_stat_string(spearman_r, spearman_p, spearman_n, kendall_r, kendall_p, kendall_n, pearson_r, pearson_p, pearson_n) -> str:
        """
        generates a string with the statistics
        """
        outstring  = f"Spearman r={spearman_r:.2f}\n"
        outstring += f"Spearman p={spearman_p:.4e}\n"
        outstring += f"Spearman n={spearman_n}\n\n"
        outstring += f"Kendall r={kendall_r:.2f}\n"
        outstring += f"Kendall p={kendall_p:.4e}\n"
        outstring += f"Kendall n={kendall_n}\n\n"
        outstring += f"Pearson (log/log) r={pearson_r:.2f}\n"
        outstring += f"Pearson (log/log) p={pearson_p:.4e}\n"
        outstring += f"Pearson (log/log) n={pearson_n}"
        return outstring

    if not suppress_plotting:
        # EXTINCTION INTENSITY
        # now make the plots
        red = "#D22C16"
        blue = "#3054A3"
        num_rows = 6
        fig, ax = plt.subplots(num_rows, 4, figsize=(20, 20))
        fontsize = 8
        for coli in [0,1,2,3]:
            ylabel1 = "Changes/branch"
            ylabel2 = "Fusions/branch"
            ylabel3 = "Losses/branch"
            ylabel4 = "Changes/million years"
            ylabel5 = "Fusions/million years"
            ylabel6 = "Losses/million years"

            yfusiontop = plotdf["fusions_ratio"]
            yfusionbot = plotdf["fusion_rate_at_this_age_mean"]
            ylosstop   = plotdf["losses_ratio"]
            ylossbot   = plotdf["loss_rate_at_this_age_mean"]
            if coli < 2:
                xlabel = "Extinction Intensity (%)"
                x      = plotdf["intensity_of_extinction"]
                title1 = "Fusions and Losses/branch vs Extinction Intensity"
                title2 = "Fusions/branch vs Extinction Intensity"
                title3 = "Losses/branch vs Extinction Intensity"
                title4 = "Fusions and Losses/million years vs Extinction Intensity"
                title5 = "Fusions/million years vs Extinction Intensity"
                title6 = "Losses/million years vs Extinction Intensity"
                stattext1 = gen_stat_string(fusion_spearman_corr,      fusion_spearman_p,      fusion_spearman_n,      fusion_kendalltau_corr,      fusion_kendalltau_p,      fusion_kendalltau_n     ,      fusion_pearson_corr,      fusion_pearson_p,      fusion_pearson_n)
                stattext2 = gen_stat_string(loss_spearman_corr,        loss_spearman_p,        loss_spearman_n,        loss_kendalltau_corr,        loss_kendalltau_p,        loss_kendalltau_n       ,        loss_pearson_corr,        loss_pearson_p,        loss_pearson_n)
                stattext3 = gen_stat_string(rate_fusion_spearman_corr, rate_fusion_spearman_p, rate_fusion_spearman_n, rate_fusion_kendalltau_corr, rate_fusion_kendalltau_p, rate_fusion_kendalltau_n, rate_fusion_pearson_corr, rate_fusion_pearson_p, rate_fusion_pearson_n)
                stattext4 = gen_stat_string(rate_loss_spearman_corr,   rate_loss_spearman_p,   rate_loss_spearman_n,   rate_loss_kendalltau_corr,   rate_loss_kendalltau_p,   rate_loss_kendalltau_n  ,   rate_loss_pearson_corr,   rate_loss_pearson_p,   rate_loss_pearson_n)
            else:
                xlabel = "Origination Intensity (%)"
                x      = plotdf["intensity_of_origination"]
                title1 = "Fusions and Losses/branch vs Origination Intensity"
                title2 = "Fusions/branch vs Origination Intensity"
                title3 = "Losses/branch vs Origination Intensity"
                title4 = "Fusions and Losses/million years vs Origination Intensity"
                title5 = "Fusions/million years vs Origination Intensity"
                title6 = "Losses/million years vs Origination Intensity"
                stattext1 = gen_stat_string(ofusion_spearman_corr,      ofusion_spearman_p,      ofusion_spearman_n,      ofusion_kendalltau_corr,      ofusion_kendalltau_p,      ofusion_kendalltau_n     ,      ofusion_pearson_corr,      ofusion_pearson_p,      ofusion_pearson_n )
                stattext2 = gen_stat_string(oloss_spearman_corr,        oloss_spearman_p,        oloss_spearman_n,        oloss_kendalltau_corr,        oloss_kendalltau_p,        oloss_kendalltau_n       ,        oloss_pearson_corr,        oloss_pearson_p,        oloss_pearson_n )
                stattext3 = gen_stat_string(rate_ofusion_spearman_corr, rate_ofusion_spearman_p, rate_ofusion_spearman_n, rate_ofusion_kendalltau_corr, rate_ofusion_kendalltau_p, rate_ofusion_kendalltau_n, rate_ofusion_pearson_corr, rate_ofusion_pearson_p, rate_ofusion_pearson_n )
                stattext4 = gen_stat_string(rate_oloss_spearman_corr,   rate_oloss_spearman_p,   rate_oloss_spearman_n,   rate_oloss_kendalltau_corr,   rate_oloss_kendalltau_p,   rate_oloss_kendalltau_n  ,   rate_oloss_pearson_corr,   rate_oloss_pearson_p,   rate_oloss_pearson_n )
            for rowi in range(num_rows):
                ax[rowi][coli].set_xlabel(xlabel, fontsize=fontsize)
            if coli % 2 == 1:
                yfusiontop = np.log2(yfusiontop)
                yfusionbot = np.log2(yfusionbot)
                ylosstop   = np.log2(ylosstop)
                ylossbot   = np.log2(ylossbot)
            # first, combined
            ax[0][coli].set_title(title1, fontsize=fontsize)
            ax[0][coli].scatter(x, yfusiontop, color=blue)
            ax[0][coli].scatter(x, ylosstop  , color=red)
            ax[0][coli].set_ylabel(ylabel1, fontsize=fontsize)
            # second, just fusions
            ax[1][coli].scatter(x, yfusiontop, color=blue)
            ax[1][coli].set_title(title2,   fontsize=fontsize)
            ax[1][coli].set_ylabel(ylabel2, fontsize=fontsize)
            ax[1][coli].text(0.5, 0.5, stattext1, fontsize=fontsize, transform=ax[1][coli].transAxes, va = "center")
            # third, just losses
            ax[2][coli].scatter(x, ylosstop, color=red)
            ax[2][coli].set_title(title3,   fontsize=fontsize)
            ax[2][coli].set_ylabel(ylabel3, fontsize=fontsize)
            ax[2][coli].text(0.5, 0.5, stattext2, fontsize=fontsize, transform=ax[2][coli].transAxes, va = "center")

            # fourth, combined, but for the rate rather than for the ratio
            ax[3][coli].scatter(x, yfusionbot, color=blue)
            ax[3][coli].scatter(x, ylossbot  , color=red)
            ax[3][coli].set_title(title4,   fontsize=fontsize)
            ax[3][coli].set_ylabel(ylabel4, fontsize=fontsize)
            # second, just fusions
            ax[4][coli].scatter(x, yfusionbot, color=blue)
            ax[4][coli].set_title(title5,   fontsize=fontsize)
            ax[4][coli].set_ylabel(ylabel5, fontsize=fontsize)
            ax[4][coli].text(0.5, 0.5, stattext3, fontsize=fontsize, transform=ax[4][coli].transAxes, va = "center")
            # third, just losses
            ax[5][coli].scatter(x, ylossbot  , color=red)
            ax[5][coli].set_title(title6,   fontsize=fontsize)
            ax[5][coli].set_ylabel(ylabel6, fontsize=fontsize)
            ax[5][coli].text(0.5, 0.5, stattext4, fontsize=fontsize, transform=ax[5][coli].transAxes, va = "center")

        # change the fontsize of the axes and the titles
        for axi in range(num_rows):
            for axj in range(4):
                ax[axi][axj].tick_params(axis='both', which='major', labelsize=fontsize)
                ax[axi][axj].set_title(ax[axi][axj].get_title(), fontsize=fontsize)

        # force change the dot size
        dotsize = 1
        for axi in range(num_rows):
            for axj in [0,1,2,3]:
                for coll in ax[axi][axj].collections:
                    coll.set_sizes([dotsize])

        # change the opacity to 0.25
        opacity = 0.25
        for axi in range(num_rows):
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

def get_edge_stats_single_taxid(taxid, edgedf) -> pd.DataFrame:
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
    columns = {"total_branches_at_this_age":      0,
               "total_fusions_at_this_age":       0,
               "total_losses_at_this_age":        0,
               "fusion_rate_at_this_age_mean":   -1,
               "fusion_rate_at_this_age_median": -1,
               "fusion_rate_at_this_age_list":   [],
               "loss_rate_at_this_age_mean":     -1,
               "loss_rate_at_this_age_median":   -1,
               "loss_rate_at_this_age_list":     [],
               }

    # get the oldest age in source_age in the edgedf
    oldest_age = int(edgedf['parent_age'].max())
    # populate the data structure with one row per integer
    count_struct = {i: copy.deepcopy(columns) for i in range(oldest_age + 1)}

    total_fusions_in_this_branch = 0
    total_losses_in_this_branch = 0
    # Now go through and record the total number of branches at each age. This is done just using the edgedf.
    # This value is what we
    for i, row in edgedf.iterrows():
        # We just need to calculate this once
        num_fusions_my_this_branch = row["num_fusions_per_my_this_branch"]
        if np.isnan(num_fusions_my_this_branch) or np.isinf(num_fusions_my_this_branch):
            num_fusions_my_this_branch = 0
        num_losses_my_this_branch = row["num_losses_per_my_this_branch"]
        if np.isnan(num_losses_my_this_branch) or np.isinf(num_losses_my_this_branch):
            num_losses_my_this_branch = 0
        # now update the count_struct for each date
        for j in range(int(row['child_age']), int(row['parent_age'])+1):
            count_struct[j]['total_branches_at_this_age'] += 1
            # if value is nan or inf, set to 0
            count_struct[j]['fusion_rate_at_this_age_list'].append(num_fusions_my_this_branch)
            count_struct[j]['loss_rate_at_this_age_list'].append(num_losses_my_this_branch)

            # increment the total number of fusions and losses
            if row["num_fusions_this_branch"] > 0:
                count_struct[j]['total_fusions_at_this_age'] += row["num_fusions_this_branch"]
                total_fusions_in_this_branch += row["num_fusions_this_branch"]
            if row["num_losses_this_branch"] > 0:
                count_struct[j]['total_losses_at_this_age'] += row["num_losses_this_branch"]
                total_losses_in_this_branch += row["num_losses_this_branch"]

    # now calculate the mean rate for each position
    for i in range(oldest_age + 1):
        # make sure that the lists are the same lengths as the total_branches_at_this_age
        #print("this age: {}".format(i))
        #print("this branch: {}".format(taxid))
        #print("Len of total branches at this age: {}".format(count_struct[i]['total_branches_at_this_age']))
        #print("Len of fusion rate at this age list: {}".format(len(count_struct[i]['fusion_rate_at_this_age_list'])))
        #print("Len of loss rate at this age list: {}".format(len(count_struct[i]['loss_rate_at_this_age_list'])))
        if not count_struct[i]['total_branches_at_this_age'] == len(count_struct[i]['fusion_rate_at_this_age_list']):
            raise ValueError("The number of branches at this age does not match the number of fusions at this age.")
        if not count_struct[i]['total_branches_at_this_age'] == len(count_struct[i]['loss_rate_at_this_age_list']):
            raise ValueError("The number of branches at this age does not match the number of losses at this age.")
        count_struct[i]['fusion_rate_at_this_age_mean']   = np.mean(count_struct[i]['fusion_rate_at_this_age_list'])
        count_struct[i]['fusion_rate_at_this_age_median'] = np.median(count_struct[i]['fusion_rate_at_this_age_list'])
        count_struct[i]['loss_rate_at_this_age_mean']     = np.mean(count_struct[i]['loss_rate_at_this_age_list'])
        count_struct[i]['loss_rate_at_this_age_median']   = np.median(count_struct[i]['loss_rate_at_this_age_list'])
        #del count_struct[i]['fusion_rate_at_this_age_list']
        #del count_struct[i]['loss_rate_at_this_age_list']

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

    # We want some clade names
    NCBI = ete3.NCBITaxa()
    for i, row in nodedf.iterrows():
        node = row["taxid"]
        clade_name = NCBI.get_taxid_translator([node])[node]
        clade_name = clade_name[0].upper() + clade_name[1:]

        if " " in clade_name:
            clade_name = "".join([x.capitalize() for x in clade_name.split(" ")])
        print("  - We are processing node {} / {}".format(i, len(nodedf)))
        resultsdf, (fusions, losses) = get_edge_stats_single_taxid(node, edgedf)
        outprefix0 = f"{clade_name}_{node}_changes_vs_age"
        resultsdf.to_csv(outprefix0 + ".tsv", sep='\t', index=False)
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