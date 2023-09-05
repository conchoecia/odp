#!/usr/bin/env python

"""
    Filename:   plot_decay_pairwise_steps_test.py
   File type:   python script (.py)
      Author:   darrin t schultz (github: @conchoecia)
Date created:   July 29th, 2023

Description:
  - The purpose of this script is to demonstrate ALG decay by performing pairwise comparisons
     over a tree. For the purposes of demonstration there will be parameters that are hardcoded
     to work with specific species combinations. The end goal is that these analyses will work
     on a phylogenetic tree.
"""

# plotting options
import matplotlib.pyplot as plt
import odp_plotting_functions as odp_plot
#
import argparse
import pandas as pd
import os
import random
import sys
import yaml
from math import ceil as ceil

# set up argparse method to get the directory of the .tsv files we want to plot
def parse_args():
    """
    For the final file we will need a tree file,
    """
    parser = argparse.ArgumentParser(description="Plot the pairwise decay of chromosomes over a tree.")
    parser.add_argument("-t", "--tree", help="The tree file (newick) to use for the analysis. The tree must have divergence dates estimated for the nodes.")
    # path to a directory containing the .rbh files with significance values calculated
    parser.add_argument("-d", "--directory", help="The directory containing the .rbh files to use for the analysis.")
    # Path to a config.yaml file that contains the parameters for the analysis we want.
    # Temporary until the tree functionality is added.
    parser.add_argument("-c", "--config", help="The config.yaml file containing the parameters for the analysis. Temporary until trees are added.")
    # add a chrom file for the species we want to plot. This will help us determine which things are chroms
    parser.add_argument("-C", "--chrom", help = "The .chrom file for the species we want to plot. This will help us determine which things are chromosomes.")
    # Target_species is the species we want to plot
    # We plot this species against everything else
    parser.add_argument("-T", "--target_species", help="The species we want to plot. This species will be plotted against all other species found in the config file.")
    args = parser.parse_args()
    # if the length of args is 0 print the help and quit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # check that the tree file actually exists
    for thisfile in [args.tree]:
        if thisfile is not None:
            if not os.path.isfile(thisfile):
                raise Exception("The file you provided does not exist: {}".format(thisfile))
        else:
            pass # no biggie if we didn't specify this file. We probably don't care.
    # check that the directory actually exists
    if not os.path.isdir(args.directory):
        raise Exception("The directory you provided does not exist.")
    return args

def read_yaml_file(file_path):
    """
    from ChatGPT
    """
    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)
    return data

def parse_config(config_file, directory_of_rbh_files, target_species):
    """
    This will be a yaml file with analyses, detailing the target species to plot,
      the other species, and the divergence times.

    """
    species_prefix_to_filename = {}
    filelist = os.listdir(directory_of_rbh_files) # get a list of files in the directory
    config = read_yaml_file(config_file)          # load the config file
    # add the target species to the config if it isn't there already
    if "target_species" not in config:
        config["target_species"] = [target_species]

    # check if the necessary field, analyses, is in the file
    if "divergence_times" not in config:
        raise IOError("The field 'divergence_times' must be present in the config file.")

    # Ensure that the target_species is in the divergence_times field
    unseen_target_species = set([x for x in config["target_species"]])
    all_species_set = set()
    for sp1 in config["divergence_times"].keys():
        if sp1 in config["target_species"]:
            # remove this from the set of unseen target species
            # don't raise an error if it isn't in the set
            unseen_target_species.discard(sp1)
        for sp2 in config["divergence_times"][sp1].keys():
            all_species_set.add(sp1)
            all_species_set.add(sp2)
            if sp2 in config["target_species"]:
                # remove this from the set of unseen target species
                unseen_target_species.discard(sp2)
    if len(unseen_target_species) > 0:
        raise IOError("The target_species [{}] is not in the divergence_times field.".format(unseen_target_species))

    # This is not efficient, but go through the loop again to determine which species pairs we need
    #  We need this because later we need to check that all of these species pairs have files
    sp_to_expected_pairs = {x: set() for x in config["target_species"]}
    if "analyses" not in config:
        config["analyses"] = {x:{} for x in config["target_species"]}

    for target in config["analyses"]:
        for sp in all_species_set:
            if sp != target:
                sptup = tuple(sorted((target, sp)))
                sp_to_expected_pairs[target].add(sptup)
                # now we also need to add the divergence time to the analyses field
                if sp not in config["analyses"][target]:
                    s12 = tuple(sorted((target, sp)))
                    config["analyses"][target][sp] = config["divergence_times"][s12[0]][s12[1]]
    #for entry in sorted(config["analyses"]["Pectenmaximus6579"]):
    #    print("{}: {}".format(entry, config["analyses"]["Pectenmaximus6579"][entry]))
    #print(len(config["analyses"]["Pectenmaximus6579"]))

    # FILE PAIRING
    # for each analysis we now need to get the file to look at for the pairwise comparison
    if "analysis_files" not in config:
        config["analysis_files"] = {x:{} for x in config["target_species"]}

    # get the pairs from the files
    pair_to_file = {}
    for thisfile in [x for x in filelist if x.endswith(".rbh")]:
        fields = thisfile.split("_")
        analysis_pair = tuple(sorted((fields[0], fields[1])))
        if analysis_pair not in pair_to_file:
            pair_to_file[analysis_pair] = []
        complete_filepath = os.path.join(directory_of_rbh_files, thisfile)
        pair_to_file[analysis_pair].append(complete_filepath)

    # go through the target species analysis pairs as inferred from the divergence_times field
    #sp_to_expected_pairs = {x: set() for x in config["target_species"]}
    for target in config["target_species"]:
        for analysis_pair in sp_to_expected_pairs[target]:
            if analysis_pair not in pair_to_file:
                raise IOError("The analysis pair {} is not in the filelist.".format(analysis_pair))
            not_target = [x for x in analysis_pair if x != target][0]
            config["analysis_files"][target][not_target] = pair_to_file[analysis_pair][0]

    return config

def decay_of_one_species_pair(rawdf, sp1, sp2, sp_to_chr_to_size):
    """
    This gets the 1:1-ish chromosomes from the vantage point of sp1.
    Returns a dictionary of sp1 chromosomes as keys, and lists of sp2
     chromosomes that are significant matches.

    Returns a decay dataframe that looks like this:

    #  sp1_scaf  sp2_scaf            sp1_scaf_genecount conserved  scattered
    #  PMA1      ['AHY11','AHY13']   590                339        251
    #  PMA10     ['AHY11', 'AHY9']   361                229        132
    #  PMA11     ['AHY9']            423                312        111

    TODO: We need to add all of the SP1 chromosomes, even if they don't have any significant matches.
    """
    # whole_FET is the column where we record the significance using FET between the two columns
    df = rawdf.copy()
    df = df[df["whole_FET"] < 0.05]
    # At this point we only need the columns that are {sp1|sp2}_scaf and whole_FET.
    #  We have the necessary information for checking for 1:1 relationships between the chromosomes
    df = df[[x for x in df.columns if x.endswith("_scaf") or x == "whole_FET"]]
    # filter the {sp1|sp2}_scaf to only contain scaffold IDs with lengths > MIN_SCAF_LEN in the sp_to_chr_to_size dict
    #  We need to do this because we don't want to plot the small scaffolds
    MIN_SCAF_LEN = 500000
    df = df[df["{}_scaf".format(sp1)].isin([x for x in sp_to_chr_to_size[sp1].keys() if sp_to_chr_to_size[sp1][x] > MIN_SCAF_LEN])]
    # do the same thing for sp2
    df = df[df["{}_scaf".format(sp2)].isin([x for x in sp_to_chr_to_size[sp2].keys() if sp_to_chr_to_size[sp2][x] > MIN_SCAF_LEN])]
    # for each chromosome in SP1, get how many significantly related chromosomes there are in SP2
    #  We need to consider 1:1 relationships
    #  We also need to consider 1:many relationships
    gb = df.groupby(["{}_scaf".format(sp1)])
    sp1_scaf_to_sp2_scaf = {x[0]:x[1]["{}_scaf".format(sp2)].unique().tolist()
                            for x in gb}
    del df
    del gb

    # For each scaf in sp1, get the percent of genes that are conserved in sp2 on the significant scaffolds from sp1_scaf_to_sp2_scaf
    sp1_scaf_to_total_genes = rawdf.groupby(["{}_scaf".format(sp1)]).size().to_dict()
    # groupby both sp1 scafs and sp2 scafs
    sp1gb = rawdf.groupby(["{}_scaf".format(sp1), "{}_scaf".format(sp2)])
    # for each group go through and figure out the number of genes that are on the chromosome pairs in sp1_scaf_to_sp2_scaf
    sp1_scaf_to_conserved_genes = {}
    for k in sp1_scaf_to_sp2_scaf.keys():
        sp1_scaf_to_conserved_genes[k] = 0
        for l in sp1_scaf_to_sp2_scaf[k]:
            if (k, l) in sp1gb.groups:
                sp1_scaf_to_conserved_genes[k] += len(sp1gb.groups[(k, l)])

    # now we can make a decay dataframe
    # The format will be like this:
    #  sp1_scaf  sp2_scaf            sp1_scaf_genecount conserved  scattered
    #  PMA1      ['AHY11','AHY13']   590                339        251
    #  PMA10     ['AHY11', 'AHY9']   361                229        132
    #  PMA11     ['AHY9']            423                312        111
    entries = []
    for sp1_scaf in sp_to_chr_to_size[sp1].keys():
        thisentry = { "sp1_scaf": sp1_scaf,
                      "sp2_scaf": [],
                      "sp1_scaf_genecount": 0,
                      "conserved":          0,
                      "scattered":          0
                      }
        if sp1_scaf in sp1_scaf_to_total_genes.keys():
            thisentry["sp1_scaf_genecount"] = sp1_scaf_to_total_genes[sp1_scaf]
            thisentry["scattered"]          = sp1_scaf_to_total_genes[sp1_scaf]
        if sp1_scaf in sp1_scaf_to_conserved_genes.keys():
            thisentry["conserved"]          = sp1_scaf_to_conserved_genes[sp1_scaf]
        if sp1_scaf in sp1_scaf_to_sp2_scaf.keys():
            thisentry["sp2_scaf"]           = sp1_scaf_to_sp2_scaf[sp1_scaf]
        if (sp1_scaf in sp1_scaf_to_total_genes) and (sp1_scaf in sp1_scaf_to_conserved_genes):
            thisentry["scattered"]          = sp1_scaf_to_total_genes[sp1_scaf] - sp1_scaf_to_conserved_genes[sp1_scaf]
        entries.append(thisentry)
    decaydf = pd.DataFrame(entries)

    return decaydf

def jitter(iterable, maxjitter):
    """
    Add a little jitter to the iterable using the maxjitter value
      as the maximum or minimum value to add.
    """
    import random
    return [max(0, x) + random.uniform(-maxjitter, maxjitter) for x in iterable]

def rbh_files_to_sp_to_chr_to_size(rbh_filelist):
    """
    Get the chromosome sizes by using the rbh files and getting the max gene indices
    """
    sp_to_chr_to_size = {}
    sp_to_scaf_to_genecount = {}
    for thisfile in rbh_filelist:
        # load as pandas df
        df = pd.read_csv(thisfile, sep="\t")
        allsp = [x.replace("_gene", "") for x in df.columns if x.endswith("_gene")]
        for sp in allsp:
            spdf = df[[x for x in df.columns if x.startswith(sp)]]
            # groupby {}_scaf
            spgrp = spdf.groupby(["{}_scaf".format(sp)])
            # make a dict of the {}_scaf as the key and the max {}_pos column as the value
            dict_of_maxes = {x[0]:x[1]["{}_pos".format(sp)].max() for x in spgrp}
            # Update sp_to_chr_to_size, making sure the add this species if it doesn't exist.
            #  We also need to make sure that the chromosome is in the dictionary.
            #  Only update if the current max is bigger than the existing max
            if sp not in sp_to_chr_to_size:
                sp_to_chr_to_size[sp] = {}
            for scaf in dict_of_maxes.keys():
                if scaf not in sp_to_chr_to_size[sp]:
                    sp_to_chr_to_size[sp][scaf] = 0
            for k in dict_of_maxes.keys():
                if dict_of_maxes[k] > sp_to_chr_to_size[sp][k]:
                    sp_to_chr_to_size[sp][k] = dict_of_maxes[k]

            # now add the gene list to the sp_to_scaf_to_genecount
            if sp not in sp_to_scaf_to_genecount:
                sp_to_scaf_to_genecount[sp] = {}
            for scaf in df["{}_scaf".format(sp)].unique().tolist():
                if scaf not in sp_to_scaf_to_genecount[sp]:
                    sp_to_scaf_to_genecount[sp][scaf] = set()
                subdf = df[df["{}_scaf".format(sp)] == scaf]
                sp_to_scaf_to_genecount[sp][scaf].update(subdf["{}_gene".format(sp)].tolist())

    # now return the size of the set for sp_to_scaf_to_genecount
    for sp in sp_to_scaf_to_genecount.keys():
        for scaf in sp_to_scaf_to_genecount[sp].keys():
            sp_to_scaf_to_genecount[sp][scaf] = len(sp_to_scaf_to_genecount[sp][scaf])

    return sp_to_chr_to_size, sp_to_scaf_to_genecount

def calculate_pairwise_decay_sp1_vs_many(sp1, config, sp_to_chr_to_size,
                                         sp_to_keepscafs, outdir="./"):
    """
    Calculates the pairwise chromosomal decay between two species.
    Saves the decay dataframes to files. Each file is sp1 vs sp2.
    """
    # just hold onto these until the end to avoid writing only some files
    sp_to_decay_df = {}
    # iterate through the pairs of species
    for sp2 in config["analyses"][sp1].keys():
        analysis_pair = tuple(sorted((sp1, sp2)))
        rbhfile = config["analysis_files"][sp1][sp2]
        # read in the rbh file as a pandas df
        rawdf = pd.read_csv(rbhfile, sep="\t")

        # get the corresponding chromosomes
        sp1_sp2_decay = decay_of_one_species_pair(rawdf, sp1, sp2, sp_to_chr_to_size)

        # Add the divergence times to the dataframe
        sp1_sp2_decay["divergence_time"] = config["analyses"][sp1][sp2]

        # Add a percent conserved column
        sp1_sp2_decay["fraction_conserved"] = sp1_sp2_decay["conserved"] / sp1_sp2_decay["sp1_scaf_genecount"]

        # only keep the scaffolds for species 1 that we know are valid
        sp1_sp2_decay = sp1_sp2_decay[sp1_sp2_decay["sp1_scaf"].isin(sp_to_keepscafs[sp1])]

        # stash this to save to a file later
        sp_to_decay_df[sp2] = sp1_sp2_decay
        #print("\n", sp2, "\n", sp1_sp2_decay, file = sys.stderr)

    # safely make the outdir if it doesn't exist
    os.makedirs(outdir, exist_ok=True)

    sp_to_file_df = {}
    # now save the files and save the paths to a structure
    for sp2 in config["analyses"][sp1].keys():
        outprefix = "{}_vs_{}_chromosomal_decay.tsv".format(sp1, sp2)
        outdir_prefix = os.path.join(outdir, outprefix)
        # save the decay dataframe
        sp_to_decay_df[sp2].to_csv(outdir_prefix, sep="\t", index=False)
        # save the path to the decay dataframe
        sp_to_file_df[sp2] = outdir_prefix

    return {sp1: sp_to_file_df}

def plot_pairwise_decay_sp1_vs_all(sp1, filestruct, outdir="./"):
    """
    This takes a list of files and plots the decay of sp1 vs all the other species
      Does this for whole chromosomes and for whole genomes (chromosomes summed).

    The left subplot will have the whole-genome conservation vs divergence time.
    The right subplot will have the per-chromosome conservation vs divergence time.

    sp1 is the focal species that appears in every pairwise comparison
    filestruct is a dictionary of dictionaries. The first key is the focal species.
      The second keys are the species to which the focal species is being compared.
      The value of the second key is the path to the tsv file to use for the comparison.

    Plot [0][0] (top-left)  is the whole-genome conservation vs divergence time.
         [0][1] (top-right) is the per-chromosome conservation vs divergence time.
         [1][0] (bottom-left)  is a violin plot of every 25 million years of divergence time-whole genomes
         [1][1] (bottom-right) is a violin plot of every 25 million years of divergence time-per-chromosome
    """
    # BIN_SIZE is the number of millions of years to bin the data
    BIN_SIZE = 50
    sp_bins  = { x:[] for x in range(BIN_SIZE, 1500, BIN_SIZE) }
    chr_bins = { x:[] for x in range(BIN_SIZE, 1500, BIN_SIZE) }
    decay_bins
    bin_to_index = { x: int((x/BIN_SIZE) - 1)+2 for x in range(BIN_SIZE, 1500, BIN_SIZE)}
    # This sort order is used to make sure that all of the chromosomes are in the same order
    # for the later plots.
    sp1_sort_order = []

    # We start by making a subplot array to make the two variants. A two-column, one-row plot.
    NUMBER_OF_ROWS = 2 + len(sp_bins.keys())
    NUMBER_OF_COLS = 2
    fig, axes = plt.subplots(NUMBER_OF_ROWS, NUMBER_OF_COLS, figsize = (7.5 * NUMBER_OF_COLS, 6 * NUMBER_OF_ROWS))
    fig.suptitle("{} decay versus divergence time".format(sp1))
    for sp2 in filestruct[sp1].keys():
        sp1_sp2_decay = pd.read_csv(filestruct[sp1][sp2], sep="\t")
        if sp1_sort_order == []:
            # we should sort this by the chromosome size, smallest to largest
            sp1_sort_order = sp1_sp2_decay.sort_values(by="sp1_scaf_genecount", ascending=True)["sp1_scaf"].tolist()

        # sort the dataframe by the sp1_sort_order
        sp1_sp2_decay = sp1_sp2_decay.set_index("sp1_scaf").reindex(sp1_sort_order).reset_index()
        # get the most abundant divergence time from sp1_sp2_decay
        # They should all be the same, but this is the most robust thing to do.
        divergence_time = sp1_sp2_decay["divergence_time"].mode()[0]
        # figure out which bin this divergence time goes into.
        # The bin will be a multiple of 25.
        # The correct bin is the closest multiple of 25 to the divergence time, rounded up
        divergence_time_bin = int(BIN_SIZE * ceil(divergence_time/BIN_SIZE))
        # Make a whole-genome version of the dataframe. Just sum up the columns and recalculate the percent conserved
        sp1_sp2_whole = sp1_sp2_decay.sum(axis=0).to_frame().transpose()
        # only keep certain columns
        sp1_sp2_whole = sp1_sp2_whole[["sp1_scaf_genecount", "conserved", "scattered"]]
        sp1_sp2_whole["divergence_time"] = divergence_time
        sp1_sp2_whole["fraction_conserved"] = sp1_sp2_whole["conserved"] / sp1_sp2_whole["sp1_scaf_genecount"]

        #on the left-plot just do a scatterplot of the fraction conserved vs divergence time
        axes[0][0].scatter(sp1_sp2_whole["divergence_time"], sp1_sp2_whole["fraction_conserved"], label = "{}".format(sp2))
        axes[0][0].set_xlabel("Divergence time (MYA)")
        axes[0][0].set_ylabel("Fraction conserved on orthologous chromosomes")
        axes[0][0].set_title("Whole-genome conservation vs divergence time")

        # the right plot is per-chromosome. Add a little jitter to the x-axis so we can see the points
        axes[0][1].scatter(jitter(sp1_sp2_decay["divergence_time"], 20), sp1_sp2_decay["fraction_conserved"],
                        label = "{}".format(sp2), alpha = 0.1, edgecolors='none')
        axes[0][1].set_xlabel("Divergence time (MYA) (+- 20 MYA jitter)")
        axes[0][1].set_ylabel("Fraction conserved on orthologous chromosomes")
        axes[0][1].set_title("Orthologous chromosome conservation vs divergence time")

        # add these values to the bins for plotting later
        sp_bins[divergence_time_bin].append( sp1_sp2_whole["fraction_conserved"].tolist()[0])
        chr_bins[divergence_time_bin].extend(sp1_sp2_decay["fraction_conserved"].tolist())

        # make the line plots of the decay per-bin
        plotting_bin = bin_to_index[divergence_time_bin]
        # A



    # remove the empty bins
    sp_bins  = {x-(BIN_SIZE/2):sp_bins[x] for x in sp_bins.keys() if len(sp_bins[x]) > 0}
    chr_bins = {x-(BIN_SIZE/2):chr_bins[x] for x in chr_bins.keys() if len(chr_bins[x]) > 0}
    # now we need to make the violin plots of all the bins
    axes[1][0].violinplot(sp_bins.values(), sp_bins.keys(),   points=10, widths=20, showmeans=False, showextrema=True, showmedians=True)
    axes[1][1].violinplot(chr_bins.values(), chr_bins.keys(), points=5, widths=20, showmeans=False, showextrema=True, showmedians=True)

    # get the xlims and ylims for the top right plot
    xlims = axes[0][1].get_xlim()
    # set the xlims of all the plots
    axes[0][0].set_xlim(xlims)
    axes[1][0].set_xlim(xlims)
    axes[1][1].set_xlim(xlims)

    # set the ylims of all the figures from -0.05 to 1
    axes[0][0].set_ylim(-0.05, 1.05)
    axes[0][1].set_ylim(-0.05, 1.05)
    axes[1][0].set_ylim(-0.05, 1.05)
    axes[1][1].set_ylim(-0.05, 1.05)

    # CALL THIS TO GET THE VISUAL STYLE WE NEED
    odp_plot.format_matplotlib()

    # safely make the output directory if it does not yet exist
    os.makedirs(outdir, exist_ok=True)
    # Save the plot as a jpeg
    outprefix = "{}_decay_plot_vs_divergence_time".format(sp1)
    outdir_prefix = os.path.join(outdir, outprefix)
    plt.savefig("{}.pdf".format(outdir_prefix), format='pdf')


def plot_decay_twospecies(sp1, sp2, path_to_tsv, scaffolds_to_keep_sp1, outdir):
    """
    This plots the decay of an ALG between number of genes in the main chromosome,
    and the number of genes in smaller chromosomes

    Parameters:
        sp1: One species that is being plotted. This text title should match the text in the tsv file.
        sp2: Another species that is being plotted. This text should match the text of the tsv file.
        path_to_tsv: The path to the tsv file that contains the data to plot. The structure is described below.
        scaffolds_to_keep_sp1: A list of scaffolds to keep for sp1. This is necessary because we don't want to plot the small scaffolds.
        outdir: The directory to which we will save the plot.

    The input is the tsv output by calculate_pairwise_decay_sp1_vs_many:
    For example, here is one df with PMA (scallop) as sp1 and PFI (sponge) as sp2:

        sp1_scaf        sp2_scaf  sp1_scaf_genecount  conserved  scattered  divergence_time  fraction_conserved
     0      PMA1    [PFI8, PFI1]                 552        322        230              800            0.583333
     1     PMA10   [PFI13, PFI1]                 310        182        128              800            0.587097
     2     PMA11          [PFI7]                 376        262        114              800            0.696809

    The output is one figure with two subplots.
    The left subplot is the ranked sizes of the chromosomes in sp1. The right subplot is the actual size of the chromosomes in sp1
    """
    BARS = True
    df = pd.read_csv(path_to_tsv, sep="\t")
    # only keep the scaffolds for species 1 that we know are valid
    df = df[df["sp1_scaf"].isin(scaffolds_to_keep_sp1)]

    # rank the chromosomes based on their size and sort by the rank
    df["sp1_ranked"] = df["sp1_scaf_genecount"].rank(ascending=True, method="first")
    df = df.sort_values(by="sp1_ranked")

    # set up the two panels of the plot
    NUMBER_OF_ROWS = 2
    NUMBER_OF_COLS = 2
    fig, axes = plt.subplots(NUMBER_OF_ROWS, NUMBER_OF_COLS,
                             figsize = (7.5 * NUMBER_OF_COLS, 7*(NUMBER_OF_ROWS - 1)),
                             gridspec_kw={'height_ratios': [1, 2.5]})

    fig.suptitle("{} and {} chromosome conservation vs {} chromosome size".format(sp1, sp2, sp1))

    if BARS:
        # on the left plot the total gene count as red bars
        axes[1][0].bar(df["sp1_ranked"], df["sp1_scaf_genecount"], color = "red", alpha = 1)
    else:
        # plot the chromosome sizes by rank on the left
        axes[1][0].plot(df["sp1_ranked"], df["sp1_scaf_genecount"], "bo")

    # plot the chromosomes by actual size on the right
    axes[1][1].plot(df["sp1_scaf_genecount"], df["sp1_scaf_genecount"], "bo")

    # add some horizontal space between axes[0] and axes[1]
    fig.subplots_adjust(wspace=0.5)

    # make vertical lines on the left and the right plot. Do it by iterating through the dataframe
    for index, row in df.iterrows():
        if BARS:
            pass
        else:
            axes[1][0].plot([row["sp1_ranked"], row["sp1_ranked"]], [0, row["sp1_scaf_genecount"]], "k-", alpha = 0.33)
        axes[1][1].plot([row["sp1_scaf_genecount"], row["sp1_scaf_genecount"]], [0, row["sp1_scaf_genecount"]], "k-", alpha = 0.33)

    # now we plot the blue points
    if BARS:
        # now we plot blue bars for the number of genes conserved
        axes[1][0].bar(df["sp1_ranked"], df["conserved"], color = "blue", alpha = 1)
    else:
        # now we plot blue points for number of genes degraded
        axes[1][0].plot(df["sp1_ranked"], df["scattered"], "ro")
    axes[1][1].plot(df["sp1_scaf_genecount"], df["scattered"], "ro")

    # add some yaxis labels. make the color red to match the dots. Then make the tick labels red too
    color = "blue"
    axes[1][0].set_ylabel("Number of orthologs on chromosome", color=color)
    axes[1][1].set_ylabel("Number of orthologs on chromosome", color=color)
    axes[1][0].tick_params(axis='y', labelcolor=color)
    axes[1][1].tick_params(axis='y', labelcolor=color)

    # add some xaxis labels
    axes[1][0].set_xlabel("Chromosome ranked by ortholog count")
    axes[1][1].set_xlabel("Number of orthologs on chromosome")

    # on the left side we will add x-axis ticks that at the sp1 chromosome names, and rotate everything 45 degrees
    axes[1][0].set_xticks(df["sp1_ranked"])
    axes[1][0].set_xticklabels(df["sp1_scaf"], rotation=90, ha="center")

    # get the y-axis limits for the left plot
    left_ylim = axes[1][0].get_ylim()
    left_xlim = axes[1][0].get_xlim()

    # get the max value of the sp1_scaf_genecount
    left_maxgene = df["sp1_scaf_genecount"].max()
    # print out the ratios of the total limits to the max gene count
    ylim_scale_factor = abs(1 - (left_ylim[1]/left_maxgene))

    # set the y-limits of the left plot to the right plot
    axes[1][0].set_ylim(axes[1][1].get_ylim())

    # now we clone the axes and plot the percent conserved on the y-axes
    axL = axes[0][0]
    axR = axes[0][1]
    ylim_scale_number = 100
    axL.set_ylim( [-5, 100] )
    axL.set_xlim( left_xlim )
    axR.set_ylim( [-5, 100] )

    color = 'black'
    axL.set_ylabel('percent conserved on ALGs', color=color)  # we already handled the x-label with ax1
    axR.set_ylabel('percent conserved on ALGs', color=color)
    # set color of axL and axR yaxis ticks to blue
    axL.tick_params(axis='y', labelcolor=color)
    axR.tick_params(axis='y', labelcolor=color)
    # change the x-axis ticks to be the same positions as the bottom left
    # turn off the tick labels for the top left
    axL.set_xticks(df["sp1_ranked"])
    axL.set_xticklabels(["" for x in df["sp1_ranked"]])
    # now plot the data
    axL.plot(df["sp1_ranked"], 100*(df["conserved"]/df["sp1_scaf_genecount"]), color = color, lw = 1)
    axR.plot(df["sp1_scaf_genecount"], 100*(df["conserved"]/df["sp1_scaf_genecount"]), color = color, lw = 1)

    # turn off the top and right bar of the frames of each axes
    for ax in axes.flatten():
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # CALL THIS TO GET THE VISUAL STYLE WE NEED
    odp_plot.format_matplotlib()

    # adjust the bounding to fit text that went outside the limit of the plot
    plt.tight_layout()
    # safe make the directory
    os.makedirs(outdir, exist_ok=True)
    outprefix = "{}_and_{}_chromosome_conservation".format(sp1, sp2)
    outdir_prefix = os.path.join(outdir, outprefix)
    plt.savefig("{}.pdf".format(outdir_prefix), format='pdf')
    plt.close()

def main():
    # parse the arguments
    args = parse_args()

    # we must parse the config file to get the analysis parameters
    config = parse_config(args.config, args.directory, args.target_species)

    # build a sp_to_chr_to_size nested dictionary. We wouldn't need this if we finished clink.
    # nested for loop to get all the files for all the species
    rbh_filelist = set()
    for sp1 in config["analysis_files"].keys():
        for sp2 in config["analysis_files"][sp1].keys():
            rbh_filelist.add(config["analysis_files"][sp1][sp2])
    rbh_filelist = list(rbh_filelist)
    sp_to_chr_to_size, sp_to_scaf_to_genecount = rbh_files_to_sp_to_chr_to_size(rbh_filelist)
    ## print the target genecounts
    #print(sp_to_scaf_to_genecount[args.target_species])

    # For the target species, get a list of scaf names if they have more
    #  than 1% of the total amount of genes in the geneome.
    #  This is to avoid plotting the small scaffolds.
    target_keep_these_scafs_gt_one_percent_genes = {x:set() for x in config["target_species"]}
    for sp in config["target_species"]:
        total_genes = sum([sp_to_scaf_to_genecount[sp][scaf] for scaf in sp_to_scaf_to_genecount[sp]])
        for scaf in sp_to_scaf_to_genecount[sp]:
            if sp_to_scaf_to_genecount[sp][scaf] >= (0.01 * total_genes):
                target_keep_these_scafs_gt_one_percent_genes[sp].add(scaf)

    # safely make the directory called 'odp_pairwise_decay'
    os.makedirs("odp_pairwise_decay", exist_ok=True)
    # make a plot using the data
    for sp1 in config["target_species"]:
        outdir = os.path.join("odp_pairwise_decay", sp1)
        outdir = os.path.join(outdir, "decay_dataframes")
        # calculate the pairwise decay in chromosomes, save the files, get the list of files
        filestruct = calculate_pairwise_decay_sp1_vs_many(sp1, config, sp_to_chr_to_size,
                                                          target_keep_these_scafs_gt_one_percent_genes, outdir)

        # make the summary plot of all the chromosomes
        outdir = os.path.join("odp_pairwise_decay", sp1)
        outdir = os.path.join(outdir, "plot_overview_sp_sp")
        plot_pairwise_decay_sp1_vs_all(sp1, filestruct, outdir=outdir)
        #
        ### make individual sp-sp scatterplots
        #outdir = os.path.join("odp_pairwise_decay", sp1)
        #outdir = os.path.join(outdir, "plot_individual_sp_sp")
        #for sp2 in sorted(filestruct[sp1].keys()):
        #    plot_decay_twospecies(sp1, sp2, filestruct[sp1][sp2],
        #                          target_keep_these_scafs_gt_one_percent_genes[sp1], outdir)


if __name__ == '__main__':
    main()