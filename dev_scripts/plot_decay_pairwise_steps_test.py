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

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import os
import random
import sys
import yaml

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
    args = parser.parse_args()

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

def parse_config(config_file, directory_of_rbh_files):
    """
    This will be a yaml file with analyses, detailing the target species to plot,
      the other species, and the divergence times.

    """
    species_prefix_to_filename = {}
    filelist = os.listdir(directory_of_rbh_files) # get a list of files in the directory
    config = read_yaml_file(config_file)          # load the config file
    # check if the necessary field, analyses, is in the file
    if "analyses" not in config:
        raise IOError("The field 'analyses' must be present in the config file.")

    # FILE PAIRING
    # for each analysis we now need to get the file to look at for the pairwise comparison
    if "analysis_files" not in config:
        config["analysis_files"] = {x:{} for x in config["analyses"].keys()}

    # get the pairs from the files
    pair_to_file = {}
    for thisfile in [x for x in filelist if x.endswith(".rbh")]:
        fields = thisfile.split("_")
        analysis_pair = tuple(sorted((fields[0], fields[1])))
        if analysis_pair not in pair_to_file:
            pair_to_file[analysis_pair] = []
        complete_filepath = os.path.join(directory_of_rbh_files, thisfile)
        pair_to_file[analysis_pair].append(complete_filepath)

    # go through the files and figure out the pair
    for sp1 in config["analyses"].keys():
        for sp2 in config["analyses"][sp1].keys():
            analysis_pair = tuple(sorted((sp1, sp2)))
            if analysis_pair not in pair_to_file:
                raise IOError("The analysis pair {} is not in the filelist.".format(analysis_pair))
            config["analysis_files"][sp1][sp2] = pair_to_file[analysis_pair][0]

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
                      "sp1_scaf_genecount": sp1_scaf_to_total_genes[sp1_scaf],
                      "conserved":          0,
                      "scattered":          sp1_scaf_to_total_genes[sp1_scaf]
                      }
        if sp1_scaf in sp1_scaf_to_sp2_scaf.keys():
            thisentry["sp2_scaf"]           = sp1_scaf_to_sp2_scaf[sp1_scaf]
            thisentry["sp1_scaf_genecount"] = sp1_scaf_to_total_genes[sp1_scaf]
            thisentry["conserved"]          = sp1_scaf_to_conserved_genes[sp1_scaf]
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
    return sp_to_chr_to_size

def plot_pairwise_decay_single_species(sp1, config, sp_to_chr_to_size):
    """
    Outputs a plot wherein each source species has its own figure.
    There will be two subplots per figures.
    The left subplot will have the whole-genome conservation vs divergence time.
    The right subplot will have the per-chromosome conservation vs divergence time.
    """
    # We start by making a subplot array to make the two variants. A two-column, one-row plot.
    NUMBER_OF_ROWS = 1
    NUMBER_OF_FIGS = 2
    fig, axes = plt.subplots(NUMBER_OF_ROWS, NUMBER_OF_FIGS, figsize = (7.5 * NUMBER_OF_FIGS, 6*NUMBER_OF_ROWS))
    fig.suptitle("{} decay versus divergence time".format(sp1))

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

        print(sp2, "\n", sp1_sp2_decay)
        print("ABOUT TO CRASH")

        # Add a percent conserved column
        sp1_sp2_decay["fraction_conserved"] = sp1_sp2_decay["conserved"] / sp1_sp2_decay["sp1_scaf_genecount"]

        # Make a whole-genome version of the dataframe. Just sum up the columns and recalculate the percent conserved
        sp1_sp2_whole = sp1_sp2_decay.sum(axis=0).to_frame().transpose()
        # only keep certain columns
        sp1_sp2_whole = sp1_sp2_whole[["sp1_scaf_genecount", "conserved", "scattered"]]
        sp1_sp2_whole["divergence_time"] = config["analyses"][sp1][sp2]
        sp1_sp2_whole["fraction_conserved"] = sp1_sp2_whole["conserved"] / sp1_sp2_whole["sp1_scaf_genecount"]

        print(sp1_sp2_decay)
        print(sp1_sp2_whole)

        # on the left-plot just do a scatterplot of the fraction conserved vs divergence time
        axes[0].scatter(sp1_sp2_whole["divergence_time"], sp1_sp2_whole["fraction_conserved"], label = "{}".format(sp2))
        axes[0].set_xlabel("Divergence time (MYA)")
        axes[0].set_ylabel("Fraction conserved on orthologous chromosomes")
        axes[0].set_title("Whole-genome conservation vs divergence time")

        # the right plot is per-chromosome. Add a little jitter to the x-axis so we can see the points
        axes[1].scatter(jitter(sp1_sp2_decay["divergence_time"], 20), sp1_sp2_decay["fraction_conserved"], label = "{}".format(sp2), alpha = 0.25)
        axes[1].set_xlabel("Divergence time (MYA) (+- 20 MYA jitter)")
        axes[1].set_ylabel("Fraction conserved on orthologous chromosomes")
        axes[1].set_title("Orthologous chromosome conservation vs divergence time")



    # Save the plot as a jpeg
    outprefix = "{}_decay_plot_vs_divergence_time".format(sp1)
    plt.savefig("{}.jpg".format(outprefix), format='jpeg')


def main():
    # parse the arguments
    args = parse_args()

    # we must parse the config file to get the analysis parameters
    config = parse_config(args.config, args.directory)

    # build a sp_to_chr_to_size nested dictionary. We wouldn't need this if we finished clink.
    # nested for loop to get all the files for all the species
    rbh_filelist = set()
    for sp1 in config["analysis_files"].keys():
        for sp2 in config["analysis_files"][sp1].keys():
            rbh_filelist.add(config["analysis_files"][sp1][sp2])
    rbh_filelist = list(rbh_filelist)
    sp_to_chr_to_size = rbh_files_to_sp_to_chr_to_size(rbh_filelist)

    # make a plot using the data
    for sp1 in config["analyses"].keys():
        plot_pairwise_decay_single_species(sp1, config, sp_to_chr_to_size)

if __name__ == '__main__':
    main()