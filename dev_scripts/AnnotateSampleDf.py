#!/usr/bin/env python

"""
This program contains the functions used by AnnotateSampleDf.snakefile
"""

# This block imports fasta-parser as fasta
from ast import literal_eval
import os
import sys
thispath = os.path.dirname(os.path.realpath(__file__))
dependencies_path = os.path.join(thispath, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

import matplotlib.pyplot as plt
import odp_plotting_functions as odp_plot
import pandas as pd
from rbh_tools import parse_rbh

def gen_rbh_stats(samplerbhfilepath, algrbhfilepath, ALGname, outfilepath):
    """
    This function generates the stats of an rbh file - namely the dispersion.
    Things that are calculated for this are:
     - the number of proteins in the rbh file
     - the number of gene groups in the rbh file
     - the number of genes for each gene group

    Input:
      - It takes in a single argument, the path to one rbh file.
    Output:
      - The output is a text file that contains the analysis as key: value pairs.
      - The fields that are output are:
        - frac_ologs:           The fraction of genes of ANY ALG that are present at all in the rbh file. len(rbhdf) / total_genes_ALGs
        - frac_ologs_sig:       The fraction of genes of ANY ALG that are significantly on any chromosome, as defined by whole_FET
        - frac_ologs_single:    The fraction of genes of ANY ALG that are significantly on the largest chromosome, as defined by whole_FET
        - frac_ologs_{ALGNAME}: The fraction of genes of INDIVIDUAL ALGs that are significantly on any chromosome
    """
    ALGrbhdf = parse_rbh(algrbhfilepath)
    total_genes_ALGs = len(ALGrbhdf)
    genes_per_ALG    = ALGrbhdf.groupby("gene_group").size().to_dict()

    # now parse the rbh file
    rbhdf = parse_rbh(samplerbhfilepath)
    sigdf = rbhdf[rbhdf["whole_FET"] < 0.05]

    # get the sample scaf column
    sample_loc_col = [col for col in rbhdf.columns if (col.endswith("_scaf")) and (ALGname not in col)][0]
    frac_ologs = len(rbhdf)/total_genes_ALGs
    frac_ologs_sig = len(rbhdf[rbhdf["whole_FET"] < 0.05])/total_genes_ALGs
    # groupby the gene_group, then get the rows with the most frequent sample_loc_col value
    entries = []
    for gene_group, groupdf in sigdf.groupby("gene_group"):
        max_single = 0 if len(groupdf) == 0 else groupdf[sample_loc_col].value_counts().max()
        entries.append({"gene_group":     gene_group,
                        "max_single":     max_single,
                        "genes_in_group": genes_per_ALG[gene_group]})
    # frac_ologs_single is the sum o
    if len(entries) == 0:
        frac_ologs_single = float(0)
    else:
        frac_ologs_single = pd.DataFrame(entries)["max_single"].sum() / total_genes_ALGs
    # print all of these fields to a text file in the format:
    # field: value
    with open(outfilepath, "w") as f:
        f.write(f"frac_ologs: {frac_ologs}\n")
        f.write(f"frac_ologs_sig: {frac_ologs_sig}\n")
        f.write(f"frac_ologs_single: {frac_ologs_single}\n")
        for ALG in genes_per_ALG:
            f.write(f"frac_ologs_{ALG}: {len(sigdf[sigdf['gene_group'] == ALG])/genes_per_ALG[ALG]}\n")

def gen_annotation_stats(sampleproteinfilepath, algrbhfilepath, outfilepath):
    """
    This generates information about the genome annotation. Specifically, it looks at the proteins in the annotation.
    Things that are calculated:
      - the number of proteins
      - the mean protein length
      - the median protein length
      - the longest protein
      - the smallest protein
      - whether the proteins are from a real annotation or from the RBH entries
    Input:
      - algrbhfile
      - protein fasta.gz file
    output:
      - A text file with the following fields:
        - num_proteins: {num_proteins}
        - mean_protein_length: {mean_protein_length}
        - median_protein_length: {median_protein_length}
        - longest_protein: {longest_protein}
        - smallest_protein: {smallest_protein}
        - from_rbh: {from_rbh}
    """
    # read in the ALG_rbh file as a pandas df
    df = parse_rbh(algrbhfilepath)
    rbh_names = list(df["rbh"])
    # read in the proteins. Make a list of putative rbh proteins. Get the other stats.
    entries = []
    for record in fasta.parse(sampleproteinfilepath):
        entries.append({"protname":           record.id,
                        "protlen" :           len(record.seq),
                        "putative_rbh_name" : "_".join(record.id.split("_")[:-1]) if "_" in record.id else record.id })
    protdf = pd.DataFrame(entries)
    num_proteins = len(protdf)
    mean_protein_length = protdf["protlen"].mean()
    median_protein_length = protdf["protlen"].median()
    longest_protein  = protdf["protlen"].max()
    smallest_protein = protdf["protlen"].min()
    # count the number of times the putative_rbh_name is in the rbh_names
    if protdf["putative_rbh_name"].isin(rbh_names).sum() > (0.25 * len(protdf)):
        from_rbh = True
    else:
        from_rbh = False
    # Print all of these fields to a text file in the format:
    # field: value
    with open(outfilepath, "w") as f:
        f.write(f"num_proteins: {num_proteins}\n")
        f.write(f"mean_protein_length: {mean_protein_length}\n")
        f.write(f"median_protein_length: {median_protein_length}\n")
        f.write(f"longest_protein: {longest_protein}\n")
        f.write(f"smallest_protein: {smallest_protein}\n")
        f.write(f"from_rbh: {from_rbh}\n")

def gen_genome_stats(genomefilepath, outfilepath):
    """
    This generates statistics about a genome assembly.
    The stats that are output are:
     - the number of scaffolds
     - the GC content
     - the genome size
     - the median scaffold length
     - the mean scaffold length
     - scaffold N50
     - longest scaffold
     - smallest scaffold
     - percent Ns

    The parameters are:
      - genomefilepath: The path to the genome fasta file
      - outfilepath: The path to the output file
    The output:
      - a key: value text file with the above fields.
    """
    entries= []
    for record in fasta.parse(genomefilepath):
        entries.append({"scafname": record.id,
                        "scaflen" : len(record.seq),
                        "gc" : (record.seq.count("G") + record.seq.count("C")) / len(record.seq),
                        "Ns" : record.seq.count("N"),
                        # gaps are the number of sequential Ns of length 10 or more
                        "num_gaps": len(record.seq.upper().split("NNNNNNNNNN")) - 1})

    # make a dataframe from the entries
    df = pd.DataFrame(entries)

    num_scaffolds = len(df)
    GC_content = df["gc"].mean()
    genome_size = df["scaflen"].sum()
    median_scaffold_length = df["scaflen"].median()
    mean_scaffold_length = df["scaflen"].mean()
    scaffold_N50 = df["scaflen"].sort_values(ascending=False).cumsum().searchsorted(genome_size/2)
    longest_scaffold  = df["scaflen"].max()
    smallest_scaffold = df["scaflen"].min()
    fraction_Ns = df["Ns"].sum() / genome_size
    number_of_gaps = df["num_gaps"].sum()
    # print all of these fields to a text file in the format:
    # field: value
    with open(outfilepath, "w") as f:
        f.write(f"num_scaffolds: {num_scaffolds}\n")
        f.write(f"GC_content: {GC_content}\n")
        f.write(f"genome_size: {genome_size}\n")
        f.write(f"median_scaffold_length: {median_scaffold_length}\n")
        f.write(f"mean_scaffold_length: {mean_scaffold_length}\n")
        f.write(f"scaffold_N50: {scaffold_N50}\n")
        f.write(f"longest_scaffold: {longest_scaffold}\n")
        f.write(f"smallest_scaffold: {smallest_scaffold}\n")
        f.write(f"fraction_Ns: {fraction_Ns}\n")
        f.write(f"number_of_gaps: {number_of_gaps}\n")

def stats_filepath_to_dict(stats_filepath):
    """
    This reads in a stats file and returns a dictionary of the stats
    """
    entries = {}
    with open(stats_filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                key, value = line.split(": ")
                # check if the value can be cast to a float
                if value.replace(".","").isdigit():
                    if "." in value:
                        value = float(value)
                    else:
                        value = int(value)
                else:
                    # check if it is a boolean
                    if value in ["True", "False"]:
                        # we have to do this because bool("False") evaluates to True
                        if value == "True":
                            value = True
                        elif value == "False":
                            value = False
                entries[key] = value
    return entries

def plot_decay(ax, df,
               alg_to_size_dict,
               x_axis_absolute_or_ranked,
               plot_title = "",
               ymin = -1, ymax = 1,
               y_axis_label = "y_axis_label_me",
               x_axis_label = "x_axis_label_me",
               x_axis_labels_on = True,
               x_axis_labels_are_ALGs = False,
               plot_bars = False):
    """
    Makes a plot of the decay values, with the values based on difference from median.

    The y-axis is linear
    """
    if x_axis_absolute_or_ranked not in ["absolute", "ranked", "boxplot"]:
        raise ValueError("x_axis_absolute_or_ranked must be either 'absolute' or 'ranked', or 'boxplot'")
    # check that alg_to_size_dict is not empty
    if len(alg_to_size_dict) == 0:
        raise IOError("alg_to_size_dict must not be empty.")

    # figure out the plot order based on size. These are the values.
    ALG_order = dict(sorted(alg_to_size_dict.items(), key=lambda item: item[1]))

    # look in the df_dict and check that there is a frac_ologs_{ALG} column for each ALG.
    missing_columns = []
    for ALGname in alg_to_size_dict:
        if f"frac_ologs_{ALGname}" not in df.columns:
            missing_columns.append(f"frac_ologs_{ALGname}")
    if len(missing_columns) > 0:
        raise ValueError("The following columns are missing from the df_dict: {}".format(missing_columns))

    # The x-axis will be shared by both axes, so just get these once
    if x_axis_absolute_or_ranked == "ranked":
        x = list(range(len(ALG_order)))
    elif x_axis_absolute_or_ranked == "absolute":
        x = [alg_to_size_dict[thisALG] for thisALG in ALG_order]
        ax.set_xlabel('ALG size (total genes)')

    # This is the same for both plots
    if plot_bars and (x_axis_absolute_or_ranked != "boxplot"):
        bar_alpha = 0.15
        # In this scenario we want to plot a bar graph as a visual aide
        ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
        bar_ymin = 0
        bar_ymax = max(alg_to_size_dict.values()) * 1.05
        ax2.set_ylim([bar_ymin, bar_ymax])

        color = 'tab:blue'
        ax2.set_ylabel('ALG size', color=color)  # we already handled the x-label with ax1
        # we already defined x,  above
        y = [alg_to_size_dict[thisALG] for thisALG in ALG_order]
        # make a plot with single unconnected blue dots
        ax2.plot(x, y, 'o', color=color, alpha = bar_alpha)
        # make a vertical line up to the dots
        for i in range(len(x)):
            ax2.plot([x[i], x[i]], [0, y[i]], color=color, alpha = bar_alpha)
        ax2.tick_params(axis='y', labelcolor=color)

    # this plots    # Loop through each line and plot
    # iterate through each sample in the df_dict row-wise
    linealpha = max(0.01, 1/len(df))
    i = 0
    if x_axis_absolute_or_ranked == "boxplot":
        for thisALG in sorted(ALG_order, key=ALG_order.get):
            # get the values for this ALG
            y = [row[f"frac_ologs_{thisALG}"] for i, row in df.iterrows()]
            # make a box plot of each line
            flierprops = dict(marker='.', markersize=5, markerfacecolor='black', markeredgecolor='black', alpha=0.05)
            ax.boxplot(y, positions=[i], widths=0.5, patch_artist=True, flierprops = flierprops)
            ax.set_ylabel('difference from median of fraction conserved')
            i += 1
    else:
        for i, row in df.iterrows():
            # we need to change the x-axis to be ranked dependent on the size of the ALG
            if x_axis_absolute_or_ranked == "ranked":
                ax.set_xlabel('ALG size (ranked)')
            elif x_axis_absolute_or_ranked == "absolute":
                ax.set_xlabel('ALG size (total genes)')
            # We have already defined x
            y = [row[f"frac_ologs_{thisALG}"] for thisALG in ALG_order]
            # make a line plot of each line
            ax.plot(x, y,
                     alpha=linealpha, color ="black", lw = 1.0)
            ax.set_ylabel('difference from median of fraction conserved')

    # turn off the text of the xaxis if we must
    if not x_axis_labels_on:
        ax.set_xticklabels([])

    # If we want to label the x-axis with the ALG ids.
    # To do this we must have the alg_to_size_dict
    if x_axis_labels_are_ALGs:
        if len(alg_to_size_dict) == 0:
            raise Exception("We must have something in alg_to_size_dict to plot as the axis labels.")
        # get the sorted x ALG sizes as a list, and also save the labels to another list
        alg_ID_smallest_to_largest    = sorted(ALG_order, key=ALG_order.get)
        alg_value_smallest_to_largest = [ALG_order[x] for x in alg_ID_smallest_to_largest]
        alg_index = range(len(alg_ID_smallest_to_largest))

        index_or_absolute = None
        if x_axis_absolute_or_ranked in ["boxplot", "ranked"]:
            index_or_absolute = "index"
        elif x_axis_absolute_or_ranked in ["absolute"]:
            index_or_absolute = "absolute"
        else:
            raise ValueError("x_axis_absolute_or_ranked must be either 'absolute' or 'ranked', or 'boxplot")
        if index_or_absolute == "index":
            ax.set_xticks(alg_index)
        elif index_or_absolute == "absolute":
            # now set the tick positions and their labels
            ax.set_xticks(alg_value_smallest_to_largest)
        else:
            raise ValueError("index_or_absolute must be 'index' or 'absolute'")

        # make the text a little smaller and the horizontal alignment to right
        ax.set_xticklabels(alg_ID_smallest_to_largest, rotation=40, ha='center')

    ax.set_ylim(ymin, ymax)

    # set the titles for this ax object
    ax.set_title(plot_title)

    return ax

def bin_and_plot_decay(ALGrbhdf_filepath: str, algrbhrbhstatsdf_filepath: str, outpdf_filepath: str,
                       ALGname: str, num_bins: int = 5):
    """
    This function is derived from "fig2" in plot_decay_many_species.py.
    This takes in a rbhstats dataframe.
    The dataframe has the columns:
      - sample: the name of the same that we are analyzing
      - frac_ologs: the fraction of orthologs in the sample.
      - frac_ologs_sig: the fraction of orthologs that are significantly on any chromosome
      - frac_ologs_single: the fraction of orthologs that are significantly on the largest chromosome
      - frac_ologs_{ALG}: the fraction of orthologs that are significantly on any chromosome for a specific ALG

    This script works by binning the samples based on frac_ologs, and then plots the results based on frac_ologs_{ALG}
    """
    for checkthis in [ALGrbhdf_filepath, algrbhrbhstatsdf_filepath]:
        if not os.path.isfile(checkthis):
            raise ValueError(f"File {checkthis} does not exist.")
    # check that the type of ALGname is a string
    if not isinstance(ALGname, str):
        raise ValueError("ALGname must be a string")
    # check that the type of num_bins is an int
    if not isinstance(num_bins, int):
        raise ValueError("num_bins must be an int")

    # generate a alg_to_size_dict dictionary from the ALGrbhdf
    ALGrbhdf         = parse_rbh(ALGrbhdf_filepath)
    alg_to_size_dict = ALGrbhdf.groupby("gene_group").size().to_dict()

    # Each column of plotting is a single bin with a minimum and maximum cutoff for the ALG conservation values
    # We must label each column with the bins.
    # COORDINATE SYSTEM:
    #  x-axis: Bins in decreasing value, can access
    print("The number of bins are: {}".format(num_bins))
    NUMBER_OF_BINS = num_bins
    BINRANGES      = [(i / NUMBER_OF_BINS, (i + 1) / NUMBER_OF_BINS) for i in range(NUMBER_OF_BINS)][::-1]
    NUMBER_OF_FIGS = 4
    fig, axes = plt.subplots(NUMBER_OF_FIGS, NUMBER_OF_BINS, figsize = (5 * NUMBER_OF_BINS, 4*NUMBER_OF_FIGS))
    fig.suptitle(f"% of orthologs conserved (median) on {ALGname} ALGs in many animals.")

    plotdf = pd.read_csv(algrbhrbhstatsdf_filepath, sep="\t")

    # now we make the individual plots for each cutoff
    for i in range(NUMBER_OF_BINS):
        # add some vertical and horizontal space between the plots
        fig.subplots_adjust(hspace=0.5)
        fig.subplots_adjust(wspace=0.5)

        print("We are in bin {}. Bin values are {}".format(i, BINRANGES))
        # Get the samples that fall into the range of the bin.
        # We use "frac_ologs_sig", because this is the fraction of orthologs that are significantly on any chromosome
        if i == NUMBER_OF_BINS-1:
            plotdfs = plotdf[(plotdf["frac_ologs_sig"] >= BINRANGES[i][0]) & (plotdf["frac_ologs_sig"] <= BINRANGES[i][1])]
        else:
            plotdfs = plotdf[(plotdf["frac_ologs_sig"] >= BINRANGES[i][0]) & (plotdf["frac_ologs_sig"] < BINRANGES[i][1])]
        print("The number of samples in this bin is: {}".format(len(plotdfs)))
        if len(plotdfs) > 0:
            axes[0, i] = plot_decay( axes[0, i], plotdfs,
                                    alg_to_size_dict,
                                    "absolute",
                                    plot_title = "{}%-{}% conserved".format(BINRANGES[i][0], BINRANGES[i][1]),
                                    ymin = -0.05, ymax = 1.05,
                                    y_axis_label = "fraction conserved on ALG",
                                    x_axis_label = "ALG size based on genes found in blast",
                                    plot_bars = True)
            axes[1, i] = plot_decay( axes[1, i], plotdfs,
                                    alg_to_size_dict,
                                    "ranked",
                                    plot_title = "{}%-{}% conserved".format(BINRANGES[i][0], BINRANGES[i][1]),
                                    ymin = -0.05, ymax = 1.05,
                                    y_axis_label = "fraction conserved on ALG",
                                    x_axis_label = "Rank of absolute ALG size, smallest to largest",
                                    plot_bars = True,
                                    x_axis_labels_on = True,
                                    x_axis_labels_are_ALGs = True)
            axes[2, i] = plot_decay( axes[2, i], plotdfs,
                                    alg_to_size_dict,
                                    "boxplot",
                                    plot_title = "{}%-{}% conserved".format(BINRANGES[i][0], BINRANGES[i][1]),
                                    ymin = -0.05, ymax = 1.05,
                                    y_axis_label = "fraction conserved on ALG",
                                    x_axis_label = "Rank of absolute ALG size, smallest to largest",
                                    plot_bars = True,
                                    x_axis_labels_on = True,
                                    x_axis_labels_are_ALGs = True)

            # for the last axis object, we want to use that space to add a legend of the file names for this bin
            #  we will use the axes[2, i] object to do this
            text_to_plot = list(sorted(plotdfs["sample"].values))
            text_ymin = 1
            for j in range(len(text_to_plot)):
                # the text position needs to change in each loop
                text_position = (0.0, 0.95 - (j * 0.06))
                text_ymin = text_position[1]
                axes[3, i].text(text_position[0], text_position[1], text_to_plot[j], transform=axes[3, i].transAxes)
            # Plot title is the species that have occurred in this bin
            axes[3, i].set_title("Species in this bin")

            # set the y-axis min and max now
            axes[3, i].set_ylim(text_ymin - 0.05, 1.05)
            ## now turn off the lines and ticks for this axis
            axes[3, i].axis('off')


    # CALL THIS TO GET THE VISUAL STYLE WE NEED
    odp_plot.format_matplotlib()
    # Save the plot as a jpeg
    plt.savefig(outpdf_filepath, format='pdf')

def taxid_list_include_exclude_boolean(taxid_list, taxids_to_include, taxids_to_exclude):
    """
    This function returns True if the taxid is to be included in the plot.
    Returns False if it is to be excluded.
    This is a helper function that will be called from plot_UMAP_highlight_subclade() in a lambda function.
    """
    # make sure the taxid_list is a list, if not, interpret it as a literal
    if not isinstance(taxid_list, list):
        taxid_list = literal_eval(taxid_list)
    # logic for checking whether to include
    taxid_list = set(taxid_list)
    include_intersection = taxid_list.intersection(set(taxids_to_include))
    exclude_intersection = taxid_list.intersection(set(taxids_to_exclude))
    return (len(include_intersection) > 0) and (len(exclude_intersection) == 0)

def plot_UMAP_highlight_subclade(df_filepath, title, taxids_to_include,
                                 taxids_to_exclude, pdfout,
                                 gray_color = "#D6D9E0") -> None:
    """
    This plot shows the points of a UMAP, and highlights species in a subclade. Non-highlighted species are plotted in gray.

    Inputs:
      - df_filepath: The path to the dataframe output from the rule make_composite_dataframe from the snakefile.
      - title: The title of the plot that will be put on the pdf.
      - taxids_to_include: A list of taxids that will be highlighted in the plot.
                            If the taxid string has this integer in it, it will be highlighted.
      - taxids_to_exclude: A list of taxids that will be excluded from the plot.
                            If the taxid string has this integer in it, it will be excluded,
                            even if it is not in the taxids_to_include list.
      - pdfout: The path to the output pdf file.
      - gray_color: The color of the points that are not in the taxids_to_include list.
    """
    # enforce the types of the inputs
    if not isinstance(df_filepath, str):
        raise ValueError("df_filepath must be a string.")
    if not isinstance(title, str):
        raise ValueError("title must be a string.")
    if not isinstance(taxids_to_include, list):
        raise ValueError("taxids_to_include must be a list.")
    if not isinstance(taxids_to_exclude, list):
        raise ValueError("taxids_to_exclude must be a list.")
    if not isinstance(pdfout, str):
        raise ValueError("pdfout must be a string.")
    # ensure that df_filepath exists
    if not os.path.isfile(df_filepath):
        raise ValueError(f"File {df_filepath} does not exist.")
    for thislist in [taxids_to_include, taxids_to_exclude]:
        if not all(isinstance(x, int) for x in thislist):
            raise ValueError("All elements of taxids_to_include and taxids_to_exclude must be integers.")

    # read in the dataframe
    df = pd.read_csv(df_filepath, sep="\t")
    print(df)
    # change the color values of the points based on the taxids
    df["color"] = df.apply(lambda x:
                    x["color"] if taxid_list_include_exclude_boolean(x["taxid_list"],
                        taxids_to_include, taxids_to_exclude)
                    else gray_color, axis=1)
    # make a plot with a single panel that is 1x1 inches.
    fig = plt.figure(figsize=(2, 2))

    # scatter the UMAP1 and UMAP2 columns of the df_embedding
    plt.scatter(df["UMAP1"], df["UMAP2"],
                c = list(df["color"]), s = 1,
                linewidths = 0, edgecolors=None)

    # Set the title
    plt.title(title, fontsize=6)

    # turn off the x and y axis ticks and tick values
    plt.xticks([])
    plt.yticks([])
    # save the figure
    #print(pdfout)
    plt.savefig(pdfout, format='pdf')
    plt.close()

