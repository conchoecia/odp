#!/usr/bin/env python

"""
This contains the helper scripts to annotate the UMAP plots with the blast results.
There is a module, because this allows us to import the functions into other programs for testing.
"""

from ast import literal_eval
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import pandas as pd
import random

def umapdf_one_species_one_query(UMAPdf, blastp, analysis, ALG, n, m, query, outputPDF, species = None):
    """
    This function takes a dataframe as input and makes a UMAP
        basedir + "/blast_filt/{query}/{sample}_results.filt.blastp",
        UMAPdf     = basedir + "/ALG_reimbedding/{analysis}/{ALG}/{analysis}.neighbors_{n}.mind_{m}.missing_large.subchrom.df",
        algrbhfile = config["ALG_rbh_file"],

    """
    if species == None:
        raise IOError("Please use this function with a species. The inputs are different from the other function and this plots something in the context of one species. Exiting")
    #         ┓   ┓•┓
    # ┏┳┓┏┓╋┏┓┃┏┓╋┃┓┣┓
    # ┛┗┗┗┻┗┣┛┗┗┛┗┗┗┗┛
    #       ┛
    UMAPempty   = False
    BLASTPempty = False
    if os.stat(UMAPdf).st_size == 0:
        UMAPempty = True
    if os.stat(blastp).st_size == 0:
        BLASTPempty = True
    if UMAPempty or BLASTPempty:
        # make a 2in x 2in plot telling the user the message
        fig = plt.figure(figsize=(2,2))
        if UMAPempty and BLASTPempty:
            plt.text(0.5, 0.5, "The input df and the blastp file were empty, so we couldn't plot anything", fontsize = 3)
        elif UMAPempty and not BLASTPempty:
            plt.text(0.5, 0.5, "The input df was empty, so we couldn't plot anything", fontsize = 3)
        elif not UMAPempty and BLASTPempty:
            plt.text(0.5, 0.5, "The blastp file was empty, so we couldn't plot anything", fontsize = 3)
        plt.text(0.5, 0.5, "The input df was empty, so we couldn't plot anything", fontsize = 3)
        # make another line telling the user which file, exactly, was empty
        plt.text(0.5, 0.4, f"The input file was {UMAPdf}", fontsize = 3)
        # turn off the axis ticks
        for ax in fig.axes:
            ax.set_xticks([])
            ax.set_yticks([])
        # turn off the axes
        plt.axis('off')
        try:
            # make the plot tight to not cut off the text
            plt.tight_layout()
        except:
            pass
        # save the figure
        plt.savefig(outputPDF)
    else:
        dot = 3
        bigger_dot = 4
        # The embedding in this case only has the UMAP1 and UMAP2 columns.
        # We need later to go through the blast results to figure out which dots to annotate
        df_embedding = pd.read_csv(UMAPdf, sep = "\t", index_col = 0)
        print(df_embedding)
        # make a matplotlib plot of the UMAP with the df_embedding, and the color_dict from SplitLossColocTree as the legend
        # make a figure that is 5x5 inches
        fig, ax = plt.subplots(figsize=(2, 2))
        # scatter the UMAP1 and UMAP2 columns of the df_embedding
        ax.scatter(df_embedding["UMAP1"], df_embedding["UMAP2"], c = df_embedding["color"], s = dot)
        # now we read in the filtered blast results for this species
        blastp = pd.read_csv(blastp, sep = "\t", header = None)
        blastp.columns = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore", "scaf_of_nearest_ALG", "nearest_ALG", "position_of_nearest_ALG", "dist_to_nearest_ALG"]
        blastp["closest_ALG"] = None
        blastp["closest_ALG_dist"] = None
        blastp["closest_ALG_position"] = None
        # get the blastp rows with nearest_ALG values in the UMAPdf.index
        for i, row in blastp.iterrows():
            thisquery = row["qseqid"]
            nearest_list = literal_eval(row["nearest_ALG"])
            for j in range(len(nearest_list)):
                if nearest_list[j] in df_embedding.index:
                    print(f"found a match for ALG {ALG}, blastp query {thisquery} {nearest_list[j]} in the UMAPdf.index")
                    nearest_pos = literal_eval(row["position_of_nearest_ALG"])[j]
                    nearest_dist = literal_eval(row["dist_to_nearest_ALG"])[j]
                    blastp.at[i, "closest_ALG"]          = nearest_list[j]
                    blastp.at[i, "closest_ALG_dist"]     = nearest_dist
                    blastp.at[i, "closest_ALG_position"] = nearest_pos
                    break
        # remove rows where closest_ALG is None
        blastp = blastp[blastp["closest_ALG"].notnull()]

        # In many cases, the closest ALG is going to be the same for many of the blastp hits, so we should group these together
        gb = blastp.groupby("closest_ALG")
        # Go through each group, and color that dot in the UMAP
        legend_patches = []
        annotate_marker = '*'
        for thisALG, group in gb:
            # generate a random color for this dot
            thiscolor = "#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
            # label is all of the qseqids that are in this group
            label = ", ".join(group["qseqid"].values)
            legend_patches.append(  mlines.Line2D([], [], color=thiscolor, marker=annotate_marker, linestyle='None', markersize=dot*1.2, label=label))
            # x is the location of this ALG in the UMAPdf
            x = df_embedding.loc[thisALG, "UMAP1"]
            y = df_embedding.loc[thisALG, "UMAP2"]
            ax.scatter(x, y, c = thiscolor, s = dot, marker = annotate_marker, alpha = 1)

        # add the entries to the legend
        ax.legend(handles=legend_patches, loc="upper right", fontsize = 3)
        # Remove the frame from the legend
        ax.get_legend().set_frame_on(False)
        # turn off the axis ticks
        ax.set_xticks([])
        ax.set_yticks([])
        # set the title based on the input
        ax.set_title(f"UMAP of {analysis}, ALG {ALG},\ngenome {species} with {query} blast results,\nmin_dist {m}, n_neighbors {n}",
                   fontsize = 3)
        # save the figure
        plt.savefig(outputPDF, bbox_inches='tight')
