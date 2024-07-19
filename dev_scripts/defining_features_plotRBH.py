#!/usr/bin/env

"""
For this plot, we will plot an oxford dot plot of the RBHs, we will also use the
 unique pairs file to plot the strongest connections between the genomes. We will
 use arcs to draw those.

# TODO - make the arcs not change in a gradient, but instead be colored by three arcs, left, right, middle(black)
"""

import argparse
import ast
import ete3
import numpy as np
import os
import pandas as pd
import random
import rbh_tools
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from   matplotlib.path import Path
import matplotlib.path as mpath
import matplotlib.patches as patches
from   matplotlib.patches import Rectangle

import odp_plotting_functions as odp_plot


# we want to parse a genome fasta file
thispath = os.path.dirname(os.path.realpath(__file__))
dependencies_path = os.path.join(thispath, "../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

def parse_args():
    """
    We need the following files:
      - a plottable rbh file - this is from a specific species rather than being from a 
      - a unique pairs tsv file
      - a genome assembly file
      - the taxid of the clade to plot
      - an optional blast file that will will use to color the arcs
      - a window range to look in, type int, unit is basepairs
    """
    parser = argparse.ArgumentParser(description="Plot RBHs")
    parser.add_argument("--rbh_file",          type=str, help="Path to the rbh file", required=True)
    parser.add_argument("--unique_pairs_path", type=str, help="Path to the unique pairs file", required=True)
    parser.add_argument("--clade_for_pairs",   type=int, help="The taxid of the clade to plot", required=True)
    parser.add_argument("--genome_file",       type=str, help="Path to the genome file", required=True)
    parser.add_argument("--blast_file",        type=str, help="Path to the blast file to color the arcs", required=False)
    parser.add_argument("--chrom_file",        type=str, help="Optional path to the chrom file", required=False)
    parser.add_argument("--window",            type=int, help="The window size to look for blast results", default=1000000)
    args = parser.parse_args()
    # check that the files exist
    if not os.path.exists(args.rbh_file):
        raise FileNotFoundError(f"Could not find the rbh file at {args.rbh_file}")
    if not os.path.exists(args.unique_pairs_path):
        raise FileNotFoundError(f"Could not find the unique pairs file at {args.unique_pairs_path}")
    # check if the blast file exists if the user has specified it
    if args.blast_file is not None:
        if not os.path.exists(args.blast_file):
            raise FileNotFoundError(f"Could not find the blast file at {args.blast_file}")
    if args.chrom_file is not None:
        if not os.path.exists(args.chrom_file):
            raise FileNotFoundError(f"Could not find the chrom file at {args.chrom_file}")
    return args

def interpolate_color(start_color, end_color, factor):
    """
    Interpolate between two colors.
    """
    start_rgba = np.array(to_rgba(start_color))
    end_rgba = np.array(to_rgba(end_color))
    return start_rgba + (end_rgba - start_rgba) * factor

def de_casteljau(P0, P1, P2, P3, t):
    """ De Casteljau's algorithm to split a cubic Bezier curve at t """
    P01 = (1 - t) * P0 + t * P1
    P12 = (1 - t) * P1 + t * P2
    P23 = (1 - t) * P2 + t * P3
    P012 = (1 - t) * P01 + t * P12
    P123 = (1 - t) * P12 + t * P23
    P0123 = (1 - t) * P012 + t * P123
    return P01, P12, P23, P012, P123, P0123

def plot_bezier_arc(panel, startx, starty, stopx, stopy, height, color, alpha,
                    color_gradient = False, lw = 0.25):
    """
    Plot bezier curves between chromosome coordinates of different species.
    - This plots a bezier between two points, centered on a midline

    The optional argument, color_gradient, will make the color of the arc change from the start to the stop
    """
    if color_gradient == False:
        start_color = color
        end_color   = color
    else:
        start_color = color_gradient[0]
        end_color = color_gradient[1]

    indent = 1.0
    # plot the indices
    diffx = abs(startx - stopx)
    middlex = min([startx, stopx]) + (diffx/2)

    second            = (startx,  height)
    second_point_five = (middlex, height)
    third             = (stopx,   height)
    path_data = [
        (Path.MOVETO, (startx, starty)),
        (Path.CURVE4, second),
        (Path.CURVE4, third),
        (Path.CURVE4, (stopx, stopy)),
        ]
    codes, verts = zip(*path_data)
    path  = mpath.Path(verts, codes)
    if start_color == end_color:
        # In this case, we plot a black line between the two points.
        # This is the easiest to work with later as a vector image
        if start_color == "#000000":
            zord = -50
        elif alpha < 0.5:
            zord = -99
        else:
            zord = 1
        patch = patches.PathPatch(path, fill = False,
                                  facecolor=start_color, lw = lw,
                                  alpha=alpha, edgecolor = start_color,
                                  zorder = zord)
        panel.add_patch(patch)
    else:
        # In this case, we plot an arc that is a gradient between the two colors
        # how many segments to use to change the color?
        # Previously we used a graduent, but now we don't like that because it makes too many lines. Now, just use 3 max
        P0 = np.array([startx, starty])
        P1 = np.array([startx, height])
        P2 = np.array([stopx,  height])
        P3 = np.array([stopx,  stopy])

        # Plot first half of the bezier arc in the start color
        # Split the bezier curve at t = 0.5
        P01, P12, P23, P012, P123, P0123 = de_casteljau(P0, P1, P2, P3, 0.5)

        # Plot the left half of the Bezier arc
        left_path_data = [
            (mpath.Path.MOVETO, P0),
            (mpath.Path.CURVE4, P01),
            (mpath.Path.CURVE4, P012),
            (mpath.Path.CURVE4, P0123),
        ]
        left_codes, left_verts = zip(*left_path_data)
        left_path = mpath.Path(left_verts, left_codes)
        left_patch = patches.PathPatch(left_path, fill=False, lw=lw, alpha=alpha, edgecolor=start_color)
        panel.add_patch(left_patch)

        # Plot the right half of the Bezier arc
        right_path_data = [
            (mpath.Path.MOVETO, P0123),
            (mpath.Path.CURVE4, P123),
            (mpath.Path.CURVE4, P23),
            (mpath.Path.CURVE4, P3),
        ]
        right_codes, right_verts = zip(*right_path_data)
        right_path = mpath.Path(right_verts, right_codes)
        right_patch = patches.PathPatch(right_path, fill=False, lw=lw, alpha=alpha, edgecolor=end_color)
        panel.add_patch(right_patch)

        # Add a thicker black "handle" in the middle of the arc
        # Split the bezier curve at t = 0.42 and t = 0.58 for the handle
        P01_42, P12_42, P23_42, P012_42, P123_42, P0123_42 = de_casteljau(P0, P1, P2, P3, 0.42)
        P01_58, P12_58, P23_58, P012_58, P123_58, P0123_58 = de_casteljau(P0, P1, P2, P3, 0.58)

        handle_path_data = [
            (mpath.Path.MOVETO, P0123_42),
            (mpath.Path.CURVE4, P123_42),
            (mpath.Path.CURVE4, P012_58),
            (mpath.Path.CURVE4, P0123_58),
        ]
        handle_codes, handle_verts = zip(*handle_path_data)
        handle_path = mpath.Path(handle_verts, handle_codes)
        handle_patch = patches.PathPatch(handle_path, fill=False, lw=lw * 2, alpha=alpha, edgecolor="#000000")
        panel.add_patch(handle_patch)

    return panel

def legal_color(color):
    """
    Check
    - the type is string
    - the first character is #
    - the length is 7
    - the rest of the characters are 0-9 or A-F
    """
    if not isinstance(color, str):
        return False
    if color[0] != "#":
        return False
    if len(color) != 7:
        return False
    for char in color[1:]:
        if char not in "0123456789ABCDEF":
            return False
    return True

def plot_arcs(arc_panel, genome1, pairtype, rbh, unique_pairs, fontsize,
              blastdf = None, blast_rbh_list = None, chromdf = None,
              scaf_order = None, scaf_to_len = None, window = 1000000,
              attempt_to_gradient_color = False):
    """
    This takes a panel
      - (probably one of the arc panel)
      - a rbh df
      - a unique pairs df
      - an optional blastdf that can be interpreted to color the arcs
    """
    if attempt_to_gradient_color is True:
        # we need to check if there is a color column in the blastdf
        # if there is, make a rbh_to_color dict
        if not "color" in rbh.columns:
            raise ValueError("You must have a color column in the blastdf to use the attempt_to_gradient_color option")
        rbh_to_color = {}
        for i, row in rbh.iterrows():
            # if the value is np.nan, make it black
            if not legal_color(row["color"]):
               rbh_to_color[row["rbh"]] = "#000000"
            else:
                rbh_to_color[row["rbh"]] = row["color"]

    # change the fontsize of the ticks on the arc plot
    arc_panel.tick_params(axis="y", labelsize=fontsize)
    # turn off the y-axis for the arc plot
    arc_panel.tick_params(axis="y", left = False, labelleft=False)
    # also turn off the left spine
    arc_panel.spines["left"].set_visible(False)

    # Now let's plot the stable pairs
    stable_pairs = unique_pairs[unique_pairs[pairtype] == 1]
    # For each stable pair, find if that stable pair is present in BCnSSimakov2022_gene,
    #  and plot an arc between the two points on the top marginal plot, the height is relative
    #  to the distance between the two points divided by the max distance of any pair
    stable_pair_entry = []
    # we do this now so we can simplify the plotting code
    for i, row in stable_pairs.iterrows():
        one_in = row["ortholog1"] in rbh["BCnSSimakov2022_gene"].values
        two_in = row["ortholog2"] in rbh["BCnSSimakov2022_gene"].values
        if one_in and two_in:
            x1 = rbh[rbh["BCnSSimakov2022_gene"] == row["ortholog1"]][f"{genome1}_plotpos"].values[0]
            x2 = rbh[rbh["BCnSSimakov2022_gene"] == row["ortholog2"]][f"{genome1}_plotpos"].values[0]
            x1_chrom = rbh[rbh["BCnSSimakov2022_gene"] == row["ortholog1"]][f"{genome1}_scaf"].values[0]
            x2_chrom = rbh[rbh["BCnSSimakov2022_gene"] == row["ortholog2"]][f"{genome1}_scaf"].values[0]
            if x1_chrom == x2_chrom:
                y1 = rbh[rbh["BCnSSimakov2022_gene"] == row["ortholog1"]]["BCnSSimakov2022_plotpos"].values[0]
                y2 = rbh[rbh["BCnSSimakov2022_gene"] == row["ortholog2"]]["BCnSSimakov2022_plotpos"].values[0]
                xD = abs(x2 - x1)
                yD = abs(y2 - y1)
                entry = {"ortholog1": row["ortholog1"],
                         "ortholog2": row["ortholog2"],
                         "x1": x1,
                         "x2": x2,
                         "x1_chrom": x1_chrom,
                         "x2_chrom": x2_chrom,
                         "y1": y1,
                         "y2": y2,
                         "xD": xD,
                         "yD": yD}
                stable_pair_entry.append(entry)
    plot_stable = pd.DataFrame(stable_pair_entry)
    # now plot arcs
    if len(plot_stable) == 0:
        max_xD = 1
    else:
        max_xD = plot_stable["xD"].max()
    rbhs_to_plot = set()
    for i, row in plot_stable.iterrows():
        modulo = 1
        color = "#000000"
        if blastdf is not None:
            if row["ortholog1"] in blast_rbh_list or row["ortholog2"] in blast_rbh_list:
                #print("These are two stable OGs that will be plotted below: ", row["ortholog1"], row["ortholog2"])
                #print(row)
                modulo = -1
                color = "#FF0000"
                if row["ortholog1"] in blast_rbh_list:
                    rbhs_to_plot.add(row["ortholog1"])
                if row["ortholog2"] in blast_rbh_list:
                    rbhs_to_plot.add(row["ortholog2"])
        # plot an arc from x1 to x2, passes through (in the middle) xD/max_xD
        if attempt_to_gradient_color:
            startcolor = rbh_to_color[row["ortholog1"]]
            stopcolor  = rbh_to_color[row["ortholog2"]]
            arc_panel = plot_bezier_arc(arc_panel, row["x1"], 0, row["x2"], 0, modulo * row["xD"]/max_xD, color, 0.6,
                                        color_gradient=(startcolor, stopcolor))
        else:
            arc_panel = plot_bezier_arc(arc_panel, row["x1"], 0, row["x2"], 0, modulo * row["xD"]/max_xD, color, 0.6)

    # plot the blast results
    if len(rbhs_to_plot) > 0:
        # figure out the offsets
        scaf_to_offset = {}
        for i, scaf in enumerate(scaf_order):
            scaflen = scaf_to_len[scaf]
            if i == 0:
                scaf_to_offset[scaf] = 0
                continue
            else:
                scaf_to_offset[scaf] = scaf_to_offset[scaf_order[i-1]] + scaf_to_len[scaf_order[i-1]]
        print(scaf_to_offset)
        if (len(chromdf) == 0) or (chromdf is None):
            raise ValueError("You must specify a chrom file to plot the blast results")
        genes_to_plot = set()
        for locus in rbhs_to_plot:
            for thisgene in blast_rbh_list[locus]:
                genes_to_plot.add(thisgene)
        # get the gene-to-sseqid mapping
        # zip together these two fields from blastdf and turn into a dict
        blast_qseqid_to_sseqid = dict(zip(blastdf["qseqid"], blastdf["sseqid"]))
        proteins_to_plot = set([blast_qseqid_to_sseqid[gene] for gene in genes_to_plot])
        protein_to_qseqid = {}
        for thisprot in proteins_to_plot:
            # get the qseqids that map to this protein from the blastdf
            protein_to_qseqid[thisprot] = set(blastdf[blastdf["sseqid"] == thisprot]["qseqid"])
        print(genes_to_plot)
        print(proteins_to_plot)
        print(protein_to_qseqid)
        # plot each protein in the panel
        protein_to_color = {x: "#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for x in protein_to_qseqid}
        yoffset = 0.05
        for thisprot in protein_to_qseqid:
            # make a patch
            scaf  = blastdf[blastdf["sseqid"] == thisprot]["scaf"].values[0]
            start = chromdf[chromdf["genome"] == thisprot]["start"].values[0]
            start = scaf_to_offset[scaf] + start
            stop  = chromdf[chromdf["genome"] == thisprot]["stop"].values[0]
            stop  = scaf_to_offset[scaf] + stop
            width = stop - start
            height = yoffset * 2
            bottom_left = (start, -yoffset)
            rect = Rectangle(bottom_left, width, height, color=protein_to_color[thisprot], alpha=1, linewidth = 0)
            arc_panel.add_patch(rect)
        # for each protein, add a legend entry
        for i, thisprot in enumerate(protein_to_qseqid):
            proteins = ",".join(list(protein_to_qseqid[thisprot]))
            # get the rightmost value of the panel x-axis
            maxx = arc_panel.get_xlim()[1]
            arc_panel.text(maxx, -0.9 - (i * yoffset), f"{thisprot}: {proteins}", fontsize=fontsize, color=protein_to_color[thisprot], ha='right', va='center')

    # change the y limit of ax_arcx_stable to 1
    Mb_window = window/1000000
    if (blastdf is not None) and (len(rbhs_to_plot) > 0):
        arc_panel.set_ylim(-1, 1)
        if pairtype == "stable_in_clade":
            arc_panel.text(0, -0.9, f"Stable in clade and within {Mb_window:.2f} MBp of blast results",   color = "#FF0000", fontsize=fontsize)
        elif pairtype == "unstable_in_clade":
            arc_panel.text(0, -0.9, f"Unstable in clade and within {Mb_window:.2f} MBp of blast results", color = "#FF0000", fontsize=fontsize)
        elif pairtype == "close_in_clade":
            arc_panel.text(0, -0.9, f"Close in clade and within {Mb_window:.2f} MBp of blast results",    color = "#FF0000", fontsize=fontsize)
        # we should also plot the blast results
    else:
        arc_panel.set_ylim(0, 1)
    # add text in the top left corner of ax_arcx_stable saying "stable pairs in clade"
    if pairtype == "stable_in_clade":
        arc_panel.text(0, 0.9, "Stable pairs in clade", fontsize=fontsize)
    elif pairtype == "unstable_in_clade":
        arc_panel.text(0, 0.9, "Unstable pairs in clade", fontsize=fontsize)
    elif pairtype == "close_in_clade":
        arc_panel.text(0, 0.9, "Pairs close in clade", fontsize=fontsize)
    # remove the bottom, right, and top spines
    arc_panel.spines["bottom"].set_visible(False)
    arc_panel.spines["right"].set_visible(False)
    arc_panel.spines["top"].set_visible(False)
    return arc_panel

def main():
    fontsize = 6
    args = parse_args()
    odp_plot.format_matplotlib()

    # get the exact scaffold_size
    scaf_to_len = {}
    for record in fasta.parse(args.genome_file):
        scaf_to_len[record.id] = len(record.seq)

    window = args.window
    # if the blast file is specified, read it in
    if args.blast_file is not None:
        blastdf = pd.read_csv(args.blast_file, sep="\t", header=None)
        colnames = ["qseqid", "sseqid", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send",
                    "evalue", "bitscore", "scaf", "rbh_list", "pos_list", "dist_list"]
        blastdf.columns = colnames
        # eval the next three columns. They are lists of strings
        blastdf["rbh_list"]   = blastdf["rbh_list"].apply( ast.literal_eval)
        blastdf["pos_list"]   = blastdf["pos_list"].apply( ast.literal_eval)
        blastdf["dist_list"]  = blastdf["dist_list"].apply(ast.literal_eval)
        # for each cell, get all the indices with values that are less than 1000000
        blastdf["index_list"] = blastdf["dist_list"].apply(lambda x: [i for i, val in enumerate(x) if val < window])
        blastdf["rbh_list"]   = blastdf.apply(lambda x: [x["rbh_list"][i] for i in x["index_list"]], axis=1)
        blastdf["pos_list"]   = blastdf.apply(lambda x: [x["pos_list"][i] for i in x["index_list"]], axis=1)
        blastdf["dist_list"]  = blastdf.apply(lambda x: [x["dist_list"][i] for i in x["index_list"]], axis=1)
        blast_rbh_list = {}
        for i, row in blastdf.iterrows():
            for rbh in row["rbh_list"]:
                if rbh not in blast_rbh_list:
                    blast_rbh_list[rbh] = set()
                blast_rbh_list[rbh].add(row["qseqid"])
        print(blast_rbh_list)
        print(blastdf)

    # if the chrom file is not none, read it in
    if args.chrom_file is not None:
        chromdf = pd.read_csv(args.chrom_file, sep="\t", header=None)
        chromdf.columns = ["genome", "scaf", "strand", "start", "stop"]
        print(chromdf)

    # figure out which unique pairs we want to plot
    unique_pairs = pd.read_csv( args.unique_pairs_path, sep="\t" )
    unique_pairs = unique_pairs[unique_pairs["taxid"] == args.clade_for_pairs]
    print(unique_pairs)
    if len(unique_pairs) == 0:
        print(f"Could not find any unique pairs for clade {args.clade_for_pairs}")
        sys.exit()

    rbh = rbh_tools.parse_rbh(args.rbh_file)

    # get the columnname that ends with "_gene" that is not BCnSSimakov2022
    # this will be the first genome
    genome1 = [col for col in rbh.columns if col.endswith("_gene") and "BCnSSimakov2022" not in col][0].replace("_gene", "")
    # make a two-panel plot using matplotlib and the gridspec
    NCBI = ete3.NCBITaxa()
    cladename = NCBI.get_taxid_translator([args.clade_for_pairs])[args.clade_for_pairs]

    fig = plt.figure(figsize=(12, 12))
    # Add a gridspec with two rows and one columns and a ratio of 1 to 4 between
    # the size of the marginal Axes and the main Axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(6, 6,  width_ratios=(16, 1, 2, 4, 4, 4), height_ratios=(4, 4, 4, 2, 1, 16),
                          left=0.2, right=0.9, bottom=0.2, top=0.9,
                          wspace=0.05, hspace=0.05)
    # set the title of the whole plot
    plottitle  = f"ODP between {genome1} and BCnSSimakov2022 ALGs"
    plottitle += f"\n Genome loci stable in {cladename} (TaxID: {args.clade_for_pairs})"
    plottitle += f"\nwindow size for pair/blast inclusion: {window} bp"
    fig.suptitle(plottitle, fontsize=fontsize)

    # force gs[1] to be square
    # Create the Axes.
    ax = fig.add_subplot(gs[5, 0])
    # change the x-axis to be f"{genome1} scaffolds}"
    ax.set_xlabel(f"{genome1} scaffolds", fontsize = fontsize)
    ax.set_ylabel("BCnSSimakov2022 scaffolds", fontsize = fontsize)
    ax_mixing_frac   = fig.add_subplot(gs[4, 0], sharex=ax)
    ax_mixing_plot   = fig.add_subplot(gs[3, 0], sharex=ax)
    ax_arcx_stable   = fig.add_subplot(gs[2, 0], sharex=ax)
    ax_arcx_unstable = fig.add_subplot(gs[1, 0], sharex=ax)
    ax_arcx_close    = fig.add_subplot(gs[0, 0], sharex=ax)
    # turn off the x-tic axis for the marginal plot
    ax_arcx_stable.tick_params(  axis="x", bottom = False, labelbottom=False)
    ax_arcx_unstable.tick_params(axis="x", bottom = False, labelbottom=False)
    ax_arcx_close.tick_params(   axis="x", bottom = False, labelbottom=False)

    # plot the rbh on the main plot
    ax.scatter(rbh[f"{genome1}_plotpos"], rbh["BCnSSimakov2022_plotpos"], s=1, c=rbh["color"])
    # change the axes to start at 0 and go to the max of the plotpos for both x and y
    ax.set_xlim(0, rbh[f"{genome1}_plotpos"].max())
    ax.set_ylim(0, rbh["BCnSSimakov2022_plotpos"].max())

    # sort rbh by plotpos for x, make vertical lines between the transitions for scaffolds
    #pd.set_option('display.max_rows', None)
    #pd.set_option('display.max_columns', None)
    #pd.set_option('display.width', None)

    # THIS MAKES THE X-AXIS TICKS
    rbh = rbh.sort_values(by=f"{genome1}_plotpos")
    rbh = rbh.reset_index(drop=True)
    transition_starts = []
    labels = []
    for i, row in rbh.iterrows():
        if i == 0:
            labels.append(row[f"{genome1}_scaf"])
            continue
        if row[f"{genome1}_scaf"] != rbh.iloc[i-1][f"{genome1}_scaf"]:
            transition_starts.append(i)
            labels.append(row[f"{genome1}_scaf"])
    midlines = []
    for i in transition_starts:
        midline = ((rbh.iloc[i][f"{genome1}_plotpos"] - rbh.iloc[i-1][f"{genome1}_plotpos"])/2) + rbh.iloc[i-1][f"{genome1}_plotpos"]
        midlines.append(midline)
        # plot the midline as a vertical line
        ax.axvline(midline, color="black", lw=0.25)
    # figure out the midpoints of the transitions
    tickpositions = []
    for i in range(len(midlines)):
        if i == 0:
            tickpositions.append(midlines[i]/2)
        else:
            tickpositions.append((midlines[i] - midlines[i-1])/2 + midlines[i-1])
    # for the last midpoint, (subtract the max of _plotpos from the last transition start)/2 + the last transition start
    tickpositions.append((rbh[f"{genome1}_plotpos"].max() - midlines[-1])/2 + midlines[-1])
    # plot the xticks and xlabels
    ax.set_xticks(tickpositions)
    ax.set_xticklabels(labels, rotation=90, fontsize = fontsize)
    xlabels = labels

    # Now we have labels for the x-axis, we can now plot the chromosome composition.
    gb = rbh.groupby(f"{genome1}_scaf")
    profiles = []
    panel_max = 0
    # binsize is the max plotpos/1000
    binsize = max(rbh[f"{genome1}_plotpos"])/1000
    for scaf, group in gb:
        group = group.sort_values(by=f"{genome1}_plotpos", ascending = True)
        # do a rolling count of 20 loci in rbh, centered on the current locus
        xstart = group[f"{genome1}_plotpos"].min()
        xstop  = group[f"{genome1}_plotpos"].max()
        # split this up into 200 bins, look in the left, middle, and right bins to get the color counts
        bins = [(x, x+binsize) for x in np.arange(xstart, xstop, binsize)]
        lookaround_width = 2
        for i in range(len(bins)):
            indices = [x for x in range(i-lookaround_width, i+lookaround_width+1) if x >= 0 and x < len(bins)]
            left = bins[indices[0]][0]
            right = bins[indices[-1]][1]
            plotleft  = bins[i][0]
            plotright = bins[i][1]
            plotmiddle = plotleft + ((plotleft + plotright)/2)
            clipped = group[(group[f"{genome1}_plotpos"] >= left) & (group[f"{genome1}_plotpos"] < right)]
            colorcounts = clipped["color"].value_counts().sort_index()
            countsum = colorcounts.sum()
            if int(countsum/2) +1 > panel_max:
                panel_max = int(countsum/2) + 1
            bottom = 0 - (countsum/2)
            bottom_frac = 0
            profiles.append(countsum/2)
            for thiscolor in colorcounts.index:
                count      = colorcounts[thiscolor]
                top = bottom + count
                rect = Rectangle((plotleft, bottom), binsize, count, color=thiscolor, alpha=1, linewidth = 0)
                ax_mixing_plot.add_patch(rect)
                bottom = top
                # now we handle fracs
                rect = Rectangle((plotleft, bottom_frac), binsize, count/countsum, color=thiscolor, alpha=1, linewidth = 0)
                ax_mixing_frac.add_patch(rect)
                bottom_frac += count/countsum

    ax_mixing_plot.set_ylim(0-panel_max, panel_max)
    # turn off the panel ticks, the spines
    for thisax in [ax_mixing_plot, ax_mixing_frac]:
        thisax.tick_params(axis="x", bottom = False, labelbottom=False)
        thisax.tick_params(axis="y", left = False, labelleft=False)
        thisax.spines["bottom"].set_visible(False)
        thisax.spines["left"].set_visible(False)
        thisax.spines["right"].set_visible(False)
        thisax.spines["top"].set_visible(False)
        for i in midlines:
            thisax.axvline(i, color="black", lw=0.25)

    # THIS MAKES THE Y-AXIS TICKS
    rbh = rbh.sort_values(by="BCnSSimakov2022_plotpos")
    rbh = rbh.reset_index(drop=True)
    transition_starts = []
    labels = []
    for i, row in rbh.iterrows():
        if i == 0:
            labels.append(row["BCnSSimakov2022_scaf"])
            continue
        if row["BCnSSimakov2022_scaf"] != rbh.iloc[i-1]["BCnSSimakov2022_scaf"]:
            transition_starts.append(i)
            labels.append(row["BCnSSimakov2022_scaf"])
    midlines = []
    for i in transition_starts:
        midline = ((rbh.iloc[i]["BCnSSimakov2022_plotpos"] - rbh.iloc[i-1]["BCnSSimakov2022_plotpos"])/2) + rbh.iloc[i-1]["BCnSSimakov2022_plotpos"]
        midlines.append(midline)
        # plot the midline as a vertical line
        ax.axhline(midline, color="black", lw=0.25)
    # figure out the midpoints of the transitions
    tickpositions = []
    for i in range(len(midlines)):
        if i == 0:
            tickpositions.append(midlines[i]/2)
        else:
            tickpositions.append((midlines[i] - midlines[i-1])/2 + midlines[i-1])
    # for the last midpoint, (subtract the max of _plotpos from the last transition start)/2 + the last transition start
    tickpositions.append((rbh["BCnSSimakov2022_plotpos"].max() - midlines[-1])/2 + midlines[-1])
    # plot the xticks and xlabels
    ax.set_yticks(tickpositions)
    ax.set_yticklabels(labels, fontsize = fontsize)

    # change the fontsize of the ticks on the arc plot
    if args.blast_file is not None:
        ax_arcx_stable = plot_arcs(ax_arcx_stable, genome1,     "stable_in_clade",
                                   rbh, unique_pairs, fontsize,
                                   blastdf = blastdf, blast_rbh_list = blast_rbh_list,
                                   chromdf = chromdf,
                                   scaf_order = xlabels, scaf_to_len = scaf_to_len,
                                   window = window, attempt_to_gradient_color=True)
        ax_arcx_unstable = plot_arcs(ax_arcx_unstable, genome1, "unstable_in_clade",
                                   rbh, unique_pairs, fontsize,
                                   blastdf = blastdf, blast_rbh_list = blast_rbh_list,
                                   chromdf = chromdf,
                                   scaf_order = xlabels, scaf_to_len = scaf_to_len,
                                   window = window, attempt_to_gradient_color=True)
        ax_arcx_close    = plot_arcs(ax_arcx_close, genome1,    "close_in_clade",
                                   rbh, unique_pairs, fontsize,
                                   blastdf = blastdf, blast_rbh_list = blast_rbh_list,
                                   chromdf = chromdf,
                                   scaf_order = xlabels, scaf_to_len = scaf_to_len,
                                   window = window, attempt_to_gradient_color=True)
    else:
        ax_arcx_stable   = plot_arcs(ax_arcx_stable, genome1,   "stable_in_clade",
                                   rbh, unique_pairs, fontsize,
                                   attempt_to_gradient_color = True)
        ax_arcx_unstable = plot_arcs(ax_arcx_unstable, genome1, "unstable_in_clade",
                                   rbh, unique_pairs, fontsize,
                                   attempt_to_gradient_color = True)
        ax_arcx_close    = plot_arcs(ax_arcx_close, genome1,    "close_in_clade",
                                   rbh, unique_pairs, fontsize,
                                   attempt_to_gradient_color = True)
    # save to a pdf
    outfile = f"{genome1}_pairsfrom_{cladename}_{args.clade_for_pairs}"
    if args.blast_file is not None:
        if len(blastdf) > 0:
            # get the first 3 blast hits, concatenate with _
            blaststring = "_blast_" + "_".join(blastdf["qseqid"].unique()[:3].tolist())
            outfile += blaststring
        else:
            "no_blast_results"
    outfile += ".pdf"
    plt.savefig(outfile)
    plt.close()

if __name__ == "__main__":
    main()