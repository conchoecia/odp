#!/usr/bin/env python
"""
This plot contains all of the functions
  related to making the ribbon plots.

These functions can be called by the other
  programs.
"""

# ODP-specific imports
import odp_plotting_functions as odp_plot

# Other standard python libraries
import logging

# non-standard dependencies
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.path as mpath
import matplotlib.patches as patches
import pandas as pd

import sys

def _quality_check_chromosome_list(sp, templist, sp_to_chr_to_size,
                                   sp_to_gene_order, sp_min_chr_size):
    """
    We perform some quality checks on the list of chromosomes to make sure that they are all valid.

    templist is the temporary list of scaffold IDs
    """
    # Now we perform some checks on the list of chromosomes and add them to the dict
    # check here to make sure that there aren't duplicate entries in the chromosomes
    if len(templist) != len(set(templist)):
        raise IOError("There are some duplicate chromosome entries for {}: {}.".format(sp, templist))
    # Make sure that all the chromosomes in the chromorder are actually in the fasta file
    checklist = [x for x in templist if x not in sp_to_chr_to_size[sp]]
    if len(checklist) > 0:
        raise IOError("These chromosomes in the chromosome order list for sp {} were not in the fasta file: {}".format(
            sp, checklist))
    # if we included a list of chromosomes below, remove those 
    #TODO DELETE
    print("sp_to_gene_order: {}".format(sp_to_gene_order))
    if sp in sp_to_gene_order:
        templist2 = [x for x in templist if x in sp_to_gene_order[sp]]
        templist = templist2
    else:
        # remove things smaller than the minimum size
        templist2 = [x for x in templist if (x in sp_to_chr_to_size[sp]) and (sp_to_chr_to_size[sp][x] >= sp_min_chr_size[sp])]
        templist = templist2
    return templist

def _optimize_spA_based_on_rbh(sp, prevsp, rbhdf, sp_to_chromorder, sp_to_gene_order):
    """
    Takes in a species, and optimizes the order of the chromosomes based on the rbh file and the order of the other species' chromosomes.

    sp is the species for which we want to optimize the order
    prevsp is the other species in the rbh dataframe

    Returns the optimized order of this species' chromosomes
    """
    # get the significance of the interactions
    rbhdf = rbhdf.copy()
    rbhdf = rbhdf[["{}_scaf".format(x) for x in [sp, prevsp]] + ["whole_FET"]]
    countsdf = rbhdf.groupby(by=["{}_scaf".format(sp), "{}_scaf".format(prevsp)]).count().reset_index()
    countsdf.columns = ["{}_scaf".format(sp), "{}_scaf".format(prevsp), "count"]
    countsdf = countsdf.sort_values(by=["count"], ascending=False).reset_index(drop=True)
    rbhdf = rbhdf.drop_duplicates().reset_index(drop=True)
    rbhdf = rbhdf.merge(countsdf, on=["{}_scaf".format(sp), "{}_scaf".format(prevsp)])
    # get rid of scaffolds that weren't in the custom list, if the custom list exists
    if sp in sp_to_gene_order:
        rbhdf = rbhdf.loc[rbhdf["{}_scaf".format(prevsp)].isin(sp_to_chromorder[prevsp])]
    # put the scaffold in the best place based on what is above it
    rbhdf["sort_index"] = rbhdf["{}_scaf".format(prevsp)].map(sp_to_chromorder[prevsp])
    rbhdf.dropna(subset=["sort_index"], inplace=True)
    rbhdf = rbhdf.sort_values(by=["{}_scaf".format(sp), "whole_FET"], ascending=[True, True])
    rbhdf.drop_duplicates(subset=["{}_scaf".format(sp)], inplace=True)
    rbhdf.sort_values(by=["sort_index", "whole_FET"], ascending=[True, True], inplace=True)
    rbhdf.reset_index(drop=True, inplace=True)
    return rbhdf["{}_scaf".format(sp)].tolist()

def plot_bezier_lines(panel, topxL, bottomxL, colors, alpha, topy, bottomy):
    """
    Plot bezier curves between chromosome coordinates of different species.
    - Returns the panel, but with the lines now

    Parameters:
      - panel:    The panel on which we want to plot the lines
      - topxL:    X-axis value where the lines start on the top of the plot
      - bottomxL: X-axis value where the lines stop on the bottom of the plot
      - colors:   The colors of each of the lines
      - alpha:    The alpha values of each of the lines
      - topy:     The y-value of the top of the lines
      - bottomy:  The y-value of the bottom of the lines
    """
    indent = 1.0
    # plot the indices
    for i in range(len(topxL)):
        topx     = topxL[i]
        bottomx  = bottomxL[i]
        diff = abs(topx - bottomx)
        middlex = min([topx, bottomx]) + (diff/2)
        leftx  = middlex - ((diff/2)*indent)
        rightx = middlex + ((diff/2)*indent)
        second = (-1,-1)
        third  = (-1,-1)
        if topx <= bottomx:
            second = (leftx,  topy+0.5)
            third  = (rightx, topy+0.5)
        else:
            second = (rightx, topy+0.5)
            third  = (leftx,  topy+0.5)
        path_data = [
            (Path.MOVETO, (topx, topy)),
            (Path.CURVE4, second),
            (Path.CURVE4, third),
            (Path.CURVE4, (bottomx, bottomy)),
            ]
        codes, verts = zip(*path_data)
        path  = mpath.Path(verts, codes)
        if colors[i] == "#000000":
            zord = -50
        elif alpha[i] < 0.5:
            zord = -99
        else:
            zord = 1
        patch = patches.PathPatch(path, fill = False,
                                  facecolor=[0,0,0,0], lw = 0.25,
                                  alpha=alpha[i], edgecolor = colors[i],
                                  zorder = zord)
        panel.add_patch(patch)
    return panel

def ribbon_plot(species_order, rbh_filelist,
                sp_to_chr_to_size,
                sp_min_chr_size, outfile,
                sp_to_gene_order = {},
                chr_sort_order   = "custom",
                plot_all = False):
    """
    Takes in a list of species as the plotting order,
     a list of rbh files, and a dict of species_to_chr_to_sizes

    In the future for the rbh file parse the headers.

    There are several ways that the chromosomes can be sorted.
      chr_sort_order < custom | optimal-top | optimal-size | optimal-random >
        custom         - use the custom sorting order for EVERY species in chromorder
        optimal-top    - use the custom order for the topmost species, then optimizes everything else
        optimal-size   - sort the top species' chromosomes by number of genes, optimize everything else
        optimal-chr-or - use `chromorder` when possible, optimize everything else
        optimal-random - randomly sort the chromosomes of the top species, optimize everything else
    """

    # check sp_to_gene_order. Reset to empty dict if None
    if type(sp_to_gene_order) == type(None):
        sp_to_gene_order = {}

    import random

    # make a list of the dataframes to open for these analyses
    rbh_df_list = [pd.read_csv(x, sep = "\t",
                   header = "infer",
                   index_col = None) for x in rbh_filelist]
    
    # This list tells the program later which index of RBH files to look at.
    # Useful for cases where we submit just one .rbh file with all of the entries
    sp_pair_to_rbh_df_list_index = {}
    for i in range(len(rbh_df_list)):
        thisdf = rbh_df_list[i]
        sp_list = [x.split("_")[0] for x in thisdf.columns if "_scaf" in x]
        for j in range(len(sp_list)-1):
            for k in range(j+1, len(sp_list)):
                sp_combo = tuple(sorted([sp_list[j], sp_list[k]]))
                if sp_combo not in sp_pair_to_rbh_df_list_index:
                    sp_pair_to_rbh_df_list_index[sp_combo] = i

    sp_to_genesdfs = {}
    # make composite gene index dataframe for each species
    #  we will use this later to plot by gene index
    #  rather than the chromosome index

    for thisrbh in rbh_df_list:
        thesesp = [x.split("_")[0] for x in thisrbh.columns if "_scaf" in x]
        for thissp in thesesp:
            if thissp not in sp_to_genesdfs:
                sp_to_genesdfs[thissp] = []
                sp_to_genesdfs[thissp].append(thisrbh[[x for x in thisrbh.columns if thissp in x]])

    # we now concat and deduplicate the gene index dfs
    for thissp in sp_to_genesdfs:
        sp_to_genesdfs[thissp] = pd.concat(sp_to_genesdfs[thissp]
                                  ).sort_values(by=["{}_gene".format(thissp)]
                                  ).drop_duplicates(subset=["{}_gene".format(thissp)]
                                  ).sort_values(by=["{}_scaf".format(thissp),
                                                    "{}_pos".format(thissp)],
                                                ascending=[True, True]
                                  ).reset_index(drop=True)
        sp_to_genesdfs[thissp] = sp_to_genesdfs[thissp][["{}_gene".format(thissp),
                                                         "{}_scaf".format(thissp),
                                                         "{}_pos".format(thissp)]]
        sp_to_genesdfs[thissp] = sp_to_genesdfs[thissp].groupby(
                                   by=["{}_scaf".format(thissp)], as_index=False
                                   ).apply(lambda x: x.reset_index(drop=True)
                                   )
        sp_to_genesdfs[thissp].index.names = ["{}_delete".format(thissp),"{}_chromIx".format(thissp)]
        sp_to_genesdfs[thissp].reset_index(inplace =True)
        sp_to_genesdfs[thissp] = sp_to_genesdfs[thissp][[x for x
                                  in sp_to_genesdfs[thissp] if "delete" not in x]]
        #print(sp_to_genesdfs[thissp])

    # This block is used for determining the chromosome order for the figure
    sp_to_chromorder     = {}
    print("species_order: {}".format(species_order))
    for i in range(0, len(species_order)):
        sp = species_order[i]

        templist = [] # templist is used to hold the chromosome order
        if i == 0:
            # This is the first species, there is a special case for it
            #      chr_sort_order < custom | optimal-top | optimal-size | optimal-random >
            if chr_sort_order in ["custom", "optimal-top", "optimal-chr-or"]:
                # in this case we use the custom order or take the order from the file
                if not sp_to_gene_order or (sp not in sp_to_gene_order):
                    templist = sorted(list(rbh_df_list[i]["{}_scaf".format(sp)].unique()))
                else:
                    templist = sp_to_gene_order[sp]
            elif chr_sort_order in ["optimal-size", "optimal-random"]:
                # make a sort of descending size of genes on the pairwise comparisons
                templist = list(rbh_df_list[i]["{}_scaf".format(sp)].value_counts().index)
                if sp in sp_to_gene_order:
                    # append the chromosomes that are in the config file but not in the rbh file
                    for x in sp_to_gene_order[sp]:
                        if x not in templist:
                            templist.append(x)
                if chr_sort_order == "optimal-random":
                    random.shuffle(templist)
        else:
            # we construct this variable to lookup the species in sp_pair_to_rbh_df_list_index
            sp_pair_lookup = tuple(sorted([species_order[i-1], sp])) 
            lookup_index = sp_pair_to_rbh_df_list_index[sp_pair_lookup]

            # this is the second or later species, change sort order depending on what was there
            if chr_sort_order == "custom":
                # in this case we use the custom order or take the order from the file
                if not sp_to_gene_order or (sp not in sp_to_gene_order):
                    # This is the case where we don't have a custom order for this species
                    #templist = sorted(list(rbh_df_list[i]["{}_scaf".format(sp)].unique()))
                    templist = sorted(list(rbh_df_list[lookup_index]["{}_scaf".format(sp)].unique()))
                else:
                    templist = sp_to_gene_order[sp]
            elif (chr_sort_order == "optimal-chr-or") and (sp in sp_to_gene_order):
                    templist = sp_to_gene_order[sp]
            else:
                # we optimize every other case
                prevsp = species_order[i-1]
                templist = _optimize_spA_based_on_rbh(sp, prevsp, rbh_df_list[lookup_index],
                                                      sp_to_chromorder, sp_to_gene_order)
        # Now we perform some checks on the list of chromosomes and add them to the dict
        # check here to make sure that there aren't duplicate entries in the chromosomes
        templist = _quality_check_chromosome_list(sp, templist, sp_to_chr_to_size,
                                                  sp_to_gene_order, sp_min_chr_size)
        sp_to_chromorder[sp] = {templist[i]: i for i in range(len(templist))}
        #print(sp, sp_to_chromorder)

    # There is an edge case where, if we used option "optimal-chr-or", but the 0th through nth species weren't
    #  in the chromosome order, we need to optimize the species in reverse order
    if chr_sort_order == "optimal-chr-or":
        # If someone uses this option when there is no chromosome order
        #  we just skip this step
        if len(sp_to_gene_order) > 0:
            first_optimized = 99999999
            for i in range(len(species_order)):
                sp = species_order[i]
                if sp in sp_to_gene_order:
                    first_optimized = i
                    break
            while first_optimized != 0:
                sp = species_order[first_optimized - 1]
                prevsp = species_order[first_optimized]
                # optimize the chromosome order
                templist = _optimize_spA_based_on_rbh(sp, prevsp, rbh_df_list[first_optimized - 1], sp_to_chromorder, sp_to_gene_order)
                templist = _quality_check_chromosome_list(sp, templist, sp_to_chr_to_size, sp_to_gene_order, sp_min_chr_size)
                sp_to_chromorder[sp] = {templist[i]: i for i in range(len(templist))}
                first_optimized -= 1

    # now construct dataframes describing how to plot the chromosomes based on
    #  gene index or on chromosome coordinate
    sp_to_chrom_to_order = {}
    sp_to_chromdf        = {} # this is for plotting by chromosome or gene ix coordinates
    for i in range(len(species_order) - 1):
        thissp = species_order[i]
        nextsp = species_order[i+1]
        # make a sp_to_chromorder dict entry
        for sp in [thissp, nextsp]:
            ########################################
            # START THE DATAFRAME FOR THIS SPECIES
            ########################################
            chromdf = pd.DataFrame([k for k,v in sorted(sp_to_chromorder[sp].items(),
                                              key=lambda item: item[1])],
                                   columns=['chrom'])

            ########################################
            # MAGIC NUMBERS
            ########################################
            ## We make a special dataframe to figure out how to plot the chromosomes in percent of chroms
            #space_between_chrom = 0.0175
            #percent_chrom_as_spaces = space_between_chrom * (len(sp_to_chromorder[sp]) - 1)
            # This doesn't work well for large numbers of chromosomes, so use a specific percent.
            # - I am finding that it helps to use different percents depending on how many chromosomes there are.
            #   This is an awkward block of code, but it is readable and is easy to change.
            num_chroms = len(sp_to_chromorder[sp])
            min_chr_size = 10
            max_chr_size = 50
            min_gap_percent = 0.12
            max_gap_percent = 0.5
            if num_chroms < min_chr_size:
                percent_chrom_as_spaces = min_gap_percent
            elif num_chroms > max_chr_size:
                percent_chrom_as_spaces = max_gap_percent
            else:
                # scale the percent_chrom_as_spaces between min_gap_percent and max_gap_percent based on the number of chromosomes
                percent_chrom_as_spaces = min_gap_percent + (max_gap_percent - min_gap_percent) * (num_chroms - min_chr_size) / (max_chr_size - min_chr_size)
            print("sp: {}, num_chroms: {}, percent_chrom_as_spaces: {}".format(sp, num_chroms, percent_chrom_as_spaces))
            percent_chrom_as_chroms = 1 - percent_chrom_as_spaces
            space_between_chrom = percent_chrom_as_spaces / (len(sp_to_chromorder[sp]) - 1)

            ######################################################
            # THIS SECTION IS FOR ABSOLUTE CHROMOSOME COORDINATES
            ######################################################
            # We make a special dataframe to figure out how to plot the chromosomes in absolute basepair coordinates.
            #  We determine what percent of the horizontal line will be occupied by gaps.
            #  This number space_between_chrom is a value between 0-1.
            chromdf["chromix"] = chromdf["chrom"].map(sp_to_chromorder[sp])
            chromdf["chrsize"] = chromdf["chrom"].map(sp_to_chr_to_size[sp])
            total_chrom_len = sum(chromdf["chrsize"])
            chromdf["chrPlotPercent"] = (chromdf["chrsize"]/total_chrom_len
                                         ) * percent_chrom_as_chroms
            chromdf["chrPlotOffset"] = chromdf["chrPlotPercent"].cumsum(
                                       ).shift(1).fillna(0) + \
                                       (space_between_chrom * chromdf["chromix"])

            ######################################################
            # THIS SECTION IS FOR RELATIVE CHROMOSOME COORDINATES
            ######################################################
            # This section does the same thing as above, but it plots the
            #  chromosomes as relative coordinates based on number of genes.
            chromdf["ixSize"] = chromdf["chrom"].map(
                sp_to_genesdfs[sp].groupby("{}_scaf".format(sp)).size().to_dict()).fillna(0)
            totalIxSize = sum(chromdf["ixSize"])
            chromdf["ixPlotPercent"] = (chromdf["ixSize"]/totalIxSize) * percent_chrom_as_chroms
            chromdf["ixPlotOffset"] = chromdf["ixPlotPercent"].cumsum(
                                       ).shift(1).fillna(0) + \
                                       (space_between_chrom * chromdf["chromix"])
            #print(chromdf)
            sp_to_chromdf[sp] = chromdf

    # CALL THIS TO GET THE VISUAL STYLE WE NEED
    odp_plot.format_matplotlib()

    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
    # Preserve the vertical order of embedded images:
    matplotlib.rcParams['image.composite_image'] = False
    # text as font in pdf
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    # first we need to figure out the dimensions of the figure.
    #  Just make it the number of samples times the space, plus the buffer
    interspeciesHeight = 0.5
    panelHeight = interspeciesHeight * len(species_order)
    panelWidth = 7.15

    #           two panels        top, bottom, middle
    bufferHeight = 1.5
    figHeight = (panelHeight*2) + (bufferHeight * 3)
    figWidth = 8

    fig = plt.figure(figsize=(figWidth,figHeight))

    #set the panel dimensions
    leftStart = (figWidth - panelWidth)/1.25
    bottomMargin = bufferHeight
    # pChr will host the chrom coordinate plots
    plt.gcf().text((leftStart + (figWidth/2))/figWidth,
                   (bottomMargin+panelHeight+0.25)/figHeight,
                   "p<=0.05 RBH results (Chr-coords)",
                   fontsize = 12, ha = "center", va = "bottom")
    pChr = plt.axes([leftStart/figWidth, #left
                   bottomMargin/figHeight,    #bottom
                   panelWidth/figWidth,   #width
                   panelHeight/figHeight])     #height
    pChr.tick_params(axis='both',which='both',
                   bottom=False, labelbottom=False,
                   left=False, labelleft=True,
                   right=False, labelright=False,
                   top=False, labeltop=False)
    # pChr will host the chrom coordinate plots
    plt.gcf().text((leftStart + (figWidth/2))/figWidth,
                   ((bottomMargin*2)+(panelHeight*2)+0.25)/figHeight,
                   "p<=0.05 RBH results (RBH-gene-coords)",
                   fontsize = 12, ha = "center", va = "bottom")
    pIx = plt.axes([leftStart/figWidth, #left
                   ((bottomMargin*2)+panelHeight)/figHeight,    #bottom
                   panelWidth/figWidth,   #width
                   panelHeight/figHeight])     #height
    pIx.tick_params(axis='both',which='both',
                   bottom=False, labelbottom=False,
                   left=False, labelleft=True,
                   right=False, labelright=False,
                   top=False, labeltop=False)
    # make a panel for the legend, too
    panellg = plt.axes([ 0/figWidth, #left
                         0/figHeight,    #bottom
                         figWidth/figWidth,   #width
                         (bottomMargin*0.75)/figHeight])     #height
    panellg.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=False, labelright=False,
                        top=False, labeltop=False)

    color_labels = set()
    # get the gene group column where the filtercol row is less than or equal to 0.05
    for thisdf in rbh_df_list:
        filtercol = "break_FET"
        if filtercol not in thisdf.columns:
            if "alpha" in thisdf.columns:
                filtercol = "alpha"
            elif "pval" in thisdf.columns:
                filtercol = "pval"
            else:
                # there is no acceptable column, raise an error
                raise ValueError("No acceptable column for filtering in colordf")
        # now that we have picked the filter column, get the color labels
        color_labels = color_labels.union( \
                        set(thisdf.loc[thisdf[filtercol] <= 0.05, "gene_group"].tolist())) 
    color_labels = sorted(list(color_labels))

    # Get a dataframe of group to color based on frequency
    colordf = pd.concat(rbh_df_list)[["gene_group", "color"]
               ].groupby(["gene_group", "color"]
               ).size(
               ).reset_index(name='Freq'
               ).sort_values(by = ["Freq"], ascending = False
               ).drop_duplicates(subset = "gene_group"
               ).sort_values(by = ["gene_group"], ascending = True
               ).reset_index(drop=True)

    # get the colors of the labels in color_labels
    # zip color_df["gene_group"] and color_df["color"] to a dict
    color_dict = dict(zip(colordf["gene_group"], colordf["color"]))
    color_colors = [color_dict[x] for x in color_labels]


    # make a legend
    # set a legend if it is prot_to_color_mode
    legend_elements = []
    for ii in range(len(color_labels)):
        legend_elements.append(
            patches.Patch(facecolor=color_colors[ii],
            edgecolor='black', lw = 0,
            label=color_labels[ii])
           )
    panellg.legend(handles=legend_elements,
                   ncol = 10,
                   fontsize = 8, loc='center left')

    # plot all the lines
    for i in range(len(species_order)-1):
        thissp = species_order[i]
        nextsp = species_order[i+1]
        thisspChrom = sp_to_chromdf[thissp]
        nextspChrom = sp_to_chromdf[nextsp]
        thisspGenes = sp_to_genesdfs[thissp]
        nextspGenes = sp_to_genesdfs[nextsp]

        sp_pair_lookup = tuple(sorted([thissp, nextsp])) 
        lookup_index = sp_pair_to_rbh_df_list_index[sp_pair_lookup]

        tr = rbh_df_list[lookup_index].copy()
        if not plot_all:
            # if we have performed FET use that, otherwise look up alpha
            filtercol = "break_FET"
            if "break_FET" not in tr.columns:
                if "alpha" in tr.columns:
                    filtercol = "alpha"
                else:
                    print("ERROR: no break_FET or alpha column in rbh_df_list[{}]".format(i))
                    sys.exit(1)
            
            tr = tr.loc[tr[filtercol] <= 0.05, ]
            tr["alpha"] = tr.apply(lambda x: 0.8
                                   if x[filtercol] <= 0.05 else 0.1,
                                   axis = 1)
        tr = tr[["{}_gene".format(thissp),
                 "{}_scaf".format(thissp),
                  "{}_pos".format(thissp),
                 "{}_gene".format(nextsp),
                 "{}_scaf".format(nextsp),
                  "{}_pos".format(nextsp),
                  "alpha",
                 "color"]]
        # handle the index plotting info
        #print(thisspChrom)
        #print()

        # get the index plot position for the top
        chrom_to_ixSize    = dict(zip(thisspChrom["chrom"], thisspChrom["ixSize"]))
        chrom_to_ixPercent = dict(zip(thisspChrom["chrom"], thisspChrom["ixPlotPercent"]))
        chrom_to_ixOffset  = dict(zip(thisspChrom["chrom"], thisspChrom["ixPlotOffset"]))
        gene_to_ix = dict(zip(thisspGenes["{}_gene".format(thissp)],
                              thisspGenes["{}_chromIx".format(thissp)]))
        tr["topIx"] = tr["{}_gene".format(thissp)].map(gene_to_ix)
        tr["topIx_ChromSize"] = tr["{}_scaf".format(thissp)].map(chrom_to_ixSize)
        tr["topIx_ChromPercent"] = tr["{}_scaf".format(thissp)].map(chrom_to_ixPercent)
        tr["topIx_geneOffset"]   = (tr["topIx"]/tr["topIx_ChromSize"]) * tr["topIx_ChromPercent"]
        tr["topIx_ChromOffset"] = tr["{}_scaf".format(thissp)].map(chrom_to_ixOffset)
        tr["topIx_finalOffset"] = tr["topIx_ChromOffset"] + tr["topIx_geneOffset"]
        delete = ["topIx", "topIx_ChromSize", "topIx_ChromPercent",
                       "topIx_geneOffset","topIx_ChromOffset"]
        tr = tr[[x for x in tr.columns if x not in delete]]
        tr = tr.sort_values(by=["topIx_finalOffset"]).reset_index(drop=True)

        # get the index plot position for the bottom
        chrom_to_ixSize = dict(zip(nextspChrom["chrom"], nextspChrom["ixSize"]))
        chrom_to_ixPercent = dict(zip(nextspChrom["chrom"], nextspChrom["ixPlotPercent"]))
        chrom_to_ixOffset = dict(zip(nextspChrom["chrom"], nextspChrom["ixPlotOffset"]))
        gene_to_ix = dict(zip(nextspGenes["{}_gene".format(nextsp)],
                              nextspGenes["{}_chromIx".format(nextsp)]))
        tr["bottomIx"] = tr["{}_gene".format(nextsp)].map(gene_to_ix)
        tr["bottomIx_ChromSize"] = tr["{}_scaf".format(nextsp)].map(chrom_to_ixSize)
        tr["bottomIx_ChromPercent"] = tr["{}_scaf".format(nextsp)].map(chrom_to_ixPercent)
        tr["bottomIx_geneOffset"]   = (tr["bottomIx"]/tr["bottomIx_ChromSize"]) * tr["bottomIx_ChromPercent"]
        tr["bottomIx_ChromOffset"] = tr["{}_scaf".format(nextsp)].map(chrom_to_ixOffset)
        tr["bottomIx_finalOffset"] = tr["bottomIx_ChromOffset"] + tr["bottomIx_geneOffset"]
        delete = ["bottomIx", "bottomIx_ChromSize", "bottomIx_ChromPercent",
                       "bottomIx_geneOffset","bottomIx_ChromOffset"]
        tr = tr[[x for x in tr.columns if x not in delete]]
        tr = tr.sort_values(by=["bottomIx_finalOffset"]).reset_index(drop=True)

        # get the chrom plot position for the top
        chrom_to_chrSize =    dict(zip(thisspChrom["chrom"], thisspChrom["chrsize"]))
        chrom_to_chrPercent = dict(zip(thisspChrom["chrom"], thisspChrom["chrPlotPercent"]))
        chrom_to_chrOffset =  dict(zip(thisspChrom["chrom"], thisspChrom["chrPlotOffset"]))
        tr["topChr_ChromSize"] = tr["{}_scaf".format(thissp)].map(chrom_to_chrSize)
        tr["topChr_ChromPercent"] = tr["{}_scaf".format(thissp)].map(chrom_to_chrPercent)
        tr["topChr_geneOffset"]   = (tr["{}_pos".format(thissp)
                                        ]/tr["topChr_ChromSize"]) * tr["topChr_ChromPercent"]
        tr["topChr_ChromOffset"] = tr["{}_scaf".format(thissp)].map(chrom_to_chrOffset)
        tr["topChr_finalOffset"] = tr["topChr_ChromOffset"] + tr["topChr_geneOffset"]
        delete = ["topChr", "topChr_ChromSize", "topChr_ChromPercent",
                       "topChr_geneOffset","topChr_ChromOffset"]
        tr = tr[[x for x in tr.columns if x not in delete]]
        tr = tr.sort_values(by=["topChr_finalOffset"]).reset_index(drop=True)

        # get the chrom plot position for the bottom
        chrom_to_chrSize =    dict(zip(nextspChrom["chrom"], nextspChrom["chrsize"]))
        chrom_to_chrPercent = dict(zip(nextspChrom["chrom"], nextspChrom["chrPlotPercent"]))
        chrom_to_chrOffset =  dict(zip(nextspChrom["chrom"], nextspChrom["chrPlotOffset"]))
        tr["bottomChr_ChromSize"] = tr["{}_scaf".format(nextsp)].map(chrom_to_chrSize)
        tr["bottomChr_ChromPercent"] = tr["{}_scaf".format(nextsp)].map(chrom_to_chrPercent)
        tr["bottomChr_geneOffset"]   = (tr["{}_pos".format(nextsp)
                                        ]/tr["bottomChr_ChromSize"]) * tr["bottomChr_ChromPercent"]
        tr["bottomChr_ChromOffset"] = tr["{}_scaf".format(nextsp)].map(chrom_to_chrOffset)
        tr["bottomChr_finalOffset"] = tr["bottomChr_ChromOffset"] + tr["bottomChr_geneOffset"]
        delete = ["bottomChr", "bottomChr_ChromSize", "bottomChr_ChromPercent",
                       "bottomChr_geneOffset","bottomChr_ChromOffset"]
        tr = tr[[x for x in tr.columns if x not in delete]]
        tr = tr.sort_values(by=["bottomChr_finalOffset"]).reset_index(drop=True)

        # plot the indices
        pIx = plot_bezier_lines(pIx,
                         tr["topIx_finalOffset"],
                         tr["bottomIx_finalOffset"],
                         tr["color"],
                         tr["alpha"],
                         i, i+1)

        pChr = plot_bezier_lines(pChr,
                         tr["topChr_finalOffset"],
                         tr["bottomChr_finalOffset"],
                         tr["color"],
                         tr["alpha"],
                         i, i+1)

    # Now we plot all the chroms
    for i in range(len(species_order)):
        sp = species_order[i]
        for index, row in sp_to_chromdf[sp].iterrows():
            # plot the line and the text
            x1 = row["chrPlotOffset"]
            x2 = row["chrPlotOffset"] + row["chrPlotPercent"]
            pChr.plot([x1,x2],[i,i],'k-')
            pChr.text(x1, i-0.03, row["chrom"],  ha="left", va = "bottom", fontsize =5)

            # plot the line and text for the index plot
            x1 = row["ixPlotOffset"]
            x2 = row["ixPlotOffset"] + row["ixPlotPercent"]
            pIx.plot([x1,x2],[i,i],'k-')
            pIx.text(x1, i-0.03, row["chrom"], ha="left", va = "bottom", fontsize = 4)

    # flip the y axes
    for p in [pChr, pIx]:
        # remove some spines
        for side in ["top", "right", "bottom", "left"]:
            p.spines[side].set_visible(False)
        # set the limits
        p.set_xlim([-0.02,1.02])
        p.set_ylim([-0.03,len(species_order)-1+0.03])
        # flip the axes
        p.set_ylim(p.get_ylim()[::-1])
        # set the axis labels
        p.set_yticks(list(range(len(species_order))))
        p.set_yticklabels(species_order)

    plt.savefig(outfile)
    return sp_to_chromorder
