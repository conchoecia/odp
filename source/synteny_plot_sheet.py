"""
This was taken out of the odp script.
"""

import os
from matplotlib              import font_manager # https://github.com/conchoecia/odp/issues/34#issuecomment-1636475654
from matplotlib.patches      import Patch
from matplotlib.lines        import Line2D
from matplotlib.font_manager import FontProperties
import sys

# ODP-specific requirements
import rbh_tools
from fishers_exact import FET_graph, FET_connected_components, FET_sample_to_chrom_to_CC
import odp_plotting_functions as odp_plot
import odp_color_manager as oc

def synteny_plot_sheet(df_file, plotdf_file, synplot,
                 xsample, ysample,
                 s2c2s,
                 xorder,  yorder,
                 xbreaks, ybreaks,
                 orthology_method,
                 sort_y_by_x = True,
                 prot_to_color_name = None,
                 prot_to_color_dir  = None,
                 prot_to_color_hmm_results = None):
    """
    If the user provided a plot order, then we should not skip any scaffolds.

    This is the main plotting script for the synteny plot.
    Parameters:
      - df_file - the file with the rbh data
      - plotdf_file - the file to which the final plotting data will be saved
      - synplot - the name of the synteny plot that will be saved as a pdf
      - xsample - the name of the x sample
      - ysample - the name of the y sample
      - s2c2s - species-to-chrom-to-size dictionary, map of the species, which chrom, and how big that chrom is
      - xorder - the order in which the x scaffolds should be plotted, if the user has a preference.
      - yorder - the order in which the y scaffolds should be plotted, if the user has a preference.
      - xbreaks - the x breaks for the plot. Could be centromeres, breaks in synteny, et cetera.
      - ybreaks - the y breaks for the plot. See above comment.
      - orthology_method - the orthology method used to generate the rbh file. Diamond or blastp.
      - sort_y_by_x - if True, then the y scaffolds will be sorted by the best match with x scaffolds.
      - prot_to_color_name - TODO finish this comment.
      - prot_to_color_dir - TODO finish this comment.
      - prot_to_color_hmm_results - TODO finish this comment.

    TODO: this function needs a parameter where the user passes a chrom-to-size file, or a dict of chroms to sizes.
          This would replace species_to_chrom_to_size() and would be much faster.
    """
    import pandas as pd
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    #import matplotlib.patches as mplpatches
    from matplotlib.ticker import StrMethodFormatter, NullFormatter
    import numpy as np

    prot_to_color_mode = False

    # check if we are in prot_to_color_mode
    if prot_to_color_name and os.path.exists(prot_to_color_dir):
        prot_to_color_mode = True

    # read in the rbh file. This will already have the Dx, Dy, and FET
    df = rbh_tools.parse_rbh(df_file)
    # Change the _ix columns to type float because later we take the mean of them.
    # We first observed this error in 2025 after recent updates to pandas.
    for thisspecies in [xsample, ysample]:
        df["{}_ix".format(thisspecies)] = df["{}_ix".format(thisspecies)].astype(float)

    # now make a df to figure out the sort order
    sort_df = df.copy()
    # We keep the scaf, the index to figure out the center of mass on the other chrom,
    #  and the FET value to figure out the order
    keep_these_columns = ["{}_scaf".format(xsample), "{}_scaf".format(ysample),
                          "{}_ix".format(xsample), "{}_ix".format(ysample),
                          "whole_FET"]
    sort_df = sort_df[keep_these_columns]

    # get the size of each group so we can add it to sort_df to sort by x-size
    size_df = sort_df.groupby(["{}_scaf".format(xsample), "{}_scaf".format(ysample)]).size().to_frame('size').reset_index()
    size_df = size_df.sort_values(by = ["size"], ascending = [False]).reset_index(drop=True)

    # get the size of the x and y chroms to add to the sort_df
    x_chrom_to_size = sort_df.groupby("{}_scaf".format(xsample)).size().to_frame('x_size').reset_index()
    y_chrom_to_size = sort_df.groupby("{}_scaf".format(ysample)).size().to_frame('y_size').reset_index()

    # groupby the scafs, figure out what the mean value of the other species' index is,
    #  then flatten back into a regular df
    # For example, if LVA is one species and BFL is another, the column "LVA_ix" will now
    #  refer to the mean index of the BFL genes that are in the same chrom pair as the LVA gene
    sort_df = sort_df.groupby(by=["{}_scaf".format(xsample), "{}_scaf".format(ysample)], dropna=False).mean().reset_index()

    # Now merge the size_df and sort_df. This only adds the column "size" to sort_df.
    # The column "size" is the number of genes in the x-chrom and y-chrom pair
    sort_df = sort_df.merge(size_df, on=["{}_scaf".format(xsample), "{}_scaf".format(ysample)], how="left").sort_values(
        by = ["size", "whole_FET"], ascending = [False, True]).reset_index(drop=True)

    # We still need to add the total number of genes on each chrom
    ysample_chrom_to_size = df.groupby("{}_scaf".format(ysample)).size().to_frame("{}_size".format(xsample)).reset_index()
    xsample_chrom_to_size = df.groupby("{}_scaf".format(xsample)).size().to_frame("{}_size".format(ysample)).reset_index()

    # merge the above two dataframes with sort_df
    sort_df = sort_df.merge(ysample_chrom_to_size, on=["{}_scaf".format(ysample)], how="left")
    sort_df = sort_df.merge(xsample_chrom_to_size, on=["{}_scaf".format(xsample)], how="left")
    # sort based on size and FET
    sort_df = sort_df.sort_values(by = ["{}_size".format(xsample), "whole_FET"], ascending = [False, True]).reset_index(drop=True)

    # now we must find the connected components to make ordering by FET clusters easier
    # This code was adapted from this SO question: https://stackoverflow.com/questions/53573865
    #  and this specific answer: https://stackoverflow.com/a/53574094/5843327
    #  by SO user: https://stackoverflow.com/users/4001592/dani-mesejo
    # first we prepend the scaf names for a fake column to get around limitations of the above code
    splitter = ":;:;:;:"
    sort_df["{}_scaf_absolute".format(xsample)] = "{}{}".format(xsample, splitter) + sort_df["{}_scaf".format(xsample)].astype(str)
    sort_df["{}_scaf_absolute".format(ysample)] = "{}{}".format(ysample, splitter) + sort_df["{}_scaf".format(ysample)].astype(str)
    G = FET_graph(sort_df.loc[sort_df["whole_FET"] <= 0.05,][["{}_scaf_absolute".format(xsample), "{}_scaf_absolute".format(ysample)]].values)
    components = [set(x) for x in FET_connected_components(G)]
    # convert components into a dictionary for the x or y-chroms lookup. We just want the CC, because we can filter on these later
    sample_to_chrom_to_CC = FET_sample_to_chrom_to_CC(components, xsample, ysample, splitter)
    # we have the connected components now, and can just add the group field with map on the xsample scafs
    sort_df['CC_group'] = sort_df["{}_scaf".format(xsample)].map(sample_to_chrom_to_CC[xsample])

    # Some of the x_samples will not have connected components, because they are not associated with any y_sample
    #  this will likely only happen with fragmented genomes, so we must add these as their own connected components
    xsamples_not_in_CC = set(sort_df["{}_scaf".format(xsample)].unique()) - set(sample_to_chrom_to_CC[xsample].keys())
    #for x in xsamples_not_in_CC:
    #    components.append(set(["{}{}{}".format(xsample, splitter, x)]))
    # These y-samples must be associated with some x-sample in order to determine how to order them
    ysamples_not_in_CC = set(sort_df["{}_scaf".format(ysample)].unique()) - set(sample_to_chrom_to_CC[ysample].keys())
    unplaced_ys = sort_df.loc[sort_df["{}_scaf".format(ysample)].isin(ysamples_not_in_CC),]
    unplaced_ys = unplaced_ys.sort_values(by = ["{}_scaf".format(ysample), "whole_FET"],
                                          ascending = [False, True]).reset_index(drop=True)
    # We just sorted these samples by y-sample scaffold, then FET. We now remove unique values
    #  to get a list of associations to add to the existing connected components
    unplaced_ys.drop_duplicates(subset = ["{}_scaf".format(ysample)], keep = "first", inplace = True)
    # groupby x-axis scaffold and do processing to add these to the connected components
    unplaced_ys_gb = unplaced_ys.groupby("{}_scaf_absolute".format(xsample))
    # - We go through each group of y scaffs associated with x-axis scafs, and add to the connected components
    # - If the connected component doesn't exist yet, we make it. 
    for name, grouped in unplaced_ys_gb:
        CC_index = grouped["CC_group"].unique()[0]
        if CC_index in list(range(len(components))):
            # convert to int now - we didn't case yet because it may have been None
            CC_index = int(CC_index)
            # this x-scaffold is already in a connected component, but add it anyway just to be safe :)
            #  and then we add the x-scaffold to the y-scaffold connected component
            newset = set.union(set([name]),
                               set(grouped["{}_scaf_absolute".format(ysample)]) )
            components[CC_index] = set.union(components[CC_index],
                                             newset)
        else:
            # this x-scaffold is not in a connected component, so we need to make a new one
            #  and then we add the x-scaffold to the y-scaffold connected component
            newset = set.union(set([name]),
                               set(grouped["{}_scaf_absolute".format(ysample)]) )
            components.append(newset)

    # recalculate this since we added some new connected components
    sample_to_chrom_to_CC = FET_sample_to_chrom_to_CC(components, xsample, ysample, splitter)
    # There might still be some unplaced scaffolds, so just add them in any order now
    xsamples_not_in_CC = set(sort_df["{}_scaf".format(xsample)].unique()) - set(sample_to_chrom_to_CC[xsample].keys())
    for x in xsamples_not_in_CC:
        components.append(set(["{}{}{}".format(xsample, splitter, x)]))

    # recalculate the connected components
    sample_to_chrom_to_CC = FET_sample_to_chrom_to_CC(components, xsample, ysample, splitter)
    # At this point all of the y-axis samples should be added to a connected component. Check to make sure.
    if set(sample_to_chrom_to_CC[ysample].keys()) != set(sort_df["{}_scaf".format(ysample)].unique()):
        print("ERROR: not all {} scaffolds have been identified in a connected component".format(ysample))
        missing_from_CC = set(sort_df["{}_scaf".format(ysample)].unique()) - set(sample_to_chrom_to_CC[ysample].keys())
        print("The missing scaffolds are:")
        print(missing_from_CC)
        sys.exit()
    else:
        print("All the {} scaffolds have been identified in a connected component, continuing".format(ysample), file = sys.stderr)
    # also check that the x-axis samples are all in a connected component
    if set(sample_to_chrom_to_CC[xsample].keys()) != set(sort_df["{}_scaf".format(xsample)].unique()):
        print("ERROR: not all {} scaffolds have been identified in a connected component".format(xsample))
        missing_from_CC = set(sort_df["{}_scaf".format(xsample)].unique()) - set(sample_to_chrom_to_CC[xsample].keys())
        print("The missing scaffolds are:")
        print(missing_from_CC)
        sys.exit()
    else:
        print("All the {} scaffolds have been identified in a connected component, continuing".format(xsample), file = sys.stderr)

    # remove the absolute fake columns - we don't need them now
    sort_df = sort_df[[x for x in sort_df.columns if "_absolute" not in x]]
    # go through again and add the CC_group to the dataframe based on the ysample scaf
    sort_df["CC_group"] = sort_df["{}_scaf".format(ysample)].apply(lambda x: sample_to_chrom_to_CC[ysample][x])
    # convert the CC_group to int
    sort_df["CC_group"] = sort_df["CC_group"].astype(int)
    # check that there are no Nones in the CC_group column
    if sort_df["CC_group"].isnull().values.any():
        print("ERROR: there are Nones in the CC_group column")
        sys.exit()

    # now we figure out the plotting order for both axes.
    # The algorithm is:
    #   - Find the largest x-scaffold that we have not yet plotted
    #   - Find the order of the y-scaffolds so that they are organized by the mean location along the x-axis
    #   - Keep track of which CCs have been added to the plotting order lists
    #   - Do this until all x-axis scaffolds have been seen. Check list of seen CCs against len of CCs to see if we missed something
    CCs_seen = set()
    x_scaffolds_by_size = df["{}_scaf".format(xsample)].value_counts().index.tolist() 
    if len(x_scaffolds_by_size) != len(df["{}_scaf".format(xsample)].unique()):
        print("ERROR: these x-axis counts should match")
        sys.exit()

    optimized_xorder = []
    optimized_yorder = []
    for thisx in x_scaffolds_by_size:
        if thisx in optimized_xorder:
            # we don't need to do anything
            pass
        else:
            # we need to add this x-scaffold to the optimized_xorder
            this_CC_index = sample_to_chrom_to_CC[xsample][thisx]
            CCs_seen.add(this_CC_index)
            x_keeps = [x for x in sample_to_chrom_to_CC[xsample] if sample_to_chrom_to_CC[xsample][x] == this_CC_index]
            x_scaffolds = [x for x in x_scaffolds_by_size if x in x_keeps]
            # only get the y-scaffolds that are in this CC
            tempdf = sort_df.loc[sort_df["CC_group"] == this_CC_index, ]
            for tempx in x_scaffolds:
                # and only get the x-scaffolds that are in this CC
                subtempdf = tempdf.loc[tempdf["{}_scaf".format(xsample)] == tempx, ]
                # sort by x-chrom scaffold and the mean position of all the genes on that scaffold
                subtempdf = subtempdf.sort_values(by = ["{}_scaf".format(xsample), "{}_ix".format(xsample)], ascending = [True, True])
                # now we add all of these things to the optimized_orders
                optimized_xorder.append(tempx)
                for tempy in subtempdf["{}_scaf".format(ysample)].unique().tolist():
                    if tempy not in optimized_yorder:
                        optimized_yorder.append(tempy)
    #print("optimized_xorder")
    #print(len(optimized_xorder), len(set(optimized_xorder)))
    #print("optimized_yorder")
    #print(len(optimized_yorder), len(set(optimized_yorder)))
    #print(optimized_yorder)

    # now we check that we have seen all the CCs
    if len(CCs_seen) != len(components):
        print("ERROR: we have not seen all the Connected Components")
        print("We have seen {} CCs, but there are {} CCs".format(len(CCs_seen), len(components)))
        print("The CCs we have seen are:")
        print(CCs_seen)
        print("The CCs we have not seen are:")
        missing_CCs = [x for x in range(len(components)) if x not in CCs_seen]
        print(missing_CCs) 
        sys.exit()
    # we check that we have seen all the x-scaffolds
    if set(optimized_xorder) != set(df["{}_scaf".format(xsample)].unique()):
        print("ERROR: we have not seen all the x-scaffolds")
        print("We have seen {} x-scaffolds, but there are {} x-scaffolds".format(len(set(optimized_xorder)), len(df["{}_scaf".format(xsample)].unique())))
        print("The x-scaffolds we have not seen are:")
        missing_x = [x for x in df["{}_scaf".format(xsample)].unique() if x not in set(optimized_xorder)]
        print(missing_x)
        sys.exit()
    # we also check that we have seen all the y-scaffolds
    if set(optimized_yorder) != set(df["{}_scaf".format(ysample)].unique()):
        print("ERROR: we have not seen all the y-scaffolds")
        print("We have seen {} y-scaffolds, but there are {} y-scaffolds".format(len(set(optimized_yorder)), len(df["{}_scaf".format(ysample)].unique())))
        print("The y-scaffolds we have not seen are:")
        missing_y = [x for x in df["{}_scaf".format(ysample)].unique() if x not in set(optimized_yorder)]
        print(missing_y)
        sys.exit()

    # at this point we now have plotting order for both the x and y scaffolds
    xorder = optimized_xorder
    yorder = optimized_yorder

    # now make a separate df for x and y for plotting order
    xdf = False
    xsorter = xorder
    if len(xorder) == 0:
        # sort by largest to smallest
        xgb = df.groupby(by="{}_scaf".format(xsample), dropna=False)
        xsorter = list(xgb.size().sort_values(ascending = False).index)
    xgb = df.groupby(by="{}_scaf".format(xsample))
    # sorts the rows by the scaf order, and within the scafs by the gene index
    xdf = pd.concat([xgb.get_group(x).sort_values(
        by = "{}_ix".format(xsample), ascending=True).reset_index(drop=True)
                           for x in xsorter]).reset_index(drop=True)
    xdf["{}_plotindex".format(xsample)] = xdf.index

    # I think it would be better to put this block of code later.
    # Specifically, I think it is important to wait until the last minute
    #  to set the plotpos for the x-axis chroms since things might need to have
    #  a more complicated sorting if there are splits or fusions
    xsorter_to_cumsum = [0]
    for i in range(1,len(xsorter)):
        xsorter_to_cumsum.append( xsorter_to_cumsum[i-1] + s2c2s[xsample][xsorter[i-1]])
    x_to_cumsum = {xsorter[i]:xsorter_to_cumsum[i] for i in range(len(xsorter_to_cumsum))}
    xdf["{}_plotpos".format(xsample)] = xdf["{}_scaf".format(xsample)].map(x_to_cumsum) + xdf["{}_pos".format(xsample)]

    # now make a separate df for x and y for plotting order
    ydf = False
    ysorter = yorder
    if len(yorder) == 0:
        # sort by largest to smallest
        ygb = df.groupby(by="{}_scaf".format(ysample), dropna=False)
        ysorter = list(ygb.size().sort_values(ascending = False).index)
    ygb = xdf.groupby(by="{}_scaf".format(ysample))
    ydf = pd.concat([ygb.get_group(x).sort_values(
        by = ["{}_ix".format(ysample)], ascending=True).reset_index(drop=True)
                           for x in ysorter]).reset_index(drop=True)
    ydf["{}_plotindex".format(ysample)] = ydf.index
    ysorter_to_cumsum = [0]
    for i in range(1,len(ysorter)):
        ysorter_to_cumsum.append( ysorter_to_cumsum[i-1] + s2c2s[ysample][ysorter[i-1]])
    y_to_cumsum = {ysorter[i]:ysorter_to_cumsum[i] for i in range(len(ysorter_to_cumsum))}
    ydf["{}_plotpos".format(ysample)] = ydf["{}_scaf".format(ysample)].map(y_to_cumsum) + ydf["{}_pos".format(ysample)]

    # this is what we plot
    plotdf = ydf
    #print("plotdf")
    #print(plotdf)

    # CALL THIS TO GET THE VISUAL STYLE WE NEED
    odp_plot.format_matplotlib()

    # now make a scatter plot
    #figDouble = 16
    figWidth  = 14
    figHeight = 18
    fig = plt.figure(figsize=(figWidth,figHeight))

    #set the panel dimensions
    panelWidth = 4
    panelHeight = 4
    dpanel_width = 0.25
    #find the margins to center the panel in figure
    leftStart = 2.5
    secondLeftStart = 8.5
    bottomMargin = 12.5
    # panel1 will host the index-based plot
    plt.gcf().text((leftStart-0.75)/figWidth,
                    (bottomMargin+panelHeight+0.25)/figHeight,
                   "a", weight = "bold",
                   fontsize = 20, ha = "left")

    panel1 = plt.axes([leftStart/figWidth, #left
                         bottomMargin/figHeight,    #bottom
                         panelWidth/figWidth,   #width
                         panelHeight/figHeight])     #height
    panelxd = plt.axes([leftStart/figWidth, #left
                         (bottomMargin+panelHeight+0.5)/figHeight,    #bottom
                         panelWidth/figWidth,   #width
                         dpanel_width/figHeight])     #height
    panelyd = plt.axes([(leftStart+panelWidth + 0.5)/figWidth, #left
                         bottomMargin/figHeight,    #bottom
                         dpanel_width/figWidth,   #width
                         panelHeight/figHeight])     #height

    # panel2 will host the position-based plot
    plt.gcf().text((secondLeftStart-0.75)/figWidth,
                    (bottomMargin+panelHeight+0.25)/figHeight,
                   "b", weight = "bold",
                   fontsize = 20, ha = "left")
    panel2 = plt.axes([secondLeftStart/figWidth, #left
                         bottomMargin/figHeight,    #bottom
                         panelWidth/figWidth,   #width
                         panelHeight/figHeight])     #height
    pospanelxd = plt.axes([secondLeftStart/figWidth, #left
                         (bottomMargin+panelHeight+0.5)/figHeight,    #bottom
                         panelWidth/figWidth,   #width
                         dpanel_width/figHeight])     #height
    pospanelyd = plt.axes([(secondLeftStart+panelWidth + 0.5)/figWidth, #left
                         bottomMargin/figHeight,    #bottom
                         dpanel_width/figWidth,   #width
                         panelHeight/figHeight])     #height


    panellg = plt.axes([ 0.4/figWidth, #left
                         4/figHeight,    #bottom
                         0.1/figWidth,   #width
                         figHeight/figHeight])     #height

    panel1.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=True,
                        left=False, labelleft=True,
                        right=False, labelright=True,
                        top=False, labeltop=False)
    panel2.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=True,
                        left=False, labelleft=True,
                        right=False, labelright=True,
                        top=False, labeltop=False)
    for xpanels in [panelxd, pospanelxd]:
        xpanels.tick_params(axis='both',which='both',
                            bottom=False, labelbottom=False,
                            left=False, labelleft=False,
                            right=False, labelright=False,
                            top=False, labeltop=False)
    for ypanels in [panelyd, pospanelyd]:
        ypanels.tick_params(axis='both',which='both',
                            bottom=False, labelbottom=False,
                            left=False, labelleft=False,
                            right=False, labelright=False,
                            top=False, labeltop=False)
    panellg.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=False, labelright=False,
                        top=False, labeltop=False)
    # set the panel linewidth thinner
    for this_panel in [panel1, panelxd, panelyd, pospanelxd, pospanelyd]:
        for axis in ['top','bottom','left','right']:
            this_panel.spines[axis].set_linewidth(0.5)
    # turn off the axis spines
    for this_panel in [panelxd, panelyd, pospanelxd, pospanelyd]:
        this_panel.spines['top'].set_visible(False)
        this_panel.spines['right'].set_visible(False)

    # set labels
    plt.gcf().text((leftStart + (panelWidth/2))/figWidth,
                   17.35/figHeight,
                   "{} vs {} synteny (gene indices)".format(xsample, ysample),
                   fontsize = 15, ha = "center")
    plt.gcf().text((secondLeftStart + (panelWidth/2))/figWidth,
                   17.35/figHeight,
                   "{} vs {} synteny (chromosome coordinates)".format(xsample, ysample),
                   fontsize = 15, ha = "center")

    color_list = []
    alpha_list = []
    legend_elements = []
    if prot_to_color_mode:
        # parse the colorset
        LG = oc.LG_db(prot_to_color_name, prot_to_color_dir,
                      prot_to_color_hmm_results)
        plotdf = LG.color_dataframe(plotdf)
        # handle the legend now
        legend_elements += [
          Patch(facecolor=LG.group_to_color[key],
                edgecolor='black', lw = 0,
                label=key) for key in LG.group_to_color
         if key not in ["none", "na", "NA"]]
    else:
        # we're not in prot_to_color_mode, so check to make sure there is a color
        if "color" not in plotdf.columns:
            plotdf["color"] = "#000000"
        else:
            # if there is a color column, make sure that it is a real column
            plotdf["color"] = plotdf["color"].apply(lambda x: x if "#" in x else "#000000")
        # We also have to check if there are groups already.
        #  There shouldn't be but we must check anyway.
        if "gene_group" not in plotdf.columns:
            plotdf["gene_group"] = None
        else:
            # there's already a column, so leave the groups alone
            pass

    legend_elements += [
      Line2D([], [],
            marker='o',
            markerfacecolor=[0,0,0,0.333],
            markeredgecolor='black',
            markeredgewidth = 0,
            markersize = 10,
            linewidth = 0,
            label="p > 0.05"),
      Line2D([], [],
            marker='o',
            markerfacecolor=[0,0,0,1],
            markeredgecolor='black',
            markeredgewidth = 0,
            markersize = 10,
            linewidth = 0,
            label="p <= 0.05")]
    panellg.legend(title="Linkage Group\nColors", handles=legend_elements,
                   fontsize = 10, loc='center left')

    [s.set_visible(False) for s in panellg.spines.values()]
    [t.set_visible(False) for t in panellg.get_xticklines()]
    [t.set_visible(False) for t in panellg.get_yticklines()]

    # first the index section
    x = np.array(plotdf["{}_plotindex".format(xsample)])
    y = np.array(plotdf["{}_plotindex".format(ysample)])
    xmax = max(plotdf["{}_plotindex".format(xsample)])
    ymax = max(plotdf["{}_plotindex".format(ysample)])

    color_list = [oc.h2r(x) for x in plotdf["color"]]
    #print(color_list)
    alpha_list = [[0.25] if x > 0.05 else [0.9] for x in plotdf["break_FET"]]
    # make the black dots a little lighter
    alpha_list = [[0.75] if (plotdf["break_FET"][i] <= 0.05 and plotdf["color"][i] == "#000000") else alpha_list[i] \
                  for i in range(len(alpha_list))]
    composite_color = np.array([color_list[i] + alpha_list[i]
                                for i in range(len(color_list))])

    posx = np.array(plotdf["{}_plotpos".format(xsample)])
    posy = np.array(plotdf["{}_plotpos".format(ysample)])

    pos_xmax = sum([s2c2s[xsample][x] for x in
                plotdf["{}_scaf".format(xsample)].unique()
                if x in s2c2s[xsample]])
    pos_ymax = sum([s2c2s[ysample][x] for x in
                plotdf["{}_scaf".format(ysample)].unique()
                if x in s2c2s[ysample]])
    # sort so the colors are on top
    # make two sets of values to use to make sure colors are on top
    zorder = [-999     if c == "#000000" else 1 for c in list(plotdf["color"])]
    order = np.argsort(zorder)

    panel1.scatter(x[order], y[order], color = composite_color[order],
                   ec = None, s=6, linewidths = 0)
    panel2.scatter(posx[order], posy[order], color = composite_color[order],
                   ec = None, s=6, linewidths = 0)
    # set mins and max for index-panel
    panel1.set_xlim([0, xmax])
    panel1.set_ylim([0, ymax])
    # for position-panel
    panel2.set_xlim([0, pos_xmax])
    panel2.set_ylim([0, pos_ymax])

    # now that we've plotted everything, save the dataframe
    # Header needs to be True, otherwise it will not be read in correctly for other tools
    plotdf.to_csv(plotdf_file, sep ="\t", header = True, index = None)

    # plot vertical x lines
    xlabelpos   = [] # x scaf name positions
    xlabels     = [] # x scaf labels
    ylabelpos   = []
    ylabels     = []
    xlinepos    = [] # delimits the scafs
    ylinepos    = []
    xbreakpos   = [] # delimits where the breaks are
    ybreakpos   = []
    xcountpos   = [] # positions of where gene counts are
    xcountlabel = [] # how many genes are in each scaffold
    ycountpos   = []
    ycountlabel = []
    pos_xlabelpos   = [] # x scaf name positions
    pos_xlabels     = [] # x scaf labels
    pos_ylabelpos   = []
    pos_ylabels     = []
    pos_xlinepos    = [] # delimits the scafs
    pos_ylinepos    = []
    pos_xbreakpos   = [] # delimits where the breaks are
    pos_ybreakpos   = []
    pos_xcountpos   = [] # positions of where gene counts are
    pos_xcountlabel = [] # how many genes are in each scaffold
    pos_ycountpos   = []
    pos_ycountlabel = []

    # first do this for the plot with the indices
    # X-AXIS LABELS AND VERTICAL LINES
    tempgb = plotdf.sort_values(by=["{}_plotindex".format(xsample)]).groupby(by="{}_scaf".format(xsample))
    # we can get these based on index
    for xscaf in xsorter:
        # sort the indices
        minpos = tempgb.get_group(xscaf)["{}_plotindex".format(xsample)].min()
        maxpos = tempgb.get_group(xscaf)["{}_plotindex".format(xsample)].max()
        midpos = ((maxpos - minpos)/2) + minpos
        xlabelpos.append(midpos)
        xlabels.append(xscaf)
        # now for the positions
        pos_xlabels.append(xscaf)
        pos_xlabelpos.append(x_to_cumsum[xscaf] + (s2c2s[xsample][xscaf]/2))
    # for the scaf labels for the chrom size

    # find vertical lines, get everything except the last, and add 0.5 to each value
    vlinesdf = tempgb.last().reset_index()
    tempxdf = vlinesdf.loc[vlinesdf["{}_scaf".format(xsample)].isin(xsorter[:-1]), ]
    xlinepos = list(tempxdf["{}_plotindex".format(xsample)] + 0.5)
    pos_xlinepos = xsorter_to_cumsum[1::]
    # get the positions and quantities of genes
    xcountpos = [0] + list(vlinesdf.sort_values(by=["{}_plotindex".format(xsample)], 
                                                ascending = True)["{}_plotindex".format(xsample)])
    xcountlabel = [str(x) for x in xcountpos]
    # now indices
    pos_xcountpos = [0] + list(vlinesdf.sort_values(by=["{}_plotindex".format(xsample)], 
                                                ascending = True)["{}_plotpos".format(xsample)])
    pos_xcountlabel = xcountlabel

    # Y-AXIS LABELS AND HORIZONTAL LINES
    tempgb = plotdf.sort_values(by=["{}_plotindex".format(ysample)]).groupby(by="{}_scaf".format(ysample))
    for yscaf in ysorter:
        minpos = tempgb.get_group(yscaf)["{}_plotindex".format(ysample)].min()
        maxpos = tempgb.get_group(yscaf)["{}_plotindex".format(ysample)].max()
        midpos = ((maxpos - minpos)/2) + minpos
        ylabelpos.append(midpos)
        ylabels.append(yscaf)
        # now for the positions
        pos_ylabels.append(yscaf)
        pos_ylabelpos.append(y_to_cumsum[yscaf] + (s2c2s[ysample][yscaf]/2))

    #find horizontal lines, get everything except the last, and add 0.5 to each value
    hlinesdf = tempgb.last().reset_index()
    tempydf = hlinesdf.loc[hlinesdf["{}_scaf".format(ysample)].isin(ysorter[:-1]), ]
    ylinepos = list(tempydf["{}_plotindex".format(ysample)] + 0.5)
    pos_ylinepos = ysorter_to_cumsum[1::]
    # get the positions and quantities of genes
    ycountpos = [0] + list(hlinesdf["{}_plotindex".format(ysample)])
    ycountlabel = ycountpos
    # now indices
    pos_ycountpos = [0] + list(hlinesdf.sort_values(by=["{}_plotindex".format(ysample)], 
                                                ascending = True)["{}_plotpos".format(ysample)])
    pos_ycountlabel = [0] + list(hlinesdf.sort_values(by=["{}_plotindex".format(ysample)], 
                                                ascending = True)["{}_plotindex".format(ysample)])

    # FIND VERTICAL BREAKS
    tempxdf = plotdf.sort_values(by=["{}_plotindex".format(xsample)]).groupby(
        by="{}_breakchrom".format(xsample)).last().reset_index()
    xbreakpos = [x for x in list(tempxdf["{}_plotindex".format(xsample)] + 0.5)
                 if ((x not in xlinepos) and (x < xmax))]
    tempxdf = tempxdf.loc[~tempxdf["{}_breakchrom".format(xsample)].str.contains("end"), ]
    pos_xbreakpos = [x for x in list(tempxdf["{}_plotpos".format(xsample)] + 0.5)
                 if ((x not in pos_xlinepos) and (x < pos_xmax))]

    # FIND HORIZONTAL BREAKS
    tempydf = plotdf.sort_values(by=["{}_plotindex".format(ysample)]).groupby(
        by=["{}_breakchrom".format(ysample)]).last().reset_index()
    ybreakpos = [y for y in list(tempydf["{}_plotindex".format(ysample)] + 0.5)
                 if ((y not in ylinepos) and (y < ymax))]
    tempydf = tempydf.loc[~tempydf["{}_breakchrom".format(ysample)].str.contains("end"), ]
    pos_ybreakpos = [y for y in list(tempydf["{}_plotindex".format(ysample)])
                 if ((y not in pos_ylinepos) and (y < pos_ymax))]

    linecolor="#8FA0B5"
    linecolor = [0, 0.169, 0.357, 0.5]
    #plot vertical lines and breaks
    for xval in xlinepos:
        panel1.axvline(x=xval, color=linecolor, lw=0.5)
    for xval in pos_xlinepos:
        panel2.axvline(x=xval, color=linecolor, lw=0.5)
    for xval in xbreakpos:
        panel1.axvline(x=xval, color=[0,0,0,0.35], lw=0.5, linestyle="dotted")
    for xval in pos_xbreakpos:
        panel2.axvline(x=xval, color=[0,0,0,0.35], lw=0.5, linestyle="dotted")
    #plot horizontal lines and breaks
    for yval in ylinepos:
        panel1.axhline(y=yval, color=linecolor, lw=0.5)
    for yval in pos_ylinepos:
        panel2.axhline(y=yval, color=linecolor, lw=0.5)
    for yval in ybreakpos:
        panel1.axhline(y=yval, color=[0,0,0,0.35], lw=0.5, linestyle="dotted")
    for yval in pos_ybreakpos:
        panel2.axhline(y=yval, color=[0,0,0,0.35], lw=0.5, linestyle="dotted")

    # plot x axis labels
    #index plotting
    panel1.tick_params(bottom=True, labelbottom = True, top = False, labeltop = False)
    panel1.set_xticks(xlabelpos)
    panel1.set_xticklabels(xlabels, fontsize=8, rotation=90)
    panel1.set_xlabel(xsample + " scaffolds")
    #position plotting
    panel2.tick_params(bottom=True, labelbottom = True, top = False, labeltop = False)
    panel2.set_xticks(pos_xlabelpos)
    panel2.set_xticklabels(pos_xlabels, fontsize=8, rotation=90)
    panel2.set_xlabel(xsample + " scaffolds")

    axT = panel1.twiny()
    axT.tick_params(top = True, labeltop = True)
    axT.set_xticks(xcountpos)
    axT.set_xticklabels(xcountlabel, fontsize = 8, rotation=90)
    axTT = panel2.twiny()
    axTT.tick_params(top = True, labeltop = True)
    axTT.set_xticks(pos_xcountpos)
    axTT.set_xticklabels(pos_xcountlabel, fontsize = 8, rotation=90)

    # plot y axis labels
    # index plotting
    panel1.tick_params(left=True, labelleft = True, right = False, labelright = False)
    panel1.set_yticks(ylabelpos)
    panel1.set_yticklabels(ylabels, fontsize=8, rotation=0)
    panel1.set_ylabel(ysample + " scaffolds")
    #position plotting
    panel2.tick_params(left=True, labelleft = True, right = False, labelright = False)
    panel2.set_yticks(pos_ylabelpos)
    panel2.set_yticklabels(pos_ylabels, fontsize=8, rotation=0)
    panel2.set_ylabel(ysample + " scaffolds")

    # for indices
    axR = panel1.twinx()
    axR.tick_params(right = True, labelright = True)
    axR.set_yticks(ycountpos)
    axR.set_yticklabels(ycountlabel, fontsize = 8)
    # for positions
    axRR = panel2.twinx()
    axRR.tick_params(right = True, labelright = True)
    axRR.set_yticks(pos_ycountpos)
    axRR.set_yticklabels(pos_ycountlabel, fontsize = 8)

    ## turn on x-axis ticks on the Dx plot
    #newarrlabels = [round(x/1000000, 1) for x in newarr]
    #panelxd.tick_params(top=True, labeltop=True)
    #panelxd.set_xticks(newarr)
    #panelxd.set_xticklabels(newarrlabels, fontsize=8, rotation=90)
    #panelxd.xaxis.set_label_position("top")
    #panelxd.set_xlabel("Mb")

    ## turn on y-axis ticks on the Dy plot
    #newarrlabels = [round(x/1000000, 1) for x in newarr]
    #panelyd.tick_params(right=True, labelright=True)
    #panelyd.set_yticks(newarr)
    #panelyd.set_yticklabels(newarrlabels, fontsize=8)
    #panelyd.yaxis.set_label_position("right")
    #panelyd.set_ylabel("Mb")

    # set the x and y labels on Dy and Dx
    # by indices
    panelxd.bar(x = plotdf["{}_plotindex".format(xsample)],
                height=plotdf["{}_D".format(xsample)],
                width = 1,
                align = "center", lw=0, color="blue", zorder = 2)
    panelxd.set_xlim([0,xmax])
    panelxd.set_ylabel('Dx', fontsize=10)

    # now by position
    plotdf = plotdf.sort_values(by=["{}_plotpos".format(xsample)]).reset_index()
    plotdf["{}_Dwidths".format(xsample)] = plotdf["{}_plotpos".format(xsample)].diff().fillna(0)

    pospanelxd.bar(x = plotdf["{}_plotpos".format(xsample)],
                height=plotdf["{}_D".format(xsample)],
                width = plotdf["{}_Dwidths".format(xsample)],
                align = "center", lw=0, color="blue", zorder = 2)
    pospanelxd.set_xlim([0,pos_xmax])
    pospanelxd.set_ylabel('Dx', fontsize=10)

    # by indices
    panelyd.barh(y = plotdf["{}_plotindex".format(ysample)],
                width=plotdf["{}_D".format(ysample)],
                height = 1,
                align = "center", lw=0, color="blue", zorder = 2)
    panelyd.set_ylim([0,ymax])
    panelyd.set_xlabel('Dy', fontsize=10)

    # now by position
    plotdf = plotdf.sort_values(by=["{}_plotpos".format(ysample)]).reset_index()
    plotdf["{}_Dwidths".format(ysample)] = plotdf["{}_plotpos".format(ysample)].diff().fillna(0)

    pospanelyd.barh(y = plotdf["{}_plotpos".format(ysample)],
                width=plotdf["{}_D".format(ysample)],
                height = plotdf["{}_Dwidths".format(ysample)],
                align = "center", lw=0, color="blue", zorder = 2)

    pospanelyd.set_ylim([0,pos_ymax])
    pospanelyd.set_xlabel('Dy', fontsize=10)

    # TABLE OF THE CHROMOSOME PIECES and significance
    pd.set_option('display.float_format', '{:.2E}'.format)
    table_title_font_size = 11
    rowheight = 0.14
    cutoff = 0.1

    # table measurements
    table_left_margin = 2.15
    left_table_width = 3.9
    middle_table_width = 2.9
    right_table_width = 3.4
    table_spacing = 0.4
    cde_h_offset = 0.15
    cde_v_offset = 0.1

    gb =  plotdf.loc[plotdf["break_FET"] <= cutoff, ].groupby(by=["{}_breakchrom".format(xsample), "{}_breakchrom".format(ysample), "break_FET"]).size().reset_index()
    gb.columns = ["{} scaf".format(xsample), "{} scaf".format(ysample), "p (FET)", "Count"]
    gb = gb.sort_values(by=["Count", "{} scaf".format(xsample), "{} scaf".format(ysample)],
                        ascending = False).reset_index(drop=True)
    gb["p (FET)"] = gb["p (FET)"].apply(lambda x: "{:.2e}".format(x))
    gb["{} scaf".format(xsample)] = gb["{} scaf".format(xsample)].apply(
        lambda x: " : ".join(str(x).split(":")))
    gb["{} scaf".format(ysample)] = gb["{} scaf".format(ysample)].apply(
        lambda x: " : ".join(str(x).split(":")))

    plt.gcf().text((table_left_margin + left_table_width + (2*table_spacing) + \
                    middle_table_width + (right_table_width/2))/figWidth,
                    (bottomMargin-1.25)/figHeight,
                   "{} vs {} FET, chrom. pieces".format(xsample, ysample),
                   fontsize = table_title_font_size, ha = "center")
    plt.gcf().text((table_left_margin + left_table_width + (2*table_spacing) + \
                    middle_table_width - cde_h_offset)/figWidth,
                    (bottomMargin-1.25 + cde_v_offset)/figHeight,
                   "e", weight = "bold",
                   fontsize = 20, ha = "left")
    tableheight = rowheight*(len(gb)+1)

    chrompiecetable = plt.axes([(table_left_margin + left_table_width + table_spacing + middle_table_width + table_spacing)/figWidth, #left
                         (bottomMargin - tableheight- 1.35)/figHeight,    #bottom
                         right_table_width/figWidth,   #width
                         tableheight/figHeight])     #height
    chrompiecetable.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=False, labelright=False,
                        top=False, labeltop=False)
    if len(gb.values) > 0:
        thistab_chrompiece = chrompiecetable.table(rowLabels=["{}".format(x) for x in gb.index],
                       cellText=gb.values,
                       colLabels=gb.columns,
                       cellLoc = "left",
                       colLoc  = "left",
                       edges = "horizontal",
                       bbox = [0,0,1,1])
        thistab_chrompiece.auto_set_column_width(col=list(range(len(gb.columns))))
        thistab_chrompiece.auto_set_font_size(False)
        thistab_chrompiece.set_fontsize(7)


    # TABLE OF WHOLE CHROMOSOME SIGNIFICANCE
    plt.gcf().text((table_left_margin + left_table_width + table_spacing + (middle_table_width*0.5))/figWidth,
                    (bottomMargin-1.25)/figHeight,
                   "{} vs {} FET, whole chroms.".format(xsample, ysample),
                   fontsize = table_title_font_size, ha = "center")
    plt.gcf().text((table_left_margin + left_table_width + (table_spacing) - cde_h_offset)/figWidth,
                    (bottomMargin-1.25 + cde_v_offset)/figHeight,
                   "d", weight = "bold",
                   fontsize = 20, ha = "left")
    gb2 = plotdf.loc[plotdf["whole_FET"] <= cutoff, ].groupby(by=["{}_scaf".format(xsample), "{}_scaf".format(ysample), "whole_FET"]).size().reset_index()
    gb2.columns = ["{} scaf".format(xsample), "{} scaf".format(ysample), "p (FET)", "Count"]
    gb2 = gb2.sort_values(by=["Count", "{} scaf".format(xsample), "{} scaf".format(ysample)],
                        ascending = False).reset_index(drop=True)
    gb2["p (FET)"] = gb2["p (FET)"].apply(lambda x: "{:.2e}".format(x))
    tableheight = rowheight*(len(gb2)+1)
    chromwholetable = plt.axes([(table_left_margin+left_table_width +table_spacing)/figWidth, #left
                         (bottomMargin - tableheight- 1.35)/figHeight,    #bottom
                         (middle_table_width)/figWidth,   #width
                         tableheight/figHeight])     #height
    chromwholetable.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=False, labelright=False,
                        top=False, labeltop=False)
    if len(gb2.values) > 0:
        thistab_chromwhole = chromwholetable.table(
                                rowLabels=["{}".format(x) for x in gb2.index],
                                cellText=gb2.values,
                                colLabels=gb2.columns,
                                cellLoc = "left",
                                colLoc  = "left",
                                edges = "horizontal",
                                bbox = [0,0,1,1])
        thistab_chromwhole.auto_set_column_width(col=list(range(len(gb2.columns))))
        thistab_chromwhole.auto_set_font_size(False)
        thistab_chromwhole.set_fontsize(7)

    # TABLE OF COLOR GROUP RESULTS
    plt.gcf().text((table_left_margin + (left_table_width/2))/figWidth,
                    (bottomMargin-1.25)/figHeight,
                   "{} vs {} groupings".format(xsample, ysample),
                   fontsize = table_title_font_size, ha = "center")
    plt.gcf().text((table_left_margin - cde_h_offset)/figWidth,
                    (bottomMargin-1.25 + cde_v_offset)/figHeight,
                   "c", weight = "bold",
                   fontsize = 20, ha = "left")
    gb3 = plotdf.loc[plotdf["whole_FET"] <= 1, ].groupby(by=["{}_scaf".format(xsample), "{}_scaf".format(ysample),"whole_FET", "gene_group"]).size().reset_index()
    gb3.columns = ["{} scaf".format(xsample), "{} scaf".format(ysample), "p (FET)", "Group", "Count"]
    gb3 = gb3.loc[gb3["Count"] >= 5, ]

    if prot_to_color_mode:
        # if we have color mappings from the gene group we're using, use those colors
        gb3["Color"] = gb3["Group"].map(LG.group_to_color)
    else:
        # otherwise, generate the group to color just for this block from our plotting df
        tempgroup_to_color_df = plotdf.groupby(
                 by=["gene_group"])["color"].value_counts(
                 ascending = False).rename("Counts").reset_index(
                 ).drop_duplicates("gene_group")
        tempgroup_to_color = dict(zip(tempgroup_to_color_df["gene_group"],
                                      tempgroup_to_color_df["color"]))
        gb3["Color"] = gb3["Group"].map(tempgroup_to_color)

    gb3 = gb3.sort_values(by=["{} scaf".format(xsample), "{} scaf".format(ysample), "p (FET)", "Count"],
                          ascending = [True, True, True, False]).reset_index(drop=True)
    gb3["p (FET)"] = gb3["p (FET)"].apply(lambda x: "{:.2e}".format(x))
    # swap the columns
    gb3 = gb3[["{} scaf".format(xsample), "{} scaf".format(ysample), "p (FET)", "Count", "Group", "Color"]]

    color_matrix = []
    for index, row in gb3.iterrows():
        # example from
        #  https://stackoverflow.com/questions/46663911
        # the_table[(1, 0)].set_facecolor("#56b5fd")
        base = ["w"] * len(gb3.columns)
        thiscolor = row["Color"] if row["Color"] != "#000000" else "w"
        base[-1] = thiscolor
        base[-2] = thiscolor
        color_matrix.append(base)

    tableheight = rowheight*(len(gb3)+1)
    colortable = plt.axes([table_left_margin/figWidth, #left
                         (bottomMargin - tableheight- 1.35)/figHeight,    #bottom
                         left_table_width/figWidth,   #width
                         tableheight/figHeight])     #height
    colortable.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=False, labelright=False,
                        top=False, labeltop=False)
    if len(gb3.values) > 0:
        thistab_color = colortable.table(
                           rowLabels=["{}".format(x) for x in gb3.index],
                           cellText=gb3.values,
                           cellColours=color_matrix,
                           colLabels=gb3.columns,
                           cellLoc = "left",
                           colLoc  = "left",
                           #edges = "horizontal", #bug where this doesn't work with cell colors https://stackoverflow.com/questions/67890401/
                           bbox = [0,0,1,1])
        thistab_color.auto_set_column_width(col=list(range(len(gb3.columns))))
        thistab_color.auto_set_font_size(False)
        thistab_color.set_fontsize(7)


    #print(thistab_color._cells)
    for rowi in range(1,len(gb3)+1):
        for coli in [len(gb3.columns)-2, len(gb3.columns)-1]:
            inv = oc.inverse_color(gb3.loc[rowi-1, "Color"])
            if inv == "#ffffff":
                thistab_color._cells[(rowi, coli)]._text.set_color(inv)

    captionax = plt.axes([1.0/figWidth, #left
                          1.0/figHeight,    #bottom
                          (figWidth-2)/figWidth,   #width
                          1/figHeight])     #height
    captionax.tick_params(axis='both',which='both',
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=False, labelright=False,
                        top=False, labeltop=False)

    # top text
    top_text = "{}: {} orthologs on {} {} scaffolds and {} {} scaffolds.".format(
                       orthology_method, len(plotdf),
                       len(plotdf["{}_scaf".format(xsample)].unique()), xsample,
                       len(plotdf["{}_scaf".format(ysample)].unique()), ysample)
    if prot_to_color_mode:
        top_text += " Colored by {}.".format(LG.color_method)
    plt.gcf().text(0.1/figHeight,
                   (figHeight-0.1)/figHeight,
                   top_text,
                   weight = "bold",
                   fontsize = 13, ha = "left", va="top")
    # caption
    if orthology_method == "diamond":
        orthology_method = "reciprocal best diamond blastp match between two species"
    elif orthology_method == "blastp":
        orthology_method = "reciprocal best blastp match between two species"
    t = [
        "(a) depicts the chromosomes/scaffolds of {} (x-axis) plotted".format(xsample),
        "against the chromosomes/scaffolds of {} (y_axis).".format(ysample),
        "Each dot in the plot represents an ortholog, specifically a {}.".format(orthology_method),
        "The unit of the x- and y-axes are the number of orthologous proteins between these two species:",
        "{} orthologs found between {} {} scaffolds and {} {} scaffolds.".format(
            len(plotdf),
            len(plotdf["{}_scaf".format(xsample)].unique()), xsample,
            len(plotdf["{}_scaf".format(ysample)].unique()), ysample),
        "If there are chromosome breaks, Fisher's exact test is used to calculate the",
        "significance of the interactions between the sub-chromosomal pieces. Otherwise",
        "Fisher's exact test is calculated on whole chromosomes.",
        "The opacity of the dots depict the significance from Fisher's exact test.",
        "Dots that are a solid color are in cells with a FET p-value less than or equal to 0.05.",
        "Dots that are translucent are in cells with a FET p-value greater than 0.05.",
        "The Dx and Dy values depict places where there may be sudden breaks in synteny.",
        "See the supplementary information the following paper for more information on Dx and Dy: Simakov, Oleg, et al.",
        "Deeply conserved synteny resolves early events in vertebrate evolution. Nature Ecology & Evolution 4.6 (2020): 820-830.",
        "(b) depicts the same information as panel a, but it is plotted in the",
        "organisms chromosome basepair coordinates rather than gene index.",
        "This is useful for visualizing gene-poor regions of the chromosomes",
        "(c) shows which gene groups are most prevalent in each chromosome pair.",
        "The FET p-value in this table corresponds to the whole-chormosome FET p-value",
        ", and is not a FET value of the correlation between the chromosome pair and the gene group.",
        "(d) shows chromosome-scale significance values.",
        "This table shows the same information as c, but does not factor in gene group information.",
        "(e) shows the FET p-values of the sub-chromosomal compartments. Useful for",
        "showing if single arms of chromosomes are correlated with other regions."
         ]
    t = " ".join(t)

    captext = captionax.text(0, 1, t, fontsize = 10, ha="left", va="top", wrap=True)

    for this_panel in [chrompiecetable, chromwholetable, colortable, captionax]:
        this_panel.spines['top'].set_visible(False)
        this_panel.spines['left'].set_visible(False)
        this_panel.spines['bottom'].set_visible(False)
        this_panel.spines['right'].set_visible(False)

    plt.savefig(synplot)
