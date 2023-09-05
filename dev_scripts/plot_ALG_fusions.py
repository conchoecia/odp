#!/usr/bin/env python

"""
author: dts
github: @conchoecia
date: September 2023

This script plots the size of ALGs (x-axis) against the number of colocalized ALGs (the fusion part of FWM) (y-axis).
"""

import argparse
import numpy as np
from sklearn.decomposition import PCA
import os
import pandas as pd
import sys


# plotting options
import matplotlib.pyplot as plt
#import odp_plotting_functions as odp_plot


def parse_args():
    """
    The args that we need to parse are the following:
      - directory: the directory where the ALG files are saved
      - ALGname: the name of the ALG database - use this to correctly identify the columns in the rbh file
      - prefix: the prefix of the rbh files. Also potentially includes the directory for the rbh files.
      - ALG_rbh: the name of the rbh file that contains the ALG information. This is used to get the ALG names, sizes, and colors.
    """
    parser = argparse.ArgumentParser(description='Plot the size of ALGs against the number of colocalized ALGs')
    parser.add_argument('-d', '--directory', type=str, required=True, help='The directory where the ALG files are saved')
    parser.add_argument('-a', '--ALGname', type=str, required=True, help='The name of the ALG database')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='The prefix of the rbh files. Also potentially includes the directory for the rbh files.')
    parser.add_argument('-r', '--ALG_rbh', type=str, required=True, help='The name of the rbh file that contains the ALG information. This is used to get the ALG names, sizes, and colors.')

    args = parser.parse_args()
    # make sure that the directory exists
    if not os.path.exists(args.directory):
        raise ValueError('The directory {} does not exist'.format(args.directory))
    return args

def parse_rbh_file(rbh_file, ALGname) -> pd.DataFrame:
    """
    Reads in an rbh file from ALGs and returns the following dataframe:


     ALGname  Color     Size
     Qb       #C72480   12
     Qc       #DCC0F3   14
     C2       #387FB2   18
     Qd       #94C47F   22
      R       #F4B93E   24
     Qa       #78A6AF   30
     A2       #8B4E67   41
     B3       #FA9A26   46
     O2       #AB5BA8   46
     Eb       #B76BED   47
    A1b       #C33D53   51
     J1       #54AB53   54
     O1       #FBD76C   55
     J2       #E64657   66
      P       #C33E51   78
     B2       #1F779A   86
      I       #3425FB   90
     B1       #2F54E3   95
      M       #A45530  102
      L       #7DC29F  104
      N       #D8BE3C  107
     Ea       #AB7E26  115
      K       #170B88  119
      H       #F04C08  135
      G       #E97B4A  138
     C1       #B07DF4  142
      F       #9B6870  145
      D       #47957F  172
    A1a       #4DB5E3  207

    The dataframe will be used later to determine plotting parameters.
    """
    df = pd.read_csv(rbh_file, sep='\t')
    # First make sure that there are columns called "gene_group", "color". We will only use these.
    if not ("gene_group" in df.columns and "color" in df.columns):
        raise IOError("The rbh file {} does not have the correct columns".format(rbh_file))
    # just subset these columns since we only need them
    df = df[["gene_group", "color"]]
    sizemap = df.groupby("gene_group").size()
    # groupby the "gene_group" column, then reconstitute with the size and the most common color for that group
    df = df.groupby("gene_group").agg({"color": lambda x: x.value_counts().index[0]})
    df = df.reset_index()
    # use the sizemap to get the size of each ALG. Use the gene_group as the key to look up the value in sizemap
    df["size"] = df["gene_group"].map(sizemap)
    # rename the columns {"gene_group": "ALGname", "color": "Color", "size": "Size"}
    df = df.rename(columns={"gene_group": "ALGname", "color": "Color", "size": "Size"})
    # sort by increasing size
    df = df.sort_values(by="Size", ascending=True)
    df = df.reset_index(drop=True)
    return df

def parse_ALG_fusions(list_of_rbh_files, ALG_df, ALGname) -> pd.DataFrame:
    """
    Takes in a list of rbh files and plots the size of ALGs (increasing) against a violin plot of the number
    of fusions that this ALG participates in for all of the species found in the RBH file.

    Outputs the following data table. The columns are redundant to facilitate plotting.
    The output will be a pandas dataframe.

    ALGname | Species | Fused_with
    G       | PMA     | []
    R       | PMA     | ["Qa", "Qb"]
    Qa      | PMA     | ["R", "Qb"]
    Qb      | PMA     | ["R", "Qa"]
    """
    # We must structure the data in a way that can be output and used for plotting later.
    #  Each entry will be formatted like this, and will be a row in a dataframe later: {"ALGname": str, "Species": str, "Fused_with": str}
    entries = []

    # iterate through the rbh files
    for thisfile in list_of_rbh_files:
        # read in the rbh file as a pandas dataframe
        rbh_df = pd.read_csv(thisfile, sep='\t')

        # First make sure that there are columns that start with ALGname and end with "_scaf", "_gene", and "_pos"
        # If not, then something is wrong with this file and we should raise an error.
        col1 = "{}_scaf".format(ALGname)
        col2 = "{}_gene".format(ALGname)
        col3 = "{}_pos".format( ALGname)
        if (col1 not in rbh_df.columns) or (col2 not in rbh_df.columns) or (col3 not in rbh_df.columns):
            raise IOError("The rbh file {} does not have the correct columns for the ALG {}".format(thisfile, ALGname))

        # now get the other species name that is not the ALG name
        species = [col.replace("_scaf","") for col in rbh_df.columns if col.endswith("_scaf") and (not col.startswith(ALGname))][0]
        # check that a _scaf, _gene, and _pos column exists in the database for this species
        col1 = "{}_scaf".format(species)
        col2 = "{}_gene".format(species)
        col3 = "{}_pos".format( species)

        if (col1 not in rbh_df.columns) or (col2 not in rbh_df.columns) or (col3 not in rbh_df.columns):
            raise IOError("The rbh file {} does not have the correct columns for the species {}".format(thisfile, species))

        # get the rows where whole_FET <= 0.05
        rbh_df = rbh_df[rbh_df['whole_FET'] <= 0.005]
        # now just keep the columns
        keep_these_columns = ["rbh", "gene_group", "color",
                              "{}_scaf".format(ALGname),
                              "{}_scaf".format(species)]
        rbh_df = rbh_df[keep_these_columns]

        # this is used to keep track of things to get rid of later
        ALGs_seen   = set()

        if len(rbh_df) > 0:
            # there are some ALGs that are significantly correlated with the chromosomes
            grouped = rbh_df.groupby("{}_scaf".format(species))
            # iterate through the groups and for each chromosome, get the other gene groups colocalized on the same chromosome
            for name, group in grouped:
                # get all of the ALGs that are significantly associated with this chromosome
                this_ALG_list = group["gene_group"].unique()
                for thisALG in this_ALG_list:
                    others = [x for x in this_ALG_list if x != thisALG]
                    thisentry = {"Species": species, "ALGname": thisALG,  "Fused_with": others}
                    ALGs_seen.add(thisALG)
                    entries.append(thisentry)

        # There now will be some ALGs that are not significantly correlated with the chromosomes
        #  These will be the ones that are not in the ALGs_seen set
        #  We need to make entries for these, too.
        for thisALG in ALG_df["ALGname"]:
            if thisALG not in ALGs_seen:
                thisentry = {"Species": species, "ALGname": thisALG,  "Fused_with": []}
                entries.append(thisentry)

        # done. next file.

    # make a df of the entries
    df = pd.DataFrame(entries)
    df["color"] = df["ALGname"].map(ALG_df.set_index("ALGname")["Color"])
    df["fused_quantity"] = df["Fused_with"].apply(lambda x: len(x))
    return df

def plot_ALG_fusions(Fusion_df, ALG_df, ALGname, outprefix=None):
    """
    produces a plot of the ALG fusions
    """
    # make a df of the entries. Groupy ALGname and make a list of the values of fused_quantity
    plotdf = Fusion_df.groupby("ALGname").agg({"fused_quantity": lambda x: list(x)})
    plotdf = plotdf.reset_index()
    plotdf["size"] = plotdf["ALGname"].map(ALG_df.set_index("ALGname")["Size"])
    plotdf["color"] = plotdf["ALGname"].map(ALG_df.set_index("ALGname")["Color"])
    plotdf = plotdf.sort_values(by="size", ascending=True)
    plotdf = plotdf.reset_index(drop=True)

    # now we plot
    NUMBER_OF_ROWS = 3
    NUMBER_OF_COLS = 2
    fig, axes = plt.subplots(NUMBER_OF_ROWS, NUMBER_OF_COLS, figsize = (7.5 * NUMBER_OF_COLS, 6 * NUMBER_OF_ROWS))

    # make a title for the whole figure
    fig.suptitle("{} ALGs sizes (x) vs. number of fusions to same chromosome (y) in {} species".format(ALGname, len(Fusion_df["Species"].unique())))

    # plot each row as a violin plot in the first figure
    for i, row in plotdf.iterrows():
        # get the x and y values
        x = i
        y = row["fused_quantity"]
        # get the color
        color = row["color"]
        # ALGsize
        ALGsize = int(row["size"])
        # plot the violin plot - TOP LEFT
        violin_parts = axes[0][0].violinplot(y,
                                             positions=[x], showmeans=False, showmedians=True, showextrema=False)
        violin_parts["cmedians"].set_color("black")
        for pc in violin_parts['bodies']:
            # change the face to be the specified color
            pc.set_facecolor(color)

        # below, plot with no zeros - BOTTOM LEFT
        violin_parts = axes[1][0].violinplot([x for x in y if x > 0],
                                             positions=[x], showmeans=False, showmedians=True, showextrema=False)
        violin_parts["cmedians"].set_color("black")
        for pc in violin_parts['bodies']:
            # change the face to be the specified color
            pc.set_facecolor(color)

        # TOP RIGHT - maxwidth is 5
        violin_parts = axes[0][1].violinplot(y,
                                             positions=[ALGsize],  widths = 5, showmeans=False, showmedians=True, showextrema=False)
        violin_parts["cmedians"].set_color("black")
        for pc in violin_parts['bodies']:
            # change the face to be the specified color
            pc.set_facecolor(color)

        # BOTTOM RIGHT
        violin_parts = axes[1][1].violinplot([x for x in y if x > 0], widths = 5,
                                             positions=[ALGsize], showmeans=False, showmedians=True, showextrema=False)
        violin_parts["cmedians"].set_color("black")
        for pc in violin_parts['bodies']:
            # change the face to be the specified color
            pc.set_facecolor(color)

    # make a linear plot of the number of fusions vs. the ALG size
    Fusion_df["size"] = Fusion_df["ALGname"].map(ALG_df.set_index("ALGname")["Size"])
    # TOP-RIGHT, we have not yet removed the zeros
    x = Fusion_df["size"]
    y = Fusion_df["fused_quantity"]
    m,b = np.polyfit(x, y, 1)
    # determine the r2 of the fit
    r2 = np.corrcoef(x, y)[0,1]**2
    for i in [0,1]:
        if i == 1:
            axes[0][i].plot(x, m*x + b, color="black", linestyle="dashed", linewidth=1)
        # add the equation to the top-right of the TOP-RIGHT plot
        axes[0][i].text(0.05, 0.95, "y = {:.8f}x + {:.2f}".format(m,b), transform=axes[0][i].transAxes, fontsize=10)
        # and add the r2 value below that
        axes[0][i].text(0.05, 0.90, "r2 = {:.8f}".format(r2), transform=axes[0][i].transAxes, fontsize=10)

    # BOTTOM-RIGHT, we have removed the zeros
    Fusion_df_NoZeros = Fusion_df[Fusion_df["fused_quantity"] > 0]
    x = Fusion_df_NoZeros["size"]
    y = Fusion_df_NoZeros["fused_quantity"]
    m,b = np.polyfit(x, y, 1)
    # determine the r2 of the fit
    r2 = np.corrcoef(x, y)[0,1]**2
    for i in [0,1]:
        if i == 1:
            axes[1][i].plot(x, m*x + b, color="black", linestyle="dashed", linewidth=1)
        # add the equation to the top-right of the BOTTOM-RIGHT plot
        axes[1][i].text(0.05, 0.95, "y = {:.8f}x + {:.2f}".format(m,b), transform=axes[1][i].transAxes, fontsize=10)
        # and add the r2 value below that
        axes[1][i].text(0.05, 0.90, "r2 = {:.8f}".format(r2), transform=axes[1][i].transAxes, fontsize=10)

    # *************************************************************
    #   CUMULATIVE FUSIONS
    # *************************************************************
    # now plot the cumulative number of fusions vs. the ALG size
    plotdf["cumulative_fusions"] = plotdf["fused_quantity"].apply(lambda x: sum(x))
    # make a bar chart base on the index for axes[2][0]
    axes[2][0].bar(plotdf.index, plotdf["cumulative_fusions"], color=plotdf["color"])

    # make a bar chart based on the ALG size for axes[2][1]
    axes[2][1].bar(plotdf["size"], plotdf["cumulative_fusions"], color=plotdf["color"])

    # for plot axes[2][1], determine the linear regression
    x = plotdf["size"]
    y = plotdf["cumulative_fusions"]
    m,b = np.polyfit(x, y, 1)
    # determine the r2 of the fit
    r2 = np.corrcoef(x, y)[0,1]**2
    for i in [0,1]:
        if i == 1:
            axes[2][i].plot(x, m*x + b, color="black", linestyle="dashed", linewidth=1)
        axes[2][i].text(0.05, 0.95, "y = {:.8f}x + {:.2f}".format(m,b), transform=axes[2][i].transAxes, fontsize=10)
        axes[2][i].text(0.05, 0.90, "r2 = {:.8f}".format(r2), transform=axes[2][i].transAxes, fontsize=10)
    print(plotdf)

    # set the x-axis labels THESE ARE THE ROWS FOR WHICH WE HAVEN"T REMOVED THE ZEROS
    for i in [0,1,2]:
        axes[i][0].set_xticks(range(len(plotdf)))
        axes[i][0].set_xticklabels(plotdf["ALGname"], rotation=90)

    # set the x-axis labels for the plots on the right, plotting by ALG size
    # THIS IS THE SET FOR WHICH WE HAVE REMOVED ZEROS
    for i in [0,1,2]:
        axes[i][1].set_xticks(plotdf["size"])
        axes[i][1].set_xticklabels(plotdf["ALGname"], rotation=90)


    # set titles for each of the panels
    axes[0][0].set_title("ALG size (ranked) vs. number of fusions (all)")
    axes[1][0].set_title("ALG size (ranked) vs. number of fusions (no zeros)")
    axes[2][0].set_title("ALG size (ranked) vs. cumulative fusions")
    axes[0][1].set_title("ALG size (increasing) vs. number of fusions (all)")
    axes[1][1].set_title("ALG size (increasing) vs. number of fusions (no zeros)")
    axes[2][1].set_title("ALG size (increasing) vs. cumulative fusions")

    # save this as a pdf
    outfile = ""
    if outprefix is not None:
        outfile = "{}_ALG_fusions.pdf".format(outprefix)
    else:
        outfile = "ALG_fusions.pdf"
    plt.savefig(outfile, bbox_inches='tight')

def main():
    # parse the args
    args = parse_args()

    # get the rbh files in the directory
    rbh_files = [os.path.join(args.directory, f)
                 for f in os.listdir(args.directory)
                 if f.endswith('.rbh')]

    # first we need to read through all of the rbh files to get all of the possible ALGs to plot
    ALGdf = parse_rbh_file(args.ALG_rbh, args.ALGname)

    ## *************************************************************
    ##  MAKING THE PLOTS ABOUT THE NUMBER OF FUSIONS
    ## *************************************************************
    ## parse the fusions
    #Fusion_df = parse_ALG_fusions(rbh_files, ALGdf, args.ALGname)
    ## plot the fusions
    #plot_ALG_fusions(Fusion_df, ALGdf, args.ALGname, outprefix=args.prefix)

    # *************************************************************
    #  MAKING THE PCA FOR PER-CHROMOSOME ALG COLOCALIZATION
    # *************************************************************
    entries = []
    for file in rbh_files:
        # read in the file as a pandas dataframe
        df = pd.read_csv(file, sep='\t')
        # get the species name
        species = [col.replace("_scaf","") for col in df.columns if col.endswith("_scaf") and (not col.startswith(args.ALGname))][0]

        # groupby the chromosome names
        grouped = df.groupby("{}_scaf".format(species))
        # iterate through the chromosomes.
        for name, group in grouped:
            # filter out the rows that are not significantly correlated
            group = group[group["whole_FET"] <= 0.005]
            thisentry = {x: 0 for x in ALGdf["ALGname"]}
            thisentry["sp.|.chromosome"] = "{}.|.{}".format(species,name)
            # if this chromosome has any significantly correlated ALGs, then add them to the dictionary
            if len(group) > 0:
                chroms = group["{}_scaf".format(args.ALGname)].unique()
                for chrom in chroms:
                    thisentry[chrom] = 1
            entries.append(thisentry)
    # convert the entries into a dataframe
    df = pd.DataFrame(entries)
    # change index to sp.|.chromosome
    df = df.set_index("sp.|.chromosome")

    # make a df of counts where two columns had values of 1 on the same row
    # counts is a dict of all the chromosomes
    ALG_list = list(sorted(ALGdf["ALGname"].tolist()))
    counts = {}
    for i in range(len(ALG_list)-1):
        for j in range(i+1, len(ALG_list)):
            sorted_list = sorted([ALG_list[i], ALG_list[j]])
            thistuple = (sorted_list[0], sorted_list[1])
            counts[thistuple] = 0
    # go through the rows of the df and count the number of times that two ALGs are colocalized
    for i, row in df.iterrows():
        # get the ALGs that are colocalized on this chromosome
        ALGs = [x for x in row.index if x in ALG_list and row[x] == 1]
        # get all of the combinations of these ALGs
        for i in range(len(ALGs)-1):
            for j in range(i+1, len(ALGs)):
                sorted_list = sorted([ALGs[i], ALGs[j]])
                thistuple = (sorted_list[0], sorted_list[1])
                counts[thistuple] += 1
    # print all the counts that are 0
    for k,v in counts.items():
        if v == 0:
            print("{}: {} occurrences".format(k, v))

    # convert counts to a dataframe where tuple 0 is a col called ALG1, tuple 1 is a cold called ALG2, and the value is the number of times they are colocalized
    countsdf = pd.DataFrame.from_dict(counts, orient='index').reset_index()
    # splot the "index" column into two columns
    countsdf[["ALG1", "ALG2"]] = pd.DataFrame(countsdf["index"].tolist(), index=countsdf.index)
    # drop the index column
    countsdf = countsdf.drop(columns=["index"])
    # rename the columns
    countsdf = countsdf.rename(columns={0: "count"})
    # unmelt the dataframe so that it is a square matrix. All of the columns will be ints
    countsdf = countsdf.pivot(index="ALG1", columns="ALG2", values="count")
    # fill the NaNs to -1
    countsdf = countsdf.fillna(-1)
    # convert all the columns to ints
    countsdf = countsdf.astype(int)
    # flip the columns and rows
    countsdf = countsdf.transpose()
    print(countsdf)

   # remove rows where all values are 0
   # df = df.loc[(df!=0).any(axis=1)]
   # pca = PCA(n_components=2)
   # pca.fit(df)
   # print(pca.components_)
   # print(pca.explained_variance_)

   # df2 = pd.DataFrame(pca.transform(df), columns = ['first', 'second'])
   # print(df)
   # print(df2)
   # # df2.plot.scatter(x = 'first', y = 'second')

   # #plt.show()

   # from mpl_toolkits import mplot3d
   # from mpl_toolkits.mplot3d import Axes3D
   # from mpl_toolkits.mplot3d import proj3d
   # from matplotlib.text import Annotation

   # x = df2['first']
   # y = df2['second']
   # # labels is a list of the df1 index values
   # labels = df.index.values.tolist()

   # # Create the scatter plot
   # fig, ax = plt.subplots()
   # scatter = ax.scatter(x, y, picker=True)
   # # plot the text labels
   # for i, txt in enumerate(labels):
   #     ax.annotate(txt, (x[i], y[i]))


   # # Show the plot
   # plt.show()


if __name__ == '__main__':
    main()