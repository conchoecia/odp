#!/usr/bin/env python

"""
author: dts
github: @conchoecia
date: September 2023

This script plots the size of ALGs (x-axis) against the number of colocalized ALGs (the fusion part of FWM) (y-axis).

It also makes a table of detectable fusion-with-mixing events in all of the rbh files. This table looks like this:
      A1a  A1b  A2  B1  B2  B3
A1a     0   85  53   0  14   1
A1b    85    0  44   1   2   7
A2     53   44   0   0   9   3
B1      0    1   0   0  19   3
B2     14    2   9  19   0   2
B3      1    7   3   3   2   0

Also makes a file that has per-species fusion with mixing. Columns include the genus, species name, NCBI Taxid, NCBI Taxid string,
 and columns where each FWM event is represented.

PREREQUISITES:
  - This script requires the ete3 toolkit to be installed. This can be done with conda:
    https://anaconda.org/conda-forge/ete3
  - You must also preload the taxonomy toolkit with the following commands:
    http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html
    ```
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    ```
"""

import argparse
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import itertools
import os
import pandas as pd
import re
import sys

from PIL import Image
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster import hierarchy

# get the warnings
import warnings
warnings.filterwarnings('error')


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
      - minsig: the minimum significance value for the whole_FET column in the rbh files. This is used to filter the rbh files.
    """
    parser = argparse.ArgumentParser(description='Plot the size of ALGs against the number of colocalized ALGs')
    parser.add_argument('-d', '--directory', type=str,   required=True,   help='The directory where the ALG files are saved')
    parser.add_argument('-a', '--ALGname',   type=str,   required=True,   help='The name of the ALG database')
    parser.add_argument('-p', '--prefix',    type=str,   required=True,   help='The prefix of the rbh files. Also potentially includes the directory for the rbh files.')
    parser.add_argument('-r', '--ALG_rbh',   type=str,   required=True,   help='The name of the rbh file that contains the ALG information. This is used to get the ALG names, sizes, and colors.')
    parser.add_argument('-m', '--minsig',    type=float, default = 0.005, help='The minimum significance value for the whole_FET column in the rbh files. This is used to filter the rbh files.')

    args = parser.parse_args()
    # make sure that the directory exists
    if not os.path.exists(args.directory):
        raise ValueError('The directory {} does not exist'.format(args.directory))
    return args

def hex_to_rgb(hex):
    """
    Converts a hex color to an rgb color. The hex color can be in the form "#FFFFFF" or "FFFFFF".
    The RGB values will be in the range of 0-255.
    """
    hex = hex.lstrip('#')
    hlen = len(hex)
    return tuple(int(hex[i:i+hlen//3], 16) for i in range(0, hlen, hlen//3))

def parse_rbh_file(rbh_file) -> pd.DataFrame:
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

def parse_ALG_fusions(list_of_rbh_files, ALG_df, ALGname, minsig) -> pd.DataFrame:
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

    Inputs:
      - list_of_rbh_files: a list of rbh files
      - ALG_df: a dataframe of the ALG names, colors, and sizes
      - ALGname: the name of the ALG database, for example, BCnSSimakov2022
      - minsig: the minimum significance value for the whole_FET column in the rbh files. The default of this program is 0.005.
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

        # get the rows where whole_FET <= minsig
        rbh_df = rbh_df[rbh_df['whole_FET'] <= minsig]
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
    plt.savefig(outfile)

def image_sp_matrix_to_lineage(taxidstring) -> Image:
    """
    Required loadings:
      - from PIL import Image
      - import numpy as np

    This takes a list called taxidstring. The taxidstrings will be taxids delimited with a semicolon.
        1;131567;2759;33154;33208;6072;33213;33317;120...
        1;131567;2759;33154;33208;6072;33213;33317;120...
        1;131567;2759;33154;33208;6072;33213;33317;120...
        1;131567;2759;33154;33208;6072;33213;33317;120...
        1;131567;2759;33154;33208;6072;33213;33317;120...
    Algorithm:
      - This will first go through all of the taxid strings and figure out the longest one.
      - Then it will construct an image where each row is a species, and each pixel is a taxid.
      - Going from the rightmost taxid, turn the pixel on if that taxid has not yet been seen.
      - Returns the img object to be handled later - probably to be added to another Image object.

    We want to color certain taxonomic units:
      - CLADE         - NCBI Taxid - COLOR
      - Ctenophores   - 10197      -  #54AB53
      - Sponges       - 6040       -  #DCC0F3
      - Cnidarians    - 6073       -  #387FB2
      - Placozoans    - 10226      -  #C72480
      - Protostomes   - 33317      -  #F4B93E
      - Deuterostomes - 33511      -  #78A6AF
    """
    split_taxids = [x.split(";") for x in taxidstring]
    # get the indices of the longest taxid strings, sorted descending
    indices_by_length = sorted(range(len(split_taxids)), key=lambda x: len(split_taxids[x]), reverse=True)
    maxindexlength = len(split_taxids[indices_by_length[0]])

    taxid_to_x_index = {}
    # go through all of the taxid strings and figure out at which x-value they should occur in the matrix
    for thisindex in indices_by_length:
        thisstring = split_taxids[thisindex]
        # if the length of thisstring = maxindexlength, then we don't need to do anything
        # if the length of the taxid string is shorter than the maxindexlength, transform the coordinates to the closest integer in the range of maxindexlength
        if len(thisstring) < maxindexlength:
            # make a dict where the index of the thisstring is transformed to an int close to the maxindexlength
            index_to_transformed_range = {i: int(np.round(np.linspace(0, maxindexlength-1, len(thisstring)))[i]) for i in range(len(thisstring))}
        elif len(thisstring) == maxindexlength:
            index_to_transformed_range = {i: i for i in range(len(thisstring))}
        for i in range(len(thisstring)):
            # Now add the taxid at its transformed index to the dictionary.
            # The key is the taxid, and the value is the transformed index.
            # Only add if we haven't seen it yet
            if thisstring[i] not in taxid_to_x_index:
                taxid_to_x_index[thisstring[i]] = index_to_transformed_range[i]

    # make a numpy matrix of the longest width (x-axis) and the length of the taxidstring length (y-axis)
    matrix = np.zeros((len(split_taxids), maxindexlength))
    seen = set()
    # On = white, black = off
    # Go through the taxid strings and turn on the pixels, starting from the rightmost taxid, if it has not yet been seen
    # Use the taxid_to_x_index dictionary to look up the x-index of the taxid
    color_dict = {10197:  "#54AB53", #Ctenophores   - 10197      -  #54AB53
                  6040 :  "#DCC0F3", #Sponges       - 6040       -  #DCC0F3
                  6073 :  "#387FB2", #Cnidarians    - 6073       -  #387FB2
                  10226:  "#C72480", #Placozoans    - 10226      -  #C72480
                  33317:  "#F4B93E", #Protostomes   - 33317      -  #F4B93E
                  33511:  "#78A6AF"  #Deuterostomes - 33511      -  #78A6AF
                  }
    row_color = []
    for i, taxidlist in enumerate(split_taxids):
        thiscolor = "#FFFFFF"
        for j in range(len(taxidlist)-1, -1, -1):
            if int(taxidlist[j]) in color_dict:
                thiscolor = color_dict[int(taxidlist[j])]
            if taxidlist[j] not in seen:
                seen.add(taxidlist[j])
                matrix[i][taxid_to_x_index[taxidlist[j]]] = 1
        row_color.append(thiscolor)
    # Convert the matrix to a Pillow Image
    # Instead of using white and black, each row gets a specific color if on, or black if off.
    # Use the color dict to look up the color
    #image = Image.fromarray((matrix * 255).astype('uint8'), 'L')
    image = Image.new('RGB', (maxindexlength, len(split_taxids)), color = (0, 0, 0))
    for i in range(len(split_taxids)):
        this_color = hex_to_rgb(row_color[i])
        for j in range(maxindexlength):
            if matrix[i][j] == 1:
                image.putpixel((j,i), this_color)
    return image
    # Save the image as PNG or JPEG
    image.save(outfile)

def _image_helper_get_ALG_columns(ALG_pres_abs_dataframe) -> list:
    """
    Takes the ALG presence/absence dataframe and returns a list of the ALG columns.
    The list will be output in the order in which the ALGs appear in the dataframe.
    """
    # get all of the columns that are a tuple
    tuple_columns = [x for x in ALG_pres_abs_dataframe.columns if isinstance(x, tuple)]
    # decompose the tuples into a set of unique entries
    unique_entries = set()
    for col in tuple_columns:
        for entry in col:
            unique_entries.add(entry)
    # change these to be the same order as the ALG_pres_abs_dataframe
    unique_entries = [x for x in ALG_pres_abs_dataframe.columns
                      if x in unique_entries]
    return unique_entries

def image_sp_matrix_to_presence_absence(ALG_pres_abs_dataframe, color_dict = None) -> Image:
    """
    Required Loadings:
      from PIL import Image
      import numpy as np

    Description:
      - Takes a df of ALG presence/absence, figures out what are the ALGS, and makes an image of it.
      - Returns the img to be handled later - probably to be added to another Image object.
      - If color_dict is not None, then it will make a key of the colors at the top of the image.
    """
    # get a list of the ALGs by parsing the column names. Returns them in the same order as the dataframe.
    unique_entries = _image_helper_get_ALG_columns(ALG_pres_abs_dataframe)
    # Just get a subdf of the columns that have the unique entries.
    # These are the ALGs.
    filtdf = ALG_pres_abs_dataframe[unique_entries]
    matrix = filtdf.values
    # make a PIL image in the same coordinates as the filtdf
    image = Image.new('RGB', (len(filtdf.columns), len(filtdf)),
                      color = (0, 0, 0))
    # iterate through the matrix i,j, and
    # if the value is 0, black, any other value, white
    for i in range(len(filtdf)):
        for j in range(len(filtdf.columns)):
            if matrix[i][j] == 0:
                #image.putpixel((j,i), (0,0,0))
                pass
            else:
                color = (255,255,255)
                #if color_dict is not None:
                #    color = hex_to_rgb(color_dict[filtdf.columns[j]])
                image.putpixel((j,i), color)

    if color_dict is not None:
        # make a ten-pixel high colorbar on the top.
        colorbar = Image.new('RGB', (len(filtdf.columns), 10), color = (0, 0, 0))
        # The colorbar encodes the ALG colocalization pairs.
        # The top 8 pixels are colored, and the bottom two pixels are black
        for i in range(len(filtdf.columns)):
            color = hex_to_rgb(color_dict[filtdf.columns[i]])
            for j in range(8):
                colorbar.putpixel((i,j), color)
        # concatenate the colorbar to the top of the image
        image = Image.fromarray(np.concatenate((np.array(colorbar), np.array(image)), axis=0))

    return image

def dict_BCnSALG_to_color() -> dict:
    """
    Makes a dictionary of the BCnS ALG strings as keys, the colors as the values.
    """
    return {
          "Qb": "#C72480",  #  12
          "Qc": "#DCC0F3",  #  14
          "C2": "#387FB2",  #  18
          "Qd": "#94C47F",  #  22
           "R": "#F4B93E",  #  24
          "Qa": "#78A6AF",  #  30
          "A2": "#8B4E67",  #  41
          "B3": "#FA9A26",  #  46
          "O2": "#AB5BA8",  #  46
          "Eb": "#B76BED",  #  47
         "A1b": "#C33D53",  #  51
          "J1": "#54AB53",  #  54
          "O1": "#FBD76C",  #  55
          "J2": "#E64657",  #  66
           "P": "#C33E51",  #  78
          "B2": "#1F779A",  #  86
           "I": "#3425FB",  #  90
          "B1": "#2F54E3",  #  95
           "M": "#A45530",  # 102
           "L": "#7DC29F",  # 104
           "N": "#D8BE3C",  # 107
          "Ea": "#AB7E26",  # 115
           "K": "#170B88",  # 119
           "H": "#F04C08",  # 135
           "G": "#E97B4A",  # 138
          "C1": "#B07DF4",  # 142
           "F": "#9B6870",  # 145
           "D": "#47957F",  # 172
         "A1a": "#4DB5E3"}  # 207


def image_colocalization_matrix(ALG_pres_abs_dataframe, color_dict = None,
                                clustering = True, missing_data_color = "#990000") -> Image:
    """
    Required Loadings:
      from PIL import Image
      import numpy as np
      from sklearn.metrics.pairwise import cosine_similarity
      from scipy.cluster import hierarchy

    Description:
      - Plots the colocalization matrix from the ALG_fusions dataframe.
      - If we choose to sort the fusions, then sort only by similarity on the x-axis. This is a type of clustering.
        We should not sort on the y-axis, because the dataframe will already have been sorted by the NCBI TaxID string.

      If there is an input dictionary of colors, makes a two-pixel high colorbar on the top.
      The colorbar encodes the ALG colocalization pairs.
    """
    # get all of the ALGs in the dataframe
    unique_entries = _image_helper_get_ALG_columns(ALG_pres_abs_dataframe)
    # Just get a subdf of the columns that have the ALG entries.
    ALG_df = ALG_pres_abs_dataframe[unique_entries]

    # make a dict of rows, and the ALGs that have a 1 in that row
    ALG_presence = {}
    for i, row in ALG_df.iterrows():
        if i not in ALG_presence:
            ALG_presence[i] = {x: int(row[x]) for x in ALG_df.columns}
    # make a subdf of the ALG combinations. These are all of the columns that are tuples.
    coloc_df = ALG_pres_abs_dataframe[[x for x in ALG_pres_abs_dataframe.columns if isinstance(x, tuple)]]

    if clustering == True:
        # CLUSTER THE COLUMNS BASED ON SIMILARITY BETWEEN SPECIES
        # Calculate similarity matrix (using correlation here)
        #similarity_matrix = coloc_df.corr()
        from sklearn.metrics.pairwise import cosine_similarity
        similarity_matrix = cosine_similarity(coloc_df.transpose())
        # Perform hierarchical clustering
        linkage_matrix = hierarchy.linkage(similarity_matrix, method='average')
        clustered_order = hierarchy.leaves_list(linkage_matrix)
        ## Dendrogram for visualization
        #dendrogram = hierarchy.dendrogram(linkage_matrix, labels=similarity_matrix.columns, orientation='top')
        #plt.show()
        # Get the clustered column order
        clustered_order = hierarchy.leaves_list(linkage_matrix)
        # Rearrange columns based on clustering
        coloc_df = coloc_df.iloc[:, clustered_order]

    # make an image that is the same size as the coloc_df
    image = Image.new('RGB', (len(coloc_df.columns), len(coloc_df)), color = (0, 0, 0))
    # Go through the pixels and turn them to white if the value is 1, black if 0.
    # Special case for 0 - if the value is 0, then we need to check if the ALG is present in the row.
    #   if either of the ALGs are not present in that row, then turn it the missing_data_color color
    # Lookup the actual index of the row from its index on the left. In otherwords, the index != the row number.
    for i, row in coloc_df.iterrows():
        # get the row index
        ri = ALG_df.index.get_loc(i)
        for j, col in enumerate(coloc_df.columns):
            if row[col] == 0:
                # check if the ALG is present in the row
                if ALG_presence[i][col[0]] == 0 or ALG_presence[i][col[1]] == 0:
                    image.putpixel((j,ri), hex_to_rgb(missing_data_color))
                    # change this value in the dataframe to be -1, to represent no potential to find something there
                    coloc_df.at[i, col] = -1
                else:
                    image.putpixel((j,ri), (0,0,0))
            else:
                image.putpixel((j,ri), (255,255,255))

    # make a colorbar if we have a color_dict
    if color_dict is not None:
        # make a new image that is 10 pixels high and the same width as the image
        # The top 8 pixels are the colors and the bottom two rows are black
        colorbar = Image.new('RGB', (len(coloc_df.columns), 10), color = (0, 0, 0))
        # go through the columns and add the color to the colorbar
        for i, col in enumerate(coloc_df.columns):
            # get the color
            colortop = color_dict[col[0]]
            colorbot = color_dict[col[1]]
            # add the color to the colorbar
            for j in range(4):
                colorbar.putpixel((i,j), hex_to_rgb(colortop))
            for j in range(4,8):
                colorbar.putpixel((i,j), hex_to_rgb(colorbot))
        # concatenate the colorbar to the image
        image = Image.fromarray(np.concatenate((np.array(colorbar), np.array(image)), axis=0))

    return image

def image_concatenate_vertically(image_objects):
    """
    DEPRECATED. Do we need this?
    Pastes images together vertically. If input is [1,2,3]
    1
    v
    2
    v
    3
    Everything will be left-aligned.
    """
    # Calculate the total height and maximum width among all images
    total_height = sum(img.height for img in image_objects)
    max_width = max(img.width for img in image_objects)

    # Initialize the combined image
    combined_image = Image.new('RGB', (max_width, total_height), (255, 255, 255))  # White background

    # Paste each image at the correct position
    current_y = 0
    for img in image_objects:
        # Calculate the x-coordinate to center the image horizontally
        x_offset = (max_width - img.width) // 2
        combined_image.paste(img, (x_offset, current_y))
        current_y += img.height
    return combined_image

def image_concatenate_horizontally(image_objects, valign="bottom"):
    """
    Takes a bunch of images and concatenates them horizontally 1>2>3.
    They will be aligned at the bottom.
    """
    # Calculate the total width and maximum height among all images
    total_width = sum(img.width for img in image_objects)
    max_height = max(img.height for img in image_objects)

    # Initialize the combined image
    combined_image = Image.new('RGB', (total_width, max_height), (0, 0, 0))  # Black background

    # Paste each image at the correct position
    current_x = 0
    for img in image_objects:
        if valign == "center":
            # The y position of everything is vertically centered
            y_offset = (max_height - img.height) // 2
        elif valign == "bottom":
            # all of the images are aligned along the bottom
            y_offset = max_height - img.height
        combined_image.paste(img, (current_x, y_offset))
        current_x += img.width

    return combined_image

def image_vertical_barrier(width, height, color = "#F0F93E") -> Image:
    """
    returns an image that is a vertical line
    of the specified width, height, and color
    """
    RGB_color = hex_to_rgb(color)
    return Image.new('RGB', (width, height),
                      color = RGB_color)

def standard_plot_out(perspchrom, outprefix)->None:
    """
    Makes a pixel-wise plot of a dataframe.
    Currently has a "phylogenetic" bar on the left side.
    Then there is an accounting of which ALGs are present in each species.
    Then there is a plot of the ALG colocalization matrix.

    At the top there is a colorbar that encodes the ALG colocalization pairs.
    """
    # convert matrix to a species image
    # make an image of a pseudo-phylogenetic tree
    tree_image    = image_sp_matrix_to_lineage(perspchrom["taxidstring"])

    bar1 = image_vertical_barrier(2, len(perspchrom), color = "#000000")
    bar2 = image_vertical_barrier(1, len(perspchrom), color = "#F0F93E")
    bar3 = image_vertical_barrier(2, len(perspchrom), color = "#000000")

    # make an image of the ALG presence/absence matrix
    presabs_image = image_sp_matrix_to_presence_absence(perspchrom, color_dict = dict_BCnSALG_to_color())

    # gap bars
    bar4 = image_vertical_barrier(2, len(perspchrom), color = "#000000")
    bar5 = image_vertical_barrier(1, len(perspchrom), color = "#F0F93E")
    bar6 = image_vertical_barrier(2, len(perspchrom), color = "#000000")

    # make an image of the colocalization matrix
    coloc_image = image_colocalization_matrix(perspchrom, clustering = True, color_dict = dict_BCnSALG_to_color())
    coloc_image_unclust = image_colocalization_matrix(perspchrom, clustering = False, color_dict = dict_BCnSALG_to_color())

    # concatenate all the images
    concatenate_these_images = [tree_image, bar1, bar2, bar3, presabs_image,
                                bar4, bar5, bar6, coloc_image]
    composite_image = image_concatenate_horizontally(concatenate_these_images)
    composite_image.save("{}_composite_image.png".format(outprefix))

    # make the same plot, but not sorted
    concatenate_these_images = [tree_image, bar1, bar2, bar3, presabs_image,
                                bar4, bar5, bar6, coloc_image_unclust]
    composite_image = image_concatenate_horizontally(concatenate_these_images)
    composite_image.save("{}_composite_image_unclustered.png".format(outprefix))

def plot_missing_vs_colocalized(perspchrom, fileprefix):
    """
    This makes a dotplot of the number of missing combinations vs the number of colocalizations.

    To find the number of missing combinations, we find the non-tuple columns with zeros, then make
     combinations of them.
    To find the number of colocalizations, we find the number of 1s in the tuple columns.
    """
    # Find all of the columns that are not tuples or tuple-like
    non_tuple_columns = [x for x in perspchrom.columns if not isinstance(x, tuple)]
    # remove the columns that are ["species", "taxid", "taxidstring"]
    remove_these = ["species", "taxid", "taxidstring"]
    non_tuple_columns = [x for x in non_tuple_columns if x not in remove_these]
    # filter the df
    entries = []
    # go through each row of the pandas dataframe perspchrom
    for i, row in perspchrom.iterrows():
        # get missing non-tuple columns
        missing_ALGs = [x for x in non_tuple_columns if row[x] == 0]
        # get the number of n choose two combinations of the missing ALGs
        missing_combinations = len(list(itertools.combinations(missing_ALGs, 2)))
        # get the number of colocalizations
        number_colocalizations = len([x for x in row.index if isinstance(x, tuple) and row[x] == 1])
        entries.append({"missing_combinations": missing_combinations,
                        "number_colocalizations": number_colocalizations})
    # Make a scatterplot of the missing combinations vs the number of colocalizations. Use matplotlib only.
    # The x-axis is the number of missing combinations
    # The y-axis is the number of colocalizations
    df = pd.DataFrame(entries)
    plt.scatter(df["missing_combinations"], df["number_colocalizations"])
    plt.xlabel("Number of missing combinations")
    plt.ylabel("Number of colocalizations")
    filename = "{}_missing_vs_colocalized.pdf".format(fileprefix)
    plt.show()

def missing_present_ALGs(df, min_for_missing = 0.8):
    """
    This takes in a dataframe and returns lists of which ALGs are missing and present.
    The min_for_missing is the minimum fraction of species in this clade that must have the ALG missing
      for this ALG to be considered missing.

    Returns the missing_ALGs and the present_ALGs as two separate lists.
    """
    # first get all of the ALGs from this dataframe by finding the non-tuple columns that are not in ["species", "taxid", "taxidstring"]
    remove_these_columns = ["species", "taxid", "taxidstring", "changestrings"]
    ALG_columns = [x for x in df.columns if not isinstance(x, tuple) and x not in remove_these_columns]
    # Now we find the ALGs that have the min_for_missing fraction of species missing
    missing_ALGs = []
    for thisALG in ALG_columns:
        # get the fraction of species that have this ALG missing
        fraction_missing = len(df[df[thisALG] == 0]) / len(df)
        if fraction_missing >= min_for_missing:
            missing_ALGs.append(thisALG)
    present_ALGs = [x for x in ALG_columns if x not in missing_ALGs]
    return missing_ALGs, present_ALGs

def separate_ALG_pairs(df, min_for_noncolocalized = 0.5):
    """
    This takes in a dataframe and returns a list of the ALGs
     that are confirmed to be separated in this clade.

    We need to identify the ALGs that are separated in this clade.
    We do not try to find ALGs that are fused in this clade, because that
      doesn't allow us to polarize against what we know about the clade in question:
      the ALGs that have been fused in this clade.
    Additionally, we do not try to find ALGs that are fused, because we may have
      lost the ability to detect fusion in this clade depending on the genome quality
      or if the ALGs in question have dispersed.

    By using the logic above, we push back the fusion event to the earliest possible
      node, rather than the latest. This is the same preference that we have given
      to detecting the node on which the ALGs are lost.
    """
    # First get all of the ALGs from this dataframe by finding the non-tuple
    #   columns that are not in ["species", "taxid", "taxidstring"]
    remove_these_columns = ["species", "taxid", "taxidstring", "changestrings"]
    ALG_columns  = [x for x in df.columns if not isinstance(x, tuple) and x not in remove_these_columns]
    pair_columns = [x for x in df.columns if isinstance(x, tuple) and x not in remove_these_columns]

    # We must count the number of detectable pairs,
    #  and we must also count the number of times the pairs are NOT colocalized.
    # Definition: a pair is detectable if both ALGs are present in that species. (1 values in the ALG columns)
    # Definition: a pair is colocalized if both ALGs are present in the same chromosome in that species.
    pair_detectability_count = {x: 0 for x in pair_columns}
    pair_nonlocalized_count  = {x: 0 for x in pair_columns}
    # go through the rows of the dataframe. For the detectability count, go through ALG column pairs to see if they are both 1.
    # For the detected count, go through the pair columns to see if they are 1.
    for i, row in df.iterrows():
        for pair in pair_columns:
            # if the pair is detectable, then add one to the detectability count
            if row[pair[0]] == 1 and row[pair[1]] == 1:
                pair_detectability_count[pair] += 1
                # In this case, as we only see if they are nonlocalized if they were both detected in the first place.
                # if the pair is NOT detected, then add one to the detected count
                if row[pair] == 0:
                    pair_nonlocalized_count[pair] += 1
    # In this dictionary, a value of 0 means that the pair was never detected.
    #  The reason it was never detected could be because it wasn't detectable (The ALGs have dispersed)
    #  Or the ALGs were detectable, but they were colocalized
    pair_nonlocalized_detectability_ratio = {}
    for thispair in pair_columns:
        if pair_detectability_count[thispair] > 0:
            pair_nonlocalized_detectability_ratio[thispair] = pair_nonlocalized_count[thispair] / pair_detectability_count[thispair]
        else:
            pair_nonlocalized_detectability_ratio[thispair] = 0
    # return the pairs that are not colocalized in at least min_for_noncolocalized fraction of the species
    return [x for x in pair_nonlocalized_detectability_ratio
            if pair_nonlocalized_detectability_ratio[x] >= min_for_noncolocalized]

def main():
    # parse the args
    args = parse_args()

    # get the rbh files in the directory
    rbh_files = list(sorted([os.path.join(args.directory, f)
                 for f in os.listdir(args.directory)
                 if f.endswith('.rbh')], reverse = True))

    print("The rbh file set length before filtering is {}".format(len(rbh_files)))
    # just pick the first rbh file to get the one species with each taxid. All of our files will have the NCBI taxid in the file name.
    seen_taxids = set()
    keep_these = []
    # save the list of selected genomes to a file:
    outhandle  = open("selected_genomes.txt", 'w')
    for thisfile in rbh_files:
        # get just the filename from the whole path
        thisfilename = os.path.basename(thisfile)
        # get the taxid. When we split on '-', it will be the 1st element, zero-based indexing.
        taxid = int(thisfilename.split('-')[1])
        if taxid not in seen_taxids:
            keep_these.append(thisfile)
            print("{}".format(thisfilename), file=outhandle)
        seen_taxids.add(taxid)
    outhandle.close()

    rbh_files = keep_these
    print("The rbh file set length after filtering is {}".format(len(rbh_files)))

    # first we need to read through all of the rbh files to get all of the possible ALGs to plot
    #ALGdf = parse_rbh_file(args.ALG_rbh, args.ALGname)
    ALGdf = parse_rbh_file(args.ALG_rbh)

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
            # filter out the rows that are not significantly correlated, using minsig
            group = group[group["whole_FET"] <= args.minsig]
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

    # *********************************************************************************************************************
    # This block makes a symmetrical table of fusion events that were significantly detected for all the species
    # *********************************************************************************************************************
    # make a df of counts where two columns had values of 1 on the same row
    # counts is a dict of all the chromosomes
    ALG_list = list(sorted(ALGdf["ALGname"].tolist()))
    counts = {}
    for i in range(len(ALG_list)):
        for j in range(i, len(ALG_list)):
            sorted_list = sorted([ALG_list[i], ALG_list[j]])
            thistuple = (sorted_list[0], sorted_list[1])
            thattuple = (sorted_list[1], sorted_list[0])
            counts[thistuple] = 0
            counts[thattuple] = 0

    # go through the rows of the df and count the number of times that two ALGs are colocalized
    for i, row in df.iterrows():
        # get the ALGs that are colocalized on this chromosome
        ALGs = [x for x in row.index if x in ALG_list and row[x] == 1]
        # get all of the combinations of these ALGs
        for i in range(len(ALGs)):
            for j in range(i, len(ALGs)):
                sorted_list = sorted([ALGs[i], ALGs[j]])
                thistuple = (sorted_list[0], sorted_list[1])
                thattuple = (sorted_list[1], sorted_list[0])
                # if we're looking at the same ALG, don't do anything
                if i != j:
                    counts[thistuple] += 1
                    counts[thattuple] += 1

    # convert counts to a dataframe where tuple 0 is a col called ALG1, tuple 1 is a col called ALG2, and the value is the number of times they are colocalized
    countsdf = pd.DataFrame.from_dict(counts, orient='index').reset_index()
    # split the "index" column into two columns
    countsdf[["ALG1", "ALG2"]] = pd.DataFrame(countsdf["index"].tolist(), index=countsdf.index)
    # drop the index column
    countsdf = countsdf.drop(columns=["index"])
    # rename the columns
    countsdf = countsdf.rename(columns={0: "count"})
    # unmelt the dataframe so that it is a square matrix. All of the columns will be ints
    countsdf = countsdf.pivot(index="ALG1", columns="ALG2", values="count")
    # fill the NaNs to -1
    # raise an error if there are any NAs. We have already made the matrix symmetric
    if countsdf.isnull().values.any():
        raise ValueError("There are NAs in the countsdf")
    # convert all the columns to ints
    countsdf = countsdf.astype(int)
    # flip the columns and rows
    countsdf = countsdf.transpose()
    # save this to a file that is called "all_ALG_fusion_table.tsv"
    countsdf.to_csv("all_ALG_fusion_table.tsv", sep='\t')

    # *********************************************************************************************************************
    #  This block prints out a species-specific table of fusion events
    #   - each row is one species
    #   - Just a reminder, the cutoff we use to detect is FET <= 0.005, which is the Bonferroni-corrected p-value
    # *********************************************************************************************************************
    # first make a blank data structure to hold all of the possible ALG combinations
    #  Don't inlcude the ALG with itself, like before.
    #  Instead, add an entry of the ALG by itself to indicate whether that ALG exists at all in this species.
    #  This will be useful for tracking whether that ALG just simply disappears.
    ALG_list = list(sorted(ALGdf["ALGname"].tolist()))
    counts = {(x): 0 for x in ALG_list}
    for i in range(len(ALG_list)-1):
        for j in range(i+1, len(ALG_list)):
            sorted_list = sorted([ALG_list[i], ALG_list[j]])
            thistuple = (sorted_list[0], sorted_list[1])
            counts[thistuple] = 0
    # go through the rows of the df and count the number of times that two ALGs are colocalized
    # the keys here are just the species names as they are in the file. Right now the taxid is encoded at the end.
    sp_df = {}
    for i, row in df.iterrows():
        thissp = i.split(".|.")[0]
        if thissp not in sp_df:
            sp_df[thissp] = counts.copy()
        # get the ALGs that are colocalized on this chromosome
        ALGs = [x for x in row.index if x in ALG_list and row[x] == 1]
        # Mark whether or not we have seen this ALG in this species.
        # We use a separate for loop here due to the fact that in the next for loop it won't reach the end
        for i in range(len(ALGs)):
            # we just mark whether or not we have seen this ALG in this species
            if sp_df[thissp][(ALGs[i])] == 0:
                sp_df[thissp][(ALGs[i])] = 1
        # get all of the combinations of these ALGs. Add a 1 if we have seen this combo in this species
        for i in range(len(ALGs)-1):
            for j in range(i+1, len(ALGs)):
                sorted_list = sorted([ALGs[i], ALGs[j]])
                thistuple = (sorted_list[0], sorted_list[1])
                # mark that this combination was seen in this species
                if sp_df[thissp][thistuple] == 0:
                    sp_df[thissp][thistuple] = 1

    from ete3 import NCBITaxa,Tree

    ncbi = NCBITaxa()
    # Now that we have looked at every scaffold in every species, add a few missing fields
    for thissp in sp_df:
        # Just get the integer at the end of this string by removing all the alpha characters
        # Note 202401 - The place we put the taxid has changed since last time. It used to be at the end of the line.
        thistaxid = int(thissp.split("-")[1])
        # Raise an error if there is something in the string that is not 0-9
        if not re.match(r"^[0-9]*$", str(thistaxid)):
            raise ValueError("There is a non-numeric character in the taxid string")
        # get the NCBI taxid string
        thistaxidstring = ";".join([str(x) for x in ncbi.get_lineage(thistaxid)])
        sp_df[thissp]["taxid"] = thistaxid
        sp_df[thissp]["taxidstring"] = thistaxidstring

    # now make a pandas dataframe
    perspchrom = pd.DataFrame.from_dict(sp_df, orient='index').reset_index()
    # move the taxid and taxidstring columns to the front
    perspchrom.insert(0, "taxidstring", perspchrom.pop("taxidstring"))
    perspchrom.insert(0, "taxid", perspchrom.pop("taxid"))
    perspchrom.insert(0, "index", perspchrom.pop("index"))
    # rename the first column to the species string
    perspchrom = perspchrom.rename(columns={"index": "species"})
    # sort by the taxidstring
    perspchrom = perspchrom.sort_values(by="taxidstring", ascending=True)
    # THIS FILE IS SAVED LATER

    # now save it in a way that can be used by R
    # get all of the columns that are not the species, taxid, or taxidstring
    drop_these_columns = ["species", "taxid", "taxidstring"]
    if "changestrings" in perspchrom.columns:
        drop_these_columns.append("changestrings")
    perspchrom_R = perspchrom.drop(columns=drop_these_columns)
    # save it in a format that can be used by R. meaning that we surround the indices with quotation marks and sep with commas
    perspchrom_R.to_csv("per_species_ALG_presence_fusions_R.tsv", sep=',', index=True, header=True)

    ## make a scatterplot of the # of lost ALG combos vs number of fused ALG combos. Color by taxid
    #plot_missing_vs_colocalized(perspchrom, "per_species_ALG_presence_colocalizations")

    # make figures of the per-species plots
    standard_plot_out(perspchrom, "perSpecies")

    ### This block makes a simple binary image of the fusion matrix
    # Create a binary matrix (0s and 1s)
    matrix = perspchrom_R.values
    # Convert the matrix to a Pillow Image
    image = Image.fromarray((matrix * 255).astype('uint8'), 'L')
    # Save the image as PNG or JPEG
    image.save("output.png")
    # print the tree to a .tre file
    tree = ncbi.get_topology([int(x) for x in perspchrom["taxid"].tolist()])

    # *********************************************************************************************************************
    #     ┏┓┓ ┏┓  ┳┓┳┏┓┏┓┏┓┳┓┏┓┳┏┓┳┓  ┏┓┳┓┳┓  ┏┓┳┳┏┓┳┏┓┳┓
    #     ┣┫┃ ┃┓  ┃┃┃┗┓┃┛┣ ┣┫┗┓┃┃┃┃┃  ┣┫┃┃┃┃  ┣ ┃┃┗┓┃┃┃┃┃
    #     ┛┗┗┛┗┛  ┻┛┻┗┛┃ ┗┛┛┗┗┛┻┗┛┛┗  ┛┗┛┗┻┛  ┻ ┗┛┗┛┻┗┛┛┗
    # *********************************************************************************************************************

    # Determine per-species number of changes
    # When I mark "changes on a node", I mean that there was a change on the branch leading up to this node.
    #  This logic applies to the species-level, but also to internal nodes.
    #  For each node, there is a possibility to lose ALGs, or to gain fusions.
    #  For the fusions, we must should consider that AxB fusing with CxD is just one event, and we should not count AxC and AxD and BxC and BxD as separate events.
    #   To satisfy this we can build a graph and count connected components.
    # Do this species-wise by iterative through the perspchrom df

    # We need cutoffs to identify if something is missing or not
    # The cutoffs for missing ALGs and missing ALG combos are different
    # For an ALG to be missing, it is likely missing in the entire clade. To give some wiggle room,
    #  we will ask if this ALG is missing in 80% of the comparison clade to determine whether it is absent or not.
    #  ----- ALG PRESENCE ABSENCE -----
    #  For each species, when we begin we must mark which ALGs are conserved in this species.
    #  These can never be lost, so these will never appear in the changeString.
    #  If the ALG is missing in 80% of the comparison clade,
    #      and the ALG is absent in the sample, this means that the loss was earlier.
    #          We don't do anything in this case.
    #      but the ALG is present in the sample, this likely means that the ALG was
    #          lost somewhere else in this clade.
    #          We don't do anything in this case.
    #  Otherwise (if the ALG is not missing in 80% of the comparison clade),
    #      and if the ALG is absent in this sample, then this suggests that the ALG was
    #          dispersed on this branch. We should mark this as a loss, and record the loss.
    #      and if the ALG is is present in the sample, then this suggests that the ALG is conserved in both clades.
    #          We don't need to do anything in this case.
    perspchrom = perspchrom.reset_index(drop=True)
    # We need to know the node for which the ALGs were inferred.
    #  Once we hit this node in the NCBI taxonomy string, we place the remaining "lost ALGs" on the branch between the previous node and this node.
    ALG_node = "1;131567;2759;33154;33208"
    # First we need to track which ALGs are present in this dataset, and what combinations are tracked in the df.
    remove_these = ["species", "taxid", "taxidstring", "changestrings"]
    ALG_columns = [x for x in perspchrom.columns if not isinstance(x, tuple) and not x in remove_these]
    ALG_combos  = [x for x in perspchrom.columns if isinstance(x, tuple) and not x in remove_these]
    perspchrom["changestrings"] = ""

    # This is a magic number used to determine whether an ALG is missing or not.
    #  If an ALG is missing in 80%, or 0.8, of the clade in question, then we consider it missing.
    min_for_missing = 0.8
    # If the other clade's ALGs are detectable, and are not colocalized in at least 50% of the clade,
    #  then we consider that this ALG pair is not colocalized in this clade.
    min_for_noncolocalized = 0.5

    for i, row in perspchrom.iterrows():
        # ----- ALG PRESENCE ABSENCE -----
        # For this species, we keep track of the changes on each branch with a structured list.
        # The format is [taxid, "(gain|loss)", "taxid", "(gain|loss)", ...]
        #  We add new entries as we go through the tree to the root.
        #  At the end, we flip the order and make a parsable string
        changeString = []
        # we have access to ["species", "taxid", "taxidstring"]
        thistaxidstring = row["taxidstring"]
        ALGsAlreadyFused = []
        # note the ALGs that have a value of 1 here. These are conserved in this species and will never be in changeString
        ALGsConservedInThisSample    = [x for x in ALG_columns if row[x] == 1]
        # Note the things specifically missing in this species. These will be added to the changeString at some point
        #  This list will not be changed at all during execution of this loop. Only the next ALGsMissingInThisSampleAccountedFor will be changed.
        ALGsMissingInThisSample      = [x for x in ALG_columns if row[x] == 0]
        # We will add to this list as we go through the tree. This will be the list of ALGs that were lost on various branches.
        ALGsMissingInThisSampleAccountedFor = [] # At the end, all of the missing ALGs must appear somewhere in ALGsMissingInThisSample

        # ----- ALG COLOCALIZATION -----
        # For each species, we similarly need to keep track of which ALG combos are gained on each branch.
        ColocalizedPairsInThisSample             = [x for x in ALG_combos if row[x] == 1]
        ColocalizedPairsInThisSampleAccountedFor = []

        #print("Starting taxidstring: {}".format(thistaxidstring))
        #print("species", row["species"])
        #print("ALGsMissingInThisSample: {}".format(ALGsMissingInThisSample))
        #print("ALGsConservedInThisSample: {}".format(ALGsConservedInThisSample))
        #print("ColocalizedPairsInThisSample: {}".format(ColocalizedPairsInThisSample))

        while thistaxidstring.count(";") > 0:
            # we need to know where we came from. Keep the previous taxidstring
            prevtaxidstring = thistaxidstring
            thisdf = perspchrom[perspchrom["taxidstring"] == thistaxidstring]

            # remove the last taxid from the string
            thistaxidstring = ";".join(thistaxidstring.split(";")[:-1])

            # these will be used to keep track of what is gained or lost on this branch
            ALG_colocalizations_on_this_branch = []
            ALG_losses_on_this_branch = []
            #print("New shorter string: {}".format(thistaxidstring))

            # Now we perform the logic of estimating whether an ALG was lost or gained on this branch,
            #  and whether an ALG combo was gained on this branch.
            if thistaxidstring == ALG_node:
                # ----- ALG PRESENCE ABSENCE -----
                # In this case we are at the node for which the ALGs were inferred.
                # For example, for the BCnS ALGs, the node is Metazoa: "1;131567;2759;33154;33208"
                #  We need to place the remaining ALGs on the branch between the previous node and this node.
                ALG_losses_on_this_branch = [x for x in ALGsMissingInThisSample if x not in ALGsMissingInThisSampleAccountedFor]
                # add all the new ALG losses to the ALGsMissingInThisSampleAccountedFor
                ALGsMissingInThisSampleAccountedFor.extend(ALG_losses_on_this_branch)

                # ----- ALG COLOCALIZATION EVENTS -----
                ALG_colocalizations_on_this_branch = [x for x in ColocalizedPairsInThisSample if x not in ColocalizedPairsInThisSampleAccountedFor]
                # add all the new ALG colocalizations to the ColocalizedPairsInThisSampleAccountedFor
                ColocalizedPairsInThisSampleAccountedFor.extend(ALG_colocalizations_on_this_branch)
            else:
                # Because we're not at the last node we need to do some more work
                # Get the subdf of the perspchrom df for rows with a taxidstring that contains the new thistaxidstring, but not the prevtaxidstring
                #  Excluding the prevtaxidstring is important, because it compares this clade only to sister clades at the same level,
                #  and removes the possibility of paraphyletic comparisons.
                subdf = perspchrom[perspchrom["taxidstring"].str.contains(thistaxidstring) & ~perspchrom["taxidstring"].str.contains(prevtaxidstring)]
                if len(subdf) == 0:
                    # If there is nothing to compare to in this clade, then we don't need to do anything and can proceed to printing the changeString
                    # In the future there could be an option to put ambiguous changes of gains or losses here, but for now
                    #  we are just going to put the changes on the oldest nodes that we can find.
                    pass
                else:
                    # ----- ALG PRESENCE ABSENCE -----
                    # If there are some species to compare here then we apply the logic detailed above in the section ALG PRESENCE ABSENCE
                    missingALGsInThisClade, presentALGsInThisClade = missing_present_ALGs(subdf, min_for_missing = min_for_missing)
                    # The only time we mark a loss of an ALG is if it is not missing in 80% of the comparison clade.
                    #  These are the things that are in the presentALGsInThisClade list.
                    #  If there is something in presentALGsInThisClade, but in ALGsMissingInThisSample, then we mark it as a loss on this branch.
                    ALGsLostOnThisBranch = [x for x in presentALGsInThisClade
                                            if x in ALGsMissingInThisSample
                                            and x not in ALGsMissingInThisSampleAccountedFor]
                    #print("Present ALGs: {}".format(presentALGsInThisClade))
                    #print("ALGs accounted for so far: {}".format(ALGsMissingInThisSampleAccountedFor))
                    #print("ALGs lost on this branch: {}", ALGsLostOnThisBranch)
                    # note that we have to print these.
                    ALG_losses_on_this_branch.extend(ALGsLostOnThisBranch)
                    # then update the dataframe saying that we know that these ALGs were lost on this branch
                    ALGsMissingInThisSampleAccountedFor.extend(ALGsLostOnThisBranch)

                    # ----- ALG COLOCALIZATION EVENTS -----
                    notLocalizedInThisClade = separate_ALG_pairs(subdf, min_for_noncolocalized = min_for_noncolocalized)
                    pairsColocalizedOnThisBranch = [x for x in notLocalizedInThisClade
                                                    if x in ColocalizedPairsInThisSample
                                                    and x not in ColocalizedPairsInThisSampleAccountedFor]
                    #print("notLocalizedInThisClade: {}".format(notLocalizedInThisClade))
                    #print("pairsColocalizedOnThisBranch: {}".format(pairsColocalizedOnThisBranch))
                    # Mark which colocalizations we have accounted for already
                    ColocalizedPairsInThisSampleAccountedFor.extend(pairsColocalizedOnThisBranch)
                    # Mark which colocalizations we have gained on this branch
                    ALG_colocalizations_on_this_branch.extend(pairsColocalizedOnThisBranch)
                #print()
            thisChange = "({}|{})".format(sorted(ALG_colocalizations_on_this_branch), sorted(ALG_losses_on_this_branch))
            changeString.append(prevtaxidstring.split(";")[-1])
            changeString.append(thisChange)
        # Now we should check that ALGsMissingInThisSampleAccountedFor is the same as ALGsMissingInThisSample.
        # This means that we have found all of the missing ALGs at some point in the tree.
        # If this is not the case, then we need to raise an error.
        if sorted(ALGsMissingInThisSampleAccountedFor) != sorted(ALGsMissingInThisSample):
            raise IOError("There is a discrepancy between the ALGsMissingInThisSampleAccountedFor and ALGsMissingInThisSample. Write more debugging code to figure out what the issue is, because I haven't worked on this yet.")
        # we should do the same for the colocalizations
        if sorted(ColocalizedPairsInThisSampleAccountedFor) != sorted(ColocalizedPairsInThisSample):
            raise IOError("There is a discrepancy between the ColocalizedPairsInThisSampleAccountedFor and ColocalizedPairsInThisSample. Write more debugging code to figure out what the issue is, because I haven't worked on this yet.")
        # We have stepped out of the for loop, now we add the last taxidstring to the changeString
        changeString.append(thistaxidstring)
        # flip the changeString and make it a parsable string
        changeString = "-".join(changeString[::-1])
        # The final string will look something like this: 1-([]|[])-131567-([]|[])-2759-([]|[])-33154-([]|[])-33208-([]|['A1b', 'A2', 'B2', 'B3', 'C1', 'C2', 'D', 'Ea', 'F', 'G', 'H', 'I', 'J1', 'K', 'L', 'M', 'N', 'O1', 'P', 'Qa', 'Qb', 'Qc', 'Qd', 'R'])-6040-([]|[])-60882-([]|[])-60883-([]|[])-60884-([]|[])-472148-([]|[])-111877-([]|[])-111878
        # We add the final string to the perspchrom df
        perspchrom.at[i, "changestrings"] = changeString
        # print a progress bar on the same line to tell the user how much longer we have to go
        print("\r{}/{}          ".format(i+1, len(perspchrom)), end="")

    # save the file to a tsv
    perspchrom.to_csv("per_species_ALG_presence_fusions.tsv", sep='\t', index=False)

    # ---------------------------------------------------------------------------------------------
    #  Move onto the per-node analysis
    # ---------------------------------------------------------------------------------------------
    # Get the labels. They have to be in the order that the leaves are returned
    leaves = [int(str(x).split("-")[-1]) for x in tree.get_leaves()]
    # Make a dict with the taxid and species cols, then make a label from the lookup with leaves
    lookup = dict(zip(perspchrom["taxid"], perspchrom["species"]))
    labels = [lookup[x] for x in leaves]
    for leaf, label in zip(tree.get_leaves(), labels):
        leaf.name = f"{label}_{leaf.name}"
    tree.write(format=1, outfile="species_tree.tre")

    # Now annotate all of the nodes.
    # Yes, this loops through the table again, but I don't have a more elegant solution
    # Right now this doesn't actually annotate any nodes, it just makes a dictionary of the annotations
    node_annotations = {}
    # Make a dict of annotations for each node
    for i, row in perspchrom.iterrows():
        spstring = row["species"]
        taxidstring = [int(x) for x in row["taxidstring"].split(";")]
        for thisid in  taxidstring:
            if int(thisid) not in node_annotations:
                node_annotations[int(thisid)] = set()
            node_annotations[int(thisid)].add(spstring)

    # we now have a dictionary with which species belong in each node
    # iterate through all of the nodes of the tree, recursively. ACTUALLY IT ISN'T DOING THAT NOW
    entries = []
    for taxid in node_annotations:
        #print(taxid, node_annotations[taxid])
        thistaxid = int(taxid)
        thisnodename = ncbi.get_taxid_translator([thistaxid]).get(thistaxid, "Unknown")
        # get the NCBI taxid lineage for this node
        thislineage = ";".join([str(x) for x in ncbi.get_lineage(thistaxid)])
        # get the set of species that belong to this node
        these_species = node_annotations[thistaxid]
        # get a sub table of the perspchrom table that only has these species in the species column
        subdf = perspchrom[perspchrom["species"].isin(these_species)]
        # sum up the dataframe to get the number of fusions for each ALG, get rid of all the other columns
        subdf = subdf.drop(columns=["species", "taxid", "taxidstring"])
        subdf = subdf.sum(axis=0)
        # now add the other information as new columns, thistaxid/thisnodename/thislineage
        subdf["taxid"] = thistaxid
        subdf["nodename"] = thisnodename
        subdf["taxidstring"] = thislineage
        subdf["spinthisclade"] = ",".join(these_species)
        entries.append(subdf.copy())
    # condense all of the entries into a single df
    per_node_df = pd.DataFrame(entries)
    # move the taxid, nodename, thislineage columns to the front
    per_node_df.insert(0, "spinthisclade", per_node_df.pop("spinthisclade") )
    per_node_df.insert(0, "taxidstring",   per_node_df.pop("taxidstring")   )
    per_node_df.insert(0, "nodename",      per_node_df.pop("nodename")      )
    per_node_df.insert(0, "taxid",         per_node_df.pop("taxid")         )
    # save this!
    per_node_df.to_csv("all_nodes_ALG_presence_fusions.tsv", sep='\t', index=False)

    # make figures of the per-species plots
    standard_plot_out(per_node_df, "perNode")
    # Nodes we want to compare:
    #  10197 - Ctenophora
    #  6040 - Sponges
    #  10226 - Placozoa
    #  6073 - Cnidaria
    #  6231 - Nematodes
    #  88770 - Panarthropoda
    #  2697495 - Spiralia
    #  7711 - Chordata
    #  7586 - Echinodermata
    #  10219 - Hemichordata
    # pull out a df of just these nodes. Use an exact match of the taxid column
    nodes_to_compare = [10197, 6040, 10226, 6073, 6231, 88770, 2697495, 7711, 7586, 10219]
    comparison_set_df = per_node_df[per_node_df["taxid"].isin(nodes_to_compare)]
    # Sort the comparison set by the taxidstring of the nodes_to_compare list.
    # Use the order of numbers in nodes_to_compare to sort the dataframe.
    comparison_set_df = comparison_set_df.sort_values(by=["taxidstring"], ascending=True,
                                                      key=lambda x: x.map(dict(zip(nodes_to_compare, range(len(nodes_to_compare))))))
    print("The species in the node comparison are:")
    print(comparison_set_df)
    standard_plot_out(comparison_set_df, "comparisonNodes")

    # print the tree to a .tre file
    tree = ncbi.get_topology([int(x) for x in perspchrom["taxid"].tolist()])

    #taxid_to_query = 2301116
    #species_under_taxid = get_species_under_taxid(taxid_to_query)
    #print(f"Species under taxid {taxid_to_query}: {species_under_taxid}")

    ## remove rows where all values are 0
    #df = df.loc[(df!=0).any(axis=1)]
    #pca = PCA(n_components=2)
    #pca.fit(df)
    #print(pca.components_)
    #print(pca.explained_variance_)

    #df2 = pd.DataFrame(pca.transform(df), columns = ['first', 'second'])
    #print(df)
    #print(df2)
    ## df2.plot.scatter(x = 'first', y = 'second')

    ##plt.show()

    #from mpl_toolkits import mplot3d
    #from mpl_toolkits.mplot3d import Axes3D
    #from mpl_toolkits.mplot3d import proj3d
    #from matplotlib.text import Annotation

    #x = df2['first']
    #y = df2['second']
    ## labels is a list of the df1 index values
    #labels = df.index.values.tolist()

    ## Create the scatter plot
    #fig, ax = plt.subplots()
    #scatter = ax.scatter(x, y, picker=True)
    ## plot the text labels
    #for i, txt in enumerate(labels):
    #    ax.annotate(txt, (x[i], y[i]))


    ## Show the plot
    #plt.show()


if __name__ == '__main__':
    main()