#!/usr/bin/env python

"""
Program  : plot_ALG_fusions_v2.py
Language : python
Date     : 2024-02-07
Author   : Darrin T. Schultz
Email    : darrin.schultz@univie.ac.at
Github   : https://github.com/conchoecia/odp
Support  : For issues or questions, please search if the topic has been discussed already
           on github and open a new issue if not: https://github.com/conchoecia/odp/issues
License  : GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007. See the LICENSE file.
Citation : If you use this software for your scientific publication, please cite:
           Schultz, DT; Haddock, SHD; Bredeson, JV; Green, RE; Simakov, O & Rokhsar, DS
           Ancient gene linkages support ctenophores as sister to other animals. Nature (2023).
           https://doi.org/10.1038/s41586-023-05936-6

Description:
  This program is an updated version of the plot_ALG_fusions.py.
  In this version, we construct event string for each species. We do not bother to make the rest of the table.
    The rest of the table is the part that contains the matrix of whether or not an ALG is colocalized on the same chromosome in that species.

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
from   sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import networkx as nx
import itertools
import os
import pandas as pd
import re
import sys
import time
import umap
#import warnings
#warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

from PIL import Image
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster import hierarchy

# import odp-specific functions
import rbh_tools

# import the stuff to work with lineages
from ete3 import NCBITaxa,Tree

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

def rgb_255_float_to_hex(rgb_floats):
    """
    Converts a single rgb 0-255 to a hex string.
    """
    return '#%02x%02x%02x' % (int(rgb_floats[0]), int(rgb_floats[1]), int(rgb_floats[2]))


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

def taxids_to_taxidstringdict(taxids) -> dict:
    """
    This function takes a list of taxids and returns a dictionary where the key is the species name and the value is the taxid string.
    The taxid string will be a string of taxids separated by a semicolon.
    """
    ncbi = NCBITaxa()

    # check that taxids is an iterable
    acceptable_iterables = [list, set, dict]
    thistype = type(taxids)
    if thistype not in acceptable_iterables:
        raise ValueError("The taxids must be an iterable of taxids that we can look through.")

    # make sure that all the taxids are interpretable as ints
    for taxid in taxids:
        if not re.match(r"^[0-9]*$", str(taxid)):
            raise ValueError("There is a non-numeric character in the taxid string")

    # This is the dict that we will return
    taxid_to_taxidstring = {}
    for taxid in taxids:
        # get the lineage of the taxid
        lineage = ncbi.get_lineage(taxid)
        # Return the complete taxid string, text delimited by a semicolon
        returnstr = ";".join([str(x) for x in lineage])
        taxid_to_taxidstring[taxid] = returnstr
    return taxid_to_taxidstring

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

    Notes:
      - 20240208: One finding was that sometimes, if a clade had a lot of losses, then the
        ALG fusions were not pushed back to the correct node. I am making modifications to
        not report something as split if it is not detectable in the first place.
    """
    # First get all of the ALGs from this dataframe by finding the non-tuple
    #   columns that are not in ["species", "taxid", "taxidstring"]
    remove_these_columns = ["species", "taxid", "taxidstring", "changestrings"]
    pair_columns = [x for x in df.columns if isinstance(x, tuple) and x not in remove_these_columns]

    pair_QC = {}
    for pair in pair_columns:
        # get the species that have both of these ALGs detectable
        both_detectable = df[(df[pair[0]] > 1) & (df[pair[1]] > 1)]
        if not both_detectable.empty:
            pair_QC[pair] = both_detectable[pair].mean()

    # return the pairs that are not colocalized in at least min_for_noncolocalized fraction of the species
    return [x for x in pair_QC if pair_QC[x] <= min_for_noncolocalized]

def unsplit_ALGs(df, max_frac_split = 0.5):
    """
    This takes in a dataframe and returns a list of ALGs that appear to be unsplit in this clade.

    When polarized against a clade for which we know the ALGs are split, we can find the branches on which
     the changes occurred.

    The max_frac_unsplit is the maximum fraction of species in this clade that can have the ALG split
     before declaring the ALG to be split in this clade - thereby not returning it as unsplit.
    """
    # First we get all the columns that are not tuples
    remove_these_columns = ["species", "taxid", "taxidstring", "changestrings"]
    ALG_columns = [x for x in df.columns if not isinstance(x, tuple)
                   and (x not in remove_these_columns)]

    ALG_qc = {}
    for ALG in ALG_columns:
        # get the fraction of species that have this ALG split across two or more.
        if len(df[df[ALG] >= 1]) > 0:
            ALG_qc[ALG] = len(df[df[ALG] > 1]) / len(df[df[ALG] >= 1])
    return [x for x in ALG_qc if ALG_qc[x] <= max_frac_split ]

def rbh_files_to_locdf_and_perspchrom(rbh_files, ALGrbhfile, minsig, ALGname) -> (pd.DataFrame, pd.DataFrame):
    """
    This takes in a list of rbh files and returns two dataframes.
    The RBH files are those that are species against ALGs.

    The locdf looks like this:
                                            sample gene_group    scaffold        pvalue  num_genes  frac_of_this_ALG_on_this_scaffold
        allooctoploidhybrid-2876849-GCA024542945.1        A1a  CM045036.1  2.713967e-05         15                           0.082873
        allooctoploidhybrid-2876849-GCA024542945.1        A1a  CM045038.1  2.440960e-02         11                           0.060773
        allooctoploidhybrid-2876849-GCA024542945.1        A1a  CM045052.1  1.984298e-04         12                           0.066298
        allooctoploidhybrid-2876849-GCA024542945.1        A1a  CM045053.1  9.969622e-13         25                           0.138122
        allooctoploidhybrid-2876849-GCA024542945.1        A1a  CM045054.1  3.496336e-08         16                           0.088398
                                               ...        ...         ...           ...        ...                                ...
                    Zeusfaber-64108-GCA960531495.1         Qd  OY482860.1  8.905043e-04          8                           0.400000
               Zeugodacustau-137263-GCA031772095.1         Ea  CM062648.1  1.620573e-03          9                           0.115385
               Zeugodacustau-137263-GCA031772095.1         Eb  CM062648.1  4.219110e-02          5                           0.161290
               Zeugodacustau-137263-GCA031772095.1          G  CM062650.1  2.983613e-02         34                           0.326923
               Zeugodacustau-137263-GCA031772095.1          I  CM062649.1  3.457063e-09         38                           0.633333
    """
    # Check that the list of rbh files is not empty
    if len(rbh_files) == 0:
        raise IOError("The list of rbh files is empty.")
    # Check that all of the files in the rbh_files list exist
    for file in rbh_files:
        if not os.path.exists(file):
            raise IOError(f"The file {file} does not exist.")
    # check that the ALGrbhfile exists
    if not os.path.exists(ALGrbhfile):
        raise IOError(f"The file {ALGrbhfile} does not exist.")
    # if minsig is greater than 0.05, then raise an error telling the user that the value is too high
    if minsig > 0.05:
        raise ValueError("The minsig value is too high. It should be less than 0.05.")
    # OK, we're done being paranoid. Let's get to work.

    # first we need to read through all of the rbh files to get all of the possible ALGs to plot
    ALGdf = rbh_tools.parse_ALG_rbh_to_colordf(ALGrbhfile)

    # First we figure out on which chromosomes the ALGs are located. Some may be split.
    sample_to_chromnum    = {}
    sample_to_taxidstring = {}
    sample_to_taxid       = {}
    # This is just used later to concatenate all of the results
    entries = []
    for i in range(len(rbh_files)):
        print(f"\r  Parsing rbh file {i+1}/{len(rbh_files)}", end = "", file = sys.stderr)
        file = rbh_files[i]
        rbhdf = rbh_tools.parse_rbh(file)
        splitdf, samplename = rbh_tools.rbhdf_to_alglocdf(rbhdf, minsig, ALGname)
        chromnum = rbh_tools.rbh_to_scafnum(rbhdf, samplename)
        sample_to_chromnum[samplename] = chromnum
        entries.append(splitdf)

        # we know where the NCBI taxid will be in the file name. Just extract it.
        taxid = samplename.split("-")[1]
        # check that the taxid is an integer
        if not re.match(r"^[0-9]*$", str(taxid)):
            raise ValueError("There is a non-numeric character in the taxid string")
        sample_to_taxid[samplename] = int(taxid)
    print()

    # convert the entries into a dataframe
    locdf = pd.concat(entries).reset_index(drop=True)

    # Now that we know the NCBI taxid for each sample, generate the taxid_to_lineagestring
    sample_to_taxidstring = taxids_to_taxidstringdict([sample_to_taxid[k] for k in sample_to_taxid])

    # now make a pandas dataframe of the sample_to_taxidstring. The columns are "sample", "taxid"
    perspchrom = pd.DataFrame.from_dict(sample_to_taxid, orient='index')
    # change the index to a column. the former index is "species", the other column is "taxid"
    perspchrom = perspchrom.reset_index()
    perspchrom = perspchrom.rename(columns={"index": "species", 0: "taxid"})
    # the taxidstring is the sample_to_taxidstring dictionary. use map.
    perspchrom["taxidstring"] = perspchrom["taxid"].map(sample_to_taxidstring)
    perspchrom = perspchrom.sort_values(by="taxidstring", ascending=True).reset_index(drop=True)

    # how many queries do we need to make?
    total_queries = len(ALGdf)
    counter = 0
    for i in range(len(ALGdf)-1):
        for ii in range(i+1, len(ALGdf)):
            total_queries += 1

    # ┏┓┓ ┏┓  ┏┓┳┓┏┓┏┓┏┓┳┓┏┓┏┓ ╻ ┏┓┳┓┏┓┏┓┳┓┏┓┏┓  ┏┓┏┓┓ ┳┳┳┳┓┳┓┏┓
    # ┣┫┃ ┃┓  ┃┃┣┫┣ ┗┓┣ ┃┃┃ ┣ ━╋━┣┫┣┫┗┓┣ ┃┃┃ ┣   ┃ ┃┃┃ ┃┃┃┃┃┃┃┗┓
    # ┛┗┗┛┗┛  ┣┛┛┗┗┛┗┛┗┛┛┗┗┛┗┛ ╹ ┛┗┻┛┗┛┗┛┛┗┗┛┗┛  ┗┛┗┛┗┛┗┛┛ ┗┛┗┗┛
    ## make a new column for each of the ALGs, sorted by largest to smallest
    ##  This will be the presence/absence of the ALGs in the sample.
    ##  The presence/absence of the ALGs will be used to make a pixel-wise plot of the ALGs.
    ##  The presence/absence of the ALGs will also be used to make a plot of the ALG colocalization matrix.
    ## If the value is 0, then the ALG is not present on any chromosomes.
    ## If the value is 1, then the ALG is present on only one chromosome.
    ## If 2, the ALG is present on two chromosomes, etc.
    sorted_ALG_list = list(ALGdf.sort_values("Size", ascending = False)["ALGname"])
    # add these columns to perspchrom, initialize with 0
    results = []
    # go through all the rows in perspchrom
    numsamp = len(perspchrom)
    for i, row in perspchrom.iterrows():
        print(f"\r  Analyzing the ALG composition of sample: {i}/{numsamp}          ", end = "", file = sys.stderr)
        thissample = row["species"]
        results.append({x: 0 for x in sorted_ALG_list})
        # update with the valuecounts for the gene_group column
        results[-1].update(locdf[locdf["sample"] == thissample]["gene_group"].value_counts().to_dict())
        results[-1]["species"] = thissample
    print()
    # make a dataframe from the results list. The columns are the ALGs and the rows are the samples. Use the order of the sorted_ALG_list
    resdf = pd.DataFrame(results, columns = ["species"] + sorted_ALG_list)
    # Merge the resdf with perspchrom.
    perspchrom = pd.merge(perspchrom, resdf, on = "species")

    # ┏┓┓ ┏┓  ┏┓┏┓┓ ┏┓┏┓┏┓┓ ┳┏┓┏┓┏┳┓┳┏┓┳┓  ┏┓┏┓┓ ┳┳┳┳┓┳┓┏┓
    # ┣┫┃ ┃┓  ┃ ┃┃┃ ┃┃┃ ┣┫┃ ┃┏┛┣┫ ┃ ┃┃┃┃┃  ┃ ┃┃┃ ┃┃┃┃┃┃┃┗┓
    # ┛┗┗┛┗┛  ┗┛┗┛┗┛┗┛┗┛┛┗┗┛┻┗┛┛┗ ┻ ┻┗┛┛┗  ┗┛┗┛┗┛┗┛┛ ┗┛┗┗┛
    # These columns mark if the ALGs are present on the same chromosomes
    columnnames = []
    # now add a column for each of the ALG pairs
    for i in range(len(sorted_ALG_list)-1):
        for ii in range(i+1, len(sorted_ALG_list)):
            Ai = sorted_ALG_list[i]
            Aii = sorted_ALG_list[ii]
            columnnames.append(tuple(sorted((Ai, Aii))))
    results = []
    # now, for each species go through and add a 1 to the column if the ALGs are present on the same chromosome
    for i, row in perspchrom.iterrows():
        print(f"\r  Analyzing the colocalizations of sample: {i+1}/{numsamp}          ", end = "", file = sys.stderr)
        results.append({x: 0 for x in columnnames})
        results[-1]["species"] = row["species"]
        gb = locdf[locdf["sample"] == row["species"]].groupby("scaffold")
        # groupby the chromosome number, if there are groups that have multiple ALGs, then add a 1 to the appropriate column
        for name, group in gb:
            # get the ALGs in the group
            ALGs_in_group = list(set(group["gene_group"]))
            # go through the combinations and add a value to the perspchrom dataframe
            for j in range(len(ALGs_in_group)-1):
                for jj in range(j+1, len(ALGs_in_group)):
                    Ai = ALGs_in_group[j]
                    Aii = ALGs_in_group[jj]
                    thiscolname = tuple(sorted((Ai, Aii)))
                    results[-1][thiscolname] += 1
    print()
    # make a df from the results dataframe. Use the columnnames list for the order.
    resdf = pd.DataFrame(results, columns = ["species"] + columnnames)
    # now merge back with the perspchrom dataframe
    perspchrom = pd.merge(perspchrom, resdf, on = "species")

    return locdf, perspchrom

def assign_colors_to_nodes(graph, root, colors):
    """
    Assigns colors to nodes in a directed acyclic graph such that colors approach the average of all colors as 
    you move towards the root.

    Parameters:
        graph (nx.DiGraph): Directed acyclic graph.
        root: Root node of the graph.
        colors (dict): Dictionary mapping node names to colors in hexadecimal notation.
    """

    def calculate_average_color(node):
        """
        Recursively calculates the average color of a node and its daughters.
        """
        daughters = list(graph.successors(node))
        if not daughters:
            return colors[node]

        daughter_colors = [calculate_average_color(daughter) for daughter in daughters]
        avg_color = np.mean(daughter_colors, axis=0)
        return avg_color

    def assign_color(node, parent_avg_color):
        """
        Assigns a color to a node based on the average color of its daughters and the parent's average color.
        """
        if node != root:
            avg_color = calculate_average_color(node)
            # Interpolate between the average color of daughters and the parent's average color
            alpha = 0.5  # Adjust this parameter for the rate of interpolation
            node_color = interpolate_color(avg_color, parent_avg_color, alpha)
            colors[node] = node_color
        else:
            # For the root node, use its average color directly
            colors[node] = parent_avg_color

        for daughter in graph.successors(node):
            assign_color(daughter, colors[node])

    def interpolate_color(color1, color2, alpha):
        """
        Interpolates between two colors in hexadecimal notation.
        """
        ## Interpolate between RGB triples
        #print("alpha ",   alpha)
        #print("color1 ", color1)
        #print("color2 ", color2)
        #print("(alpha * color1) ", (alpha * color1))
        #print("(1 - alpha) * color2) ", (1 - alpha) * color2)

        interpolated_rgb = (alpha * color1) + ((1 - alpha) * color2)
        return interpolated_rgb

    # Start assigning colors from the root
    root_color = calculate_average_color(root)
    assign_color(root, root_color)

def save_UMAP_plotly(tree_df, outprefix):
    import plotly.express as px
    # only get the rows where the '-' characters is in the index
    for mode in ["withnodes", "withoutnodes"]:
        tempdf = tree_df.copy()
        outhtml = None
        if mode == "withoutnodes":
            tempdf = tempdf[tempdf.index.str.contains("-")]
            outhtml = f"{outprefix}.withoutnodes.html"
        elif mode == "withnodes":
            outhtml = f"{outprefix}.withnodes.html"
        else:
            raise ValueError("The mode must be either 'withnodes' or 'withoutnodes'")
        plotdf = tempdf.drop(columns = ["color"])
        X = plotdf.values
        print("This is X")
        print(X)

        # Add random noise to the data
        noise_level = 0.1  # Adjust the noise level as needed
        X_noisy = X + np.random.normal(0, noise_level, X.shape)

        # Apply UMAP to reduce dimensionality
        reducer = umap.UMAP()
        embedding = reducer.fit_transform(X_noisy)

        # Create a DataFrame with the embedding
        df_embedding = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'])

        # Add the indices as labels
        df_embedding['label'] = plotdf.index

        # Add colors to the plot
        # Assuming you have a 'color' column in your DataFrame indicating the color of each point
        color_column = tempdf['color']  # Replace 'color' with the actual column name containing colors
        fig = px.scatter(df_embedding,
                         x='UMAP1', y='UMAP2',
                         hover_name='label', color=color_column)
        # Show the plot
        fig.write_html(outhtml)

def save_UMAP(tree_df):
    """
    This method saves the UMAP of the tree.
    """
    # make a UMAP of the self.tree_df, each point on the UMAP is one row of the dataframe, which is one node in the tree
    # make a UMAP of the self.tree_df
    # remove the colors dataframe
    tempdf = tree_df.drop(columns = ["color"])
    X = tempdf.values
    print("This is X")
    print(X)

    # Add random noise to the data
    noise_level = 0.1  # Adjust the noise level as needed
    X_noisy = X + np.random.normal(0, noise_level, X.shape)

    # Apply UMAP to reduce dimensionality
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(X_noisy)

    # Plot the reduced data, add the colors from the dataframe
    #plt.scatter(embedding[:, 0], embedding[:, 1], s=5)  # Adjust 's' for the size of points
    plt.scatter(embedding[:, 0], embedding[:, 1], s=5, c=tree_df["color"])  # Adjust 's' for the size of points
    plt.gca().set_aspect('equal', 'datalim')  # Set equal aspect ratio
    plt.title('UMAP visualization')
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    plt.show()

class SplitLossColocTree:
    """
    This class is used to store a phylogenetic tree.
    The tree is implemented as a directional graph.
    The tree can be constructed by lists of edges. The nodes are inferred from the edges.
    The lineage is input as a string of taxids, separated by semicolons, plus the sample name.

    The point of this class;
      - Each leave will have a state of how the ALGs are configured in that sample.
        Using this information, we will infer the state of the ALGs on the internal nodes.
      - With the internal information, we will determine on which branches different events happened.

    # Every node has the following properties:
      - completed: a boolean that is True if we know the state of the ALGs at this node.
      -
    """
    color_dict_top = {
                  33317:     "#F4B93E", # Protostomes
                   1215728:  "#2ECAC8", # Scalidophora
                   6231:     "#9F2ECA", # Nematodes
                   88770:    "#BF2424", # Panarthropoda
                    32341:   "#F93AD5", # Sophophora - these are the fruit flies
                    41084:   "#F1890D", # Polyphaga - these are the beetles
                    3042114: "#EBB836", # Anthophila - these are the bees
                   2697495:  "#2432BF", # Spiralia
                    6178:    "#24ED54", # Trematoda
                    6447:    "#0CD6F2", # Mollusca
                     6606:   "#F2B70C", # Coleoidea
                    6340:    "#0CF247", # Annelida
                   2697496:  "#7624BF", # Gnathifera
                  33511:   "#78A6AF", # Deuterostomes
                   # Chordates
                    # Vertebrates
                     # Cyclostomes
                     1476529: "#F23341", # Jawless fishes - hagfish lampreys
                     # Gnathostomes
                      # Cartilaginous fishes - Condrichthyes
                       7777: "#78E576", # Cartilaginous fishes - sharks, rays, chimaeras
                      # Bony fishes - Osteichthyes
                       # Actinopterygii - Ray-finned fishes
                        1338366: "#DEC06A", # Cladistia - Bichirs, Reedfish
                        32440:   "#636E42", # Chondrostei - Sturgeons, Paddlefish
                        # Teleost fish
                        41665:   "#BD5812", # Neopterygii - Teleost fish, >32k living species
                          186628:  "#CF99FF", # Characiphysae - some group of fish, including Astanyx, the blind cave fish
                          7952:  "#4BE2BA", # Cypriniformes - these are the minnows, carps, loaches, et cetera
                       # Sarcopterygii - Lobe-fins (fish and tetrapods)
                        118072: "#445FCF", # Actinistia - Coelacanths
                        7878: "#DEC06A", # Dipnoi - Lungfish
                        # Tetrapods
                         # Amniota
                          9443:  "#F58E8C",   # Primates
                          33554: "#8A0303",   # Carnivora
                          1579337: "#2C4BD8", # Durocryptodira - turtles
                          9263:  "#5FA4E2",   # Metatheria - marsupials
                          9126:  "#242B7D",   # Passeriformes - more than half of all bird species
                          8826:  "#16A75C",   # Anseriformes - waterfowl
                          35500: "#4E4EC4",   # Pecorids - these are sheep, bovines, deer, muntjac deer, et cetera
                         # Amphibia
                          8342:  "#AFC73E",
                  10197:   "#54AB53", # Ctenophores
                  6040 :   "#DCC0F3", # Sponges
                  6073 :   "#387FB2", # Cnidarians
                  10226:   "#C72480", # Placozoans
                  }
    def __init__(self, perspchrom) -> None:
        # initialize the bidirectional graph using networkx
        self.G = nx.DiGraph()
        self.perspchrom = perspchrom
        # Before we delete anything, make a dict of samples to the taxidstring
        self.sample_to_taxidlist = self.perspchrom.set_index("species")["taxidstring"].to_dict()
        self.sample_to_taxidlist = {k: [int(x) for x in v.split(";")] for k, v in self.sample_to_taxidlist.items()}
        self._build_tree_from_perspchrom()

        # define the ALG nodes and the combination nodes
        # Just define the ALG columns now
        reject_these = ['species', 'taxid', 'taxidstring']
        self.ALGcols = [ x for x in self.perspchrom if not isinstance(x, tuple) and x not in reject_these]
        self.TUPcols = [ x for x in self.perspchrom if     isinstance(x, tuple) and x not in reject_these]
        # add the completed property to all of the nodes
        for node in self.G.nodes():
            if node in self.perspchrom["species"].values:
                self.G.nodes[node]["completed"] = True
            else:
                self.G.nodes[node]["completed"] = False
        # keep track of all the leaves as places to start
        self.leaves = [x for x in self.G.nodes() if self.G.out_degree(x) == 0]

        # reindex the perspchrom dataframe with the species column
        self.perspchrom = self.perspchrom.set_index("species")
        # drop the columns taxid and taxidstring
        self.perspchrom = self.perspchrom.drop(columns = ["taxid", "taxidstring"])
        # copy the first row of perspchrom to be the empty predecessor. Get it as a df, not a series.
        self.empty_predecessor = self.perspchrom.iloc[[0]].copy().reset_index(drop=True)

        # set all the values to -1
        for col in self.empty_predecessor.columns:
            self.empty_predecessor[col] = -1
        print(self.empty_predecessor)

        # for all of the leaves, add the datatypes
        for thisnode in self.leaves:
            # make sure it returns a dataframe and not a series
            self.G.nodes[thisnode]["dataframe"] = self.perspchrom.loc[[thisnode]].copy()

        # assign the colors for all the nodes
        self.node_to_color = {}
        self.assign_colors()

        # Final df for the tree
        self.tree_df = None

    def assign_colors(self):
        """
        Assigns colors to the nodes based on some preferences.
        """
        node_colors = {}
        # go through the class variable, color_dict, and find all the leaves, then assign the colors
        for sample in self.sample_to_taxidlist:
            # first do the top-level colors
            for thistop in self.color_dict_top:
                if thistop in self.sample_to_taxidlist[sample]:
                    node_colors[sample] = self.color_dict_top[thistop]
            # then the subclade colors
            for thissub in self.color_dict:
                if thissub in self.sample_to_taxidlist[sample]:
                    node_colors[sample] = self.color_dict[thissub]
        # convert the node_colors to np arrays
        node_colors = {node: np.array(hex_to_rgb(color)) for node, color in node_colors.items()}
        # go through the leaves, and if the color is not assigned, give it a non-offensive blue "#3f3f7f"
        for leaf in self.leaves:
            if leaf not in node_colors:
                node_colors[leaf] = np.array(hex_to_rgb("#3f3f7f"))

        # Assign colors to nodes
        root = [n for n,d in self.G.in_degree() if d==0][0]
        assign_colors_to_nodes(self.G, root, node_colors)

        ## Print the assigned colors for each node
        #print("Node Colors:")
        ## results
        ##                       #E: #2faf1f
        ##           #B: #5f5f3f #D: #af2f1f
        ## A: #3f3f7f
        ##           #C: #1f1fbf #F: #0f0fdf
        #for node, color in node_colors.items():
        #    newcolor = rgb_255_float_to_hex(color)
        #    print(f"{node}: {newcolor}")

        # go through the graph and add a color to each node
        for node in self.G.nodes():
            if node not in node_colors:
                raise IOError(f"The node {node} does not have a color assigned to it.")
            else:
                self.G.nodes[node]["color"] = rgb_255_float_to_hex(node_colors[node])
                self.node_to_color[node]    = rgb_255_float_to_hex(node_colors[node])

    def save_tree_to_df(self, filename):
        """
        This method saves the rows of the tree as a dataframe
        """
        # get all the nodes in the tree
        all_nodes = list(self.G.nodes())
        all_dfs = [self.G.nodes[x]["dataframe"] for x in all_nodes]
        all_dfs = pd.concat(all_dfs)
        # change the columns of self.ALGcols to be ints
        for col in self.ALGcols:
            all_dfs[col] = all_dfs[col].astype(int)

        # Add the colors column to the dataframe. Use concat because pandas
        #  complains about performance if I try to add a column using
        #  older spells: a la all_dfs["color"] = colors
        colors = [self.node_to_color[x] for x in all_dfs.index]
        all_dfs = pd.concat([pd.DataFrame({'color': colors}, index=all_dfs.index),
                             all_dfs],
                             axis=1)

        # save to self.tree_df
        self.tree_df = all_dfs
        # save as a tsv
        # When the floats were allowed to be any length, the tree size was 17 MB.
        # When the floats were limited to 3 decimal places.
        self.tree_df.to_csv(filename, sep = "\t")

    def _get_predecessor(self, node):
        """
        Returns the predecessor of the node. If the node is the root, then returns None.
        """
        predecessors = list(self.G.predecessors(node))
        if len(predecessors) == 0:
            # we don't need to do anything because we're at the end of the tree.
            return None
        elif len(predecessors) > 1:
            raise IOError(f"There should only be one predecessor in a phylogenetic tree. Found {predecessors}")
        elif len(predecessors) == 1:
            # there is only one parent
            parent_node = list(self.G.predecessors(node))[0]
            return parent_node

    def percolate(self):
        """
        This function goes through the leaves in a breadth-first search.
        Eventually, randomnly choose a starting order as another source of stochasticity in the algorithm.

        To stochasitcally sample the trees, we must choose probabilities of changes in states.

        Probabilities we must choose:
          - ALG losses:
            - When an ALG is absent in node A, but is present in sister node(s) B, we must decide whether the ALG was present
              or not in the parent of A and B. Given that we know that these ALGs were present in the ancestor of animals, the most
              likely scenario is that the ALG was present in the parent of A and B. The directionality of ALGs changing over time
              is simply loss. We express this probability as (pALGeP|!A and eB) = 0.9999.
              In other words, we choose that the ALG will appear in the parent of A and B in 99.99% of the cases.
                    |
                ----P----      (pALG eP | eA and eB) = 0.99999  (variable name prob_eP_eAeB)
                |   |   |  AND
                |   |   |      (pALG eP | !A and |eB| ) = 1 - ((1 - 0.99999) ** |eB|)
                xA  eB  eB      In other words, the probability of the ALG being present in the parent increases with every observation.


            - When the ALG is not present in A or B, the probability that it is present in the parent is negligibly low.
              We express this probability as (pALGeP| !A and !B) = 0.00001. This means that 0.001% of the time, the ALG will be present
              in the parent.
                     |
                   --P--  (pALG eP | !A and !B) = 0.00001  (variable name prob_eP_eAeB)
                   |   |
                   |   |
                   xA  xB

            - If there are no siblings to check, we can do nothing but to inherit the state of the one existing child node, A.
              There is no probability to calculate in this case. This is the only option when there are no siblings.
                     |
                   --P--  (pP | A ) = 1  (The parent just inherits the state of the child. Because B is missing.)
                   |
                   |
                   xA
        """
        prob_eP_xAeB = 0.99999
        prob_eP_xAxB = 0.00001

        # start at the leaves
        queue = list(self.leaves)
        while len(queue) > 0:
            #print(f"Queue length: {len(queue)}, queue: {queue}", file = sys.stderr)
            thisnode = queue.pop(0)

            # ┏┓┏┓┏┳┓  ┏┓┏┓┳┓┏┓┳┓┏┳┓ - We need to get the parent node of thisnode.
            # ┃┓┣  ┃   ┃┃┣┫┣┫┣ ┃┃ ┃    Then we must determine whether the ancestral state of the parent node
            # ┗┛┗┛ ┻   ┣┛┛┗┛┗┗┛┛┗ ┻    has already been determined.
            # get sister nodes. Do this by getting the node from the in edges
            predecessor = self._get_predecessor(thisnode)
            # If the predecessor is None, then this means we're at the end of the tree, and we can't do anything.
            # If we have already completed the predecessor node, then that means it was visited from another sibling node.
            # Therefore, we don't need to analyze this node.
            if (predecessor is not None) and (self.G.nodes[predecessor]["completed"] is False):
                # We should check that the predecessor node does not have a dataframe. If it does, but was not marked as completed,
                #  then this means that there was some problem with the algorithm. It is important that we mark the predecessor node as
                #  completed if we modify the dataframe.
                if "dataframe" in self.G.nodes[predecessor]:
                    raise IOError(f"The predecessor node {predecessor} has a dataframe, but was not marked as completed.")
                #DEBUG
                ## Now we know that we will modify the parent node. We can make a copy of the dataframe of thisnode and set all the
                ##  the values to -1. this way, we will know if we have modified the value or not.
                ## convert it to a series since it is just one row.
                #if predecessor == "390379":
                #    print(f"{thisnode} is the node and the predecessor is {predecessor}", file = sys.stderr)
                #    print(f"The dataframe of this node is:")
                #    print(self.G.nodes[thisnode]["dataframe"], file = sys.stderr)
                #    sys.exit()
                # make an empty dataframe
                # set the index to the predecessor id
                predecessordf = self.empty_predecessor.copy()
                predecessordf.index = [predecessor]

                # ┏┓┳┳┓┓ ┳┳┓┏┓┏┓ - We have determined there are multiple sibling clades. We get them
                # ┗┓┃┣┫┃ ┃┃┃┃┓┗┓   and then we can determine the state of the parent node.
                # ┗┛┻┻┛┗┛┻┛┗┗┛┗┛
                # get the other direct descendants of the parent
                siblings = [x for x in self.G.successors(predecessor) if x != thisnode]
                # If there are no siblings, then the parent node inherits the state of this node.
                # We then add the parent to the queue and mark it as completed.
                if len(siblings) == 0:
                    self.G.nodes[predecessor]["dataframe"] = self.G.nodes[thisnode]["dataframe"].copy()
                    # set the index
                    self.G.nodes[predecessor]["dataframe"].index = [predecessor]
                    self.G.nodes[predecessor]["completed"] = True
                    queue.append(predecessor)
                # Now we filter to just get the siblings that are completed
                siblings = [x for x in siblings if self.G.nodes[x]["completed"]]
                if len(siblings) == 0:
                    # For this node, there are no siblings that are completed, so we can't do anything.
                    # We have to wait until the siblings are completed. We will come to this node later.
                    # Add it to the back of the queue.
                    queue.append(thisnode)
                else:
                    # let's make a tempdf of this dataframe and the sibling dataframes
                    siblingdf = pd.concat([self.G.nodes[thisnode]["dataframe"]] + [self.G.nodes[x]["dataframe"] for x in siblings])

                    values = {}
                    # We first address the presence/absence of the ALGs in the parent node.
                    # ┏┓┓ ┏┓  ┏┓┳┓┏┓┏┓┏┓┳┓┏┓┏┓  ┏┓┳┓┳┓  ┏┓┳┓┏┓┏┓┳┓┏┓┏┓
                    # ┣┫┃ ┃┓  ┃┃┣┫┣ ┗┓┣ ┃┃┃ ┣   ┣┫┃┃┃┃  ┣┫┣┫┗┓┣ ┃┃┃ ┣
                    # ┛┗┗┛┗┛  ┣┛┛┗┗┛┗┛┗┛┛┗┗┛┗┛  ┛┗┛┗┻┛  ┛┗┻┛┗┛┗┛┛┗┗┛┗┛
                    # For evey ALG, we apply the logic above. Use the helper method to do this.
                    pdf = self._determine_parental_ALG_PresAbs(predecessordf, siblingdf,
                                                               prob_eP_xAeB = prob_eP_xAeB,
                                                               prob_eP_xAxB = prob_eP_xAxB)

                    # ┏┓┓ ┏┓ ┏┓┏┓┓ ┳┏┳┓┏┓
                    # ┣┫┃ ┃┓ ┗┓┃┃┃ ┃ ┃ ┗┓ - We now infer what the number of ALGs was at each node.
                    # ┛┗┗┛┗┛ ┗┛┣┛┗┛┻ ┻ ┗┛
                    # For evey ALG, we apply the logic above. Use the helper method to do this.
                    pdf = self._determine_parental_ALG_Splits(pdf, siblingdf)

                    # ┏┓┓ ┏┓  ┏┓┏┓┓ ┏┓┏┓┏┓
                    # ┣┫┃ ┃┓  ┃ ┃┃┃ ┃┃┃ ┗┓
                    # ┛┗┗┛┗┛  ┗┛┗┛┗┛┗┛┗┛┗┛
                    # now we should determine the state of the colocalized ALGs
                    pdf = self._determine_ALG_colocalization(pdf, siblingdf)

                    # add the parent to the queue, mark it as completed
                    self.G.nodes[predecessor]["dataframe"] = pdf
                    self.G.nodes[predecessor]["completed"] = True
                    queue.append(predecessor)

    def _conservation_of_colocalizations(self, df):
        """
        Takes in a dataframe and returns a dataframe of how often each tuple is conserved
        """
        results = {}
        for thistup in self.TUPcols:
            ALG1 = thistup[0]
            ALG2 = thistup[1]
            # The subdf the rows in which both ALG1 and ALG2 are present, get the value counts of the tuple
            subdf = df[(df[ALG1] >= 1) & (df[ALG2] >= 1)][thistup]
            if len(subdf) > 0:
                numconserved = len(subdf[subdf == 1])
                # round the next value to 3 decimal places
                results[thistup] = round(numconserved / len(subdf), 3)
            else:
                # We can't make an inference, so this colocalization gets a value of 0
                results[thistup] = 0
        assert len(results) == len(self.TUPcols)
        return results

    def _determine_ALG_colocalization(self, pdf, sdf):
        """
        Use the leaves to determine the state of the parent.
        """
        # get all the leaves from the parent
        predecessor_node = pdf.index[0]
        leaves = [x for x in nx.descendants(self.G, predecessor_node) if self.G.out_degree(x) == 0]
        # It is problematic if there are no leaves. Every predecessor's existence is predicated on a leaf's existence.
        if len(leaves) == 0:
            raise IOError(f"There are no leaves in the predecessor {predecessor_node}. This should not happen.")
        # get the dataframes of the leaves
        ldf = self.perspchrom.loc[leaves]
        results = self._conservation_of_colocalizations(ldf)
        # update the parent dataframe with the results. The results are tuple columns
        for thistup in self.TUPcols:
            pdf[thistup] = results[thistup]
        return pdf

    def _determine_parental_ALG_Splits(self, pdf, sdf):
        """
        Infers how many ALGs were at each node.

        This rule determines what the colocalization state of the parent nodes are.
        Like in the _determine_parental_ALG_Splits rule, we will do some filtering
          to pick out the higher-quality genomes.

        There are a few cases to handle.
        - If there is only one genome in the sdf, then the parent inherits this state.
        - If both of the values for the tuple are the same,
          it is easy to determine what the value of the parent should be.
        - If the values are different, look at the genomes of the sister clade.
        """
        sdf = self._filter_sdf_for_high_quality(sdf)

        # If the sdf dataframe has a length of 1, then we just multiply the values by the current values in the pdf.
        # Let's check quickly that the ALG columns in the pdf do not have any values that are -1. If they do, this means
        #  that we didn't finish assigning the values during the ALG presence/absence step.
        if (pdf[self.ALGcols] == -1).any().any():
            raise IOError(f"The pdf has -1 values in the ALG columns. This means that we didn't finish assigning the values during the ALG presence/absence step.")
        # Now we continue. If the sdf dataframe has a length of 1, then we just multiply the values by the current values in the pdf.
        if len(sdf) == 0:
            raise IOError(f"The sdf dataframe has a length of 0. This should not happen.")
        elif len(sdf) == 1:
            pdf[self.ALGcols] = pdf[self.ALGcols].multiply(sdf.iloc[0][self.ALGcols], axis = 1)
        elif len(sdf) > 1:
            # We must pick one of these numbers. For consistency, just pick an entire row from sdf and update the pdf.
            randindex = np.random.choice(sdf.index)
            pdf[self.ALGcols] = pdf[self.ALGcols].multiply(sdf.loc[randindex][self.ALGcols], axis = 1)

        return pdf

    def _filter_sdf_for_high_quality(self, sdf) -> pd.DataFrame:
        """
        This method is used to pick out the high-quality genomes from the dataframe.
        Returns a filtered dataframe.
        """
        # First we must infer if we are looking at leaves or not.
        # To check if we are looking at leaves, check if all the nodes from the index of sdf have any descendants in the graph.
        leaves = []
        for i in range(len(sdf)):
            thisnode = sdf.index[i]
            if thisnode in self.G.nodes():
                if self.G.out_degree(thisnode) == 0:
                    leaves.append(True)
                else:
                    leaves.append(False)
        if all(leaves):
            # There are many differences between the RefSeq and GenBank versions of the genome in the case of chromosome fusions.
            # Because in theory we should trust the RefSeq version more, we will look for cases where there are both versions,
            #  and we will only further consider the RefSeq version.
            # To find it. we we change one character in the accession number. GCF=RefSeq, GCA=GenBank.
            #
            # Here is an example from Takifugu flavidus, where there are differences between the GCA and GCF versions of the genome.
            #                                          A1a  D  F  C1  G  H  Ea  N  L  M  B1  I  B2  O1  A1b  Eb  O2  (A1a, A1b)  (D, O1)  (D, O2)  (Ea, O1)  (Ea, Eb)  (L, M)  (B1, B2)  (O1, O2)
            #  Takifuguflavidus-433684-GCF003711565.1    3  3  3   2  4  3   4  2  3  3   2  2   2   3    2   3   2           2        2        2         2         3       3         2         2
            #  Takifuguflavidus-433684-GCA003711565.2    3  3  2   3  2  3   4  3  3  2   2  2   1   3    2   3   2           2        2        2         2         3       2         1         2
            stripped_rows     = [".".join(x.split(".")[:-1]) for x in sdf.index]
            remove_these_rows = [i for i in sdf.index
                                 if (i.split("-")[-1].startswith("GCA"))
                                 and (".".join(i.replace("-GCA", "-GCF").split(".")[:-1]) in stripped_rows)]
            sdf = sdf.drop(index = remove_these_rows)
            # Now, there is an issue where there may be many poor-quality assemblies.
            # If there are both RefSeq (GCF) and GenBank (GCA) versions of the genome, then we will only consider the RefSeq version.
            #
            # For example, look at all of these pig assemblies.
            #                                A1a  D  F  C1  G  H  K  Ea  N  L  M  B1  I  O1  A1b  Eb  O2  Qa  (C1, M)  (G, H)  (Ea, Eb)  (O1, O2)
            # Susscrofa-9823-GCA031225015.1    2  3  2   2  2  3  2   3  3  2  1   2  1   2    2   1   2   1        1       1         1         2
            # Susscrofa-9823-GCF000003025.6<   2  3  2   3  3  3  2   3  3  3  3   2  1   2    2   2   2   2        2       2         2         2
            # Susscrofa-9823-GCA031306245.1    2  3  2   2  2  3  2   2  3  2  1   2  1   2    2   1   2   1        1       1         1         2
            # Susscrofa-9823-GCA002844635.1    2  3  2   2  2  3  2   3  2  3  1   2  1   2    2   1   2   1        1       1         1         2
            # Susscrofa-9823-GCA023065335.1    0  0  0   0  0  0  0   0  0  0  0   0  0   0    0   0   0   0        0       0         0         0
            # Susscrofa-9823-GCA015776825.1    2  3  2   2  2  3  2   2  3  2  1   3  2   2    2   1   2   1        1       1         1         2
            # Susscrofa-9823-GCA023065355.1    0  0  0   0  0  0  0   0  0  0  0   0  0   0    0   0   0   0        0       0         0         0
            # Susscrofa-9823-GCA030704935.1    2  3  2   2  2  3  2   2  3  2  1   2  1   2    2   1   2   2        1       1         1         2
            # Susscrofa-9823-GCA007644095.1    2  3  2   2  2  3  1   2  2  2  1   1  1   2    2   1   2   1        1       1         1         2
            # Susscrofa-9823-GCA024718415.1    2  3  2   2  2  3  2   3  3  2  1   1  1   2    2   1   2   1        1       1         1         2
            # Susscrofa-9823-GCA900119615.2    0  0  0   0  0  0  0   0  0  0  0   0  0   0    0   0   0   0        0       0         0         0

            # If there are any rows that are GCF, then keep only those rows.
            if any([x.split("-")[-1].startswith("GCF") for x in sdf.index]):
                keep_rows = [x for x in sdf.index if x.split("-")[-1].startswith("GCF")]
                sdf = sdf.loc[keep_rows]
            # if there was a value in sdf > 1, print
        else:
            # We currently have no special rules for internal nodes. We made most of the inferences based on the leaves.
            pass
        return sdf

    def _parental_probability_log(self, count,
                                  prob_eP_xAeB,
                                  prob_eP_xAxB):
        """
        Returns the logarithm of the probability.
        """
        if count == 0:
            return np.log(prob_eP_xAxB)
        if count == 1:
            return np.log(prob_eP_xAeB)
        else:
            base_probability_log = np.log(1 - prob_eP_xAeB)
            return np.log1p(-np.exp(count * base_probability_log))

    def _count_values_ge_1(self, column):
        """
        Counts the number of values greater than or equal to 1 in a column.
        """
        return (column >= 1).sum()

    def _determine_parental_ALG_PresAbs(self, pdf, sdf,
                                        prob_eP_xAeB = 0.99999,
                                        prob_eP_xAxB = 0.00001):
        """
        This is a helper method for self.percolate().
          - It uses only the tempdf.
          - It modifies a dataframe provided for the parent.
          - the default probabilities are coded into the parameters

        The input parameters are:
          - pdf - the parental df that we will be modifying
          - sdf - the dataframe of the sibling nodes
          - prob_eP_xAeB = 0.99999
          - prob_eP_eAeB = 0.99999
          - prob_eP_xAxB = 0.00001
        """
        # for each of the ALGs, check the condition and update the pDF
        # just get the sum of the columns
        ALGtemp = sdf[self.ALGcols]
        probabilities_log = ALGtemp.apply(lambda col: self._parental_probability_log(
                                          self._count_values_ge_1(col), prob_eP_xAeB, prob_eP_xAxB))
        # generate random numbers between 0 and 1. If the value is less than the probability, then we set the value to 1.
        #  Otherwise, we set the value to 0. Check it in log space.
        results = probabilities_log.apply(lambda x: 1 if np.log(np.random.random()) < x else 0)
        # print the results sideways, so we can see the results of the random number generation
        # The pdf has the same colnames as results. Use the results dataframe to update these values in the pdf.
        pdf[results.index] = results
        return pdf

    def add_taxname_to_all_nodes(self):
        """
        This function adds the taxname to a single node. Uses ete3.
        """
        # use ete3 to get the names of the taxids
        ncbi = NCBITaxa()
        for node in self.G.nodes():
            self.G.nodes[node]["taxname"] = ncbi.get_taxid_translator([node])[node].replace(" ", "-")

    def add_lineage_string_sample(self, lineage_string, samplename) -> int:
        """
        The lineage strings look like this:
          - 1;131567;2759;33154;33208;6040;6042;1779146;1779162;6060;6061;1937956;6063

        Then, there is a sample name too for the last node.

        Notes:
          - The edges from this string will be (1, 131567), (131567, 2759), (2759, 33154), etc.
        """
        fields = [int(x) for x in lineage_string.split(";")]
        for i in range(len(fields)-1):
            self.G.add_edge(fields[i], fields[i+1])
        # add the final edge
        self.G.add_edge(fields[-1], samplename)
        return 0

    def _build_tree_from_perspchrom(self) -> int:
        """
        This function takes in a per_sp_chrom_df and builds a tree from it.
        """
        # add each lineage string to the tree
        for i, row in self.perspchrom.iterrows():
            self.add_lineage_string_sample(row["taxidstring"], row["species"])
        return 0

    def _get_edges_in_clade_helper(self, node):
        """
        This is the recursive case for the get_edges_in_clade function.
        """
        # get the outgoing edges from this node.
        out_edges = list(self.G.out_edges(node))
        # recursive break condition - if there are no outgoing edges, then return an empty list
        if len(out_edges) == 0:
            return []
        out_nodes = [x[1] for x in out_edges]
        for thisnode in out_nodes:
            out_edges += self._get_edges_in_clade_helper(thisnode)
        return out_edges

    def get_edges_in_clade(self, node) -> list:
        """
        This function takes in a node ID (clade and returns a recursive list of all
          the outgoing edges, and the single incoming edge.
        """
        if not isinstance(node, int):
            node = int(node)

        # get the single incoming edge. Make sure it is a tuple
        in_edges = list(self.G.in_edges(node))
        if len(in_edges) > 1:
            raise Exception("There should only be one incoming edge. We don't allow reticulate phylogenetic trees. Found {}".format(in_edges))

        return in_edges + self._get_edges_in_clade_helper(node)

def main():
    # parse the args
    args = parse_args()

    # get the rbh files in the directory
    rbh_files = list(sorted([os.path.join(args.directory, f)
                 for f in os.listdir(args.directory)
                 if f.endswith('.rbh')], reverse = True))
    #rbh_files = rbh_files[:500]

    # There are two files that we only need to calculate once.
    # The locdf is the dataframe that contains the location of the ALGs on the chromosomes.
    # The perspchrom is the dataframe that contains the presence/absence of the ALGs in the species.
    # The perspchrom is a derivative of the locdf, and takes a while to calculate, so it is better to just
    #  do it once and save it to a file.
    locdf      = None
    perspchrom = None
    overwrite = False
    if os.path.exists("locdf.tsv") and os.path.exists("perspchrom.tsv"):
        # These files both exist, so we can just read them in.
        print("Reading in the locdf and perspchrom from file.", file = sys.stderr)
        locdf      = pd.read_csv("locdf.tsv", sep = "\t")
        perspchrom = pd.read_csv("perspchrom.tsv", sep = "\t")
        # if there is a '(' or ')' in the column names, then we need to convert them to tuples
        perspchrom.columns = [tuple(eval(x)) if isinstance(x, str) and "(" in x else x for x in perspchrom.columns]
    else:
        # The files do not yet exist, so we need to calculate them
        print("Calculating the locdf and perspchrom df from the input files.", file = sys.stderr)
        locdf, perspchrom = rbh_files_to_locdf_and_perspchrom(rbh_files, args.ALG_rbh,
                                                              args.minsig, args.ALGname)
        # save the locdf and perspchrom to a file
        locdf.to_csv("locdf.tsv", sep = "\t", index = False)
        perspchrom.to_csv("perspchrom.tsv", sep = "\t", index = False)
        overwrite = True

    ## filters for testing just get the chordates. Check if ";7711;" is in the taxidstring
    #perspchrom = perspchrom[perspchrom['taxidstring'].str.contains(';7711;')]

    resultstsv = "tree1.tsv.gz"
    resultsdf = None
    if overwrite:
        # now that we have the perspchrom, we should construct the tree structure
        T = SplitLossColocTree(perspchrom)
        print("Percolating", file = sys.stderr)
        T.percolate()
        print("Saving the tree to tsv", file = sys.stderr)
        T.save_tree_to_df(resultstsv)
        resultsdf = T.tree_df
    elif not overwrite:
        # if we don't want to overwrite the results (tree1.tsv)
        # then we just read it in as a pdf
        if os.path.exists(resultstsv):
            resultsdf = pd.read_csv(resultstsv, sep = "\t", index_col = 0)
        else:
            raise IOError(f"The file {resultstsv} does not exist. We need to calculate the tree.")
    #
    print("Making a UMAP with matplotlib", file = sys.stderr)
    #save_UMAP(resultsdf)
    print("Making a UMAP with plotly", file = sys.stderr)
    save_UMAP_plotly(resultsdf, "tree1")

    sys.exit()
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
        # The format is [taxid, "(gain|loss|split)", "taxid", "(gain|loss|split)", ...]
        #  We add new entries as we go through the tree to the root.
        #  At the end, we flip the order and make a parsable string
        changeString = []
        # we have access to ["species", "taxid", "taxidstring"]
        thistaxidstring = row["taxidstring"]

        # --------- CONSERVED ----------- #
        # note the ALGs that have a value of 1 here. These are conserved in this species and will never be in changeString
        ALGsConservedInThisSample = [x for x in ALG_columns if row[x] >= 1]
        # Note the things specifically missing in this species. These will be added to the changeString at some point

        # ---------- MISSING ------------ #
        #  This list will not be changed at all during execution of this loop. Only the next ALGsMissingInThisSampleAccountedFor will be changed.
        ALGsMissingInThisSample   = [x for x in ALG_columns if row[x] == 0]
        # We will add to this list as we go through the tree. This will be the list of ALGs that were lost on various branches.
        ALGsMissingInThisSampleAccountedFor = [] # At the end, all of the missing ALGs must appear somewhere in ALGsMissingInThisSample

        # ----------- SPLITS ------------ #
        ALGsSplitInThisSample = [x for x in ALG_columns if row[x] > 1]
        ALGsSplitInThisSampleAccountedFor = []

        # ----- ALG COLOCALIZATION ------ #
        # For each species, we similarly need to keep track of which ALG combos are gained on each branch.
        ColocalizedPairsInThisSample             = [x for x in ALG_combos if row[x] == 1]
        ColocalizedPairsInThisSampleAccountedFor = []

        ## DEBUG
        #print()
        #print("The sample is {}".format(row["species"]))
        #print("ALGsConservedInThisSample: {}".format(ALGsConservedInThisSample))
        #print("ALGsMissingInThisSample: {}".format(ALGsMissingInThisSample))
        #print("ALGsSplitInThisSample: {}".format(ALGsSplitInThisSample))
        #print("num ColocalizedPairsInThisSample: {}".format(len(ColocalizedPairsInThisSample)))
        #print()

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
            ALG_splits_on_this_branch = []

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

                # ----- ALG SPLITS -----
                ALG_splits_on_this_branch = [x for x in ALGsSplitInThisSample if x not in ALGsMissingInThisSampleAccountedFor]
                # add all the new ALG splits to the ALGsMissingInThisSampleAccountedFor
                ALGsSplitInThisSampleAccountedFor.extend(ALG_splits_on_this_branch)
            else:
                # Because we're not at the last node we need to do some more work
                # Get the subdf of the perspchrom df for rows with a taxidstring that contains the new thistaxidstring, but not the prevtaxidstring
                #  Excluding the prevtaxidstring is important, because it compares this clade only to sister clades at the same level,
                #  and removes the possibility of paraphyletic comparisons.
                subdf = perspchrom[perspchrom["taxidstring"].str.contains(thistaxidstring) &
                                   ~perspchrom["taxidstring"].str.contains(prevtaxidstring)]
                if len(subdf) == 0:
                    # If there is nothing to compare to in this clade, then we don't need to do anything and can proceed to printing the changeString.
                    #   - 20240207: The fundamental reason that we can't do anything here is that we don't have any information about the ALGs in this clade.
                    # Deprecated:
                    #   - 20240207: In the future there could be an option to put ambiguous changes of gains or losses here, but for now
                    #     we are just going to put the changes on the oldest nodes that we can find.
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

                    # ----- ALG SPLITS -----
                    notSplitInThisClade = unsplit_ALGs(subdf,
                                                       max_frac_split = min_for_noncolocalized)
                    ALGsplitsThisBranch = [x for x in notSplitInThisClade
                                           if x in ALGsSplitInThisSample
                                           and x not in ALGsSplitInThisSampleAccountedFor]
                    ALG_splits_on_this_branch.extend(ALGsplitsThisBranch)
                #print()
            thisChange = "({}|{}|{})".format(
                sorted(ALG_colocalizations_on_this_branch),
                sorted(ALG_losses_on_this_branch),
                sorted(ALG_splits_on_this_branch))
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
        # we should do this for splits
        if sorted(ALGsSplitInThisSampleAccountedFor) != sorted(ALGsSplitInThisSample):
            raise IOError("There is a discrepancy between the ALGsSplitInThisSampleAccountedFor and ALGsSplitInThisSample. Write more debugging code to figure out what the issue is, because I haven't worked on this yet.")
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

    # all of this is for doing the clade-level analysis
    ## ---------------------------------------------------------------------------------------------
    ##  Move onto the per-node analysis
    ## ---------------------------------------------------------------------------------------------
    ## Get the labels. They have to be in the order that the leaves are returned
    #leaves = [int(str(x).split("-")[-1]) for x in tree.get_leaves()]
    ## Make a dict with the taxid and species cols, then make a label from the lookup with leaves
    #lookup = dict(zip(perspchrom["taxid"], perspchrom["species"]))
    #labels = [lookup[x] for x in leaves]
    #for leaf, label in zip(tree.get_leaves(), labels):
    #    leaf.name = f"{label}_{leaf.name}"
    #tree.write(format=1, outfile="species_tree.tre")

    ## Now annotate all of the nodes.
    ## Yes, this loops through the table again, but I don't have a more elegant solution
    ## Right now this doesn't actually annotate any nodes, it just makes a dictionary of the annotations
    #node_annotations = {}
    ## Make a dict of annotations for each node
    #for i, row in perspchrom.iterrows():
    #    spstring = row["species"]
    #    taxidstring = [int(x) for x in row["taxidstring"].split(";")]
    #    for thisid in  taxidstring:
    #        if int(thisid) not in node_annotations:
    #            node_annotations[int(thisid)] = set()
    #        node_annotations[int(thisid)].add(spstring)

    ## we now have a dictionary with which species belong in each node
    ## iterate through all of the nodes of the tree, recursively. ACTUALLY IT ISN'T DOING THAT NOW
    #entries = []
    #for taxid in node_annotations:
    #    #print(taxid, node_annotations[taxid])
    #    thistaxid = int(taxid)
    #    thisnodename = ncbi.get_taxid_translator([thistaxid]).get(thistaxid, "Unknown")
    #    # get the NCBI taxid lineage for this node
    #    thislineage = ";".join([str(x) for x in ncbi.get_lineage(thistaxid)])
    #    # get the set of species that belong to this node
    #    these_species = node_annotations[thistaxid]
    #    # get a sub table of the perspchrom table that only has these species in the species column
    #    subdf = perspchrom[perspchrom["species"].isin(these_species)]
    #    # sum up the dataframe to get the number of fusions for each ALG, get rid of all the other columns
    #    subdf = subdf.drop(columns=["species", "taxid", "taxidstring"])
    #    subdf = subdf.sum(axis=0)
    #    # now add the other information as new columns, thistaxid/thisnodename/thislineage
    #    subdf["taxid"] = thistaxid
    #    subdf["nodename"] = thisnodename
    #    subdf["taxidstring"] = thislineage
    #    subdf["spinthisclade"] = ",".join(these_species)
    #    entries.append(subdf.copy())
    ## condense all of the entries into a single df
    #per_node_df = pd.DataFrame(entries)
    ## move the taxid, nodename, thislineage columns to the front
    #per_node_df.insert(0, "spinthisclade", per_node_df.pop("spinthisclade") )
    #per_node_df.insert(0, "taxidstring",   per_node_df.pop("taxidstring")   )
    #per_node_df.insert(0, "nodename",      per_node_df.pop("nodename")      )
    #per_node_df.insert(0, "taxid",         per_node_df.pop("taxid")         )
    ## save this!
    #per_node_df.to_csv("all_nodes_ALG_presence_fusions.tsv", sep='\t', index=False)

    ## make figures of the per-species plots
    #standard_plot_out(per_node_df, "perNode")
    ## Nodes we want to compare:
    ##  10197 - Ctenophora
    ##  6040 - Sponges
    ##  10226 - Placozoa
    ##  6073 - Cnidaria
    ##  6231 - Nematodes
    ##  88770 - Panarthropoda
    ##  2697495 - Spiralia
    ##  7711 - Chordata
    ##  7586 - Echinodermata
    ##  10219 - Hemichordata
    ## pull out a df of just these nodes. Use an exact match of the taxid column
    #nodes_to_compare = [10197, 6040, 10226, 6073, 6231, 88770, 2697495, 7711, 7586, 10219]
    #comparison_set_df = per_node_df[per_node_df["taxid"].isin(nodes_to_compare)]
    ## Sort the comparison set by the taxidstring of the nodes_to_compare list.
    ## Use the order of numbers in nodes_to_compare to sort the dataframe.
    #comparison_set_df = comparison_set_df.sort_values(by=["taxidstring"], ascending=True,
    #                                                  key=lambda x: x.map(dict(zip(nodes_to_compare, range(len(nodes_to_compare))))))
    #print("The species in the node comparison are:")
    #print(comparison_set_df)
    #standard_plot_out(comparison_set_df, "comparisonNodes")

    ## print the tree to a .tre file
    #tree = ncbi.get_topology([int(x) for x in perspchrom["taxid"].tolist()])

    ##taxid_to_query = 2301116
    ##species_under_taxid = get_species_under_taxid(taxid_to_query)
    ##print(f"Species under taxid {taxid_to_query}: {species_under_taxid}")

    ### remove rows where all values are 0
    ##df = df.loc[(df!=0).any(axis=1)]
    ##pca = PCA(n_components=2)
    ##pca.fit(df)
    ##print(pca.components_)
    ##print(pca.explained_variance_)

    ##df2 = pd.DataFrame(pca.transform(df), columns = ['first', 'second'])
    ##print(df)
    ##print(df2)
    ### df2.plot.scatter(x = 'first', y = 'second')

    ###plt.show()

    ##from mpl_toolkits import mplot3d
    ##from mpl_toolkits.mplot3d import Axes3D
    ##from mpl_toolkits.mplot3d import proj3d
    ##from matplotlib.text import Annotation

    ##x = df2['first']
    ##y = df2['second']
    ### labels is a list of the df1 index values
    ##labels = df.index.values.tolist()

    ### Create the scatter plot
    ##fig, ax = plt.subplots()
    ##scatter = ax.scatter(x, y, picker=True)
    ### plot the text labels
    ##for i, txt in enumerate(labels):
    ##    ax.annotate(txt, (x[i], y[i]))


    ### Show the plot
    ##plt.show()


if __name__ == '__main__':
    main()