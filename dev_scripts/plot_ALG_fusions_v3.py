#!/usr/bin/env python

"""
Program  : plot_ALG_fusions_v3.py
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

DESCRIPTION:
============
This program analyzes ancestral linkage group (ALG) fusions, losses, and splits across species
by constructing phylogenetic event strings for each species. It tracks chromosomal rearrangements
across evolutionary time and generates interactive visualizations.

This version (v3) is an updated version of plot_ALG_fusions.py that:
- Constructs detailed event strings for each species showing gains/losses/splits
- Maps ALG colocalizations onto phylogenetic trees
- Generates UMAP visualizations of chromosomal architecture patterns
- Creates per-species and per-node analyses of chromosomal changes

KEY FEATURES:
=============
1. ALG Colocalization Analysis: Identifies which ALGs are found on the same chromosome
2. Phylogenetic Event Mapping: Tracks when fusions, losses, and splits occurred
3. Interactive Visualizations: Creates UMAP plots showing species clustering by chromosomal patterns
4. Automated Caching: Saves intermediate results to speed up re-runs

INPUT FILES:
============
1. Directory of RBH files (-d, --directory):
   - One .rbh file per species comparison against ALGs
   - Each file contains reciprocal best hits with columns:
     * {ALGname}_scaf, {ALGname}_gene, {ALGname}_pos
     * {species}_scaf, {species}_gene, {species}_pos
     * whole_FET: Fisher's Exact Test p-value for colocalization
     * rbh, gene_group, color

2. ALG RBH reference file (-r, --ALG_rbh):
   - Master RBH file defining ALG properties (names, sizes, colors)
   - Used as the reference database for all ALGs
   - Typically named something like "BCnS.rbh" or similar

3. ALG database name (-a, --ALGname):
   - Name of the ALG database (e.g., "BCnSSimakov2022", "BCnS")
   - Must match the column prefix in the RBH files

4. Calibrated tree file (-t, --tree_info) [OPTIONAL]:
   - Path to node_information.tsv from Newick_to_common_ancestors.py
   - Contains custom topology (e.g., Ctenophora placement) and calibrated divergence times
   - If not provided, falls back to NCBI taxonomy with apply_custom_phylogeny()
   - Recommended for consistent phylogenetic framework across analyses

REQUIRED PARAMETERS:
====================
-d, --directory    : Directory containing species .rbh files
-a, --ALGname      : Name of ALG database (e.g., "BCnSSimakov2022")
-r, --ALG_rbh      : Path to master ALG RBH file

OPTIONAL PARAMETERS:
====================
-t, --tree_info    : Path to node_information.tsv from Newick_to_common_ancestors.py
                     Uses pre-built calibrated tree with custom topology and TimeTree ages
                     If not provided, builds tree from NCBI + custom phylogeny
-m, --minsig       : Minimum significance for whole_FET (default: 0.005)
                     Lower values = more stringent colocalization requirement

OUTPUT FILES:
=============
1. locdf.tsv
   - Location dataframe showing which ALGs are on which scaffolds
   - Columns: sample, gene_group, scaffold, pvalue, num_genes, frac_of_this_ALG_on_this_scaffold

2. perspchrom.tsv
   - "Perspective chromosome" dataframe with presence/absence and colocalization
   - One row per species
   - Columns for each ALG (binary: present=1, absent=0, split>1)
   - Columns for each ALG pair (binary: colocalized=1, separate=0)
   - Contains 'changestrings' encoding evolutionary events

3. per_species_ALG_presence_fusions.tsv
   - Extended perspchrom with detailed changestrings
   - Changestring format: taxid-([coloc_gains]|[losses]|[splits])-taxid-...
   - Example: "1-([]|[]|[])-131567-([('A1b','B3')]|[]|[])-2759-..."

4. tree1.tsv.gz
   - Phylogenetic tree structure with all evolutionary events
   - Used for downstream plotting and analysis

5. tree1_umap.html
   - Interactive UMAP visualization (Plotly)
   - Shows species clustering based on chromosomal architecture
   - Color-coded by taxonomy

CHANGESTRING FORMAT:
====================
The changestring encodes events on branches leading TO each node:
Format: taxid-([colocalizations]|[losses]|[splits])-taxid-(...)
- [colocalizations]: List of ALG pairs that fused on this branch
- [losses]: List of ALGs that were lost on this branch
- [splits]: List of ALGs that split on this branch

Example:
"1-([]|[]|[])-131567-([]|[]|[])-2759-([('A','B')]|['C']|[])-33208-..."
This shows:
- At taxid 2759: ALGs A and B fused, ALG C was lost
- No events at taxids 1, 131567

ALGORITHM OVERVIEW:
===================
1. Parse all RBH files to identify ALG-chromosome associations
2. Build presence/absence matrix (perspchrom) for all species
3. For each species, walk up phylogenetic tree (NCBI taxonomy)
4. At each node, compare to sister clades to infer:
   - New fusions (ALG pairs colocalized here but not in sister clades)
   - Losses (ALGs present in sister but absent here)
   - Splits (ALGs on 1 chromosome in sister, >1 here)
5. Construct phylogenetic tree object with all events
6. Generate UMAP projection for visualization

THRESHOLDS:
===========
- min_for_missing = 0.8
  An ALG is considered "missing" in a clade if absent in 80%+ of species
  
- min_for_noncolocalized = 0.5
  ALG pairs are "separated" if colocalized in <50% of sister clade

USAGE EXAMPLES:
===============
Basic usage:
    python plot_ALG_fusions_v3.py \\
        -d /path/to/rbh_files/ \\
        -a BCnSSimakov2022 \\
        -r /path/to/BCnS.rbh

With calibrated tree (recommended):
    python plot_ALG_fusions_v3.py \\
        -d /path/to/rbh_files/ \\
        -a BCnSSimakov2022 \\
        -r /path/to/BCnS.rbh \\
        -t output_prefix.node_information.tsv

With custom significance threshold:
    python plot_ALG_fusions_v3.py \\
        -d ./rbh_output/ \\
        -a BCnSSimakov2022 \\
        -r BCnS.rbh \\
        -t calibrated_tree.node_information.tsv \\
        -m 0.001

PERFORMANCE NOTES:
==================
- First run calculates locdf and perspchrom (slow, ~minutes to hours depending on dataset size)
- Subsequent runs use cached files (fast, seconds)
- Delete locdf.tsv and perspchrom.tsv to force recalculation
- Tree construction is also cached in tree1.tsv.gz

PREREQUISITES:
==============
Required Python packages:
- pandas, numpy
- ete4 (phylogenetic trees)
- matplotlib
- umap-learn
- plotly (for interactive plots)
- scikit-learn
- networkx
- PIL

NCBI Taxonomy Database:
You must initialize the NCBI taxonomy database before first use:
    ```python
    from ete4 import NCBITaxa
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    ```
    
See: http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html

NOTES:
======
- The script automatically caches intermediate results to speed up re-runs
- ALG inference node is hardcoded as "1;131567;2759;33154;33208" (Metazoa)
- Modify this value in the code if analyzing different taxonomic groups
- The script uses Fisher's Exact Test to determine significant colocalizations
- Connected components are used to group multi-way fusions (AxBxC counts as 1 event)

AUTHORS: Darrin T. Schultz
VERSION: 3.0
DATE: 2024
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
from multiprocessing import Pool
import threading
#import warnings
#warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

from PIL import Image
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster import hierarchy

# import odp-specific functions
thisfile_path = os.path.dirname(os.path.realpath(__file__))
source_path = os.path.join(thisfile_path, "../source")
sys.path.insert(1, source_path)
import rbh_tools

from taxid_tools import NCBI_taxid_to_taxdict

# import the stuff to work with lineages
from ete4 import NCBITaxa,Tree

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
      - ALG_rbh: the name of the rbh file that contains the ALG information. This is used to get the ALG names, sizes, and colors.
      - tree_info: (optional) path to node_information.tsv from Newick_to_common_ancestors.py with calibrated tree
      - minsig: the minimum significance value for the whole_FET column in the rbh files. This is used to filter the rbh files.
    """
    parser = argparse.ArgumentParser(description='Plot the size of ALGs against the number of colocalized ALGs')
    parser.add_argument('-d', '--directory', type=str,   required=True,   help='The directory where the ALG files are saved')
    parser.add_argument('-a', '--ALGname',   type=str,   required=True,   help='The name of the ALG database')
    parser.add_argument('-r', '--ALG_rbh',   type=str,   required=True,   help='The name of the rbh file that contains the ALG information. This is used to get the ALG names, sizes, and colors.')
    parser.add_argument('-t', '--tree_info', type=str,   required=False,  help='Path to node_information.tsv from Newick_to_common_ancestors.py with calibrated tree topology. If not provided, falls back to NCBI taxonomy with custom phylogeny.')
    parser.add_argument('-m', '--minsig',    type=float, default = 0.005, help='The minimum significance value for the whole_FET column in the rbh files. This is used to filter the rbh files.')
    parser.add_argument('--parallel', action='store_true', help='Use multiprocessing for changestring generation (recommended for large datasets)')
    parser.add_argument('--ncores', type=int, default=20, help='Number of CPU cores to use for parallel processing (default: 20)')

    args = parser.parse_args()
    # Check if tree_info file exists if provided
    if args.tree_info and not os.path.exists(args.tree_info):
        raise ValueError(f'The tree info file {args.tree_info} does not exist')
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

def apply_custom_phylogeny(lineage, taxid, ncbi):
    """
    Apply custom phylogenetic placement that differs from NCBI taxonomy.
    
    CUSTOM TOPOLOGY:
    ================
    This function implements the following alternative animal phylogeny:
    
    Metazoa (33208)
    ├─ Myriazoa (-67) [CUSTOM NODE - not in NCBI]
    │  ├─ Porifera (6040)
    │  └─ Eumetazoa (6072)
    │     ├─ Cnidaria (6073)
    │     └─ [other Eumetazoa]
    └─ Ctenophora (10197) [MOVED - sister to Myriazoa]
    
    In NCBI, Ctenophora is nested within Eumetazoa (6072), but phylogenomic
    evidence suggests they are sister to all other animals (Schultz et al. 2023).
    
    Parameters:
    -----------
    lineage : list
        Original NCBI lineage as list of taxids
    taxid : int
        The taxid being processed
    ncbi : NCBITaxa
        NCBI taxonomy database object
        
    Returns:
    --------
    list : Modified lineage with custom phylogenetic placement
    
    Notes:
    ------
    - Myriazoa (-67) is a fake taxid representing Porifera + Eumetazoa
    - Negative taxids will not conflict with real NCBI taxids (all positive)
    - The function checks if taxid is Ctenophora or descendant thereof
    - If not Ctenophora-related, checks if taxid is Porifera/Eumetazoa/descendants
    """
    CTENOPHORA_TAXID = 10197
    EUMETAZOA_TAXID = 6072
    PORIFERA_TAXID = 6040
    METAZOA_TAXID = 33208
    MYRIAZOA_TAXID = -67  # Custom fake taxid for Porifera + Eumetazoa clade
    
    # Check if this taxid is Ctenophora or a descendant
    if CTENOPHORA_TAXID in lineage:
        # Find where Metazoa is in the lineage
        if METAZOA_TAXID in lineage:
            metazoa_index = lineage.index(METAZOA_TAXID)
            # Build custom lineage: [root...Metazoa, Ctenophora, ...descendants]
            # Remove any intermediate nodes between Metazoa and Ctenophora
            new_lineage = lineage[:metazoa_index + 1]  # Everything up to and including Metazoa
            
            # Find the position of Ctenophora in original lineage
            cteno_index = lineage.index(CTENOPHORA_TAXID)
            # Add everything from Ctenophora onwards
            new_lineage.extend(lineage[cteno_index:])
            return new_lineage
    
    # Check if this taxid is Porifera, Eumetazoa, or descendants thereof
    elif PORIFERA_TAXID in lineage or EUMETAZOA_TAXID in lineage:
        # Find where Metazoa is in the lineage
        if METAZOA_TAXID in lineage:
            metazoa_index = lineage.index(METAZOA_TAXID)
            # Build custom lineage: [root...Metazoa, Myriazoa, Porifera/Eumetazoa, ...descendants]
            new_lineage = lineage[:metazoa_index + 1]  # Everything up to and including Metazoa
            new_lineage.append(MYRIAZOA_TAXID)  # Insert fake Myriazoa node
            
            # Add the rest of the lineage after Metazoa
            new_lineage.extend(lineage[metazoa_index + 1:])
            return new_lineage
    
    # For all other taxa, return original lineage unchanged
    return lineage


def taxids_to_taxidstringdict(taxids, use_custom_phylogeny=True) -> dict:
    """
    This function takes a list of taxids and returns a dictionary where the key is the species name and the value is the taxid string.
    The taxid string will be a string of taxids separated by a semicolon.
    
    Parameters:
    -----------
    taxids : list/set/dict
        Iterable of taxids to process
    use_custom_phylogeny : bool, default=True
        If True, applies custom phylogenetic corrections (e.g., Ctenophora placement)
        If False, uses NCBI taxonomy exactly as-is
        
    Notes:
    ------
    When use_custom_phylogeny=True:
    - Ctenophora (10197) is placed as sister to Myriazoa (Porifera + Eumetazoa)
    - Myriazoa is represented by fake taxid -67
    - This reflects the phylogenomic hypothesis from Schultz et al. (2023) Nature
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
        
        # Apply custom phylogeny if requested
        if use_custom_phylogeny:
            lineage = apply_custom_phylogeny(lineage, taxid, ncbi)
        
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

def standard_plot_out(perspchrom, outprefix, taxid_order = None, safe = False)->None:
    """
    Makes a pixel-wise plot of a dataframe.
    Currently has a "phylogenetic" bar on the left side.
    Then there is an accounting of which ALGs are present in each species.
    Then there is a plot of the ALG colocalization matrix.

    At the top there is a colorbar that encodes the ALG colocalization pairs.
    """
    # first, convert the colnames that have "(" as the first character to a tuple. just use eval
    for i in range(len(perspchrom.columns)):
        col = perspchrom.columns[i]
        if type(col) == str:
            if col[0] == "(" and col[-1] == ")":
                perspchrom.rename(columns={col: eval(col)}, inplace=True)

    concatenate_these_images = []
    # convert matrix to a species image
    # make an image of a pseudo-phylogenetic tree
    if taxid_order is None:
        # In this case we use whatever sorting the dataframe has
        tree_image    = image_sp_matrix_to_lineage(perspchrom["taxidstring"])
        bar1 = image_vertical_barrier(2, len(perspchrom), color = "#000000")
        bar2 = image_vertical_barrier(1, len(perspchrom), color = "#F0F93E")
        bar3 = image_vertical_barrier(2, len(perspchrom), color = "#000000")
        concatenate_these_images += [tree_image, bar1, bar2, bar3]
    else:
        print(taxid_order)

        # In this case, get the first taxid that matches the one in the taxid_order
        # Then sort the dataframe by this order
        keep_indices = []
        for taxid in taxid_order:
            # get the first taxid that matches the one in the taxid_order
            if len(perspchrom[perspchrom["taxid"] == taxid]) == 0:
                if safe:
                    raise ValueError("The taxid {} is not in the dataframe.".format(taxid))
            else:
                keep_indices.append(perspchrom[perspchrom["taxid"] == taxid].index[0])
        print(keep_indices)
        # now it is sorted by the taxid_order
        perspchrom = perspchrom.loc[keep_indices]
        perspchrom = perspchrom.reset_index(drop=True)

    # make an image of the ALG presence/absence matrix
    presabs_image = image_sp_matrix_to_presence_absence(perspchrom, color_dict = dict_BCnSALG_to_color())
    # gap bars
    bar4 = image_vertical_barrier(2, len(perspchrom), color = "#000000")
    bar5 = image_vertical_barrier(1, len(perspchrom), color = "#F0F93E")
    bar6 = image_vertical_barrier(2, len(perspchrom), color = "#000000")
    concatenate_these_images += [presabs_image, bar4, bar5, bar6]

    # make an image of the colocalization matrix
    coloc_image = image_colocalization_matrix(perspchrom, clustering = True, color_dict = dict_BCnSALG_to_color())
    coloc_image_unclust = image_colocalization_matrix(perspchrom, clustering = False, color_dict = dict_BCnSALG_to_color())

    # concatenate all the images
    concatenate_clustered = concatenate_these_images + [coloc_image]
    composite_image = image_concatenate_horizontally(concatenate_clustered)
    composite_image.save("{}_composite_image.png".format(outprefix))

    # make the same plot, but not sorted
    concatenate_unclustered = concatenate_these_images + [coloc_image_unclust]
    composite_image = image_concatenate_horizontally(concatenate_unclustered)
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

def load_calibrated_tree(node_info_file):
    """
    Load calibrated phylogenetic tree from Newick_to_common_ancestors.py output.
    
    This function reads the node_information.tsv file which contains the custom topology
    and calibrated divergence times from TimeTree.org.
    
    Parameters:
    -----------
    node_info_file : str
        Path to {prefix}.node_information.tsv from Newick_to_common_ancestors.py
        
    Returns:
    --------
    dict : Dictionary with two keys:
        'lineages' : dict mapping taxid (int) -> lineage_string (str)
            Where lineage_string is semicolon-delimited taxid path from root to tip
            Example: "1;131567;2759;33154;33208;-67;6072;33213;7711;9606"
        'ages' : dict mapping taxid (int) -> nodeage (float)
            Where nodeage is the divergence time in millions of years ago (MYA)
            Only includes internal nodes that have calibrated ages
        
    Notes:
    ------
    - Preserves custom topology (e.g., Ctenophora placement, Myriazoa node)
    - Uses lineage_string column which has the full taxid path
    - Only includes leaf nodes (species) in lineages, internal nodes in ages
    - Ages are from TimeTree.org molecular clock calibration
    """
    print(f"Loading calibrated tree from: {node_info_file}", file=sys.stderr)
    
    # Read the node information file
    node_df = pd.read_csv(node_info_file, sep="\t")
    
    # Filter to only leaf nodes (species) - these have lineage_string values
    # Internal nodes might have different structure
    leaf_nodes = node_df[node_df['lineage_string'].notna()].copy()
    
    # Build mapping of taxid -> lineage_string
    taxid_to_lineage = {}
    for idx, row in leaf_nodes.iterrows():
        taxid = int(row['taxid'])
        lineage_str = str(row['lineage_string'])
        taxid_to_lineage[taxid] = lineage_str
    
    # Build mapping of taxid -> nodeage for all nodes with calibrated ages
    # This includes internal nodes that have divergence times
    taxid_to_age = {}
    age_nodes = node_df[node_df['nodeage'].notna()].copy()
    for idx, row in age_nodes.iterrows():
        taxid = int(row['taxid'])
        nodeage = float(row['nodeage'])
        taxid_to_age[taxid] = nodeage
    
    print(f"  Loaded lineages for {len(taxid_to_lineage)} species from calibrated tree", file=sys.stderr)
    print(f"  Loaded node ages for {len(taxid_to_age)} nodes from calibrated tree", file=sys.stderr)
    
    return {'lineages': taxid_to_lineage, 'ages': taxid_to_age}

def rbh_files_to_locdf_and_perspchrom(rbh_files, ALGrbhfile, minsig, ALGname, calibrated_tree_data=None) -> (pd.DataFrame, pd.DataFrame):
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
    # Use calibrated tree if provided, otherwise fall back to NCBI + custom phylogeny
    sample_to_taxidstring = {}
    
    if calibrated_tree_data is not None:
        # Use pre-built lineages from calibrated tree
        # calibrated_tree_data is now a dict with 'lineages' and 'ages' keys
        lineages_dict = calibrated_tree_data['lineages']
        print("  Using lineages from calibrated tree", file=sys.stderr)
        missing_taxa = []
        for k in sample_to_taxid:
            taxid = sample_to_taxid[k]
            if taxid in lineages_dict:
                sample_to_taxidstring[k] = lineages_dict[taxid]
            else:
                missing_taxa.append((k, taxid))
        
        if missing_taxa:
            print(f"\n  WARNING: {len(missing_taxa)} taxa not found in calibrated tree", file=sys.stderr)
            print(f"  These samples will be EXCLUDED from analysis to avoid topology conflicts\n", file=sys.stderr)
            
            # Write full list to a file for later review
            missing_file = "missing_taxa_from_calibrated_tree.txt"
            with open(missing_file, 'w') as f:
                f.write(f"# {len(missing_taxa)} taxa found in RBH files but missing from calibrated tree\n")
                f.write(f"# These samples are EXCLUDED to avoid topology conflicts\n")
                f.write(f"# The calibrated tree has different custom nodes (-67, -68) than apply_custom_phylogeny()\n")
                f.write(f"#\n")
                f.write("# Format: sample_name\ttaxid\n")
                for k, taxid in sorted(missing_taxa, key=lambda x: x[0]):
                    f.write(f"{k}\t{taxid}\n")
            
            print(f"  Full list saved to: {missing_file}", file=sys.stderr)
            print(f"  First 20 missing samples:", file=sys.stderr)
            for k, taxid in sorted(missing_taxa, key=lambda x: x[0])[:20]:
                print(f"    - {k} (taxid: {taxid})", file=sys.stderr)
            if len(missing_taxa) > 20:
                print(f"    ... and {len(missing_taxa) - 20} more (see {missing_file})", file=sys.stderr)
            
            # Remove missing samples from sample_to_taxid dict
            # This prevents them from being added to perspchrom
            for k, taxid in missing_taxa:
                del sample_to_taxid[k]
            
            print(f"\n  Continuing with {len(sample_to_taxid)} samples that are in calibrated tree\n", file=sys.stderr)
            
            # Filter locdf to only include samples that have lineages
            excluded_samples = [k for k, taxid in missing_taxa]
            original_len = len(locdf)
            locdf = locdf[~locdf["sample"].isin(excluded_samples)].reset_index(drop=True)
            print(f"  Filtered locdf from {original_len} to {len(locdf)} entries\n", file=sys.stderr)
    else:
        # Fall back to NCBI taxonomy with custom phylogeny correction
        print("  Building lineages from NCBI taxonomy with custom phylogeny", file=sys.stderr)
        NCBI = NCBITaxa()
        for k in sample_to_taxid:
            taxid = sample_to_taxid[k]
            taxdict = NCBI_taxid_to_taxdict(NCBI, taxid)
            lineage = taxdict["taxid_list"]
            
            # Apply custom phylogeny correction (Ctenophora as sister to Myriazoa)
            lineage = apply_custom_phylogeny(lineage, taxid, NCBI)
            
            lineage_string = ";".join([str(x) for x in lineage])
            sample_to_taxidstring[k] = lineage_string

    # now make a pandas dataframe of the sample_to_taxidstring. The columns are "sample", "taxid"
    perspchrom = pd.DataFrame.from_dict(sample_to_taxid, orient='index')
    # change the index to a column. the former index is "species", the other column is "taxid"
    perspchrom = perspchrom.reset_index()
    perspchrom = perspchrom.rename(columns={"index": "species", 0: "taxid"})
    # the taxidstring is the sample_to_taxidstring dictionary. Map by species name, not taxid.
    perspchrom["taxidstring"] = perspchrom["species"].map(sample_to_taxidstring)
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
                    43845:   "#F93AD5", # Drosophilinae - these are the fruit flies
                    41084:   "#F1890D", # Polyphaga - these are the beetles
                    3042114: "#EBB836", # Anthophila - these are the bees
                    7147:    "#C02424", # Diptera
                      43741:  "#7373E6", # Acalyptratae - acalpytrate muscoid flies
                      43742:  "#314006", # Calyptratae - calyptrate muscoid flies
                      43838:  "#7EBB27", # Syrphidae-Syrphinae   (hoverflies)
                      115244: "#2654BB", # Syrphidae-Eristalinae (hoverflies)
                      52735:  "#35B2ED", # Tipulinae (craneflies)
                      7157:   "#E2BE64", # Culicidae
                      7158:   "#A8FF2E", # genus Aedes
                      7164:   "#FAA42C", # Anopheles
                      7174:   "#FA2213", # Culex
                      46207:  "#FBEF2F", # Toxorhynchites
                      53549:  "#AB9C01", # Tribe Sabethini
                      139055: "#C231FF", # Uranotaenia
                   2697495:  "#2432BF", # Spiralia
                    6178:    "#24ED54", # Trematoda
                    6447:    "#0CD6F2", # Mollusca
                      6563:       "#2451B8", # Ostreidae (true oysters)
                      32584:      "#338FEF", # Scaphopoda (tusk shells)
                      6544:     "#3B9110", # Bivalvia
                        6547:     "#304C1F", # Mytilidae
                        278205:   "#57D076", # Myida
                        6592:     "#41BF97", # Veneridae
                      6448:     "#BD47EA", # Gastropoda
                        6463:     "#691CAF", # Patella
                        6579:     "#1FE21D", # Pecten maximus
                      6605:   "#F2B70C", # Cephalopod
                        215450:  "#FA9420", # Decapodiformes
                        215451:  "#B2991D", # Octopodiformes
                    6340:    "#0CF247", # Annelida
                   2697496:  "#7624BF", # Gnathifera
                  33511:   "#78A6AF", # Deuterostomes
                   # Chordates 7711
                    # Vertebrates 7742
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
                        7878: "#DEBA97", # Dipnoi - Lungfish
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
                  6040 :   "#682CAE", # Sponges
                  6073 :   "#387FB2", # Cnidarians
                  10226:   "#C72480", # Placozoans
                  }
    # Secondary color dict for more specific subclades (empty by default, can be populated)
    color_dict = {}
    
    def __init__(self, perspchrom, calibrated_ages=None) -> None:
        # initialize the bidirectional graph using networkx
        self.G = nx.DiGraph()
        self.perspchrom = perspchrom
        # Before we delete anything, make a dict of samples to the taxidstring
        self.sample_to_taxidlist = self.perspchrom.set_index("species")["taxidstring"].to_dict()
        # Handle cases where taxidstring might be NaN (float) - report and filter those out
        invalid_samples = []
        valid_sample_to_taxidlist = {}
        for k, v in self.sample_to_taxidlist.items():
            if isinstance(v, str):
                valid_sample_to_taxidlist[k] = [int(x) for x in v.split(";")]
            else:
                invalid_samples.append((k, v))
        
        if invalid_samples:
            print(f"WARNING: {len(invalid_samples)} samples have invalid taxidstring (NaN or non-string):", file=sys.stderr)
            for sample, value in invalid_samples:
                print(f"  - {sample}: {value}", file=sys.stderr)
        
        self.sample_to_taxidlist = valid_sample_to_taxidlist
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
        
        # Rebuild lineage strings from actual topology to handle pruned nodes
        self._rebuild_lineage_strings_from_topology()

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

        # Store calibrated ages if provided and add them to graph nodes
        self.calibrated_ages = calibrated_ages
        if self.calibrated_ages is not None:
            self._add_node_ages()

        # assign the colors for all the nodes
        self.node_to_color = {}
        self.assign_colors()

        # Final df for the tree
        self.tree_df = None

    def _rebuild_lineage_strings_from_topology(self):
        """
        Build lineage strings from actual graph topology.
        This ensures consistency when nodes have been pruned from the calibrated tree.
        """
        print("  Rebuilding lineage strings from tree topology", file=sys.stderr)
        
        updated_count = 0
        conflicts = []  # Collect all conflicts for detailed error reporting
        cycles = []     # Collect cycle information
        
        for sample in self.leaves:
            # Traverse from leaf to root, collecting taxid nodes
            lineage = []
            current = sample
            visited = set()  # Prevent infinite loops
            
            while current is not None:
                if current in visited:
                    cycles.append((sample, current))
                    break
                visited.add(current)
                
                # Only include integer taxid nodes in lineage
                if isinstance(current, int):
                    lineage.append(current)
                
                # Get parent
                predecessors = list(self.G.predecessors(current))
                if len(predecessors) == 0:
                    break
                elif len(predecessors) == 1:
                    current = predecessors[0]
                else:
                    # Record conflict: node has multiple predecessors (reticulation)
                    conflicts.append((sample, current, predecessors))
                    break
            
            # Reverse to get root → leaf order
            lineage.reverse()
            
            # Update if different from original
            if sample in self.sample_to_taxidlist:
                if self.sample_to_taxidlist[sample] != lineage:
                    updated_count += 1
            self.sample_to_taxidlist[sample] = lineage
        
        # Check for topology conflicts and crash with detailed diagnostics
        if conflicts or cycles:
            error_msg = "\n" + "="*80 + "\n"
            error_msg += "ERROR: Graph topology conflicts detected!\n"
            error_msg += "="*80 + "\n\n"
            
            if conflicts:
                error_msg += f"Found {len(conflicts)} node(s) with MULTIPLE PREDECESSORS:\n"
                error_msg += "This indicates the graph contains references to non-existent internal nodes.\n\n"
                error_msg += "CAUSE: perspchrom.tsv contains lineage strings with taxids that were\n"
                error_msg += "pruned from the calibrated tree (internal nodes with only one child).\n"
                error_msg += "These pruned nodes create edges to nodes that already have parents in\n"
                error_msg += "the calibrated tree, resulting in multiple predecessors.\n\n"
                error_msg += "SOLUTION: Delete cache files and rebuild from calibrated tree:\n"
                error_msg += "  rm locdf.tsv perspchrom.tsv\n"
                error_msg += "  bash run_persp_chr.sh\n\n"
                error_msg += "First 10 conflicts (sample → node → predecessors):\n"
                for i, (sample, node, preds) in enumerate(conflicts[:10]):
                    error_msg += f"  {i+1}. Sample '{sample}' → Node {node} has predecessors: {preds}\n"
                if len(conflicts) > 10:
                    error_msg += f"  ... and {len(conflicts) - 10} more conflicts\n"
            
            if cycles:
                error_msg += f"\nFound {len(cycles)} cycle(s) in graph:\n"
                for sample, node in cycles[:10]:
                    error_msg += f"  Sample '{sample}' → Cycle at node {node}\n"
            
            error_msg += "\n" + "="*80 + "\n"
            
            print(error_msg, file=sys.stderr)
            raise IOError(error_msg)
        
        print(f"    Updated {updated_count} lineage strings from topology", file=sys.stderr)

    def _add_node_ages(self):
        """
        Add node ages from calibrated tree to graph nodes.
        
        Node naming convention:
        - Internal nodes: integers (taxids)
        - Leaf nodes: strings (sample names)
        
        For leaf nodes, we look up the terminal taxid from their lineage.
        For internal nodes, we use the node name directly as taxid.
        """
        print("  Adding node ages from calibrated tree to graph nodes", file=sys.stderr)
        nodes_with_ages = 0
        
        for node in self.G.nodes():
            node_taxid = None
            
            # Determine node type and extract taxid
            if isinstance(node, int):
                # Internal node - node name IS the taxid
                node_taxid = node
            elif isinstance(node, str):
                # Leaf node - sample name, get taxid from lineage
                if node in self.sample_to_taxidlist:
                    taxidlist = self.sample_to_taxidlist[node]
                    if len(taxidlist) > 0:
                        # Last taxid in lineage is the terminal taxid for this species
                        node_taxid = taxidlist[-1]
            
            # Look up age if we found a taxid
            if node_taxid is not None and node_taxid in self.calibrated_ages:
                self.G.nodes[node]['nodeage'] = self.calibrated_ages[node_taxid]
                nodes_with_ages += 1
        
        print(f"    Added ages to {nodes_with_ages} nodes", file=sys.stderr)

    def calculate_branch_lengths(self):
        """
        Calculate branch lengths in millions of years for all edges.
        Branch length = parent.nodeage - child.nodeage
        
        Stores 'branch_length_mya' attribute on each edge.
        Only calculates for nodes that have nodeage attributes.
        """
        print("  Calculating branch lengths from node ages", file=sys.stderr)
        edges_with_lengths = 0
        
        for parent, child in self.G.edges():
            parent_age = self.G.nodes[parent].get('nodeage', None)
            child_age = self.G.nodes[child].get('nodeage', None)
            
            if parent_age is not None and child_age is not None:
                branch_length = parent_age - child_age
                if branch_length < 0:
                    print(f"    WARNING: Negative branch length {branch_length:.2f} for edge {parent} -> {child}", file=sys.stderr)
                    branch_length = 0.0
                self.G.edges[parent, child]['branch_length_mya'] = branch_length
                edges_with_lengths += 1
        
        print(f"    Calculated branch lengths for {edges_with_lengths} edges", file=sys.stderr)

    def calculate_event_rates(self):
        """
        Calculate fusion and loss rates per million years for all edges.
        Rate = number_of_events / branch_length_mya
        
        Stores 'fusions_per_my' and 'losses_per_my' attributes on each edge.
        Only calculates for edges that have branch_length_mya attributes.
        Handles zero branch length by setting rate to NaN.
        """
        print("  Calculating event rates (fusions and losses per MY)", file=sys.stderr)
        edges_with_rates = 0
        
        for parent, child in self.G.edges():
            branch_length = self.G.edges[parent, child].get('branch_length_mya', None)
            
            if branch_length is not None:
                # Get fusion and loss counts from edge attributes
                num_fusions = self.G.edges[parent, child].get('num_fusions', 0)
                num_losses = self.G.edges[parent, child].get('num_losses', 0)
                
                # Calculate rates, handle zero branch length
                if branch_length > 0:
                    fusions_per_my = num_fusions / branch_length
                    losses_per_my = num_losses / branch_length
                else:
                    # Zero branch length - set to NaN
                    fusions_per_my = float('nan')
                    losses_per_my = float('nan')
                
                self.G.edges[parent, child]['fusions_per_my'] = fusions_per_my
                self.G.edges[parent, child]['losses_per_my'] = losses_per_my
                edges_with_rates += 1
        
        print(f"    Calculated rates for {edges_with_rates} edges", file=sys.stderr)

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

        # Add node name column at the beginning for easier grepping
        node_names = [str(x) for x in all_dfs.index]
        all_dfs = pd.concat([pd.DataFrame({'node': node_names}, index=all_dfs.index),
                             all_dfs],
                             axis=1)

        # Add the colors column to the dataframe. Use concat because pandas
        #  complains about performance if I try to add a column using
        #  older spells: a la all_dfs["color"] = colors
        colors = [self.node_to_color[x] for x in all_dfs.index]
        all_dfs = pd.concat([pd.DataFrame({'color': colors}, index=all_dfs.index),
                             all_dfs],
                             axis=1)

        # Add node ages if available (from calibrated tree)
        if self.calibrated_ages is not None:
            nodeages = []
            for node in all_dfs.index:
                node_age = self.G.nodes[node].get('nodeage', float('nan'))
                nodeages.append(node_age)
            all_dfs = pd.concat([all_dfs,
                                 pd.DataFrame({'nodeage_mya': nodeages}, index=all_dfs.index)],
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
            # Add debugging information about the problematic node
            print(f"\nERROR: Multiple predecessors found!", file=sys.stderr)
            print(f"  Current node: {node}", file=sys.stderr)
            print(f"  Type: {type(node)}", file=sys.stderr)
            print(f"  Predecessors: {predecessors}", file=sys.stderr)
            
            # If it's a leaf node (string), show its lineage
            if isinstance(node, str) and node in self.sample_to_taxidlist:
                print(f"  Lineage: {';'.join(map(str, self.sample_to_taxidlist[node]))}", file=sys.stderr)
            
            # Try to identify what taxa these are
            try:
                from ete4 import NCBITaxa
                ncbi = NCBITaxa()
                for pred in predecessors:
                    if isinstance(pred, int):
                        name = ncbi.get_taxid_translator([pred]).get(pred, "Unknown")
                        print(f"    {pred}: {name}", file=sys.stderr)
            except:
                pass
                
            raise IOError(f"There should only be one predecessor in a phylogenetic tree. Node: {node}, Found predecessors: {predecessors}")
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
        total_nodes = len(self.G.nodes())
        processed_nodes = len([n for n in self.G.nodes() if self.G.nodes[n]["completed"]])
        last_report = 0
        
        while len(queue) > 0:
            thisnode = queue.pop(0)
            
            # Report progress every 100 nodes
            current_processed = len([n for n in self.G.nodes() if self.G.nodes[n]["completed"]])
            if current_processed - last_report >= 100:
                print(f"\r  Percolating: {current_processed}/{total_nodes} nodes completed, queue: {len(queue)}", 
                      end="", file=sys.stderr)
                last_report = current_processed

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
        
        # Print final completion message
        final_completed = len([n for n in self.G.nodes() if self.G.nodes[n]["completed"]])
        print(f"\r  Percolating: {final_completed}/{total_nodes} nodes completed - Done!        ", file=sys.stderr)
        print("", file=sys.stderr)  # New line after progress
        print("", file=sys.stderr)  # Additional blank line for spacing

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
        pdf.loc[:, results.index] = results.values
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

def compute_changestring_for_species(species_row, perspchrom, ALG_columns, ALG_combos,
                                       min_for_missing, min_for_noncolocalized, checkpoint_dir, ALG_node):
    """
    Compute changestring for a single species.
    
    Parameters:
    -----------
    species_row : pandas.Series
        Row from perspchrom containing species data
    perspchrom : pandas.DataFrame
        Full perspchrom dataframe for comparison queries
    ALG_columns : list
        List of ALG column names
    ALG_combos : list
        List of ALG combination tuple column names
    min_for_missing : float
        Threshold for considering ALG missing in clade (default 0.8)
    min_for_noncolocalized : float
        Threshold for considering ALG pair non-colocalized in clade (default 0.5)
    checkpoint_dir : str
        Directory to save/load checkpoints
    ALG_node : str
        Node where ALGs were inferred
        
    Returns:
    --------
    tuple : (changestring, was_loaded)
        changestring: str - the computed or loaded changestring
        was_loaded: bool - True if loaded from checkpoint, False if computed
    """
    # Check for existing checkpoint first
    checkpoint_file = os.path.join(checkpoint_dir, f"{species_row['species']}.txt")
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'r') as f:
            return (f.read().strip(), True)
    
    # Initialize changeString for this species
    changeString = []
    thistaxidstring = species_row["taxidstring"]
    
    # --------- CONSERVED ----------- #
    ALGsConservedInThisSample = [x for x in ALG_columns if species_row[x] >= 1]
    
    # ---------- MISSING ------------ #
    ALGsMissingInThisSample = [x for x in ALG_columns if species_row[x] == 0]
    ALGsMissingInThisSampleAccountedFor = []
    
    # ---------- COLOCALIZED -------- #
    ColocalizedPairsInThisSample = [x for x in ALG_combos if species_row[x] >= 1]
    ColocalizedPairsInThisSampleAccountedFor = []
    
    # ---------- SPLIT -------------- #
    ALGsSplitInThisSample = [x for x in ALG_columns if species_row[x] > 1]
    ALGsSplitInThisSampleAccountedFor = []
    
    # Walk through lineage from tip to root
    while thistaxidstring.count(";") > 0:
        prevtaxidstring = thistaxidstring
        thisdf = perspchrom[perspchrom["taxidstring"] == thistaxidstring]
        
        # Move up one node
        thistaxidstring = ";".join(thistaxidstring.split(";")[:-1])
        
        # Initialize event lists for this branch
        ALG_colocalizations_on_this_branch = []
        ALG_losses_on_this_branch = []
        ALG_splits_on_this_branch = []
        
        # At the ALG inference node, all remaining events are placed here
        if thistaxidstring == ALG_node:
            ALG_losses_on_this_branch = [x for x in ALGsMissingInThisSample if x not in ALGsMissingInThisSampleAccountedFor]
            ALGsMissingInThisSampleAccountedFor.extend(ALG_losses_on_this_branch)
            
            ALG_colocalizations_on_this_branch = [x for x in ColocalizedPairsInThisSample if x not in ColocalizedPairsInThisSampleAccountedFor]
            ColocalizedPairsInThisSampleAccountedFor.extend(ALG_colocalizations_on_this_branch)
            
            ALG_splits_on_this_branch = [x for x in ALGsSplitInThisSample if x not in ALGsSplitInThisSampleAccountedFor]
            ALGsSplitInThisSampleAccountedFor.extend(ALG_splits_on_this_branch)
        else:
            # Compare to other clades at this level
            subdf = perspchrom[perspchrom["taxidstring"].str.contains(thistaxidstring) &
                               ~perspchrom["taxidstring"].str.contains(prevtaxidstring)]
            
            if len(subdf) > 0:
                # ----- ALG LOSSES -----
                missingALGsInThisClade, presentALGsInThisClade = missing_present_ALGs(subdf, min_for_missing=min_for_missing)
                ALGsLostOnThisBranch = [x for x in presentALGsInThisClade
                                         if x in ALGsMissingInThisSample
                                         and x not in ALGsMissingInThisSampleAccountedFor]
                ALGsMissingInThisSampleAccountedFor.extend(ALGsLostOnThisBranch)
                ALG_losses_on_this_branch.extend(ALGsLostOnThisBranch)
                
                # ----- ALG COLOCALIZATION EVENTS -----
                notLocalizedInThisClade = separate_ALG_pairs(subdf, min_for_noncolocalized=min_for_noncolocalized)
                pairsColocalizedOnThisBranch = [x for x in notLocalizedInThisClade
                                                 if x in ColocalizedPairsInThisSample
                                                 and x not in ColocalizedPairsInThisSampleAccountedFor]
                ColocalizedPairsInThisSampleAccountedFor.extend(pairsColocalizedOnThisBranch)
                ALG_colocalizations_on_this_branch.extend(pairsColocalizedOnThisBranch)
                
                # ----- ALG SPLITS -----
                notSplitInThisClade = unsplit_ALGs(subdf, max_frac_split=min_for_noncolocalized)
                ALGsplitsThisBranch = [x for x in notSplitInThisClade
                                        if x in ALGsSplitInThisSample
                                        and x not in ALGsSplitInThisSampleAccountedFor]
                ALGsSplitInThisSampleAccountedFor.extend(ALGsplitsThisBranch)
                ALG_splits_on_this_branch.extend(ALGsplitsThisBranch)
        
        thisChange = "({}|{}|{})".format(
            sorted(ALG_colocalizations_on_this_branch),
            sorted(ALG_losses_on_this_branch),
            sorted(ALG_splits_on_this_branch))
        changeString.append(prevtaxidstring.split(";")[-1])
        changeString.append(thisChange)
    
    # Handle unaccounted ALG losses (widespread ancient losses)
    if sorted(ALGsMissingInThisSampleAccountedFor) != sorted(ALGsMissingInThisSample):
        unaccounted_losses = [x for x in ALGsMissingInThisSample 
                              if x not in ALGsMissingInThisSampleAccountedFor]
        if unaccounted_losses:
            ALGsMissingInThisSampleAccountedFor.extend(unaccounted_losses)
    
    # Verify accounting for losses
    if sorted(ALGsMissingInThisSampleAccountedFor) != sorted(ALGsMissingInThisSample):
        raise IOError(f"ALG loss accounting failed for {species_row['species']}")
    
    # Handle unaccounted colocalizations (widespread)
    if sorted(ColocalizedPairsInThisSampleAccountedFor) != sorted(ColocalizedPairsInThisSample):
        unaccounted_colocs = [x for x in ColocalizedPairsInThisSample 
                              if x not in ColocalizedPairsInThisSampleAccountedFor]
        if unaccounted_colocs:
            ColocalizedPairsInThisSampleAccountedFor.extend(unaccounted_colocs)
    
    if sorted(ColocalizedPairsInThisSampleAccountedFor) != sorted(ColocalizedPairsInThisSample):
        raise IOError(f"ALG colocalization accounting failed for {species_row['species']}")
    
    # Handle unaccounted splits (widespread)
    if sorted(ALGsSplitInThisSampleAccountedFor) != sorted(ALGsSplitInThisSample):
        unaccounted_splits = [x for x in ALGsSplitInThisSample 
                              if x not in ALGsSplitInThisSampleAccountedFor]
        if unaccounted_splits:
            ALGsSplitInThisSampleAccountedFor.extend(unaccounted_splits)
    
    if sorted(ALGsSplitInThisSampleAccountedFor) != sorted(ALGsSplitInThisSample):
        raise IOError(f"ALG split accounting failed for {species_row['species']}")
    
    # Finalize changestring - add the root node
    changeString.append(thistaxidstring)
    changeString = "-".join(changeString[::-1])
    
    # Save checkpoint
    with open(checkpoint_file, 'w') as f:
        f.write(changeString)
    
    return (changeString, False)

def process_species_batch(args):
    """
    Process a batch of species indices in parallel.
    Worker function for multiprocessing.Pool.
    
    Parameters:
    -----------
    args : tuple
        (indices, perspchrom, ALG_columns, ALG_combos, min_for_missing, min_for_noncolocalized, checkpoint_dir, ALG_node)
        
    Returns:
    --------
    list : List of (index, changestring, was_loaded) tuples
    """
    indices, perspchrom, ALG_columns, ALG_combos, min_for_missing, min_for_noncolocalized, checkpoint_dir, ALG_node = args
    
    results = []
    for i in indices:
        row = perspchrom.iloc[i]
        changestring, was_loaded = compute_changestring_for_species(
            row, perspchrom, ALG_columns, ALG_combos,
            min_for_missing, min_for_noncolocalized, checkpoint_dir, ALG_node
        )
        results.append((i, changestring, was_loaded))
    
    return results

def monitor_progress(checkpoint_dir, total, stop_event):
    """
    Monitor progress by counting completed checkpoint files.
    Runs in background thread during parallel processing.
    """
    while not stop_event.is_set():
        try:
            completed = len([f for f in os.listdir(checkpoint_dir) if f.endswith('.txt')])
            print(f"\r  Progress: {completed}/{total} species completed          ", end="", file=sys.stderr)
        except:
            pass
        time.sleep(2)

def main():
    """
    Main execution function for ALG fusion analysis.
    
    WORKFLOW:
    =========
    1. Parse command line arguments
    2. Load or calculate location dataframe (locdf) and perspective chromosome dataframe (perspchrom)
       - locdf: Which ALGs are on which scaffolds in each species
       - perspchrom: Presence/absence matrix + colocalization data
    3. Build phylogenetic tree structure with evolutionary events
    4. Generate UMAP visualization of chromosomal architecture patterns
    
    CACHING BEHAVIOR:
    =================
    The script uses intelligent caching to speed up re-runs:
    
    - If locdf.tsv AND perspchrom.tsv exist:
      * Reads cached files (fast)
      * Sets overwrite=False
      
    - If either file is missing:
      * Calculates both from scratch (slow, could take hours)
      * Saves to locdf.tsv and perspchrom.tsv
      * Sets overwrite=True
      
    - If overwrite=True:
      * Constructs full phylogenetic tree with event strings
      * Saves to tree1.tsv.gz
      
    - If overwrite=False:
      * Reads cached tree1.tsv.gz if it exists
      * Raises error if tree1.tsv.gz doesn't exist
    
    OUTPUTS GENERATED:
    ==================
    Always generated:
    - tree1_umap.html: Interactive UMAP visualization (Plotly)
    
    Generated on first run (then cached):
    - locdf.tsv: ALG locations on scaffolds
    - perspchrom.tsv: Presence/absence and colocalization matrix
    - per_species_ALG_presence_fusions.tsv: Extended perspchrom with changestrings
    - tree1.tsv.gz: Complete phylogenetic tree with events
    
    FORCE RECALCULATION:
    ====================
    To force complete recalculation, delete:
    - locdf.tsv
    - perspchrom.tsv
    - tree1.tsv.gz (if you want to rebuild the tree too)
    
    NOTES:
    ======
    - Progress is printed to stderr
    - Uses NCBI taxonomy database (must be initialized first)
    - Exits after generating UMAP (sys.exit() before deprecated code sections)
    """
    # parse the args
    args = parse_args()

    # Load calibrated tree if provided
    calibrated_tree_data = None
    if args.tree_info:
        calibrated_tree_data = load_calibrated_tree(args.tree_info)
    else:
        print("No calibrated tree provided (--tree_info), will use NCBI taxonomy with custom phylogeny", file=sys.stderr)

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
        print("  locdf contains ALG localization data: which ALGs appear on which", file = sys.stderr)
        print("    scaffolds in each sample, with significance values, gene counts,", file = sys.stderr)
        print("    and fraction of each ALG on each scaffold.", file = sys.stderr)
        print("  perspchrom contains species-level summary data: taxonomic lineages,", file = sys.stderr)
        print("    chromosome counts, ALG presence/absence (number of scaffolds per", file = sys.stderr)
        print("    ALG), and ALG colocalization patterns (pairs of ALGs on same scaffold).", file = sys.stderr)
        locdf, perspchrom = rbh_files_to_locdf_and_perspchrom(rbh_files, args.ALG_rbh,
                                                              args.minsig, args.ALGname,
                                                              calibrated_tree_data=calibrated_tree_data)
        # save the locdf and perspchrom to a file
        locdf.to_csv("locdf.tsv", sep = "\t", index = False)
        perspchrom.to_csv("perspchrom.tsv", sep = "\t", index = False)
        overwrite = True

    ## filters for testing just get the chordates. Check if ";7711;" is in the taxidstring
    #perspchrom = perspchrom[perspchrom['taxidstring'].str.contains(';7711;')]

    resultstsv = "tree1.tsv.gz"
    resultsdf = None
    
    # If tree file doesn't exist, we need to build it regardless of overwrite flag
    if not os.path.exists(resultstsv):
        overwrite = True
    
    if overwrite:
        # now that we have the perspchrom, we should construct the tree structure
        # Pass calibrated ages if available
        calibrated_ages = None
        if calibrated_tree_data is not None:
            calibrated_ages = calibrated_tree_data['ages']
        
        T = SplitLossColocTree(perspchrom, calibrated_ages=calibrated_ages)
        print("Percolating", file = sys.stderr)
        T.percolate()
        
        # Calculate branch lengths and rates if we have calibrated ages
        if calibrated_ages is not None:
            T.calculate_branch_lengths()
            T.calculate_event_rates()
        
        print("Saving the tree to tsv", file = sys.stderr)
        T.save_tree_to_df(resultstsv)
        resultsdf = T.tree_df
    else:
        # tree1.tsv.gz exists, just read it
        resultsdf = pd.read_csv(resultstsv, sep = "\t", index_col = 0)
    #
    print("Making a UMAP with matplotlib", file = sys.stderr)
    #save_UMAP(resultsdf)
    print("Making a UMAP with plotly", file = sys.stderr)
    try:
        save_UMAP_plotly(resultsdf, "tree1")
    except ModuleNotFoundError:
        print("  Skipping plotly UMAP (plotly not installed)", file = sys.stderr)

    # Continue to changestring generation
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

    print("\nGenerating changestrings for each species", file=sys.stderr)
    print(f"  Processing {len(perspchrom)} species to identify ALG fusions, losses, and splits\n", file=sys.stderr)

    # Create checkpoint directory for crash recovery and future parallelization
    checkpoint_dir = "changestring_checkpoints"
    os.makedirs(checkpoint_dir, exist_ok=True)
    checkpoints_loaded = 0
    checkpoints_computed = 0

    # This is a magic number used to determine whether an ALG is missing or not.
    #  If an ALG is missing in 80%, or 0.8, of the clade in question, then we consider it missing.
    min_for_missing = 0.8
    # If the other clade's ALGs are detectable, and are not colocalized in at least 50% of the clade,
    #  then we consider that this ALG pair is not colocalized in this clade.
    min_for_noncolocalized = 0.5

    if args.parallel:
        # ============= PARALLEL MODE =============
        print(f"  Using {args.ncores} parallel processes for changestring generation", file=sys.stderr)
        
        # Split indices into batches for workers
        all_indices = list(range(len(perspchrom)))
        batch_size = max(1, len(perspchrom) // args.ncores)
        batches = [all_indices[i:i+batch_size] 
                   for i in range(0, len(all_indices), batch_size)]
        
        print(f"  Processing {len(perspchrom)} species in {len(batches)} batches", file=sys.stderr)
        
        # Prepare worker arguments (all batches get same parameters)
        worker_args = [(batch, perspchrom, ALG_columns, ALG_combos,
                        min_for_missing, min_for_noncolocalized, checkpoint_dir, ALG_node)
                       for batch in batches]
        
        # Start progress monitoring thread
        stop_event = threading.Event()
        monitor_thread = threading.Thread(target=monitor_progress,
                                           args=(checkpoint_dir, len(perspchrom), stop_event))
        monitor_thread.start()
        
        # Process in parallel using multiprocessing Pool with immediate error handling
        checkpoints_loaded = 0
        checkpoints_computed = 0
        try:
            with Pool(args.ncores) as pool:
                # Use imap_unordered to get results as they complete (fail fast on errors)
                for batch_results in pool.imap_unordered(process_species_batch, worker_args):
                    # Collect results from this batch
                    for i, changestring, was_loaded in batch_results:
                        perspchrom.at[i, "changestrings"] = changestring
                        if was_loaded:
                            checkpoints_loaded += 1
                        else:
                            checkpoints_computed += 1
        except Exception as e:
            # If any worker fails, stop progress monitoring and terminate pool
            stop_event.set()
            monitor_thread.join()
            print(f"\n\nERROR: Worker process failed. Shutting down all workers.", file=sys.stderr)
            raise
        
        # Stop the progress monitoring thread
        stop_event.set()
        monitor_thread.join()
        
        print(f"\r  Changestrings complete (parallel): {checkpoints_loaded} loaded from checkpoints, {checkpoints_computed} computed - Done!     ", file=sys.stderr)
    
    else:
        # ============= SEQUENTIAL MODE (ORIGINAL) =============
        print("  Using sequential processing for changestring generation", file=sys.stderr)
        
        for i, row in perspchrom.iterrows():
            # Check for existing checkpoint
            checkpoint_file = os.path.join(checkpoint_dir, f"{row['species']}.txt")
            if os.path.exists(checkpoint_file):
                with open(checkpoint_file, 'r') as f:
                    changeString = f.read().strip()
                perspchrom.at[i, "changestrings"] = changeString
                checkpoints_loaded += 1
                if (i+1) % 100 == 0 or i == 0:
                    print(f"\r  Loaded from checkpoints: {checkpoints_loaded}, Computing: {checkpoints_computed}, Progress: {i+1}/{len(perspchrom)}          ", end="", file=sys.stderr)
                continue
            
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
            
            # Handle colocalizations that couldn't be placed anywhere in phylogeny
            # This occurs when ALG pairs are colocalized in >50% of both the focal clade and sister clades
            if sorted(ColocalizedPairsInThisSampleAccountedFor) != sorted(ColocalizedPairsInThisSample):
                unaccounted_colocs = [x for x in ColocalizedPairsInThisSample 
                                      if x not in ColocalizedPairsInThisSampleAccountedFor]
                if unaccounted_colocs:
                    # Mark these as accounted (ancient/widespread colocalizations)
                    ColocalizedPairsInThisSampleAccountedFor.extend(unaccounted_colocs)
            
            # Final check for colocalizations - should now pass
            if sorted(ColocalizedPairsInThisSampleAccountedFor) != sorted(ColocalizedPairsInThisSample):
                raise IOError("There is a discrepancy between the ColocalizedPairsInThisSampleAccountedFor and ColocalizedPairsInThisSample. Write more debugging code to figure out what the issue is, because I haven't worked on this yet.")
            
            # Handle splits that couldn't be placed anywhere in phylogeny
            # This occurs when ALGs are split in >50% of both the focal clade and sister clades
            # (widespread splits that predate observable divergences)
            if sorted(ALGsSplitInThisSampleAccountedFor) != sorted(ALGsSplitInThisSample):
                unaccounted_splits = [x for x in ALGsSplitInThisSample 
                                      if x not in ALGsSplitInThisSampleAccountedFor]
                if unaccounted_splits:
                    # Mark these as accounted (ancient/widespread splits)
                    ALGsSplitInThisSampleAccountedFor.extend(unaccounted_splits)
            
            # Final check - should now pass
            if sorted(ALGsSplitInThisSampleAccountedFor) != sorted(ALGsSplitInThisSample):
                error_msg = "\n" + "="*80 + "\n"
                error_msg += "ERROR: ALG split accounting discrepancy detected\n"
                error_msg += "="*80 + "\n\n"
                error_msg += f"Species: {row['species']} (taxid: {row['taxid']})\n\n"
                error_msg += f"ALGs split in this sample (expected): {len(ALGsSplitInThisSample)}\n"
                error_msg += f"  {sorted(ALGsSplitInThisSample)}\n\n"
                error_msg += f"ALGs split accounted for: {len(ALGsSplitInThisSampleAccountedFor)}\n"
                error_msg += f"  {sorted(ALGsSplitInThisSampleAccountedFor)}\n\n"
                
                missing_from_accounted = set(ALGsSplitInThisSample) - set(ALGsSplitInThisSampleAccountedFor)
                extra_in_accounted = set(ALGsSplitInThisSampleAccountedFor) - set(ALGsSplitInThisSample)
                
                if missing_from_accounted:
                    error_msg += f"Missing from accounted: {sorted(missing_from_accounted)}\n"
                if extra_in_accounted:
                    error_msg += f"Extra in accounted: {sorted(extra_in_accounted)}\n"
                
                error_msg += "\nThis may indicate an issue with the split detection logic or\n"
                error_msg += "with the clade comparison thresholds (min_for_noncolocalized).\n"
                error_msg += "="*80 + "\n"
                print(error_msg, file=sys.stderr)
                raise IOError(error_msg)
            # We have stepped out of the for loop, now we add the last taxidstring to the changeString
            changeString.append(thistaxidstring)
            # flip the changeString and make it a parsable string
            changeString = "-".join(changeString[::-1])
            # The final string will look something like this: 1-([]|[])-131567-([]|[])-2759-([]|[])-33154-([]|[])-33208-([]|['A1b', 'A2', 'B2', 'B3', 'C1', 'C2', 'D', 'Ea', 'F', 'G', 'H', 'I', 'J1', 'K', 'L', 'M', 'N', 'O1', 'P', 'Qa', 'Qb', 'Qc', 'Qd', 'R'])-6040-([]|[])-60882-([]|[])-60883-([]|[])-60884-([]|[])-472148-([]|[])-111877-([]|[])-111878
            # We add the final string to the perspchrom df
            perspchrom.at[i, "changestrings"] = changeString
            
            # Save checkpoint immediately
            with open(checkpoint_file, 'w') as f:
                f.write(changeString)
            checkpoints_computed += 1
            
            # print a progress bar on the same line to tell the user how much longer we have to go
            print(f"\r  Loaded from checkpoints: {checkpoints_loaded}, Computing: {checkpoints_computed}, Progress: {i+1}/{len(perspchrom)}          ", end="", file=sys.stderr)
    
        print(f"\r  Changestrings complete (sequential): {checkpoints_loaded} loaded from checkpoints, {checkpoints_computed} computed - Done!     ", file=sys.stderr)
    
    # Common final steps (both parallel and sequential)
    print(f"  Checkpoints saved to: {checkpoint_dir}/", file=sys.stderr)
    print()

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