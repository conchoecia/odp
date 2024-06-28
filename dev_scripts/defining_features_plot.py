#!/usr/bin/env python

"""
The point of this script is to make plots to visualize the genomic colocalization features that
  we found from the colocalization matrix. We will use the output of defining_features.py to
  identify the top defining features for each clade, and then will plot them in the context of
  the ODOG UMAP. Later we may implement plotting of these in genomes, or in ODOL UMAPs.
"""

import argparse
import numpy as np
import pandas as pd
import ete3

import odp_plotting_functions as odp_plot
from defining_features import load_coo

from PhyloTreeUMAP import algcomboix_file_to_dict

import rbh_tools

#plotting stuff
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.cm as cm
import sys


def parse_args():
    """
    We need to load the following files:
      - unique_pairs.tsv.gz  - This file is where we get the best feature from each clade
      - a path to a coo file - This file contains all the distances for each pair of genomes
      - a path to a sample file that matches the coo file
      - a path to a sample df with the umap coordinates that we will plot
      - a path to the coo combination file. We will use this to get the index to the feature
      - an rbh file. We will use this to annotate the pdfs with some relevant information
      - the min number of genomes that must be in a clade to plot it
    """
    parser = argparse.ArgumentParser(description="Plot defining features" )
    parser.add_argument("--unique_pairs_path", type=str, help="Path to the unique pairs file", required=True)
    parser.add_argument("--coo_path", type=str, help="Path to the coo file", required=True)
    parser.add_argument("--sample_df_path", type=str, help="Path to the sample df", required=True)
    parser.add_argument("--umapdf_path", type=str, help="Path to the umap df", required=True)
    parser.add_argument("--coo_combination_path", type=str, help="Path to the coo combination file", required=True)
    parser.add_argument("--rbh_file", type=str, help="Path to the rbh file", required=True)
    parser.add_argument("--min_num_samples", type=int, help="The minimum number of samples in a clade to plot it", default=10)
    return parser.parse_args()

def main():
    # first load the data
    odp_plot.format_matplotlib()

    # get the args
    args = parse_args()

    algrbh = rbh_tools.parse_rbh(args.rbh_file)
    print(algrbh)

    # ALGcomboix
    ALGcomboix = algcomboix_file_to_dict(args.coo_combination_path)
    ALGcomboix_reverse = {v:k for k,v in ALGcomboix.items()}

    umapdf   =  pd.read_csv(args.umapdf_path, sep="\t", index_col=0)
    print(umapdf)
    print(umapdf.columns)
    cdf    = pd.read_csv(args.sample_df_path, sep="\t", index_col=0)
    cdf["taxid_list"] = cdf["taxid_list"].apply(eval)
    print(cdf.columns)
    matrix = load_coo(cdf, args.coo_path, ALGcomboix, missing_value_as = np.nan)
    # make sure that the sampledf len is the same as the matrix len
    assert len(cdf) == matrix.shape[0]
    # print the column at index 658688. this is a numpy ndarray
    mdf = pd.DataFrame(matrix, index = cdf.index,
                      columns = [ALGcomboix_reverse[i] for i in range(matrix.shape[1])])
    mdf["sample"] = cdf["sample"]
    print(mdf)
    print(mdf.columns)

    NCBI = ete3.NCBITaxa()
    pairsdf = pd.read_csv(args.unique_pairs_path, sep="\t")
    all_cdf_taxids = set()
    for i, row in cdf.iterrows():
        all_cdf_taxids.update(row["taxid_list"])
    # we first should check that EVERY single taxid in pairsdf is in the cdf
    for i, row in pairsdf.iterrows():
        taxid = row["taxid"]
        if taxid not in all_cdf_taxids:
            # print the full lineage of the taxids, using the translated names from NCBI
            lineage = NCBI.get_lineage(taxid)
            translated = NCBI.get_taxid_translator(lineage)
            raise ValueError(f"Taxid {taxid} ({translated}) is not in the cdf - in other words maybe your pairs file is from a different set of organisms than your sample df")

    pairsdf["nodename"] = pairsdf["taxid"].apply(lambda x: NCBI.get_taxid_translator([int(x)])[int(x)])
    # sort by taxid
    pairsdf = pairsdf.sort_values("nodename")
    print(pairsdf)
    pdf_path = "definitive_colocalizations.pdf"

    # setup the colormap
    dotgray = "#D5D7DF"
    dotred  = "#EB1F22"
    fontsize = 8
    # Create the colormap
    cmap = LinearSegmentedColormap.from_list("gray_red", [dotred, dotgray])

    ## open a pdf, where we will write one page per NCBI taxid that we iterate through in pairsdf
    with PdfPages(pdf_path) as pdf:
        counter = 0
        for i, row in pairsdf.iterrows():
            if row["num_samples_in_taxid"] < args.min_num_samples:
                continue
            cladename = row["nodename"]
            unique_pairs = eval(row["unique_pairs"])
            if len(unique_pairs) == 0:
                # there's nothing in here, just skip it.
                continue
            zeroth_entry = eval(row["unique_pairs"])[0][0]
            pairname     = ALGcomboix_reverse[zeroth_entry]
            taxid = row["taxid"]
            #print(f"The cladename is {cladename},")
            #print(f"  and the row is {row}")
            #print(f"  and the zeroth entry is {zeroth_entry}")
            #print(f"  and the pairname is {pairname}")
            #print(f"  and this column in mdf is {mdf[pairname]}")
            # get the clade
            x = umapdf["UMAP1"].to_list()
            y = umapdf["UMAP2"].to_list()
            values = []
            # not every value in the umapdf will be present in the matrix df (mdf), so we have to default to nan
            for i, row in umapdf.iterrows():
                umapdf_sample = row["sample"]
                if umapdf_sample in mdf["sample"].values:
                    values.append(mdf[mdf["sample"] == umapdf_sample][pairname].values[0] + 1)
                else:
                    values.append(np.nan)
            # Normalize the data
            #minval = 1/9999999999 # the minval is the value that we use for high-distance encoding. This is like a 10Gb chromosome
            #maxval = np.max([x for x in values if not np.isnan(x)])
            minval = 1000 # the minval is the value that we use for high-distance encoding. This is like a 10Gb chromosome
            maxval = np.max([x for x in values if not np.isnan(x)])
            print("the min of values is {}".format(minval))
            print("the max of values is {}".format(maxval))
            norm = Normalize(vmin=minval, vmax=maxval)
            colors = []
            for v in values:
                if np.isnan(v):
                    colors.append(cmap(norm(maxval)))
                else:
                    colors.append(cmap(norm(v)))
            # get the indices where the values are not nan
            non_nan_indices = [i for i, v in enumerate(values) if not np.isnan(v)]
            print([x[i] for i in non_nan_indices])
            print([y[i] for i in non_nan_indices])
            print([values[i] for i in non_nan_indices])
            print([colors[i] for i in non_nan_indices])
            # now we just plot a scatter. Enforce that the aspect ratio is 1
            fig,ax = plt.subplots(1,1, figsize=(5,5))
            plt.gca().set_aspect('equal', adjustable='box')
            scatter = ax.scatter(x, y, c=colors, lw = 0, s = 3)
            # set the title
            title = f"Top pair for {cladename}, TaxID {taxid}\nis this pair:"
            pairname0 = pairname[0]
            pairname1 = pairname[1]
            ALG0 = algrbh[algrbh["rbh"] == pairname0]["gene_group"].values[0]
            ALG1 = algrbh[algrbh["rbh"] == pairname1]["gene_group"].values[0]
            title += f"\n{pairname0} - {ALG0}\n{pairname1} - {ALG1}"
            plt.title(title, fontsize = fontsize * 0.75)
            # set the x-axis to say "UMAP"
            plt.xlabel("UMAP1", fontsize = fontsize)
            # set the y-axis to say "UMAP"
            plt.ylabel("UMAP2", fontsize = fontsize)
            # Turn off the axis ticks
            plt.xticks([])
            plt.yticks([])

            # Create a scalar mappable object with the colormap and normalization
            norm = Normalize(vmin=maxval, vmax=minval)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])  # This is needed for the colorbar to work properly

            # Add a colorbar to the current figure
            cbar = plt.colorbar(sm, ax=ax, shrink = 0.5)
            #cbar.set_ticks([0.0, 1.0])
            #cbar.set_ticklabels([1/maxval, 1000])

            ## Customize the min and max values of the colorbar
            #cbar.set_clim(0.2, 0.8)  # Adjust these values as needed
            cbar.set_label('bp distance', fontsize = fontsize * 0.75)  # Optional: add a label
            cbar.ax.yaxis.get_offset_text().set_fontsize(fontsize * 0.75)
            cbar.ax.tick_params(labelsize=fontsize *0.75)  # Adjust the font size of the colorbar ticks

            # Save the current figure to the PDF
            pdf.savefig()
            # Close the current figure
            plt.close()
            counter += 1

if __name__ == "__main__":
    main()