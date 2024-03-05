#!/usr/bin/env python
"""
# Take in a list of datafraes from samples and constructs a comparison of the UMAP plots.
"""
import argparse
import numpy as np
import os
import pandas as pd
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import matplotlib.colors as mcolors

import odp_plotting_functions as odpf

def parse_args():
    """
    The thing we need to read in now is a list of dataframes.
    With this list of dataframes we will infer the parameters used from the file names.
    Then just make a plot in a grid from n_neighbors and min_dist. Make one plot each for small/large

    Args:
      -d --directory: The directory to read in the dataframes from.
      -f --filelist: The list of dataframes to read in. Space separated. We will infer the parameters from the filenames.
      -o --outdir: The directory to which we will write the plots.
      --plot_features: Looks in the DF for features to plot. Plots everything on the same plot. Only takes in one dataframe.
    """
    parser = argparse.ArgumentParser(description = "Take in a list of datafraes from samples and constructs a comparison of the UMAP plots.")
    parser.add_argument("-d", "--directory",help = "The directory to read in the dataframes from.")
    parser.add_argument("-f", "--filelist", help = "The list of dataframes to read in. Space separated. We will infer the parameters from the filenames.")
    parser.add_argument("-o", "--outpdf",   help = "The pdf file to which we want to save our results.", required = True)
    parser.add_argument("--plot_features", action = "store_true", help = "Looks in the DF for features to plot. Plots everything on the same plot. Only takes in one dataframe.")
    args = parser.parse_args()

    # Make sure that both directory and filelist are not specified.
    # If they are both specified, we don't know which one to use.
    if args.directory and args.filelist:
        raise ValueError("Both directory and filelist are specified. We don't know which one to use. Please just use one.")

    # If we have turned on the plot_features flag, then we need to make sure that we only have one file in the filelist.
    if args.plot_features:
        if args.filelist:
            if len(args.filelist.split(" ")) > 1:
                raise ValueError("You have turned on the plot_features flag, but you have more than one file in the filelist. We can only plot one file at a time with this flag.")
        elif args.directory:
            if len([x for x in os.listdir(args.directory) if x.endswith(".df")]) > 1:
                raise ValueError("You have turned on the plot_features flag, but you have more than one file in the directory. We can only plot one file at a time with this flag.")

    return args

def plot_paramsweep(args):
    """
    Makes the plot for the parameter sweep plot when we provide multiple dataframes.
    """
    df_filelist = []
    if args.directory:
        # Get all the dataframes from the directory
        df_filelist = [x for x in os.listdir(args.directory)
                       if x.endswith(".df")]
    elif args.filelist:
        # Get all the dataframes from the filelist
        df_filelist = args.filelist.split(" ")

    print("df_filelist is {}".format(df_filelist))

    # There filenames are formatted like this:
    #   Metazoa_33208.neighbors_20.mind_0.0.missing_large.df
    samplename = set([os.path.basename(x).split(".neighbors_")[0] for x in df_filelist])
    if len(samplename) > 1:
        raise ValueError("More than one sample name found in the file list. We only allow one for now.")
    samplename = list(samplename)[0]

    print("samplename is {}".format(samplename))

    # get the number of neighbors
    num_neighbors = [os.path.basename(x).split(".neighbors_")[1].split(".")[0] for x in df_filelist]
    # check that everything can be cast to an int
    for x in num_neighbors:
        try:
            int(x)
        except ValueError:
            raise ValueError(f"Could not cast {x} to an int.")
    num_neighbors = [int(x) for x in num_neighbors]
    print("num_neighbors is {}".format(num_neighbors))
    # get the min_dist
    min_dist = [os.path.basename(x).split(".mind_")[1].split(".missing_")[0] for x in df_filelist]
    # check that all min_dist can be cast to a floar
    for x in min_dist:
        try:
            float(x)
        except ValueError:
            raise ValueError(f"Could not cast {x} to a float.")
    min_dist = [float(x) for x in min_dist]
    print("min_dist is {}".format(min_dist))
    # get whether it is from the small or large dataset
    miss_size = [os.path.basename(x).split(".missing_")[1].split(".")[0] for x in df_filelist]
    # make sure that there is only one value for small or large here. We can't deal with both.
    if len(set(miss_size)) > 1:
        raise ValueError("More than one size found in the file list. We only allow one for now: small or large")
    print("Missing size is {}".format(miss_size))

    # collate together the html_filelist, num_neighbors, min_dist, and size
    df = pd.DataFrame({"dffile":        df_filelist,
                       "num_neighbors": num_neighbors,
                       "min_dist":      min_dist,
                       "size":          miss_size})
    # sort by num_neighbors and min_dist, ascending both
    df = df.sort_values(["num_neighbors", "min_dist"], ascending=[True, True])
    # groupby size and make one plot for each size
    gb = df.groupby("size")
    for name, group in gb:
        # make subplots.
        # The number of unique things in num_neighbors is the number of rows
        # the numer of unique things in min_dist is the number of columns
        num_rows = len(group["num_neighbors"].unique())
        num_cols = len(group["min_dist"].unique())
        # make a grid of squares to plot each of these on
        # reduce the space between all the plots
        # make the figure size such that all the plots are square
        square_size = 1.5
        fig, axes = plt.subplots(num_rows, num_cols,
                                 figsize=(num_cols*square_size, num_rows*square_size))
        plt.subplots_adjust(wspace=0.1, hspace=0.1)

        # turn off all the ticks
        # Turn off the axes
        for i in range(num_rows):
            for j in range(num_cols):
                axes[i, j].set_xticks([])
                axes[i, j].set_yticks([])
                # Turn off the lines around the plot
                for spine in ['top', 'right', 'bottom', 'left']:
                    axes[i, j].spines[spine].set_visible(False)

        # set the title as samplename and the whether it is small or large
        fig.suptitle(f"{samplename} {name} NaNs")
        # set absolute left label as the number of neighbors
        fig.text(0.06, 0.5, 'Number of Neighbors', va='center', rotation='vertical')
        fig.text(0.5, 0.92, 'Min Distance', ha='center')
        for i, row in enumerate(sorted(group["num_neighbors"].unique(), reverse = False)):
            # for the left-most plot in each row, set the ylabel as the number of neighbors
            axes[i, 0].set_ylabel(f"{row}", rotation=0, ha='right')
            for j, col in enumerate(sorted(group["min_dist"].unique(), reverse = False)):
                #for the top-most plot in each column, set the xlabel as the min_dist
                # put the xlabel on the top
                axes[0, j].xaxis.set_label_position('top')
                axes[0, j].set_xlabel(f"{col}")
                # get the df file for this row and column from the dffile
                dffile = group[(group["num_neighbors"] == row) & (group["min_dist"] == col)]["dffile"].values[0]
                # if the type of dffile is NoneType, then we didn't find the file
                if type(dffile) == type(None):
                    # write into the ax[i, j] that we didn't find the file
                    axes[i, j].text(0.5, 0.5, f"Missing file", fontsize=6, ha='center')
                elif type(dffile) == str:
                    # now we try to find the file
                    if os.path.exists(dffile):
                        # if the file is empty, then we write into the ax[i, j] that the file is empty
                        if os.path.getsize(dffile) == 0:
                            axes[i, j].text(0.5, 0.5, f"Empty file", fontsize=6, ha='center')
                        else:
                            tempdf = pd.read_csv(dffile, sep="\t", index_col=0)
                            # Plot the UMAP1 UMAP2 with the color column as the color.
                            # Make the dot size small
                            axes[i, j].scatter(tempdf["UMAP1"], tempdf["UMAP2"],
                                             s=0.5, lw = 0, alpha=0.5,
                                             color=list(tempdf["color"]))
                    else:
                        # The user provided the path to this file, but it doesn't exist.
                        # This means that the user made a mistake in writing the file name.
                        raise ValueError(f"The file {dffile} does not exist.")
                else:
                    raise ValueError(f"Type of dffile is not a string or NoneType. It is {type(dffile)}")
        # make sure that the aspect ratio for all of these is the same
        # Iterate over each axis and set aspect ratio to 'equal'
        for row in axes:
            for col in row:
                col.set_aspect('equal', adjustable='box')

        # Now make vertical and horizontal lines to separate the plots. Make them medium gray.
        # The lines will be on the figure, and not in the plots.
        # Draw horizontal lines between rows
        # Get the bounding boxes of the axes including text decorations

        # Get the bounding boxes of the axes including text decorations
        r = fig.canvas.get_renderer()
        get_bbox = lambda ax: ax.get_tightbbox(r).transformed(fig.transFigure.inverted())
        bbox_list = [get_bbox(ax) for ax in axes.flat]

        # Create an empty array with the correct shape and dtype
        bboxes = np.empty(axes.shape, dtype=object)

        # Fill the array with the bounding boxes
        for idx, bbox in np.ndenumerate(bboxes):
            bboxes[idx] = bbox_list[idx[0] * axes.shape[1] + idx[1]]

        # Get the minimum and maximum extent, get the coordinate half-way between those
        ymax = np.array(list(map(lambda b: b.y1, bboxes.flat))).reshape(axes.shape).max(axis=1)
        ymin = np.array(list(map(lambda b: b.y0, bboxes.flat))).reshape(axes.shape).min(axis=1)
        ys = np.c_[ymax[1:], ymin[:-1]].mean(axis=1)

        # Draw horizontal lines at those coordinates
        for y in ys:
            line = plt.Line2D([0.125, 0.9], [y, y], transform=fig.transFigure, color="#BBBBBB")
            fig.add_artist(line)

        # Get the minimum and maximum extent, get the coordinate half-way between those for vertical lines
        xmax = np.array(list(map(lambda b: b.x1, bboxes.flat))).reshape(axes.shape).max(axis=0)
        xmin = np.array(list(map(lambda b: b.x0, bboxes.flat))).reshape(axes.shape).min(axis=0)
        xs = np.c_[xmax[1:], xmin[:-1]].mean(axis=1)

        # Draw vertical lines at those coordinates
        for xi in range(len(xs)):
            x = xs[xi]
            if xi == 0:
                x = x + 0.0125
            line = plt.Line2D([x, x], [0.1, 0.875], transform=fig.transFigure, color="#BBBBBB")
            fig.add_artist(line)

        # save the figure as f"{samplename}_{name}.pdf"
        # name is just the size, small or large
        print("saving the file to {}".format(args.outpdf))
        plt.savefig(args.outpdf)
        # close the figure
        plt.close(fig)

def interpolate_color(value, vmin, vmax, start_color, end_color):
    """
    Interpolates between two colors based on a given value and a range.

    Parameters:
        value (float): The value to map to a color.
        vmin (float): The minimum value of the range.
        vmax (float): The maximum value of the range.
        start_color (str): Hexadecimal color string for the start color.
        end_color (str): Hexadecimal color string for the end color.

    Returns:
        str: Hexadecimal color string interpolated between start and end colors.
    """
    # Convert hex color strings to RGB tuples
    start_rgb = mcolors.hex2color(start_color)
    end_rgb = mcolors.hex2color(end_color)

    # Normalize value to range [0, 1]
    normalized_value = (value - vmin) / (vmax - vmin)

    # Interpolate RGB values
    interpolated_rgb = [
        start_rgb[channel] + normalized_value * (end_rgb[channel] - start_rgb[channel])
        for channel in range(3)
    ]

    # Convert interpolated RGB values back to hexadecimal color string
    interpolated_color = mcolors.rgb2hex(interpolated_rgb)

    return interpolated_color

def plot_features(args):
    """
    This makes a plot of the features of a single dataframe.
    All of the plots will be plotted along the axes of the UMAP1 and UMAP2.
    """
    # get the dataframe to load in
    df = pd.read_csv(args.filelist, sep="\t", index_col=0)
    # remove the rows that have a value of NaN in the "smallest_protein" column
    df = df[~df["smallest_protein"].isna()]
    print(df.columns)
    print(df)
    # the columns we wa:nt to annotate are defined in the AnnotateSampleDf pipeline
    #  - things from the genome fasta file
    #    - num_scaffolds
    #    - GC_content
    #    - genome_size
    #    - median_scaffold_length
    #    - mean_scaffold_length
    #    - scaffold_N50
    #    - longest_scaffold
    #    - smallest_scaffold
    #    - fraction_Ns
    #    - number_of_gaps
    #  - things from the protein file:
    #    - num_proteins
    #    - mean_protein_length
    #    - median_protein_length
    #    - longest_protein
    #    - smallest_protein
    #    - from_rbh
    #  - things from the rbh file:
    #    - frac_ologs:           The fraction of genes of ANY ALG that are present at all in the rbh file. len(rbhdf) / total_genes_ALGs
    #    - frac_ologs_sig:       The fraction of genes of ANY ALG that are significantly on any chromosome, as defined by whole_FET
    #    - frac_ologs_single:    The fraction of genes of ANY ALG that are significantly on the largest chromosome, as defined by whole_FET
    #    - frac_ologs_{ALGNAME}: The fraction of genes of INDIVIDUAL ALGs that are significantly on any chromosome

    regular_columns_to_plot = ["num_scaffolds", "GC_content", "genome_size", "median_scaffold_length",
                               "mean_scaffold_length", "scaffold_N50", "longest_scaffold", "smallest_scaffold",
                               "fraction_Ns", "number_of_gaps", "num_proteins", "mean_protein_length",
                               "median_protein_length", "longest_protein", "smallest_protein", "from_rbh",
                               "frac_ologs", "frac_ologs_sig", "frac_ologs_single"]
    olog_columns_to_plot = [x for x in df.columns if x.startswith("frac_ologs_") and (x != "frac_ologs_single") and (x != "frac_ologs_sig")]
    all_columns_to_plot = ["color"] + regular_columns_to_plot + olog_columns_to_plot
    total_num_cols_to_plot = len(all_columns_to_plot)
    # figure out what number we should pick to make the shape closest to a square, when plotting N x N plots
    plots_per_row = int(np.ceil(np.sqrt(total_num_cols_to_plot)))

    num_rows = plots_per_row
    num_cols = plots_per_row
    # make a grid of squares to plot each of these on
    # reduce the space between all the plots
    # make the figure size such that all the plots are square
    square_size = 1.5
    fig, axes = plt.subplots(num_rows, num_cols,
                             figsize=(num_cols*square_size, num_rows*square_size))
    plt.subplots_adjust(wspace=0.1, hspace=0.1)

    # turn off all the ticks
    # Turn off the axes
    for i in range(num_rows):
        for j in range(num_cols):
            axes[i, j].set_xticks([])
            axes[i, j].set_yticks([])
            # Turn off the lines around the plot
            for spine in ['top', 'right', 'bottom', 'left']:
                axes[i, j].spines[spine].set_visible(False)

    # set the title as samplename and the whether it is small or large
    fig.suptitle(f"Paramplot for {args.filelist}")
    # set absolute left label as the number of neighbors
    #fig.text(0.06, 0.5, 'Number of Neighbors', va='center', rotation='vertical')
    #fig.text(0.5, 0.92, 'Min Distance', ha='center')
    i = 0
    j = 0
    for thiscol in all_columns_to_plot:
        # for the left-most plot in each row, set the ylabel as the number of neighbors
        #axes[i, j].set_ylabel(f"{row}", rotation=0, ha='right')
        #for the top-most plot in each column, set the xlabel as the min_dist
        # put the xlabel on the top
        axes[i, j].xaxis.set_label_position('top')
        if thiscol == "color":
            axes[i, j].set_xlabel(f"Clade color", fontsize = 5)
        else:
            axes[i, j].set_xlabel(f"{thiscol}", fontsize = 5)
        # figure out the type of the column
        coltype = df[thiscol].dtype
        #if the type of the column is an object, then plot True as #074FF7 and False as #FD6117
        # If True/False not in that column , then randomly select a color
        colors = []
        print('coltype of {} is {}'.format(thiscol, coltype))
        if thiscol == "color":
            colors = list(df[thiscol])
        else:
            if ((True in df[thiscol].unique()) or (False in df[thiscol].unique())) and (coltype == "object"):
                # print the counts of True and False
                print(df[thiscol].value_counts())
                colordict = {True: "#074FF7", False: "#FD6117"}
                colors = [colordict[x] for x in df[thiscol]]
            else:
                # make an object to interpolate color between #DCDEE3 and #FF2608
                # if the max value is above 1, then we will interpolate between 0 and the max of the column
                if df[thiscol].max() > 1:
                    colors = [interpolate_color(x, 0, df[thiscol].max(), "#DCDEE3", "#FF2608") for x in df[thiscol]]
                else:
                    # otherwise go between 0 and 1
                    colors = [interpolate_color(x, 0, 1, "#DCDEE3", "#FF2608") for x in df[thiscol]]

        # get the df file for this row and column from the dffile
        axes[i, j].scatter(df["UMAP1"], df["UMAP2"],
                         s=0.5, lw = 0, alpha=0.5,
                         color=colors)
        # iterate i and j
        j += 1
        if j == num_cols:
            j = 0
            i += 1

    # make sure that the aspect ratio for all of these is the same
    # Iterate over each axis and set aspect ratio to 'equal'
    for row in axes:
        for col in row:
            col.set_aspect('equal', adjustable='box')

    # Now make vertical and horizontal lines to separate the plots. Make them medium gray.
    # The lines will be on the figure, and not in the plots.
    # Draw horizontal lines between rows
    # Get the bounding boxes of the axes including text decorations

    # Get the bounding boxes of the axes including text decorations
    r = fig.canvas.get_renderer()
    get_bbox = lambda ax: ax.get_tightbbox(r).transformed(fig.transFigure.inverted())
    bbox_list = [get_bbox(ax) for ax in axes.flat]

    # Create an empty array with the correct shape and dtype
    bboxes = np.empty(axes.shape, dtype=object)

    # Fill the array with the bounding boxes
    for idx, bbox in np.ndenumerate(bboxes):
        bboxes[idx] = bbox_list[idx[0] * axes.shape[1] + idx[1]]

    # Get the minimum and maximum extent, get the coordinate half-way between those
    ymax = np.array(list(map(lambda b: b.y1, bboxes.flat))).reshape(axes.shape).max(axis=1)
    ymin = np.array(list(map(lambda b: b.y0, bboxes.flat))).reshape(axes.shape).min(axis=1)
    ys = np.c_[ymax[1:], ymin[:-1]].mean(axis=1)

    # Draw horizontal lines at those coordinates
    for y in ys:
        line = plt.Line2D([0.125, 0.9], [y, y], transform=fig.transFigure, color="#BBBBBB")
        fig.add_artist(line)

    # Get the minimum and maximum extent, get the coordinate half-way between those for vertical lines
    xmax = np.array(list(map(lambda b: b.x1, bboxes.flat))).reshape(axes.shape).max(axis=0)
    xmin = np.array(list(map(lambda b: b.x0, bboxes.flat))).reshape(axes.shape).min(axis=0)
    xs = np.c_[xmax[1:], xmin[:-1]].mean(axis=1)

    # Draw vertical lines at those coordinates
    for xi in range(len(xs)):
        x = xs[xi]
        if xi == 0:
            x = x + 0.0125
        line = plt.Line2D([x, x], [0.1, 0.875], transform=fig.transFigure, color="#BBBBBB")
        fig.add_artist(line)

    # save the figure as f"{samplename}_{name}.pdf"
    # name is just the size, small or large
    print("saving the file to {}".format(args.outpdf))
    plt.savefig(args.outpdf)
    # close the figure
    plt.close(fig)

def main():
    odpf.format_matplotlib()
    args = parse_args()
    if args.plot_features:
        plot_features(args)
    else:
        plot_paramsweep(args)

if __name__ == "__main__":
    main()