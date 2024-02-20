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

def hex_to_rgb(value):
    """
    converts hex colors to rgb colors
    """

def parse_args():
    """
    The thing we need to read in now is a list of dataframes.
    With this list of dataframes we will infer the parameters used from the file names.
    Then just make a plot in a grid from n_neighbors and min_dist. Make one plot each for small/large

    Args:
      -d --directory: The directory to read in the dataframes from.
    """
    parser = argparse.ArgumentParser(description="Take in a list of datafraes from samples and constructs a comparison of the UMAP plots.")
    parser.add_argument("-d", "--directory", help="The directory to read in the dataframes from.", required=True)

    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    # Get all the dataframes from the directory
    html_filelist = [x for x in os.listdir(args.directory)
                   if x.endswith(".html")]
    # There filenames are formatted like this:
    #   Metazoa_33208.neighbors_20.mind_0.0.missing_large.df
    samplename = set([f"{x.split('_')[0]}_{x.split('_')[1].split('.')[0]}" for x in html_filelist])
    if len(samplename) > 1:
        raise ValueError("More than one sample name found in the file list. We only allow one for now.")
    samplename = list(samplename)[0]
    print("samplename is {}".format(samplename))

    # get the number of neighbors
    num_neighbors = [int(x.split(".")[1].split("_")[1]) for x in html_filelist]
    # get the min_dist
    min_dist = [x.split(".mind_")[1].split(".missing")[0] for x in html_filelist]
    # get whether it is from the small or large dataset
    size     = [x.split(".missing_")[1].split(".bokeh")[0] for x in html_filelist]

    # collate together the html_filelist, num_neighbors, min_dist, and size
    df = pd.DataFrame({"html_filelist": html_filelist,
                       "num_neighbors": num_neighbors,
                       "min_dist": min_dist, "size": size})
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
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols*square_size, num_rows*square_size))
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
        for i, row in enumerate(group["num_neighbors"].unique()):
            # for the left-most plot in each row, set the ylabel as the number of neighbors
            axes[i, 0].set_ylabel(f"{row}", rotation=0, ha='right')
            for j, col in enumerate(group["min_dist"].unique()):
                #for the top-most plot in each column, set the xlabel as the min_dist
                # put the xlabel on the top
                axes[0, j].xaxis.set_label_position('top')
                axes[0, j].set_xlabel(f"{col}")
                # get the df file for this row and column
                #Metazoa_33208.neighbors_20.mind_0.0.missing_large.df
                dffile = f"{samplename}.neighbors_{row}.mind_{col}.missing_{name}.df"
                # read in the html file
                fullpath = os.path.join(args.directory, dffile)
                if os.path.exists(fullpath):
                    print(f"We are looking at this file: {fullpath}")
                    tempdf = pd.read_csv(fullpath, sep="\t", index_col=0)
                    # Plot the UMAP1 UMAP2 with the color column as the color.
                    # Make the dot size small
                    axes[i, j].scatter(tempdf["UMAP1"], tempdf["UMAP2"],
                                     s=0.5, lw = 0, alpha=0.5,
                                     color=list(tempdf["color"]))
                else:
                    # write into the ax[i, j] that we didn't find the file
                    axes[i, j].text(0.5, 0.5, f"Missing file", fontsize=6, ha='center')

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
        plt.savefig(f"{samplename}_{name}.pdf")
        # close the figure
        plt.close(fig)

if __name__ == "__main__":
    main()