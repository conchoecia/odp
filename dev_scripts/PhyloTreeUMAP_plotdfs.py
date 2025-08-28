#!/usr/bin/env python
"""
# Take in a list of datafraes from samples and constructs a comparison of the UMAP plots.
"""
import argparse
import numpy as np
import os
import pandas as pd
import random
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import odp_plotting_functions as odpf

# for the html version
from bokeh.plotting import figure, output_file, save
from bokeh.layouts import gridplot, column
from bokeh.models import ColumnDataSource, Div, TapTool, OpenURL

import colorsys
def generate_distinct_colors(n, saturation=0.65, lightness=0.5):
    """
    Generate `n` visually distinct colors using evenly spaced HSL hue values.
    Returns a list of hex color strings.
    """
    colors = []
    for i in range(n):
        hue = i / n
        rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
        rgb_scaled = tuple(int(x * 255) for x in rgb)
        hex_color = '#{:02x}{:02x}{:02x}'.format(*rgb_scaled)
        colors.append(hex_color)
    return colors

# Official Benedictus stops (alpha 'FF' trimmed)
BENEDICTUS_HEX = [
    "#9A133D", "#B93961", "#D8527C", "#F28AAA", "#F9B4C9", "#F9E0E8",
    "#FFFFFF",
    "#EAF3FF", "#C5DAF6", "#A1C2ED", "#6996E3", "#4060C8", "#1A318B",
]

def benedictus_cmap(name="Benedictus", reverse=False, N=256):
    """Continuous diverging colormap from Benedictus stops."""
    stops = BENEDICTUS_HEX[::-1] if reverse else BENEDICTUS_HEX
    return LinearSegmentedColormap.from_list(name, stops, N=N)

def benedictus_listed(name="Benedictus_discrete", reverse=False):
    """Discrete colormap using the 13 official stops."""
    import matplotlib.colors as mcolors
    stops = BENEDICTUS_HEX[::-1] if reverse else BENEDICTUS_HEX
    return mcolors.ListedColormap(stops, name=name)

def benedictus_n(n, reverse=False):
    """Get n evenly spaced colors sampled from the continuous Benedictus map."""
    cmap = benedictus_cmap(reverse=reverse)
    return [mcolors.to_hex(cmap(t)) for t in np.linspace(0, 1, n)]

def parse_args():
    """
    The thing we need to read in now is a list of dataframes.
    With this list of dataframes we will infer the parameters used from the file names.
    Then just make a plot in a grid from n_neighbors and min_dist. Make one plot each for small/large

    Args:
      -d --directory: The directory to read in the dataframes from.
      -f --filelist: The list of dataframes to read in. Space separated. We will infer the parameters from the filenames.
      -p --prefix:   The files will be saved to this + ".pdf" or ".html"
      --genome-min-bp: Minimum genome size (bp). Values <= this are shown as --genome-min-color (grey).
      --genome-max-bp: Maximum genome size (bp). Values >= this are shown as --genome-max-color (red).
      --legend-scale: Scale factor for legend (colorbar) size and font. 1.0 = original; 0.5 = half size.
      --genome-min-color: Hex color for genome sizes <= genome-min-bp (default grey).
      --genome-max-color: Hex color for genome sizes >= genome-max-bp (if not set, uses cmap endpoint color).
      --benedictus: Use the Benedictus diverging colormap for genome size panels,
                    ignoring genome size min/max thresholds.
      --metadata: space-separated list of metadata files to join against the main dataframe. This will be type list [str] of filenames.
      --pdf: save a {prefix}.pdf file
      --html: save a {prefix}.html file
      --plot_features: Looks in the DF for features to plot. Plots everything on the same plot. Only takes in one dataframe.
    """
    parser = argparse.ArgumentParser(description = "Take in a list of datafraes from samples and constructs a comparison of the UMAP plots.")
    parser.add_argument("-d", "--directory",help = "The directory to read in the dataframes from.")
    flstr  = "The list of dataframes to read in. Space separated. We will infer the parameters from the filenames.\n"
    flstr += "  This cannot be used in combination with the --directory flag.\n"
    flstr += "  If you use this in combination with the --plot_features flag, this must only be one file."
    parser.add_argument("-f", "--filelist", help = flstr)
    parser.add_argument("-p", "--prefix",   help = "The pdf file to which we want to save our results.", required = True)
    mdstr  = "Optional metadata files with one 'rbh' column to join against the main dataframe, and other columns to annotate the plot."
    mdstr += "  If the metadata column does not contain an additional column called *_color, the dots will be assigned colors."
    mdstr += "  The colors will be assigned a gradient if numeric, or a random color if the column if categorical."
    parser.add_argument("--metadata", help = mdstr)
    parser.add_argument("--pdf", action = "store_true", help = "Save a {prefix}.pdf file")
    parser.add_argument("--html", action = "store_true", help = "Save a {prefix}.html file")
    parser.add_argument("--plot_features", action = "store_true", help = "Looks in the DF for features to plot. Plots everything on the same plot. Only takes in one dataframe.")
    parser.add_argument("--genome-min-bp", type=float, default=None,
                        help="Minimum genome size (bp). Values <= this are shown as --genome-min-color (grey).")
    parser.add_argument("--genome-max-bp", type=float, default=None,
                        help="Maximum genome size (bp). Values >= this are shown as --genome-max-color (default: cmap endpoint).")
    parser.add_argument("--legend-scale", type=float, default=0.5,
                        help="Scale factor for legend (colorbar) size and font. 1.0 = original; 0.5 = half size.")
    parser.add_argument("--genome-min-color", type=str, default="#DCDEE3",
                        help="Hex color for genome sizes <= genome-min-bp (default grey).")
    parser.add_argument("--genome-max-color", type=str, default="#FF2608",
                        help="Hex color for genome sizes >= genome-max-bp (if not set, uses cmap endpoint color).")
    parser.add_argument("--benedictus", action="store_true",
                        help="Use Benedictus color scheme for genome size panels (overrides genome min/max options).")

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

    # Make sure all the files exist
    if args.directory:
        df_filelist = [os.path.join(args.directory, f) for f in os.listdir(args.directory) if f.endswith(".df")]
        if not df_filelist:
            raise ValueError(f"No .df files found in directory: {args.directory}")
    elif args.filelist:
        df_filelist = args.filelist.split(" ")
        for filepath in df_filelist:
            if not os.path.exists(filepath):
                raise ValueError(f"File does not exist: {filepath}")

    # Check that the metadata file(s) exist(s) if specified. They are space-separated, and the output will be a list of files, even if just one file is specified.
    if args.metadata:
        metadata_files = args.metadata.split(" ")
        for metadata_file in metadata_files:
            if not os.path.exists(metadata_file):
                raise ValueError(f"Metadata file does not exist: {metadata_file}")
        args.metadata = metadata_files

    # just make sure the max isn't less than the min
    if args.genome_min_bp is not None and args.genome_max_bp is not None:
        if args.genome_min_bp >= args.genome_max_bp:
            raise ValueError("genome-min-bp must be < genome-max-bp")

    return args

def generate_df_dict(args):
    """
    Reads .df files from a directory or file list, extracts parameters from filenames,
    and returns a dictionary where keys are (num_neighbors, min_dist) and values are DataFrames.

    Parameters:
    - args: Argument object with `directory` or `filelist` attributes.

    Returns:
    - df_dict: {(num_neighbors, min_dist): pd.DataFrame}
    """
    df_filelist = []
    if args.directory:
        df_filelist = [x for x in os.listdir(args.directory) if x.endswith(".df")]
        df_filelist = [os.path.join(args.directory, f) for f in df_filelist]
    elif args.filelist:
        df_filelist = args.filelist.split(" ")

    df_dict = {}

    for filepath in df_filelist:
        filename = os.path.basename(filepath)

        try:
            # Decapodiformes_215450_without_None.method_phylogenetic.neighbors_50.mind_0.9.missing_large.subchrom.df
            samplename  =       filename.split(".")[0]
            avgmethod   =       filename.split(".method_")[0].split(".")[0]
            miss_size   =       filename.split(".missing_")[1].split(".")[0]
            num_neighbors = int(filename.split(".neighbors_")[1].split(".")[0])
            min_dist    = float(filename.split(".mind_")[1].split(".missing_")[0])
        except (IndexError, ValueError):
            raise ValueError(f"Invalid filename: {filename}")

        # Read DataFrame
        df = pd.read_csv(filepath, sep="\t", index_col=0, header=0)
        results = {"df": df, "samplename": samplename,
                   "filepath": filepath,
                   "num_neighbors": num_neighbors,
                   "min_dist": min_dist,
                   "size": miss_size,
                   "method": avgmethod}
        df_dict[(num_neighbors, min_dist)] = results

    return df_dict

def plot_paramsweep(df_dict, outpdf):
    """
    Makes the plot for the parameter sweep plot when we provide multiple dataframes.

    Uses the df_dict as input.
    """
    # Extract sorted unique values for num_neighbors & min_dist
    num_neighbors_list = sorted(set(k[0] for k in df_dict.keys()))
    min_dist_list = sorted(set(k[1] for k in df_dict.keys()))

    # The rest of the numbers are calculated based on these two
    #     x   0 1 2 3 4
    #  y +-----------------+
    #    |
    #  0 |    o o o o o
    #  1 |    o o o o o
    #  2 |    o o o o o
    #  3 |    o o o o o
    #    |
    #    +-----------------+
    #
    # Determine Figure size
    # setup the plot based on what we know the parameters will be
    # These are the the magic numbers! We only need to adjust the size of each panel and how big the margins will be
    text_size = 10
    panel_width = 1
    margin = 0.25
    panel_height = panel_width
    # the width will have 4 margins, plot, margin, plot... 4 margins
    fig_width  = (margin * 4) + (panel_width * len(min_dist_list))       + (margin * (len(min_dist_list) - 1)) + (margin * 4)
    # the height will have the same thing, but for the number of neighbors
    fig_height = (margin * 4) + (panel_height * len(num_neighbors_list)) + (margin * (len(num_neighbors_list) - 1)) + (margin * 4)
    fig = plt.figure(figsize=(fig_width, fig_height))

    # Determine figure size
    text_size = 10
    panel_width, margin = 1, 0.25
    fig_width = (margin * 4) + (panel_width * len(min_dist_list)) + (margin * (len(min_dist_list) - 1)) + (margin * 4)
    fig_height = (margin * 4) + (panel_width * len(num_neighbors_list)) + (margin * (len(num_neighbors_list) - 1)) + (margin * 4)
    fig = plt.figure(figsize=(fig_width, fig_height))

    # Create axes grid with correct dimensions
    axes = [[None for _ in min_dist_list] for _ in num_neighbors_list]

    # This creates the grid of plots and removes the ticks and spines
    for y_idx, num_neighbors in enumerate(num_neighbors_list):
        for x_idx, min_dist in enumerate(min_dist_list):
            left = (4 * margin) + (x_idx * panel_width) + (x_idx * margin)
            bottom = fig_height - ((4 * margin) + ((y_idx + 1) * panel_width) + (y_idx * margin))
            axes[y_idx][x_idx] = fig.add_axes([
                left / fig_width,
                bottom / fig_height,
                panel_width / fig_width,
                panel_width / fig_height
            ])
            axes[y_idx][x_idx].set_xticks([])
            axes[y_idx][x_idx].set_yticks([])
            for spine in ['top', 'right', 'bottom', 'left']:
                axes[y_idx][x_idx].spines[spine].set_visible(False)

    # Plot data from df_dict
    for (num_neighbors, min_dist), data in df_dict.items():
        y_idx, x_idx = num_neighbors_list.index(num_neighbors), min_dist_list.index(min_dist)
        ax = axes[y_idx][x_idx]
        df = data["df"]

        if df.empty:
            ax.text(0.5, 0.5, "Empty file", fontsize=3, ha='center')
        else:
            ax.scatter(df["UMAP1"], df["UMAP2"], s=0.5, lw=0, alpha=0.5, color=df["color"])
            minval, maxval = df[["UMAP1", "UMAP2"]].min().min(), df[["UMAP1", "UMAP2"]].max().max()
            ax.set_xlim([minval - abs(.05 * minval), 1.05 * maxval])
            ax.set_ylim([minval - abs(.05 * minval), 1.05 * maxval])

        # If we're at the absolute left (first column), add a Y-axis label
        if x_idx == 0:
            ax.yaxis.set_label_position("left")
            ax.set_ylabel(num_neighbors, rotation=0, ha="right", fontsize=text_size)

        # If we're at the absolute top (first row), add an X-axis label
        if y_idx == 0:
            ax.xaxis.set_label_position("top")
            ax.set_xlabel(min_dist, fontsize=text_size)

    # Add titles and labels
    fig.suptitle(f"{data['samplename']}, {data['size']} values for non-colocalized loci,\n {data['method']} method for averaging",
                 fontsize=text_size)

    fig.text(0.5, (fig_height - (margin * 3)) / fig_height, "Min Distance", ha="center", fontsize=text_size)
    fig.text((margin * 2) / fig_width, 0.5, "Number of Neighbors", va="center", rotation="vertical", fontsize=text_size)

    # Add grid dividers to separate plots
    for x_idx in range(1, len(min_dist_list)):
        x1 = ((4 * margin) + (x_idx * panel_width) + (x_idx * margin) - (margin / 2)) / fig_width
        y1, y2 = ((fig_height - (4 * margin)) / fig_height, (4 * margin) / fig_height)
        fig.add_artist(plt.Line2D([x1, x1], [y1, y2], transform=fig.transFigure, color="#BBBBBB"))

    for y_idx in range(1, len(num_neighbors_list)):
        y1 = ((fig_height - ((4 * margin) + (y_idx * panel_width) + (y_idx * margin) - (margin / 2))) / fig_height)
        x1, x2 = (4 * margin) / fig_width, (fig_width - (4 * margin)) / fig_width
        fig.add_artist(plt.Line2D([x1, x2], [y1, y1], transform=fig.transFigure, color="#BBBBBB"))

    # Save and close
    print(f"Saving file to {outpdf}")
    plt.savefig(outpdf)
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

def plot_features(args, outpdf, metadata_df=None, legend_scale=0.5,
                  genome_min_bp=None, genome_max_bp=None,
                  genome_min_color="#DCDEE3", genome_max_color="#FF2608",
                  use_benedictus=False):
    """
    Make a grid of UMAP scatter panels colored by many features.
    Adds colorbars for genome_size, genome_size_log2, genome_size_log10 with:
      - legend_scale to shrink/grow the colorbars and fonts,
      - genome_min_bp/genome_max_bp thresholds (raw bp) that clamp colors:
          <= min -> genome_min_color, >= max -> genome_max_color.
    """
    # ---------- helpers ----------
    def human_readable_bp(n):
        try:
            n = float(n)
        except Exception:
            return str(n)
        if n >= 1e9:
            return f"{n/1e9:.2f} Gb"
        if n >= 1e6:
            return f"{n/1e6:.2f} Mb"
        if n >= 1e3:
            return f"{n/1e3:.1f} kb"
        return f"{int(n)} bp"

    # local fallback if odpf.interpolate_color isn't available
    def _interp_color(x, vmin, vmax, c0="#DCDEE3", c1="#FF2608"):
        if x is None or not np.isfinite(x):
            return c0
        if vmax == vmin:
            t = 0.0
        else:
            t = (float(x) - float(vmin)) / (float(vmax) - float(vmin))
            t = 0.0 if t < 0 else (1.0 if t > 1.0 else t)
        r0, g0, b0 = mcolors.to_rgb(c0)
        r1, g1, b1 = mcolors.to_rgb(c1)
        r = r0 + t * (r1 - r0)
        g = g0 + t * (g1 - g0)
        b = b0 + t * (b1 - b0)
        return mcolors.to_hex((r, g, b))

    # pick an interpolate fn: prefer your odpf implementation if present
    if hasattr(odpf, "interpolate_color"):
        interpolate_color = odpf.interpolate_color
    else:
        interpolate_color = _interp_color

    # ---------- load DF ----------
    df = pd.read_csv(args.filelist, sep="\t", index_col=0)

    if "smallest_protein" in df.columns:
        df = df[~df["smallest_protein"].isna()]

    # Merge metadata if provided (on 'rbh' = index of df)
    if metadata_df is not None:
        metadata_df = metadata_df.set_index("rbh")
        metadata_df.index = metadata_df.index.astype(str).str.strip()
        df.set_index("rbh", inplace=True, drop=False)  # ensure 'rbh' is the index
        df.index = df.index.astype(str).str.strip()
        matched = metadata_df.index.intersection(df.index)
        print(f"Metadata merge: matched {len(matched)} of {len(df)} UMAP RBH entries")
        df = df.join(metadata_df, how="left")  # join by index

    # ---------- columns to plot ----------
    regular_columns_to_plot = [
        "num_scaffolds", "GC_content",
        "genome_size", "genome_size_log2", "genome_size_log10",
        "median_scaffold_length", "mean_scaffold_length", "scaffold_N50",
        "longest_scaffold", "smallest_scaffold", "fraction_Ns",
        "number_of_gaps", "num_proteins", "mean_protein_length",
        "median_protein_length", "longest_protein", "smallest_protein",
        "from_rbh", "frac_ologs", "frac_ologs_sig", "frac_ologs_single"
    ]

    if "genome_size" in df.columns:
        if "genome_size_log2" not in df.columns:
            df["genome_size_log2"] = np.log2(df["genome_size"] + 1)
        if "genome_size_log10" not in df.columns:
            df["genome_size_log10"] = np.log10(df["genome_size"] + 1)

    regular_columns_to_plot = [c for c in regular_columns_to_plot if c in df.columns]
    olog_columns_to_plot = [
        x for x in df.columns
        if x.startswith("frac_ologs_") and x not in ("frac_ologs_sig", "frac_ologs_single")
    ]

    known_umap_cols = {"UMAP1", "UMAP2", "color"}
    metadata_columns = [col for col in df.columns
                        if col not in known_umap_cols and (
                            col.endswith("_color") or
                            (not col.endswith("_color") and f"{col}_color" in df.columns)
                        )]

    all_columns_to_plot = ["color"] + regular_columns_to_plot + olog_columns_to_plot
    all_columns_to_plot = list(dict.fromkeys(all_columns_to_plot))  # drop dups, keep order

    total_num_cols_to_plot = len(all_columns_to_plot) + int(len(metadata_columns) / 2)
    plots_per_row = int(np.ceil(np.sqrt(total_num_cols_to_plot)))
    num_rows = plots_per_row
    num_cols = plots_per_row

    # ---------- figure & axes grid ----------
    margin = 0.25
    panel_width = 1.5
    fig_width = (margin * 4) + (panel_width * num_cols) + (margin * (num_cols - 1)) + (margin * 4)
    fig_height = (margin * 4) + (panel_width * num_rows) + (margin * (num_rows - 1)) + (margin * 4)
    fig = plt.figure(figsize=(fig_width, fig_height))

    axes = [[None for _ in range(num_cols)] for _ in range(num_rows)]
    for ii in range(num_rows):
        for jj in range(num_cols):
            left = (4 * margin) + (jj * panel_width) + (jj * margin)
            bottom = fig_height - ((4 * margin) + ((ii + 1) * panel_width) + (ii * margin))
            ax = fig.add_axes([
                left / fig_width,
                bottom / fig_height,
                panel_width / fig_width,
                panel_width / fig_height
            ])
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ["top", "right", "bottom", "left"]:
                ax.spines[spine].set_visible(False)
            axes[ii][jj] = ax

    fig.suptitle(f"Paramplot for {args.filelist}", fontsize=4)

    # Colormap for genome panels uses your chosen endpoint colors
    custom_cmap = (benedictus_cmap()
                   if use_benedictus
                   else LinearSegmentedColormap.from_list("genome_cmap",
                                                         [genome_min_color, genome_max_color]))

    # pull user params (already in function signature)
    # legend_scale, genome_min_bp, genome_max_bp, genome_min_color, genome_max_color

    # ---------- plot the main columns ----------
    i = 0
    j = 0
    for thiscol in all_columns_to_plot:
        ax = axes[i][j]
        ax.xaxis.set_label_position("top")
        ax.set_xlabel("Clade color" if thiscol == "color" else f"{thiscol}", fontsize=5)

        # Special handling for genome size panels (add colorbars, clamp colors)
        if thiscol in ("genome_size", "genome_size_log2", "genome_size_log10"):
            vals = df[thiscol].to_numpy(dtype=float)

            if use_benedictus:
                use_vmin = np.nanmin(vals)
                use_vmax = np.nanmax(vals)
            else:
                if thiscol == "genome_size":
                    use_vmin = genome_min_bp if genome_min_bp is not None else np.nanmin(vals)
                    use_vmax = genome_max_bp if genome_max_bp is not None else np.nanmax(vals)
                elif thiscol == "genome_size_log2":
                    use_vmin = np.log2(genome_min_bp + 1) if genome_min_bp is not None else np.nanmin(vals)
                    use_vmax = np.log2(genome_max_bp + 1) if genome_max_bp is not None else np.nanmax(vals)
                else:  # log10
                    use_vmin = np.log10(genome_min_bp + 1) if genome_min_bp is not None else np.nanmin(vals)
                    use_vmax = np.log10(genome_max_bp + 1) if genome_max_bp is not None else np.nanmax(vals)

            # safety
            if not np.isfinite(use_vmin):
                use_vmin = np.nanmin(vals)
            if not np.isfinite(use_vmax):
                use_vmax = np.nanmax(vals)
            if use_vmin == use_vmax:
                use_vmax = use_vmin + 1.0

            norm = Normalize(vmin=use_vmin, vmax=use_vmax)
            cmap = custom_cmap

            mapped_rgba = cmap(norm(vals))

            # Convert panel values back to *raw bp* for threshold checks / labels
            if thiscol == "genome_size":
                raw_vals = vals
            elif thiscol == "genome_size_log2":
                raw_vals = (2.0 ** vals) - 1.0
            else:
                raw_vals = (10.0 ** vals) - 1.0

            if not use_benedictus:
                if genome_min_bp is not None:
                    mask_min = np.isfinite(raw_vals) & (raw_vals <= genome_min_bp)
                    if mask_min.any():
                        mapped_rgba[mask_min] = mcolors.to_rgba(genome_min_color)
                if genome_max_bp is not None:
                    mask_max = np.isfinite(raw_vals) & (raw_vals >= genome_max_bp)
                    if mask_max.any():
                        mapped_rgba[mask_max] = mcolors.to_rgba(genome_max_color)

            nan_mask = ~np.isfinite(raw_vals)
            if nan_mask.any():
                mapped_rgba[nan_mask] = mcolors.to_rgba("#DDDDDD")

            ax.scatter(df["UMAP1"], df["UMAP2"], s=0.5, lw=0, alpha=0.5, color=mapped_rgba)

            sm = ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            cbar_fraction = max(0.005, 0.046 * float(legend_scale))
            cbar = fig.colorbar(sm, ax=ax, fraction=cbar_fraction, pad=0.02)
            cbar.ax.tick_params(labelsize=max(1, int(4 * float(legend_scale))))

            ticks = np.linspace(use_vmin, use_vmax, 5)
            cbar.set_ticks(ticks)
            if thiscol == "genome_size":
                labels = [human_readable_bp(t) for t in ticks]
                cbar.set_label("Genome size", fontsize=max(2, int(5 * float(legend_scale))))
            elif thiscol == "genome_size_log2":
                labels = [human_readable_bp((2.0 ** t) - 1.0) for t in ticks]
                cbar.set_label("Genome size (log2)", fontsize=max(2, int(5 * float(legend_scale))))
            else:
                labels = [human_readable_bp((10.0 ** t) - 1.0) for t in ticks]
                cbar.set_label("Genome size (log10)", fontsize=max(2, int(5 * float(legend_scale))))
            cbar.set_ticklabels(labels)

        else:
            # Non-genome-size panels
            coltype = df[thiscol].dtype
            if thiscol == "color":
                colors = list(df[thiscol])
                ax.scatter(df["UMAP1"], df["UMAP2"], s=0.5, lw=0, alpha=0.5, color=colors)
            else:
                # boolean-like in object dtype -> map two fixed colors
                uniques = df[thiscol].unique()
                if ((True in uniques) or (False in uniques)) and (coltype == "object"):
                    colordict = {True: "#074FF7", False: "#FD6117"}
                    colors = [colordict.get(x, "#999999") for x in df[thiscol]]
                    ax.scatter(df["UMAP1"], df["UMAP2"], s=0.5, lw=0, alpha=0.5, color=colors)
                else:
                    # numeric gradient using your interpolation helper
                    arr = pd.to_numeric(df[thiscol], errors="coerce")
                    maxval = np.nanmax(arr.to_numpy(dtype=float))
                    if np.isfinite(maxval) and maxval > 1:
                        colors = [interpolate_color(x, 0, maxval, "#DCDEE3", "#FF2608") for x in arr.fillna(0)]
                    else:
                        colors = [interpolate_color(x, 0, 1, "#DCDEE3", "#FF2608") for x in arr.fillna(0)]
                    ax.scatter(df["UMAP1"], df["UMAP2"], s=0.5, lw=0, alpha=0.5, color=colors)

        # advance grid position
        j += 1
        if j == num_cols:
            j = 0
            i += 1
            if i >= num_rows:
                break  # safety

    # ---------- plot metadata columns (pre-colored) ----------
    for thiscol in metadata_columns:
        if thiscol.endswith("_color"):
            continue
        if i >= num_rows:
            break
        ax = axes[i][j]
        ax.xaxis.set_label_position("top")
        ax.set_xlabel(f"{thiscol}", fontsize=5)

        color_col = f"{thiscol}_color"
        if color_col not in df.columns:
            raise ValueError(f"Missing expected color column {color_col} for metadata.")
        colors = df[color_col]
        ax.scatter(df["UMAP1"], df["UMAP2"], s=0.5, lw=0, alpha=0.5, color=colors)

        # add legend for categorical metadata
        if pd.api.types.is_object_dtype(df[thiscol]) or pd.api.types.is_categorical_dtype(df[thiscol]):
            unique_vals = df[[thiscol, color_col]].dropna().drop_duplicates()
            handles = [
                plt.Line2D([0], [0], marker="o", color="none",
                           label=str(row[thiscol]),
                           markerfacecolor=row[color_col],
                           markersize=2, markeredgewidth=0, markeredgecolor="none")
                for _, row in unique_vals.iterrows()
            ]
            ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
                      fontsize=2, frameon=False)

        # advance grid
        j += 1
        if j == num_cols:
            j = 0
            i += 1
            if i >= num_rows:
                break

    # ---------- aspect ratio ----------
    for row in axes:
        for col in row:
            if col is not None:
                col.set_aspect("equal", adjustable="box")

    # ---------- figure grid lines ----------
    xs = []
    for jj in range(1, num_cols):
        left_current = (4 * margin) + (jj * panel_width) + ((jj - 1) * margin)
        xs.append(left_current / fig_width)
    for x in xs:
        line = plt.Line2D([x, x], [0, 1], transform=fig.transFigure, color="#BBBBBB", lw=0.5)
        fig.add_artist(line)

    ys = []
    for ii in range(1, num_rows):
        bottom_current = (4 * margin) + (ii * panel_width) + ((ii - 1) * margin)
        y = 1 - (bottom_current / fig_height)
        ys.append(y)
    for y in ys:
        line = plt.Line2D([0, 1], [y, y], transform=fig.transFigure, color="#BBBBBB", lw=0.5)
        fig.add_artist(line)

    # ---------- save ----------
    print(f"saving the file to {outpdf}")
    plt.savefig(outpdf)
    plt.close(fig)


#def plot_features(args, outpdf, metadata_df=None):
#    """
#    This makes a plot of the features of a single dataframe.
#    All of the plots will be plotted along the axes of the UMAP1 and UMAP2.
#    """
#    # get the dataframe to load in
#    df = pd.read_csv(args.filelist, sep="\t", index_col=0)
#    # remove the rows that have a value of NaN in the "smallest_protein" column
#    if "smallest_protein" in df.columns:
#        df = df[~df["smallest_protein"].isna()]
#
#    # Merge metadata if provided (on 'rbh' = index of df)
#    if metadata_df is not None:
#        metadata_df = metadata_df.set_index("rbh")
#        metadata_df.index = metadata_df.index.astype(str).str.strip()
#        df.set_index("rbh", inplace=True, drop=False)  # ensure 'rbh' is the index
#        df.index = df.index.astype(str).str.strip()
#
#        matched = metadata_df.index.intersection(df.index)
#        print(f"Metadata merge: matched {len(matched)} of {len(df)} UMAP RBH entries")
#
#        df = df.join(metadata_df, how="left")  # join by index
#
#    print(df.columns)
#    print(df)
#    # the columns we want to annotate are defined in the AnnotateSampleDf pipeline
#    #  - things from the genome fasta file
#    #    - num_scaffolds
#    #    - GC_content
#    #    - genome_size
#    #    - median_scaffold_length
#    #    - mean_scaffold_length
#    #    - scaffold_N50
#    #    - longest_scaffold
#    #    - smallest_scaffold
#    #    - fraction_Ns
#    #    - number_of_gaps
#    #  - things from the protein file:
#    #    - num_proteins
#    #    - mean_protein_length
#    #    - median_protein_length
#    #    - longest_protein
#    #    - smallest_protein
#    #    - from_rbh
#    #  - things from the rbh file:
#    #    - frac_ologs:           The fraction of genes of ANY ALG that are present at all in the rbh file. len(rbhdf) / total_genes_ALGs
#    #    - frac_ologs_sig:       The fraction of genes of ANY ALG that are significantly on any chromosome, as defined by whole_FET
#    #    - frac_ologs_single:    The fraction of genes of ANY ALG that are significantly on the largest chromosome, as defined by whole_FET
#    #    - frac_ologs_{ALGNAME}: The fraction of genes of INDIVIDUAL ALGs that are significantly on any chromosome
#
#    regular_columns_to_plot = ["num_scaffolds", "GC_content", "genome_size", "genome_size_log2", "genome_size_log10",
#                               "median_scaffold_length", "mean_scaffold_length", "scaffold_N50", "longest_scaffold",
#                               "smallest_scaffold", "fraction_Ns", "number_of_gaps", "num_proteins", "mean_protein_length",
#                               "median_protein_length", "longest_protein", "smallest_protein", "from_rbh",
#                               "frac_ologs", "frac_ologs_sig", "frac_ologs_single"]
#    # if there is no log2 genome size plot, make one
#    if "genome_size_log2" not in df.columns:
#        df["genome_size_log2"]  = np.log2(df["genome_size"] + 1)
#    if "genome_size_log10" not in df.columns:
#        df["genome_size_log10"] = np.log10(df["genome_size"] + 1)
#    # If we're plotting metadata, we might not actually have these columns in the dataframe.
#    regular_columns_to_plot = [col for col in regular_columns_to_plot if col in df.columns]
#    olog_columns_to_plot = [x for x in df.columns if x.startswith("frac_ologs_") and x not in ("frac_ologs_sig", "frac_ologs_single")]
#
#    known_umap_cols = {"UMAP1", "UMAP2", "color"}
#    metadata_columns = [col for col in df.columns
#                        if col not in known_umap_cols and (
#                            col.endswith("_color") or 
#                            (not col.endswith("_color") and f"{col}_color" in df.columns)
#                        )
#                       ]
#    all_columns_to_plot = ["color"] + regular_columns_to_plot + olog_columns_to_plot
#    all_columns_to_plot = list(dict.fromkeys(all_columns_to_plot))  # remove duplicates while preserving order
#    total_num_cols_to_plot = len(all_columns_to_plot) + int(len(metadata_columns)/2) # divide by 2 because we include both the column and its colors
#    # figure out what number we should pick to make the shape closest to a square, when plotting N x N plots
#    plots_per_row = int(np.ceil(np.sqrt(total_num_cols_to_plot)))
#
#    num_rows = plots_per_row
#    num_cols = plots_per_row
#    # make a grid of squares to plot each of these on
#    # reduce the space between all the plots
#    # make the figure size such that all the plots are square
#
#    # make the axes
#    margin = 0.25
#    panel_width = 1.5
#    fig_width = (margin * 4) + (panel_width * num_cols) + (margin * (num_cols - 1)) + (margin * 4)
#    fig_height = (margin * 4) + (panel_width * num_rows) + (margin * (num_rows - 1)) + (margin * 4)
#    fig = plt.figure(figsize=(fig_width, fig_height))
#
#    axes = [[None for _ in range(num_cols)] for _ in range(num_rows)]
#    for i in range(num_rows):
#        for j in range(num_cols):
#            left = (4 * margin) + (j * panel_width) + (j * margin)
#            bottom = fig_height - ((4 * margin) + ((i + 1) * panel_width) + (i * margin))
#            ax = fig.add_axes([
#                left / fig_width,
#                bottom / fig_height,
#                panel_width / fig_width,
#                panel_width / fig_height
#            ])
#            ax.set_xticks([])
#            ax.set_yticks([])
#            for spine in ['top', 'right', 'bottom', 'left']:
#                ax.spines[spine].set_visible(False)
#            axes[i][j] = ax
#
#    # set the title as samplename and the whether it is small or large
#    fig.suptitle(f"Paramplot for {args.filelist}", fontsize = 4)
#    # set absolute left label as the number of neighbors
#    #fig.text(0.06, 0.5, 'Number of Neighbors', va='center', rotation='vertical')
#    #fig.text(0.5, 0.92, 'Min Distance', ha='center')
#    i = 0
#    j = 0
#    for thiscol in all_columns_to_plot:
#        # for the left-most plot in each row, set the ylabel as the number of neighbors
#        #axes[i][j].set_ylabel(f"{row}", rotation=0, ha='right')
#        #for the top-most plot in each column, set the xlabel as the min_dist
#        # put the xlabel on the top
#        axes[i][j].xaxis.set_label_position('top')
#        if thiscol == "color":
#            axes[i][j].set_xlabel(f"Clade color", fontsize = 5)
#        else:
#            axes[i][j].set_xlabel(f"{thiscol}", fontsize = 5)
#        # figure out the type of the column
#        coltype = df[thiscol].dtype
#        #if the type of the column is an object, then plot True as #074FF7 and False as #FD6117
#        # If True/False not in that column , then randomly select a color
#        colors = []
#        print('coltype of {} is {}'.format(thiscol, coltype))
#        if thiscol == "color":
#            colors = list(df[thiscol])
#        else:
#            if ((True in df[thiscol].unique()) or (False in df[thiscol].unique())) and (coltype == "object"):
#                # print the counts of True and False
#                print(df[thiscol].value_counts())
#                colordict = {True: "#074FF7", False: "#FD6117"}
#                colors = [colordict[x] for x in df[thiscol]]
#            else:
#                # make an object to interpolate color between #DCDEE3 and #FF2608
#                # if the max value is above 1, then we will interpolate between 0 and the max of the column
#                if df[thiscol].max() > 1:
#                    colors = [interpolate_color(x, 0, df[thiscol].max(), "#DCDEE3", "#FF2608") for x in df[thiscol]]
#                else:
#                    # otherwise go between 0 and 1
#                    colors = [interpolate_color(x, 0, 1, "#DCDEE3", "#FF2608") for x in df[thiscol]]
#
#        # get the df file for this row and column from the dffile
#        axes[i][j].scatter(df["UMAP1"], df["UMAP2"],
#                         s=0.5, lw = 0, alpha=0.5,
#                         color=colors)
#        # iterate i and j
#        j += 1
#        if j == num_cols:
#            j = 0
#            i += 1
#
#    # Now we plot the metadata columns. This is less complicated
#    #  since we already calculated the colors for everything
#    for thiscol in metadata_columns:
#        if thiscol.endswith("_color"):
#            # if the column is not a color column, then we will plot the column and its color
#            continue
#        # get the color column
#        color_col = f"{thiscol}_color"
#        if color_col not in df.columns:
#            # We should have generated this in the other function, so raise an error here
#            raise ValueError(f"The column {color_col} is not present in the dataframe. This should have been generated in the parse_metadata_dfs() function.")
#        axes[i][j].xaxis.set_label_position('top')
#        axes[i][j].set_xlabel(f"{thiscol}", fontsize = 5)
#        colors = df[color_col]
#        axes[i][j].scatter(df["UMAP1"], df["UMAP2"],
#                            s=0.5, lw = 0, alpha=0.5,
#                            color=colors)
#        # If the column is categorical, make a legend
#        if pd.api.types.is_object_dtype(df[thiscol]) or pd.api.types.is_categorical_dtype(df[thiscol]):
#            unique_vals = df[[thiscol, color_col]].dropna().drop_duplicates()
#
#            handles = [
#                plt.Line2D(
#                    [0], [0],
#                    marker='o',
#                    color='none',
#                    label=str(row[thiscol]),
#                    markerfacecolor=row[color_col],
#                    markersize=2,
#                    markeredgewidth=0,
#                    markeredgecolor='none'
#                )
#                for _, row in unique_vals.iterrows()
#            ]
#
#            axes[i][j].legend(
#                handles=handles,
#                loc='center left',
#                bbox_to_anchor=(1.0, 0.5),
#                fontsize=2,
#                frameon=False
#            )
#
#        # iterate i and j
#        j += 1
#        if j == num_cols:
#            j = 0
#            i += 1
#
#    # make sure that the aspect ratio for all of these is the same
#    # Iterate over each axis and set aspect ratio to 'equal'
#    for row in axes:
#        for col in row:
#            col.set_aspect('equal', adjustable='box')
#
#    # Now make vertical and horizontal lines to separate the plots. Make them medium gray.
#    xs = []
#    for j in range(1, num_cols):
#        left_current = (4 * margin) + (j * panel_width) + ((j - 1) * margin)
#        xs.append(left_current / fig_width)
#    for x in xs:
#        line = plt.Line2D([x, x], [0, 1], transform=fig.transFigure, color="#BBBBBB", lw=0.5)
#        fig.add_artist(line)
#
#    ys = []
#    for i in range(1, num_rows):
#        bottom_current = (4 * margin) + (i * panel_width) + ((i - 1) * margin)
#        y = 1 - (bottom_current / fig_height)
#        ys.append(y)
#    for y in ys:
#        line = plt.Line2D([0, 1], [y, y], transform=fig.transFigure, color="#BBBBBB", lw=0.5)
#        fig.add_artist(line)
#
#    # save the figure as f"{samplename}_{name}.pdf"
#    # name is just the size, small or large
#    print("saving the file to {}".format(outpdf))
#    plt.savefig(outpdf)
#    # close the figure
#    plt.close(fig)

def generate_umap_grid_bokeh(df_dict, output_html):
    """
    Takes a dictionary of UMAP DataFrames and generates a Bokeh grid plot.
    - Rows represent different `num_neighbors` values.
    - Columns represent different `min_dist` values.
    - Saves as an interactive HTML file.

    Parameters:
    - df_dict: Dict where keys are (num_neighbors, min_dist) and values are dictionaries with:
      - "df": DataFrame with UMAP coordinates
    - output_html: Path to the output HTML file.
    """

    # Extract unique num_neighbors (rows) and min_dist (columns) values
    num_neighbors_list = sorted(set(k[0] for k in df_dict.keys()))
    min_dist_list = sorted(set(k[1] for k in df_dict.keys()))

    # Adjust plot sizing to reduce overall width and make each plot smaller
    plot_width = 150  # Reduced plot width
    plot_height = 150  # Reduced plot height

    # Retrieve general metadata for the title (from any entry)
    sample_info = next(iter(df_dict.values()))  # Get any sample metadata
    title_text = f"{sample_info['samplename']}, {sample_info['size']} values for non-colocalized loci, {sample_info['method']} method for averaging"

    # Create a title for the entire plot
    title_div = Div(text=f"<h3>{title_text}</h3>", width=plot_width * len(min_dist_list), height=40, styles={"text-align": "center"})

    # Create a dictionary to store plots in a grid layout
    plot_grid = [[None for _ in min_dist_list] for _ in num_neighbors_list]

    # Generate scatter plots for each dataset
    for (num_neighbors, min_dist), data in df_dict.items():
        df = data["df"]  # Extract the DataFrame

        # Ensure only valid DataFrame columns are passed to Bokeh
        valid_columns = ["UMAP1", "UMAP2", "color"]  # Keep only numeric/iterable columns
        df_filtered = df[valid_columns] if set(valid_columns).issubset(df.columns) else df

        source = ColumnDataSource(df_filtered)

        # Create Bokeh scatter plot with adjusted transparency, smaller dots, and no grid lines
        p = figure(width=plot_width, height=plot_height, tools="", toolbar_location=None)
        p.scatter(x="UMAP1", y="UMAP2", source=source, size=1, color="color", alpha=0.3, line_color=None)  # Smaller, transparent dots with no outlines

        # Remove grid lines, axis labels, and ticks
        p.xgrid.visible = False
        p.ygrid.visible = False
        p.outline_line_color = None  # Removes outer plot border
        p.xaxis.visible = False
        p.yaxis.visible = False

        # Determine grid position
        row_idx = num_neighbors_list.index(num_neighbors)
        col_idx = min_dist_list.index(min_dist)
        plot_grid[row_idx][col_idx] = p

    # Create centered labels for rows (num_neighbors) and columns (min_dist)
    row_labels = [Div(text=f"<b>{n}</b>", width=30, height=plot_height, styles={"text-align": "center", "display": "flex", "align-items": "center", "justify-content": "center"}) for n in num_neighbors_list]
    col_labels = [Div(text=f"<b>{d}</b>", width=plot_width, height=30, styles={"text-align": "center"}) for d in min_dist_list]

    # Arrange plots into a grid layout with labels
    full_grid = [[Div(text="", width=30, height=30)] + col_labels]  # Top row with column labels
    for row_label, plots in zip(row_labels, plot_grid):
        full_grid.append([row_label] + plots)

    # Full layout including title
    layout = column(title_div, gridplot(full_grid))

    # Output to an HTML file
    # Print confirmation message
    print(f"Saving Bokeh grid plot to {output_html}")
    output_file(output_html)
    save(layout)

def parse_metadata_dfs(df_filelist: list):
    """
    This function reads in a series of files that contain at least two columns:
    -                      rbh: The rbh identifier, which is the same as the index of the main dataframe.
    -       <your_column_name>: The column that you want to use to annotate the plot. This
                                  could be a categorical or numeric column. If there is no additional
                                  column called <your_column_name>_color, then the colors will be assigned
                                  a gradient if numeric, or a random color if the column is categorical.
    - <your_column_name>_color: The column that contains the colors to use for the plot.
                                  If this column is not present, then the colors will be assigned
                                  a gradient if numeric, or a random color if the column is categorical.

    Notes:
      - The argument df_filelist is a list of files that contain the metadata. This should be a list even if there is only one file.
      - There can be multiple columns with <your_column_name>_color, and they will be used to color the points in the plot. These
        will get merged together into a single dataframe, with a series of columns that contain the original data, and corresponding columns
        that contain the colors to use for the plot.
    """
    # first enforce that the type of df_filelist is a list
    if not isinstance(df_filelist, list):
        raise ValueError("The df_filelist argument must be a list of files.")
    # ensure that all of the files exist
    for df_file in df_filelist:
        if not os.path.exists(df_file):
            raise ValueError(f"Metadata file does not exist: {df_file}")
    list_of_dfs = [] # the dfs pre-merge will be added into here
    for df_file in df_filelist:
        # read in the dataframe
        df = pd.read_csv(df_file, sep="\t", header=0)
        # check that the rbh column is present
        if "rbh" not in df.columns:
            raise ValueError(f"The metadata file {df_file} does not contain a 'rbh' column.")
        # check that there is at least one column that is not rbh
        if len(df.columns) < 2:
            raise ValueError(f"The metadata file {df_file} does not contain any columns other than 'rbh'.")
        # Raise an error if there is a column called "rbh_color", this conflicts with the column we will merge against
        if "rbh_color" in df.columns:
            raise ValueError(f"The metadata file {df_file} contains a column called 'rbh_color', which conflicts with the column we will merge against.")
        # Get a list of the columns ending in _color. If there is not a corresponding column without _color, then we will raise an error.
        color_columns = [col for col in df.columns if col.endswith("_color")]
        for color_col in color_columns:
            non_color_col = color_col[:-6]  # Remove the '_color' suffix
            if non_color_col not in df.columns:
                raise ValueError(f"The metadata file {df_file} contains a column called '{color_col}', but does not contain a corresponding column called '{non_color_col}'.")
        # For all the color columns, ensure that the colors are valid hex colors. Tell the user that we only allow hex colors and not RGB or other formats.
        for color_col in color_columns:
            if not all(df[color_col].apply(lambda x: isinstance(x, str) and x.startswith("#") and len(x) == 7)):
                raise ValueError(f"The metadata file {df_file} contains a column called '{color_col}', but it does not contain valid hex colors. We only allow hex colors in the format #RRGGBB.")

        # NOW WE COLOR THE COLUMNS IF THERE IS NO _color COLUMN
        # For all the columns that do not have a _color column, we will assign a color based on the values in the column.
        non_color_non_rbh_columns = [col for col in df.columns if col != "rbh" and not col.endswith("_color")]
        for thiscol in non_color_non_rbh_columns:
            # If the column is numeric, we will assign a gradient color
            if pd.api.types.is_numeric_dtype(df[thiscol]):
                # Get the min and max values of the column
                vmin = df[thiscol].min()
                vmax = df[thiscol].max()
                # Assign a color based on the values in the column - this is a bit of a magic "color" since we use it in the other function, but leaving it here for now
                df[f"{thiscol}_color"] = df[thiscol].apply(lambda x: interpolate_color(x, vmin, vmax, "#DCDEE3", "#FF2608"))
            else:
                unique_values = df[thiscol].unique()
                tcolor = "#074FF7"
                fcolor = "#FD6117"

                # there there is only one value, then raise a warning. There isn't a point for the user to color a column with only one value. Just color everything black.
                if len(unique_values) == 1:
                    print(f"Warning: The column '{thiscol}' in the metadata file '{df_file}' contains only one unique value. All points will be colored black.")
                    df[f"{thiscol}_color"] = "#000000"
                elif len(unique_values) == 2:
                    # When there are only two things to color, blue/orange is a good contrasting choice that is colorblind-friendly.
                    # If the column contains True/Fales values, color them with binary. In fact, if the column is boolean at all, just assign the two True/False colors randomly
                    if set(unique_values).issubset({True, False}):
                        # If the column is boolean, we will assign a color based on the True/False values
                        df[f"{thiscol}_color"] = df[thiscol].map({True: tcolor, False: fcolor})
                    else:
                        # If there are only two unique values, we will assign tcolor or fcolor randomly
                        df[f"{thiscol}_color"] = df[thiscol].apply(lambda x: tcolor if x == unique_values[0] else fcolor)
                else:
                    # This column is not numeric and not boolean, so we will assign a random color to each unique value
                    # Generate a random color for each unique value. Use generate_random_color()
                    unique_values = pd.Series(df[thiscol].dropna().unique())
                    distinct_colors = generate_distinct_colors(len(unique_values))
                    color_map = dict(zip(unique_values, distinct_colors))
                    df[f"{thiscol}_color"] = df[thiscol].map(color_map).fillna("#000000")
        # append each of these to the list of dataframes
        list_of_dfs.append(df)
    # Now we do some checks before merging to ensure their safety
    # Check for column name conflicts (excluding 'rbh')
    all_columns = [c for df in list_of_dfs for c in df.columns if c != "rbh"]
    if len(all_columns) != len(set(all_columns)):
        dupes = [c for c in set(all_columns) if all_columns.count(c) > 1]
        raise ValueError(f"Conflicting column names across metadata files: {dupes}. "
                         f"Each metadata file must have unique annotation column names.")

    # Now we merge all the dataframes together on the "rbh" column
    # Set index on rbh for merging
    for i in range(len(list_of_dfs)):
        list_of_dfs[i] = list_of_dfs[i].set_index("rbh")

    # Merge all metadata dataframes on 'rbh'
    from functools import reduce
    merged_df = reduce(lambda left, right: left.join(right, how="outer"), list_of_dfs)

    # Restore rbh as a column (optional, depending on your downstream needs)
    merged_df.reset_index(inplace=True)

    return merged_df

def main():
    odpf.format_matplotlib()
    args = parse_args()

    if args.plot_features:
        metadatadf = parse_metadata_dfs(args.metadata) if args.metadata else None
        print(f"Metadata DataFrame:\n {metadatadf.head() if metadatadf is not None else 'No metadata provided'}")
        outpdf = args.prefix + ".features.pdf"
        plot_features(args, outpdf, metadata_df = metadatadf, legend_scale = args.legend_scale,
                      genome_min_bp = args.genome_min_bp, genome_min_color = args.genome_min_color,
                      genome_max_bp = args.genome_max_bp, genome_max_color = args.genome_max_color,
                      use_benedictus = args.benedictus)
    else:
        # In this scenario, we will simply plot the dataframes in a grid based on the parameters.
        df_dict = generate_df_dict(args) # a special function that reads in the dataframes and extracts the parameters from the filenames
        ## print the keys of the df_dict
        #print(df_dict.keys())

        # plot the parameter sweep plot, as a pdf or html or both
        if args.pdf:
            outpdf = args.prefix + ".pdf"
            plot_paramsweep(df_dict, outpdf)
        if args.html:
            outhtml = args.prefix + ".html"
            generate_umap_grid_bokeh(df_dict, outhtml)

if __name__ == "__main__":
    main()