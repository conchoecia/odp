#!/usr/bin/env python
"""
# Take in a list of datafraes from samples and constructs a comparison of the UMAP plots.
"""
import argparse
import numpy as np
import os
import pandas as pd
import random
import re
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

# Phylum dictionary for phyla-based plotting (sorted alphabetically)
phylum_to_taxid = {
    "Annelida":         {"taxid": 6340,    "index": 1,  "size": "macroscopic", "length": "0.5mm-3m"},
    "Arthropoda":       {"taxid": 6656,    "index": 2,  "size": "macroscopic", "length": "0.08mm-1m"},
    "Brachiopoda":      {"taxid": 7568,    "index": 3,  "size": "macroscopic", "length": "1mm-100mm"},
    "Bryozoa":          {"taxid": 10205,   "index": 4,  "size": "microscopic", "length": "0.5mm"},
    "Chaetognatha":     {"taxid": 10229,   "index": 5,  "size": "macroscopic", "length": "2mm-120mm"},
    "Chordata":         {"taxid": 7711,    "index": 6,  "size": "macroscopic", "length": "0.5mm-30m"},
    "Cnidaria":         {"taxid": 6073,    "index": 7,  "size": "macroscopic", "length": "0.006mm-2m"},
    "Ctenophora":       {"taxid": 10197,   "index": 8,  "size": "macroscopic", "length": "10mm-15cm"},
    "Cycliophora":      {"taxid": 69815,   "index": 9,  "size": "microscopic", "length": "0.1-0.5mm"},
    "Echinodermata":    {"taxid": 7586,    "index": 10, "size": "macroscopic", "length": "4mm-3m"},
    "Entoprocta":       {"taxid": 43120,   "index": 11, "size": "macroscopic", "length": "0.1mm-7mm"},
    "Gastrotricha":     {"taxid": 33313,   "index": 12, "size": "microscopic", "length": "0.06mm-3mm"},
    "Gnathostomulida":  {"taxid": 66780,   "index": 13, "size": "microscopic", "length": "0.5mm-1mm"},
    "Hemichordata":     {"taxid": 10219,   "index": 14, "size": "macroscopic", "length": "0.6mm-2.5m"},
    "Kinorhyncha":      {"taxid": 51516,   "index": 15, "size": "microscopic", "length": "0.1mm-1mm"},
    "Loricifera":       {"taxid": 310840,  "index": 16, "size": "microscopic", "length": "0.1mm-1mm"},
    "Mollusca":         {"taxid": 6447,    "index": 17, "size": "macroscopic", "length": "0.1mm-1m"},
    "Nematoda":         {"taxid": 6231,    "index": 18, "size": "microscopic", "length": "1mm-7mm"},
    "Nematomorpha":     {"taxid": 33310,   "index": 19, "size": "macroscopic", "length": "5cm-10cm"},
    "Nemertea":         {"taxid": 6217,    "index": 20, "size": "macroscopic", "length": "3mm-54m"},
    "Onychophora":      {"taxid": 27563,   "index": 21, "size": "macroscopic", "length": "0.1cm-22cm"},
    "Orthonectida":     {"taxid": 33209,   "index": 22, "size": "microscopic", "length": "0.035mm-0.3mm"},
    "Phoronida":        {"taxid": 120557,  "index": 23, "size": "macroscopic", "length": "1mm-50cm"},
    "Placozoa":         {"taxid": 10226,   "index": 24, "size": "microscopic", "length": "0.1mm-0.75mm"},
    "Platyhelminthes":  {"taxid": 6157,    "index": 25, "size": "macroscopic", "length": "1mm-60cm"},
    "Porifera":         {"taxid": 6040,    "index": 26, "size": "macroscopic", "length": "10mm-3.5m"},
    "Priapulida":       {"taxid": 33467,   "index": 27, "size": "macroscopic", "length": "0.2mm-39cm"},
    "Rotifera":         {"taxid": 10190,   "index": 28, "size": "microscopic", "length": "0.1mm-2mm"},
    "Sipuncula":        {"taxid": 6433,    "index": 29, "size": "macroscopic", "length": "0.5cm-10cm"},
    "Tardigrada":       {"taxid": 42241,   "index": 30, "size": "microscopic", "length": "0.1mm-1mm"},
    "Xenacoelomorpha":  {"taxid": 1312402, "index": 31, "size": "macroscopic", "length": "2mm-50cm"}
}


def return_kingdom_full_sort_order():
    """Return a list of the sort order for the taxonomic rankings."""
    return ["superkingdom",
            "kingdom",
            "subkingdom",
            "infrakingdom",
            "superphylum",
            "phylum",
            "subphylum",
            "infraphylum",
            "superclass",
            "class",
            "subclass",
            "infraclass",
            "parvclass",
            "cohort",
            "subcohort",
            "superorder",
            "order",
            "suborder",
            "infraorder",
            "parvorder",
            "superfamily",
            "family",
            "subfamily",
            "tribe",
            "subtribe",
            "genus",
            "subgenus",
            "section",
            "subsection",
            "series",
            "subseries",
            "species group",
            "species subgroup",
            "species",
            "subspecies",
            "allsamples"]

# Build a fast lookup once
_CANONICAL_RANKS = return_kingdom_full_sort_order()
# token form: spaces->underscore, lowercase
_RANK_TOKEN_MAP = {
    rk.lower().replace(" ", "_"): rk
    for rk in _CANONICAL_RANKS
}

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
                    still respecting genome size min/max thresholds but ignoring
                    the custom min/max colors.
      --metadata: space-separated list of metadata files to join against the main dataframe. This will be type list [str] of filenames.
      --phylo-map
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
    parser.add_argument("--plot-phyla", action = "store_true", help = "Plot UMAP with each phylum highlighted individually. Requires a 'phylum' column in the dataframe.")
    parser.add_argument("--phyla-rotation", type=str, default="maximize_square",
                        choices=["maximize_square", "minimize_vertical"],
                        help="Rotation strategy for --plot-phyla: 'maximize_square' fits data optimally in square panels, 'minimize_vertical' minimizes vertical extent for condensed figures (default: maximize_square)")
    parser.add_argument("--phyla-clean-output", action="store_true",
                        help="For --plot-phyla, also output a clean version without labels, vectors, or grid lines (saved as {prefix}_clean.pdf)")
    parser.add_argument("--phyla-order", type=str, default=None,
                        help="Space-delimited list of phyla names to plot in custom order (e.g., 'Chordata Arthropoda Mollusca'). Only specified phyla will be plotted. The 2x2 'All Phyla' panel will still be shown.")
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
    parser.add_argument("--threecolor", action="store_true",
                        help="Use Benedictus three-color scheme for genome size panels (ignores custom min/max colors but still applies genome size thresholds).")
                        # benedictus is from here: https://emilhvitfeldt.github.io/r-color-palettes/discrete/MetBrewer/Benedictus/index.html
                        #TODO This should also be implemented with other column types, like ALG % retention
    parser.add_argument("--phylolist", nargs= "+",
                        help=("space-separated .df files. Rank is inferred from each filename. "
                          "Enables a vertical phylo-resampling grid (rows=ranks, "
                          "cols=(n_neighbors,min_dist) inferred from filenames)."))
    parser.add_argument("--num-cols", type=int, default=None,
                        help="Number of columns for --plot_features or --plot-phyla grid layout. If not specified, uses sqrt of total panels (square grid).")

    args = parser.parse_args()

    # --- normalize --phylolist to a list of paths ---
    files = args.phylolist
    if files is None:
        files = []
    elif isinstance(files, list) and len(files) == 1 and (" " in files[0] or "\n" in files[0]):
        # User provided one big quoted string → split it into tokens
        import shlex
        files = [p for p in shlex.split(files[0]) if p]
    # overwrite with standardized list
    args.phylolist = files

    # Make sure that both directory and filelist are not specified.
    # If they are both specified, we don't know which one to use.
    if args.directory and args.filelist:
        raise ValueError("Both directory and filelist are specified. We don't know which one to use. Please just use one.")

    # conflicts
    if args.phylolist and (args.directory or args.filelist):
        raise ValueError("--phylolist cannot be used with --directory/--filelist.")

    # existence check (optional but nice)
    for f in args.phylolist:
        if not os.path.exists(f):
            raise ValueError(f"File does not exist: {f}")

    # If we have turned on the plot_features flag, then we need to make sure that we only have one file in the filelist.
    if args.plot_features:
        if args.filelist:
            if len(args.filelist.split(" ")) > 1:
                raise ValueError("You have turned on the plot_features flag, but you have more than one file in the filelist. We can only plot one file at a time with this flag.")
        elif args.directory:
            if len([x for x in os.listdir(args.directory) if x.endswith(".df")]) > 1:
                raise ValueError("You have turned on the plot_features flag, but you have more than one file in the directory. We can only plot one file at a time with this flag.")

    # If we have turned on the plot_phyla flag, then we need to make sure that we only have one file in the filelist.
    if args.plot_phyla:
        if args.filelist:
            if len(args.filelist.split(" ")) > 1:
                raise ValueError("You have turned on the plot_phyla flag, but you have more than one file in the filelist. We can only plot one file at a time with this flag.")
        elif args.directory:
            if len([x for x in os.listdir(args.directory) if x.endswith(".df")]) > 1:
                raise ValueError("You have turned on the plot_phyla flag, but you have more than one file in the directory. We can only plot one file at a time with this flag.")

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

    args.benedictus = args.threecolor

    return args


def generate_df_dict(args):
    if args.directory:
        df_filelist = [os.path.join(args.directory, f)
                       for f in os.listdir(args.directory) if f.endswith(".df")]
    else:
        df_filelist = args.filelist.split(" ")

    df_dict = {}

    for filepath in df_filelist:
        filename = os.path.basename(filepath)

        # Defaults (optional fields)
        samplename = filename.split(".neighbors_")[0] if ".neighbors_" in filename else filename.split(".")[0]
        avgmethod  = None
        miss_size  = None
        metric     = None

        # Optional: method
        m = re.search(r"\.method_([^.]+)", filename)
        if m: avgmethod = m.group(1)

        # neighbors (required)
        m = re.search(r"\.neighbors_(\d+)", filename)
        if not m:
            raise ValueError(f"Invalid filename (neighbors): {filename}")
        num_neighbors = int(m.group(1))

        # min_dist (required)
        m = re.search(r"\.mind_([0-9]*\.?[0-9]+)", filename)
        if not m:
            raise ValueError(f"Invalid filename (mind): {filename}")
        min_dist = float(m.group(1))

        # Optional: missing size
        m = re.search(r"\.missing_([^.]+)", filename)
        if m: miss_size = m.group(1)

        # Optional: metric before .df (euclidean/cosine/etc.)
        m = re.search(r"\.(euclidean|cosine|manhattan|chebyshev|minkowski)\.df$", filename)
        if m: metric = m.group(1)

        # Load DF
        df = pd.read_csv(filepath, sep="\t", index_col=0, header=0)

        results = {
            "df": df,
            "samplename": samplename,
            "filepath": filepath,
            "num_neighbors": num_neighbors,
            "min_dist": min_dist,
            "size": miss_size,          # may be None
            "method": avgmethod,        # may be None
            "metric": metric            # may be None
        }
        df_dict[(num_neighbors, min_dist)] = results

    return df_dict

def set_square_limits(ax, x, y, q=(0.002, 0.998), pad=0.03):
    """
    Set per-axis limits using optional quantile clipping, then expand to a
    square view with a small padding. Keeps aspect='equal'.
    """
    import numpy as np
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)

    if q is not None:
        x0, x1 = np.quantile(x, q); y0, y1 = np.quantile(y, q)
    else:
        x0, x1 = x.min(), x.max(); y0, y1 = y.min(), y.max()

    if x0 == x1: x0 -= 0.5; x1 += 0.5
    if y0 == y1: y0 -= 0.5; y1 += 0.5

    px = (x1 - x0) * pad
    py = (y1 - y0) * pad
    cx, cy = 0.5*(x0 + x1), 0.5*(y0 + y1)
    side = max((x1 - x0) + 2*px, (y1 - y0) + 2*py)

    ax.set_xlim(cx - side/2, cx + side/2)
    ax.set_ylim(cy - side/2, cy + side/2)
    ax.set_aspect("equal", adjustable="box")

def auto_point_size(
    n_points: int,
    ax=None,
    *,
    # Area-based mode (when ax is given)
    target_fill: float = 0.2,   # fraction of axes area to cover with markers (0.03–0.06 is typical)
    min_size: float = 0.25,       # clamp (pt²)
    max_size: float = 6.0,        # clamp (pt²)
    # Count-only fallback (when ax is None)
    n_ref: int = 20000,           # reference point count
    s_ref: float = 1.0,           # size (pt²) to use at n_ref
    power: float = 0.5            # s ~ (n_ref / n_points)^power
) -> float:
    """
    Returns a scatter size 's' in pt². If 'ax' is provided, uses axes area so the
    total marker area ~= target_fill * axes area; otherwise uses a count-only power law.
    """
    if n_points <= 0:
        return min_size

    # --- Area-based sizing (uses figure + axes geometry; no renderer needed) ---
    if ax is not None and getattr(ax, "figure", None) is not None:
        fig = ax.figure
        # Axes bbox in figure fraction → convert to inches
        bbox = ax.get_position()                   # [0..1] figure fraction
        fig_w_in, fig_h_in = fig.get_size_inches() # inches
        ax_w_in = bbox.width  * fig_w_in
        ax_h_in = bbox.height * fig_h_in
        ax_area_in2 = ax_w_in * ax_h_in
        ax_area_pt2 = ax_area_in2 * (72.0 ** 2)    # 1 in = 72 pt

        # Share 'target_fill' across all points
        s = (target_fill * ax_area_pt2) / float(n_points)

    else:
        # --- Count-only fallback (simple power-law) ---
        s = s_ref * (float(n_ref) / float(n_points)) ** power

    # Clamp to sane bounds
    if s < min_size:
        s = min_size
    elif s > max_size:
        s = max_size
    return float(s)

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
            s = auto_point_size(len(df), ax=ax)  # dynamic point size (pt²)
            colors = df["color"] if "color" in df.columns else None
            ax.scatter(df["UMAP1"], df["UMAP2"], s=s, lw=0, alpha=0.5, color=colors)

            # Quantile-based, per-axis limits; square view so UMAP isn’t distorted
            set_square_limits(ax, df["UMAP1"].values, df["UMAP2"].values,
                              q=(0.002, 0.998),  # set to None to disable clipping
                              pad=0.03)

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
    Adds colorbars for genome_size, genome_size_log10 with:
      - legend_scale to shrink/grow the colorbars and fonts,
      - genome_min_bp/genome_max_bp thresholds (raw bp) that clamp colors:
          <= min -> genome_min_color, >= max -> genome_max_color.
    When `use_benedictus` is True, the Benedictus diverging colormap is used
    instead of the custom min/max colors but the genome size thresholds are
    still respected.
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
        "genome_size", "genome_size_log10",
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
    
    # Use user-specified number of columns, or calculate square grid
    if args.num_cols is not None:
        num_cols = args.num_cols
    else:
        num_cols = int(np.ceil(np.sqrt(total_num_cols_to_plot)))
    
    num_rows = int(np.ceil(total_num_cols_to_plot / num_cols))

    # ---------- figure & axes grid ----------
    margin = 0.25
    panel_width = 1.5
    colorbar_space = 1.5  # Extra space on the right for colorbars
    fig_width = (margin * 4) + (panel_width * num_cols) + (margin * (num_cols - 1)) + (margin * 4) + colorbar_space
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

    # Track colorbars to draw on the right side
    colorbar_info = []
    seen_colorbar_types = set()  # Track which types we've already added

    # ---------- plot the main columns ----------
    i = 0
    j = 0
    for thiscol in all_columns_to_plot:
        ax = axes[i][j]
        # Put label at the very top edge of the plot area
        ax.text(0.5, 1.0, "Clade color" if thiscol == "color" else f"{thiscol}", 
                transform=ax.transAxes, fontsize=10, ha='center', va='top')

        # Special handling for genome size panels (add colorbars, clamp colors)
        if thiscol in ("genome_size", "genome_size_log2", "genome_size_log10"):
            vals = df[thiscol].to_numpy(dtype=float)

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

            norm = Normalize(vmin=use_vmin, vmax=use_vmax, clip=True)
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

            # Set square limits to stop lines cleanly at edges with minimal padding
            set_square_limits(ax, df["UMAP1"].values, df["UMAP2"].values,
                              q=(0.002, 0.998), pad=0.001)

            # Store colorbar info for later drawing on the right side
            if thiscol == "genome_size":
                label = "Genome size"
                labels = [human_readable_bp(t) for t in np.linspace(use_vmin, use_vmax, 5)]
            elif thiscol == "genome_size_log2":
                label = "Genome size (log2)"
                labels = [human_readable_bp((2.0 ** t) - 1.0) for t in np.linspace(use_vmin, use_vmax, 5)]
            else:
                label = "Genome size (log10)"
                labels = [human_readable_bp((10.0 ** t) - 1.0) for t in np.linspace(use_vmin, use_vmax, 5)]
            
            colorbar_info.append({
                'norm': norm,
                'cmap': cmap,
                'label': label,
                'ticklabels': labels,
                'ticks': np.linspace(use_vmin, use_vmax, 5)
            })

        else:
            # Non-genome-size panels
            coltype = df[thiscol].dtype
            if thiscol == "color":
                colors = list(df[thiscol])
                ax.scatter(df["UMAP1"], df["UMAP2"], s=0.5, lw=0, alpha=0.5, color=colors)
                # Set square limits to stop lines cleanly at edges with minimal padding
                set_square_limits(ax, df["UMAP1"].values, df["UMAP2"].values,
                                  q=(0.002, 0.998), pad=0.001)
            else:
                # boolean-like in object dtype or actual boolean dtype -> map two fixed colors
                uniques = df[thiscol].unique()
                if ((True in uniques) or (False in uniques)) and (coltype == "object" or coltype == "bool"):
                    colordict = {True: "#074FF7", False: "#FD6117"}
                    colors = [colordict.get(x, "#999999") for x in df[thiscol]]
                    ax.scatter(df["UMAP1"], df["UMAP2"], s=0.5, lw=0, alpha=0.5, color=colors)
                    # Set square limits to stop lines cleanly at edges with minimal padding
                    set_square_limits(ax, df["UMAP1"].values, df["UMAP2"].values,
                                      q=(0.002, 0.998), pad=0.001)
                    
                    # Store binary colorbar info (only once)
                    if 'binary' not in seen_colorbar_types:
                        colorbar_info.append({
                            'type': 'binary',
                            'label': 'from_rbh',
                            'colors': colordict
                        })
                        seen_colorbar_types.add('binary')
                else:
                    # numeric gradient using your interpolation helper
                    arr = pd.to_numeric(df[thiscol], errors="coerce")
                    maxval = np.nanmax(arr.to_numpy(dtype=float))
                    if np.isfinite(maxval) and maxval > 1:
                        colors = [interpolate_color(x, 0, maxval, "#DCDEE3", "#FF2608") for x in arr.fillna(0)]
                        vmin, vmax = 0, maxval
                    else:
                        colors = [interpolate_color(x, 0, 1, "#DCDEE3", "#FF2608") for x in arr.fillna(0)]
                        vmin, vmax = 0, 1
                    ax.scatter(df["UMAP1"], df["UMAP2"], s=0.5, lw=0, alpha=0.5, color=colors)
                    # Set square limits to stop lines cleanly at edges with minimal padding
                    set_square_limits(ax, df["UMAP1"].values, df["UMAP2"].values,
                                      q=(0.002, 0.998), pad=0.001)
                    
                    # Store numeric gradient colorbar info (only once, for [0,1] range)
                    if 'gradient_0_1' not in seen_colorbar_types:
                        gradient_cmap = LinearSegmentedColormap.from_list("gradient", ["#DCDEE3", "#FF2608"])
                        norm = Normalize(vmin=0, vmax=1)
                        colorbar_info.append({
                            'type': 'gradient',
                            'norm': norm,
                            'cmap': gradient_cmap,
                            'label': 'Fraction [0-1]',
                            'ticks': np.linspace(0, 1, 5),
                            'ticklabels': [f"{t:.2f}" for t in np.linspace(0, 1, 5)]
                        })
                        seen_colorbar_types.add('gradient_0_1')

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
        # Put label at the very top edge of the plot area
        ax.text(0.5, 1.0, f"{thiscol}", 
                transform=ax.transAxes, fontsize=10, ha='center', va='top')

        color_col = f"{thiscol}_color"
        if color_col not in df.columns:
            raise ValueError(f"Missing expected color column {color_col} for metadata.")
        colors = df[color_col]
        ax.scatter(df["UMAP1"], df["UMAP2"], s=0.5, lw=0, alpha=0.5, color=colors)
        
        # Set square limits to stop lines cleanly at edges with minimal padding
        set_square_limits(ax, df["UMAP1"].values, df["UMAP2"].values,
                          q=(0.002, 0.998), pad=0.001)

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
    # Vertical lines (between columns) - centered between panels like in plot_paramsweep
    for jj in range(1, num_cols):
        x_pos = ((4 * margin) + (jj * panel_width) + (jj * margin) - (margin / 2)) / fig_width
        y1 = ((fig_height - (4 * margin)) / fig_height)
        y2 = ((4 * margin) / fig_height)
        line = plt.Line2D([x_pos, x_pos], [y1, y2], transform=fig.transFigure, color="#BBBBBB", lw=1.0)
        fig.add_artist(line)

    # Horizontal lines (between rows) - centered between panels like in plot_paramsweep
    # Note: x2 should stop at the plot area, not extend into colorbar space
    plot_area_right_for_lines = ((4 * margin) + (num_cols * panel_width) + ((num_cols - 1) * margin)) / fig_width
    for ii in range(1, num_rows):
        y_pos = ((fig_height - ((4 * margin) + (ii * panel_width) + (ii * margin) - (margin / 2))) / fig_height)
        x1 = ((4 * margin) / fig_width)
        x2 = plot_area_right_for_lines
        line = plt.Line2D([x1, x2], [y_pos, y_pos], transform=fig.transFigure, color="#BBBBBB", lw=1.0)
        fig.add_artist(line)

    # ---------- draw colorbars on the right side of the figure ----------
    if colorbar_info:
        # Calculate where the plot area ends
        plot_area_right = ((4 * margin) + (num_cols * panel_width) + ((num_cols - 1) * margin)) / fig_width
        
        # Position colorbars in the right margin, after the plot area
        cbar_width = 0.02  # width of each colorbar in figure coordinates
        cbar_left = plot_area_right + (2 * margin / fig_width)  # Start after plot area + some margin
        cbar_spacing = 0.03  # vertical spacing between colorbars
        
        # Calculate heights (binary colorbars are shorter)
        heights = []
        for cbar_data in colorbar_info:
            if cbar_data.get('type') == 'binary':
                heights.append(0.06)  # Short for binary
            else:
                heights.append(0.15)  # Standard height for continuous
        
        # Calculate total height needed and starting position
        total_cbar_height = sum(heights) + (len(colorbar_info) - 1) * cbar_spacing
        cbar_top = 0.5 + total_cbar_height / 2  # center vertically
        
        current_y = cbar_top
        for idx, cbar_data in enumerate(colorbar_info):
            cbar_height = heights[idx]
            cbar_bottom = current_y - cbar_height
            
            # Create a new axes for the colorbar
            cax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])
            
            # Create the colorbar based on type
            if cbar_data.get('type') == 'binary':
                # Binary colorbar (True/False)
                from matplotlib.patches import Rectangle
                cax.add_patch(Rectangle((0, 0), 1, 0.5, facecolor=cbar_data['colors'][False], edgecolor='black', linewidth=0.5))
                cax.add_patch(Rectangle((0, 0.5), 1, 0.5, facecolor=cbar_data['colors'][True], edgecolor='black', linewidth=0.5))
                cax.set_xlim(0, 1)
                cax.set_ylim(0, 1)
                cax.set_yticks([0.25, 0.75])
                cax.set_yticklabels(['False', 'True'])
                cax.set_xticks([])
                cax.tick_params(labelsize=max(5, int(6 * float(legend_scale))))
                cax.set_ylabel(cbar_data['label'], fontsize=max(6, int(7 * float(legend_scale))))
                for spine in ['top', 'right', 'bottom', 'left']:
                    cax.spines[spine].set_visible(True)
                    cax.spines[spine].set_linewidth(0.5)
            else:
                # Continuous colorbar (genome size or gradient)
                sm = ScalarMappable(norm=cbar_data['norm'], cmap=cbar_data['cmap'])
                sm.set_array([])
                cbar = plt.colorbar(sm, cax=cax)
                cbar.set_ticks(cbar_data['ticks'])
                cbar.set_ticklabels(cbar_data['ticklabels'])
                cbar.ax.tick_params(labelsize=max(5, int(6 * float(legend_scale))))
                cbar.set_label(cbar_data['label'], fontsize=max(6, int(7 * float(legend_scale))))
            
            # Move to next position
            current_y = cbar_bottom - cbar_spacing

    # ---------- save ----------
    print(f"saving the file to {outpdf}")
    plt.savefig(outpdf)
    plt.close(fig)

def plot_phyla(args, outpdf, metadata_df=None):
    """
    Plot UMAP coordinates with each phylum highlighted individually.
    Creates a multi-page PDF with:
      - First page: all species colored by their phylum
      - Subsequent pages: each phylum highlighted with its color, others greyed out
    
    Parses taxid information from 'taxid_list' or 'taxid_list_str' columns.
    Uses the phylum_to_taxid dictionary to map taxids to phyla.
    """
    from matplotlib.backends.backend_pdf import PdfPages
    import ast
    
    GREY_COLOR = "#D5D7DF"
    
    # Create reverse mapping: taxid -> phylum name
    taxid_to_phylum = {}
    for phylum_name, info in phylum_to_taxid.items():
        taxid_to_phylum[info['taxid']] = phylum_name
    
    # Load dataframe
    df = pd.read_csv(args.filelist, sep="\t", index_col=0)
    
    # Merge metadata if provided
    if metadata_df is not None:
        metadata_df = metadata_df.set_index("rbh")
        metadata_df.index = metadata_df.index.astype(str).str.strip()
        df.set_index("rbh", inplace=True, drop=False)
        df.index = df.index.astype(str).str.strip()
        matched = metadata_df.index.intersection(df.index)
        print(f"Metadata merge: matched {len(matched)} of {len(df)} UMAP RBH entries")
        df = df.join(metadata_df, how="left")
    
    # Check for required columns
    if "UMAP1" not in df.columns or "UMAP2" not in df.columns:
        raise ValueError("DataFrame must contain 'UMAP1' and 'UMAP2' columns")
    
    # Parse taxid information to determine phylum
    def parse_taxids(row):
        """Parse taxid_list or taxid_list_str to extract phylum.
        Returns the most specific (last/deepest) phylum match in the lineage."""
        taxids = []
        
        # Try taxid_list first (format: "[1, 131567, 2759, ...]")
        if "taxid_list" in row.index and pd.notna(row["taxid_list"]):
            try:
                if isinstance(row["taxid_list"], str):
                    taxids = ast.literal_eval(row["taxid_list"])
                elif isinstance(row["taxid_list"], list):
                    taxids = row["taxid_list"]
            except Exception as e:
                pass  # silently skip parse errors
        
        # Try taxid_list_str if taxid_list didn't work (format: "1;131567;2759;...")
        if not taxids and "taxid_list_str" in row.index and pd.notna(row["taxid_list_str"]):
            try:
                taxids = [int(x) for x in str(row["taxid_list_str"]).split(";") if x.strip()]
            except Exception as e:
                pass  # silently skip parse errors
        
        # Special case: Acanthocephala (taxid 10232) is now considered part of Rotifera
        # Check for this taxid first before doing normal phylum lookup
        if 10232 in taxids:
            return "Rotifera"
        
        # Find the most specific (last) phylum taxid in the list
        # This handles cases where one phylum is nested within another (e.g., Sipuncula within Annelida)
        last_phylum = None
        for taxid in taxids:
            if taxid in taxid_to_phylum:
                last_phylum = taxid_to_phylum[taxid]
        
        return last_phylum
    
    # Check if we have the required columns
    if "taxid_list" not in df.columns and "taxid_list_str" not in df.columns:
        raise ValueError("DataFrame must contain either 'taxid_list' or 'taxid_list_str' column for --plot-phyla mode")
    
    # Parse phylum for each row
    print("Parsing taxid information to determine phyla...")
    df["phylum"] = df.apply(parse_taxids, axis=1)
    
    # Get unique phyla present in the data
    phyla_in_data_original = df["phylum"].dropna().unique()
    phyla_in_data_sorted = sorted(phyla_in_data_original)
    
    # Handle custom phyla order if provided
    phyla_to_plot = phyla_in_data_sorted  # Default: plot all phyla in sorted order
    if args.phyla_order is not None:
        # Parse the space-delimited phyla names
        custom_phyla_list = args.phyla_order.split()
        
        # Validate that all requested phyla are in the data
        requested_phyla = set(custom_phyla_list)
        found_phyla_set = set(phyla_in_data_original)
        missing_requested = requested_phyla - found_phyla_set
        
        if missing_requested:
            print(f"\nWarning: Requested phyla not found in data: {sorted(missing_requested)}")
        
        # Filter to only include phyla that are both requested and found in data
        phyla_to_plot = [p for p in custom_phyla_list if p in found_phyla_set]
        
        if not phyla_to_plot:
            raise ValueError("None of the requested phyla were found in the data")
        
        print(f"\nUsing custom phyla order ({len(phyla_to_plot)} phyla): {phyla_to_plot}")
    
    # phyla_in_data is used for individual panels (respects custom order)
    phyla_in_data = phyla_to_plot
    # phyla_in_data_all is used for the "All Phyla" panel (always shows everything)
    phyla_in_data_all = list(phyla_in_data_original)
    
    # Report which phyla were found and which were not
    all_known_phyla = set(phylum_to_taxid.keys())
    found_phyla = set(phyla_in_data_original)
    missing_phyla = all_known_phyla - found_phyla
    
    print(f"\nFound {len(phyla_in_data_all)} phyla in data:")
    total_phylum_count = 0
    for phylum in sorted(phyla_in_data_all):
        count = (df["phylum"] == phylum).sum()
        print(f"  - {phylum}: {count} samples")
        total_phylum_count += count
    print(f"\nTotal samples assigned to phyla: {total_phylum_count}")
    print(f"Total samples in dataframe: {len(df)}")
    
    # If custom order was specified, also report which phyla will be plotted
    if args.phyla_order is not None:
        print(f"\nPlotting {len(phyla_in_data)} phyla in custom order:")
        for phylum in phyla_in_data:
            count = (df["phylum"] == phylum).sum()
            print(f"  - {phylum}: {count} samples")
    
    if missing_phyla:
        print(f"\nPhyla not found in data ({len(missing_phyla)}):")
        for phylum in sorted(missing_phyla):
            print(f"  - {phylum}")
    else:
        print("\nAll known phyla were found in the data!")
    
    # Report samples without phylum assignment
    unassigned_count = df["phylum"].isna().sum()
    if unassigned_count > 0:
        print(f"\nWarning: {unassigned_count} samples could not be assigned to a known phylum")
        # Get unassigned samples dataframe
        unassigned_df = df[df["phylum"].isna()]
        
        # Print detailed information for each unassigned sample
        print("\nUnassigned samples details:")
        for idx, row in unassigned_df.iterrows():
            # Get sample identifier
            sample_id = idx
            if "rbh" in df.columns:
                sample_id = row["rbh"]
            
            # Get taxid information
            taxid_info = "No taxid information"
            if "taxid_list" in row.index and pd.notna(row["taxid_list"]):
                taxid_info = f"taxid_list: {row['taxid_list']}"
            elif "taxid_list_str" in row.index and pd.notna(row["taxid_list_str"]):
                taxid_info = f"taxid_list_str: {row['taxid_list_str']}"
            
            print(f"  - {sample_id}: {taxid_info}")
    
    print(f"Found {len(phyla_in_data_all)} unique phyla in data: {sorted(phyla_in_data_all)}")
    
    # Generate colors for each phylum
    # If the dataframe has a 'color' column, use those colors directly for plotting
    # Otherwise, generate distinct colors for each phylum
    if "color" in df.columns:
        # Use the existing color column directly - each sample keeps its own color
        # This preserves the original taxonomic-specific colors
        use_existing_colors = True
        phylum_colors = {}  # We'll still need this for the legend/reference
        for phylum in phyla_in_data_all:
            phylum_mask = df["phylum"] == phylum
            if phylum_mask.any():
                # Get the most common color for this phylum (for legend purposes)
                colors = df.loc[phylum_mask, "color"].mode()
                if len(colors) > 0:
                    phylum_colors[phylum] = colors.iloc[0]
                else:
                    phylum_colors[phylum] = "#000000"
    else:
        # Generate distinct colors for each phylum
        use_existing_colors = False
        colors_list = generate_distinct_colors(len(phyla_in_data_all))
        phylum_colors = {phylum: colors_list[i] for i, phylum in enumerate(sorted(phyla_in_data_all))}
        # Assign colors to dataframe based on phylum
        df["phylum_color"] = df["phylum"].map(phylum_colors)
        df["phylum_color"] = df["phylum_color"].fillna(GREY_COLOR)
    
    # Create a color column based on phylum if it doesn't exist
    df["phylum_color"] = df["phylum"].map(phylum_colors)
    df["phylum_color"] = df["phylum_color"].fillna(GREY_COLOR)
    
    # Rotate UMAP coordinates based on the selected rotation strategy
    rotation_mode = args.phyla_rotation
    print(f"\nRotating UMAP coordinates using '{rotation_mode}' strategy...")
    umap_coords = df[["UMAP1", "UMAP2"]].values
    
    # Center the data
    mean_coords = np.mean(umap_coords, axis=0)
    centered_coords = umap_coords - mean_coords
    
    # Try multiple rotation angles to find the optimal one
    best_angle = 0
    best_metric = float('inf')
    
    for angle_deg in range(0, 180, 1):  # Try every degree from 0 to 180
        angle_rad = np.deg2rad(angle_deg)
        cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
        rotation_matrix = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
        
        rotated = centered_coords @ rotation_matrix.T
        
        # Calculate bounding box (using full range, no outlier filtering)
        x_range = rotated[:, 0].max() - rotated[:, 0].min()
        y_range = rotated[:, 1].max() - rotated[:, 1].min()
        
        if rotation_mode == "maximize_square":
            # Minimize the larger dimension (best fit for a square)
            metric = max(x_range, y_range) ** 2
        elif rotation_mode == "minimize_vertical":
            # Minimize the vertical extent
            metric = y_range
        
        if metric < best_metric:
            best_metric = metric
            best_angle = angle_deg
    
    # Apply the best rotation
    angle_rad = np.deg2rad(best_angle)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    rotation_matrix = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
    
    rotated_coords = centered_coords @ rotation_matrix.T
    
    # Now find the best reflection to minimize density in top-left quadrant
    # Test 4 orientations: original, flip-x, flip-y, flip-both
    print("Finding optimal reflection to minimize top-left density...")
    
    best_reflection = None
    best_topleft_density = float('inf')
    reflection_options = [
        ("none", np.array([[1, 0], [0, 1]])),
        ("flip_x", np.array([[-1, 0], [0, 1]])),
        ("flip_y", np.array([[1, 0], [0, -1]])),
        ("flip_both", np.array([[-1, 0], [0, -1]]))
    ]
    
    for refl_name, refl_matrix in reflection_options:
        reflected = rotated_coords @ refl_matrix.T
        
        # Get quantile-based bounds for defining quadrants
        x_median = np.median(reflected[:, 0])
        y_median = np.median(reflected[:, 1])
        y_75 = np.quantile(reflected[:, 1], 0.75)
        
        # Count points in top-left quadrant (x < median, y > 75th percentile)
        top_left_mask = (reflected[:, 0] < x_median) & (reflected[:, 1] > y_75)
        top_left_count = np.sum(top_left_mask)
        
        if top_left_count < best_topleft_density:
            best_topleft_density = top_left_count
            best_reflection = (refl_name, refl_matrix)
    
    # Apply the best reflection
    refl_name, refl_matrix = best_reflection
    rotated_coords = rotated_coords @ refl_matrix.T
    
    # Update rotation matrix to include reflection
    rotation_matrix = rotation_matrix @ refl_matrix.T
    
    # Store original axis directions (for legend)
    original_x_axis = rotation_matrix @ np.array([1, 0])
    original_y_axis = rotation_matrix @ np.array([0, 1])
    
    # Update dataframe with rotated and reflected coordinates
    df["UMAP1"] = rotated_coords[:, 0]
    df["UMAP2"] = rotated_coords[:, 1]
    
    print(f"Applied {best_angle}° rotation + '{refl_name}' reflection")
    print(f"  Top-left quadrant density: {best_topleft_density} points")
    
    # Report rotation results
    final_x_range = np.quantile(df["UMAP1"], 0.998) - np.quantile(df["UMAP1"], 0.002)
    final_y_range = np.quantile(df["UMAP2"], 0.998) - np.quantile(df["UMAP2"], 0.002)
    aspect_ratio = final_x_range / final_y_range
    print(f"Applied {best_angle}° rotation ({rotation_mode})")
    print(f"  Final aspect ratio (width/height): {aspect_ratio:.2f}")
    print(f"  X range: {final_x_range:.2f}, Y range: {final_y_range:.2f}")
    
    # Determine panel aspect ratio based on rotation mode
    if rotation_mode == "minimize_vertical":
        # Use rectangular panels that match the data aspect ratio
        use_square_panels = False
        panel_aspect_ratio = aspect_ratio
        print(f"  Using rectangular panels with aspect ratio {aspect_ratio:.2f}")
    else:
        # Use square panels (default for maximize_square)
        use_square_panels = True
        panel_aspect_ratio = 1.0
        print(f"  Using square panels")
    
    # Calculate grid layout - "All Phyla" will be 2x2, others are 1x1
    # The "All Phyla" panel takes up 4 slots (positions 0,0 to 1,1)
    num_phyla_panels = len(phyla_in_data)
    
    # Calculate grid size: we need space for num_phyla_panels regular panels 
    # plus the 2x2 "All Phyla" panel (which uses 4 grid positions but counts as 1 panel)
    # We'll add 3 to account for the 4 slots used by the large panel
    total_grid_slots = num_phyla_panels + 3
    
    # Respect --num-cols if provided, otherwise calculate from sqrt
    if args.num_cols is not None:
        num_cols = args.num_cols
    else:
        num_cols = int(np.ceil(np.sqrt(total_grid_slots)))
    
    num_rows = int(np.ceil(total_grid_slots / num_cols))
    
    # Ensure we have at least 2 rows and 2 cols for the large panel
    if num_cols < 2:
        num_cols = 2
    if num_rows < 2:
        num_rows = 2
    
    # Figure dimensions - adjust based on panel aspect ratio
    margin = 0.25
    if use_square_panels:
        panel_width = 1.5
        panel_height = 1.5
    else:
        # For minimize_vertical mode, use wider panels
        panel_height = 1.0
        panel_width = panel_height * panel_aspect_ratio
    
    fig_width = (margin * 4) + (panel_width * num_cols) + (margin * (num_cols - 1)) + (margin * 4)
    fig_height = (margin * 4) + (panel_height * num_rows) + (margin * (num_rows - 1)) + (margin * 4)
    
    # Helper function to create a phyla plot figure with configurable elements
    def create_phyla_figure(include_points=True, include_labels=True, include_grid=True, include_vectors=True):
        """
        Create a phyla plot figure with configurable elements.
        
        Args:
            include_points: Whether to plot scatter points
            include_labels: Whether to include text labels (title, phylum names)
            include_grid: Whether to include grid separator lines
            include_vectors: Whether to include UMAP axis vector legend
        
        Returns:
            matplotlib figure object
        """
        fig_new = plt.figure(figsize=(fig_width, fig_height))
        
        # Create axes grid - regular sized panels
        axes_new = [[None for _ in range(num_cols)] for _ in range(num_rows)]
        for ii in range(num_rows):
            for jj in range(num_cols):
                if ii < 2 and jj < 2:
                    continue
                left = (4 * margin) + (jj * panel_width) + (jj * margin)
                bottom = fig_height - ((4 * margin) + ((ii + 1) * panel_height) + (ii * margin))
                ax = fig_new.add_axes([
                    left / fig_width,
                    bottom / fig_height,
                    panel_width / fig_width,
                    panel_height / fig_height
                ])
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ["top", "right", "bottom", "left"]:
                    ax.spines[spine].set_visible(False)
                axes_new[ii][jj] = ax
        
        # Create the large "All Phyla" panel spanning 2x2 in top-left
        left_large = (4 * margin)
        bottom_large = fig_height - ((4 * margin) + (2 * panel_height) + (1 * margin))
        width_large = (2 * panel_width) + margin
        height_large = (2 * panel_height) + margin
        ax_large_new = fig_new.add_axes([
            left_large / fig_width,
            bottom_large / fig_height,
            width_large / fig_width,
            height_large / fig_height
        ])
        ax_large_new.set_xticks([])
        ax_large_new.set_yticks([])
        for spine in ["top", "right", "bottom", "left"]:
            ax_large_new.spines[spine].set_visible(False)
        
        # Add figure title if requested
        if include_labels:
            fig_new.suptitle(f"Phyla plot for {args.filelist}", fontsize=10)
        
        # Plot the large "All Phyla" panel
        if include_points:
            if use_existing_colors:
                # Use phyla_in_data_all to show ALL phyla, not just the custom selection
                for phylum in sorted(phyla_in_data_all):
                    phylum_mask = df["phylum"] == phylum
                    phylum_df = df[phylum_mask]
                    ax_large_new.scatter(phylum_df["UMAP1"], phylum_df["UMAP2"],
                              s=3.0, lw=0, alpha=0.6,
                              color=phylum_df["color"])
            else:
                # Use phyla_in_data_all to show ALL phyla, not just the custom selection
                for phylum in sorted(phyla_in_data_all):
                    phylum_mask = df["phylum"] == phylum
                    phylum_df = df[phylum_mask]
                    ax_large_new.scatter(phylum_df["UMAP1"], phylum_df["UMAP2"],
                              s=3.0, lw=0, alpha=0.6,
                              color=phylum_colors[phylum])
        
        if include_labels:
            ax_large_new.text(0.5, 1.0, "All Phyla", 
                    transform=ax_large_new.transAxes, fontsize=12, ha='center', va='top', fontweight='bold')
        
        # Set limits
        if use_square_panels:
            set_square_limits(ax_large_new, df["UMAP1"], df["UMAP2"], q=None, pad=0.01)
        else:
            x = df["UMAP1"].values
            y = df["UMAP2"].values
            x0, x1 = x.min(), x.max()
            y0, y1 = y.min(), y.max()
            pad_val = 0.01
            px = (x1 - x0) * pad_val
            py = (y1 - y0) * pad_val
            ax_large_new.set_xlim(x0 - px, x1 + px)
            ax_large_new.set_ylim(y0 - py, y1 + py)
            ax_large_new.set_aspect("equal", adjustable="box")
        
        # Add axis direction legend if requested
        if include_vectors:
            xlim = ax_large_new.get_xlim()
            ylim = ax_large_new.get_ylim()
            arrow_scale = 0.15 * min(xlim[1] - xlim[0], ylim[1] - ylim[0])
            arrow_origin_x = xlim[0] + 0.1 * (xlim[1] - xlim[0])
            arrow_origin_y = ylim[0] + 0.1 * (ylim[1] - ylim[0])
            
            ax_large_new.arrow(arrow_origin_x, arrow_origin_y,
                           original_x_axis[0] * arrow_scale, original_x_axis[1] * arrow_scale,
                           head_width=arrow_scale*0.15, head_length=arrow_scale*0.2,
                           fc='red', ec='red', alpha=0.7, lw=1.5)
            ax_large_new.text(arrow_origin_x + original_x_axis[0] * arrow_scale * 1.3,
                          arrow_origin_y + original_x_axis[1] * arrow_scale * 1.3,
                          'UMAP1', fontsize=7, color='red', ha='center', va='center', fontweight='bold')
            
            ax_large_new.arrow(arrow_origin_x, arrow_origin_y,
                           original_y_axis[0] * arrow_scale, original_y_axis[1] * arrow_scale,
                           head_width=arrow_scale*0.15, head_length=arrow_scale*0.2,
                           fc='blue', ec='blue', alpha=0.7, lw=1.5)
            ax_large_new.text(arrow_origin_x + original_y_axis[0] * arrow_scale * 1.3,
                          arrow_origin_y + original_y_axis[1] * arrow_scale * 1.3,
                          'UMAP2', fontsize=7, color='blue', ha='center', va='center', fontweight='bold')
        
        # Plot individual phylum panels
        phylum_idx_new = 0
        for ii in range(num_rows):
            for jj in range(num_cols):
                if ii < 2 and jj < 2:
                    continue
                if phylum_idx_new >= len(phyla_in_data):
                    break
                
                ax = axes_new[ii][jj]
                if ax is None:
                    continue
                    
                # phyla_in_data is already in the correct order (custom or sorted)
                phylum = phyla_in_data[phylum_idx_new]
                
                # Plot points if requested
                if include_points:
                    # Plot grey background
                    # Exclude both the current phylum AND Sipuncula from Annelida (since Sipuncula is nested)
                    if phylum == "Annelida":
                        # For Annelida, exclude both Annelida AND Sipuncula from background
                        other_mask = (df["phylum"] != phylum) & (df["phylum"] != "Sipuncula")
                    else:
                        other_mask = df["phylum"] != phylum
                    other_df = df[other_mask]
                    if len(other_df) > 0:
                        ax.scatter(other_df["UMAP1"], other_df["UMAP2"],
                                  s=1.0, lw=0, alpha=0.3,
                                  color=GREY_COLOR, zorder=1)
                    
                    # Plot highlighted phylum
                    # For Annelida, explicitly exclude Sipuncula since it's nested within Annelida taxonomically
                    if phylum == "Annelida":
                        phylum_mask = (df["phylum"] == phylum)
                        # The parse_taxids function already assigns Sipuncula samples to "Sipuncula", not "Annelida"
                        # So this mask will naturally exclude Sipuncula, but we're being explicit here for clarity
                    else:
                        phylum_mask = df["phylum"] == phylum
                    phylum_df = df[phylum_mask]
                    if len(phylum_df) > 0:
                        if use_existing_colors:
                            ax.scatter(phylum_df["UMAP1"], phylum_df["UMAP2"],
                                      s=1.5, lw=0, alpha=0.8,
                                      color=phylum_df["color"], zorder=2)
                        else:
                            ax.scatter(phylum_df["UMAP1"], phylum_df["UMAP2"],
                                      s=1.5, lw=0, alpha=0.8,
                                      color=phylum_colors[phylum], zorder=2)
                
                # Add label if requested
                if include_labels:
                    phylum_count = (df["phylum"] == phylum).sum()
                    ax.text(0.5, 1.0, f"{phylum}\nn = {phylum_count}", 
                            transform=ax.transAxes, fontsize=7, ha='center', va='top',
                            multialignment='center')
                
                # Set limits
                if use_square_panels:
                    set_square_limits(ax, df["UMAP1"], df["UMAP2"], q=None, pad=0.01)
                else:
                    x = df["UMAP1"].values
                    y = df["UMAP2"].values
                    x0, x1 = x.min(), x.max()
                    y0, y1 = y.min(), y.max()
                    pad_val = 0.01
                    px = (x1 - x0) * pad_val
                    py = (y1 - y0) * pad_val
                    ax.set_xlim(x0 - px, x1 + px)
                    ax.set_ylim(y0 - py, y1 + py)
                    ax.set_aspect("equal", adjustable="box")
                phylum_idx_new += 1
            
            if phylum_idx_new >= len(phyla_in_data):
                break
        
        # Add grid separator lines if requested
        if include_grid:
            for jj in range(1, num_cols):
                x_pos = ((4 * margin) + (jj * panel_width) + (jj * margin) - (margin / 2)) / fig_width
                
                if jj < 2:
                    y_top = ((fig_height - ((4 * margin) + (2 * panel_height) + (2 * margin))) / fig_height)
                    y_bottom = ((4 * margin) / fig_height)
                    line = plt.Line2D([x_pos, x_pos], [y_top, y_bottom], 
                                    transform=fig_new.transFigure, color="#BBBBBB", lw=1.0)
                    fig_new.add_artist(line)
                elif jj == 2:
                    y_top = ((fig_height - (4 * margin)) / fig_height)
                    y_bottom = ((4 * margin) / fig_height)
                    line = plt.Line2D([x_pos, x_pos], [y_top, y_bottom], 
                                    transform=fig_new.transFigure, color="#BBBBBB", lw=1.0)
                    fig_new.add_artist(line)
                else:
                    y_top = ((fig_height - (4 * margin)) / fig_height)
                    y_bottom = ((4 * margin) / fig_height)
                    line = plt.Line2D([x_pos, x_pos], [y_top, y_bottom], 
                                    transform=fig_new.transFigure, color="#BBBBBB", lw=1.0)
                    fig_new.add_artist(line)
            
            for ii in range(1, num_rows):
                y_pos = ((fig_height - ((4 * margin) + (ii * panel_height) + (ii * margin) - (margin / 2))) / fig_height)
                
                if ii < 2:
                    x_left = ((4 * margin) + (2 * panel_width) + (2 * margin)) / fig_width
                    x_right = ((fig_width - (4 * margin)) / fig_width)
                    line = plt.Line2D([x_left, x_right], [y_pos, y_pos], 
                                    transform=fig_new.transFigure, color="#BBBBBB", lw=1.0)
                    fig_new.add_artist(line)
                elif ii == 2:
                    x_left = ((4 * margin) / fig_width)
                    x_right = ((fig_width - (4 * margin)) / fig_width)
                    line = plt.Line2D([x_left, x_right], [y_pos, y_pos], 
                                    transform=fig_new.transFigure, color="#BBBBBB", lw=1.0)
                    fig_new.add_artist(line)
                else:
                    x_left = ((4 * margin) / fig_width)
                    x_right = ((fig_width - (4 * margin)) / fig_width)
                    line = plt.Line2D([x_left, x_right], [y_pos, y_pos], 
                                    transform=fig_new.transFigure, color="#BBBBBB", lw=1.0)
                    fig_new.add_artist(line)
        
        return fig_new
    
    # Create and save the main figure with all elements
    fig = create_phyla_figure(include_points=True, include_labels=True, include_grid=True, include_vectors=True)
    print(f"Saving phyla plot to {outpdf}")
    plt.savefig(outpdf)
    print(f"Saved phyla plot to {outpdf}")
    
    # If clean output requested, create additional versions
    if args.phyla_clean_output:
        # Clean version: only data points, no annotations
        clean_outpdf = outpdf.replace('.pdf', '_clean.pdf')
        print(f"\nCreating clean version without annotations...")
        fig_clean = create_phyla_figure(include_points=True, include_labels=False, include_grid=False, include_vectors=False)
        print(f"Saving clean phyla plot to {clean_outpdf}")
        plt.savefig(clean_outpdf)
        print(f"Saved clean phyla plot to {clean_outpdf}")
        plt.close(fig_clean)
        
        # Annotations-only version: all annotations, no data points
        annotations_outpdf = outpdf.replace('.pdf', '_annotations.pdf')
        print(f"\nCreating annotations-only version (no data points)...")
        fig_annot = create_phyla_figure(include_points=False, include_labels=True, include_grid=True, include_vectors=True)
        print(f"Saving annotations-only phyla plot to {annotations_outpdf}")
        plt.savefig(annotations_outpdf)
        print(f"Saved annotations-only phyla plot to {annotations_outpdf}")
        plt.close(fig_annot)
    
    plt.close(fig)

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

def infer_rank_from_subsample_filename(filepath: str) -> str:
    """
    Expect filenames like: subsample_{rank}.missing_{...}.<whatever>.{df|pdf}
    Returns the canonical rank string from return_kingdom_full_sort_order().

    Examples that match:
      subsample_phylum.neighbors_15.mind_0.1.df
      subsample_species_group.missing_large.paramsweep.pdf
    """
    fname = os.path.basename(filepath).lower()
    m = re.search(r"subsample_([^.]+)\.", fname)   # capture until first dot
    if not m:
        raise ValueError(f"Could not find 'subsample_{{rank}}.' pattern in: {filepath}")
    token = m.group(1).replace("-", "_")          # allow dashes too

    # Direct token match
    if token in _RANK_TOKEN_MAP:
        return _RANK_TOKEN_MAP[token]

    # A little forgiveness: collapse multiple underscores, strip stray suffixes
    token2 = re.sub(r"_+", "_", token).strip("_")
    if token2 in _RANK_TOKEN_MAP:
        return _RANK_TOKEN_MAP[token2]

    # Last chance: turn underscores back to spaces and check raw list
    as_words = token.replace("_", " ")
    for rk in _CANONICAL_RANKS:
        if as_words == rk.lower():
            return rk

    allowed = ", ".join(sorted(_RANK_TOKEN_MAP.keys()))
    raise ValueError(f"Unrecognized rank token '{token}' in: {filepath}\n"
                     f"Allowed tokens are: {allowed}")


def _parse_df_filename(filepath):
    """
    Parse your standard filename:
      subsample_<rank>.neighbors_<int>.mind_<float>[.missing_<...>][.method_<...>][.<metric>].df
    Returns: (num_neighbors:int, min_dist:float, method:str|None, size:str|None, metric:str|None, samplename:str)
    """
    filename = os.path.basename(filepath)

    # sample name (before .neighbors_)
    samplename = filename.split(".neighbors_")[0] if ".neighbors_" in filename else filename.split(".")[0]

    m = re.search(r"\.neighbors_(\d+)", filename)
    if not m:
        raise ValueError(f"Invalid filename (neighbors): {filename}")
    num_neighbors = int(m.group(1))

    m = re.search(r"\.mind_([0-9]*\.?[0-9]+)", filename)
    if not m:
        raise ValueError(f"Invalid filename (mind): {filename}")
    min_dist = float(m.group(1))

    # optional pieces
    method = None
    m = re.search(r"\.method_([^.]+)", filename)
    if m: method = m.group(1)

    size = None
    m = re.search(r"\.missing_([^.]+)", filename)
    if m: size = m.group(1)

    metric = None
    m = re.search(r"\.(euclidean|cosine|manhattan|chebyshev|minkowski)\.df$", filename)
    if m: metric = m.group(1)

    return num_neighbors, min_dist, method, size, metric, samplename

def load_phylo_df_by_rank_from_phylolist(files: list):
    df_by_rank = {}
    all_params = set()
    row_labels = {}

    for path in files:
        rank = infer_rank_from_subsample_filename(path)
        n, md, method, size, metric, samplename = _parse_df_filename(path)
        d = pd.read_csv(path, sep="\t", index_col=0, header=0)

        df_by_rank.setdefault(rank, {})
        df_by_rank[rank][(n, md)] = {
            "df": d, "samplename": samplename, "filepath": path,
            "num_neighbors": n, "min_dist": md, "size": size,
            "method": method, "metric": metric
        }
        row_labels[rank] = rank
        all_params.add((n, md))

    return df_by_rank, sorted(all_params), row_labels


def plot_phylo_resampling_grid(
    df_by_rank, all_params, row_labels, outpdf,
    panel=2.0,                # size of each panel (inches)
    inner=0.18,               # gap between panels (inches)
    outer=0.50,               # top/bottom and RIGHT margin (inches)
    left_gutter=1.20,         # reserved space for row labels (inches)
    sep_color="#BBBBBB",
    sep_lw=0.6
):
    import numpy as np
    import matplotlib.pyplot as plt

    RANK_ORDER = return_kingdom_full_sort_order()
    ranks_known   = [rk for rk in RANK_ORDER if rk in df_by_rank]
    ranks_unknown = sorted([rk for rk in df_by_rank.keys() if rk not in RANK_ORDER])
    ranks = ranks_known + ranks_unknown
    if not ranks:
        raise ValueError("No ranks to plot (df_by_rank is empty).")

    num_rows, num_cols = len(ranks), len(all_params)

    # Figure geometry: left gutter + right outer
    fig_w = left_gutter + (num_cols * panel + (num_cols - 1) * inner) + outer
    fig_h = (2 * outer) + (num_rows * panel + (num_rows - 1) * inner)
    fig = plt.figure(figsize=(fig_w, fig_h))

    # --- axes placement: start at LEFT_GUTTER ---
    def ax_rect(i, j):
        # i=row (0=top), j=col (0=left)
        left   = left_gutter + j * (panel + inner)     # <— was outer + ...
        bottom = outer + (num_rows - 1 - i) * (panel + inner)
        return [left / fig_w, bottom / fig_h, panel / fig_w, panel / fig_h]

    axes = [[None for _ in range(num_cols)] for _ in range(num_rows)]
    for i in range(num_rows):
        for j in range(num_cols):
            ax = fig.add_axes(ax_rect(i, j))
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ("top", "right", "bottom", "left"):
                ax.spines[sp].set_visible(False)
            axes[i][j] = ax

    # column headers
    for j, (n, md) in enumerate(all_params):
        ax = axes[0][j]
        ax.xaxis.set_label_position("top")
        ax.set_xlabel(f"n={n}, md={md}", fontsize=8)

    # row labels (now the gutter actually exists left of the panels)
    for i, rank in enumerate(ranks):
        ax = axes[i][0]
        ax.yaxis.set_label_position("left")
        ax.set_ylabel(row_labels.get(rank, rank),
                      rotation=0, ha="right", va="center",
                      fontsize=10, labelpad=10)

    # plot cells (unchanged except your auto_point_size & square limits)
    for i, rank in enumerate(ranks):
        by_param = df_by_rank.get(rank, {})
        for j, key in enumerate(all_params):
            ax = axes[i][j]
            info = by_param.get(key)
            if info is None:
                ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=8); continue
            d = info["df"]
            s = auto_point_size(len(d), ax=ax)
            if d.empty:
                ax.text(0.5, 0.5, "Empty", ha="center", va="center", fontsize=6); continue
            if "color" in d.columns:
                ax.scatter(d["UMAP1"], d["UMAP2"], s=s, lw=0, alpha=0.5, color=d["color"])
            else:
                ax.scatter(d["UMAP1"], d["UMAP2"], s=s, lw=0, alpha=0.5)

            # square, per-axis limits
            # quantile-based, per-axis limits; square view to not distort the umap
            set_square_limits(ax, d["UMAP1"].values, d["UMAP2"].values,
                  q=(0.002, 0.998),  # set to None to disable clipping
                  pad=0.03)

    # --- separator lines that span EXACTLY the panels area ---
    x_axes_left   = left_gutter
    x_axes_right  = left_gutter + num_cols * panel + (num_cols - 1) * inner
    y_axes_bottom = outer
    y_axes_top    = outer + num_rows * panel + (num_rows - 1) * inner

    xmin, xmax = x_axes_left / fig_w,  x_axes_right / fig_w
    ymin, ymax = y_axes_bottom / fig_h, y_axes_top   / fig_h

    # vertical separators (between columns)
    for j in range(num_cols - 1):
        x_mid = (left_gutter + (j + 1) * panel + j * inner + inner / 2) / fig_w
        fig.add_artist(plt.Line2D([x_mid, x_mid], [ymin, ymax],
                                  transform=fig.transFigure, color=sep_color, lw=sep_lw))
    # horizontal separators (between rows)
    for i in range(num_rows - 1):
        y_mid = (outer + (i + 1) * panel + i * inner + inner / 2) / fig_h
        fig.add_artist(plt.Line2D([xmin, xmax], [y_mid, y_mid],
                                  transform=fig.transFigure, color=sep_color, lw=sep_lw))

    print(f"saving the file to {outpdf}")
    # IMPORTANT: don't use bbox_inches='tight' or it will clip your gutter
    plt.savefig(outpdf)
    plt.close(fig)


def main():
    odpf.format_matplotlib()
    args = parse_args()

    if args.phylolist:
        df_by_rank, all_params, row_labels = load_phylo_df_by_rank_from_phylolist(args.phylolist)
        outpdf = args.prefix + ".phyloresample.pdf"
        plot_phylo_resampling_grid(df_by_rank, all_params, row_labels, outpdf)
        return

    if args.plot_phyla:
        metadatadf = parse_metadata_dfs(args.metadata) if args.metadata else None
        print(f"Metadata DataFrame:\n {metadatadf.head() if metadatadf is not None else 'No metadata provided'}")
        outpdf = args.prefix + ".phyla.pdf"
        plot_phyla(args, outpdf, metadata_df = metadatadf)
        return

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