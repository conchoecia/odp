#!/usr/bin/env python

"""
Program  : PhyloTreeUMAP.py
Language : python
Date     : 2024-02-08
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
  This program takes in a list of RBH files.
  It constructs a phylogenetic tree with those files, and then uses UMAP to visualize the tree based on the
    distance of ALG ortholog pairs from each other.

Usage instructions:
  - None yet.
"""

import argparse
from  ast import literal_eval as aliteraleval
import bokeh           # bokeh is used to visualize and save the UMAP
import networkx as nx
import numpy as np
import os
import pandas as pd
import pickle
import re
import scipy.sparse
from scipy.sparse import coo_matrix, lil_matrix, save_npz, load_npz
import sys
import time
import umap
import umap.plot
import warnings
warnings.filterwarnings("ignore", message="Graph is not fully connected", category=UserWarning)
warnings.filterwarnings("ignore", message="Hammer edge bundling is expensive for large graphs!")

# stuff for taxonomy
from ete3 import NCBITaxa,Tree

# matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from plot_ALG_fusions_v3 import assign_colors_to_nodes, SplitLossColocTree, hex_to_rgb, rgb_255_float_to_hex

# odp-specific imports
from rbh_tools import parse_rbh
from plot_ALG_fusions_v2 import taxids_to_taxidstringdict

from itertools import combinations

class PhyloTree:
    """
    This class is used to store a phylogenetic tree.
    The tree is implemented as a directional graph.
    The tree can be constructed by lists of edges. The nodes are inferred from the edges.
    """
    def __init__(self) -> None:
        # initialize the graph using networkx
        self.G = nx.DiGraph()
        self.sample_to_locdict = {}
        self.sample_to_taxidlist = {}
        self.locdf = pd.DataFrame(columns=["rbh1", "rbh2", "distance"])
        # this is a dictionary that will be used to store all the distance matrices for each sample before merging
        self.algrbhdf = None
        self.algname = None
        self.alg_combo_to_ix = None
        self.num_features = 0
        self.num_plotlevel_rows = 10

    def add_taxname_to_all_nodes(self):
        """
        This function adds the taxname to a single node. Uses ete3.
        """
        # use ete3 to get the names of the taxids
        ncbi = NCBITaxa()
        for node in self.G.nodes():
            taxid = None
            # If the node is a leaf, we need to get the taxid from the node name.
            # check if leaf if there are no descendants
            if nx.descendants(self.G, node) == set():
                # get the taxid from the node name
                taxid = int(node.split("-")[1])
            else:
                taxid = int(node)
            # for each node, make the full lineage string, in this form "Metazoa;Bilateria;Protostomes"
            lineage = ncbi.get_lineage(taxid)
            names   = ncbi.get_taxid_translator(lineage)
            self.G.nodes[node]["taxonomy_list"] = [names[taxid] for taxid in lineage]
            self.G.nodes[node]["taxid_list"]    = [taxid for taxid in lineage]
            self.G.nodes[node]["taxname"] = names[taxid]
            self.G.nodes[node]["plot_string"] = ";".join([f" {names[taxid]} ({taxid})" for taxid in lineage])
            npl = 4 # number of taxonomic units per level
            for i in range(1,self.num_plotlevel_rows+1):
                thislevel = f"level_{i}"
                j = (i-1)*npl
                self.G.nodes[node][thislevel] = ";".join([f" {names[taxid]} ({taxid})" for taxid in lineage[j:j+npl]])

        # assign the colors to the nodes
        self._assign_colors()

    def _assign_colors(self):
        """
        Assigns colors to the nodes based on some preferences.
        """
        node_colors = {}
        # go through the leaves, then assign the colors
        for thisnode in self.G.nodes():
            if nx.descendants(self.G, thisnode) == set():
                # first do the top-level colors
                for thistop in SplitLossColocTree.color_dict_top:
                    if thistop in self.G.nodes[thisnode]["taxid_list"]:
                        node_colors[thisnode] = SplitLossColocTree.color_dict_top[thistop]

        # convert the node_colors to np arrays
        node_colors = {node: np.array(hex_to_rgb(color))
                       for node, color in node_colors.items()}
        # go through the leaves, and if the color is not assigned, give it a non-offensive blue "#3f3f7f"
        # go through the leaves, then assign the colors
        for thisnode in self.G.nodes():
            if nx.descendants(self.G, thisnode) == set():
                if thisnode not in node_colors:
                    node_colors[thisnode] = np.array(hex_to_rgb("#3f3f7f"))

        # Assign colors to nodes
        root = [n for n,d in self.G.in_degree() if d==0][0]
        assign_colors_to_nodes(self.G, root, node_colors)

        # go through the graph and add a color to each node
        for node in self.G.nodes():
            if node not in node_colors:
                raise IOError(f"The node {node} does not have a color assigned to it.")
            else:
                self.G.nodes[node]["color"] = rgb_255_float_to_hex(node_colors[node])

    def ingest_ALG_rbh(self, ALGname, rbhfile) -> int:
        """
        Takes in an ALG rbh file and stores it as a dataframe.
        Safely read it in with rbh_tools.
        """
        self.algname  = ALGname
        # first check that the rbhfilepath exists
        if not os.path.exists(rbhfile):
            raise IOError(f"The file {rbhfile} does not exist.")
        self.algrbhdf = parse_rbh(rbhfile)
        # for all the values in the rbh column, get all the possible combinations of the values, and assign them an index 
        self.alg_combo_to_ix = {tuple(sorted(x)): i
                                for i, x in enumerate(list(combinations(
                                    self.algrbhdf["rbh"], 2)))}
        self.num_features = len(self.alg_combo_to_ix)
        # for all the values in the self.alg_combo_to_ix, reverse the order and add it to the dict
        for key in list(self.alg_combo_to_ix.keys()):
            self.alg_combo_to_ix[tuple(reversed(key))] = self.alg_combo_to_ix[key]
        # just return 0 in case nothing went wrong
        return 0

    def add_lineage_string_sample_distances(self, lineage_string, sample, ALGname, distdf) -> int:
        """
        The lineage strings look like this:
          - 1;131567;2759;33154;33208;6040;6042;1779146;1779162;6060;6061;1937956;6063
        The samples look like this:
          Xiphophorushellerii-8084-GCF003331165.1
        The ALGname looks like this:
          BCnSSimakov2022
        The distdf looks like this:
          rbh1                              rbh2                              distance
          Simakov2022BCnS_genefamily_12988  Simakov2022BCnS_genefamily_6122   10885675
          Simakov2022BCnS_genefamily_12988  Simakov2022BCnS_genefamily_7201   10538458
          Simakov2022BCnS_genefamily_6122   Simakov2022BCnS_genefamily_7201   347217
          Simakov2022BCnS_genefamily_12988  Simakov2022BCnS_genefamily_7465   8881006
          Simakov2022BCnS_genefamily_6122   Simakov2022BCnS_genefamily_7465   2004669
          Simakov2022BCnS_genefamily_7201   Simakov2022BCnS_genefamily_7465   1657452
          Simakov2022BCnS_genefamily_12988  Simakov2022BCnS_genefamily_9113   7935905
          Simakov2022BCnS_genefamily_6122   Simakov2022BCnS_genefamily_9113   2949770
          Simakov2022BCnS_genefamily_7201   Simakov2022BCnS_genefamily_9113   2602553

        Notes:
          - The edges from this string will be (1, 131567), (131567, 2759), ... ,
                                               (6063, Xiphophorushellerii-8084-GCF003331165.1), etc.
        """
        not_yet_seen = set()
        fields = [x for x in lineage_string.split(";")] + [sample]
        for i in range(len(fields)-1):
            self.G.add_edge(fields[i], fields[i+1])
        self.sample_to_locdict[sample] = distdf
        return 0

    def merge_sampledistances_to_locdf(self):
        """
        All this does is merge all of the samples in the sample_to_locdict to the locdf.
        Then we will modify the locdf to have the correct format, in which we can link the individual samples to the distances.
        """
        # since we're plotting everything, now is a good time to add the extra information to the nodes for plotting
        # add the lineage information to all the nodes
        self.add_taxname_to_all_nodes()

        # now we check that every node has a tstring and a taxname
        for thisnode in self.G.nodes():
            if "taxonomy_list" not in self.G.nodes[thisnode]:
                raise ValueError(f"The node {thisnode} does not have a taxonomy_string.")
            if "taxid_list" not in self.G.nodes[thisnode]:
                raise ValueError(f"The node {thisnode} does not have a taxid_list.")
            if "taxname" not in self.G.nodes[thisnode]:
                raise ValueError(f"The node {thisnode} does not have a taxname.")

        import plotly.express as px
        # make a concatdf of the dfs in the sample_to_locdict
        for key in self.sample_to_locdict:
            self.sample_to_locdict[key]["sample"] = key

        # assign each sample a row
        sample_to_row = {sample_id: i for i, sample_id in enumerate(self.sample_to_locdict.keys())}
        row_to_sample = {i: sample_id for i, sample_id in enumerate(self.sample_to_locdict.keys())}

        concatdf = pd.concat([x for x in self.sample_to_locdict.values()])
        start = time.time()
        concatdf["pair"] = concatdf.apply(lambda x: (x["rbh1"], x["rbh2"]), axis = 1)
        stop = time.time()
        print ("It took {} seconds to add the pair column with apply".format(stop - start))
        start = time.time()
        concatdf["col_indices"] = concatdf["pair"].map(self.alg_combo_to_ix)
        stop = time.time()
        print("It took {} seconds to add the col_indices column with map".format(stop - start))
        concatdf["row_indices"] = concatdf["sample"].map(sample_to_row)

        # DIMENSIONALITY REDUCTION
        # count the number of times each pair occurs
        pair_counts = concatdf["pair"].value_counts()
        # get the most frequent 1000 combinations
        pair_counts = pair_counts.nlargest(10000)
        print("The largest 1000 pairs are: ", pair_counts)
        pair_to_ix = {pair: i for i, pair in enumerate(pair_counts.index)}
        # filter the concatdf to only contain the most frequent 1000 pairs
        concatdf = concatdf[concatdf["pair"].isin(pair_counts.index)]
        # reset the values of the col_indices column of concatdf to match the row index in pair_counts
        concatdf["col_indices"] = concatdf["pair"].map(pair_to_ix)
        print(concatdf)

        values       = np.array(concatdf["distance"]   )
        row_indices  = np.array(concatdf["row_indices"])
        col_indices  = np.array(concatdf["col_indices"])
        num_features = len(pair_counts)
        num_samples  = len(sample_to_row)

        # construct the COO matrix
        #sparse_matrix = coo_matrix((concatdf["distance"],
        #                            (concatdf["row_indices"], concatdf["col_indices"])),
        #                            shape = (len(sample_to_row), self.num_features))
        sparse_matrix = coo_matrix(( values,
                                    (row_indices, col_indices)),
                                    shape = (num_samples, num_features))
        del concatdf

        sparse_matrix = sparse_matrix.tolil()
        # set the missing values of the sparse matrix to 999999999999
        sparse_matrix.data[sparse_matrix.data == 0] = 999999999999

        print("Fitting the UMAP")
        reducer = umap.UMAP(low_memory=True)
        start = time.time()
        mapper = reducer.fit(sparse_matrix)
        stop = time.time()
        print("It took {} seconds to fit_transform the UMAP".format(stop - start))
        del sparse_matrix
        #              ┓    •
        # ┓┏┏┳┓┏┓┏┓  ┏┓┃┏┓╋╋┓┏┓┏┓
        # ┗┻┛┗┗┗┻┣┛  ┣┛┗┗┛┗┗┗┛┗┗┫
        #        ┛   ┛          ┛
        color_dict = {i: self.G.nodes[row_to_sample[i]]["color"]
                      for i in sorted(row_to_sample.keys())}
        hover_data = pd.DataFrame({
                                   "label":   [row_to_sample[i]                          for i in sorted(row_to_sample.keys())],
                                   "taxname": [self.G.nodes[row_to_sample[i]]["taxname"] for i in sorted(row_to_sample.keys())],
                                   "color":   [self.G.nodes[row_to_sample[i]]["color"]   for i in sorted(row_to_sample.keys())]
                                   })
        for i in range(1,self.num_plotlevel_rows+1):
            thislevel = f"level_{i}"
            hover_data[thislevel] = [self.G.nodes[row_to_sample[i]][thislevel] for i in sorted(row_to_sample.keys())]

        print(hover_data)
        plot = umap.plot.interactive(mapper,
                                     color_key = color_dict,
                                     labels = [row_to_sample[i] for i in sorted(row_to_sample.keys())],
                                     hover_data = hover_data,
                                     point_size = 4
                                     )
        # output to an HTML file
        bokeh.io.output_file("distances_UMAP_sparse_bokeh.html")
        # Save the plot to an HTML file
        bokeh.io.save(plot)

        # ┏┓┓ ┏┓┏┳┓┓ ┓┏
        # ┃┃┃ ┃┃ ┃ ┃ ┗┫
        # ┣┛┗┛┗┛ ┻ ┗┛┗┛
        # get the coordinates of the UMAP
        df_embedding = pd.DataFrame(mapper.embedding_, columns=['UMAP1', 'UMAP2'])
        # Add the indices as labels. Use the keys of the sample_to_locdict, sorted by the values
        df_embedding['label'] = df_embedding.index.map(row_to_sample)
        df_embedding["color"] = df_embedding.index.map(color_dict)
        print(df_embedding)
        # Add colors to the plot
        # Assuming you have a 'color' column in your DataFrame indicating the color of each point
        #fig = px.scatter(df_embedding,
        #                 x='UMAP1', y='UMAP2',
        #                 color='color',
        #                 hover_name='label')
        fig = px.scatter()
        for color, data in df_embedding.groupby('color'):
            fig.add_scatter(
                x=data['UMAP1'],
                y=data['UMAP2'],
                mode='markers',
                marker=dict(color=color),
                text=data['label'],
                name=color  # Optional: Set the legend name
            )
        # Show the plot
        outhtml = "distances_UMAP_sparse_plotly.html"
        fig.write_html(outhtml)
        # clear the figure
        plt.clf()
        #         ┓   ┓•┓
        # ┏┳┓┏┓╋┏┓┃┏┓╋┃┓┣┓
        # ┛┗┗┗┻┗┣┛┗┗┛┗┗┗┗┛
        #       ┛
        # make a matplotlib plot of the UMAP with the df_embedding, and the color_dict from SplitLossColocTree as the legend
        # make a figure that is 5x5 inches
        fig = plt.figure(figsize=(5,5))
        # scatter the UMAP1 and UMAP2 columns of the df_embedding
        fig = plt.scatter(df_embedding["UMAP1"], df_embedding["UMAP2"], c = df_embedding["color"])
        # make a legend with the color_dict from SplitLossColocTree
        nbci = NCBITaxa()
        # get the name of the ncbi taxid from the SplitLossColocTree color_dict
        legend_dict = {}
        for key in SplitLossColocTree.color_dict_top:
            taxid = int(key)
            taxname = nbci.get_taxid_translator([taxid])[taxid]
            legend_dict[taxname] = SplitLossColocTree.color_dict_top[key]
        print("This is the legend dict")
        print(legend_dict)
        legend_patches = [mpatches.Patch(color=color, label=label)
                          for label, color in legend_dict.items()]
        # add the entries to the legend
        fig = plt.legend(handles=legend_patches, loc="upper right")
        # save the figure
        plt.savefig("distances_UMAP_sparse_matplotlib.pdf")
        sys.exit()

def umap_mapper_to_QC_plots(mapper, outfile, title = "UMAP Connectivity"):
    """
    This makes all the QC plots that the UMAP program can make.
    """
    umap.plot.diagnostic(mapper, diagnostic_type='pca')
    # add the title to the plot
    ax.set_title(title)
    # save the plot to a file
    # if the output type is a raster, change the output resolution to 600 dpi
    raster_formats = [".png", ".jpg", ".jpeg", ".tiff", ".tif"]
    if any([outfile.endswith(x) for x in raster_formats]):
        plt.savefig(outfile, dpi = 900)
    else:
        plt.savefig(outfile)
    plt.savefig(outfile)

def umap_mapper_to_connectivity(mapper, outfile, bundled = False, title = "UMAP Connectivity"):
    """
    This makes connectivity plots of the UMAP to visualize the data.
    """
    if bundled:
        # ignore the UserWarning: Hammer edge bundling is expensive for large graphs!
        ax = umap.plot.connectivity(mapper, show_points = True, edge_bundling='hammer')
    else:
        ax = umap.plot.connectivity(mapper, show_points = True)
    # add the title to the plot
    ax.set_title(title)
    # save the plot to a file
    # if the output type is a raster, change the output resolution to 600 dpi
    raster_formats = [".png", ".jpg", ".jpeg", ".tiff", ".tif"]
    if any([outfile.endswith(x) for x in raster_formats]):
        plt.savefig(outfile, dpi = 900)
    else:
        plt.savefig(outfile)
    plt.savefig(outfile)

def umap_mapper_to_df(mapper, cdf):
    """
    This function takes a UMAP mapper and a dataframe of the distances and returns a dataframe with the UMAP coordinates.
    """
    # get the coordinates of the UMAP
    df_embedding = pd.DataFrame(mapper.embedding_, columns=['UMAP1', 'UMAP2'])
    # This is all for debugging
    #print(dir(mapper))
    #print("This is densmap")
    #print(mapper.densmap)
    #print("This is the embedding")
    #print(mapper.embedding_)
    #print("This is the fit")
    #print(mapper.fit)
    #print("This is the fit_transform")
    #print(mapper.fit_transform)
    #print("This is knn_dists")
    #print(mapper.knn_dists)
    #print("THis is local_connectivity")
    #print(mapper.local_connectivity)
    #print("This is output metric")
    #print(mapper.output_metric)
    #print("This is precomputed_knn")
    #print(mapper.precomputed_knn)
    #print("This is random_state")
    #print(mapper.random_state)
    #print("This is repulsion_strength")
    #print(mapper.repulsion_strength)
    #print("This is tqdm_kwds")
    #print(mapper.tqdm_kwds)
    # add those two columns to the cdf
    return pd.concat([cdf, df_embedding], axis = 1)

def rbh_to_gb(sample, rbhdf, outfile):
    """
    Converts the rbh dataframe to a groupby object. This is saved to a file for later consumption.
    """
    # Initialize an empty DataFrame to store the results
    result_df = pd.DataFrame(columns=["rbh1", "rbh2", "distance"])

    # Group by scaffold and iterate over each group
    gb = rbhdf.groupby(f"{sample}_scaf")
    for name, group in gb:
        # Get all combinations of rbh pairs within the group
        combos = pd.DataFrame(list(combinations(group["rbh"], 2)), columns=["rbh1", "rbh2"])

        # Merge combinations with group to get position information
        merged = pd.merge(combos, group[[f"{sample}_pos", "rbh"]], left_on="rbh1", right_on="rbh")
        merged = pd.merge(merged, group[[f"{sample}_pos", "rbh"]], left_on="rbh2", right_on="rbh")

        # Calculate absolute distance and add to result DataFrame
        merged["distance"] = abs(merged[f"{sample}_pos_x"] - merged[f"{sample}_pos_y"])
        result_df = pd.concat([result_df, merged[["rbh1", "rbh2", "distance"]]], ignore_index=True)

    # Swap the values of the rbh1 and rbh2 columns if they are not in alphabetical order.
    result_df["rbh1"], result_df["rbh2"] = np.where(result_df["rbh1"] < result_df["rbh2"],
                                                    (result_df["rbh1"], result_df["rbh2"]),
                                                    (result_df["rbh2"], result_df["rbh1"]))
    # verify that all of the values in the rbh1 column are lexicographically less than the values in the rbh2 column
    if not all(result_df["rbh1"] < result_df["rbh2"]):
        raise IOError("The values in the rbh1 column are not lexicographically less than the values in the rbh2 column. These need to be sorted.")
    # DO NOT sort by the rbh1 and rbh2 columns. This does not help compression

    # Save the result DataFrame to a tsv file
    result_df.to_csv(outfile, sep="\t", index=False, compression="gzip")

def NCBI_taxid_to_taxdict(ncbi, taxid):
    """
    Takes a single NCBI taxid as input and returns a dictionary with useful information:

    Input:
      - ncbi:  The NCBITaxa object
      - taxid: The NCBI taxid
    Output:
      - A dictionary with the following
        taxid: The taxid, same as the input
        taxname: The name of this specific taxid
        taxname_list: A list of the taxonomy names, like ["root", "cellular organisms", "Eukaryota", "Opisthokonta"]
        taxid_list: A list of the taxids, like [1, 131567, 2759, 33154]
        level_1: The first level of the taxid, like "root (1); cellular organisms (131567); Eukaryota (2759); Opisthokonta (33154)"
        ... up to level_10
        printstring: The printstring of the taxid, like "root (1); cellular organisms (131567); Eukaryota (2759); Opisthokonta (33154)"
    """
    if isinstance(taxid, str):
        # first check that the taxid is an integer
        if not re.match(r"^[0-9]*$", taxid):
            raise ValueError(f"There is a non-numeric character in the taxid string, {taxid}, for file {thisfile}. Exiting.")
    elif isinstance(taxid, int):
        pass
    else:
        raise ValueError(f"The taxid is not a string or an integer. It is a {type(taxid)}. Exiting.")

    # safe, get the lineage
    entry = {"taxid": taxid}
    # for each node, make the full lineage string, in this form "Metazoa;Bilateria;Protostomes"
    lineage = ncbi.get_lineage(taxid)
    names   = ncbi.get_taxid_translator(lineage)
    # make sure that the lineage and the names are not empty
    if len(lineage) == 0:
        raise ValueError(f"The lineage is empty for the taxid {taxid}. Exiting.")
    if len(names) == 0:
        raise ValueError(f"The names are empty for the taxid {taxid}. Exiting.")
    entry["taxname"]          = names[taxid]
    entry["taxid_list"]       = [taxid for taxid in lineage]
    entry["taxid_list_str"]   = ";".join([str(taxid) for taxid in lineage])
    entry["taxname_list"]     = [names[taxid] for taxid in lineage]
    entry["taxname_list_str"] = ";".join([names[taxid] for taxid in lineage])

    npl = 4 # number of taxonomic units per level
    num_rows = 10 # what level do we want to go to
    for i in range(1, num_rows+1):
        thislevel = f"level_{i}"
        j = (i-1)*npl
        entry[thislevel] = ";".join([f" {names[taxid]} ({taxid})" for taxid in lineage[j:j+npl]])

    entry["printstring"]  = ";".join([f" {names[taxid]} ({taxid})" for taxid in lineage])
    return entry

def create_directories_if_not_exist(file_path):
    """
    This takes an iterable of file paths or directories for which we want to safely create the directories.
    """
    print("requested file path: ", file_path)
    # Split the path into directories
    directories = file_path.split(os.sep)

    # Iterate over each directory and create if it doesn't exist
    path_so_far = ''
    for directory in directories:
        path_so_far = os.path.join(path_so_far, directory)
        if not os.path.exists(path_so_far):
            os.makedirs(path_so_far)

def rbh_to_distance_gbgz(rbhfile, outfile, ALGname):
    """
    This takes a single rbh file and converts it to a distance matrix.
    It compresses the distance matrix and saves it to a file.
    This program does two things.
      1. It takes all of the rbh files in the directory and calculates the distance matrix of all of the BCnS ALGs.
        - The distance matrix is saved in a file called GTUMAP/distance_matrices/{sample}.gb.gz
        - The columns of the distance matrix are: rbh1, rbh2, distance. See the example below
        ```
        rbh1                              rbh2                              distance
        Simakov2022BCnS_genefamily_10008  Simakov2022BCnS_genefamily_1023   98045821
        Simakov2022BCnS_genefamily_10008  Simakov2022BCnS_genefamily_10663  3425056
        Simakov2022BCnS_genefamily_1023   Simakov2022BCnS_genefamily_10663  101470877
        Simakov2022BCnS_genefamily_10008  Simakov2022BCnS_genefamily_10751  86114004
        Simakov2022BCnS_genefamily_1023   Simakov2022BCnS_genefamily_10751  11931817
        Simakov2022BCnS_genefamily_10663  Simakov2022BCnS_genefamily_10751  89539060
        Simakov2022BCnS_genefamily_10008  Simakov2022BCnS_genefamily_10927  42360830
        ```
    """
    # We must check that all of the rbh files, when split on '-', have an integer as the 2nd element.
    # If not, the filename needs to be changed. Right now we parse the taxid from the filename.
    thisfilename = os.path.basename(rbhfile)
    taxid = thisfilename.split('-')[1]
    if not re.match(r"^[0-9]*$", str(taxid)):
        raise ValueError(f"There is a non-numeric character in the taxid string for file {rbhfile}. Exiting.")

    df = parse_rbh(rbhfile)
    # make sure that ALG "_scaf", "_gene", and "_pos" are in the columns.
    for col in ["_scaf", "_gene", "_pos"]:
        thiscol = f"{ALGname}{col}"
        if thiscol not in df.columns:
            raise IOError(f"The column {thiscol} is not in the rbh file {rbhfile}. Exiting.")
    thissample = [x for x in df.columns
                  if "_scaf" in x
                  and ALGname not in x][0].split("_")[0]
    # check that the second field when splitting on '-' is an integer
    if not re.match( r"^[0-9]*$", thissample.split("-")[1] ):
        raise ValueError( f"There is a non-numeric character in the taxid string for the sample {thissample} when split with '-'. The file was {rbhfile} Exiting." )
    gb_filepath = outfile
    # check that the file ends in .gb.gz
    if not gb_filepath.endswith(".gb.gz"):
        raise IOError(f"The file {gb_filepath} does not end in .gb.gz. Exiting.")
    if not os.path.exists(gb_filepath):
        rbh_to_gb(thissample, df, gb_filepath)

def sampleToRbhFileDict_to_sample_matrix(sampleToRbhFileDict, ALGname,
                                         gbgz_directory,
                                         outtsv, unannotated_color = "#3f3f7f"):
    """
    This is similar to the rbh_directory_to_distance_matrix,
    but it does not calculate the distance matrix.
    """

    # We must check that all of the rbh files, when split on '-', have an integer as the 2nd element.
    # If not, the filename needs to be changed. Right now we parse the taxid from the filename.
    for thissample in sampleToRbhFileDict:
        taxid = thissample.split('-')[1]
        if not re.match(r"^[0-9]*$", str(taxid)):
            raise ValueError(f"There is a non-numeric character in the taxid string for file {rbhfile}. Exiting.")

    # the entries dict will contain the sample information before concatenating to a new df.
    entries = []
    ncbi = NCBITaxa() # set this up, as we will use this tool once for each sample
    # print the rbh files
    i = 1
    for key in sampleToRbhFileDict:
        print("\r   Parsing the rbh file: {}/{}   ".format(i+1, len(sampleToRbhFileDict)), end="", file = sys.stdout)
        rbhfile = sampleToRbhFileDict[key]
        # get the taxid from the filename
        thisfilename = os.path.basename(rbhfile)
        # get the taxid. When we split on '-', it will be the 1st element, zero-based indexing.
        taxid = key.split('-')[1]
        # check that the taxid is an integer
        if not re.match(r"^[0-9]*$", str(taxid)):
            raise ValueError("There is a non-numeric character in the taxid string")
        taxid = int(taxid)

        df = parse_rbh(rbhfile)
        # make sure that ALG "_scaf", "_gene", and "_pos" are in the columns.
        for col in ["_scaf", "_gene", "_pos"]:
            thiscol = f"{ALGname}{col}"
            if thiscol not in df.columns:
                raise IOError(f"The column {thiscol} is not in the rbh file {rbhfile}. Exiting.")
        thissample = [x for x in df.columns
                      if "_scaf" in x
                      and ALGname not in x][0].split("_")[0]
        if thissample != key:
            raise ValueError(f"The sample {thissample} is not the same as the key {key}. Exiting.")
        # check that the second field when splitting on '-' is an integer
        if not re.match( r"^[0-9]*$", thissample.split("-")[1] ):
            raise ValueError( f"There is a non-numeric character in the taxid string for the sample {thissample} when split with '-'. The file was {rbhfile} Exiting." )
        # This is where we skip the gb.gz distance matrix creation.

        gb_filepath = os.path.join(gbgz_directory, f"{thissample}.gb.gz")
        # now we add the remaining necessary information to the entries dict
        taxid_dict = {"sample": thissample}
        # add all the outputs of NCBITaxa to the taxid_dict
        taxid_dict.update(NCBI_taxid_to_taxdict(ncbi, taxid))
        # now we need to add the genome size, number of chromosomes, and the filename
        # The genome size is the maximum value of each of the summed {thissample}_pos columns when grouped by {thissample}_scaf
        taxid_dict["genome_size"] = df.groupby(f"{thissample}_scaf").max().sum()[f"{thissample}_pos"]
        # The number of chromosomes is the number of unique {thissample}_scaf values
        taxid_dict["number_of_chromosomes"] = df[f"{thissample}_scaf"].nunique()
        taxid_dict["rbh_filepath"] = rbhfile
        taxid_dict["rbh_filepath_abs"] = os.path.abspath(rbhfile)
        taxid_dict["rbh_filename"] = os.path.basename(rbhfile)
        taxid_dict["dis_filepath"] = gb_filepath
        taxid_dict["dis_filepath_abs"] = os.path.abspath(gb_filepath)
        taxid_dict["dis_filename"] = os.path.basename(gb_filepath)
        taxid_dict["color"] = unannotated_color
        # now see if we should update the color further
        for thistaxid in taxid_dict["taxid_list"][::-1]:
            if int(thistaxid) in SplitLossColocTree.color_dict_top:
                taxid_dict["color"] = SplitLossColocTree.color_dict_top[thistaxid]
                break
        entries.append(taxid_dict)
        i = i+1

    # make a dataframe from the entries dict. Save it as a tsv to the outtsv file
    sampledf = pd.DataFrame(entries)
    # sort by the taxid_list_str column, reset index
    sampledf = sampledf.sort_values(by = "taxid_list_str").reset_index(drop = True)
    # move the color column to right after the sample column
    sampledf = sampledf[["sample", "color"] + [col for col in sampledf.columns if col not in ["sample", "color"]]]
    sampledf.to_csv(outtsv, sep = "\t", index = True)
    return sampledf

def rbh_directory_to_distance_matrix(rbh_directory, ALGname, unannotated_color = "#3f3f7f",
                                     outtsv = "GTUMAP/sampledf.tsv",
                                     outputdir = "GTUMAP/distance_matrices/") -> pd.DataFrame:
    """
    This program does two things.
      1. It takes all of the rbh files in the directory and calculates the distance matrix of all of the BCnS ALGs.
        - The distance matrix is saved in a file called GTUMAP/distance_matrices/{sample}.gb.gz
        - The columns of the distance matrix are: rbh1, rbh2, distance. See the example below
        ```
        rbh1                              rbh2                              distance
        Simakov2022BCnS_genefamily_10008  Simakov2022BCnS_genefamily_1023   98045821
        Simakov2022BCnS_genefamily_10008  Simakov2022BCnS_genefamily_10663  3425056
        Simakov2022BCnS_genefamily_1023   Simakov2022BCnS_genefamily_10663  101470877
        Simakov2022BCnS_genefamily_10008  Simakov2022BCnS_genefamily_10751  86114004
        Simakov2022BCnS_genefamily_1023   Simakov2022BCnS_genefamily_10751  11931817
        Simakov2022BCnS_genefamily_10663  Simakov2022BCnS_genefamily_10751  89539060
        Simakov2022BCnS_genefamily_10008  Simakov2022BCnS_genefamily_10927  42360830
        ```
      2. It collates information from the rbh files into a sampledf.
        - This is the dataframe with the following columns:
          - index: The index of the dataframe is important because this will be the order of the samples in the distance matrix. Everything will be ordered by this index.
          - sample: This is the sample name. This is the same sample information that will be in the rbh file columns, and in the distance matrix.
          - taxid: This is the NCBI TAXid of the sample.
          - taxname: This is the name of the taxid. For example 7777 is "Condrichthyes".
          - taxid_list: This is a list of all the taxids in the lineage of the sample from closest to root to furthest.
          - taxid_list_str: A string version of taxid_list joined together with ';' characters
          - taxname_list: This is a list of all the taxnames in the lineage of the sample from closest to root to furthest. Matches the indices of taxid_list.
          - taxname_list_str: A string version of taxname_list joined together with ';' characters
          - level_1: Here, the level_1, level_2, et cetera are splits of the NCBI taxid for easy plotting. These go up to level_10.
              - Some examples are:
                ```
                level_1: root (1); cellular organisms (131567); Eukaryota (2759); Opisthokonta (33154)
                level_2: Metazoa (2); Eumetazoa (6072); Bilateria (33213); Protostomia (33317)
                et cetera
                ```
          - printstring: This is the printstring of the taxid, like "root (1); cellular organisms (131567); Eukaryota (2759); Opisthokonta (33154)"
          - #NOT YET annotation_method: This shows what method was used to annotate the sample, for example "BCnSSimakov2022_protmap" or "NCBI/Existing"
          - genome_size: This is the size of the genome in base pairs.
          - #NOT YET gc_content: This is the GC content of the genome. This is expressed as a value from 0 to 1.
          - number_of_chromosomes: The haploid number of chromosomes that this species has.
          - rbh_filepath: The filepath, as provided, of the rbh file that was used to generate this information.
          - rbh_filepath_abs: The filepath, resolved, of the rbh file that was used to generate this information.
          - rbh_filename: The filename of the rbh file that was used to generate this information.
          - dis_filepath: The filepath, as provided, of the distance .gb.gz file that was generated.
          - dis_filepath_abs: The filepath, resolved, of the distance .gb.gz file that was generated.
          - dis_filename: The filename of the gb.gz file that was generated.

      The input:
       - rbh_directory: The directory that contains the RBH files. These RBH files are from odp.
       - ALGname: The name of the ALG that we are looking at. This is required to help us parse the columns of the rbh files.
       - outputdir

      The output:
        - One file per sample in this format:
          `GTUMAP/distance_matrices/{sample}.gb.gz`
        - One file called `GTUMAP/sampledf.tsv` that contains the sampledf.
        - Returns the sampledf as a pandas dataframe.

    Saves all of them in a directory called results. They
    """
    # safely create the required output directories.
    create_directories_if_not_exist(outtsv)
    create_directories_if_not_exist(outputdir)

    # get the rbh files in the directory
    rbh_files = list(sorted([os.path.join(rbh_directory, f)
                 for f in os.listdir(rbh_directory)
                 if f.endswith('.rbh')], reverse = True))

    # We must check that all of the rbh files, when split on '-', have an integer as the 2nd element.
    # If not, the filename needs to be changed. Right now we parse the taxid from the filename.
    for rbhfile in rbh_files:
        thisfilename = os.path.basename(rbhfile)
        taxid = thisfilename.split('-')[1]
        if not re.match(r"^[0-9]*$", str(taxid)):
            raise ValueError(f"There is a non-numeric character in the taxid string for file {rbhfile}. Exiting.")

    # the entries dict will contain the sample information before concatenating to a new df.
    entries = []
    ncbi = NCBITaxa() # set this up, as we will use this tool once for each sample
    # print the rbh files
    for i in range(len(rbh_files)):
        print("\r   Parsing the rbh file: {}/{}   ".format(i+1, len(rbh_files)), end="", file = sys.stdout)
        rbhfile = rbh_files[i]
        # get the taxid from the filename
        thisfilename = os.path.basename(rbhfile)
        # get the taxid. When we split on '-', it will be the 1st element, zero-based indexing.
        taxid = thisfilename.split('-')[1]
        # check that the taxid is an integer
        if not re.match(r"^[0-9]*$", str(taxid)):
            raise ValueError("There is a non-numeric character in the taxid string")
        taxid = int(taxid)

        df = parse_rbh(rbhfile)
        # make sure that ALG "_scaf", "_gene", and "_pos" are in the columns.
        for col in ["_scaf", "_gene", "_pos"]:
            thiscol = f"{ALGname}{col}"
            if thiscol not in df.columns:
                raise IOError(f"The column {thiscol} is not in the rbh file {rbhfile}. Exiting.")
        thissample = [x for x in df.columns
                      if "_scaf" in x
                      and ALGname not in x][0].split("_")[0]
        # check that the second field when splitting on '-' is an integer
        if not re.match( r"^[0-9]*$", thissample.split("-")[1] ):
            raise ValueError( f"There is a non-numeric character in the taxid string for the sample {thissample} when split with '-'. The file was {rbhfile} Exiting." )
        gb_filepath = os.path.join(outputdir, f"{thissample}.gb.gz")
        if not os.path.exists(gb_filepath):
            rbh_to_gb(thissample, df, gb_filepath)

        # now we add the remaining necessary information to the entries dict
        taxid_dict = {"sample": thissample}
        # add all the outputs of NCBITaxa to the taxid_dict
        taxid_dict.update(NCBI_taxid_to_taxdict(ncbi, taxid))
        # now we need to add the genome size, number of chromosomes, and the filename
        # The genome size is the maximum value of each of the summed {thissample}_pos columns when grouped by {thissample}_scaf
        taxid_dict["genome_size"] = df.groupby(f"{thissample}_scaf").max().sum()[f"{thissample}_pos"]
        # The number of chromosomes is the number of unique {thissample}_scaf values
        taxid_dict["number_of_chromosomes"] = df[f"{thissample}_scaf"].nunique()
        taxid_dict["rbh_filepath"] = rbhfile
        taxid_dict["rbh_filepath_abs"] = os.path.abspath(rbhfile)
        taxid_dict["rbh_filename"] = os.path.basename(rbhfile)
        taxid_dict["dis_filepath"] = gb_filepath
        taxid_dict["dis_filepath_abs"] = os.path.abspath(gb_filepath)
        taxid_dict["dis_filename"] = os.path.basename(gb_filepath)
        taxid_dict["color"] = unannotated_color
        # now see if we should update the color further
        for thistaxid in taxid_dict["taxid_list"][::-1]:
            if int(thistaxid) in SplitLossColocTree.color_dict_top:
                taxid_dict["color"] = SplitLossColocTree.color_dict_top[thistaxid]
                break
        entries.append(taxid_dict)

    # Now we add the color information to each row.
    # At this point, we are not factoring in information about internal nodes,
    # so we can add the colors given the annotations

    # make a dataframe from the entries dict. Save it as a tsv to the outtsv file
    sampledf = pd.DataFrame(entries)
    # sort by the taxid_list_str column, reset index
    sampledf = sampledf.sort_values(by = "taxid_list_str").reset_index(drop = True)
    # move the color column to right after the sample column
    sampledf = sampledf[["sample", "color"] + [col for col in sampledf.columns if col not in ["sample", "color"]]]
    sampledf.to_csv(outtsv, sep = "\t", index = True)
    print()
    print("Done parsing the rbh files")
    return sampledf

def parse_args():
    """
    This has all of the args that we need to parse.
    The args that we need to parse are:
      -d --directory : The directory that contains the RBH files. These RBH files
      -V --overwrite : If this is set, then we will overwrite the files in the output directory.
                       Otherwise, try to load the existing files. This is not required.
      -a --ALGname   : The name of the ALG that we are looking at. This is required.
      -r --rbhfile   : The ALG rbh file. This is required, as we will use this information later.
      -c --cladestoplot : This is a list of the clades that we want to plot. This is not required. Use comma-separated values.
      -q --qualitycontrol : If this is set, then we will run the quality control on the distance matrices. This is not required.
    """
    parser = argparse.ArgumentParser(description='This program takes in a list of RBH files. It constructs a phylogenetic tree with those files, and then uses UMAP to visualize the tree based on the distance of ALG ortholog pairs from each other.')
    parser.add_argument('-d', '--directory',
                        required = True,
                        type=str,
                        help='The directory that contains the RBH files. These RBH files')
    parser.add_argument('-V', '--overwrite', action='store_true',
                        help='If this is set, then we will overwrite the files in the output directory. Otherwise, try to load the existing files. This is not required.')
    parser.add_argument('-a', '--ALGname', type=str,
                        required = True,
                        help='The name of the ALG that we are looking at. This is required.')
    parser.add_argument('-r', '--rbhfile', type=str,
                        required = True,
                        help='The ALG rbh file. This is required, as we will use this information later.')
    parser.add_argument('-c', '--cladestoplot', type=str, default = "33208,6605", # the default is Metazoa and Cephalopoda
                        help='This is a list of the clades that we want to plot. This is not required. Use comma-separated values.')
    parser.add_argument('-q', '--qualitycontrol', action='store_true',
                        help='If this is set, then we will generate quality control plots of the data. This is not required.')
    args = parser.parse_args()

    # check that the directory exists
    if not os.path.exists(args.directory):
        raise IOError(f"The directory {args.directory} does not exist. Exiting.")
    # go through clades to plot, split on commas, check that everything can be parsed as an int, and reassign to args
    fields = args.cladestoplot.split(",")
    for thisfield in fields:
        if not re.match(r"^[0-9]*$", thisfield):
            raise ValueError(f"The field {thisfield} is not a number. Exiting.")
    args.cladestoplot = [int(x) for x in fields]

    return args

def umap_mapper_to_bokeh(mapper, sampledf, outhtml, plot_title = "UMAP"):
    """
    This takes a UMAP mapper and a sampledf and returns a bokeh plot.
    """
    if not outhtml.endswith(".html"):
        raise ValueError(f"The output file {outhtml} does not end with '.html'. Exiting.")

    #              ┓    •
    # ┓┏┏┳┓┏┓┏┓  ┏┓┃┏┓╋╋┓┏┓┏┓
    # ┗┻┛┗┗┗┻┣┛  ┣┛┗┗┛┗┗┗┛┗┗┫
    #        ┛   ┛          ┛
    hover_data = pd.DataFrame({
                               "label":   sampledf["sample"],
                               "taxname": sampledf["taxname"],
                               "color":   sampledf["color"]
                               })
    # get all the level columns from the sampledf
    level_cols = [x for x in sampledf.columns if "level_" in x]
    # replace the missing values with the string ""
    sampledf[level_cols] = sampledf[level_cols].fillna("")
    for thiscol in level_cols:
        hover_data[thiscol] = sampledf[thiscol]
    hover_data = hover_data.fillna("")

    color_dict = {i: sampledf["color"][i] for i in sampledf.index}

    plot = umap.plot.interactive(mapper,
                                 color_key = color_dict,
                                 labels = sampledf["sample"],
                                 hover_data = hover_data,
                                 tools=[],
                                 point_size = 4
                                 )
    # add a title to the plot
    plot.title.text = plot_title
    # output to an HTML file
    bokeh.io.output_file(outhtml)
    # Save the plot to an HTML file
    bokeh.io.save(plot)

def filter_sample_df_by_clades(sampledf, taxids_to_keep, taxids_to_remove) -> pd.DataFrame:
    """
    This takes, as input, a sampledf and a list of taxids to keep. It returns a filtered sampledf.
    """
    # There is a column in the sample df called taxid_list.
    # This is a list of all the taxids in the lineage of the sample from closest to root to furthest.
    # Check to see if any of the taxid_to_keep are in the taxid_list. Return a df of the rows that match.
    sampledf = sampledf[sampledf["taxid_list"].apply(lambda x: any([y in taxids_to_keep for y in aliteraleval(x)]))]
    # We now remove the taxids that we know we don't want to keep.
    return sampledf[sampledf["taxid_list"].apply(lambda x: not any([y in taxids_to_remove for y in aliteraleval(x)]))]

def ALGrbh_to_algcomboix(rbhfile) -> dict:
    """
    Returns a dictionary of the unique ALG combinations to an index.
    """
    df = parse_rbh(rbhfile)
    alg_combo_to_ix = {tuple(sorted(x)): i
                            for i, x in enumerate(list(combinations(
                                df["rbh"], 2)))}
    return alg_combo_to_ix

def algcomboix_file_to_dict(ALGcomboixfile) -> dict:
    """
    This takes in a file that contains a dictionary of the unique ALG combinations to an index.
    It returns a dictionary.
    """
    if not os.path.exists(ALGcomboixfile):
        raise IOError(f"The file {ALGcomboixfile} does not exist. Exiting.")
    alg_combo_to_ix = {}
    with open(ALGcomboixfile, "r") as infile:
        for line in infile:
            line = line.strip()
            if line:
                key, value = line.split("\t")
                rbh1 = key.replace("(", "").replace(")", "").replace("'", "").replace(" ", "").split(",")[0]
                rbh2 = key.replace("(", "").replace(")", "").replace("'", "").replace(" ", "").split(",")[1]
                alg_combo_to_ix[tuple([rbh1, rbh2])] = int(value)
    return alg_combo_to_ix

def construct_coo_matrix_from_sampledf(sampledf, alg_combo_to_ix, print_prefix = ""):
    """
    This takes in a sampledf, and constructs a coo matrix of all of the distance matrices.
    """
    # take the first couple of keys from the alg_combo_to_ix and check that they are type tuple with two type strings
    counter = 0
    for key in alg_combo_to_ix:
        if not type(key) == tuple:
            raise ValueError(f"The key {key} is not a tuple. Exiting.")
        if not len(key) == 2:
            raise ValueError(f"The key {key} is not of length 2. Exiting.")
        if not all([type(x) == str for x in key]):
            raise ValueError(f"The key {key} is not of type string. Exiting.")
        counter += 1
        if counter == 5:
            break

    # check if the max index is greater than the length of the sampledf -1
    if max(sampledf.index) > len(sampledf) - 1:
        raise ValueError(f"The maximum index of the sampledf is greater than the length of the sampledf. Exiting.")
    # This is annoying to create this temporary data structure, but it helps us pinpoint broken distance matrix .gb.gz files.
    tempdfs = []
    for i, row in sampledf.iterrows():
        thisfile = row["dis_filepath_abs"]
        try:
            tempdfs.append(pd.read_csv(thisfile, sep = "\t", compression = "gzip"))
            # add a column with the index of the sampledf
            tempdfs[-1]["row_indices"] = i
            # assert that all of the values of rbh1 are less than the values of rbh2
            assert all(tempdfs[-1]["rbh1"] < tempdfs[-1]["rbh2"])
        except:
            raise IOError(f"The file {thisfile} could not be read in with pandas. There probably was something wrong wth the compression. Try deleteing this file. It will be regenerated. Exiting.")
    concatdf = pd.concat(tempdfs)
    start = time.time()
    concatdf["pair"] = concatdf.apply(lambda x: (x["rbh1"], x["rbh2"]), axis = 1)
    stop = time.time()
    print ("{}It took {} seconds to add the pair column with apply".format(print_prefix, stop - start))
    start = time.time()

    concatdf["col_indices"] = concatdf["pair"].map(alg_combo_to_ix)
    stop = time.time()
    print("{}It took {} seconds to add the col_indices column with map".format(print_prefix, stop - start))

    sparse_matrix = coo_matrix(( concatdf["distance"],
                                (concatdf["row_indices"], concatdf["col_indices"])),
                                shape = (len(sampledf), len(alg_combo_to_ix)))
    del concatdf
    return sparse_matrix

def construct_lil_matrix_from_sampledf(sampledf, alg_combo_to_ix, print_prefix = "") -> lil_matrix:
    """
    This takes a sampledf, and a directory of distance matrices, and returns a lil matrix.
    We return a lil matrix because it is a sparse representation of the matrix.

    rbh1                              rbh2                              distance
    Simakov2022BCnS_genefamily_12988  Simakov2022BCnS_genefamily_6122   10885675
    Simakov2022BCnS_genefamily_12988  Simakov2022BCnS_genefamily_7201   10538458
    Simakov2022BCnS_genefamily_6122   Simakov2022BCnS_genefamily_7201   347217
    Simakov2022BCnS_genefamily_12988  Simakov2022BCnS_genefamily_7465   8881006
    Simakov2022BCnS_genefamily_6122   Simakov2022BCnS_genefamily_7465   2004669
    Simakov2022BCnS_genefamily_7201   Simakov2022BCnS_genefamily_7465   1657452

    The way this works is that we take in the sampledf, reads in the distance matrices as pandas dataframes,
        and then constructs a lil matrix from the distance matrices. Just use the file path from the sampledf.
    """
    # check if the max index is greater than the length of the sampledf -1
    if max(sampledf.index) > len(sampledf) - 1:
        raise ValueError(f"The maximum index of the sampledf is greater than the length of the sampledf. Exiting.")
    # This is annoying to create this temporary data structure, but it helps us pinpoint broken distance matrix .gb.gz files.
    tempdfs = []
    for i, row in sampledf.iterrows():
        thisfile = row["dis_filepath_abs"]
        try:
            tempdfs.append(pd.read_csv(thisfile, sep = "\t", compression = "gzip"))
            # add a column with the index of the sampledf
            tempdfs[-1]["row_indices"] = i
            # assert that all of the values of rbh1 are less than the values of rbh2
            assert all(tempdfs[-1]["rbh1"] < tempdfs[-1]["rbh2"])
        except:
            raise IOError(f"The file {thisfile} could not be read in with pandas. There probably was something wrong wth the compression. Try deleteing this file. It will be regenerated. Exiting.")
    concatdf = pd.concat(tempdfs)
    start = time.time()
    # I require that all of the input
    concatdf["pair"] = concatdf.apply(lambda x: (x["rbh1"], x["rbh2"]), axis = 1)
    stop = time.time()
    print ("{}It took {} seconds to add the pair column with apply".format(print_prefix, stop - start))
    start = time.time()
    concatdf["col_indices"] = concatdf["pair"].map(alg_combo_to_ix)
    stop = time.time()
    print("{}It took {} seconds to add the col_indices column with map".format(print_prefix, stop - start))

    sparse_matrix = coo_matrix(( concatdf["distance"],
                                (concatdf["row_indices"], concatdf["col_indices"])),
                                shape = (len(sampledf), len(alg_combo_to_ix)))
    del concatdf
    sparse_matrix = sparse_matrix.tolil()
    return sparse_matrix

def rbh_to_samplename(rbhfile, ALGname) -> str:
    """
    This function takes an rbh filename and an ALG name and returns the sample name.
    # All of the filenames look like this:
      - BCnSSimakov2022_Zonotrichialeucophrys-44393-GCA028769735.1_xy_reciprocal_best_hits.plotted.rbh
    """
    filename = os.path.basename(rbhfile)
    # strip the ALGname and first _ from the front of the rbhfile
    # check that the filename starts with the ALGname
    if not filename.startswith(f"{ALGname}_"):
        raise ValueError(f"The filename {filename} does not start with the ALGname {ALGname}_. Exiting.")
    # split on the first ALGname and _
    filename = filename.split(f"{ALGname}_")[1]
    # try to split on _ and return the first element
    filename = filename.split("_")[0]
    # make sure that there are three fields
    splits = filename.split("-")
    if len(splits) != 3:
        raise ValueError(f"The filename {filename} does not have three fields when split on '-'. It splits to {splits}.")

    # split on - and check that the second element is an integer
    if not re.match(r"^[0-9]*$", splits[1]):
        raise ValueError(f"There is a non-numeric character in the taxid string, {taxid}, for file {thisfile}. Exiting.")
    # I haven't makde a unit test. Not working on ssh. Oh well.
    return filename

def plot_topoumap_from_files(sampledffile, ALGcomboixfile, coofile,
                             outdir, sample, smalllargeNaN, n_neighbors, min_dist,
                             taxids_to_keep, taxids_to_remove):
    """
    This function makes a UMAP plot where the points are inverted.
    The points for this are the distances between the pairs.
    The colors are the colors of the taxids.
    """
    # make sure that taxids_to_keep and taxids_to_remove are lists
    if not type(taxids_to_keep) == list:
        raise ValueError(f"The taxids_to_keep {taxids_to_keep} is not a list. Exiting.")
    if not type(taxids_to_remove) == list:
        raise ValueError(f"The taxids_to_remove {taxids_to_remove} is not a list. Exiting.")

    # make sure that everything in taxids_to_keep and taxids_to_remove are integers
    for entry in taxids_to_keep + taxids_to_remove:
        if not re.match(r"^[0-9]*$", str(entry)):
            raise ValueError(f"The taxid {entry} is not an integer. Exiting.")

    cdf = pd.read_csv(sampledffile, sep = "\t", index_col = 0)
    cdf2 = filter_sample_df_by_clades(cdf, taxids_to_keep, taxids_to_remove)
    # Get a list of the indices that are in cdf that are not in cdf2.
    ixnotin = [x for x in cdf.index if x not in cdf2.index]
    # These are the indices that we want to remove from the lil matrix.
    lil = load_npz(coofile).tolil()
    # we are removing the row indices that are not in cdf2
    lil = lil[[x for x in range(lil.shape[0]) if x not in ixnotin]]



def plot_umap_from_files(sampledffile, ALGcomboixfile, coofile,
                         outdir, sample, smalllargeNaN, n_neighbors, min_dist):
    """
    This is an all-in-one plotting method to make UMAP plots from the files.

    Later, I could add a feature to plot a specific clade if I remove the entries from the lil matrix that are
     not the samples that we want.
    """
    # read in the sample dataframe. We will need this later
    cdf = pd.read_csv(sampledffile, sep = "\t", index_col = 0)
    # Read in the ALGcomboixfile
    ALGcomboix = algcomboix_file_to_dict(ALGcomboixfile)
    lil = load_npz(coofile).tolil()

    # check that the largest row index of the lil matrix is less than the largest index of cdf - 1
    if lil.shape[0] > max(cdf.index) + 1:
        raise ValueError(f"The largest row index of the lil matrix, {lil.shape[0]}, is greater than the largest index of cdf, {max(cdf.index)}. Exiting.")
    # check that the largest value of the ALGcomboix is less than the number of columns of the lil matrix - 1
    if max(ALGcomboix.values()) > lil.shape[1] - 1:
        raise ValueError(f"The largest value of the ALGcomboix, {max(ALGcomboix.values())}, is greater than the number of columns of the lil matrix, {lil.shape[1]}. Exiting.")

    # If we pass these checks, we should be fine

    # check that the smalllargeNaN is either small or large
    if smalllargeNaN not in ["small", "large"]:
        raise ValueError(f"The smalllargeNaN {smalllargeNaN} is not 'small' or 'large'. Exiting.")
    if smalllargeNaN == "large":
        # we have to flip the values of the lil matrix
        lil.data[lil.data == 0] = 999999999999
    # check that min_dist is between 0 and 1
    if min_dist < 0 or min_dist > 1:
        raise IOError(f"The min_dist {min_dist} is not between 0 and 1. Exiting.")

    # We need a unique set of files for each of these
    # In every case, we must produce a .df file and a .bokeh.html file
    UMAPdf    = f"{outdir}/{sample}.neighbors_{n_neighbors}.mind_{min_dist}.missing_{smalllargeNaN}.df"
    UMAPbokeh = f"{outdir}/{sample}.neighbors_{n_neighbors}.mind_{min_dist}.missing_{smalllargeNaN}.bokeh.html"
    if len(cdf) <= n_neighbors:
        print(f"    The number of samples, {len(cdf)}, is less than the number of neighbors, {n_neighbors}. Skipping.")
        # write to empty UMAPdf and UMAPbokeh files
        with open(UMAPdf, "w") as f:
            f.write("")
        with open(UMAPbokeh, "w") as f:
            f.write("")
    elif len(cdf) > n_neighbors: # we have this condition for smaller datasets
        try:
            print(f"    PLOTTING - UMAP with {smalllargeNaN} missing vals, with n_neighbors = {n_neighbors}, and min_dist = {min_dist}")
            reducer = umap.UMAP(low_memory=True, n_neighbors = n_neighbors, min_dist = min_dist)
            start = time.time()
            mapper = reducer.fit(lil)
            stop = time.time()
            print("   - It took {} seconds to fit_transform the UMAP".format(stop - start))
            # save the UMAP as a bokeh plot
            umap_mapper_to_bokeh(mapper, cdf, UMAPbokeh,
                plot_title = f"UMAP of {sample} with {smalllargeNaN} missing vals, n_neighbors = {n_neighbors}, min_dist = {min_dist}")
            umap_df = umap_mapper_to_df(mapper, cdf)
            umap_df.to_csv(UMAPdf, sep = "\t", index = True)
            # save the connectivity figure
            UMAPconnectivity = f"{outdir}/{sample}.neighbors_{n_neighbors}.mind_{min_dist}.missing_{smalllargeNaN}.connectivity.jpeg"
            UMAPconnectivit2 = f"{outdir}/{sample}.neighbors_{n_neighbors}.mind_{min_dist}.missing_{smalllargeNaN}.connectivity2.jpeg"
            try:
                umap_mapper_to_connectivity(mapper, UMAPconnectivity,
                                            title = f"UMAP of {sample} with {smalllargeNaN} missing vals, n_neighbors = {n_neighbors}, min_dist = {min_dist}")
            except:
                print(f"    Warning: Could not make the connectivity plot for {UMAPconnectivity}")
            try:
                umap_mapper_to_connectivity(mapper, UMAPconnectivit2, bundled = True,
                                            title = f"UMAP of {sample} with {smalllargeNaN} missing vals, n_neighbors = {n_neighbors}, min_dist = {min_dist}")
            except:
                print(f"    Warning: Could not make the connectivity plot for {UMAPconnectivit2}")
        except UserWarning as e:
            # Catch the specific warning about graph not being fully connected
            if "Graph is not fully connected" in str(e):
                print("    Warning: Graph is not fully connected. Can't run UMAP with these parameters.")
                # we check for file UMAPbokeh, so write this message to it
                with open(UMAPbokeh, "w") as f:
                    f.write("The graph is not fully connected. Can't run UMAP with these parameters.")
                # write an empty .df file
                with open(UMAPdf, "w") as f:
                    f.write("")
            else:
                # If it's a different warning, re-raise it
                raise e

def main():
    """
    This program interprets many RBH files, and handles the way that they are visualized.

    Steps:
      1. Generate a sampledf from the RBH files. This dataframe contains important information.
        - Columns:
          - index: The index of the dataframe is important because this will be the order of the samples in the distance matrix. Everything will be ordered by this index.
          - sample: This is the sample name. This is the same sample information that will be in the rbh file columns, and in the distance matrix.
          - color: This is the color of the sample. This is based on a custom annotation from plot_ALG_fusions_v3.py
          - taxid: This is the NCBI TAXid of the sample.
          - taxname: This is the name of the taxid. For example 7777 is "Condrichthyes".
          - taxid_list: This is a list of all the taxids in the lineage of the sample from closest to root to furthest.
          - taxid_list_str: A string version of taxid_list joined together with ';' characters
          - taxname_list: This is a list of all the taxnames in the lineage of the sample from closest to root to furthest. Matches the indices of taxid_list.
          - taxname_list_str: A string version of taxname_list joined together with ';' characters
          - level_1: Here, the level_1, level_2, et cetera are splits of the NCBI taxid for easy plotting. These go up to level_10.
              - Some examples are:
                ```
                level_1: root (1); cellular organisms (131567); Eukaryota (2759); Opisthokonta (33154)
                level_2: Metazoa (2); Eumetazoa (6072); Bilateria (33213); Protostomia (33317)
                et cetera
                ```
          - printstring: This is the printstring of the taxid, like "root (1); cellular organisms (131567); Eukaryota (2759); Opisthokonta (33154)"
          - #NOT YET annotation_method: This shows what method was used to annotate the sample, for example "BCnSSimakov2022_protmap" or "NCBI/Existing"
          - genome_size: This is the size of the genome in base pairs.
          - #NOT YET gc_content: This is the GC content of the genome. This is expressed as a value from 0 to 1.
          - number_of_chromosomes: The haploid number of chromosomes that this species has.
          - filepath: The filepath, as provided, of the rbh file that was used to generate this information.
          - filepath_abs: The filepath, resolved, of the rbh file that was used to generate this information.
          - filename: The filename of the rbh file that was used to generate this information.
      2. Read in the RBH files and convert them to distance matrix files
      (Steps 1 and 2 are combined in the function rbh_directory_to_distance_matrix)
      3. Generate the LIL matrix of the data. Save this to disk.
      4. Generate UMAP projections of the data. In this step we use different parameters to test how they affect the UMAP.
      5. Visualize the data with plotly, bokeh, and matplotlib. For these plots, also generate QC plots.
    """
    args = parse_args()

    results_base_directory = "GTUMAP"
    # ┏┓┓ ┏┓  ┏┓┏┓┏┓┏┓┳┏┓┏┓  •  ┏
    # ┣┫┃ ┃┓━━┗┓┃┃┣ ┃ ┃┣ ┗┓  ┓┏┓╋┏┓
    # ┛┗┗┛┗┛  ┗┛┣┛┗┛┗┛┻┗┛┗┛  ┗┛┗┛┗┛
    rbh_files = list(sorted([os.path.join(args.directory, f)
                 for f in os.listdir(args.directory)
                 if f.endswith('.rbh')], reverse = True))
    # Generate the distance matrices. This is the part that saves the .gb.gz files.
    outtsv = f"{results_base_directory}/sampledf.tsv"
    outputdir = f"{results_base_directory}/distance_matrices/"
    sampledf = None
    if not os.path.exists(outtsv):
        print("generating the distance matrices")
        sampledf = rbh_directory_to_distance_matrix(args.directory, args.ALGname, outtsv = outtsv, outputdir = outputdir)
    else:
        print("The sampledf file already exists. Not overwriting it. Instead, we will read it in.")
        print("If you would rather overwrite it, delete it in the file system and run this program again.")
        sampledf = pd.read_csv(outtsv, sep = "\t", index_col = 0)

    # We will need to calculate the rbh combo to index dictionary.
    # DO NOT bother reading in the existing file. It takes 5x longer
    #  to read in the file than it does to generate it.
    outfile = f"{results_base_directory}/combo_to_index.txt"
    alg_combo_to_ix = ALGrbh_to_algcomboix(args.rbhfile)
    # save the dictionary pairs
    with open(outfile, "w") as f:
        for key, value in alg_combo_to_ix.items():
            f.write(f"{key}\t{value}\n")

    ncbi = NCBITaxa()
    # Now we must construct a lil matrix from the distance matrices. This takes up significant disk space.
    # We now parse through different sets of species and parameters to see how the UMAP is affected.
    # Usage notes:
    # For about 100 samples for molluscs, took up about 29Gb virtual memory and about 8Gb of RAM.
    for thissp in args.cladestoplot:
        thiscladename = ncbi.get_taxid_translator([thissp])[thissp]
        sample_outfix = f"{thiscladename}_{thissp}"
        sampleoutdir = f"{results_base_directory}/analyses/{sample_outfix}"
        # safely make the dirs for the output files
        create_directories_if_not_exist(sampleoutdir)
        print(f"We are making a plot for this NCBI taxid: {thissp}")


        # LIL FILE
        # This whole thing could be sped up for repeated runs if I allowed the LIL matrix to be saved to disk.
        # Right now it is not saved to disk, because it is a large file.
        # This wouldn't be such a big deal if I compressed the output. I need to optimize this to make
        #  sure that I understand the inputs for constructing the LIL matrix object, though.
        print(f"  - Constructing the LIL matrix")
        samplecdf1 = f"{sampleoutdir}/{sample_outfix}.sampledf.tsv"
        samplecdf2 = f"{sampleoutdir}/{sample_outfix}.sampledf.matrixindices.tsv"
        cdf = filter_sample_df_by_clades(sampledf, [thissp])
        # Save the cdf as a tsv
        cdf.to_csv(samplecdf1, sep = "\t", index = True)
        # reset the indices, because that is what we need for the lil matrix
        cdf = cdf.reset_index(drop = True)
        cdf.to_csv(samplecdf2, sep = "\t", index = True)
        # now get the lil matrix
        lil = construct_lil_matrix_from_sampledf(cdf, alg_combo_to_ix, print_prefix = "  - ")

        # UMAP section
        # Levels for plotting text:
        # We are making a plot for this NCBI taxid: 33208
        #  - Constructing the LIL matrix for {sample_outfix}
        for missing in ["small", "large"]:
            if missing == "large":
                # we have to flip the values of the lil matrix
                lil.data[lil.data == 0] = 999999999999
            for n in [2, 5, 10, 20, 50, 100, 250]: # this is the number of neighbors
                for min_dist in [0.0, 0.1, 0.2, 0.5, 0.75, 0.9]:
                    # We need a unique set of files for each of these
                    if len(cdf) <= n:
                        print(f"    The number of samples, {len(cdf)}, is less than the number of neighbors, {n}. Skipping.")
                    if len(cdf) > n: # we have this condition for smaller datasets
                        # First check if the UMAP exists
                        UMAPdf           = f"{sampleoutdir}/{sample_outfix}.neighbors_{n}.mind_{min_dist}.missing_{missing}.df"
                        UMAPbokeh        = f"{sampleoutdir}/{sample_outfix}.neighbors_{n}.mind_{min_dist}.missing_{missing}.bokeh.html"
                        UMAPfit = None
                        if os.path.exists(UMAPbokeh):
                            print(f"    FOUND - UMAP with {missing} missing vals, with n_neighbors = {n}, and min_dist = {min_dist}")
                        else:
                            try:
                                print(f"    PLOTTING - UMAP with {missing} missing vals, with n_neighbors = {n}, and min_dist = {min_dist}")
                                reducer = umap.UMAP(low_memory=True, n_neighbors = n, min_dist = min_dist)
                                start = time.time()
                                mapper = reducer.fit(lil)
                                stop = time.time()
                                print("   - It took {} seconds to fit_transform the UMAP".format(stop - start))
                                # save the UMAP as a bokeh plot
                                umap_mapper_to_bokeh(mapper, cdf, UMAPbokeh,
                                    plot_title = f"UMAP of {sample_outfix} with {missing} missing vals, n_neighbors = {n}, min_dist = {min_dist}")
                                umap_df = umap_mapper_to_df(mapper, cdf)
                                umap_df.to_csv(UMAPdf, sep = "\t", index = True)
                                # save the connectivity figure
                                if args.qualitycontrol:
                                    UMAPconnectivity = f"{sampleoutdir}/{sample_outfix}.neighbors_{n}.mind_{min_dist}.missing_{missing}.connectivity.jpeg"
                                    UMAPconnectivit2 = f"{sampleoutdir}/{sample_outfix}.neighbors_{n}.mind_{min_dist}.missing_{missing}.connectivity2.jpeg"
                                    try:
                                        umap_mapper_to_connectivity(mapper, UMAPconnectivity,
                                                                    title = f"UMAP of {sample_outfix} with {missing} missing vals, n_neighbors = {n}, min_dist = {min_dist}")
                                    except:
                                        print(f"    Warning: Could not make the connectivity plot for {UMAPconnectivity}")
                                    try:
                                        umap_mapper_to_connectivity(mapper, UMAPconnectivit2, bundled = True,
                                                                    title = f"UMAP of {sample_outfix} with {missing} missing vals, n_neighbors = {n}, min_dist = {min_dist}")
                                    except:
                                        print(f"    Warning: Could not make the connectivity plot for {UMAPconnectivit2}")
                                    #QCplot = f"{sampleoutdir}/{sample_outfix}.neighbors_{n}.mind_{min_dist}.missing_{missing}.QC.jpeg"
                                    ##try:
                                    #umap_mapper_to_QC_plots(mapper, QCplot,
                                    #                title = f"UMAP of {sample_outfix} with {missing} missing vals, n_neighbors = {n}, min_dist = {min_dist}")
                                    ##except:
                                    ##    print(f"    Warning: failed to make the QC plot for {QCplot}.")
                            except UserWarning as e:
                                # Catch the specific warning about graph not being fully connected
                                if "Graph is not fully connected" in str(e):
                                    print("    Warning: Graph is not fully connected. Can't run UMAP with these parameters.")
                                    # we check for file UMAPbokeh, so write this message to it
                                    with open(UMAPbokeh, "w") as f:
                                        f.write("The graph is not fully connected. Can't run UMAP with these parameters.")
                                else:
                                    # If it's a different warning, re-raise it
                                    raise e


    ## ┏┓┓   ┓           •   ┏┳┓       ┏  ┏┓┓ ┏┓  ┳┓┳┓┓┏  ┏•┓
    ## ┃┃┣┓┓┏┃┏┓┏┓┏┓┏┓┏┓╋┓┏   ┃ ┏┓┏┓┏┓ ┃  ┣┫┃ ┃┓  ┣┫┣┫┣┫  ╋┓┃┏┓
    ## ┣┛┛┗┗┫┗┗┛┗┫┗ ┛┗┗ ┗┗┗   ┻ ┛ ┗ ┗  ┛  ┛┗┗┛┗┛  ┛┗┻┛┛┗  ┛┗┗┗
    ##      ┛    ┛
    ## construct our dataframe
    #T = PhyloTree()
    ## Add in the rbh information
    #T.ingest_ALG_rbh(args.ALGname, args.rbhfile)
    ## for all the files in the results directory, we will add the distances to the PhyloTree
    ##for thisfile in [x for x in os.listdir("results") if x.endswith(".gb.gz")][:100]:
    #counter = 0
    ## This part uses 16.5GB of memory for 3600 files
    #gbgzfiles = [x for x in os.listdir("results") if x.endswith(".gb.gz")]
    ## check that for every .gz.gz file, the second things is an integer
    #for thisfile in gbgzfiles:
    #    thissample = thisfile.replace(".gb.gz", "")
    #    taxid = thissample.split("-")[1]
    #    if not re.match(r"^[0-9]*$", taxid):
    #        raise ValueError(f"There is a non-numeric character in the taxid string, {taxid}, for file {thisfile}. Exiting.")
    #    if int(taxid) not in taxid_to_taxidstring:
    #        raise ValueError(f"The taxid {taxid} is not in the taxid_to_taxidstring dictionary. The file was {thisfile}. Exiting.")
    ## Process the .gb.gz files and add them to the graph
    #for thisfile in gbgzfiles:
    #    print(f"\r    Adding the file {counter}/{len(gbgzfiles)}", end="", file = sys.stdout)
    #    thissample = thisfile.replace(".gb.gz", "")
    #    taxid = int(thissample.split("-")[1])
    #    try:
    #        distdf = pd.read_csv(f"results/{thisfile}", sep = "\t", compression = "gzip")
    #    except:
    #        raise IOError(f"The file {thisfile} could not be read in with pandas. There probably was something wrong wth the compression. Try deleteing this file. It will be regenerated. Exiting.")
    #    T.add_lineage_string_sample_distances(taxid_to_taxidstring[taxid],
    #                                          thissample, args.ALGname,
    #                                          distdf)
    #    counter += 1

    #print()
    #print("Done adding the files")
    #T.merge_sampledistances_to_locdf()
    #print("done")
    ##print(T.G)

if __name__ == "__main__":
    main()