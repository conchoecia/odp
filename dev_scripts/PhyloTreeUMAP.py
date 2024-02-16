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
import bokeh           # bokeh is used to visualize and save the UMAP
import networkx as nx
import numpy as np
import os
import pandas as pd
import pickle
import re
from scipy.sparse import lil_matrix
from scipy.sparse import coo_matrix
import sys
import time
import umap
import umap.plot
import warnings
warnings.filterwarnings("ignore", message="Graph is not fully connected", category=UserWarning)

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

        sparse_matrix = sparse_matrix.tolil()
        # set the missing values of the sparse matrix to 999999999999
        sparse_matrix.data[sparse_matrix.data == 0] = 999999999999

        print("The type of the sparse matrix is ", type(sparse_matrix))
        print("Fitting the UMAP")
        reducer = umap.UMAP(low_memory=True)
        start = time.time()
        mapper = reducer.fit(sparse_matrix)
        stop = time.time()
        print("It took {} seconds to fit_transform the UMAP".format(stop - start))

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

def rbh_directory_to_distance_matrix(rbh_directory, ALGname):
    """
    Takes all of the rbh files in the directory and calculates the distance matrix of all of the BCnS ALGs.
    Saves all of them in a directory called results. They
    """
    # make sure the results directory exists
    # if the directory doesn't exist yet, make it
    if not os.path.exists("results"):
        os.makedirs("results")

    # get the rbh files in the directory
    rbh_files = list(sorted([os.path.join(rbh_directory, f)
                 for f in os.listdir(rbh_directory)
                 if f.endswith('.rbh')], reverse = True))
    #rbh_files = rbh_files[:100]

    taxid_to_taxidstring = taxids_to_taxidstringdict(
        [int(os.path.basename(x).split('-')[1]) for x in rbh_files])
#
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
        taxidstring = taxid_to_taxidstring[taxid]
        df = parse_rbh(rbhfile)
        # make sure that ALG "_scaf", "_gene", and "_pos" are in the columns.
        for col in ["_scaf", "_gene", "_pos"]:
            thiscol = f"{ALGname}{col}"
            if thiscol not in df.columns:
                raise IOError(f"The column {thiscol} is not in the rbh file {rbhfile}. Exiting.")
        thissample = [x for x in df.columns
                      if "_scaf" in x
                      and ALGname not in x][0].split("_")[0]
        gb_filepath = f"results/{thissample}.gb.gz"
        if not os.path.exists(gb_filepath):
            rbh_to_gb(thissample, df, gb_filepath)
    print("Done parsing the rbh files")

def parse_args():
    """
    This has all of the args that we need to parse.
    The args that we need to parse are:
      -d --directory : The directory that contains the RBH files. These RBH files
      -V --overwrite : If this is set, then we will overwrite the files in the output directory.
                       Otherwise, try to load the existing files. This is not required.
      -a --ALGname   : The name of the ALG that we are looking at. This is required.
      -r --rbhfile   : The ALG rbh file. This is required, as we will use this information later.
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
    args = parser.parse_args()

    # check that the directory exists
    if not os.path.exists(args.directory):
        raise IOError(f"The directory {args.directory} does not exist. Exiting.")
    return args

def main():
    args = parse_args()

    rbh_files = list(sorted([os.path.join(args.directory, f)
                 for f in os.listdir(args.directory)
                 if f.endswith('.rbh')], reverse = True))
    taxid_to_taxidstring = taxids_to_taxidstringdict(
        [int(os.path.basename(x).split('-')[1]) for x in rbh_files])

    # generate the distance matrices
    print("generating the distance matrices")
    rbh_directory_to_distance_matrix(args.directory, args.ALGname)

    # construct our dataframe
    T = PhyloTree()
    # Add in the rbh information
    T.ingest_ALG_rbh(args.ALGname, args.rbhfile)
    # for all the files in the results directory, we will add the distances to the PhyloTree
    #for thisfile in [x for x in os.listdir("results") if x.endswith(".gb.gz")][:100]:
    counter = 0
    gbgzfiles = [x for x in os.listdir("results") if x.endswith(".gb.gz")]
    for thisfile in gbgzfiles:
        print(f"\r    Adding the file {counter}/{len(gbgzfiles)}", end="", file = sys.stdout)
        thissample = thisfile.replace(".gb.gz", "")
        taxid = int(thissample.split("-")[1])
        distdf = pd.read_csv(f"results/{thisfile}", sep = "\t", compression = "gzip")
        T.add_lineage_string_sample_distances(taxid_to_taxidstring[taxid],
                                              thissample, args.ALGname,
                                              distdf)
        counter += 1
    print()
    print("Done adding the files")
    T.merge_sampledistances_to_locdf()
    print("done")
    #print(T.G)

if __name__ == "__main__":
    main()