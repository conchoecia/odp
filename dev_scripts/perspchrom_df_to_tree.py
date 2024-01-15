#!/usr/bin/env python3
# the file path is the first positional arg. use sys.
import ast
import os
from typing import Any
import numpy as np
import pandas as pd
import random
import sys
from ete3 import NCBITaxa

# import the parse_rbh_file from the plot_ALG_fusions.py script
# this is a function that parses the RBH file into a dataframe
from plot_ALG_fusions import parse_rbh_file
import matplotlib.pyplot as plt
# import patches
import matplotlib.patches as mpatches
# use twoslope norm to make a diverging color map
from matplotlib.colors import TwoSlopeNorm

# Filter out the DeprecationWarning related to the py23 module
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning, module="fontTools.misc.py23")
# import pdfpages
from matplotlib.backends.backend_pdf import PdfPages


# use networkx to make graphs for the lineage-specific fusion/losses
import networkx as nx

def parse_gain_loss_string(GL_string, samplename) -> pd.DataFrame:
    """
    This function parses a single gain/loss string.
    This is one example gain/loss string from the cnidarian Xenia:
      1-([]|[])-131567-([]|[])-2759-([]|[])-33154-([]|[])-33208-([]|[])-6072-([('A1b', 'B3'), ('A2', 'N'), ('B1', 'B2'), ('Eb', 'Qb'), ('F', 'Qb'), ('O2', 'Qd'), ('Qc', 'R')]|[])-6073-([]|[])-6101-([('A1a', 'K'), ('A1b', 'Qc'), ('A1b', 'R'), ('A2', 'Eb'), ('A2', 'F'), ('A2', 'Qb'), ('B2', 'J2'), ('B2', 'L'), ('B3', 'Qc'), ('B3', 'R'), ('C1', 'J1'), ('C1', 'Qa'), ('C2', 'J2'), ('Ea', 'G'), ('Eb', 'F'), ('Eb', 'N'), ('F', 'N'), ('I', 'P'), ('J1', 'Qa'), ('J2', 'L'), ('M', 'O2'), ('M', 'Qd'), ('N', 'Qb')]|[])-6132-([]|[])-3028843-([]|[])-86538-([]|[])-86539-([]|[])-2653460-([]|[])-2897299
    In the above example, if we ignore the colocalizations/losses, we get the following taxid string:
      1-131567-2759-33154-33208-6072-6073-6101-6132-3028843-86538-86539-2653460-2897299

    The string is formatted like so:
        taxid-([colocalization list]|[ALG dispersion list])-taxid-([colocalization list].... etc.

    The colocalizations and ALG dispersions are placed between the ALG taxids, because these chages
     happen on the phylogenetic branches between the NCBI taxids.

    Args:
      - GL_string: a string that contains the taxids and the changes that happen on the branches between them.
          Look above for the format
      - samplename: the name of the sample that this string came from.


    Return:
      We will return a dataframe that contains the edges on which the changes happened (taxid, taxid),
         and the types of change on the branch. Literally, the columns of the return dataframe are:
           - source_taxid
           - target_taxid
           - colocalizations
           - losses
    """
    splitstring = GL_string.split("-")
    # Go every two elements in the list.
    entries = []
    # These will be the taxids, and the things between them are the changes.
    for i in range(0,len(splitstring)-1, 2):
        taxid1 = int(splitstring[i])
        taxid2 = int(splitstring[i+2])
        colocloss = splitstring[i+1]
        # interpret the next two strings as lists
        colocs = eval(      colocloss.lstrip("(").split("]|[")[0] + "]")
        losses = eval("[" + colocloss.rstrip(")").split("]|[")[1]      )
        entry = {"source_taxid": taxid1,    # already an int, we cast it earlier
                 "target_taxid": taxid2,    # already an int, we cast it earlier
                 "colocalizations": colocs, # this is a list
                 "losses": losses,          # this is a list
                 "samplename": samplename,  # this is a string
                 "sample_taxid": int(GL_string.split("-")[-1])
                 } # this is an int
        entries.append(entry)
    # convert the list of dicts to a dataframe
    df = pd.DataFrame(entries)
    return df

def parse_gain_loss_from_perspchrom_df(perspchromdf) -> pd.DataFrame:
    """
    This handles parsing the gain/loss strings from the perspchrom dataframe.
    It handles the whole dataframe, and in the end outputs a dataframe that notes
    all of the chromosomal changes from the ALGs.

    Args:
      - perspchromdf: a dataframe that contains the ALG presence and colocalization
        information for each species. This dataframe is the output of plot_ALG_fusions.py

    Return:
      - a dataframe that contains the chromosomal changes between the ALGs. This dataframe
        has the following columns:
        - source_taxid
        - target_taxid
        - colocalizations
        - losses
        - samplename
        - sample_taxid
    """
    df_list = []
    # for the changestrings column, parse the coloc/loss string
    for i, row in perspchromdf.iterrows():
        spdf = parse_gain_loss_string(row["changestrings"], row["species"])
        df_list.append(spdf)
    # concatenate the list of dataframes into one big dataframe
    changedf = pd.concat(df_list)
    return changedf

def stats_on_changedf(sampledf, changedf) -> pd.DataFrame:
    """
    This function performs some basic stats on the collated list of changes from all the samples, called changedf,
      and returns a dataframe that contains information about where different changes occurred:

    Inputs:
      - sampledf: A dataframe that contains the ALG presence and colocalization information,
                  as well as the taxidstring for each sample. We need this to parse the fraction
                  of samples within a clade that provided evidence for a specific change.
      - changedf: A dataframe that contains the changes that happened on each branch. This is the output
                  of the function parse_gain_loss_from_perspchrom_df.
    Outputs:
      - A pandas
    """
    # we want a column in the dataframe that contains the taxidstring parsed out to a list of ints
    sampledf["taxidstring"] = sampledf["taxidstring"].apply(lambda x: [int(taxid) for taxid in x.split(";")])
    # For each change on each branch, we want to know how many samples with that NCBI taxid have that change recorded.
    #  We count how many samples there are for each unique taxid in our dataframe, then use this later.
    #  Just do a for loop because it is n complexity.
    taxid_to_sample_count = {}
    for i, row in sampledf.iterrows():
        for this_taxid in row["taxidstring"]:
            if this_taxid not in taxid_to_sample_count:
                taxid_to_sample_count[this_taxid] = 0
            taxid_to_sample_count[this_taxid] += 1

    # Right now colocalizations and losses are lists. We want to count the number of colocalizations and losses.
    # So we will make a new column called 'change' and another called 'change_type', then we will unwrap
    # the lists into the new columns. This will make it easier to do stats on the changes.
    entries = []
    for i, row in changedf.iterrows():
        # colocalizations
        for colocalization in row["colocalizations"]:
            entry = {"source_taxid": row["source_taxid"],
                     "target_taxid": row["target_taxid"],
                     "change":       colocalization,
                     "change_type":  "colocalization",
                     "samplename":   row["samplename"],
                     "sample_taxid": row["sample_taxid"]}
            entries.append(entry)
        # losses
        for loss in row["losses"]:
            entry = {"source_taxid": row["source_taxid"],
                     "target_taxid": row["target_taxid"],
                     "change":       loss,
                     "change_type":  "loss",
                     "samplename":   row["samplename"],
                     "sample_taxid": row["sample_taxid"]}
            entries.append(entry)
    # convert the list of dicts to a dataframe
    changedf = pd.DataFrame(entries)

    # Because this is a big structured N-sat problem to figure out the exact branch on which the change happened,
    #  we will now look at whether the change happened right after the source_taxid or right before the
    #  target_taxid.
    # To do this we will groupby on source_taxid, then count the number of changes. We will then groupby on
    #  target_taxid, then count the number of changes. We will then groupby on both source_taxid and target_taxid,
    #  then count the number of changes. For now all we can do is compare these numbers.
    # For each case, we want to independently count the losses or colocalizations.
    # Just keep the changes column in the groupby after count, then change colocalization and loss to their own columns.
    # get right of target_taxid, samplename, sample_taxid
    groupby_target = changedf.groupby(["target_taxid", "change"])
    # sort by, then print, the most common changes
    groupby_target = groupby_target.count().sort_values(by="change_type", ascending=False)
    # change the the source_taxid column to a "counts" column
    groupby_target = groupby_target.rename(columns={"source_taxid": "counts"})
    # remove the change_type, samplename, sample_taxid columns
    groupby_target = groupby_target.drop(["change_type", "samplename", "sample_taxid"], axis=1)
    # unfold this into a dataframe
    groupby_target = groupby_target.reset_index()
    # In the change column, the same object can occur multiple times. Make a frac_total that counts how large this finding is as a fraction of the total for that change.
    # get the total number of each change type, as a dict
    change_to_total_counts = groupby_target.groupby("change")["counts"].sum().to_dict()
    # make a column num_samples_in_taxid by maping the target_taxid to the total number of samples with that taxid with the dict taxid_to_sample_count
    groupby_target["number_of_samples_in_taxid"] = groupby_target["target_taxid"].map(taxid_to_sample_count)
    groupby_target["change_frac_total"] = groupby_target.apply(lambda row: row["counts"]/change_to_total_counts[row["change"]], axis=1)
    groupby_target["frac_samples_w_this_taxid_w_this_change"] = groupby_target.apply(lambda row: row["counts"]/taxid_to_sample_count[row["target_taxid"]], axis=1)
    # use ete3 to get the names of the taxids
    ncbi = NCBITaxa()
    # get the names of the target taxids
    target_taxid_names = ncbi.get_taxid_translator(groupby_target["target_taxid"].tolist())
    # add a column to the dataframe that contains the names of the target taxids
    groupby_target["target_taxid_name"] = groupby_target["target_taxid"].map(target_taxid_names)
    return groupby_target

def delete_node_resize(G, node):
    """
    This function deletes a node from a graph.
    The way that it does this is by:
      - getting the node's size, then ranomly adding that many counts to other nodes.
    We ignore the color, and now it is just deleted

    Returns: a graph with the node deleted, and with the counts added to other nodes.
    """
    # get the size of the one node
    nodesize = G.nodes[node]["size"]
    # remove the node from the graph
    G.remove_node(node)
    # now randomly add the counts to the other nodes one at a time.
    node_list = list(G.nodes)
    if len(node_list) > 0:
        for i in range(nodesize):
            # randomly choose a node
            randnode = random.choice(node_list)
            # add one to the size of that node
            G.nodes[randnode]["size"] += 1
    return G

def node_size_fraction_of_total_size(G, node):
    """
    Returns the node size as a fraction of the total size of the graph.
    """
    return G.nodes[node]["size"]/sum([G.nodes[node]["size"] for node in G.nodes])

def node_size_fraction_of_largest(G, node):
    """
    Returns the node size as a fraction of the largest node in the graph.
    """
    return G.nodes[node]["size"]/max([G.nodes[node]["size"] for node in G.nodes])

def node_size_CC(G, node):
    """
    Returns the size of the CC in which the node is.
    """
    # get the CCs of the graph
    CCs = nx.connected_components(G)
    CC_size = -1
    # get the size of the CC in which the node is
    for CC in CCs:
        if node in CC:
            return sum([G.nodes[node]["size"] for node in CC])

def node_size_CC_fraction_of_total_size(G, node):
    """
    Get the CCs of the graph, figure out the size of the CC in which the node is,
    report the size of that CC as a fraction of the total size of the graph.
    """
    # get the CCs of the graph
    CCs = nx.connected_components(G)
    CC_size = -1
    # get the size of the CC in which the node is
    for CC in CCs:
        if node in CC:
            CC_size = sum([G.nodes[node]["size"] for node in CC])
    # return the size of the CC as a fraction of the total size of the graph
    return CC_size/sum([G.nodes[node]["size"] for node in G.nodes])

def node_size_CC_fraction_of_largest(G, node):
    """
    Get the CCs of the graph, get their sizes, find the largest.
    Report the size of the CC in which this node exists as a fraction of the largest CC.
    """
    # get the CCs of the graph
    CCs = nx.connected_components(G)
    node_cc_size = -1
    largest_cc_size = -1
    # get the size of each CC
    for CC in CCs:
        this_cc_size = sum([G.nodes[node]["size"] for node in CC])
        if this_cc_size > largest_cc_size:
            largest_cc_size = this_cc_size
        if node in CC:
            node_cc_size = this_cc_size
    return node_cc_size/largest_cc_size

def nodes_in_same_CC(G, node_iterable):
    """
    This returns True if the nodes in the iterable are all in the same CC in node G.
      Returns False otherwise.
    """
    # get the CCs of the graph
    CCs = nx.connected_components(G)
    for CC in CCs:
        # if all of the nodes are in this CC, return True
        if all([node in CC for node in node_iterable]):
            return True
    return False

def colocalize_these_nodes(G, node_iterable):
    """
    This function colocalizes the nodes in the iterable in the graph G.
    At the end of this function, every connected component will be a k-complete component.
      - Every node in the CC will be connected to every other node in the CC.
    """
    # First make the connections between all of the nodes in the iterable
    for node1 in node_iterable:
        for node2 in node_iterable:
            if node1 != node2:
                G.add_edge(node1, node2)
    # Now go through each connected component, and make it a k-complete component.
    #  This means that every node in the CC will be connected to every other node in the CC.
    CCs = nx.connected_components(G)
    for CC in CCs:
        # get all of the pairs of nodes in the CC
        node_pairs = [(node1, node2) for node1 in CC for node2 in CC if node1 != node2]
        # Add edges between all of the nodes in the CC
        # Fine if the edge already exists
        G.add_edges_from(node_pairs)
    return G

def stats_df_to_loss_fusion_dfs(perspchromdf, ALGdf, randomize_ALGs = False):
    """
    Use the perspchromdf, the statsdf, and an RBH file to make loss/fusion plots.
    The question is whether larger or smaller ALGs tend to fuse or be lost more often.

    The way that we do this is by going through each species, and analyzing the changes in
        ALGs on each lineage. If a change has been observed already on a lineage, then that
        count is not double counted. TODO This may be from weighting the changes, or from
        not incrementing the counter.

    The final matrix will show on one axis, the smaller ALG contributing to a fusion, and on
        the other axis, the larger ALG contributing to a fusion. Doing this in a per-species
        fashion will allow us to track the existing fusions to more accurately consider the
        changes at each evolutionary timepoint.

    To show what is the most prevalent pattern, we will use Pearson correlation of the matrix,
        compared to obs/exp.

    Input:
      - perspchromdf:   A dataframe that contains the ALG presence and colocalization. This is just the
                          input to the whole script.
      - rbhfile:        The rbh file that contains all of the info about the ALGs that will be studied
                          for this analysis.
      - randomize_ALGs: A boolean that determines whether the sizes of the ALGs will be randomized or not.
                          Running many randomizations will allow us to get a null distribution of the
                          fusion stats on small vs large ALGs.
    Output:
      - dispersion_df:  A dataframe that contains the information about the ALG losses.
      - coloc_df:       A dataframe that contains the information about the ALG colocalizations.
    """
    ALG_info = ALGdf.copy()
    if randomize_ALGs:
        # randomize the ALG sizes
        ALG_info["Size"] = ALG_info["Size"].sample(frac=1).reset_index(drop=True)
    # Make a starting graph where each "ALGname" in ALG_info is a node,
    #  each node has a size attribute (integer), and a color attribute (string)
    G = nx.Graph()
    # add the nodes
    for i, row in ALG_info.iterrows():
        G.add_node(row["ALGname"], size=int(row["Size"]), color=row["Color"])

    # Every time we run this program we will randomize the row order. Doing this will sample the different
    #  predicted events depending on the order in which we add fusions or losses to the graph.
    perspchromdf = perspchromdf.sample(frac=1)

    already_counted        = {}
    dispersion_entries     = []
    colocalization_entries = []
    for i, row in perspchromdf.iterrows():
        # For this species, make a copy of the starting ALG graph. We will manipulate this graph. to track the changes.
        thisG = G.copy()
        # parse the gain/loss string
        changedf = parse_gain_loss_string(row["changestrings"], row["species"])

        for j, row_change in changedf.iterrows():
            thisedge = (row_change["source_taxid"], row_change["target_taxid"])
            # first check if there has been either a fusion or a loss here. If so, add an entry to the already_counted dict.
            if row_change["colocalizations"] != [] or row_change["losses"] != []:
                if thisedge not in already_counted:
                    already_counted[thisedge] = {"colocalizations": [], "losses": []}
            # first we check if there have been any ALG dispersions
            if row_change["losses"] != []:
                # iterate through each loss event
                for thisloss in row_change["losses"]:
                    # check if this loss has already been counted
                    if thisloss not in already_counted[thisedge]["losses"]:
                        # add this loss to the set of already counted losses
                        already_counted[thisedge]["losses"].append(thisloss)
                        # add this loss to the fusion entries
                        thisentry = {"thisedge": thisedge,
                                     "thisloss": thisloss,
                                     "loss_size": thisG.nodes[thisloss]["size"],
                                     "loss_percent_of_total": node_size_fraction_of_total_size(thisG, thisloss),
                                     "loss_percent_of_largest": node_size_fraction_of_largest(thisG, thisloss),
                                     "loss_CC_size": node_size_CC(thisG, thisloss),
                                     "loss_CC_percent_of_total": node_size_CC_fraction_of_total_size(thisG, thisloss),
                                     "loss_CC_percent_of_largest": node_size_CC_fraction_of_largest(thisG, thisloss)
                                     }
                        dispersion_entries.append(thisentry)
                        #print("removing node {}".format(thisloss), " nodes are: {}".format(thisG.nodes))
                        thisG = delete_node_resize(thisG, thisloss)
            # now we check if there is a colocalization here
            if row_change["colocalizations"] != []:
                # iterate through each colocalization event
                for thiscoloc in row_change["colocalizations"]:
                    if thiscoloc not in already_counted[thisedge]["colocalizations"]:
                        already_counted[thisedge]["colocalizations"].append(thiscoloc)
                        # SPECIAL CASE FOR COLOCALIZATIONS
                        #  We must now determine if this fusion is already in the graph.
                        #  It may already be in the graph if these two ALGs are already colocalized.
                        if nodes_in_same_CC(thisG, thiscoloc):
                            # We have already counted this fusion already indirectly.
                            #  This is the transitive propery of colocalization.
                            #  If A has fused with B (AxB), and B has fused with C (BxC), then A has fused with C (AxC).
                            pass
                        else:
                            # This colocalization has not already been counted.
                            # Now we count it directly in the context of the existing colocalizations.
                            # add this entry to the colocalization entries
                            thisentry = {"thisedge": thisedge,
                                         "thiscoloc": thiscoloc,
                                         "coloc0_size": thisG.nodes[thiscoloc[0]]["size"],
                                         "coloc1_size": thisG.nodes[thiscoloc[1]]["size"],
                                         "coloc0_percent_of_total": node_size_fraction_of_total_size(thisG, thiscoloc[0]),
                                         "coloc1_percent_of_total": node_size_fraction_of_total_size(thisG, thiscoloc[1]),
                                         "coloc0_percent_of_largest": node_size_fraction_of_largest(thisG, thiscoloc[0]),
                                         "coloc1_percent_of_largest": node_size_fraction_of_largest(thisG, thiscoloc[1]),
                                         "coloc0_CC_size": node_size_CC(thisG, thiscoloc[0]),
                                         "coloc1_CC_size": node_size_CC(thisG, thiscoloc[1]),
                                         "coloc0_CC_percent_of_total": node_size_CC_fraction_of_total_size(thisG, thiscoloc[0]),
                                         "coloc1_CC_percent_of_total": node_size_CC_fraction_of_total_size(thisG, thiscoloc[1]),
                                         "coloc0_CC_percent_of_largest": node_size_CC_fraction_of_largest(thisG, thiscoloc[0]),
                                         "coloc1_CC_percent_of_largest": node_size_CC_fraction_of_largest(thisG, thiscoloc[1])
                                        }
                            colocalization_entries.append(thisentry)
                            # fuse these ALGs
                            thisG = colocalize_these_nodes(thisG, thiscoloc)
    # change the dispersion entries to a df
    dispersion_df  = pd.DataFrame(dispersion_entries)
    # change the colocalization entries to a df
    coloc_df = pd.DataFrame(colocalization_entries)

    return dispersion_df, coloc_df

    # save the dispersion_unique_df

    print("this is the dispersion df")
    print(dispersion_df)
    print("this is the colocalization df")
    print(coloc_df)
    sys.exit()

    # THIS SECTION CONTAINS CODE FOR PLOTTING THAT NEEDS TO BE PARTITIONED INTO FUNCTIONS
    ## first do the colocalizations:
    #x = list(coloc_df["coloc0_CC_size"]) + list(coloc_df["coloc1_CC_size"])
    #y = list(coloc_df["coloc1_CC_size"]) + list(coloc_df["coloc0_CC_size"])
    ## make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    ## get the ALG sizes from the rbhdf
    ## Create the scatterplot. These dots will be blue.
    ##  make the dots 0.1 transparency
    #plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    #plt.xlabel('ALG Size')
    #plt.ylabel('ALG Size')
    #plt.title('Component fusion sizes')
    ## Display the plot
    #plt.show()
    #plt.close()

    ## Plot as the percent of the size of the longest
    #x = list(coloc_df["coloc0_CC_percent_of_largest"]) + list(coloc_df["coloc1_CC_percent_of_largest"])
    #y = list(coloc_df["coloc1_CC_percent_of_largest"]) + list(coloc_df["coloc0_CC_percent_of_largest"])
    ## make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    ## get the ALG sizes from the rbhdf
    ## Create the scatterplot. These dots will be blue.
    ##  make the dots 0.1 transparency
    #plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    ##  make the dots 0.1 transparency
    #plt.xlabel('CC Size as a fraction')
    #plt.ylabel('CC Size as a fraction')
    #plt.title('CC_size as fraction of largest in genome')
    ## Display the plot
    #plt.show()
    #plt.close()

    ## Plot as the percent of the size of the longest
    #x = list(coloc_df["coloc0_CC_percent_of_total"]) + list(coloc_df["coloc1_CC_percent_of_total"])
    #y = list(coloc_df["coloc1_CC_percent_of_total"]) + list(coloc_df["coloc0_CC_percent_of_total"])
    ## make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    ## get the ALG sizes from the rbhdf
    ## Create the scatterplot. These dots will be blue.
    ##  make the dots 0.1 transparency
    #plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    #plt.xlabel('CC Size as a fraction')
    #plt.ylabel('CC Size as a fraction')
    #plt.title('CC_size as fraction of the total genome size')
    ## Display the plot
    #plt.show()
    #plt.close()

    ## Plot as the percent of the size of the longest
    #x = list(coloc_df["coloc0_size"]) + list(coloc_df["coloc1_size"])
    #y = list(coloc_df["coloc1_size"]) + list(coloc_df["coloc0_size"])
    ## make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    ## get the ALG sizes from the rbhdf
    ## Create the scatterplot. These dots will be blue.
    ##  make the dots 0.1 transparency
    #plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    #plt.xlabel('CC Size, absolute')
    #plt.ylabel('CC Size, absolute')
    #plt.title('CC_size absolute')
    ## Display the plot
    #plt.show()
    #plt.close()

    ## Plot as the percent of the size of the longest
    #x = list(coloc_df["coloc0_percent_of_total"]) + list(coloc_df["coloc1_percent_of_total"])
    #y = list(coloc_df["coloc1_percent_of_total"]) + list(coloc_df["coloc0_percent_of_total"])
    ## make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    ## get the ALG sizes from the rbhdf
    ## Create the scatterplot. These dots will be blue.
    ##  make the dots 0.1 transparency
    #plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    #plt.xlabel('CC Size, absolute')
    #plt.ylabel('CC Size, absolute')
    #plt.title('CC_size colocalizations as percent of total')
    ## Display the plot
    #plt.show()
    #plt.close()

    ## Plot as the percent of the size of the longest
    #x = list(coloc_df["coloc0_percent_of_largest"]) + list(coloc_df["coloc1_percent_of_largest"])
    #y = list(coloc_df["coloc1_percent_of_largest"]) + list(coloc_df["coloc0_percent_of_largest"])
    ## make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    ## get the ALG sizes from the rbhdf
    ## Create the scatterplot. These dots will be blue.
    ##  make the dots 0.1 transparency
    #plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    #plt.xlabel('CC Size, absolute')
    #plt.ylabel('CC Size, absolute')
    #plt.title('CC_size colocalizations as percent of largest')
    ## Display the plot
    #plt.show()
    #plt.close()

class coloc_array():
    """
    In this class we have a object array that contains the colocalization information.
    We will keep track of both the observed matricies and the expected matricies.
    The observed and expected matrices will have the same dimensions to later test the
      observed versus expected.

    Input:
      - abs_bin_size:  The size of the bins for the absolute size of the ALGs. For ref,
                        the size of the BCnS ALGs go from 12-207. Default value = 25
                        For the default bin size of 25, there will be a bin [0, 25), [25, 50), etc.
      - frac_bin_size: The size of the bins for the fraction of the largest ALG size.
                        For example, if the largest ALG is 200, and the ALG is 150, then
                        the fraction is 0.75. Default value = 0.05
                        The default of 0.05 means that there will be a bin:
                            - [0, 0.05), [0.05, 0.1) . . . [0.95, 1.0]
    """
    # initialize
    def __init__(self, abs_bin_size=25, frac_bin_size=0.05):
        # assert that the abs_bin_size is an integer
        if not isinstance(abs_bin_size, int):
            raise Exception("abs_bin_size must be an integer")
        # assert that the frac_bin_size is a float
        if not isinstance(frac_bin_size, float):
            raise Exception("frac_bin_size must be a float")
        # set the bin sizes
        self.abs_bin_size  = abs_bin_size
        self.frac_bin_size = frac_bin_size

        # The number of times we ran the forward loss/fusion analysis.
        self.num_observed_observations = 0
        # The number of times we ran the randomized ALG loss/fusion analysis.
        self.num_expected_observations = 0
        # initialize the observed and expected matrices. Just an empty dict for now
        self.observed_matrix  = {}
        self.expected_matrix  = {}

        # initialize the observed and expected matrices. Use a dict to store the observed and expected matrix
        ffrange = list(np.arange(0, 1+self.frac_bin_size, self.frac_bin_size))
        self.observed_matrix_frac  = {}
        self.expected_matrix_frac  = {}
        for i in range(len(ffrange)):
            for j in range(i, len(ffrange)):
                self.observed_matrix_frac[(ffrange[i], ffrange[j])] = 0
                self.expected_matrix_frac[(ffrange[i], ffrange[j])] = 0

        # Make structures to record the traces of the bins over time to check for convergence
        # For now store as a list of dicts. Convert to a dataframe on print
        self.ove_trace_size = []
        self.ove_trace_frac = []

    def _size_to_bin(self, size):
        """
        uses the bin size information to figure our in which bin this size belongs.
        We are using left exclusive right inclusive for the bin sizes
        (0, 25], (25, 50], etc.
        """
        return int(size/self.abs_bin_size) * self.abs_bin_size

    def _frac_to_bin(self, size):
        """
        uses the bin size information to figure our in which bin this size belongs.
        We are using left exclusive right inclusive for the bin sizes
        (0, 0.05], (0.05, 0.1], etc.
        """
        return int(size/self.frac_bin_size) * self.frac_bin_size

    def add_matrix(self, coloc_array, ooe):
        """
        This function takes in a coloc_df and adds it to the observed and expected matrices.
        """
        dicts    = {"frac": {}, "abs": {}}
        opp_dict = {"frac": {}, "abs": {}}
        if ooe == "observed":
            self.num_observed_observations += 1
            dicts    = {"frac": self.observed_matrix_frac, "abs": self.observed_matrix}
        elif ooe == "expected":
            self.num_expected_observations += 1
            dicts = {"frac": self.expected_matrix_frac, "abs": self.expected_matrix}
        else:
            print("ooe was {}".format(ooe))
            raise Exception("ooe must be either observed or expected")
        # iterate through the rows
        for i, row in coloc_array.iterrows():
            # get the size of the ALGs
            bin0 = self._size_to_bin(row["coloc0_size"])
            bin1 = self._size_to_bin(row["coloc1_size"])
            size_key = tuple(sorted([bin0,bin1]))
            for thisdict in [self.observed_matrix, self.expected_matrix]:
                if size_key not in thisdict:
                    thisdict[size_key] = 0

            dicts["abs"][size_key] += 1

            # get the fraction of the largest ALG
            frac0 = self._frac_to_bin(row["coloc0_percent_of_largest"])
            frac1 = self._frac_to_bin(row["coloc1_percent_of_largest"])
            frac_key = tuple(sorted([frac0, frac1]))
            dicts["frac"][frac_key] += 1

        # Don't generate traces for individual runs. Just do it for the whole thing
        #  at the end when we combine all the MC runs.
        #trace_granularity = 1
        ## now we go through the matrix and add the traces
        #if self.num_expected_observations%trace_granularity == 0:
        #    for thisbin in dicts["abs"]:
        #        self.ove_trace_size.append({"bin": thisbin,
        #                                    "ove": 0 if self.expected_matrix[thisbin] == 0 else self.observed_matrix[thisbin]/self.expected_matrix[thisbin],
        #                                    "obs_round": self.num_observed_observations,
        #                                    "exp_round": self.num_expected_observations})
        #    for thisbin in dicts["frac"]:
        #        self.ove_trace_size.append({"bin": thisbin,
        #                                    "ove": 0 if self.expected_matrix_frac[thisbin] == 0 else self.observed_matrix_frac[thisbin]/self.expected_matrix_frac[thisbin],
        #                                    "obs_round": self.num_observed_observations,
        #                                    "exp_round": self.num_expected_observations})

    def save_obs_expected_file(self, filename):
        """
        Saves the observed and expected matrices, abs and frac, to a single file.
        Later this can be read in and added to the datastrtuctures in this object.
        """
        # Saving to pwd doesn't break anything. Just check the dirname if it is nested inside pwd or elsewhere.
        if not os.path.dirname(filename) == "":
            # check that the directory in which filename is exists
            if not os.path.exists(os.path.dirname(filename)):
                raise Exception("Directory does not exist: {}".format(os.path.dirname(filename)))
        # check that the file does not already exist. We won't overwrite.
        if os.path.exists(filename):
            raise Exception("File already exists: {}".format(filename))
        # convert the traces to a dataframe
        entries = []
        for key in self.observed_matrix:
            entries.append({"bin":       key,
                            "ob_ex":     "observed",
                            "size_frac": "size",
                            "counts":    self.observed_matrix[key],
                            "obs_count": self.num_observed_observations})
        for key in self.expected_matrix:
            entries.append({"bin":       key,
                            "ob_ex":     "expected",
                            "size_frac": "size",
                            "counts":    self.expected_matrix[key],
                            "obs_count": self.num_expected_observations})
        for key in self.observed_matrix_frac:
            entries.append({"bin":       key,
                            "ob_ex":     "observed",
                            "size_frac": "frac",
                            "counts":    self.observed_matrix_frac[key],
                            "obs_count": self.num_observed_observations})
        for key in self.expected_matrix_frac:
            entries.append({"bin":       key,
                            "ob_ex":     "expected",
                            "size_frac": "frac",
                            "counts":    self.expected_matrix_frac[key],
                            "obs_count": self.num_expected_observations})
        # make a df of the entries
        df = pd.DataFrame(entries)
        print(df)
        # save it to a file, with headers, no indices
        df.to_csv(filename, sep="\t", index=False)

    def ret_obs_expected(self):
        """
        Perform the math on the whole matrices and return
        """
        ove_size = {}
        ove_frac = {}

        for thisbin in self.observed_matrix:
            ove_size[thisbin] = 0 if self.expected_matrix[thisbin] == 0 else self.observed_matrix[thisbin]/self.expected_matrix[thisbin]
        for thisbin in self.observed_matrix_frac:
            ove_frac[thisbin] = 0 if self.expected_matrix_frac[thisbin] == 0 else self.observed_matrix_frac[thisbin]/self.expected_matrix_frac[thisbin]

        return ove_size, ove_frac

def run_n_simulations_save_results(sampledfpath, algdfpath, filename,
                                   num_sims=10, abs_bin_size=25, frac_bin_size=0.05,
                                   verbose = False):
    """
    This function runs one Monte Carlo simulation unit, and should be used for parallelization.
    This function runs an even number of simulations of the loss/fusion analysis with real and randomized ALG sizes.
    Then, it saves the results to a .tsv file.
    These TSV files can be collated later with a coloc_array object to compile a larger dataset.
    """
    sampledf = pd.read_csv(sampledfpath, sep="\t")
    algdf = parse_rbh_file(algdfpath)

    counter = 0
    c = coloc_array(abs_bin_size=abs_bin_size, frac_bin_size=frac_bin_size)
    while counter < num_sims:
        if verbose:
            print("   - Running simulation {}".format(counter + 1), end="\r")
        dispersion_df, coloc_df  = stats_df_to_loss_fusion_dfs(sampledf, algdf, randomize_ALGs=False)
        c.add_matrix(  coloc_df, "observed")
        dispersion_df, coloc_df  = stats_df_to_loss_fusion_dfs(sampledf, algdf, randomize_ALGs=True)
        c.add_matrix(  coloc_df, "expected")
        counter +=1
    print("   - Running simulation {}".format(counter + 1))
    # save the results to a file
    c.save_obs_expected_file(filename)
    # safe return value
    return 0

def generate_stats_df(sample_df_filepath, outfilename) -> int:
    """
    Reads in the sample df, gets fusion stats along the tree, and saves the results as a tsv file.
    """
    # Check that the file exists
    if not os.path.exists(sample_df_filepath):
        raise Exception("File does not exist: {}".format(sample_df_filepath))

    # Read the TSV file into a pandas dataframe
    sampledf = pd.read_csv(sample_df_filepath, sep="\t")

    changedf = parse_gain_loss_from_perspchrom_df(sampledf)
    statsdf  = stats_on_changedf(sampledf, changedf)

    # check that the directory in which the outfilename exists
    # Saving to pwd doesn't break anything. Just check the dirname if it is nested inside pwd or elsewhere.
    if not os.path.dirname(outfilename) == "":
        # check that the directory in which filename is exists
        if not os.path.exists(os.path.dirname(outfilename)):
            raise Exception("Directory does not exist: {}".format(os.path.dirname(outfilename)))

    # overwriting is ok. Don't bother to check if the file exists
    # save the stats df to a tsv file with headers
    statsdf.to_csv(outfilename, sep="\t", index=False)
    return 0

def read_simulations_and_make_heatmaps(simulation_filepaths, outfilename):
    """
    This function reads in a file list of simulation files, sums up all the data,
     and makes a heatmap of the results. We are able to do this as it is a Monte Carlo
     simulation, so the results are independent.
    """
    all_dfs = []
    for thisfile in simulation_filepaths:
        # make sure the file exists
        if not os.path.exists(thisfile):
            raise Exception("File does not exist: {}".format(thisfile))
        # read in the file as a pandas df
        all_dfs.append(pd.read_csv(thisfile, sep="\t"))
    # concatenate the dfs
    df = pd.concat(all_dfs)
    # delete all the dfs after concat
    del all_dfs

    # open a pdf to save plots independently to it
    with PdfPages(outfilename) as pdf:
        # For size, first just get the rows that are size_frac.
        # After we get the df_size_obs and df_size_exp, we need to sum up all
        #  the rows that have the same value for bin. We can ignore obs_count
        #  because it simply tells us how many times the simulation was run in
        #  this instance.
        for size_frac in ["size", "frac"]:
            # get the rows that are size_frac
            df_size = df[df["size_frac"] == size_frac]
            df_size = df_size.sort_values(by=["bin"])
            # get the observed and expected dfs
            df_size_obs = df_size[df_size["ob_ex"] == "observed"]
            df_size_exp = df_size[df_size["ob_ex"] == "expected"]
            # sum up the counts for each bin. Add a pseudo-count of 1 to avoid divide by zero errors.
            df_size_obs = df_size_obs.groupby("bin")["counts"].sum().to_dict()
            df_size_exp = df_size_exp.groupby("bin")["counts"].sum().to_dict()
            df_size_obs = {thisbin: df_size_obs[thisbin] + 1 for thisbin in df_size_obs}
            df_size_exp = {thisbin: df_size_exp[thisbin] + 1 for thisbin in df_size_exp}

            # Get the obs/exp matrix. Use log2(obs/exp). We don't need to worry about divide
            #  by zero errors because we added a pseudo-count of 1 to each bin.
            ove_size = {}
            for thisbin in df_size_obs:
                ove_size[thisbin] = np.log2(df_size_obs[thisbin]/df_size_exp[thisbin])

            # infer the step size from the bin sizes in the "bin" column. Get it from the step size
            step_size = sorted(set([ast.literal_eval(x)[0] for x in df_size["bin"].unique()]))[1] - sorted(set([ast.literal_eval(x)[0] for x in df_size["bin"].unique()]))[0]
            # make a heatmap using the ove_size.
            # We want a custom heatmap where anything above 1 is red, and anything below 1 is blue. White is 1.
            #   The colormap should fade to white only exactly at 1. Anything above that quickly goes to red. Anything below that quickly goes to blue.
            #   The keys are the (x,y) for the heatmap and the values are the ove.
            # Use matplotlib patches to make the heatmap
            # get the x and y values
            x      = [ast.literal_eval(thisbin)[0] for thisbin in ove_size]
            y      = [ast.literal_eval(thisbin)[1] for thisbin in ove_size]
            values = [ove_size[thisbin] for thisbin in ove_size]
            colormap = plt.cm.RdBu_r
            # get the min and max values
            vmin = min(values)
            vmax = max(values)
            absmax = max([abs(vmin), abs(vmax)])
            center = 0  # Center value for the colormap
            norm = TwoSlopeNorm(vcenter=center, vmin=absmax*-1, vmax=absmax)

            # Make a plot with two subplots. One for the heatmap, and one for the colorbar.
            # The heatmap will be on the left and will be a square.
            # The colorbar will be on the right and will be thin vertical rectangle.
            # The heatmap will be 90% of the width of the figure.
            # The colorbar will be 10% of the width of the figure.
            # The colorbar will be centered vertically in the figure.
            # The heatmap will be centered vertically in the figure.
            # The heatmap will be centered horizontally in the left 90% of the figure.
            # The colorbar will be centered horizontally in the right 10% of the figure.
            # The heatmap will have a title that says "observed/expected"
            # The colorbar will have a title that says "log2(observed/expected)"

            # Make the figure
            fig = plt.figure(figsize=(6, 6))
            # Make the heatmap axes. The left 90% of the figure.
            ax = fig.add_axes([0.1, 0.1, 0.75, 0.75])
            # Make the colorbar axes. The right 10% of the figure.
            ax2 = fig.add_axes([0.9, 0.1, 0.075, 0.9])

            # use matplotlib patches to make squares for the heatmap.
            for i in range(len(x)):
                # get the color of the square. The colormap should center around 0
                color = colormap(norm(values[i]))
                # print the color as rgb. The value to 2 decimal places.
                #print(x[i], y[i], "{:.2f}".format(values[i]), rgb_float_to_hex(color))
                # make the square. The x and y are the bottom left corner of the square. The width and height are the step size.
                rect = mpatches.Rectangle((x[i],y[i]), step_size, step_size,
                                          linewidth=0, edgecolor='none',
                                          facecolor=color)
                ax.add_patch(rect)
            # set the x and y limits
            ax.set_xlim(min(x), max(x)+step_size)
            ax.set_ylim(min(y), max(y)+step_size)
            # set the x and y ticks
            ax.set_xticks(np.arange(min(x), max(x)+step_size, step_size))
            ax.set_yticks(np.arange(min(y), max(y)+step_size, step_size))
            # set the x and y tick labels.
            # If the values are floats, only keep to 2 decimal places. For both x and y.
            # For x, rotate the labels 90 degrees.
            if isinstance(min(x), float):
                ax.set_xticklabels(["{:.2f}".format(xtick) for xtick in np.arange(min(x), max(x)+step_size, step_size)], rotation=90)
                ax.set_yticklabels(["{:.2f}".format(ytick) for ytick in np.arange(min(y), max(y)+step_size, step_size)])
            else:
                ax.set_xticklabels(np.arange(min(x), max(x)+step_size, step_size), rotation=90)
                ax.set_yticklabels(np.arange(min(y), max(y)+step_size, step_size))
            # set the x and y labels
            ax.set_xlabel("Smaller ALG size")
            ax.set_ylabel("Larger ALG size")
            # set the title
            ax.set_title("ALG fusion size, {}".format(size_frac))
            # 500 units from absmin and absmax, make the colorbar
            thisrange = np.linspace(absmax * -1, absmax, 500)
            stepsize = thisrange[1] - thisrange[0]
            for i in thisrange:
                ax2.add_patch(mpatches.Rectangle((0, i), 1, stepsize, color=colormap(norm(i))))
            # set the limits of ax2
            ax2.set_ylim(absmax * -1, absmax)
            # set the title of ax2
            ax2.set_title("log2(observed/expected)")
            # turn off the x axis and turn on the y-axis. Just plot integers on y -axis
            ax2.set_xticks([])
            # get the smallest integer
            smallest_int = int(absmax * -1)
            # get the largest integer
            largest_int = int(absmax)
            # y labes on the right side
            ax2.yaxis.tick_right()
            # only plot integers on the y-axis
            ax2.set_yticks(     np.arange(smallest_int, largest_int+1, 1))
            # set the y tick labels
            ax2.set_yticklabels(np.arange(smallest_int, largest_int+1, 1))
            # set the y label
            ax2.set_ylabel("log2(observed/expected)")
            # set the x label
            ax2.set_xlabel("Colorbar")
            # save the figure to a pdf.
            pdf.savefig(bbox_inches="tight")
            plt.close()

def main():
    # Specify the file path of the TSV file. Use sys.argv[1]. The file will be called something like per_species_ALG_presence_fusions.tsv
    generate_stats_df(sys.argv[1], "statsdf.tsv")

    num_simulations = 1000
    sims_per_run    = 20
    for i in range(int(num_simulations/sims_per_run)):
        print("running simulation {}/{}".format(i+1, int(num_simulations/sims_per_run)))
        outname = "sim_results_{}_{}.tsv".format(i, sims_per_run)
        run_n_simulations_save_results(sys.argv[1],
                                       sys.argv[2],
                                       outname,
                                       num_sims=20,
                                       abs_bin_size=25,
                                       frac_bin_size=0.05,
                                       verbose = True)

    #simulation_filepaths = ["sim_results_20.tsv"]
    #read_simulations_and_make_heatmaps(simulation_filepaths, "simulations.pdf")


    ## find all the
    #ove_size, ove_frac = c.ret_obs_expected()
    #print(ove_size)
    ## make a heatmap using the ove_size.
    ## the keys are (x,y) tuples. The values are the ove. Use a red-white-blue color scale. Red means over 1, blue means under 1
    ## get the x and y values
    #x = [thisbin[0] for thisbin in ove_size]
    #y = [thisbin[1] for thisbin in ove_size]
    ## get the ove values
    #ove = [ove_size[thisbin] for thisbin in ove_size]
    ## make the heatmap
    #plt.scatter(x,y, c=ove, cmap="RdBu_r")
    #plt.xlabel('ALG Size')
    #plt.ylabel('ALG Size')
    #plt.title('Component fusion sizes')
    ## Display the plot
    #plt.show()


    #print()
    ## convert the ove_trace_frac to a df and print
    #print(c.ove_trace_frac)

def rgb_float_to_hex(list_of_rgb_floats):
    """
    Converts a list of rgb floats to a hex string.
    """
    return '#%02x%02x%02x' % (int(list_of_rgb_floats[0]*255), int(list_of_rgb_floats[1]*255), int(list_of_rgb_floats[2]*255))

if __name__ == "__main__":
    main()