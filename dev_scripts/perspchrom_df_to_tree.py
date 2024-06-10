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
import glob

# import the parse_rbh_file from the plot_ALG_fusions.py script
# this is a function that parses the RBH file into a dataframe
import rbh_tools

import matplotlib.pyplot as plt
# import patches
import matplotlib.patches as mpatches
# use twoslope norm to make a diverging color map
from matplotlib.colors import Normalize, TwoSlopeNorm

# Filter out the DeprecationWarning related to the py23 module
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning, module="fontTools.misc.py23")
# import pdfpages
from matplotlib.backends.backend_pdf import PdfPages

# use networkx to make graphs for the lineage-specific fusion/losses
import networkx as nx

# odp stuff to format the plot
# ODP-specific imports
thisfile_path = os.path.dirname(os.path.realpath(__file__))
scripts_path = os.path.join(thisfile_path, "../scripts")
sys.path.insert(1, scripts_path)
import odp_plotting_functions as odp_plot


def rgb_float_to_hex(list_of_rgb_floats):
    """
    Converts a list of rgb floats to a hex string.
    """
    return '#%02x%02x%02x' % (int(list_of_rgb_floats[0]*255), int(list_of_rgb_floats[1]*255), int(list_of_rgb_floats[2]*255))

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
           - source_taxid    - the older node in the phylogenetic tree.
           - target_taxid    - the younger node in the phylogenetic tree. Together the source and target form a branch.
           - colocalizations - ALGs that came together onto the same chromosome on this branch.
           - losses          - ALGs that were lost/dispersed/became untraceable on this branch.
    """
    splitstring = GL_string.split("-")
    # Go every two elements in the list.
    entries = []
    # These will be the taxids, and the things between them are the changes.
    for i in range(0,len(splitstring)-1, 2):
        taxid1 = int(splitstring[i])
        taxid2 = int(splitstring[i+2])
        colocloss = splitstring[i+1].lstrip("(").rstrip(")").split("]|[")
        # interpret the next two strings as lists
        colocs = eval(      colocloss[0] + "]")
        losses = eval("[" + colocloss[1] + "]")
        splits = eval("[" + colocloss[2]      )
        print("These are splits: {}".format(splits))
        entry = {"source_taxid": taxid1,    # already an int, we cast it earlier
                 "target_taxid": taxid2,    # already an int, we cast it earlier
                 "colocalizations": colocs, # this is a list
                 "losses": losses,          # this is a list
                 "splits": splits,          # this is a list
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
      - A pandas dataframe with stats on the changes that happened on each branch.
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
    # The change_frac_total is, for any given event, the percent of times we inferred that event on that branch.
    #  For example, the change ('A1a', 'A1b') was observed 1230 times on the branch leading to 33213,
    #                                                         5 times on the branch leading to 41711,
    #                                                         2 times on the branch leading to 42450,
    #                                                         2 times on the branch leading to 9871,
    #                                                         and a few times in other places.
    #  To calculate the change_frac_total, we sum the observations for this change type. 1230 + 5 + 2 + 2 = 1239.
    #  Then, these are the change_frac_total values for each of the branches:
    #                                                         target_taxid  change          counts  change_frac_total
    #                                                         33213         ('A1a', 'A1b')  1230    1230/1239 = 0.993
    #                                                         41711         ('A1a', 'A1b')  5       5/1239   = 0.004
    #                                                         42450         ('A1a', 'A1b')  2       2/1239   = 0.002
    #                                                          9871         ('A1a', 'A1b')  2       2/1239   = 0.002
    #  This calculation is therefore useful for QC'ing how well the algorithm is working on known fusion or loss events.
    #   For example, 33213 is the clade leading up to bilateria, and this is the branch on which we know those changes occurred.
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

def stats_df_to_loss_fusion_dfs(perspchromdf, ALGdf,
                                obs_seed = 0,
                                randomize_ALGs = False):
    """
    Use the perspchromdf, the statsdf, and an RBH file to make loss/fusion plots.
    The question is whether larger or smaller ALGs, or specific ALG combos,tend to
        fuse or be lost more often.

    The way that we do this is by going through each species, and analyzing the changes in
        ALGs on each lineage. If a change has been observed already on a lineage, then that
        count is not double counted.

    The final matrix will show on one axis, the smaller ALG contributing to a fusion, and on
        the other axis, the larger ALG contributing to a fusion. Doing this in a per-branch
        fashion will allow us to track the existing fusions to more accurately consider the
        changes at each evolutionary timepoint.

    To show what is the most prevalent pattern, we will use a Monte Carlo simulation to determine
      the observed/expected ratio given all of the possible trees.

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
        - thisedge:     The edge on which the fusion/loss happened.
        - thiscoloc:    The colocalization that happened on this edge.
        - coloc0_size:  The size of the first ALG in the colocalization.
        - coloc1_size:  The size of the second ALG in the colocalization.
        - coloc0_percent_of_total:   The size of the first ALG in the colocalization as a fraction of the total size of the graph.
        - coloc1_percent_of_total:   The size of the second ALG in the colocalization as a fraction of the total size of the graph.
        - coloc0_percent_of_largest: The size of the first ALG in the colocalization as a fraction of the largest ALG currently in the graph.
        - coloc1_percent_of_largest: The size of the second ALG in the colocalization as a fraction of the largest ALG currently in the graph.
        - coloc0_CC_size:            The size of all the ALGs in the first connected component.
        - coloc1_CC_size:            The size of all the ALGs in the second connected component.
        - coloc0_CC_percent_of_total:   The size of the first connected component as a fraction of the total size of the graph.
        - coloc1_CC_percent_of_total:   The size of the second connected component as a fraction of the total size of the graph.
        - coloc0_CC_percent_of_largest: The size of the first connected component as a fraction of the largest ALG currently in the graph.
        - coloc1_CC_percent_of_largest: The size of the second connected component as a fraction of the largest ALG currently in the graph.
    """
    ALG_info = ALGdf.copy()
    ALG_random_lookup = {}
    if randomize_ALGs:
        # randomize the ALG sizes
        ALG_info["Size"] = ALG_info["Size"].sample(frac=1).reset_index(drop=True)
        # add a random mapping column for the ALG names
        ALG_info["random"] = ALG_info["ALGname"].sample(frac=1).reset_index(drop=True)
        # This is only used in the random case
        ALG_random_lookup = ALG_info.set_index("ALGname")["random"].to_dict()

    ALG_info = ALG_info.sort_values(by="ALGname").reset_index(drop=True)

    # Make a starting graph where each "ALGname" in ALG_info is a node,
    #  each node has a size attribute (integer), and a color attribute (string)
    G = nx.Graph()
    # add the nodes
    for i, row in ALG_info.iterrows():
        G.add_node(row["ALGname"], size=int(row["Size"]), color=row["Color"])

    # Every time we run this program we will randomize the row order. Doing this will sample the different
    #  predicted events depending on the order in which we add fusions or losses to the graph.
    perspchromdf = perspchromdf.sample(frac=1, random_state=obs_seed)

    already_counted        = {}
    dispersion_entries     = []
    colocalization_entries = []
    ALG_coloc_entries      = []
    for i, row in perspchromdf.iterrows():
        # For this species, make a copy of the starting ALG graph. We will manipulate this graph. to track the changes.
        thisG = G.copy()
        # parse the gain/loss string
        changedf = parse_gain_loss_string(row["changestrings"], row["species"])

        for j, row_change in changedf.iterrows():
            thisedge = (row_change["source_taxid"], row_change["target_taxid"])
            # first check if there has been either a fusion or a loss here.
            # If there is a fusion or loss, but we haven't seen something here before, accounted for them yet, add an entry to the already_counted dict.
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
                                     "loss_percent_of_total":      node_size_fraction_of_total_size(thisG, thisloss),
                                     "loss_percent_of_largest":    node_size_fraction_of_largest(thisG, thisloss),
                                     "loss_CC_size":               node_size_CC(thisG, thisloss),
                                     "loss_CC_percent_of_total":   node_size_CC_fraction_of_total_size(thisG, thisloss),
                                     "loss_CC_percent_of_largest": node_size_CC_fraction_of_largest(thisG, thisloss)
                                     }
                        dispersion_entries.append(thisentry)
                        #print("removing node {}".format(thisloss), " nodes are: {}".format(thisG.nodes))
                    # This must be outside the if statement, because we want to remove the node even if it has already been counted.
                    # This is the only way in which we can accurately model the loss of ALGs.
                    thisG = delete_node_resize(thisG, thisloss)
            # now we check if there is a colocalization here
            if row_change["colocalizations"] != []:
                # iterate through each colocalization event
                for thiscoloc in row_change["colocalizations"]:
                    # SECTION FOR SIZE, we use the graph structure so we have special considerations.
                    #   if the colocalization has not already been counted here, we count it.
                    if thiscoloc not in already_counted[thisedge]["colocalizations"]:
                        already_counted[thisedge]["colocalizations"].append(thiscoloc)

                        # SECTION FOR LOOKING AT ALG COLOCALIZATIONS
                        # For the ALG_coloc_entries, we want to keep track of every single coloc, even if redundant through transitive property.
                        #   In the case of randomize_ALGs, we will map the tuples to the random ALG names.
                        if randomize_ALGs:
                            ALG1 = ALG_random_lookup[thiscoloc[0]]
                            ALG2 = ALG_random_lookup[thiscoloc[1]]
                            ALG_coloc_entries.append({"thisedge": thisedge, "thiscolor": (ALG1, ALG2)})
                        else:
                            ALG_coloc_entries.append({"thisedge": thisedge, "thiscolor": thiscoloc})

                        # SPECIAL CASE FOR COLOCALIZATIONS
                        #  We must now determine if this fusion is already in the graph.
                        #  It may already be in the graph if these two ALGs are already colocalized.
                        if nodes_in_same_CC(thisG, thiscoloc):
                            # We have already counted this fusion already indirectly.
                            #  This is the transitive property of colocalization.
                            #  If A has fused with B (AxB), and B has fused with C (BxC), then A has fused with C (AxC).
                            # We don't do anything now, because we only end up in this case if we have already counted this fusion.
                            #  and added it to the colocalization_entries.
                            pass
                        else:
                            # This colocalization has not already been counted.
                            # Now we count it directly in the context of the existing colocalizations.
                            # add this entry to the colocalization entries
                            thisentry = {"thisedge": thisedge,
                                         "thiscoloc": thiscoloc,
                                         "coloc0_size": thisG.nodes[thiscoloc[0]]["size"],
                                         "coloc1_size": thisG.nodes[thiscoloc[1]]["size"],
                                         "coloc0_percent_of_total":      node_size_fraction_of_total_size(thisG, thiscoloc[0]),
                                         "coloc1_percent_of_total":      node_size_fraction_of_total_size(thisG, thiscoloc[1]),
                                         "coloc0_percent_of_largest":     node_size_fraction_of_largest(   thisG, thiscoloc[0]),
                                         "coloc1_percent_of_largest":     node_size_fraction_of_largest(   thisG, thiscoloc[1]),
                                         "coloc0_CC_size":               node_size_CC(thisG, thiscoloc[0]),
                                         "coloc1_CC_size":               node_size_CC(thisG, thiscoloc[1]),
                                         "coloc0_CC_percent_of_total":    node_size_CC_fraction_of_total_size(thisG, thiscoloc[0]),
                                         "coloc1_CC_percent_of_total":    node_size_CC_fraction_of_total_size(thisG, thiscoloc[1]),
                                         "coloc0_CC_percent_of_largest": node_size_CC_fraction_of_largest(thisG, thiscoloc[0]),
                                         "coloc1_CC_percent_of_largest": node_size_CC_fraction_of_largest(thisG, thiscoloc[1])
                                        }
                            colocalization_entries.append(thisentry)
                            # fuse these ALGs
                            thisG = colocalize_these_nodes(thisG, thiscoloc)
                    # If it has been counted already, we don't do anything.
                    else:
                        pass
    # change the dispersion entries to a df
    dispersion_df  = pd.DataFrame(dispersion_entries)
    # change the colocalization entries to a df
    coloc_df = pd.DataFrame(colocalization_entries)
    # ALG_coloc_entries is a list of dicts. We want to make this into a dataframe.
    ALG_coloc_df = pd.DataFrame(ALG_coloc_entries)
    return dispersion_df, coloc_df, ALG_coloc_df

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

        # NUMBER OF OBSERVATIONS
        # The number of times we ran the forward loss/fusion analysis.
        self.num_observed_observations = 0
        # The number of times we ran the randomized ALG loss/fusion analysis.
        self.num_expected_observations = 0

        # initialize the observed and expected matrices. Just an empty dict for now
        # The structure of these dicts is: dict[(source, target)][(bin1, bin2)] = count
        self.observed_matrix_size_abs  = {}
        self.expected_matrix_size_abs  = {}
        self.observed_matrix_size_CC  = {}
        self.expected_matrix_size_CC  = {}

        # Initialize the observed and expected matrices for the fractionals.
        #  Same as the above, we don't do initialize anything because we don't know which phylogenetic branches
        #   we will have.
        self.observed_matrix_frac_abs  = {}
        self.expected_matrix_frac_abs  = {}
        self.observed_matrix_frac_CC  = {}
        self.expected_matrix_frac_CC  = {}

        # Initialize the observed and expected matrices for the ALGs
        self.observed_matrix_ALG  = {}
        self.expected_matrix_ALG  = {}
        # The number of times we ran the forward loss/fusion analysis for ALGs.
        self.num_observed_observations_ALGs = 0
        # The number of times we ran the randomized ALG loss/fusion analysis for ALGs.
        self.num_expected_observations_ALGs = 0

        self.plotmatrix_concatdf = None  # This is the concatenation of the observed and expected matrices from many simulations
        self.plotmatrix_sumdf    = None  # This is the sum of the observed and expected matrices

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

    def plotmatrix_listoffiles_to_plotmatrix(self, simulation_filepaths):
        """
        This function takes in a list of files that contain the observed and expected matrices.
        It summarizes all of the dataframes into the summed matrix and the plotdf.

        Updates self.plotmatrix_concatdf and self.plotmatrix_sumdf
        """
        all_dfs = []
        for thisfile in simulation_filepaths:
            # make sure the file exists
            if not os.path.exists(thisfile):
                raise Exception("File does not exist: {}".format(thisfile))
            # read in the file as a pandas df
            tempdf = pd.read_csv(thisfile, sep="\t")
            all_dfs.append(tempdf)

            # From the last df, get the unique obs_count for observed.
            # If the length of this is not 1, then something went wrong with the simulation
            observed_obs_count = tempdf[tempdf["ob_ex"] == "observed"]["obs_count"].unique()
            if len(observed_obs_count) != 1:
                raise Exception("The length of the observed_obs_count is not 1. Something went wrong with the simulation.")
            self.num_observed_observations += observed_obs_count[0]
            # Also check this same thing for the expected rows
            expected_obs_count = tempdf[tempdf["ob_ex"] == "expected"]["obs_count"].unique()
            if len(expected_obs_count) != 1:
                raise Exception("The length of the expected_obs_count is not 1. Something went wrong with the simulation.")
            self.num_expected_observations += expected_obs_count[0]

        # CONCATENATED DF FOR THE TRACES
        # Concatenate the dfs. This will have all the entries from all the dfs.
        #  Many rows may be duplicated. This is only if two simulation files had the same
        #  number of counts.
        self.plotmatrix_concatdf = pd.concat(all_dfs).reset_index(drop = True)

        # check that the len of df is not zero. If it is zero, we have not successully read in any of the files
        if len(self.plotmatrix_concatdf) == 0:
            raise Exception("No data was read in. Check that the simulation filepaths are correct.")

        # SUMDF FOR FINAL PLOTTING. Sum up both the obs_count and the counts
        sumdf  = self.plotmatrix_concatdf.groupby(["ALG_num", "branch", "bin", "ob_ex", "size_frac", "abs_CC"])["counts"].sum().reset_index()
        sumdf2 = self.plotmatrix_concatdf.groupby(["ALG_num", "branch", "bin", "ob_ex", "size_frac", "abs_CC"])["obs_count"].sum().reset_index()
        # merge these so that the sumdf has both the counts and the obs_count
        self.plotmatrix_sumdf = pd.merge(sumdf, sumdf2, on=["ALG_num", "branch", "bin", "ob_ex", "size_frac", "abs_CC"])
        self.plotmatrix_sumdf["count_per_sim"] = self.plotmatrix_sumdf["counts"]/self.plotmatrix_sumdf["obs_count"]

        # make sure that ALG_num.unique: ['num', 'ALG']
        if not sorted(list(self.plotmatrix_sumdf["ALG_num"].unique())) == sorted(list(['num', 'ALG'])):
            raise Exception("self.plotmatrix_sumdf['ALG_num'].unique() must be ['num', 'ALG']")
        # make sure that ob_ex.unique: ['observed' 'expected'
        if not sorted(list(self.plotmatrix_sumdf["ob_ex"].unique())) == sorted(list(['observed', 'expected'])):
            raise Exception("self.plotmatrix_sumdf['ob_ex'].unique() must be ['observed' 'expected']")
        # make sure that size_frac.unique: ['size' 'frac']
        if not sorted(list(self.plotmatrix_sumdf["size_frac"].unique())) == sorted(list(['size', 'frac'])):
            raise Exception("self.plotmatrix_sumdf['size_frac'].unique() must be ['size' 'frac']")
        # make sure that abs_CC.unique: ['CC' 'abs']
        if not sorted(list(self.plotmatrix_sumdf["abs_CC"].unique())) == sorted(list(['CC', 'abs'])):
            raise Exception("self.plotmatrix_sumdf['abs_CC'].unique() must be ['CC' 'abs']")

        # return a safe value
        return 0

    def add_matrix_ALGs(self, ALG_coloc_df, ooe):
        """
        adds the ALG_coloc_df to the observed or expected matrix.
        """
        # assert that ooe is either observed or expected
        if not ooe in ["observed", "expected"]:
            raise Exception("ooe must be either 'observed' or 'expected'")

        dict = self.observed_matrix_ALG if ooe == "observed" else self.expected_matrix_ALG
        opp_dict =  self.expected_matrix_ALG if ooe == "observed" else self.observed_matrix_ALG

        if ooe == "observed":
            self.num_observed_observations_ALGs += 1
        elif ooe == "expected":
            self.num_expected_observations_ALGs += 1

        # now we add the entries to the dict
        for i, row in ALG_coloc_df.iterrows():
            thisbin = tuple(sorted([row["thiscolor"][0], row["thiscolor"][1]]))
            oppbin = tuple(sorted([row["thiscolor"][1], row["thiscolor"][0]]))
            for thisdict in [dict, opp_dict]:
                for bin in [thisbin, oppbin]:
                    if row["thisedge"] not in thisdict:
                        thisdict[row["thisedge"]] = {}
                    if bin not in thisdict[row["thisedge"]]:
                        thisdict[row["thisedge"]][bin] = 0
            dict[row["thisedge"]][thisbin] += 1
            dict[row["thisedge"]][oppbin]  += 1

    def _safe_add_branch_bin_to_dicts(self, dicts, thisedge, bin):
        """
        Helper function for add_matrix. For any number of dicts, add the branch and bin to the dicts.
        The structure of all the dicts will be dict[branch][bin] = 0
        """
        for thisdict in dicts:
            # first check if the phylogenetic branch is in the dict
            if thisedge not in thisdict:
                thisdict[thisedge] = {}
            if bin not in thisdict[thisedge]:
                thisdict[thisedge][bin] = 0

    def add_matrix(self, coloc_array, ooe):
        """
        This function takes in a coloc_df and adds it to the observed and expected matrices.
        This is useful for adding the data from a single sample to the matrices for simulations.
        """
        dicts     = {}
        opp_dicts = {}
        if ooe == "observed":
            self.num_observed_observations += 1
            dicts     = {"frac_abs": self.observed_matrix_frac_abs,
                         "frac_CC":  self.observed_matrix_frac_CC,
                         "size_abs": self.observed_matrix_size_abs,
                         "size_CC":  self.observed_matrix_size_CC}
            opp_dicts = {"frac_abs": self.expected_matrix_frac_abs,
                         "frac_CC":  self.expected_matrix_frac_CC,
                         "size_abs": self.expected_matrix_size_abs,
                         "size_CC":  self.expected_matrix_size_CC}
        elif ooe == "expected":
            self.num_expected_observations += 1
            dicts     = {"frac_abs": self.expected_matrix_frac_abs,
                         "frac_CC":  self.expected_matrix_frac_CC,
                         "size_abs": self.expected_matrix_size_abs,
                         "size_CC":  self.expected_matrix_size_CC}
            opp_dicts = {"frac_abs": self.observed_matrix_frac_abs,
                         "frac_CC":  self.observed_matrix_frac_CC,
                         "size_abs": self.observed_matrix_size_abs,
                         "size_CC":  self.observed_matrix_size_CC}
        else:
            print("ooe was {}".format(ooe))
            raise Exception("ooe must be either observed or expected")
        # iterate through the rows
        for i, row in coloc_array.iterrows():
            # ****************** SIZE OF ALGs ***************************
            # Keep track of the sizes of the single ALGs, not the CCs
            bin0 = self._size_to_bin(row["coloc0_size"])
            bin1 = self._size_to_bin(row["coloc1_size"])
            size_key = tuple(sorted([bin0,bin1]))
            self._safe_add_branch_bin_to_dicts([dicts["size_abs"], opp_dicts["size_abs"]],
                                               row["thisedge"],
                                               size_key)
            # In the above for loop, we have added the phylogenetic branch to the dict if it was not already there.
            # Here, we actually add the value. The above code ensures we will protect against divide by zero errors later.
            dicts["size_abs"][row["thisedge"]][size_key] += 1

            # ********************** SIZE OF CCs ***************************
            # Keep track of the CC sizes of the ALG fusions, not the size of the individual ALGs
            bin0 = self._size_to_bin(row["coloc0_CC_size"])
            bin1 = self._size_to_bin(row["coloc1_CC_size"])
            size_key = tuple(sorted([bin0,bin1]))
            self._safe_add_branch_bin_to_dicts([dicts["size_CC"], opp_dicts["size_CC"]],
                                               row["thisedge"],
                                               size_key)
            dicts["size_CC"][row["thisedge"]][size_key] += 1

            # ********************** FRACTION OF ALGs ***************************
            # get the fraction of the single ALG size versus the largest ALG size.
            frac0 = self._frac_to_bin(row["coloc0_percent_of_largest"])
            frac1 = self._frac_to_bin(row["coloc1_percent_of_largest"])
            frac_key = tuple(sorted([frac0, frac1]))
            self._safe_add_branch_bin_to_dicts([dicts["frac_abs"], opp_dicts["frac_abs"]],
                                               row["thisedge"],
                                               frac_key)
            dicts["frac_abs"][row["thisedge"]][frac_key] += 1

            #  ********************** FRACTION OF CCs ***************************
            # get the fraction of the CC size versus the largest ALG size.
            frac0 = self._frac_to_bin(row["coloc0_CC_percent_of_largest"])
            frac1 = self._frac_to_bin(row["coloc1_CC_percent_of_largest"])
            frac_key = tuple(sorted([frac0, frac1]))
            self._safe_add_branch_bin_to_dicts([dicts["frac_CC"], opp_dicts["frac_CC"]],
                                               row["thisedge"],
                                               frac_key)
            dicts["frac_CC"][row["thisedge"]][frac_key] += 1

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
        iteration_dict = {"observed_size_abs": self.observed_matrix_size_abs, "expected_size_abs": self.expected_matrix_size_abs,
                          "observed_size_CC":  self.observed_matrix_size_CC,  "expected_size_CC":  self.expected_matrix_size_CC,
                          "observed_frac_abs": self.observed_matrix_frac_abs, "expected_frac_abs": self.expected_matrix_frac_abs,
                          "observed_frac_CC":  self.observed_matrix_frac_CC,  "expected_frac_CC":  self.expected_matrix_frac_CC}

        # enforce that the number of observed observations and the number of expected observations are the same
        # This means that the values will not be over- or under-estimated.
        if not self.num_observed_observations == self.num_expected_observations:
            raise Exception("self.num_observed_observations must equal self.num_expected_observations")

        for dictkey in iteration_dict:
            for branch in iteration_dict[dictkey]:
                for bin in iteration_dict[dictkey][branch]:
                    # use the name observed if the string observed is in the variable name of thisdict
                    obex     = "observed" if "observed" in dictkey else "expected"
                    sizefrac = "size"     if "size"     in dictkey else "frac"
                    abscc    = "abs"      if "abs"      in dictkey else "CC"
                    entries.append({"ALG_num": "num",
                                    "branch": branch,
                                    "bin": bin,
                                    "ob_ex": obex,
                                    "size_frac": sizefrac,
                                    "abs_CC": abscc,
                                    "counts": iteration_dict[dictkey][branch][bin],
                                    "obs_count": self.num_observed_observations})

        # Enforce that the number of observed observations and the number of expected observations are the same
        # This means that the values will not be over- or under-estimated.
        if not self.num_observed_observations_ALGs == self.num_expected_observations_ALGs:
            raise Exception("self.num_observed_observations_ALGs must equal self.num_expected_observations_ALGs")
        # now we add the ALG entries
        iteration_dict = {"observed_ALG": self.observed_matrix_ALG, "expected_ALG": self.expected_matrix_ALG}
        for dictkey in iteration_dict:
            for branch in iteration_dict[dictkey]:
                for bin in iteration_dict[dictkey][branch]:
                    # use the name observed if the string observed is in the variable name of thisdict
                    obex     = "observed" if "observed" in dictkey else "expected"
                    entries.append({"ALG_num": "ALG",
                                    "branch": branch,
                                    "bin": bin,
                                    "ob_ex": obex,
                                    "size_frac": "size",
                                    "abs_CC": "abs",
                                    "counts": iteration_dict[dictkey][branch][bin],
                                    "obs_count": self.num_observed_observations_ALGs})

        # make a df of the entries
        df = pd.DataFrame(entries)
        # save it to a file, with headers, no indices
        df.to_csv(filename, sep="\t", index=False)

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
    algdf    = rbh_tools.parse_ALG_rbh_to_colordf(algdfpath)

    print("This is the sampledf")
    print(sampledf)
    print("This is the algdf")
    print(algdf)

    counter = 0
    c = coloc_array(abs_bin_size=abs_bin_size, frac_bin_size=frac_bin_size)
    while counter < num_sims:
        if verbose:
            print("   - Running simulation {}".format(counter + 1), end="\r")
        # get a random number. must be between 0 and 2^32 - 1
        random_integer = random.randint(0,4294967295)
        dispersion_df, coloc_df, ALG_coloc_df  = stats_df_to_loss_fusion_dfs(sampledf, algdf,
                                    obs_seed = random_integer, randomize_ALGs=False)
        c.add_matrix(          coloc_df, "observed")
        c.add_matrix_ALGs( ALG_coloc_df, "observed")
        dispersion_df, coloc_df, ALG_coloc_df  = stats_df_to_loss_fusion_dfs(sampledf, algdf,
                                    obs_seed = random_integer, randomize_ALGs=True)
        c.add_matrix(          coloc_df, "expected")
        c.add_matrix_ALGs( ALG_coloc_df, "expected")
        counter +=1
    print("   - Running simulation {}".format(counter))
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

def i2f(i, dimension):
    """
    converts the index to the fraction of the dimension.
    """
    return i/dimension

def gen_square_ax(in_offset_left, in_offset_bottom,
                  fw, fh, inches):
    """
    Makes an axis to make a square positioned at the given offsets,
     with the given number of inches in both dimensions
    """
    return [in_offset_left/fw,
            in_offset_bottom/fh,
            inches/fw, inches/fh]

def gen_square_ax_and_colorbar(
        in_offset_left, in_offset_bottom,
        fw, fh, inches):
    """
    Makes an axis to make a square positioned at the given offsets,
     with the given number of inches in both dimensions
    Also, adds a colorbar to the right as a second axis.
    """
    plot_params = [in_offset_left/fw,   # left offset
                   in_offset_bottom/fh, # bottom offset
                   inches/fw,           # width
                   inches/fh]           # height
    cbar_params = [in_offset_left/fw + inches/fw + 0.025, # left offset
                   in_offset_bottom/fh,                  # bottom offset
                   inches*0.05/fw,                        # width
                   inches/fh]                            # height
    return plot_params, cbar_params

def generate_ALG_obs_exp_counts_panel(ax, cax, sumdf, ALGdf):
    """
    Makes the observed/expected counts panel for the ALGs.
    """
    # make sure that the algdf is sorted by ALGsize
    algdf = ALGdf.copy()
    algdf = algdf.sort_values(by="Size").reset_index(drop=True)
    algdf["index"] = algdf.index
    ALG_to_index = dict(zip(algdf["ALGname"], algdf["index"]))
    alglist = list(algdf["ALGname"])

    # get the rows that are ALG
    df_ALG = sumdf[sumdf["ALG_num"] == "ALG"].copy()
    # get the first value of the bin to get the ALG name. This is easier than looking for a tuple
    df_ALG["bin"]  = df_ALG["bin"].apply(ast.literal_eval)
    df_ALG["ALG1"] = df_ALG["bin"].apply(lambda x: x[0])
    df_ALG["ALG2"] = df_ALG["bin"].apply(lambda x: x[1])

    # this makes a dict of the observed and expected values
    ove_size = df_to_obs_exp_dict(df_ALG)

    absmax = 6
    colormap = plt.cm.RdBu_r
    # we hard-code the colors so that anything above these values don't weight more
    center = 0  # Center value for the colormap
    norm = TwoSlopeNorm(vcenter=center, vmin=absmax*-1, vmax=absmax)

    # use matplotlib patches to make squares for the heatmap.
    seen_already = set()
    for thiskey in ove_size:
        ALG1, ALG2 = thiskey
        ALG1_val = ALG_to_index[ALG1]
        ALG2_val = ALG_to_index[ALG2]
        sorted_bin_tuple = tuple(sorted([ALG1, ALG2]))
        if sorted_bin_tuple in seen_already:
            raise Exception("This bin has already been seen. This should not happen. The bin is: {}".format(sorted_bin_tuple))
        seen_already.add(sorted_bin_tuple)
        # plot now
        x = -1
        y = -1
        # we want the smaller value first
        if ALG1_val < ALG2_val:
            x = ALG1_val
            y = ALG2_val
        elif ALG1_val > ALG2_val:
            x = ALG2_val
            y = ALG1_val
        value = ove_size[thiskey]
        color = colormap(norm(value))
        rect = mpatches.Rectangle((x, y), 1, 1,
                                  linewidth=0,
                                  edgecolor='none',
                                  facecolor=color)
        ax.add_patch(rect)

    # set the x and y limits
    ax.set_xlim(0, len(alglist))
    ax.set_ylim(0, len(alglist))
    # set the x and y ticks
    ax.set_xticks(np.arange(0.5, len(alglist) + 0.5, 1))
    ax.set_yticks(np.arange(0.5, len(alglist) + 0.5, 1))
    ax.set_xticklabels(alglist, rotation=90)
    ax.set_yticklabels(alglist)

    # COLORBAR
    # 500 units from absmin and absmax, make the colorbar
    thisrange = np.linspace(absmax * -1, absmax, 500)
    stepsize = thisrange[1] - thisrange[0]
    # we hard-code the colors so that anything above these values don't weight more
    center = 0  # Center value for the colormap
    norm = TwoSlopeNorm(vcenter=center, vmin=absmax*-1, vmax=absmax)
    for i in thisrange:
        cax.add_patch(mpatches.Rectangle((0, i), 1, stepsize, color=colormap(norm(i))))
    # set the limits of ax2
    cax.set_ylim(absmax * -1, absmax)
    # turn off the x axis and turn on the y-axis. Just plot integers on y -axis
    cax.set_xticks([])
    # get the smallest integer
    smallest_int = int(absmax * -1)
    # get the largest integer
    largest_int = int(absmax)
    # y labes on the right side
    cax.yaxis.tick_right()
    # only plot integers on the y-axis
    cax.set_yticks(     np.arange(smallest_int, largest_int+1, 1))
    # set the y tick labels
    cax.set_yticklabels(np.arange(smallest_int, largest_int+1, 1))
    # set the y label
    cax.set_ylabel("log2(observed/expected)")
    # set the x label
    cax.set_xlabel("Colorbar")

    ## now we make the specific labels depending on whether we're looking at the absolute size or the fraction of the largest size
    #if abs_CC == "abs":
    #    # set the x and y labels
    #    ax.set_xlabel("Smaller ALG size when fused. Not CC size.")
    #    ax.set_ylabel("Larger ALG size when fused. Not CC size")
    #    # set the title
    #    ax.set_title("Log2(obs/exp) value of fusion events for dif. sizes of ALGs, not CCs, *{}*".format(size_frac.upper()))

    #elif abs_CC == "CC":
    #    # set the x and y labels
    #    ax.set_xlabel("Smaller CC size when fused. Not ALG size.")
    #    ax.set_ylabel("Larger CC size when fused. Not ALG size")
    #    # set the title
    #    ax.set_title("Log2(obs/exp) value of fusion events for dif. sizes of CCs, not ALGs, *{}*".format(size_frac.upper()))

    return ax, cax

def generate_ALG_mean_counts_panel(ax, cax, sumdf, ALGdf):
    """
    Makes a heatmap of the mean counts per cell for the ALGs.
    """
    # make sure that the algdf is sorted by ALGsize
    algdf = ALGdf.copy()
    algdf = algdf.sort_values(by="Size").reset_index(drop=True)
    algdf["index"] = algdf.index
    ALG_to_index = dict(zip(algdf["ALGname"], algdf["index"]))
    ALG_to_size  = dict(zip(algdf["ALGname"], algdf["Size"]))
    alglist = list(algdf["ALGname"])

    # get the rows that are ALG
    df_ALG = sumdf[sumdf["ALG_num"] == "ALG"].copy()
    # get only the rows that are observed
    df_ALG = df_ALG[df_ALG["ob_ex"] == "observed"]
    # get the first value of the bin to get the ALG name. This is easier than looking for a tuple
    df_ALG["bin"]  = df_ALG["bin"].apply(ast.literal_eval)
    df_ALG["ALG1"] = df_ALG["bin"].apply(lambda x: x[0])
    df_ALG["ALG2"] = df_ALG["bin"].apply(lambda x: x[1])

    # This was for debugging - can be deleted
    # # print all the rows with the value R in col ALG1 or ALG2
    # print(df_ALG[(df_ALG["ALG1"] == "R") | (df_ALG["ALG2"] == "R")])

    # make sure that the length of the df is not zero
    if len(df_ALG) == 0:
        raise Exception("len(df_ALG) == 0. This means that, for some reason, the ALG colocalizations were not saved to the file.")

    # We want a custom heatmap where we scale from white (0) to red (the max)
    colormap = plt.cm.Reds
    norm = Normalize(vmin=0, vmax=int(max(df_ALG["count_per_sim"])) + 1)

    seen_already = set()
    for i, row in df_ALG.iterrows():
        sorted_bin_tuple = tuple(sorted(row["bin"]))
        if sorted_bin_tuple in seen_already:
            raise Exception("This bin has already been seen. This should not happen. The bin is: {}".format(sorted_bin_tuple))
        seen_already.add(sorted_bin_tuple)
        # plot now
        x = -1
        y = -1
        # we want the smaller value first
        if ALG_to_size[row["ALG1"]] < ALG_to_size[row["ALG2"]]:
            x = ALG_to_index[row["ALG1"]]
            y = ALG_to_index[row["ALG2"]]
        elif ALG_to_size[row["ALG1"]] > ALG_to_size[row["ALG2"]]:
            x = ALG_to_index[row["ALG2"]]
            y = ALG_to_index[row["ALG1"]]
        value = row["count_per_sim"]
        color = colormap(norm(value))
        rect = mpatches.Rectangle((x, y), 1, 1,
                                  linewidth=0,
                                  edgecolor='none',
                                  facecolor=color)
        ax.add_patch(rect)

    # set the x and y limits
    ax.set_xlim(0, len(alglist))
    ax.set_ylim(0, len(alglist))
    # set the x and y ticks
    ax.set_xticks(np.arange(0.5, len(alglist)+0.5, 1))
    ax.set_yticks(np.arange(0.5, len(alglist)+0.5, 1))
    # set the x and y tick labels.
    # For x, rotate the labels 90 degrees.
    ax.set_xticklabels(algdf["ALGname"], rotation=90)
    ax.set_yticklabels(algdf["ALGname"])

    # COLORBAR SECTION
    # 500 units from absmin and absmax, make the colorbar
    thisrange = np.linspace(0, int(max(df_ALG["count_per_sim"])) + 1, 500)
    stepsize = thisrange[1] - thisrange[0]
    for i in thisrange:
        cax.add_patch(mpatches.Rectangle((0, i), 1, stepsize, color=colormap(norm(i))))
    # set the limits of ax2
    cax.set_ylim(0, int(max(df_ALG["count_per_sim"])) + 1)
    # set the title of ax2
    # turn off the x axis and turn on the y-axis. Just plot integers on y -axis
    cax.set_xticks([])
    # get the smallest integer
    smallest_int = 0
    # get the largest integer
    largest_int = int(max(df_ALG["count_per_sim"])) + 1
    # y labes on the right side
    cax.yaxis.tick_right()
    ## only plot integers on the y-axis
    #cax.set_yticks(     np.arange(smallest_int, largest_int+1, 1))
    ## set the y tick labels
    #cax.set_yticklabels(np.arange(smallest_int, largest_int+1, 1))
    # set the x label
    cax.set_xlabel("Colorbar")

    # set the x and y labels
    ax.set_xlabel("Smaller ALG participating in fusion")
    ax.set_ylabel("Larger ALG participating in fusion")
    # set the title
    ax.set_title( "Mean count of fusion events for different ALG pairs")
    # set the y label
    cax.set_ylabel("Mean number of ALG fusions per topology")
    return ax, cax

def generate_mean_counts_panel(ax, cax, sumdf, size_frac, abs_CC):
    """
    Takes in the sumdf and makes a heatmap of the mean counts per cell.

    Input:
      - ax is the axis on which the heatmap will be plotted.
      - cax is the axis on which the colorbar will be plotted.
      - sumdf is the dataframe that contains all of the simulations added.
      - size_frac is either "frac" or "size", depending on what we want to plot.
      - abs_CC is either "abs" or "CC", depending on what we want to plot.
    """
    # we must assert that size_frac is either "frac" or "size"
    if size_frac not in ["frac", "size"]:
        raise Exception("size_frac must be either frac or size")
    # we must assert that abs_CC is either "abs" or "CC"
    if abs_CC not in ["abs", "CC"]:
        raise Exception("abs_CC must be either abs or CC")

    # get the rows that are size_frac
    df_size = sumdf[sumdf["size_frac"] == size_frac]
    # get the rows that match the abs_CC state
    df_size = df_size[df_size["abs_CC"] == abs_CC]
    # get the rows that are observed
    df_size = df_size[df_size["ob_ex"] == "observed"]
    # remove the rows that are ALG_num ALG
    df_size = df_size[df_size["ALG_num"] == "num"]

    x = []
    y = []
    values = []
    for i, row in df_size.iterrows():
        # get the bin
        thisbin = ast.literal_eval(row["bin"])
        # get the x and y
        x.append(thisbin[0])
        y.append(thisbin[1])
        # get the value
        values.append(row["count_per_sim"])

    # infer the step size from the bin sizes in the "bin" column. Get it from the step size
    step_size = sorted(set([ast.literal_eval(x)[0] for x in df_size["bin"].unique()]))[1] - sorted(set([ast.literal_eval(x)[0] for x in df_size["bin"].unique()]))[0]
    # make a heatmap using the count_per_sim
    # We want a custom heatmap where we scale from white (0) to red (the max)
    colormap = plt.cm.Reds
    norm = Normalize(vmin=0, vmax=int(max(values)) + 1)

    # use matplotlib patches to make squares for the heatmap.
    for i in range(len(x)):
        # only plot if the value is greater than 0
        if values[i] > 0:
            # get the color of the square. The colormap should center around 0
            color = colormap(norm(values[i]))
            # print the color as rgb. The value to 2 decimal places.
            # make the square. The x and y are the bottom left corner of the square. The width and height are the step size.
            rect = mpatches.Rectangle((x[i],y[i]), step_size, step_size,
                                      linewidth=0, edgecolor='none',
                                      facecolor=color)
            ax.add_patch(rect)
    # set the x and y limits
    ax.set_xlim(0, max(max(x), max(y))+step_size)
    ax.set_ylim(0, max(max(x), max(y))+step_size)
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

    # COLORBAR SECTION
    # 500 units from absmin and absmax, make the colorbar
    thisrange = np.linspace(0, int(max(values)) + 1, 500)
    stepsize = thisrange[1] - thisrange[0]
    for i in thisrange:
        cax.add_patch(mpatches.Rectangle((0, i), 1, stepsize, color=colormap(norm(i))))
    # set the limits of ax2
    cax.set_ylim(0, int(max(values)) + 1)
    # set the title of ax2
    # turn off the x axis and turn on the y-axis. Just plot integers on y -axis
    cax.set_xticks([])
    # get the smallest integer
    smallest_int = 0
    # get the largest integer
    largest_int = int(max(values)) + 1
    # y labes on the right side
    cax.yaxis.tick_right()
    ## only plot integers on the y-axis
    #cax.set_yticks(     np.arange(smallest_int, largest_int+1, 1))
    ## set the y tick labels
    #cax.set_yticklabels(np.arange(smallest_int, largest_int+1, 1))
    # set the x label
    cax.set_xlabel("Colorbar")

    # now we make the specific labels depending on whether we're looking at the absolute size or the fraction of the largest size
    if abs_CC == "abs":
        # set the x and y labels
        ax.set_xlabel("Smaller ALG size when fused. Not CC size.")
        ax.set_ylabel("Larger ALG size when fused. Not CC size")
        # set the title
        ax.set_title("Mean count of fusion events for different sizes of ALGs, not CC, *{}*".format(size_frac.upper()))
        # set the y label
        cax.set_ylabel("Mean number of ALG fusions (not CCs) per topology")

    elif abs_CC == "CC":
        # set the x and y labels
        ax.set_xlabel("Smaller CC size when fused. Not ALG size.")
        ax.set_ylabel("Larger CC size when fused. Not ALG size")
        # set the title
        ax.set_title("Mean count of fusion events for different sizes of CCs, not ALGs, *{}*".format(size_frac.upper()))
        # set the y label
        cax.set_ylabel("Mean number of CC fusions (not ALGs) per topology")
    return ax, cax

def generate_obs_exp_panel(ax, cax, sumdf, size_frac, abs_CC, absmax = 6):
    """
    This generates the observed/expexted panel and puts it in an existing axis.
    This uses a dataframe that has been summed up from many files

    Inputs:
      - ax:    The axis in which to put the heatmap.
      - sumdf: A dataframe that has been summed up from many files.
      - size_frac: Either "frac" or "size". This determines which heatmap to plot.
      - abs_CC: Either "abs" or "CC". This determines which heatmap to plot.
    """
    # we must assert that size_frac is either "frac" or "size"
    if size_frac not in ["frac", "size"]:
        raise Exception("size_frac must be either frac or size")
    # we must assert that abs_CC is either "abs" or "CC"
    if abs_CC not in ["abs", "CC"]:
        raise Exception("abs_CC must be either abs or CC")

    # get the rows that are size_frac
    df_size = sumdf[sumdf["size_frac"] == size_frac]
    # get the rows that match the abs_CC state
    df_size = df_size[df_size["abs_CC"] == abs_CC]
    # remove the rows in ALG_num that are ALG
    df_size = df_size[df_size["ALG_num"] == "num"]
    # this makes a dict of the observed and expected values
    ove_size = df_to_obs_exp_dict(df_size)

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
    # we hard-code the colors so that anything above these values don't weight more
    center = 0  # Center value for the colormap
    norm = TwoSlopeNorm(vcenter=center, vmin=absmax*-1, vmax=absmax)

    # use matplotlib patches to make squares for the heatmap.
    for i in range(len(x)):
        # get the color of the square. The colormap should center around 0
        color = colormap(norm(values[i]))
        # print the color as rgb. The value to 2 decimal places.
        # make the square. The x and y are the bottom left corner of the square. The width and height are the step size.
        rect = mpatches.Rectangle((x[i],y[i]), step_size, step_size,
                                  linewidth=0, edgecolor='none',
                                  facecolor=color)
        ax.add_patch(rect)
    # set the x and y limits
    ax.set_xlim(0, max(max(x), max(y))+step_size)
    ax.set_ylim(0, max(max(x), max(y))+step_size)
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

    # 500 units from absmin and absmax, make the colorbar
    thisrange = np.linspace(absmax * -1, absmax, 500)
    stepsize = thisrange[1] - thisrange[0]
    # we hard-code the colors so that anything above these values don't weight more
    center = 0  # Center value for the colormap
    norm = TwoSlopeNorm(vcenter=center, vmin=absmax*-1, vmax=absmax)
    for i in thisrange:
        cax.add_patch(mpatches.Rectangle((0, i), 1, stepsize, color=colormap(norm(i))))
    # set the limits of ax2
    cax.set_ylim(absmax * -1, absmax)
    # turn off the x axis and turn on the y-axis. Just plot integers on y -axis
    cax.set_xticks([])
    # get the smallest integer
    smallest_int = int(absmax * -1)
    # get the largest integer
    largest_int = int(absmax)
    # y labes on the right side
    cax.yaxis.tick_right()
    # only plot integers on the y-axis
    cax.set_yticks(     np.arange(smallest_int, largest_int+1, 1))
    # set the y tick labels
    cax.set_yticklabels(np.arange(smallest_int, largest_int+1, 1))
    # set the y label
    cax.set_ylabel("log2(observed/expected)")
    # set the x label
    cax.set_xlabel("Colorbar")

    # now we make the specific labels depending on whether we're looking at the absolute size or the fraction of the largest size
    if abs_CC == "abs":
        # set the x and y labels
        ax.set_xlabel("Smaller ALG size when fused. Not CC size.")
        ax.set_ylabel("Larger ALG size when fused. Not CC size")
        # set the title
        ax.set_title("Log2(obs/exp) value of fusion events for dif. sizes of ALGs, not CCs, *{}*".format(size_frac.upper()))

    elif abs_CC == "CC":
        # set the x and y labels
        ax.set_xlabel("Smaller CC size when fused. Not ALG size.")
        ax.set_ylabel("Larger CC size when fused. Not ALG size")
        # set the title
        ax.set_title("Log2(obs/exp) value of fusion events for dif. sizes of CCs, not ALGs, *{}*".format(size_frac.upper()))

    return ax, cax

def generate_trace_panel(ax, branches_in_clade, simulation_filepaths, ALG_num, frac_or_size, abs_CC):
    """
    This generates the traces for the observed/expected matrices.
    We will use this to check for convergence.
    Instead of plotting some normalized values, we will just plot the raw values.
    """
    # make sure that ALG_num is either ALG or num
    if not ALG_num in ["ALG", "num"]:
        raise Exception("ALG_num must be either ALG or num")
    # make sure that frac or size is either frac or size
    if not frac_or_size in ["frac", "size"]:
        raise Exception("frac_or_size must be either frac or size")
    # we must assert that abs_CC is either "abs" or "CC"
    if abs_CC not in ["abs", "CC"]:
        raise Exception("abs_CC must be either abs or CC")

    lines = {}
    # using dictionaries to keep track of these values, because pandas was failing sometimes
    dicts = {"observed": {},
             "expected": {}}
    num_sims = 0
    # this will be updated each time we go through
    for thisfile in simulation_filepaths:
        # if the df is None, then read in the file
        tempdf = pd.read_csv(thisfile, sep="\t")
        # filter on size_frac column
        tempdf = tempdf[tempdf["ALG_num"]   == ALG_num      ] # get the rows that are or are not ALGs
        tempdf = tempdf[tempdf["size_frac"] == frac_or_size ] # get the rows that match the method variable
        tempdf = tempdf[tempdf["abs_CC"]    == abs_CC       ] # get the rows that match the abs_CC variable
        # get only the rows that have the branches_in_clade
        tempdf = tempdf[tempdf["branch"].isin(branches_in_clade)]

        # assert that we are only dealing with a single type of size_frac and abs_CC
        assert_single_size_abs_condition(tempdf)
        num_sims += int(tempdf["obs_count"].max())
        # go through the rows and add values. Yes, using a for loop :)
        for i, row in tempdf.iterrows():
            # add the value to the dict
            if row["bin"] not in dicts[row["ob_ex"]]:
                dicts[row["ob_ex"]][   row["bin"]  ] = 0
            dicts[row["ob_ex"]][       row["bin"]  ] += row["counts"]

       # now that we have updated the database (this dicts object), we generate an ove dict and add those values to a
        ove_size = df_to_obs_exp_dict(dicts)
        # now we go through the matrix and add the traces for this value
        for thisbin in ove_size:
            if thisbin not in lines:
                lines[thisbin] = {"num_sim": [], "value": []}
            lines[thisbin]["num_sim"].append( num_sims         )
            lines[thisbin]["value"].append(   ove_size[thisbin])

    # now we plot each of the lines
    for thisbin in lines:
        #index_to_print = 2400
        ## get the index of num_sim that equals 2400
        ## TODO there is a bug with samples (100,250), (175, 275), and (150, 275)
        ##      For the bug, their bvalues go to zero in the obs/expected value after a certain number of simulations
        ##      that changes every time I run the program. I need to figure out why this is happening.
        #if index_to_print in lines[thisbin]["num_sim"]:
        #    # if the last 10 values are 0, then we can print the bin number
        #    if all([x == 0 for x in lines[thisbin]["value"][-10:]]):
        #        # get the index of the value that is closest to 2400
        #        thisix = lines[thisbin]["num_sim"].index(index_to_print)
        #        # plot the bin number above the line for their value around 2500
        #        ax.text(index_to_print, lines[thisbin]["value"][thisix], thisbin, fontsize=5)
        ax.plot(lines[thisbin]["num_sim"], lines[thisbin]["value"], color="black", alpha=0.1, lw=0.5)
    # change the xaxis to go from 0 to the num_sims
    ax.set_xlim(0, num_sims)
    # set the title
    return ax

def assert_single_size_abs_condition(df):
    """
    We have to do a lot of checks to make sure that the df is in the correct format.
     The checks make sure that we will not accidentally count more than one field at a time.

    Returns True if the df is in the correct format.
    Returns False if the df is not in the correct format, however an error will come up first
    """
    # We have to do a lot of checks to make sure that the df is in the correct format.
    #  The checks make sure that we will not accidentally count more than one field at a time.
    # We must assert that the ob_ex column is either observed or expected.
    target = list(sorted(["observed", "expected"]))
    for entry in list(sorted(df["ob_ex"].unique())):
        if not entry in target:
            raise Exception("df['ob_ex'].unique() must have contents of either observed or expected. Found {}".format(df["ob_ex"].unique()))

    # we must only find one value for size_frac. It must be either frac or size.
    if not len(df["size_frac"].unique()) == 1:
        raise Exception("df['size_frac'].unique() must be of length 1. Should be either frac or size. Found {}".format(df["size_frac"].unique()) )
    if not df["size_frac"].unique()[0] in ["frac", "size"]:
        raise Exception("df['size_frac'].unique() must be either frac or size")

    # we must only find one value for abs_CC. It must be either abs or CC.
    if not len(df["abs_CC"].unique()) == 1:
        raise Exception("df['abs_CC'].unique() must be of length 1. Should be either abs or CC")
    if not df["abs_CC"].unique()[0] in ["abs", "CC"]:
        raise Exception("df['abs_CC'].unique() must be either abs or CC")
    return True

def df_to_obs_exp_dict(df):
    """
    Takes in a df and returns the observed and expected dicts.
      - the df can be a pandas dataframe or a special dictionary structure.
      - the dictionary structure looks like this:

    dicts = {"observed": {bin: count, bin2: count},
             "expected": {bin: count, bin2: count} }
    """
    # This is the dict that we will return
    ove_size = {}

    # If the type is a pandas df.
    # This is deprecated since we are likely not using pandas dfs anymore for this function.
    if isinstance(df, pd.DataFrame):
        # we must check that the df is in the correct format. Only single fields for size_frac, and abs_CC
        assert_single_size_abs_condition(df)
        # We are done with the checks. Now we can get the observed and expected dicts.
        # get the observed and expected dfs
        df_size_obs = df[df["ob_ex"] == "observed"]
        df_size_exp = df[df["ob_ex"] == "expected"]
        # sum up the counts for each bin. Add a pseudo-count of 1 to avoid divide by zero errors.
        df_size_obs = df_size_obs.groupby("bin")["counts"].sum().to_dict()
        df_size_exp = df_size_exp.groupby("bin")["counts"].sum().to_dict()
        df_size_obs = {thisbin: df_size_obs[thisbin] + 1 for thisbin in df_size_obs}
        df_size_exp = {thisbin: df_size_exp[thisbin] + 1 for thisbin in df_size_exp}

        # Get the obs/exp matrix. Use log2(obs/exp). We don't need to worry about divide
        #  by zero errors because we added a pseudo-count of 1 to each bin.
        for thisbin in df_size_obs:
            ove_size[thisbin] = np.log2(df_size_obs[thisbin]/df_size_exp[thisbin])
    elif isinstance(df, dict) and ("observed" in df) and ("expected" in df):
        # allbins is a set of all of the bins in the observed and expected dicts
        allbins = set(list(df["observed"].keys()) + list(df["expected"].keys()))
        for thisbin in allbins:
            # add pseudocounts to avoid any divide by zero errors or log errors
            obsval = 1 if thisbin not in df["observed"] else df["observed"][thisbin] + 1
            expval = 1 if thisbin not in df["expected"] else df["expected"][thisbin] + 1
            obsexp = np.log2(obsval/expval)
            if obsexp < -10:
                obsexp = -10
            if obsexp > 10:
                obsexp = 10
            ove_size[thisbin] = obsexp
    return ove_size

class PhyloTree:
    """
    This class is used to store a phylogenetic tree.
    The tree is implemented as a directional graph.
    The tree can be constructed by lists of edges. The nodes are inferred from the edges.
    """
    def __init__(self) -> None:
        # initialize the graph using networkx
        self.G = nx.DiGraph()

    def add_taxname_to_all_nodes(self):
        """
        This function adds the taxname to a single node. Uses ete3.
        """
        # use ete3 to get the names of the taxids
        ncbi = NCBITaxa()

        for node in self.G.nodes():
            self.G.nodes[node]["taxname"] = ncbi.get_taxid_translator([node])[node].replace(" ", "-")

    def add_lineage_string(self, lineage_string) -> int:
        """
        The lineage strings look like this:
          - 1;131567;2759;33154;33208;6040;6042;1779146;1779162;6060;6061;1937956;6063

        Notes:
          - The edges from this string will be (1, 131567), (131567, 2759), (2759, 33154), etc.
        """
        fields = [int(x) for x in lineage_string.split(";")]
        for i in range(len(fields)-1):
            self.G.add_edge(fields[i], fields[i+1])
        return 0

    def build_tree_from_per_sp_chrom_df(self, per_sp_chrom_df) -> int:
        """
        This function takes in a per_sp_chrom_df and builds a tree from it.
        """
        # if this is a string look for a file
        if isinstance(per_sp_chrom_df, str):
            per_sp_chrom_df = pd.read_csv(per_sp_chrom_df, sep="\t")
        elif not isinstance(per_sp_chrom_df, pd.DataFrame):
            raise Exception("per_sp_chrom_df must be either a string or a pandas dataframe")
        # get the lineage strings
        lineage_strings = per_sp_chrom_df["taxidstring"].unique()
        # add each lineage string to the tree
        for thislineage in lineage_strings:
            self.add_lineage_string(thislineage)
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

def generate_node_taxid_file_from_per_sp_chrom_df(per_sp_chrom_df, outfile):
    """
    Builds a tree from the per_sp_chrom_df and writes a file that contains the taxid and the node ID
    """
    T = PhyloTree()
    T.build_tree_from_per_sp_chrom_df(per_sp_chrom_df)
    T.add_taxname_to_all_nodes()
    with open(outfile, "w") as f:
        for node in T.G.nodes():
            f.write("{}\t{}\n".format(node, T.G.nodes[node]["taxname"]))

def _make_one_simulation_plot(algdf, c, T, taxid, simulation_filepaths,
                              outfileprefix, absmax = 6):
    """
    makes a single plot for a single NCBI taxid

    Inputs:
      - c is the coloc_array object
      - T is the PhyloTree object
      - taxid is the NCBI taxid that we are plotting
    """
    plotting = True
    # At this point we're not even sure if this species is in the graph. We must first check.
    if not taxid in T.G.nodes():
        # it isn't in the graph, we can't plot it.
        plotting = False
        # use ete3 ncbi to get the taxon name
        ncbi = NCBITaxa()
        taxon_name = ncbi.get_taxid_translator([taxid])[taxid].replace(" ", "-")
    else:
        # There is a chance we can plot it still
        # First we figure out which branches are in the clade we care about.
        # we cast the tuples to strings because of how they're stored in the df
        branches_in_clade = [str(x) for x in T.get_edges_in_clade(taxid)]
        taxon_name = T.G.nodes[taxid]["taxname"]
        ## this is useful for debugging
        #print("We're plotting the taxid {} ({})".format(taxid, taxon_name))
        #print("branches in this clade are")
        #for thisbranch in branches_in_clade:
        #    taxid1 = int(thisbranch.split(",")[0].replace("(", ""))
        #    taxid2 = int(thisbranch.split(",")[1].replace(")", ""))
        #    print("  {}:  {}-{}".format(thisbranch, T.G.nodes[taxid1]["taxname"], T.G.nodes[taxid2]["taxname"]))
        # raise an error if ther eis nothing in this clade
        if len(branches_in_clade) == 0:
            raise Exception("There are no branches in the clade {}".format(taxid))
        # make a filtered df that only contains the branches in the clade
        plotdf = c.plotmatrix_sumdf[c.plotmatrix_sumdf["branch"].isin(branches_in_clade)]

        # now we don't care about the branch. Groupby the other columns and sum them up. They are all counts
        plotdf  = plotdf.groupby(["ALG_num", "bin", "ob_ex", "size_frac", "abs_CC" ])["counts"].sum().reset_index()
        # if ob_ex is expected, "obs_count" will be c.expected_obs_count, otherwise it will be c.observed_obs_count
        plotdf["obs_count"] = [c.num_expected_observations if x == "expected" else c.num_observed_observations for x in plotdf["ob_ex"]]
        # merge these so that the sumdf has both the counts and the obs_count
        plotdf["count_per_sim"] = plotdf["counts"] / plotdf["obs_count"]

        print("the plotdf is:\n{}".format(plotdf))
        plotting = True if len(plotdf) > 0 else False

    if not plotting:
        fig = plt.figure(figsize=(6, 2))
        # turn off all the axes
        plt.axis('off')
        figtitle = "Sorry, we did not observe and fusions on the branch leading up to or within the clade: {} - {}".format(
            taxid, taxon_name)
        plt.text(0.5, 0.5, figtitle, fontsize = 30, horizontalalignment='center',
                 verticalalignment='center', transform=plt.gca().transAxes)
        outfilename = "{}_{}_{}.pdf".format(outfileprefix, taxid, taxon_name)
        # save the plot as a pdf
        fig.savefig(outfilename, bbox_inches="tight")
        plt.close()
    elif plotting:
        # make one big plot with all of the data
        # Make the figure
        fw = 20
        fh = 32
        fig = plt.figure(figsize=(fw, fh))
        axes = []

        #for aligning all the panels
        left1 = 0.6
        left2 = 7.5
        left3 = 14

        # CC, size section below
        bottom = 0.6
        panelheight = 5
        for thisabs in ["abs", "CC"]:
            for thissize in ["frac", "size"]:
                # make a plot of the mean counts
                counts, cbar = gen_square_ax_and_colorbar(left1, bottom, fw, fh, panelheight)
                axes.append(fig.add_axes(counts))
                axes.append(fig.add_axes(cbar))
                temp1, temp2 = generate_mean_counts_panel(
                    axes[-2], axes[-1], plotdf, thissize, thisabs)
                axes[-2] = temp1
                axes[-1] = temp2

                # make the heatmap axes
                counts, cbar = gen_square_ax_and_colorbar(left2, bottom, fw, fh, panelheight)
                axes.append(fig.add_axes(counts))
                axes.append(fig.add_axes(cbar))
                temp1, temp2 = generate_obs_exp_panel(
                    axes[-2], axes[-1], plotdf, thissize, thisabs)

                # generate the square for the trace.
                axes.append(fig.add_axes(gen_square_ax(left3, bottom, fw, fh, panelheight)))
                # add the trace to the trace panel. This is sort of complicated, so add a function just to modify this panel
                axes[-1] = generate_trace_panel(axes[-1], branches_in_clade,
                                                simulation_filepaths, "num", thissize, thisabs)

                # update the bottom
                bottom = bottom + panelheight + 1.1

        # make a plot of the mean ALG counts
        counts, cbar = gen_square_ax_and_colorbar(left1, bottom, fw, fh, panelheight)
        axes.append(fig.add_axes(counts))
        axes.append(fig.add_axes(cbar))
        temp1, temp2 = generate_ALG_mean_counts_panel(
            axes[-2], axes[-1], plotdf, algdf)
        axes[-2] = temp1
        axes[-1] = temp2

        # make a plot of the mean ALG counts
        counts, cbar = gen_square_ax_and_colorbar(left2, bottom, fw, fh, panelheight)
        axes.append(fig.add_axes(counts))
        axes.append(fig.add_axes(cbar))
        temp1, temp2 = generate_ALG_obs_exp_counts_panel(
            axes[-2], axes[-1], plotdf, algdf)
        axes[-2] = temp1
        axes[-1] = temp2

        # generate the square for the trace.
        axes.append(fig.add_axes(gen_square_ax(left3, bottom, fw, fh, panelheight)))
        # add the trace to the trace panel. This is sort of complicated, so add a function just to modify this panel
        axes[-1] = generate_trace_panel(axes[-1], branches_in_clade,
                                        simulation_filepaths, "ALG", "size", "abs")

        # Add a title to the plot to indicate the NCBI taxid and the taxon name
        # Move it to the middle of the plot
        fig_title = "NCBI taxid: {}, taxon name: {}".format(taxid, taxon_name)
        plt.text(-0.75, 1.15, fig_title, fontsize = 30, horizontalalignment='right', verticalalignment='center', transform=plt.gca().transAxes)

        outfilename = "{}_{}_{}.pdf".format(outfileprefix, taxid, taxon_name)
        # save the plot as a pdf
        fig.savefig(outfilename, bbox_inches="tight")
        plt.close()

def read_simulations_and_make_heatmaps(simulation_filepaths, per_sp_df, algdfpath, outfileprefix,
                                       clades_of_interest, absmax = 6):
    """
    This function reads in a file list of simulation files, sums up all the data,
     and makes a heatmap of the results.

    Makes on plot per NCBI taxid that we care about.
    We are able to do this as it is a Monte Carlo
     simulation, so the results are independent.
    """
    # CALL THIS TO GET THE VISUAL STYLE WE NEED
    odp_plot.format_matplotlib()

    # Try to build a tree and print a clade
    T = PhyloTree()
    T.build_tree_from_per_sp_chrom_df(per_sp_df)
    T.add_taxname_to_all_nodes()

    # read in the algdf
    algdf = rbh_tools.parse_ALG_rbh_to_colordf(algdfpath)
    algdf = algdf.sort_values(by=["Size"]).reset_index(drop=True)

    c = coloc_array()
    c.plotmatrix_listoffiles_to_plotmatrix(simulation_filepaths)

    for taxid in clades_of_interest:
        _make_one_simulation_plot(algdf, c, T, taxid, simulation_filepaths,
                                  outfileprefix, absmax = absmax)

def unit_test_coloc_array_identical():
    """
    Does a test to check whether the output of the fusion df and the coloc dfs
    are exactly the same depending on the order in which they were processed.

    20240118 - right now, it appears that the disappearance dataframes are exactly the same.
             - the colocalization dataframes are not exactly the same, but they are very close.
               One difference that I have seen in the color dataframes is in the cnidarians, this
               row for example
               9           6103         86626  [(A1b, A2), (A1b, N), (A2, B1), (A2, B2), (A2,...     []  Metridiumsenile-6116-GCA949775045.1          6116
    """
    # test if these are the same output with different random seeds. Should be exactly the same.
    sampledf = pd.read_csv(sys.argv[1], sep="\t")
    algdf = rbh_tools.parse_ALG_rbh_to_colordf(sys.argv[2])

    dispersion_df, coloc_df, ALG_coloc_df = stats_df_to_loss_fusion_dfs(sampledf, algdf,
                                   obs_seed = 10,
                                   randomize_ALGs=False)

    dispersion_df2, coloc_df2, ALG_coloc_df = stats_df_to_loss_fusion_dfs(sampledf, algdf,
                                   obs_seed = 500,
                                   randomize_ALGs=False)

    ## CHECK THE DISPERSION DFs
    ## sort by thisedge, this loss, loss size. Just keep thisedge, thisloss
    #dispersion_df = dispersion_df.sort_values(by=["thisedge", "thisloss", "loss_size"])[["thisedge", "thisloss"]].reset_index(drop=True)
    #dispersion_df2 = dispersion_df2.sort_values(by=["thisedge", "thisloss", "loss_size"])[["thisedge", "thisloss"]].reset_index(drop=True)
    ##check that the two dfs are the same
    #if not dispersion_df.equals(dispersion_df2):
    #    raise Exception("dispersion_dfs are not the same")

    # CHECK THE COLOCALIZATION DFs
    coloc_df  = coloc_df.sort_values( by=["thisedge", "thiscoloc"])[["thisedge", "thiscoloc"]].reset_index(drop=True)
    coloc_df2 = coloc_df2.sort_values(by=["thisedge", "thiscoloc"])[["thisedge", "thiscoloc"]].reset_index(drop=True)
    # check that the two dfs are the same
    print(coloc_df)
    print(coloc_df2)
    if not coloc_df.equals(coloc_df2):
        # Get the unique rows in coloc_df that are not in coloc_df2.
        # Also get the unique rows in coloc_df2 that are not in coloc_df.
        # Use merge to accomplish this.
        # For example, if these occur in both dfs, remove them from both.
        # df1 7660  (3073854, 30446)    (G, H)
        # df2 7649  (3073854, 30446)    (G, H)
        coloc_unique  = pd.merge(coloc_df, coloc_df2, how="outer", indicator=True).query('_merge=="left_only"').drop('_merge', axis=1)
        coloc2_unique = pd.merge(coloc_df2, coloc_df, how="outer", indicator=True).query('_merge=="left_only"').drop('_merge', axis=1)
        print("coloc_df unique")
        print(coloc_unique)
        print("coloc2 unique")
        print(coloc2_unique)
        raise Exception("coloc_dfs are not the same")


def main():
    #generate_node_taxid_file_from_per_sp_chrom_df(sys.argv[1], "node_taxid.tsv")

    # how many BCnS ALGs correspond to two or more chromosomes
    # obs/exp value for different ALG fusion combinations
    # obs/exp value for number of losses
    # x-axis number of chromosomes, in that species, y-axis is number of fusions leading
    # x-axis number of chromosomes, in that speciers, y-axis is the number of losses
    #  ... and some variation on that.
    clades_of_interest = [6340,     # Annelida
                          6447,     # Mollusca
                          6606,     # Coleoidea
                          33511,    # Deuterostomia
                          33317,    # Protostomia
                          10197,   # Ctenophora
                          33213,   # Bilateria
                          31265,   # Acoela
                          33317,   # Protostomia
                          1206794, # Ecdysozoa
                          6231,    # Nematoda
                          88770,   # Panarthropoda
                          7088,    # Lepidoptera
                          2697495, # Spiralia
                          6544,    # Bivalvia
                          7586,    # Echinodermata
                          6073,    # Cnidaria
                          7711,    # Chordata
                          7742,    # Vertebrata
                          6605,    # cephalopods
                          33208,   # Metazoa
                          6040,    # Porifera
                          6042,    # Demospongiae
                          ]

    # Specify the file path of the TSV file. Use sys.argv[1]. The file will be called something like per_species_ALG_presence_fusions.tsv
    generate_stats_df(sys.argv[1], "statsdf.tsv")

    #run_n_simulations_save_results(sys.argv[1]       ,
    #                               sys.argv[2]       ,
    #                               "testfile.tsv"    ,
    #                               num_sims=10       ,
    #                               abs_bin_size=25   ,
    #                               frac_bin_size=0.05,
    #                               verbose = True    )
    #read_simulations_and_make_heatmaps(["testfile.tsv"], sys.argv[1], sys.argv[2], "simulations", clades_of_interest)

    #sampledf = pd.read_csv(sys.argv[1], sep="\t")
    #algdf = rbh_tools.parse_ALG_rbh_to_colordf(sys.argv[2])
    #dispersion_df, coloc_df, ALG_coloc_df  = stats_df_to_loss_fusion_dfs(
    #                                            sampledf, algdf,
    #                                            obs_seed = 1,
    #                                            randomize_ALGs=False)
    #print("This is a dispersion df")
    #print(dispersion_df)
    #print("This is a coloc df")
    #print(coloc_df)
    #print("This is an ALG coloc df")
    #print(ALG_coloc_df)
    #sys.exit()

    num_simulations = 60
    sims_per_run    = 2
    for i in range(int(num_simulations/sims_per_run)):
        outname = "simulations/sim_results_{}_{}.tsv".format(i, sims_per_run)
        if not os.path.exists(outname):
            print("running simulation {}/{}".format(i+1, int(num_simulations/sims_per_run)))
            run_n_simulations_save_results(sys.argv[1],
                                           sys.argv[2],
                                           outname,
                                           num_sims=sims_per_run,
                                           abs_bin_size=25,
                                           frac_bin_size=0.05,
                                           verbose = True)

    # find all the files in this directory that start with sim_results_ or dfsim_run_
    simulation_filepaths =  list(set(glob.glob("./simulations/sim_results_*.tsv")) | set(glob.glob("./simulations/dfsim_run_*.tsv")))
    print(simulation_filepaths)
    if len(simulation_filepaths) == 0:
        raise Exception("No simulation files found in ./simulations/")
    read_simulations_and_make_heatmaps(simulation_filepaths,
                                       sys.argv[1], sys.argv[2],
                                       "simulations", clades_of_interest)

if __name__ == "__main__":
    main()