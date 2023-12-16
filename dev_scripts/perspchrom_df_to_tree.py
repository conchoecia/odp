#!/usr/bin/env python3
# the file path is the first positional arg. use sys.
import os
import pandas as pd
import random
import sys
from ete3 import NCBITaxa

# import the parse_rbh_file from the plot_ALG_fusions.py script
# this is a function that parses the RBH file into a dataframe
from plot_ALG_fusions import parse_rbh_file
import matplotlib.pyplot as plt

# use networkx to make graphs for the lineage-specific fusion/losses
import networkx as nx

def parse_gain_loss_string(GL_string, samplename):
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
                 "sample_taxid": int(GL_string.split("-")[-1])} # this is an int
        entries.append(entry)
    # convert the list of dicts to a dataframe
    df = pd.DataFrame(entries)
    return df

def parse_gain_loss_from_perspchrom_df(perspchromdf):
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

def stats_on_changedf(changedf):
    """
    TODO
    """
    print(changedf)
    # Right now colocalizations and losses are lists. We want to count the number of colocalizations and losses.
    # So we will make a new column called 'change' and another called 'change_type', then we will unwrap
    # the lists into the new columns. This will make it easier to do stats on the changes.
    entries = []
    for i, row in changedf.iterrows():
        # colocalizations
        for colocalization in row["colocalizations"]:
            entry = {"source_taxid": row["source_taxid"],
                     "target_taxid": row["target_taxid"],
                     "change": colocalization,
                     "change_type": "colocalization",
                     "samplename": row["samplename"],
                     "sample_taxid": row["sample_taxid"]}
            entries.append(entry)
        # losses
        for loss in row["losses"]:
            entry = {"source_taxid": row["source_taxid"],
                     "target_taxid": row["target_taxid"],
                     "change": loss,
                     "change_type": "loss",
                     "samplename": row["samplename"],
                     "sample_taxid": row["sample_taxid"]}
            entries.append(entry)
    # convert the list of dicts to a dataframe
    changedf = pd.DataFrame(entries)
    print(changedf)

    # Because this is a big N-sat problem to figure out the exact branch on which the change happened,
    #  we will now look at whether the change happened right after the source_taxid or right before the
    #  target_taxid.
    # To do this we will groupby on source_taxid, then count the number of changes. We will then groupby on
    #  target_taxid, then count the number of changes. We will then groupby on both source_taxid and target_taxid,
    #  then count the number of changes. For now all we can do is compare these numbers.
    # For each case, we want to independently count the losses or colocalizations.
    # Just keep the changes column in the groupby after count, then change colocalization and loss to their own columns.
    # get right of target_taxid, samplename, sample_taxid
    #groupby_source = changedf.groupby(["source_taxid", "change_type"]).count()
    #groupby_source = groupby_source.drop(["target_taxid", "samplename", "sample_taxid"], axis=1)
    ## turn the groupby columns into a df
    #groupby_source = groupby_source.reset_index()
    #print(groupby_source)
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
    groupby_target["frac_total"] = groupby_target.apply(lambda row: row["counts"]/change_to_total_counts[row["change"]], axis=1)
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

def stats_df_to_loss_fusion_plots(perspchromdf, statsdf, rbhfile):
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
    """
    ALG_info = parse_rbh_file(rbhfile)
    # Make a starting graph where each "ALGname" in ALG_info is a node,
    #  each node has a size attribute (integer), and a color attribute (string)
    G = nx.Graph()
    # add the nodes
    for i, row in ALG_info.iterrows():
        G.add_node(row["ALGname"], size=int(row["Size"]), color=row["Color"])

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

    print(dispersion_df)
    print(coloc_df)

    # first do the colocalizations:
    x = list(coloc_df["coloc0_CC_size"]) + list(coloc_df["coloc1_CC_size"])
    y = list(coloc_df["coloc1_CC_size"]) + list(coloc_df["coloc0_CC_size"])
    # make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    # get the ALG sizes from the rbhdf
    # Create the scatterplot. These dots will be blue.
    #  make the dots 0.1 transparency
    plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    plt.xlabel('ALG Size')
    plt.ylabel('ALG Size')
    plt.title('Component fusion sizes')
    # Display the plot
    plt.show()
    plt.close()

    # Plot as the percent of the size of the longest
    x = list(coloc_df["coloc0_CC_percent_of_largest"]) + list(coloc_df["coloc1_CC_percent_of_largest"])
    y = list(coloc_df["coloc1_CC_percent_of_largest"]) + list(coloc_df["coloc0_CC_percent_of_largest"])
    # make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    # get the ALG sizes from the rbhdf
    # Create the scatterplot. These dots will be blue.
    #  make the dots 0.1 transparency
    plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    #  make the dots 0.1 transparency
    plt.xlabel('CC Size as a fraction')
    plt.ylabel('CC Size as a fraction')
    plt.title('CC_size as fraction of largest in genome')
    # Display the plot
    plt.show()
    plt.close()

    # Plot as the percent of the size of the longest
    x = list(coloc_df["coloc0_CC_percent_of_total"]) + list(coloc_df["coloc1_CC_percent_of_total"])
    y = list(coloc_df["coloc1_CC_percent_of_total"]) + list(coloc_df["coloc0_CC_percent_of_total"])
    # make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    # get the ALG sizes from the rbhdf
    # Create the scatterplot. These dots will be blue.
    #  make the dots 0.1 transparency
    plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    plt.xlabel('CC Size as a fraction')
    plt.ylabel('CC Size as a fraction')
    plt.title('CC_size as fraction of the total genome size')
    # Display the plot
    plt.show()
    plt.close()

    # Plot as the percent of the size of the longest
    x = list(coloc_df["coloc0_size"]) + list(coloc_df["coloc1_size"])
    y = list(coloc_df["coloc1_size"]) + list(coloc_df["coloc0_size"])
    # make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    # get the ALG sizes from the rbhdf
    # Create the scatterplot. These dots will be blue.
    #  make the dots 0.1 transparency
    plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    plt.xlabel('CC Size, absolute')
    plt.ylabel('CC Size, absolute')
    plt.title('CC_size absolute')
    # Display the plot
    plt.show()
    plt.close()

    # Plot as the percent of the size of the longest
    x = list(coloc_df["coloc0_percent_of_total"]) + list(coloc_df["coloc1_percent_of_total"])
    y = list(coloc_df["coloc1_percent_of_total"]) + list(coloc_df["coloc0_percent_of_total"])
    # make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    # get the ALG sizes from the rbhdf
    # Create the scatterplot. These dots will be blue.
    #  make the dots 0.1 transparency
    plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    plt.xlabel('CC Size, absolute')
    plt.ylabel('CC Size, absolute')
    plt.title('CC_size colocalizations as percent of total')
    # Display the plot
    plt.show()
    plt.close()

    # Plot as the percent of the size of the longest
    x = list(coloc_df["coloc0_percent_of_largest"]) + list(coloc_df["coloc1_percent_of_largest"])
    y = list(coloc_df["coloc1_percent_of_largest"]) + list(coloc_df["coloc0_percent_of_largest"])
    # make a scatterplot of the ALG size vs. the number of fusions. These dots will be blue
    # get the ALG sizes from the rbhdf
    # Create the scatterplot. These dots will be blue.
    #  make the dots 0.1 transparency
    plt.scatter(x,y, color="blue", alpha=0.1, lw=0)
    plt.xlabel('CC Size, absolute')
    plt.ylabel('CC Size, absolute')
    plt.title('CC_size colocalizations as percent of largest')
    # Display the plot
    plt.show()
    plt.close()

def main():
    # Specify the file path of the TSV file. Use sys.argv[1]
    file_path = sys.argv[1]
    # Check that the file exists
    if not os.path.exists(file_path):
        raise Exception("File does not exist: {}".format(file_path))

    # Read the TSV file into a pandas dataframe
    df = pd.read_csv(file_path, sep="\t")

    changedf = parse_gain_loss_from_perspchrom_df(df)
    statsdf = stats_on_changedf(changedf)
    # save the stats df to a tsv file with headers
    statsdf.to_csv("statsdf.tsv", sep="\t", index=False)

    stats_df_to_loss_fusion_plots(df, statsdf, sys.argv[2])

if __name__ == "__main__":
    main()