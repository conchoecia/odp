#!/usr/bin/env python
"""
This script takes a newick tree and identifies the divergence time of various nodes in the tree.

"""

import argparse
from collections import Counter,deque
import ete3
import networkx as nx
import os
import pandas as pd
import random
import sys
import yaml
import time
# we need the special newick package to parse this file type
from newick import read

def create_directories_recursive_notouch(path):
    """
    Unlike os.makedirs, this function will not touch a directory if it already exists.
    This is useful for snakemake because it will not re-run a rule if the output already exists.
    """
    parts = os.path.normpath(path).split(os.path.sep)
    # Determine whether to remove the last part of the path.
    # Basically we have to determine if the last part is intended to be a path or a file.
    # If it is a file, then we need to remove it.
    file_endings = [".txt", ".tsv", ".csv", ".yaml"]
    end_is_file = False
    for ending in file_endings:
        if parts[-1].endswith(ending):
            end_is_file = True
            break
    if end_is_file:
        parts = parts[:-1]

    current_path = ""
    for part in parts:
        current_path = os.path.join(current_path, part)
        if not os.path.exists(current_path):
            os.mkdir(current_path)
    # safe return code if done
    return 0

def parse_args():
    """
    The args that we need to parse are:
      - the file path of the newick tree
      - the output prefix for the file, including a prepended path
      - the original config file from which the tree was derived.
      - the taxids of interest, if we want to plot the divergence relative to a specific species. Accepts a list of taxids
      - the email address to use for programmatic access to NCBI
      - file with samples and chromosome sizes.
      - a flag -s that just tells the program to make the species list
    """
    parser = argparse.ArgumentParser(description="This script takes a newick tree and identifies the divergence time of various nodes in the tree.")
    parser.add_argument("-n", "--newick", help="The path to the newick tree file.", required=True)
    parser.add_argument("-C", "--chromosome_sizes", help="The path to the chromosome sizes")
    parser.add_argument("-p", "--prefix", help="The output prefix for the file, including a prepended path if you want another directory.", required=True)
    parser.add_argument("-c", "--config", help="The original config file from which the tree was derived.", required=False)
    parser.add_argument("-t", "--taxids", help="The taxids of interest, if we want to plot the divergence relative to a specific species. Accepts a list of taxids.", required=False)
    parser.add_argument("-e", "--email", help="The email address to use for programmatic access to NCBI.", required=False)
    parser.add_argument("-s", "--species_list", help="Just make the species list.", action="store_true")

    args = parser.parse_args()

    if args.species_list:
        # generate a species list from the config file if it does not yet exist
        # read in the config file using the yaml package
        fig = open(args.config, "r")
        config = yaml.safe_load(fig)
        binomials = set()
        for sp in config["species"].keys():
            # the field is config["species"][sp]["genus"] + " " + config["species"][sp]["species"]
            genus   = config["species"][sp]["genus"]
            species = config["species"][sp]["species"]
            binomials.add(f"{genus} {species}")

        # make a file called species_list.txt if it does not yet exist. Do it in alphabetical order.
        with open("species_list.txt", "w") as f:
            for binomial in sorted(binomials):
                f.write(f"{binomial}\n")

    # check that the newick and config files actually exist
    if not os.path.exists(args.newick):
        # raise an IO error
        raise IOError("The newick file does not exist: {}".format(args.newick))

    # make sure the chromosome_sizes file exists
    if args.chromosome_sizes:
        if not os.path.exists(args.chromosome_sizes):
            raise IOError("The chromosome sizes file does not exist: {}".format(args.chromosome_sizes))

    # optional args
    if args.config:
        # make sure that the config file exists
        if not os.path.exists(args.config):
            # raise an IO error
            raise IOError("The config file does not exist: {}".format(args.config))
        # make sure that the user has also provided an email address
        if not args.email:
            # raise an error
            raise ValueError("If you provide a config file, you must also provide an email address to programmatically access the NCBI taxonomy database.")

    ## If the type of taxids is None, return an empty list. If it is an int, return a list with that int.
    #if

    return parser.parse_args()


def get_lineage(tree, node, lineage = []):
    """
    This function returns the lineage path of this species.
    Uses recursion
    """
    if len(lineage) == 0:
        lineage = [node]
    # first we need to get the lineage path of this species
    ancestor = node.ancestor
    # break condition
    if ancestor == tree:
        return lineage
    # recursion condition
    else:
        lineage = [ancestor] + lineage
        return get_lineage(tree, ancestor, lineage)

def get_all_lineages(tree):
    """
    For any tree, gets the lineage path for all species.
    Returns a dict of the lineages.
    """
    # first we need to get all of the leaves
    leaves = tree.get_leaves()
    # now we need to get the lineage path for each of these leaves
    lineages = {leaf.name: get_lineage(tree, leaf) for leaf in leaves}
    return lineages

def find_common_ancestor_age(sp1_lineage, sp2_lineage):
    """
    Takes two lineages and finds the common ancestor.
    Does this with a Newick lineage that was extracted from two species.

    The lineage data structure for Newick looks like this:
       [Node("'414'"), Node("'421'"), Node("'327'"), Node("'358'"), Node("'329'"), Node("'349'"), Node("'330'"), Node("'286'"),
        Node("'320'"), Node("'321'"), Node("'323'"), Node("'287'"), Node("'297'"), Node("'294'"), Node("'143'"), Node("'144'"),
        Node("'146'"), Node("'22'"), Node("'23'"), Node("Abscondita_terminalis")]

    These lineage data structures are extracted from the get_all_lineages() function.

    The return type of this function is a tuple of the common ancestor and the age of the species.
    """

    # first we need to find the common ancestor.
    # just compare the two lists until they don't match anymore
    # the most recent match is the common ancestor
    common_ancestor = None
    for i in range(len(sp1_lineage)):
        if sp1_lineage[i] != sp2_lineage[i]:
            common_ancestor = sp1_lineage[i-1]
            break
    shared_species = set(sp1_lineage).intersection(set(sp2_lineage))
    unique_sp1 = [x for x in sp1_lineage if x not in shared_species]
    unique_sp2 = [x for x in sp2_lineage if x not in shared_species]
    sp1_age = sum([x.length for x in unique_sp1] ) #+ [common_ancestor.length])
    sp2_age = sum([x.length for x in unique_sp2] ) #+ [common_ancestor.length])
    # the ages should be the same, so check
    # Sometimes when one of the species has a really long branch, the ages are not exactly the same.
    # Just check that they are within 0.05 of each other in terms of precent.
    percent_diff = 0 if abs(sp1_age - sp2_age) == 0 else (abs(sp1_age - sp2_age)/sp1_age)

    # There is a weird behavior where, if percent_diff is 0, then the equality statement doesn't work as predicted.
    #  So we need to handle that case separately

    ## Getting rid of this and just assuming that the tree is correct
    #sp_1_2_within_0_0_5 = True if (percent_diff == 0) else (percent_diff < (0.15 * sp1_age))
    #if not sp_1_2_within_0_0_5:
    #    print("The two species are: {} and {}".format(sp1_lineage[-1].name, sp2_lineage[-1].name))
    #    print("The lineage of sp1 is: {}".format([x.name for x in sp1_lineage]))
    #    print("The lineage of sp2 is: {}".format([x.name for x in sp2_lineage]))
    #    print("The common ancestor is: {}".format(common_ancestor.name))
    #    print("The percent difference is: {}".format(percent_diff))
    #    raise ValueError("The ages of the two species are not the same: {} vs {}".format(sp1_age, sp2_age))
    return common_ancestor, sp1_age

def get_divergence_time_all_vs_all(tree):
    """
    Takes a tree and gets the divergence times for all species.
    Returns this information as a dictionary.
    """
    lineages = get_all_lineages(tree)

    sp_list_sorted = list(sorted(lineages.keys()))
    for i in range(len(sp_list_sorted)-1):
        for j in range(i+1, len(sp_list_sorted)):
            sp1 = sp_list_sorted[i]
            sp2 = sp_list_sorted[j]
            common_ancestor, age = find_common_ancestor_age(lineages[sp1], lineages[sp2])
            # round this to 5 decimal places
            age_print = round(age, 5)
            yield sp1, sp2, age_print

def report_divergence_time_all_vs_all(tree, output_prefix):
    """
    This method gets the divergence times and writes them to a file with the prefix.
    Also safely makes a directory if it does not yet exist.

    Makes a dictionary of the divergence times for all the species
    """
    # first come up with the outfile path
    outfile_path = "{}.divergence_times.txt".format(output_prefix)
    # data structure to save the divergence times
    divergence_times = {}
    # safely make the directories if they don't exist
    create_directories_recursive_notouch(outfile_path)
    # open the outfile for writing
    with open(outfile_path, "w") as f:
        for sp1, sp2, age in get_divergence_time_all_vs_all(tree):
            entry = (sp1, sp2)
            if entry not in divergence_times:
                divergence_times[entry] = age
            f.write("{}\t{}\t{}\n".format(sp1, sp2, age))
    return divergence_times

def convert_ncbi_entry_to_dict(ncbi_entry):
    entries = []
    for entry in ncbi_entry["LineageEx"]:
        tempdict = {}
        tempdict["TaxID"] =          int(entry["TaxId"])
        tempdict["ScientificName"] = str(entry["ScientificName"])
        tempdict["Rank"] =           str(entry["Rank"])
        entries.append(tempdict)
    new_dict = {"TaxID":          int(ncbi_entry["TaxId"]),
                "ScientificName": str(ncbi_entry["ScientificName"]),
                "Lineage":        str(ncbi_entry["Lineage"]),
                "LineageEx":      entries,
                }
    return new_dict

def get_taxonomy_info(binomial_name):
    handle = Entrez.esearch(db="taxonomy", term=binomial_name, retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    if records["Count"] == "0":
        return "Species not found."

    taxon_id = records["IdList"][0]
    handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
    record = Entrez.read(handle)
    handle.close()
    resultsdict = convert_ncbi_entry_to_dict(record[0])
    return resultsdict

def taxinfo_download_or_load(binomial_name, taxinfo_filepath):
    """
    This looks to see if a yaml file exists with the taxinfo for this species.
    If it does not, it will download the taxinfo from NCBI and save it to that yaml file.

    Sometimes the download from NCBI doesn't work, so we need to allow for failures.

    If it doesn't work, returns a 1.
    If the file exists, returns a 0.
    """
    if not os.path.exists(taxinfo_filepath):
        # safely make the directory if it doesn't exist
        create_directories_recursive_notouch(taxinfo_filepath)
        try:
            sp_tax_info = get_taxonomy_info(binomial_name)
            # now we need to write this to a yaml file
            with open(taxinfo_filepath, "w") as f:
                yaml.dump(sp_tax_info, f)
            # we need to pause if we had a successful download to avoid overloading the NCBI servers
            time.sleep(3)
            # return success
            return 0
        except:
            # return failure
            print("           ^^^ THE DOWNLOAD FOR THIS SPECIES DIDN'T WORK IN THIS ROUND.")
            return 1
    else:
        # read in the file and check if it has any contents.
        # If not, this hasn't worked, we delete the file, then return 1
        with open(taxinfo_filepath, "r") as f:
            contents = f.read()
            if len(contents) == 0:
                os.remove(taxinfo_filepath)
                return 1
            else:
                # in theory the file should be good, so return success
                return 0

def yaml_file_legal(filepath):
    """
    Returns True if the yaml file exists and has contents.
    Returns False otherwise.
    """
    with open(filepath, "r") as f:
        contents = f.read()
        if len(contents) == 0:
            return False
        else:
            # in theory the file should be good, so return success
            return True

class TaxNode:
    """
    one node of the taxonomy tree
    """
    __slots__ = ['taxid', 'parent', 'children',
                 'name', 'nodeages', 'nodeage',
                 'nodeageinterpolated', 'lineage', 'lineage_string',
                 'sort_order', 'dist_crown', 'dist_crown_plus_root',
                 "chromsize_median", "chromsize_mean", "chromsize_list"]
    def __init__(self, taxid) -> None:
        self.taxid = taxid
        self.parent = None
        self.children = set()
        self.name = None
        # node age estimates, in millions of years ago
        self.nodeages = Counter()
        # This is the singular node age estimate, once we are happy with one of them
        self.nodeage = None
        self.nodeageinterpolated = False
        # Get the lineage
        self.lineage = None
        self.lineage_string = ""
        self.sort_order = None
        # path information
        # This is the distance of all the edges from this node to the tips
        self.dist_crown = None
        # This is the distance of all the edges from this node to the tips,
        #  plus the distance of the edge leading up to this node.
        self.dist_crown_plus_root = None
        # chromsize info
        self.chromsize_median = -1
        self.chromsize_mean   = -1
        self.chromsize_list   = []

    def __str__(self) -> str:
        outstring  = "TaxNode:\n"
        outstring += "  - taxis: {}\n".format(self.taxid)
        outstring += "  - parent: {}\n".format(self.parent)
        outstring += "  - children: {}\n".format(self.children)
        outstring += "  - name: {}\n".format(self.name)
        outstring += "  - nodeages: {}\n".format(self.nodeages)
        outstring += "  - nodeage: {}\n".format(self.nodeage)
        outstring += "  - nodeageinterpolated: {}\n".format(self.nodeageinterpolated)
        outstring += "  - lineage: {}\n".format(self.lineage)
        outstring += "  - lineage_string: {}\n".format(self.lineage_string)
        outstring += "  - sort_order: {}\n".format(self.sort_order)
        outstring += "  - dist_crown: {}\n".format(self.dist_crown)
        outstring += "  - dist_crown_plus_root: {}\n".format(self.dist_crown_plus_root)
        return outstring

class TaxEdge:
    """
    One edge of the taxonomy tree.
    Useful for recording path lengths and other properties.
    Don't use it for navigating the graph. Just use the edges
    """
    __slots__ = ['parent_taxid', 'child_taxid',
                 'parent_age', 'child_age',
                 'branch_length',
                 'dist_crown_plus_this_edge',
                 'parent_lineage', 'child_lineage']
    def __init__(self, parent_taxid, child_taxid) -> None:
        self.parent_taxid = parent_taxid
        self.child_taxid = child_taxid
        self.parent_age = None
        self.child_age = None
        self.branch_length = None
        self.dist_crown_plus_this_edge = None
        self.parent_lineage = None
        self.child_lineage  = None

    def __str__(self) -> str:
        outstring  = "TaxEdge:\n"
        outstring += "  - parent_taxid: {}\n".format(self.parent_taxid)
        outstring += "  - child_taxid: {}\n".format(self.child_taxid)
        outstring += "  - branch_length: {}\n".format(self.branch_length)
        outstring += "  - parent_age: {}\n".format(self.parent_age)
        outstring += "  - child_age: {}\n".format(self.child_age)
        return outstring

class TaxIDtree:
    """
    This is a datastructure to quickly search for the most closely related species
      given search species 1 and a tree of species 2...N.
    """
    __slots__ = ["nodes", "edges", "root", "leaf_order"]
    def __init__(self) -> None:
        self.nodes = {}
        self.edges = {}
        self.root = None
        self.leaf_order = []

    def find_root(self) -> int:
        """
        Finds the root and returns its int. The root is defined as something without a parent.
        """
        potential_roots = []
        for node in self.nodes:
            if self.nodes[node].parent is None:
                potential_roots.append(node)
        if len(potential_roots) == 0:
            raise ValueError("There is no root in this tree.")
        elif len(potential_roots) > 1:
            raise ValueError("There is more than one root in this tree.")
        else:
            return potential_roots[0]

    def __str__(self) -> str:
        """
        Just the print method. Makes kind of a tree structure.
        """
        # find the root
        root = self.find_root()
        self.root = root
        newoutstring = "- root\n"
        # now make a recursive algorithm to print the whole tree
        def print_tree(node, outstring, level):
            for child in node.children:
                outstring += "{}|_ {}\n".format("  "*level, child)
                outstring = print_tree(self.nodes[child], outstring, level+1)
            return outstring
        outstring = print_tree(self.nodes[root], newoutstring, 1)
        return outstring

    def add_node(self, taxid) -> TaxNode:
        """
        Takes a taxid, name, and lineage and adds it to the tree.
        """
        if taxid not in self.nodes:
            self.nodes[taxid] = TaxNode(taxid)
        return self.nodes[taxid]

    def add_edge(self, parent_taxid, child_taxid) -> int:
        """
        Adds an edge between two nodes.
        """
        if not parent_taxid in self.nodes:
            self.nodes[parent_taxid] = TaxNode(parent_taxid)
        if not child_taxid in self.nodes:
            self.nodes[child_taxid] = TaxNode(child_taxid)
        # From the parent, add the child
        self.nodes[parent_taxid].children.add(child_taxid)
        # from the child, add the parent. If there is no parent, fine. If there is a parent already, must match the existing taxid.
        if self.nodes[child_taxid].parent is not None:
            if self.nodes[child_taxid].parent != parent_taxid:
                raise ValueError("The child already has a parent, and it is different from what we thought.")
        else:
            self.nodes[child_taxid].parent = parent_taxid

        edgekey = (parent_taxid, child_taxid)
        if edgekey not in self.edges:
            self.edges[edgekey] = TaxEdge(parent_taxid, child_taxid)
        return 0

    def add_chromosome_info_file(self, chrominfo_file):
        """
        Adds chromosome info to the nodes from a file.
        The format of the file is:
        f"%s\t%d" % (sample_string, num_chromosomes)

        The algorithm is: Load the dataframe
        """
        chromdf = pd.read_csv(chrominfo_file, sep="\t", header=None)
        # the headers are sample_string, num_chromosomes
        chromdf.columns = ["sample_string", "num_chromosomes"]
        # only get the sample_strings that have exactly two '-' characters
        chromdf = chromdf[chromdf["sample_string"].str.count("-") == 2]
        chromdf["taxid"] = chromdf["sample_string"].str.split("-").str[1]
        missing = 0
        for i, row in chromdf.iterrows():
            taxid = int(row["taxid"])
            if taxid in self.nodes:
                self.nodes[taxid].chromsize_list.append(row["num_chromosomes"])
            else:
                missing += 1
        print("There were {} missing taxids in the graph taht were present in the genomes".format(missing))

        # now do a reverse BFS to sum the chromosome lists at each node
        counter = 1
        for node in self.nodes:
            print("  - Setting up the chromsizes in the node: {} of {}".format(counter, len(self.nodes)), end="\r")
            newlist = []
            # do a dfs to get all the tips in this clade
            stack = [node]
            while len(stack) > 0:
                current = stack.pop()
                if len(self.nodes[current].children) == 0:
                    newlist.extend(self.nodes[current].chromsize_list)
                else:
                    stack.extend(self.nodes[current].children)
            self.nodes[node].chromsize_list = newlist
            counter += 1
        print()

        # Go through all the nodes and calculate the mean and median chromsize
        for node in self.nodes:
            if len(self.nodes[node].chromsize_list) > 0:
                self.nodes[node].chromsize_mean = sum(self.nodes[node].chromsize_list)/len(self.nodes[node].chromsize_list)
                self.nodes[node].chromsize_median = sorted(self.nodes[node].chromsize_list)[len(self.nodes[node].chromsize_list)//2]
            else:
                self.nodes[node].chromsize_mean = -1
                self.nodes[node].chromsize_median = -1

    def set_leaf_ages_to_zero(self) -> None:
        """
        Sets the leaf ages to zero
        """
        for node in self.nodes:
            if len(self.nodes[node].children) == 0:
                self.nodes[node].nodeages.update([0])
                self.nodes[node].nodeage = 0

    def find_closest_relative(self, NCBI_object, query_taxid) -> int:
        """
        Given any NCBI taxid - find the most closely related species in this tree.
        Do this by going up the lineage of the NCBI taxid until we find a common ancestor in the tree.
        This strictly uses NCBI taxids.

        Inputs:
        - NCBI_object - is the output of ete3.NCBITaxa(). Use this to get the lineage information without having to download it.
        - query_taxid - the taxid of the species we are looking for.
        """
        query_lineage = NCBI_object.get_lineage(query_taxid)[::-1]
        #print("The query lineage is: {}".format(query_lineage))
        # go through the query_lineage and return the first taxid that is in the tree
        target_node = None
        for i in range(len(query_lineage)):
            if query_lineage[i] in self.nodes:
                target_node = query_lineage[i]
                break
        # Now that we found something in common, we go to the tips until we hit the end. We return the id of the first tip we find.
        #  We can make a best guess, so we shuffle the children.
        while len(self.nodes[target_node].children) > 0:
            children = list(self.nodes[target_node].children)
            random.shuffle(children)
            target_node = children[0]
        return target_node

    def get_lineage(self, taxid) -> list:
        """
        traversal to the root. Returns a list of taxids from root-> tip.
        """
        lineage = []
        current_taxid = taxid
        while current_taxid is not None:
            lineage.append(current_taxid)
            current_taxid = self.nodes[current_taxid].parent
        return lineage[::-1]

    def find_LCA(self, taxid1, taxid2) -> int:
        """
        Finds the lowest common ancestor of two species in the tree.
        """
        if taxid1 == taxid2:
            return taxid1
        lineage1 = self.get_lineage(taxid1)
        lineage2 = self.get_lineage(taxid2)
        # find the common ancestor
        common_ancestor = None
        for i in range(len(lineage1)):
            if lineage1[i] != lineage2[i]:
                common_ancestor = lineage1[i-1]
                break
        return common_ancestor

    def percolate_acceptable_ages(self) -> None:
        """
        This function goes through and finds acceptable ages for each node.
        The value for each node we need to optimize is self.nodeage.
        The value of each node must come from self.nodeages, must be less than the parent self.nodeage,
            and must be greater than whatever value the childrens' self.nodeage values have.

        Some algorithm ideas are to start from the root and BFS
        """
        # find the root and set it for the object
        self.root = self.find_root()
        # first, find all the nodes that are missing ages, that have no children. These are leaves and they should have a time of 0
        for node in self.nodes:
            if len(self.nodes[node].nodeages) == 0 and len(self.nodes[node].children) == 0:
                self.nodes[node].nodeages.update([0])
        # There are a few nodes that we force the dates of, because they are the root of the tree.
        # 1       - the root of the tree - 3.7   billion years ago
        # 131567  - cellular organisms   - 3.48  billion years ago
        # 2759    - eukaryotes           - 1.898 billion years ago
        # 33154   - opisthokonta         - 1.010 billion years ago
        changes = {1: 3700, 131567: 3480, 2759: 1898, 33154: 1010}
        for node in changes:
            if node in self.nodes:
                self.nodes[node].nodeages.update([changes[node]])

        # First, just set every node to the most common value.
        for node in self.nodes:
            if len(self.nodes[node].nodeages) > 0:
                self.nodes[node].nodeage = self.nodes[node].nodeages.most_common(1)[0][0]

        # There are now a lot of nodes that have ages, but we need to make sure that the children are younger than the parents.
        # Do a BFS from the root.
        iteration = 0
        done = False
        while not done:
            # We do one of the descending routines
            if iteration % 2 == 0:
                print("We are percolating from root to tip")
                self.percolate_descending_root_to_tip()
            else:
                raise ValueError("We have not implemented this yet.")
            descending = self.check_nodeages_descending()
            if descending:
                done = True
            print("We are on iteration {}. The descending state is {}".format(iteration, descending))
            iteration += 1

    def percolate_descending_root_to_tip(self) -> None:
        """
        This does a pass through the tree and attempts to sort things in descending order from the root to the tip.
        Uses a BFS.
        """
        queue = [self.root]
        counter = 0
        while len(queue) > 0:
            counter += 1
            parent = queue.pop(0)
            # now add the kids to the queue
            for child in self.nodes[parent].children:
                queue.append(child)
            # ** Updating the actual selected age! **
            # for now, just pick a parent that is older than all of the children and younger than the parent
            # If there is no parent, just get the oldest root age.
            grandparent = self.nodes[parent].parent
            if grandparent == None:
                continue
            else:
                # enforce that we already set an age for the grandparent if it is available
                if len(self.nodes[grandparent].nodeages) != 0:
                    if self.nodes[grandparent].nodeage is None:
                        raise ValueError("The grandparent nodeage is None. This is supposed to be set in percolate_acceptable_ages()")
                if len(self.nodes[parent].nodeages) != 0:
                    if self.nodes[parent].nodeage is None:
                        raise ValueError("The parent nodeage is None. This is supposed to be set in percolate_acceptable_ages()")
                for child in self.nodes[parent].children:
                    if len(self.nodes[child].nodeages) != 0:
                        if self.nodes[child].nodeage is None:
                            raise ValueError("The child nodeage is None. This is supposed to be set in percolate_acceptable_ages()")

                # Now we can continue to see if we can adjust the current value
                # short for "grandparent has ages", "parent has ages", "children have ages"
                gha = len(self.nodes[grandparent].nodeages) > 0
                pha = len(self.nodes[parent].nodeages) > 0
                cha = any([len(self.nodes[child].nodeages) > 0 for child in self.nodes[parent].children])
                # There are a few cases to deal with
                # No Falses
                if [gha, pha, cha] == [True,  True,  True ]:
                    # All three have node ages. There is a chance to fix the parent node age.
                    maxchild = max([x for x in [self.nodes[child].nodeage for child in self.nodes[parent].children] if x is not None] + [0])
                    if (self.nodes[parent].nodeage <= self.nodes[grandparent].nodeage) and (self.nodes[parent].nodeage >= maxchild):
                        # in this case, everything is good. the parent is older than the children and younger than the grandparent.
                        pass
                    else:
                        for k in self.nodes[parent].nodeages:
                            if (k <= self.nodes[grandparent].nodeage) and (k >= maxchild):
                                self.nodes[parent].nodeage = k
                                break
                # One False
                if [gha, pha, cha] == [False, True,  True ]:
                    maxchild = max([x for x in [self.nodes[child].nodeage for child in self.nodes[parent].children] if x is not None] + [0])
                    if self.nodes[parent].nodeage >= maxchild:
                        # in this case we don't care about the grandparent, pass
                        pass
                    else:
                        for k in self.nodes[parent].nodeages:
                            if k >= maxchild:
                                self.nodes[parent].nodeage = k
                                break
                if [gha, pha, cha] == [True,  False, True ]:
                    # we don't do anything here, because we don't have any values to pick from
                    pass
                if [gha, pha, cha] == [True,  True,  False]:
                    # in this case we don't care about the children, pass
                    if self.nodes[parent].nodeage <= self.nodes[grandparent].nodeage:
                        pass
                    else:
                        for k in self.nodes[parent].nodeages:
                            if k <= self.nodes[grandparent].nodeage:
                                self.nodes[parent].nodeage = k
                                break
                # Two Falses
                if [gha, pha, cha] == [False, False, True ]:
                    # we don't do anything here, because we don't have any values to pick from
                    pass
                if [gha, pha, cha] == [True,  False, False]:
                    # we don't do anything here, because we don't have any values to pick from
                    pass
                if [gha, pha, cha] == [False, True, False ]:
                    # we don't do anything here either, because there are no limiting bounds above or below
                    pass
                # Three Falses
                if [gha, pha, cha] == [False, False, False]:
                    # we don't do anything here
                    pass
        if counter != len(self.nodes):
            raise ValueError("The number of nodes visited is not the same as the number of nodes in the tree. Counter = {}, NumNodes = {}".format(
                counter, len(self.nodes)))

    def check_nodeages_descending(self, enforce_nomissing = False) -> bool:
        """
        This returns True if the nodeages are descending from the root to the tips.

        There is an option to enforce that there are no missing nodeages.
        Does not allow for identical node ages

        Returns True if all of the node ages are in descending order from the root to the tip.
         I.E. - None of the children's nodeages can be larger than the parent's nodeage.
        """
        queue = [self.root]
        while len(queue) > 0:
            parent = queue.pop(0)
            if enforce_nomissing:
                if self.nodes[parent].nodeage is None:
                    return False
            for child in self.nodes[parent].children:
                queue.append(child)
                #print("The parent is: {}. The parent ages are {}.".format(parent, self.nodes[parent].nodeages))
                #print("  - the child is {}. The child ages are {}".format(child, self.nodes[child].nodeages))
                # Here we actually check if the nodeages are in descending order.
                if (self.nodes[child].nodeage is None) or (self.nodes[parent].nodeage is None):
                    if enforce_nomissing:
                        raise ValueError("The parent or child nodeage is None.\n  - Parent nodeage is {}\n  - Child nodeage is {}".format(
                            self.nodes[parent].nodeage, self.nodes[child].nodeage))
                    else:
                        pass
                else:
                    # We are checking descending
                    if self.nodes[child].nodeage > self.nodes[parent].nodeage:
                        print("  ** - Parent-child is: {}-{}, ages are: {}-{}".format(
                            parent,                    child,
                            self.nodes[parent], self.nodes[child]))
                        print("     - The parent is: {}. The parent ages are {}.".format(parent, self.nodes[parent].nodeages))
                        print("     - The child is {}. The child ages are {}".format(child, self.nodes[child].nodeages))
                        return False
        return True

    def interpolate_nodes(self) -> None:
        """
        Interpolates nodes with missing values
        """
        # now, find all of the nodes that are missing ages, but have children. These are the ones we need to interpolate.
        # Use BFS to find all of the missing nodes.
        queue = [self.root]
        while len(queue) > 0:
            parent = queue.pop(0)
            for child in self.nodes[parent].children:
                queue.append(child)
            if len(self.nodes[parent].nodeages) == 0:
                # This is missing, we need to figure out how to interpolate it.
                missing_node = parent
                sublineage = [missing_node]
                sublineage = [self.nodes[sublineage[0]].parent] + sublineage
                while len(self.nodes[sublineage[0]].nodeages) == 0:
                    sublineage = [self.nodes[sublineage[0]].parent] + sublineage
                # we have now found a parent with an age. We can now recursively DFS until we find a child with an age.
                # The return is a list of the children, and the last one will have a node age
                sublineage = sublineage + self.find_children_with_ages(sublineage[-1])
                print("We looked at the missing node {}".format(sublineage[0]))
                print("  - the sublineage is: {}".format(sublineage))
                print("  - The nodeages are {}".format([self.nodes[x].nodeages for x in sublineage]))
                # Now we need to interpolate the ages of the missing nodes
                #  We subtract the age of sublineage[-1] from the sublineage[0], then step linearly for the missing nodes.
                print("The sublineage is: {}".format(sublineage))
                print("  - The nodeages of sublineage[0] are: {}".format(self.nodes[sublineage[0]].nodeages))
                print("  - The nodeages of sublineage[-1] are: {}".format(self.nodes[sublineage[-1]].nodeages))
                oldage   = self.nodes[sublineage[0]].nodeages.most_common(1)[0][0]
                youngage = self.nodes[sublineage[-1]].nodeages.most_common(1)[0][0]
                agedif  = oldage - youngage
                steps_to_take = len(sublineage)
                age_per_step  = agedif / steps_to_take
                # we only modify the parent node... the one we are currently looking at
                self.nodes[parent].nodeages.update([oldage - age_per_step])
                self.nodes[parent].nodeageinterpolated = True

        # make sure that none of the nodes are missing ages
        missing_nodes = [node for node in self.nodes if len(self.nodes[node].nodeages) == 0]
        if len(missing_nodes) > 0:
            raise ValueError("There are still missing nodes in the tree.")

        # Now set all the values for the interpolated
        for node in self.nodes:
            if self.nodes[node].nodeageinterpolated:
                self.nodes[node].nodeage = self.nodes[node].nodeages.most_common(1)[0][0]

        done = False
        while not done:
            counts = self.fix_broken_interpolated_entries()
            print("We fixed {} broken interpolated entries.".format(counts))
            if counts == 0:
                done = True

    def fix_broken_interpolated_entries(self) -> int:
        """
        Fixes broken interpolated entries. Returns the number of changes
        """
        # Now we should check that any pairs that have interpolated values are appropriately set,
        #   where the parent is older than the child.
        queue = [self.root]
        counter = 0
        while len(queue) > 0:
            parent = queue.pop(0)
            if parent == 6381:
                print("We're in here")
                print("grandparent", self.nodes[self.nodes[parent].parent])
                print("parent", self.nodes[parent])
            for child in self.nodes[parent].children:
                queue.append(child)
                if self.nodes[parent].nodeageinterpolated or self.nodes[child].nodeageinterpolated:
                    if parent == 6381:
                        print(self.nodes[child])
                    # check if the parent is older than the child
                    parent_nodeage = self.nodes[parent].nodeage
                    child_nodeage  = self.nodes[child].nodeage
                    if parent_nodeage < child_nodeage:
                        # two options, if there is a grantparent, we make sure the new value is less
                        if self.nodes[parent].parent is not None:
                            grandparent_nodeage = self.nodes[self.nodes[parent].parent].nodeage
                            if child_nodeage + 1 < grandparent_nodeage:
                                self.nodes[parent].nodeage = child_nodeage + 1
                                self.nodes[parent].nodeages = Counter({child_nodeage + 1: 1})
                                counter += 1
                            else:
                                self.nodes[parent].nodeage = child_nodeage
                                self.nodes[parent].nodeages = Counter({child_nodeage: 1})
                                counter += 1
                        else:
                            self.nodes[parent].nodeage = child_nodeage + 1
                            self.nodes[parent].nodeages = Counter({child_nodeage + 1: 1})
                            counter += 1
                    if parent == 6381:
                        print("The parent should be updated")
                        print(self.nodes[parent])
        return counter

    def correct_missing_nodes(self) -> None:
        """
        Goes through the tree and interpolates the node ages based on the ages of nodes above and below it.
        This only assigns ages to the nodes that are missing info.
        """
        # find the root and set it for the object
        self.root = self.find_root()
        # first, find all the nodes that are missing ages, that have no children. These are leaves and they should have a time of 0
        for node in self.nodes:
            if len(self.nodes[node].nodeages) == 0 and len(self.nodes[node].children) == 0:
                self.nodes[node].nodeages.update([0])
        # There are a few nodes that we force the dates of, because they are the root of the tree.
        # 1       - the root of the tree - 3.7   billion years ago
        # 131567  - cellular organisms   - 3.48  billion years ago
        # 2759    - eukaryotes           - 1.898 billion years ago
        # 33154   - opisthokonta         - 1.010 billion years ago
        changes = {1: 3700, 131567: 3480, 2759: 1898, 33154: 1010}
        for node in changes:
            if node in self.nodes:
                self.nodes[node].nodeages.update([changes[node]])

        # There are now a lot of nodes that have ages, but we need to make sure that the children are older than the parents.
        # Do a BFS from the root.
        queue = [self.root]
        while len(queue) > 0:
            parent = queue.pop(0)
            acceptable_parent_ages = [x for x in self.nodes[parent].nodeages]
            all_child_ages = []
            #print("The parent is: {}. The parent ages are {}.".format(parent, self.nodes[parent].nodeages))
            for child in self.nodes[parent].children:
                if len(acceptable_parent_ages) > 0:
                    # This only works if there are some parent ages to compare to
                    # if the child has an age, we need to make sure that it is older than the parent.
                    child_ages = [x for x in self.nodes[child].nodeages]
                    all_child_ages.append(child_ages)
                    # remove the acceptable parent ages that are smaller than all of the child ages
                    if len(child_ages) > 0:
                        # we need this condition because the logic doesn't work if the child ages list is empty
                        remove_parent_age = [p for p in acceptable_parent_ages if all([p < c for c in child_ages])]
                    else:
                        remove_parent_age = []
                    acceptable_parent_ages = [p for p in acceptable_parent_ages if p not in remove_parent_age]
                    #print("  - the child is {}. The child ages are {}".format(
                    #  child, self.nodes[child].nodeages))
                queue.append(child)
            # If there were some ages to start with, but we find that none of them are acceptable, raise an error.
            if (len(acceptable_parent_ages) == 0) and len(self.nodes[parent].nodeages) > 0:
                raise ValueError("The parent {} has no acceptable ages.\nThe original parent ages were {}.\nThe original child ages were {}.".format(
                    parent, self.nodes[parent].nodeages, all_child_ages))
            # If we didn't hit an error, actually take the acceptable parent ages and update the parent node. We still want a counter object
            self.nodes[parent].nodeages = Counter({age: count for age, count in self.nodes[parent].nodeages.items() if age in acceptable_parent_ages})

        # Set the leave ages to zero
        self.set_leaf_ages_to_zero()

        # We need to adjust the ages of the children to make sure that they are older than the parent.
        self.percolate_acceptable_ages()
        print("We found that there are options where all children are younger than the parent. This is good.")

        # Set the leave ages to zero
        self.set_leaf_ages_to_zero()

        # now, find all of the nodes that are missing ages, but have children. These are the ones we need to interpolate.
        self.interpolate_nodes()
        # Set the leave ages to zero
        self.set_leaf_ages_to_zero()

        # We need to adjust the ages of the children to make sure that they are older than the parent.
        self.percolate_acceptable_ages()
        print("We found that there are options where all children are younger than the parent. This is good.")
        # Set the leave ages to zero
        self.set_leaf_ages_to_zero()

        # We should check that the path of all the node ages from the root to the tip is the same.
        # If there is not one, raise an error for now. Everything should be the same.
        # We can do a DFS to check this. We want the whole lineage path recorded.
        lineage_to_length = []
        stack = [[self.root]]
        while len(stack) > 0:
            path = stack.pop()
            node = path[-1]
            if len(self.nodes[node].children) == 0:
                # this is a tip, so we need to record the length of the paths
                lineage_to_length.append(([path], int(self.get_lineage_length(path))))
            else:
                for child in self.nodes[node].children:
                    stack.append(path + [child])
        print("The lineage to length is: {}".format(Counter([x[1] for x in lineage_to_length])))

    def get_lineage_length(self, lineage) -> float:
        """
        Given a list of nodes, this function returns the length of the lineage.
        It does so by getting the edge length between each node and summing them.
        """
        totlen = 0
        for i in range(len(lineage)-1):
            parent = lineage[i]
            child  = lineage[i+1]
            # if the child is not a child of the parent, raise an error
            if child not in self.nodes[parent].children:
                raise ValueError("The child {} is not a child of the parent {}".format(child, parent))
            distance = self.nodes[parent].nodeage - self.nodes[child].nodeage
            # check that the distance is positive
            if distance < 0:
                raise ValueError("The distance between {} and {} is negative: {}".format(parent, child, distance))
            totlen += distance
        return totlen

    def get_dominant_age(self,node) -> float:
        """
        This function returns the most common age of the node.
        The node ages are stored in a Counter object.
        """
        return self.nodes[node].nodeages.most_common(1)[0][0]

    def find_children_with_ages(self, node) -> list:
        """
        This is a recursive function that finds the first child with an age.
        """
        if len(self.nodes[node].nodeages) > 0:
            return [node]
        else:
            for child in self.nodes[node].children:
                return [node] + self.find_children_with_ages(child)

    def ensure_all_leaves_have_age_zero(self) -> None:
        """
        Make sure that all of the leaves have an age of 0.
        """
        leaves_without_zero = [node for node in self.nodes if len(self.nodes[node].children) == 0 and 0 not in self.nodes[node].nodeages]
        print("These are the leaves without zeros:")
        for node in leaves_without_zero:
            print("  - {}".format(self.nodes[node]))
        if len(leaves_without_zero) > 0:
            raise ValueError("The node {} is a leaf, but it does not have an age of 0.".format(node))

    def calc_edges(self) -> None:
        """
        From the graph of nodes, calculates the edges. Completely resets the object.
        """
        self.edges = {}
        self.root = self.find_root()
        queue = [self.root]
        while len(queue) > 0:
            parent = queue.pop(0)
            for child in self.nodes[parent].children:
                queue.append(child)
                self.edges[(parent, child)] = TaxEdge(parent, child)
                self.edges[(parent, child)].parent_age = self.nodes[parent].nodeage
                self.edges[(parent, child)].child_age  = self.nodes[child].nodeage
                parent_nodeage = self.nodes[parent].nodeage
                child_nodeage  = self.nodes[child].nodeage
                if (parent_nodeage is None) or (child_nodeage is None):
                    raise ValueError("The parent {} or child {} nodeage is None. This is supposed to be set in percolate_acceptable_ages()".format(
                        parent_nodeage, child_nodeage))
                self.edges[(parent, child)].branch_length = self.nodes[parent].nodeage - self.nodes[child].nodeage

    def calc_dist_crown(self) -> None:
        """
        For each node, calculate the distance of all the edges from this node to the tips.
        """
        # we have to calculate the edges
        self.calc_edges()

        visited = set()
        queue = deque([node for node in self.nodes if len(self.nodes[node].children) == 0])

        while len(queue) > 0:
            # print the results to make sure we're progressing
            print("    The queue length is: {}. First value is {}.".format(len(queue), queue[0]), end="\r")
            current = queue.popleft()
            if current in visited:
                continue
            parent = self.nodes[current].parent
            if parent is not None:
                thisedge = (parent, current)
            if len(self.nodes[current].children) == 0:
                # Leaf node
                self.nodes[current].dist_crown = 0
                if parent is not None:
                    self.nodes[current].dist_crown_plus_root = self.edges[thisedge].branch_length
                    self.edges[thisedge].dist_crown_plus_this_edge = self.edges[thisedge].branch_length
            else:
                if all([self.nodes[child].dist_crown is not None for child in self.nodes[current].children]):
                    self.nodes[current].dist_crown = sum(self.nodes[child].dist_crown_plus_root
                                                         for child in self.nodes[current].children)
                    if parent is not None:
                        self.nodes[current].dist_crown_plus_root = self.nodes[current].dist_crown + self.edges[thisedge].branch_length
                        self.edges[thisedge].dist_crown_plus_this_edge = self.nodes[current].dist_crown + self.edges[thisedge].branch_length
                else:
                    # Re-add the current node to the queue if not all children are processed
                    queue.append(current)
                    continue
            # Mark the current node as visited
            # This only is activated if the current node is not re-added
            visited.add(current)
            # Add the parent node to the queue for further processing
            if parent is not None and parent not in visited:
                queue.append(parent)
        print()
        # set the values for the root
        self.root = self.find_root()
        self.nodes[self.root].dist_crown = sum([self.nodes[child].dist_crown_plus_root
                                                for child in self.nodes[self.root].children])
        self.nodes[self.root].dist_crown_plus_root = self.nodes[self.root].dist_crown

        # now go through and make sure that all of the nodes have a distance to the crown
        for node in self.nodes:
            if self.nodes[node].dist_crown is None:
                raise ValueError("The node {} does not have a distance to the crown.\n  - {}".format(
                    node, self.nodes[node]))
            if self.nodes[node].dist_crown_plus_root is None:
                raise ValueError("The node {} does not have a distance to the crown plus root.\n  - {}".format(
                    node, self.nodes[node]))

    def add_lineage_info(self) -> None:
        """
        Uses BFS to add lineage information to the nodes.
        """
        root = self.find_root()
        queue = [root]
        while len(queue) > 0:
            current = queue.pop(0)
            if current == root:
                self.nodes[current].lineage = [root]
            else:
                parent = self.nodes[current].parent
                self.nodes[current].lineage = self.nodes[parent].lineage + [current]
            # add the kids to the queue
            for child in self.nodes[current].children:
                queue.append(child)
        # now format the lineages as strings
        for node in self.nodes:
            self.nodes[node].lineage_string = ";".join([str(x) for x in self.nodes[node].lineage])

        # If there are any edges, we can add the lineage information to the edges
        for edge in self.edges:
            self.edges[edge].parent_lineage = self.nodes[edge[0]].lineage
            self.edges[edge].child_lineage  = self.nodes[edge[1]].lineage

    def print_edge_information(self, outfile) -> None:
        """
        Prints the edge information of the tree. Use __slots__ to determine what to print.
        """
        # get the fields
        fields = self.edges[(self.root, list(self.nodes[self.root].children)[0])].__slots__
        with open(outfile, "w") as outhandle:
            print("\t".join(fields), file = outhandle)
            for edge in self.edges:
                print("\t".join([str(getattr(self.edges[edge], x))
                                 for x in self.edges[edge].__slots__]), file = outhandle)

    def print_node_information(self, outfile) -> None:
        """
        prints all the fields of the nodes to a file. Use the slots to determine what to print.
        """
        fields = self.nodes[self.root].__slots__
        with open(outfile, "w") as outhandle:
            print("\t".join(self.nodes[self.root].__slots__), file = outhandle)
            for node in self.nodes:
                print("\t".join([str(getattr(self.nodes[node], x))
                                 for x in self.nodes[self.root].__slots__]), file = outhandle)

def main():
    # first we need to parse the arguments from the comand line
    args = parse_args()
    print(args)

    # We will use taxnames for many things
    NCBI = ete3.NCBITaxa()
    ## now we need to load in the newick file
    tree = read(args.newick)[0]
    entries = []
    # iterate through the leaves
    for thisleaf in tree.get_leaves():
        taxname = thisleaf.name.replace("_", " ")
        taxid   = NCBI.get_name_translator([taxname])[taxname][0]
        entry = {"node": thisleaf,
                 "nodename": thisleaf.name,
                 "taxname": taxname,
                 "taxid": taxid,
                 "lineage": NCBI.get_lineage(taxid)}
        entries.append(entry)

    # the ttdf is the timetree df
    ttdf = pd.DataFrame(entries)
    ttdf_taxid_to_sp = {row["taxid"]: row["nodename"] for i, row in ttdf.iterrows()}
    print(ttdf)

    # make a tttree
    TTtree = TaxIDtree()
    for i, row in ttdf.iterrows():
        if 346063 in row["lineage"]:
            print("The lineage of {} is: {}".format(row["nodename"], row["lineage"]))
            sys.exit()
        for j in range(0, len(row["lineage"])-1):
            TTtree.add_edge(row["lineage"][j], row["lineage"][j+1])

    # get the all vs all time divergence
    # print getting the divergence times
    print("Getting the divergence times for all species.")
    divergence_times = report_divergence_time_all_vs_all(tree, args.prefix)
    # print the first few divergence times
    print("The first few divergence times are:")
    for i, (sp1, sp2) in enumerate(divergence_times):
        if i > 10:
            break
        print("{}\t{}\t{}".format(sp1, sp2, divergence_times[(sp1, sp2)]))

    # If there is an odp config file, we will try to link the species of the config file to
    #   the species in the tree. Timetree.org obviously doesn't have all of the species, so
    #   sometimes we need to just find the closest species. In this case, we will have to
    #   use something to interact with the NCBI taxonomy database to find closest species pairs.
    #   This will take something like all-v-all comparisons of the lineages to find something
    #   that is close.
    # check if the prefix exists in the config file
    # make a df, similar to that of the ttdf, but for the config file
    entries = []
    config = yaml.safe_load(open(args.config, "r"))
    for sp in config["species"]:
        taxid = config["species"][sp]["taxid"]
        # get the taxname using NCBI
        spname = NCBI.get_taxid_translator([taxid])[taxid]
        lineage = NCBI.get_lineage(config["species"][sp]["taxid"])
        # get the closest species in timetree
        entry = {"node": None,
                 "nodename": sp,
                 "taxname": spname,
                 "taxid": taxid,
                 "taxid_current": lineage[-1],
                 "taxid_in_timetree": TTtree.find_closest_relative(NCBI, taxid),
                 "lineage": lineage}
        entries.append(entry)
    configdf = pd.DataFrame(entries)

    # now make a tree of the species from the config file
    CFtree = TaxIDtree()
    for i, row in configdf.iterrows():
        for j in range(0, len(row["lineage"])-1):
            CFtree.add_edge(row["lineage"][j], row["lineage"][j+1])
    # Add the same edges to a networkX graph
    nxCF = nx.DiGraph()
    for i, row in configdf.iterrows():
        for j in range(0, len(row["lineage"])-1):
            nxCF.add_edge(row["lineage"][j], row["lineage"][j+1])
    # verify whether it is a single CC
    if not nx.is_weakly_connected(nxCF):
        raise ValueError("The config file tree is not a single connected component.")
    # set the leaf ages to zero
    CFtree.set_leaf_ages_to_zero()

    outfile = f"{args.prefix}.node_ages_for_config.tsv"
    if not os.path.exists(outfile):
        # if the weights file does not yet exist, we need to make it.
        # This step takes the longest of anything in the paper, so it is better to load if possible.
        # Now go through and add all the ages to the config file tree
        for i in range(len(configdf["taxid"])-1):
            # make a print statement that stays on one line
            print("    - Iterating through species, on sample {}/{}".format(i, len(configdf["taxid"])), end="\r")
            for j in range(i+1, len(configdf["taxid"])):
                taxid1_config = configdf["taxid_current"][i]
                taxid2_config = configdf["taxid_current"][j]
                taxid1_tt     = configdf["taxid_in_timetree"][i]
                taxid2_tt     = configdf["taxid_in_timetree"][j]
                #print("Comparing {} and {}".format(taxid1_tt, taxid2_tt))
                #print("  - config taxids are {} and {}".format(taxid1_config, taxid2_config))
                if taxid1_tt != taxid2_tt:
                    sp1_tt        = ttdf_taxid_to_sp[taxid1_tt]
                    sp2_tt        = ttdf_taxid_to_sp[taxid2_tt]
                    div_time      = divergence_times[tuple(sorted((sp1_tt, sp2_tt)))]
                    common_ancestor = CFtree.find_LCA(taxid1_config, taxid2_config)
                    # add the node age to the node in CFtree
                    CFtree.nodes[common_ancestor].nodeages.update([div_time])

        # print out the node ages for all the nodes
        with open(outfile, "w") as f:
            f.write("TaxID\tNodeAges\n")
            for thisnode in CFtree.nodes:
                f.write("{}\t{}\n".format(thisnode, CFtree.nodes[thisnode].nodeages))
    else:
        # The file is already there, so we just need to load its values into the graph.
        with open(outfile, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                taxid, ages = line.strip().split("\t")
                CFtree.nodes[int(taxid)].nodeages = eval(ages)

    # now we correct the missing nodes
    print("Now we are correcting the missing nodes")
    CFtree.correct_missing_nodes()
    CFtree.ensure_all_leaves_have_age_zero()
    # calculate the dist crown
    CFtree.calc_dist_crown()
    # calculate the lineage info
    CFtree.add_lineage_info()
    # add chromosome information
    CFtree.add_chromosome_info_file(args.chromosome_sizes)

    # Print out the info of the first 5 leaves
    print("These are the first five leaves")
    counter = 0
    for thisnode in CFtree.nodes:
        if len(CFtree.nodes[thisnode].children) == 0:
            print (CFtree.nodes[thisnode])
            counter += 1
        if counter >= 5:
            break
    # print out the first five non-leaves
    print("These are the first five non-leaves")
    counter = 0
    for thisnode in CFtree.nodes:
        if len(CFtree.nodes[thisnode].children) > 0:
            print (CFtree.nodes[thisnode])
            counter += 1
        if counter >= 5:
            break
    # now print how many nodes are missing ages
    missing_nodes = [node for node in CFtree.nodes if len(CFtree.nodes[node].nodeages) == 0]
    print("There are {} nodes missing ages.".format(len(missing_nodes)))

    # now we print the edge information
    CFtree.print_edge_information(f"{args.prefix}.edge_information.tsv")
    CFtree.print_node_information(f"{args.prefix}.node_information.tsv")
    print("success")

if __name__ == "__main__":
    main()