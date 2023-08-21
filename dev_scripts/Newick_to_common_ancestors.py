#!/usr/bin/env python

"""
This script takes a newick tree and identifies the divergence time of various nodes in the tree.
"""

import argparse
import os
import sys
# we need the special newick package to parse this file type
from newick import read

def parse_args():
    """
    The args that we need to parse are:
      - the file path of the newick tree
      - the species of interest  
    """
    parser = argparse.ArgumentParser(description="This script takes a newick tree and identifies the divergence time of various nodes in the tree.")
    parser.add_argument("-n", "--newick", help="The path to the newick tree file", required=True)
    parser.add_argument("-s", "--species", help="The species of interest", required=True)

    # check that the newick file actually exists
    args = parser.parse_args()
    if not os.path.exists(args.newick):
        # raise an IO error
        raise IOError("The newick file does not exist: {}".format(args.newick))
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
    print(leaves)
    # now we need to get the lineage path for each of these leaves
    lineages = {leaf.name: get_lineage(tree, leaf) for leaf in leaves}
    return lineages

def find_common_ancestor_age(sp1_lineage, sp2_lineage):
    """
    Takes two lineages and finds the common ancestor
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
    sp1_age = sum([x.length for x in unique_sp1]) 
    sp2_age = sum([x.length for x in unique_sp2])
    # the ages should be the same, so check
    # Sometimes when one of the species has a really long branch, the ages are not exactly the same.
    # Just check that they are within 0.05 of each other in terms of precent.
    sp_1_2_within_0_1 = (abs(sp1_age - sp2_age)/sp1_age) < (0.05 * sp1_age)
    if not sp_1_2_within_0_1:
        raise ValueError("The ages of the two species are not the same: {} vs {}".format(sp1_age, sp2_age))
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
            print("{}\t{}\t{}".format(sp1, sp2, age))

def main():
    # first we need to parse the arguments from the comand line
    args = parse_args()

    ## now we need to load in the newick file   
    tree = read(args.newick)[0]
    print(tree.ascii_art())
    print(dir(tree))

    # get the all vs all time divergence
    get_divergence_time_all_vs_all(tree)
    #print(dir(tree))
    #print("descendants", tree.descendants)
    #print("capsaspora_name", tree.descendants[0].name)
    #print("capsaspora_length", tree.descendants[0].length)
    #print("tree", tree)
    #print("comment", tree.comment)
    #print("length", tree.length)
    #print("root name", tree.name)
    #print(tree.get_leaves())
    #print(tree.get_leaves()[-1].ancestor)
    #print(tree.get_leaves()[0].ancestor)
    #print(tree.get_leaves()[0].ancestor == tree)
    #print(get_lineage(tree, tree.get_leaves()[4]))





    ## get all of the node labels
    #node_labels = [node.name for node in tree.get_leaves()]

    ## make sure that the species of interest is in the nodeL-labels
    #if args.species not in node_labels:
    #    raise ValueError("The species of interest is not in the tree: {}".format(args.species))

    ## now we get the node age of each species combination
    #for sp2 in node_labels:
    #    if sp2 != args.species:
    #        common_ancestor = find_common_ancestor(tree, [args.species, sp2])
    #        print("{}\t{}\t{}".format(args.species, sp2, common_ancestor.dist))
    
if __name__ == "__main__":
    main()