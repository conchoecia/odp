#!/usr/bin/env python

"""
This script takes a newick tree and identifies the divergence time of various nodes in the tree.
"""

import argparse
import os
import sys
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
    file_endings = [".txt", ".tsv", ".csv"]
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
      - the email address to use for programmatic access to NCBI
    """
    parser = argparse.ArgumentParser(description="This script takes a newick tree and identifies the divergence time of various nodes in the tree.")
    parser.add_argument("-n", "--newick", help="The path to the newick tree file.", required=True)
    parser.add_argument("-p", "--prefix", help="The output prefix for the file, including a prepended path if you want another directory.", required=True)
    parser.add_argument("-c", "--config", help="The original config file from which the tree was derived.", required=False)
    parser.add_argument("-e", "--email", help="The email address to use for programmatic access to NCBI.", required=False)

    # check that the newick and config files actually exist
    args = parser.parse_args()
    if not os.path.exists(args.newick):
        # raise an IO error
        raise IOError("The newick file does not exist: {}".format(args.newick))
    
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
    sp1_age = sum([x.length for x in unique_sp1] ) #+ [common_ancestor.length]) 
    sp2_age = sum([x.length for x in unique_sp2] ) #+ [common_ancestor.length])
    # the ages should be the same, so check
    # Sometimes when one of the species has a really long branch, the ages are not exactly the same.
    # Just check that they are within 0.05 of each other in terms of precent.
    percent_diff = 0 if abs(sp1_age - sp2_age) == 0 else (abs(sp1_age - sp2_age)/sp1_age)

    # There is a weird behavior where, if percent_diff is 0, then the equality statement doesn't work as predicted. 
    #  So we need to handle that case separately
    sp_1_2_within_0_0_5 = True if (percent_diff == 0) else (percent_diff < (0.05 * sp1_age))
    if not sp_1_2_within_0_0_5:
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
            # round this to 5 decimal places
            age_print = round(age, 5)
            yield sp1, sp2, age_print

def report_divergence_time_all_vs_all(tree, output_prefix):
    """
    This method gets the divergence times and writes them to a file with the prefix.
    Also safely makes a directory if it does not yet exist.
    """
    # first come up with the outfile path
    outfile_path = "{}.divergence_times.txt".format(output_prefix)
    # safely make the directories if they don't exist
    create_directories_recursive_notouch(outfile_path)
    # open the outfile for writing
    with open(outfile_path, "w") as f:
        for sp1, sp2, age in get_divergence_time_all_vs_all(tree):
            f.write("{}\t{}\t{}\n".format(sp1, sp2, age))


# Replace with your own email address
Entrez.email = "your.email@example.com"

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

    return record[0]

def main():
    # first we need to parse the arguments from the comand line
    args = parse_args()
    print(args)

    ## now we need to load in the newick file   
    tree = read(args.newick)[0]

    ## get the all vs all time divergence
    #report_divergence_time_all_vs_all(tree, args.prefix)

    # If there is an odp config file, we will try to link the species of the config file to
    #   the species in the tree. Timetree.org obviously doesn't have all of the species, so
    #   sometimes we need to just find the closest species. In this case, we will have to
    #   use something to interact with the NCBI taxonomy database to find closest species pairs.
    #   This will take something like all-v-all comparisons of the lineages to find something
    #   that is close.
    # check if the prefix exists in the config file
    if "config" in args:
        import yaml
        from Bio import Entrez
        with open(args.config, 'r') as file:
            config = yaml.safe_load(file)

        print(config)

    # This is all debug code
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
    
if __name__ == "__main__":
    main()