#!/usr/bin/env python

"""
This script takes a newick tree and identifies the divergence time of various nodes in the tree.
"""

import argparse
from Bio import Entrez
import os
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

def download_all_taxinfo(config, output_prefix, email):
    """
    This controls a loop that handles downloading all of the taxinfo
    for all of the species in the config file. Returns a dict of the taxinfo when
    done, and a path to the file where the taxinfo is stored.
    """
    # set up email for Entrez
    Entrez.email = email
    # for now we just want to investigate the best way to get the lineage information
    taxinfo_yaml = {"taxinfo": {}}
    taxinfo_yaml_filepath = "{}.taxid.yaml".format(output_prefix)

    # We will make a temporary folder to store the taxid information for each species.
    #  Maybe there are files there already, so don't overwrite them.
    tempdir = "{}.taxid.temp".format(output_prefix)
    create_directories_recursive_notouch(tempdir)

    species_remaining = set(config["species"].keys())
    species_completed_this_round = set()
    downloading_round = 0
    # now we need to loop through all of the species in the config file
    print("DOWNLOADING ROUND {}".format(downloading_round), file = sys.stderr)
    while len(species_remaining) > 0:
        for thissp in config["species"]:
            binomial = "{} {}".format(config["species"][thissp]["genus"], config["species"][thissp]["species"])
            sp_remaining = len(species_remaining) - len(species_completed_this_round)
            print("downloading", binomial, "- {} species remaining".format(sp_remaining))
            taxinfo_filepath = os.path.join(tempdir, "{}.taxinfo.yaml".format(thissp))
            success_value = taxinfo_download_or_load(binomial, taxinfo_filepath)
            if success_value == 0:
                species_completed_this_round.add(thissp)
        species_remaining = species_remaining - species_completed_this_round
        species_completed_this_round = set()
    # Now that we know that all the yaml files exist, just run that round of checks again
    for thissp in config["species"]:
        taxinfo_filepath = os.path.join(tempdir, "{}.taxinfo.yaml".format(thissp))
        if not yaml_file_legal(taxinfo_filepath):
            # There is a problem with this file, so we should make the game through an error
            raise ValueError("The yaml file for {} is not legal.".format(thissp))
        # if the file is legal, read it in and add it to the dict 
        with open(taxinfo_filepath, 'r') as file:
            thissptax = yaml.load(file, Loader=yaml.Loader )
            taxinfo_yaml["taxinfo"][thissp] = thissptax

    # save all of the taxinfo to a yaml file
    # Fine to overwrite
    create_directories_recursive_notouch(taxinfo_yaml_filepath)
    with open(taxinfo_yaml_filepath, "w") as f:
        yaml.dump(taxinfo_yaml, f)

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
        with open(args.config, 'r') as file:
            config = yaml.safe_load(file)
        
        # we may have already saved a yaml file with the taxonomy information, so check for that
        taxinfo_filepath = "{}.taxinfo.yaml".format(args.prefix)
        # safely make the directories if they don't exist
        create_directories_recursive_notouch(taxinfo_filepath)
        taxinfo_yaml = {}
        # open the taxinfo file for writing if it doesn't exist
        if os.path.exists(taxinfo_filepath):
            with open(taxinfo_filepath, "r") as f:
                taxinfo_yaml = yaml.safe_load(f)
        # if the file doesn't exist yet we have to parse the info from NCBI
        else:
            download_all_taxinfo(config, args.prefix, args.email)
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