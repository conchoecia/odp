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

def download_all_taxinfo(sp_binomials, output_prefix, email):
    """
    Inputs:
      - sp_binomials: a dict of keys for the species names
        - the keys are the special string for that species
        - the values are the binomial "Genus species" string
      - The output prefix is what the files will be saved as
      - the email is the email address to use for programmatic access to NCBI

    This controls a loop that handles downloading all of the taxinfo
    for all of the species in the sp_binomials object. Returns a dict of the taxinfo when
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

    species_remaining = set(sp_binomials.keys())
    species_completed_this_round = set()
    downloading_round = 0
    # now we need to loop through all of the species in the sp_binomials_object
    print("DOWNLOADING ROUND {}".format(downloading_round), file = sys.stderr)
    while len(species_remaining) > 0:
        for thissp in sp_binomials.keys():
            binomial = sp_binomials[thissp]
            sp_remaining = len(species_remaining) - len(species_completed_this_round)
            print("downloading", binomial, "- {} species remaining".format(sp_remaining))
            taxinfo_filepath = os.path.join(tempdir, "{}.taxinfo.yaml".format(thissp))
            success_value = taxinfo_download_or_load(binomial, taxinfo_filepath)
            if success_value == 0:
                species_completed_this_round.add(thissp)
        species_remaining = species_remaining - species_completed_this_round
        species_completed_this_round = set()
    # Now that we know that all the yaml files exist, just run that round of checks again
    for thissp in sp_binomials.keys():
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
    
    return taxinfo_yaml, taxinfo_yaml_filepath

def download_config_data(config_filepath, prefix, email):
    """
    Handles the download of the taxinfo in the config file.

    Returns a yaml file of the taxonomy info,
      and the filepath that also has all of the taxonomy information.
    """
    with open(config_filepath, 'r') as file:
        config = yaml.safe_load(file)
    
    # we may have already saved a yaml file with the taxonomy information, so check for that
    taxinfo_filepath = "{}.taxinfo.yaml".format(prefix)
    # safely make the directories if they don't exist
    create_directories_recursive_notouch(taxinfo_filepath)
    taxinfo_yaml = {}
    # open the taxinfo file for writing if it doesn't exist
    if os.path.exists(taxinfo_filepath):
        with open(taxinfo_filepath, "r") as f:
            taxinfo_yaml = yaml.safe_load(f)
    # if the file doesn't exist yet we have to parse the info from NCBI
    else:
        binomial_dict = {k: "{} {}".format(config["species"][k]["genus"],
                                                config["species"][k]["species"]) \
                                                    for k in config["species"]}

        taxinfo_yaml, taxinfo_filepath = download_all_taxinfo(binomial_dict, prefix, email)
    return taxinfo_yaml, taxinfo_filepath

class TaxNode:
    """
    one node of the taxonomy tree
    """
    def __init__(self, taxid, name, rank) -> None:
        self.taxid = taxid
        self.name = name
        self.rank = rank
        self.children = {}

class TaxIDtree:
    """
    This is a datastructure to quickly search for the most closely related species
      given search species 1 and a tree of species 2...N.
    """
    def __init__(self, taxinfoyaml = None) -> None:
        self.tree = TaxNode(taxid = -1, name = "root", rank = "root")
        if taxinfoyaml:
            self.taxinfoyaml = taxinfoyaml["taxinfo"] 
            self.build_tree_from_yaml(self.taxinfoyaml)

    def __str__(self) -> str:
        newoutstring = "- root\n"
        # now make a recursive algorithm to print the whole tree
        def print_tree(node, outstring, level):
            for child in node.children:
                outstring += "{}|_ {}: {}\n".format("  "*level, child, node.children[child].name)
                outstring = print_tree(node.children[child], outstring, level+1)
            return outstring
        outstring = print_tree(self.tree, newoutstring, 1)
        return outstring

    def build_tree_from_yaml(self, taxinfoyaml):
        """
        Takes a yaml file of the taxinfo and builds a tree from it.
        """
        # go through all of the entries in the yaml
        for entry in taxinfoyaml:
            self.add_to_tree_from_yaml_entry(taxinfoyaml[entry])

    def add_to_tree_from_yaml_entry(self, entry):
        """
        Takes a yaml entry and builds a tree from it.

        The entries look like this:
        {'Lineage': 'cellular organisms; Eukaryota; Opisthokonta; Metazoa; Porifera; Hexactinellida; Hexasterophora; Lyssacinosida; Leucopsacidae; Oopsacas',
         'LineageEx': [{'Rank': 'no rank', 'ScientificName': 'cellular organisms', 'TaxID': 131567}, {'Rank': 'superkingdom', 'ScientificName': 'Eukaryota', 'TaxID': 2759}, {'Rank': 'clade', 'ScientificName': 'Opisthokonta', 'TaxID': 33154}, {'Rank': 'kingdom', 'ScientificName': 'Metazoa', 'TaxID': 33208}, {'Rank': 'phylum', 'ScientificName': 'Porifera', 'TaxID': 6040}, {'Rank': 'class', 'ScientificName': 'Hexactinellida', 'TaxID': 60882}, {'Rank': 'subclass', 'ScientificName': 'Hexasterophora', 'TaxID': 60883}, {'Rank': 'order', 'ScientificName': 'Lyssacinosida', 'TaxID': 60884}, {'Rank': 'family', 'ScientificName': 'Leucopsacidae', 'TaxID': 472148}, {'Rank': 'genus', 'ScientificName': 'Oopsacas', 'TaxID': 111877}],
         'ScientificName': 'Oopsacas minuta', 'TaxID': 111878
         }
        """
        # Get the longform taxonomic info of this entry
        sciname = entry["ScientificName"]
        taxid   = entry["TaxID"]
        lineage = entry["LineageEx"]
        # now we iterate through the longform lineage to add to the tree
        # start ancestor at -1 for the root
        ancestor = self.tree
        for i in range(len(lineage)):
            thistaxid = lineage[i]["TaxID"]
            thisname = None
            if "ScientificName" in lineage[i]:
                thisname  = lineage[i]["ScientificName"]
            thisrank =  "None"
            if "Rank" in lineage[i]:
                thisrank = lineage[i]["Rank"]
            # now we add this to the tree
            # first we check if this is in the tree
            if thistaxid in ancestor.children:
                # this taxid is already in the tree, so we just need to move on
                ancestor = ancestor.children[thistaxid]
            else:
                # this taxid isn't in the tree, so we need to add it.
                ancestor = self.add_to_tree(thistaxid, thisname, thisrank, ancestor)
        # now we add the final entry to the tree
        ancestor = self.add_to_tree(taxid, sciname, "species", ancestor)

    def add_to_tree(self, taxid, name, rank, ancestor) -> TaxNode:
        """
        Takes a taxid, name, and lineage and adds it to the tree.
        """
        if taxid in ancestor.children:
            raise IOError("This taxid is already in the tree: {}".format(taxid)) 
        ancestor.children[taxid] = TaxNode(taxid, name, rank)
        # return the new node
        return ancestor.children[taxid]

    def find_closest_relative(self, yaml_entry) -> TaxNode:
        """
        Uses a taxonomy entry to get the closest lineage. Returns a TaxNode object

        {'Lineage': 'cellular organisms; Eukaryota; Opisthokonta; Metazoa; Porifera; Hexactinellida; Hexasterophora; Lyssacinosida; Leucopsacidae; Oopsacas',
         'LineageEx': [{'Rank': 'no rank', 'ScientificName': 'cellular organisms', 'TaxID': 131567}, {'Rank': 'superkingdom', 'ScientificName': 'Eukaryota', 'TaxID': 2759}, {'Rank': 'clade', 'ScientificName': 'Opisthokonta', 'TaxID': 33154}, {'Rank': 'kingdom', 'ScientificName': 'Metazoa', 'TaxID': 33208}, {'Rank': 'phylum', 'ScientificName': 'Porifera', 'TaxID': 6040}, {'Rank': 'class', 'ScientificName': 'Hexactinellida', 'TaxID': 60882}, {'Rank': 'subclass', 'ScientificName': 'Hexasterophora', 'TaxID': 60883}, {'Rank': 'order', 'ScientificName': 'Lyssacinosida', 'TaxID': 60884}, {'Rank': 'family', 'ScientificName': 'Leucopsacidae', 'TaxID': 472148}, {'Rank': 'genus', 'ScientificName': 'Oopsacas', 'TaxID': 111877}],
         'ScientificName': 'Oopsacas minuta', 'TaxID': 111878
         }
        """
        # first we need to get the lineage of this entry
        lineage = yaml_entry["LineageEx"]
        # now we need to iterate through the lineage to find the closest relative
        # start ancestor at -1 for the root
        ancestor = self.tree
        for i in range(len(lineage)):
            thistaxid = lineage[i]["TaxID"]
            # now we check if this taxid is in the tree
            if thistaxid in ancestor.children:
                # this taxid is already in the tree, so keep going down the tree
                ancestor = ancestor.children[thistaxid]
            else:
                # We have reached a point where this taxid is not in the tree.
                break
        # We check if the ancestor node has children.
        # If it does, that means we haven't reached the target species or a closely related species yet.
        done = False
        if len(ancestor.children) == 0:
            done = True
        else:
            # we need to find the closest relative in the tree
            while done == False:
                # greedily take the first child until we get to a leaf
                # this is a closest relative
                if len(ancestor.children) == 0:
                    done = True
                else:
                    ancestor = ancestor.children[list(ancestor.children.keys())[0]]
        return ancestor

def main():
    # first we need to parse the arguments from the comand line
    args = parse_args()
    print(args)

    ## now we need to load in the newick file
    tree = read(args.newick)[0]

    # get the all vs all time divergence
    divergence_times = report_divergence_time_all_vs_all(tree, args.prefix)

    # If there is an odp config file, we will try to link the species of the config file to
    #   the species in the tree. Timetree.org obviously doesn't have all of the species, so
    #   sometimes we need to just find the closest species. In this case, we will have to
    #   use something to interact with the NCBI taxonomy database to find closest species pairs.
    #   This will take something like all-v-all comparisons of the lineages to find something
    #   that is close.
    # check if the prefix exists in the config file
    if "config" in args:
        # download the taxinfo from the config file
        taxinfo_yaml, taxinfo_filepath = download_config_data(args.config,
                                                              args.prefix,
                                                              args.email)
        # now we download the taxinfo for the species in the tree
        # first we need to get all of the leaves
        leaves = tree.get_leaves()
        print(dir(leaves[0]))
        binomial_dict = {k.name: "{} {}".format(k.name.split("_")[0],
                                                k.name.split("_")[1]
                                                ) \
                                                for k in leaves}
        treetax_yaml, treetax_filepath = download_all_taxinfo(binomial_dict,
                                                              args.prefix + ".tree",
                                                              args.email)

        # For each species in the config file, we need to find the closest species in the tree
        # We only need to map config -> tree for these reasons:
        #  - The point of the program is to get a list of genomes from the config file.
        #  - That list of genomes is fed to timetree.org
        #  - Timetree.org doesn't have all species, so it finds the closest species.
        #  - The tree we get from Timetree.org has some species that are not in the config file.
        #  - We need to go back and figure out what species in the config file are closest to those in the tree.
        TreeTaxStruct = TaxIDtree(taxinfoyaml = treetax_yaml)

        # go through each species in the config file and find the closest species from the timetree
        closest_relative_dict = {}
        for sp in taxinfo_yaml["taxinfo"]:
            # get the closest relative to this species
            closest_relative = TreeTaxStruct.find_closest_relative(taxinfo_yaml["taxinfo"][sp])
            closest_relative_dict[sp] = closest_relative.name
            #print("{}: {}".format(sp, closest_relative.name))

        # Now construct the divergence_times_config dict with the species translated 
        divergence_times_config = {"divergence_times": {}}
        sp_names = list(sorted(closest_relative_dict.keys()))
        for i in range(len(sp_names)-1):
            for j in range(i+1, len(sp_names)):
                sp1 = sp_names[i]
                sp2 = sp_names[j]
                tree_sp1 = "_".join(closest_relative_dict[sp1].split(" "))
                tree_sp2 = "_".join(closest_relative_dict[sp2].split(" "))
                if tree_sp1 == tree_sp2:
                    age = 0
                else:
                    age = divergence_times[tuple(sorted((tree_sp1, tree_sp2)))]
                if sp1 not in divergence_times_config["divergence_times"]:
                    divergence_times_config["divergence_times"][sp1] = {}
                if sp2 not in divergence_times_config["divergence_times"][sp1]:
                    divergence_times_config["divergence_times"][sp1][sp2] = age
                #print("{}\t{}\t{}".format(sp1, sp2, divergence_times_config[(sp1, sp2)]))

        # make a new config file where we will write the divergence times
        outconfig = args.prefix + ".divTimeAdded.yaml"
        with open(outconfig, "w") as f:
            # copy the yaml from args to this new file
            with open(args.config, "r") as f2:
                config = yaml.safe_load(f2)
                yaml.dump(config, f)
            yaml.dump(divergence_times_config, f)

if __name__ == "__main__":
    main()