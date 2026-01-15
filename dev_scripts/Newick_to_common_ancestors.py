#!/usr/bin/env python
"""
Newick Tree to Common Ancestors - Divergence Time Analysis

PURPOSE:
========
This script combines two phylogenetic trees to produce a fully-annotated phylogeny:
1. A custom topology tree (from taxids_to_newick.py) with your desired tree structure and all species
2. A time-calibrated tree (from TimeTree.org) with divergence times for a subset of species

The result is a phylogeny that uses YOUR topology while incorporating TimeTree divergence times
where available and interpolating missing times. This enables custom phylogenetic hypotheses
(like Ctenophora as sister to all animals) while leveraging TimeTree's molecular clock calibration.

USAGE:
======
Standard usage:
    python Newick_to_common_ancestors.py --time_newick timetree.nwk --topology_newick custom_topology.nwk -p output_prefix -c config.yaml

With chromosome sizes:
    python Newick_to_common_ancestors.py --time_newick timetree.nwk --topology_newick custom_topology.nwk -p output_prefix -c config.yaml -C chromosome_sizes.tsv

Generate species list for TimeTree.org:
    python Newick_to_common_ancestors.py --time_newick dummy.nwk --topology_newick dummy.nwk -p dummy -c config.yaml -s

INPUT:
======
- Time newick file (.nwk) - Time-calibrated tree from TimeTree.org with branch lengths in millions of years (subset of species)
- Topology newick (.nwk) - Full topology tree from taxids_to_newick.py with all species (defines tree structure)
- Config file (YAML) - Contains species information including taxids, genus, and species names
- Chromosome sizes file (optional) - TSV with sample names and chromosome size information

Note: Species names in newick files should match "Genus_species" format (underscores, not spaces).
      Taxids are obtained from ete4 NCBITaxa database using species names - no network queries needed.
      The topology_newick defines the tree structure; time_newick provides calibration for subset of species.

OUTPUT:
=======
- {prefix}.divergence_times.txt - Pairwise divergence times for all species (TSV format)
- {prefix}.edge_information.tsv - Information about each edge/branch in the phylogeny
- {prefix}.node_information.tsv - Information about each node including ages and lineages
- species_list.txt - List of binomial names for TimeTree.org upload (if -s flag used)

WORKFLOW:
=========
Complete workflow for custom phylogenies (e.g., Ctenophora placement):

1. Generate custom topology tree:
   python taxids_to_newick.py -c config.yaml -o custom_topology.nwk --custom_phylogeny

2. Export species list for TimeTree.org:
   python taxids_to_newick.py -c config.yaml -o custom_topology.nwk --custom_phylogeny --timetree_list species_for_timetree.txt

3. Upload to TimeTree.org:
   - Go to https://timetree.org
   - Upload species_for_timetree.txt
   - Download calibrated tree (subset of species with divergence times)
   - Save as timetree_calibrated.nwk

4. Run this script with both trees:
   python Newick_to_common_ancestors.py --time_newick timetree_calibrated.nwk --topology_newick custom_topology.nwk -p output -c config.yaml

5. Result:
   - Uses YOUR custom topology (all species, custom Ctenophora placement)
   - Annotated with TimeTree divergence times where available
   - Missing times interpolated based on tree structure
   - Outputs: divergence times, node ages, edge information

NOTES:
======
- Branch lengths in input Newick must be in millions of years (MYA)
- Script handles species name matching between config and TimeTree output
- Uses ete4 NCBITaxa local database to map species names to taxids (no network access required)
- Interpolates missing node ages when necessary
- Ensure ete4 NCBI taxonomy database is initialized before running (one-time setup)

REQUIREMENTS:
=============
- ete4 (updated from ete3)
- newick - Newick format parser
- pandas, numpy, matplotlib, networkx
- PyYAML

AUTHOR: Darrin T. Schultz
DATE: 2023-2025
"""

import argparse
from collections import Counter,deque
from ete4 import NCBITaxa
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
import numpy as np
import os
import pandas as pd
import random
import sys
import yaml
import time
# we need the PhyloTree parser from ete4 to parse newick files
from ete4 import PhyloTree
# Import subspecies conversion functions from taxids_to_newick.py
from taxids_to_newick import is_subspecies_or_below, get_species_level_taxid

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
      - the file path of the time-calibrated newick tree (from TimeTree.org)
      - the file path of the topology newick tree (optional, from taxids_to_newick.py)
      - the output prefix for the file, including a prepended path
      - the original config file from which the tree was derived
      - the taxids of interest, if we want to plot the divergence relative to a specific species. Accepts a list of taxids
      - file with samples and chromosome sizes
      - a flag -s that just tells the program to make the species list
    """
    parser = argparse.ArgumentParser(description="This script takes a newick tree and identifies the divergence time of various nodes in the tree.")
    parser.add_argument("-n", "--time_newick", "--newick", dest="time_newick", help="The path to the TimeTree.org calibrated newick tree file (with divergence times in MYA).", required=True)
    parser.add_argument("--topology_newick", help="The path to the custom topology newick file from taxids_to_newick.py (with all species). This topology will be used and annotated with times from --time_newick.", required=True)
    parser.add_argument("-C", "--chromosome_sizes", help="The path to the chromosome sizes")
    parser.add_argument("-p", "--prefix", help="The output prefix for the file, including a prepended path if you want another directory.", required=True)
    parser.add_argument("-c", "--config", help="The original config file from which the tree was derived.", required=False)
    parser.add_argument("-t", "--taxids", help="The taxids of interest, if we want to plot the divergence relative to a specific species. Accepts a list of taxids.", required=False)
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

    # check that the newick files actually exist
    if not os.path.exists(args.time_newick):
        raise IOError("The time newick file does not exist: {}".format(args.time_newick))
    
    if not os.path.exists(args.topology_newick):
        raise IOError("The topology newick file does not exist: {}".format(args.topology_newick))

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
    ancestor = node.up
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
    leaves = list(tree.leaves())
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
    sp1_age = sum([x.dist for x in unique_sp1] ) #+ [common_ancestor.dist])
    sp2_age = sum([x.dist for x in unique_sp2] ) #+ [common_ancestor.dist])
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

def annotate_custom_tree_with_timetree_ages(custom_tree, timetree_tree, ncbi):
    """
    Transfers divergence time information from TimeTree.org calibrated tree to custom topology tree.
    
    This function:
    1. Extracts pairwise divergence times from the TimeTree.org tree (fewer species, has times)
    2. Maps those times to equivalent nodes in the custom topology tree (all species, no times)
    3. Uses NCBI taxids to match species between trees
    
    Parameters:
    -----------
    custom_tree : newick.Tree
        The custom topology tree from taxids_to_newick.py (all species, no branch lengths)
    timetree_tree : newick.Tree
        The TimeTree.org calibrated tree (subset of species, has branch lengths in MYA)
    ncbi : NCBITaxa
        NCBI taxonomy database object for species name/taxid mapping
        
    Returns:
    --------
    dict : Mapping of (taxid1, taxid2) -> divergence_time_mya for species in custom tree
    
    Notes:
    ------
    - Species in custom tree but not in TimeTree will have interpolated ages
    - Uses get_lineage() to find common ancestors and calculate divergence times
    """
    print("\nAnnotating custom topology tree with TimeTree divergence times...")
    
    # Get divergence times from TimeTree
    timetree_lineages = get_all_lineages(timetree_tree)
    timetree_divergences = {}
    
    # Build map of species names to taxids for TimeTree species
    timetree_name_to_taxid = {}
    for leaf_name in timetree_lineages.keys():
        # Skip leaves with None names
        if leaf_name is None:
            print(f"  Warning: Skipping TimeTree leaf with None name")
            continue
        # Parse taxid from node name if it's in brackets, e.g., "Homo_sapiens[9606]"
        if '[' in leaf_name and leaf_name.endswith(']'):
            taxid = int(leaf_name[leaf_name.rfind('[')+1:-1])
            timetree_name_to_taxid[leaf_name] = taxid
        else:
            # Fallback to old method if no brackets found
            taxname = leaf_name.replace("_", " ")
            try:
                taxid = ncbi.get_name_translator([taxname])[taxname][0]
                timetree_name_to_taxid[leaf_name] = taxid
            except:
                print(f"  Warning: Could not find taxid for TimeTree species: {leaf_name}")
                continue
    
    print(f"  TimeTree has {len(timetree_name_to_taxid)} species with known taxids")
    
    # Calculate all pairwise divergence times from TimeTree
    timetree_sp_list = list(timetree_lineages.keys())
    for i in range(len(timetree_sp_list)-1):
        for j in range(i+1, len(timetree_sp_list)):
            sp1 = timetree_sp_list[i]
            sp2 = timetree_sp_list[j]
            if sp1 not in timetree_name_to_taxid or sp2 not in timetree_name_to_taxid:
                continue
            
            common_ancestor, age = find_common_ancestor_age(timetree_lineages[sp1], timetree_lineages[sp2])
            taxid1 = timetree_name_to_taxid[sp1]
            taxid2 = timetree_name_to_taxid[sp2]
            # Store both directions for easy lookup
            timetree_divergences[(taxid1, taxid2)] = age
            timetree_divergences[(taxid2, taxid1)] = age
    
    print(f"  Extracted {len(timetree_divergences)//2} pairwise divergence times from TimeTree")
    
    # Get all species from custom tree
    custom_leaves = list(custom_tree.leaves())
    custom_name_to_taxid = {}
    
    for leaf in custom_leaves:
        # Skip leaves with None names
        if leaf.name is None:
            print(f"  Warning: Skipping custom tree leaf with None name")
            continue
        # Parse taxid from node name if it's in brackets, e.g., "Homo_sapiens[9606]"
        if '[' in leaf.name and leaf.name.endswith(']'):
            taxid = int(leaf.name[leaf.name.rfind('[')+1:-1])
            custom_name_to_taxid[leaf.name] = taxid
        else:
            # Fallback to old method if no brackets found
            taxname = leaf.name.replace("_", " ")
            try:
                taxid = ncbi.get_name_translator([taxname])[taxname][0]
                custom_name_to_taxid[leaf.name] = taxid
            except:
                print(f"  Warning: Could not find taxid for custom tree species: {leaf.name}")
                continue
    
    print(f"  Custom tree has {len(custom_name_to_taxid)} species with known taxids")
    
    # Now we return the TimeTree divergences - they'll be used to annotate nodes
    # The TaxIDtree structure will handle interpolation for missing species
    return timetree_divergences, timetree_name_to_taxid, custom_name_to_taxid

def extract_timetree_root_age_as_metazoa(timetree):
    """
    Extract the root age from a TimeTree.org calibrated newick tree.
    
    TimeTree newick files have branch lengths in millions of years (MYA) but do NOT
    have internal node labels like "Metazoa[33208]". The root node represents the
    most ancient divergence in the tree (e.g., Metazoa for metazoan-only trees).
    
    To get the root age, we walk from any leaf to the root, summing branch lengths.
    This gives us the divergence time from present to the root.
    
    Args:
        timetree: PhyloTree object from TimeTree.org newick file
        
    Returns:
        float: Age of the root node in millions of years (MYA)
    """
    # Get any leaf from the tree
    leaves = list(timetree.leaves())
    if not leaves:
        print("  Warning: TimeTree has no leaves, cannot extract root age")
        return None
    
    # Walk from leaf to root, summing branch lengths
    current_node = leaves[0]
    root_age = 0.0
    
    while current_node.up:  # While we haven't reached the root
        root_age += current_node.dist
        current_node = current_node.up
    
    print(f"  TimeTree root age (Metazoa): {root_age:.2f} MYA")
    return root_age

def get_divergence_time_all_vs_all_taxidtree(tree):
    """
    Takes a TaxIDtree and gets the divergence times for all species pairs.
    Returns this information as tuples of (taxid1, taxid2, divergence_time).
    """
    # Get all leaf nodes (species)
    leaf_taxids = [taxid for taxid in tree.nodes if len(tree.nodes[taxid].children) == 0]
    leaf_taxids_sorted = sorted(leaf_taxids)
    
    print(f"  Calculating divergence times for {len(leaf_taxids)} species ({len(leaf_taxids)*(len(leaf_taxids)-1)//2} pairs)...")
    
    # For each pair of species, find their common ancestor and calculate divergence time
    for i in range(len(leaf_taxids_sorted)-1):
        for j in range(i+1, len(leaf_taxids_sorted)):
            taxid1 = leaf_taxids_sorted[i]
            taxid2 = leaf_taxids_sorted[j]
            
            # Find the least common ancestor
            lca_taxid = tree.find_LCA(taxid1, taxid2)
            
            # The divergence time is the age of the LCA node
            if lca_taxid in tree.nodes and tree.nodes[lca_taxid].nodeage is not None:
                divergence_time = tree.nodes[lca_taxid].nodeage
                # Round to 5 decimal places
                divergence_time = round(divergence_time, 5)
                yield taxid1, taxid2, divergence_time
            else:
                print(f"  Warning: Could not find divergence time for {taxid1} and {taxid2} (LCA: {lca_taxid})")

def report_divergence_time_all_vs_all(tree, output_prefix):
    """
    This method gets the divergence times and writes them to a file with the prefix.
    Works with TaxIDtree objects.
    
    Output format: taxid1<TAB>taxid2<TAB>divergence_time_mya
    """
    # first come up with the outfile path
    outfile_path = "{}.divergence_times.txt".format(output_prefix)
    # data structure to save the divergence times
    divergence_times = {}
    # safely make the directories if they don't exist
    create_directories_recursive_notouch(outfile_path)
    
    # open the outfile for writing
    num_pairs = 0
    with open(outfile_path, "w") as f:
        for taxid1, taxid2, age in get_divergence_time_all_vs_all_taxidtree(tree):
            entry = (taxid1, taxid2)
            if entry not in divergence_times:
                divergence_times[entry] = age
            f.write("{}\t{}\t{}\n".format(taxid1, taxid2, age))
            num_pairs += 1
    
    print(f"  Wrote {num_pairs} pairwise divergence times to {outfile_path}")
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
                 'nodeageinterpolated', 'lock_age', 'lineage', 'lineage_string',
                 'sort_order', 'x', 'y',
                 'dist_crown', 'dist_crown_plus_root',
                 "chromsize_median", "chromsize_mean",
                 "num_genomes", "chromsize_list",
                 "fusions_in_this_clade",
                 "losses_in_this_clade",
                 "extinction_fusion_spearman_r",
                 "extinction_fusion_spearman_p",
                 "extinction_fusion_spearman_n",
                 "extinction_fusion_kendalltau_r",
                 "extinction_fusion_kendalltau_p",
                 "extinction_fusion_kendalltau_n",
                 "extinction_losses_spearman_r",
                 "extinction_losses_spearman_p",
                 "extinction_losses_spearman_n",
                 "extinction_losses_kendalltau_r",
                 "extinction_losses_kendalltau_p",
                 "extinction_losses_kendalltau_n",
                 "origination_fusion_spearman_r",
                 "origination_fusion_spearman_p",
                 "origination_fusion_spearman_n",
                 "origination_fusion_kendalltau_r",
                 "origination_fusion_kendalltau_p",
                 "origination_fusion_kendalltau_n",
                 "origination_losses_spearman_r",
                 "origination_losses_spearman_p",
                 "origination_losses_spearman_n",
                 "origination_losses_kendalltau_r",
                 "origination_losses_kendalltau_p",
                 "origination_losses_kendalltau_n"]
    def __init__(self, taxid, name = None) -> None:
        self.taxid = taxid
        self.parent = None
        self.children = set()
        self.name = name
        # node age estimates, in millions of years ago
        self.nodeages = Counter()
        # This is the singular node age estimate, once we are happy with one of them
        self.nodeage = None
        self.nodeageinterpolated = False
        self.lock_age = False  # True when age should not be modified (priority calibrations)
        # Get the lineage
        self.lineage = None
        self.lineage_string = ""
        self.sort_order = None
        self.x = None
        self.y = None
        # path information
        # This is the distance of all the edges from this node to the tips
        self.dist_crown = None
        # This is the distance of all the edges from this node to the tips,
        #  plus the distance of the edge leading up to this node.
        self.dist_crown_plus_root = None
        # chromsize info
        self.chromsize_median = -1
        self.chromsize_mean   = -1
        self.num_genomes      = -1
        self.chromsize_list   = []
        # fusion/loss info
        self.fusions_in_this_clade = None
        self.losses_in_this_clade = None
        # extinction correlation statistics
        self.extinction_fusion_spearman_r = None
        self.extinction_fusion_spearman_p = None
        self.extinction_fusion_spearman_n = None
        self.extinction_fusion_kendalltau_r = None
        self.extinction_fusion_kendalltau_p = None
        self.extinction_fusion_kendalltau_n = None
        self.extinction_losses_spearman_r = None
        self.extinction_losses_spearman_p = None
        self.extinction_losses_spearman_n = None
        self.extinction_losses_kendalltau_r = None
        self.extinction_losses_kendalltau_p = None
        self.extinction_losses_kendalltau_n = None
        # origination correlation statistics
        self.origination_fusion_spearman_r = None
        self.origination_fusion_spearman_p = None
        self.origination_fusion_spearman_n = None
        self.origination_fusion_kendalltau_r = None
        self.origination_fusion_kendalltau_p = None
        self.origination_fusion_kendalltau_n = None
        self.origination_losses_spearman_r = None
        self.origination_losses_spearman_p = None
        self.origination_losses_spearman_n = None
        self.origination_losses_kendalltau_r = None
        self.origination_losses_kendalltau_p = None
        self.origination_losses_kendalltau_n = None

    def __str__(self) -> str:
        outstring  = "TaxNode:\n"
        # now print out all of the fields that are currently in the object
        for slot in self.__slots__:
            if hasattr(self, slot):
                outstring += "  - {}: {}\n".format(slot, getattr(self, slot))
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
                 'parent_lineage', 'child_lineage',
                 'num_fusions_this_branch', 'num_losses_this_branch',
                 'num_fusions_per_my_this_branch', 'num_losses_per_my_this_branch',
                 'fusions', 'losses']
    def __init__(self, parent_taxid, child_taxid) -> None:
        self.parent_taxid = parent_taxid
        self.child_taxid = child_taxid
        self.parent_age = None
        self.child_age = None
        self.branch_length = None
        self.dist_crown_plus_this_edge = None
        self.parent_lineage = None
        self.child_lineage  = None
        self.num_fusions_this_branch = None
        self.num_losses_this_branch = None
        self.num_fusions_per_my_this_branch = None
        self.num_losses_per_my_this_branch = None
        self.fusions = Counter()
        self.losses  = Counter()

    def __str__(self) -> str:
        outstring  = "TaxEdge:\n"
        for slot in self.__slots__:
            #if hasattr(self, slot) and slot not in ["fusions", "losses"]:
            if hasattr(self, slot):
                outstring += "  - {}: {}\n".format(slot, getattr(self, slot))
        return outstring

class TaxIDtree:
    """
    This is a datastructure to quickly search for the most closely related species
      given search species 1 and a tree of species 2...N.
    """
    __slots__ = ["nodes", "edges", "root", "leaf_order", "NCBI", "_2d_color_params"]
    def __init__(self) -> None:
        self.nodes = {}
        self.edges = {}
        self.root = None
        self.leaf_order = []
        self.NCBI = NCBITaxa()
        self._2d_color_params = None

    def ingest_node_edge(self, node_file, edge_file):
        # Read in existing node and edge files and construct the tree.
        # Accept either filenames (str) or DataFrames
        if isinstance(node_file, str):
            nodedf = pd.read_csv(node_file, sep="\t")
        else:
            nodedf = node_file  # Already a DataFrame
        # ingest the nodes
        for idx, row in nodedf.iterrows():
            taxid = row["taxid"]
            self.add_node(taxid)
            # now go through all of the columns and add them to self.colname items
            for col in nodedf.columns:
                if col == "taxid":
                    continue
                elif col == "nodeages":
                    self.nodes[taxid].__setattr__(col, eval(row[col]))
                elif col == "lineage":
                    self.nodes[taxid].__setattr__(col, eval(row[col]))
                elif col == "chromsize_list":
                    self.nodes[taxid].__setattr__(col, eval(row[col]))
                elif col == "children":
                    self.nodes[taxid].children = set(eval(row[col]))
                else:
                    self.nodes[taxid].__setattr__(col, row[col])
            if self.nodes[taxid].num_genomes == -1:
                self.nodes[taxid].num_genomes = len(self.nodes[taxid].chromsize_list)
        # Now add the edges
        # Accept either filenames (str) or DataFrames
        if isinstance(edge_file, str):
            edgedf = pd.read_csv(edge_file, sep="\t")
        else:
            edgedf = edge_file  # Already a DataFrame
        counter = 0
        for idx, row in edgedf.iterrows():
            p = row["parent_taxid"]
            c = row["child_taxid"]
            # check that the parent and child are in the nodes
            if not p in self.nodes:
                raise ValueError("The parent taxid is not in the nodes.")
            if not c in self.nodes:
                raise ValueError("The child taxid is not in the nodes.")
            self.add_edge(p, c)
            for col in edgedf.columns:
                if col == "parent" or col == "child":
                    continue
                
                # Map new dispersals column names to old losses attribute names for backward compatibility
                attr_name = col
                if col == "num_dispersals_this_branch":
                    attr_name = "num_losses_this_branch"
                elif col == "num_dispersals_per_my_this_branch":
                    attr_name = "num_losses_per_my_this_branch"
                
                # if the column is not in self, check if it starts with num_ and ends with _this_branch
                if col.startswith("num_") and col.endswith("_this_branch") and (col not in ['num_fusions_this_branch', 'num_losses_this_branch', 'num_dispersals_this_branch', 'num_fusions_per_my_this_branch', 'num_losses_per_my_this_branch', 'num_dispersals_per_my_this_branch']):
                    # determine if it is a fusion or a loss
                    checkfield = col.split("_")[1]
                    if "+" in checkfield:
                        # it is a fusion
                        self.edges[(p, c)].fusions[col] = row[col]
                    else:
                        # it is a loss
                        self.edges[(p, c)].losses[col]  = row[col]
                else:
                    self.edges[(p, c)].__setattr__(attr_name, row[col])

    def sort_nodes(self, sort = "lineage"):
        """
        All this does is figure out the sort order for the leaves.
        It does this by sorting based on the self.lineage_string.
        It updates the self.leaf_order list.
        Then it updates the sort_order for each node.
        """
        if sort not in ["lineage", "ascending", "descending"]:
            raise ValueError("The sort parameter is not recognized.")
        if sort == "lineage":
            tempdict = {x: self.nodes[x].lineage for x in self.nodes if len(self.nodes[x].children) == 0}
            self.leaf_order = [k for k, v in sorted(tempdict.items(), key=lambda item: item[1])]
            # now update the sort order of each node
        elif sort in ["ascending", "descending"]:
            # Sorts the nodes in descending order using the least dist_crown field
            # needs to use a DFS to put things into the right order
            self.leaf_order = []
            if self.root is None:
                self.root = self.find_root()
            stack = [self.root]
            while stack:
                node = stack.pop()
                if node is None:
                    raise ValueError("The node is None. The stack is: {}".format(stack))
                if len(self.nodes[node].children) == 0:
                    self.leaf_order.append(node)
                else:
                    # We have to make some decision about which to look at first
                    #  We will look at the one with the least dist_crown first
                    children = list(self.nodes[node].children)
                    if sort == "ascending":
                        children.sort(key=lambda x: self.nodes[x].dist_crown)
                    elif sort == "descending":
                        children.sort(key=lambda x: self.nodes[x].dist_crown, reverse=True)
                    stack += children
        for i in range(len(self.leaf_order)):
            node = self.leaf_order[i]
            self.nodes[node].sort_order = i

    def plot_tree(self, ax, sort = "ascending", variable = None, lw_standard = 0.25, text_older_than = 999999999999, draw_horizontal_bars = True, split_side_mode = False, randomize_order = False, split_draw_zero_branches = False, split_buffer_pct = 0.2, split_width_scale = 2.0, magenta_color = None, bar_width_base = 1.0):
        """
        Returns an axis object with the tree plotted.
        For now, just sorts based on the ncbi lineage.
        Eventually, will plot the line thickness and color based on the variable.
        
        Parameters:
            variable: Either a single variable name (str) or a tuple of two variable names for 2D color mapping.
                     If tuple, first variable maps to red, second to blue.
            draw_horizontal_bars: If False, skips drawing horizontal lines connecting child clades
            split_side_mode: If True with 2-variable tuple, plots each variable as a half-width line on opposite sides
            randomize_order: If True, randomizes the order of edge plotting to avoid overlap artifacts
            split_draw_zero_branches: If True in split_side_mode, draws thin gray lines for branches with both values == 0
            split_buffer_pct: Percentage (0-1) of max tip spacing to use as buffer/whitespace (default 0.2 = 20%)
            split_width_scale: Multiplier for how much width increases from tips to root (default 2.0)
        """
        self.sort_nodes(sort)
        for i in range(len(self.leaf_order)):
            node = self.leaf_order[i]
            self.nodes[node].x = i
        
        # Check if we're doing 2D color mapping
        is_2d_color = isinstance(variable, tuple) and len(variable) == 2
        
        # Calculate spacing and max age for variable-width in 2D modes
        if is_2d_color:
            # Leaf spacing is 1.0 in x-axis units (since leaves are at integer positions)
            leaf_spacing = 1.0
            # Max width at tips = leaf_spacing * (1 - buffer_pct) for all 2D modes
            max_tip_width = leaf_spacing * (1.0 - split_buffer_pct)
            # Get max age (root age) for scaling
            max_age = max([self.nodes[node].nodeage for node in self.nodes if self.nodes[node].nodeage is not None])
        else:
            max_tip_width = 0.3  # Default value for non-split mode
            max_age = 1.0

        # Check if we're doing 2D color mapping
        is_2d_color = isinstance(variable, tuple) and len(variable) == 2
        
        # Store 2D color mapping parameters for legend creation
        self._2d_color_params = None
        
        if variable is not None:
            if is_2d_color:
                # Handle 2D color mapping
                var1, var2 = variable
                minv1 = min([y for y in [self.edges[x].__getattribute__(var1) for x in self.edges] if y != 0]) * 0.1
                maxv1 = max([y for y in [self.edges[x].__getattribute__(var1) for x in self.edges] if y != float("inf")]) + minv1
                minv2 = min([y for y in [self.edges[x].__getattribute__(var2) for x in self.edges] if y != 0]) * 0.1
                maxv2 = max([y for y in [self.edges[x].__getattribute__(var2) for x in self.edges] if y != float("inf")]) + minv2
                print(f"2D color mapping: {var1} [{minv1:.6f}, {maxv1:.2f}] vs {var2} [{minv2:.6f}, {maxv2:.2f}]")
                # Store for legend creation
                self._2d_color_params = (minv1, maxv1, minv2, maxv2)
            else:
                # get the min and max of the variable
                minv = min([y for y in [self.edges[x].__getattribute__(variable) for x in self.edges] if y != 0]) * 0.1
                # get the max that isn't infinity
                maxv = max([y for y in [self.edges[x].__getattribute__(variable) for x in self.edges] if y != float("inf")]) + minv
                print("We're trying variable: {} with min: {} and max: {}".format(variable, minv, maxv))

        def _log_normalize(value, min_val_linear, max_val_linear):
            # Step 1: Convert the min and max values to log space
            min_val_log = np.log(min_val_linear)
            max_val_log = np.log(max_val_linear)

            # Step 2: Apply logarithmic transformation to the input values
            log_transformed = np.log(value)

            # Step 3: Min-Max scaling to [0, 1] using the log-transformed min and max values
            scaled_values = (log_transformed - min_val_log) / (max_val_log - min_val_log)
            # check to make sure that there are no NaNs
            if np.isnan(scaled_values):
                raise ValueError("The scaled values are NaN. The input value was: {}, input min was {}, input max was {}".format(value, min_val_linear, max_val_linear))

            return scaled_values


        # now, do a reverse BFS from the tips and plot everything.
        # Things that are plotted are the edges. This can also just be inferred from the nodes.
        plotted = set()
        edge_drawn = set()
        if randomize_order:
            import random
            leaf_list = list(self.leaf_order)
            random.shuffle(leaf_list)
            queue = deque(leaf_list)
        else:
            queue = deque(self.leaf_order)
        print()
        while queue:
            node = queue.popleft()
            print(f"  - There are {len(queue)} nodes left in the queue.", end="\r")
            #print("  - self.nodes[x].x for the children is: {}".format([self.nodes[x].x for x in self.nodes[node].children]))
            # If there are children, first check if all of them are plotted, if not,
            #  we must skip this node for now and add it to the end of the queue.
            if len(self.nodes[node].children) > 0:
                if not all([child in plotted for child in self.nodes[node].children]):
                    queue.append(node)
                    continue

            # THIS ONLY HAPPENS IF WE HAVE NOT CONTINUED
            # DRAWING THE clade hat (These are drawn with a constant y, changing x)
            # we need to draw the line from the node to the parent.
            # if there are no children, we do not need to draw the vertical line
            if len(self.nodes[node].children) == 0:
                # we already set the x-value to be the sort order outside of this while loop
                pass
            elif len(self.nodes[node].children) == 1:
                # In this case, there are no children, or just one. In this case we draw no line
                #  connecting children. We just set x to whatever the child's x is
                self.nodes[node].x = self.nodes[list(self.nodes[node].children)[0]].x
            else:
                # We have to figure out how to draw the vertical line first0.15
                #  Again, y is constant, changing x
                y = self.nodes[node].nodeage
                x_min = min([self.nodes[x].x for x in self.nodes[node].children])
                if type(x_min) not in [int, float]:
                    raise ValueError("The x_min is not a number.")
                x_max = max([self.nodes[x].x for x in self.nodes[node].children])
                if type(x_max) not in [int, float]:
                    raise ValueError("The x_max is not a number.")
                if draw_horizontal_bars:
                    if split_side_mode:
                        # Draw horizontal bars first with low zorder so colored bars overlay on top
                        ax.plot([x_min, x_max], [y, y], color="grey", linewidth=lw_standard, zorder=1)
                    else:
                        ax.plot([x_min, x_max], [y, y], color="grey", linewidth=lw_standard, zorder=1)
                # We now must set the x value of the node to be the average of the
                self.nodes[node].x = (x_min + x_max)/2

            # DRAWING THE BRANCHES (Drawn with a constant x, changing y)
            if (self.nodes[node].parent is not None) and (self.nodes[node].parent != -1):
                x       = self.nodes[node].x
                y_start = self.nodes[node].nodeage
                p       = self.nodes[node].parent
                y_stop  = self.nodes[p].nodeage
                e       = (p, node)
                if e not in edge_drawn:
                    if variable is not None:
                        if is_2d_color:
                            # Get both variable values
                            value1 = self.edges[(p, node)].__getattribute__(var1)
                            value2 = self.edges[(p, node)].__getattribute__(var2)
                            # Handle NaN and infinity
                            if np.isnan(value1) or value1 == float("inf"):
                                value1 = 0
                            if np.isnan(value2) or value2 == float("inf"):
                                value2 = 0
                            
                            if split_side_mode:
                                # Import Rectangle for split-side mode
                                from matplotlib.patches import Rectangle
                                
                                # Calculate variable bar width based on child node age
                                # Width increases linearly from tips (age=0) to root (age=max_age)
                                child_age = self.nodes[node].nodeage
                                age_fraction = child_age / max_age if max_age > 0 else 0
                                # Width = max_tip_width * (bar_width_base + age_fraction * split_width_scale)
                                bar_width = max_tip_width * (bar_width_base + age_fraction * split_width_scale) / 2.0  # Divide by 2 for half-width per side
                                
                                # Check if both values are zero
                                if value1 == 0 and value2 == 0:
                                    # Draw thin light grey line for zero branches if requested (zorder=1, colored bars have higher zorder)
                                    if split_draw_zero_branches:
                                        ax.plot([x, x], [y_start, y_stop], color="lightgrey", linewidth=0.2, solid_capstyle='butt', zorder=1)
                                else:
                                    # Plot dispersals (var1, red) on the left side as a rectangle
                                    if value1 != 0:
                                        value1 += minv1
                                        norm1 = _log_normalize(value1, minv1, maxv1)
                                        # Create gradient from white to more saturated red (hex color)
                                        red_max = mcolors.to_rgb("#D22C16")
                                        cmap_red = mcolors.LinearSegmentedColormap.from_list("", ["white", red_max])
                                        color1 = cmap_red(norm1)
                                        # Rectangle: (x, y, width, height)
                                        rect1 = Rectangle((x - bar_width, y_start), bar_width, y_stop - y_start, 
                                                         facecolor=color1, edgecolor='none', zorder=5)
                                        ax.add_patch(rect1)
                                    
                                    # Plot fusions (var2, blue) on the right side as a rectangle
                                    if value2 != 0:
                                        value2 += minv2
                                        norm2 = _log_normalize(value2, minv2, maxv2)
                                        # Create gradient from white to more saturated blue (hex color)
                                        blue_max = mcolors.to_rgb("#3054A3")
                                        cmap_blue = mcolors.LinearSegmentedColormap.from_list("", ["white", blue_max])
                                        color2 = cmap_blue(norm2)
                                        # Rectangle: (x, y, width, height)
                                        rect2 = Rectangle((x, y_start), bar_width, y_stop - y_start, 
                                                         facecolor=color2, edgecolor='none', zorder=5)
                                        ax.add_patch(rect2)
                            else:
                                # Standard bivariate blended color mode with variable line width
                                # If both values are 0 and draw_horizontal_bars is True, draw as grey structure line
                                if value1 == 0 and value2 == 0:
                                    if draw_horizontal_bars:
                                        # Draw zero branch as light grey structure line
                                        ax.plot([x, x], [y_start, y_stop], color="lightgrey", linewidth=0.2, zorder=1)
                                else:
                                    # Import Rectangle for bivariate mode
                                    from matplotlib.patches import Rectangle
                                    
                                    value1 += minv1
                                    value2 += minv2
                                    # Normalize both values
                                    norm1 = _log_normalize(value1, minv1, maxv1)
                                    norm2 = _log_normalize(value2, minv2, maxv2)
                                    # Create 2D color using bilinear interpolation between corner colors
                                    # Define corner colors using more saturated hex colors
                                    white = np.array([1.0, 1.0, 1.0])
                                    blue_max = np.array(mcolors.to_rgb("#3054A3"))   # More saturated blue at (1,0)
                                    red_max = np.array(mcolors.to_rgb("#D22C16"))     # More saturated red at (0,1)
                                    # Use provided magenta_color or default to softer magenta
                                    magenta = np.array(magenta_color if magenta_color is not None else [0.6, 0.3, 0.6])  # Magenta at (1,1)
                                    # Bilinear interpolation
                                    color_rgb = ((1-norm1)*(1-norm2)*white + 
                                                norm1*(1-norm2)*blue_max + 
                                                (1-norm1)*norm2*red_max + 
                                                norm1*norm2*magenta)
                                    color = tuple(color_rgb)
                                    # Calculate variable bar width based on child node age
                                    child_age = self.nodes[node].nodeage
                                    age_fraction = child_age / max_age if max_age > 0 else 0
                                    # Width increases linearly from tips (max_tip_width) to root (max_tip_width * (bar_width_base + split_width_scale))
                                    bar_width = max_tip_width * (bar_width_base + age_fraction * split_width_scale)
                                    # Draw as centered rectangle (not split left/right like split_side_mode)
                                    rect = Rectangle((x - bar_width/2, y_start), bar_width, y_stop - y_start, 
                                                     facecolor=color, edgecolor='none', zorder=5)
                                    ax.add_patch(rect)
                        else:
                            # get the value of the variable for this edge
                            value = self.edges[(p, node)].__getattribute__(variable)
                            # if the value is nan, convert to zero
                            if np.isnan(value):
                                value = 0
                            # if the value is infinity, we plot it as the max value
                            if value == float("inf"):
                                value = 0
                            # Skip plotting if value is 0 (would be white)
                            if value != 0:
                                value += minv
                                # now we need to normalize the value
                                value = _log_normalize(value, minv, maxv)
                                # Use constant linewidth to avoid overlap issues - rely on color for encoding
                                lw = 2.5
                                # Use white-to-more-saturated-color gradients (hex colors)
                                if "fusion" in variable:
                                    blue_max = mcolors.to_rgb("#3054A3")
                                    cmap = mcolors.LinearSegmentedColormap.from_list("", ["white", blue_max])
                                elif "loss" in variable or "dispersal" in variable:
                                    red_max = mcolors.to_rgb("#D22C16")
                                    cmap = mcolors.LinearSegmentedColormap.from_list("", ["white", red_max])
                                color = cmap(value)
                                ax.plot([x, x], [y_start, y_stop], color=color, linewidth=lw, alpha=1.0, solid_capstyle='butt')
                    else:
                        # Only draw grey tree structure if not in split_side_mode
                        if not split_side_mode:
                            ax.plot([x, x], [y_start, y_stop], color="grey", linewidth=lw_standard, solid_capstyle='butt', zorder=1)
                    # now add the parent node
                    if p not in plotted:
                        queue.append(p)
                    edge_drawn.add(e)
            # we now plotted the node
            plotted.add(node)
        print()
        
        # Set explicit axis limits to show full tree extent
        # X-axis: based on number of leaves (sorted order)
        # Y-axis: based on node ages
        # Required for modes using Rectangle patches (split_side_mode and bivariate with is_2d_color)
        if split_side_mode or is_2d_color:
            # Calculate tree extent
            x_positions = [self.nodes[node].x for node in self.nodes if hasattr(self.nodes[node], 'x')]
            y_positions = [self.nodes[node].nodeage for node in self.nodes if self.nodes[node].nodeage is not None]
            
            if x_positions and y_positions:
                x_min, x_max = min(x_positions), max(x_positions)
                y_min, y_max = min(y_positions), max(y_positions)
                # Add small padding
                x_padding = (x_max - x_min) * 0.02
                y_padding = (y_max - y_min) * 0.02
                ax.set_xlim(x_min - x_padding, x_max + x_padding)
                ax.set_ylim(y_min - y_padding, y_max + y_padding)

        # Let's add text to the nodes of the plot.
        #  If the node is older than text_older_than, we will add the text.
        #  The alignment will different depending on the sort
        for node in self.nodes:
            if self.nodes[node].nodeage > text_older_than:
                nodename = self.nodes[node].name
                # check if none or nan (use pd.isna which works for all types)
                if pd.isna(nodename) or nodename == "":
                    nodename = f"Node{node}"  # Use taxid as fallback name
                if sort == "ascending":
                    ax.text(self.nodes[node].x, self.nodes[node].nodeage, nodename, ha="right", va="bottom", rotation=-45)
                elif sort == "descending":
                    ax.text(self.nodes[node].x, self.nodes[node].nodeage, nodename, ha="left", va="top", rotation=45)
                else:
                    ax.text(self.nodes[node].x, self.nodes[node].nodeage, nodename, ha="center", va="bottom", rotation=0)
        return ax
    
    @staticmethod
    def create_2d_colorbar_legend(ax, var1_name="Dispersals", var2_name="Fusions", 
                                  var1_range=None, var2_range=None, resolution=100, magenta_color=None):
        """
        Create a 2D color legend showing the bivariate color mapping.
        
        Parameters:
            ax: matplotlib axis to draw on
            var1_name: label for first variable (vertical axis, red channel)
            var2_name: label for second variable (horizontal axis, blue channel)
            var1_range: tuple of (min, max) for var1, if None uses (0, 1)
            var2_range: tuple of (min, max) for var2, if None uses (0, 1)
            resolution: number of grid points in each dimension
        """
        # Get the actual data ranges
        if var1_range is None:
            minv1, maxv1 = 0.001, 1.0  # Default log-safe range
        else:
            minv1, maxv1 = var1_range
            
        if var2_range is None:
            minv2, maxv2 = 0.001, 1.0  # Default log-safe range
        else:
            minv2, maxv2 = var2_range
        
        # Create log-spaced grids in the actual data space
        x_data = np.logspace(np.log10(minv2), np.log10(maxv2), resolution)
        y_data = np.logspace(np.log10(minv1), np.log10(maxv1), resolution)
        X_data, Y_data = np.meshgrid(x_data, y_data)
        
        # Normalize using log normalization (matching plot_tree)
        def log_normalize(value, min_val, max_val):
            min_log = np.log(min_val)
            max_log = np.log(max_val)
            log_val = np.log(value)
            return (log_val - min_log) / (max_log - min_log)
        
        # Apply log normalization to get color coordinates
        X_norm = log_normalize(X_data, minv2, maxv2)
        Y_norm = log_normalize(Y_data, minv1, maxv1)
        
        # Use proper white-to-color colormaps
        cmap_red = plt.cm.Reds
        cmap_blue = plt.cm.Blues
        
        # Define corner colors for bilinear interpolation using saturated hex colors
        white = np.array([1.0, 1.0, 1.0])
        blue_max = np.array(mcolors.to_rgb("#3054A3"))  # More saturated blue
        red_max = np.array(mcolors.to_rgb("#D22C16"))  # More saturated red
        # Use provided magenta_color or default to softer magenta
        magenta = np.array(magenta_color if magenta_color is not None else [0.6, 0.3, 0.6])
        
        # Create RGB color array using bilinear interpolation
        colors = np.zeros((resolution, resolution, 3))
        for i in range(resolution):
            for j in range(resolution):
                # Bilinear interpolation using normalized coordinates
                color = ((1-X_norm[i,j])*(1-Y_norm[i,j])*white + 
                        X_norm[i,j]*(1-Y_norm[i,j])*blue_max + 
                        (1-X_norm[i,j])*Y_norm[i,j]*red_max + 
                        X_norm[i,j]*Y_norm[i,j]*magenta)
                colors[i, j] = color
        
        # Display the color grid with square aspect ratio in log space
        ax.imshow(colors, origin='lower', extent=[np.log10(minv2), np.log10(maxv2), 
                                                   np.log10(minv1), np.log10(maxv1)], 
                  aspect='equal')
        ax.set_xlabel(f'{var2_name} (per MY) →', fontsize=12)
        ax.set_ylabel(f'{var1_name} (per MY) →', fontsize=12)
        ax.set_title('2D Color Scale (White = Low, Purple = Both High)', fontsize=14)
        
        # Set log-spaced ticks with actual values
        # For x-axis (var2/fusions)
        x_ticks_log = np.linspace(np.log10(minv2), np.log10(maxv2), 5)
        x_ticks_actual = 10**x_ticks_log
        ax.set_xticks(x_ticks_log)
        ax.set_xticklabels([f'{v:.3f}' if v < 1 else f'{v:.1f}' for v in x_ticks_actual])
        
        # For y-axis (var1/dispersals)
        y_ticks_log = np.linspace(np.log10(minv1), np.log10(maxv1), 5)
        y_ticks_actual = 10**y_ticks_log
        ax.set_yticks(y_ticks_log)
        ax.set_yticklabels([f'{v:.3f}' if v < 1 else f'{v:.1f}' for v in y_ticks_actual])
        
        return ax

    def find_root(self) -> int:
        """
        Finds the root and returns its int. The root is defined as something without a parent.
        """
        potential_roots = []
        for node in self.nodes:
            if node == 1:
                print(self.nodes[node])
            if (self.nodes[node].parent is None) or (self.nodes[node].parent == -1) or (self.nodes[node].parent == np.nan):
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
        Safely adds, so we can call this multiple times.
        Doesn't fail if the node already exists.
        """
        if taxid not in self.nodes:
            self.nodes[taxid] = TaxNode(taxid)
            # For negative taxids (custom nodes like Myriazoa=-67), use a generic name
            if taxid < 0:
                clade_name = f"CustomNode{abs(taxid)}"
            else:
                # For positive taxids, use NCBI to get the name
                clade_name = self.NCBI.get_taxid_translator([taxid])[taxid]
                # Replace spaces with underscores to preserve all fields (e.g., "Homo sapiens neanderthalensis" -> "Homo_sapiens_neanderthalensis")
                if " " in clade_name:
                    clade_name = "_".join(clade_name.split(" "))
            self.nodes[taxid].name = clade_name
        return self.nodes[taxid]

    def add_edge(self, parent_taxid, child_taxid) -> int:
        """
        Adds an edge between two nodes.
        """
        # add the nodes. The addition is safe, so we can call this multiple times.
        self.add_node(parent_taxid)
        self.add_node(child_taxid)
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
    
    def build_from_newick_tree(self, newick_tree, ncbi):
        """
        Builds the TaxIDtree structure from a newick tree topology.
        
        This preserves the custom topology from the newick file instead of 
        rebuilding from NCBI lineages.
        
        Parameters:
        -----------
        newick_tree : newick.Node
            The newick tree object to build from
        ncbi : NCBITaxa
            NCBI taxonomy database for species name to taxid conversion
        """
        print("\nBuilding TaxIDtree from newick topology...")
        
        # Map species names to taxids
        leaf_name_to_taxid = {}
        for leaf in list(newick_tree.leaves()):
            # Skip leaves with None names
            if leaf.name is None:
                print(f"  Warning: Skipping newick tree leaf with None name")
                continue
            # Parse taxid from node name if it's in brackets, e.g., "Homo_sapiens[9606]"
            if '[' in leaf.name and leaf.name.endswith(']'):
                taxid = int(leaf.name[leaf.name.rfind('[')+1:-1])
                leaf_name_to_taxid[leaf.name] = taxid
            else:
                # Fallback to old method if no brackets found
                taxname = leaf.name.replace("_", " ")
                try:
                    taxid = ncbi.get_name_translator([taxname])[taxname][0]
                    leaf_name_to_taxid[leaf.name] = taxid
                except:
                    print(f"  Warning: Could not find taxid for: {leaf.name}")
                    continue
        
        print(f"  Mapped {len(leaf_name_to_taxid)} species to taxids")
        
        # Build tree structure by traversing the newick tree
        # We'll assign internal node IDs based on a simple counter for custom nodes
        # or use actual taxids when we can infer them from NCBI lineages
        node_to_id = {}
        next_internal_id = -1000  # Use negative IDs for internal nodes without bracket notation
        custom_taxids_found = []  # Track custom taxids we extract
        
        def traverse_and_build(node, parent_id=None):
            nonlocal next_internal_id
            
            # Initialize custom_name for all node types
            custom_name = None
            
            # Assign an ID to this node
            if node.is_leaf:
                # Leaf nodes get their taxid
                if node.name in leaf_name_to_taxid:
                    node_id = leaf_name_to_taxid[node.name]
                else:
                    return None  # Skip nodes we couldn't map
            else:
                # Internal nodes - extract taxid from bracket notation if present
                if node.name and '[' in node.name and node.name.endswith(']'):
                    try:
                        # Extract taxid from brackets, e.g., "Myriazoa[-67]" -> -67
                        node_id = int(node.name[node.name.rfind('[')+1:-1])
                        # Extract the name part before the bracket
                        custom_name = node.name[:node.name.rfind('[')]
                        # Track custom taxids (negative ones)
                        if node_id < 0:
                            custom_taxids_found.append((node.name, node_id))
                    except (ValueError, IndexError):
                        # If parsing fails, use sequential negative IDs
                        node_id = next_internal_id
                        next_internal_id -= 1
                else:
                    # No bracket notation, use sequential negative IDs
                    node_id = next_internal_id
                    next_internal_id -= 1
            
            node_to_id[node] = node_id
            
            # Add edge from parent if this isn't the root
            if parent_id is not None:
                self.add_edge(parent_id, node_id)
            else:
                # Root node - add it explicitly since it has no parent
                self.add_node(node_id)
            
            # Override the name with the extracted name from newick file (preserves custom names like "Myriazoa")
            if custom_name:
                self.nodes[node_id].name = custom_name
            
            # Recursively process children
            for child in node.children:
                traverse_and_build(child, node_id)
            
            return node_id
        
        # Start traversal from root
        root_id = traverse_and_build(newick_tree)
        print(f"  Built tree with {len(self.nodes)} nodes and {len(self.edges)} edges")
        print(f"  Root node ID: {root_id}")
        
        # Report custom taxids found
        if custom_taxids_found:
            print(f"  Found {len(custom_taxids_found)} custom taxids from bracket notation:")
            for name, taxid in custom_taxids_found:
                print(f"    {name} -> taxid {taxid}")
        
        return root_id

    def add_chromosome_info_file(self, chrominfo_file):
        """
        Adds chromosome info to the nodes from a file.
        The format of the file is:
        f"%s\t%d" % (sample_string, num_chromosomes)

        The algorithm is: Load the dataframe
        """
        chromdf = pd.read_csv(chrominfo_file, sep="\t", header=0)
        # Rename columns to standardized names for internal use
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
            self.nodes[node].num_genomes = len(self.nodes[node].chromsize_list)

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
        - NCBI_object - is the output of ete4's NCBITaxa(). Use this to get the lineage information without having to download it.
        - query_taxid - the taxid of the species we are looking for.
        
        Returns:
        - The taxid of the closest relative in this tree, or None if no match found.
        """
        query_lineage = NCBI_object.get_lineage(query_taxid)[::-1]
        #print("The query lineage is: {}".format(query_lineage))
        # go through the query_lineage and return the first taxid that is in the tree
        target_node = None
        for i in range(len(query_lineage)):
            if query_lineage[i] in self.nodes:
                target_node = query_lineage[i]
                break
        
        # If no matching taxid found in tree, return None
        if target_node is None:
            return None
        
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
        # BUT: Only apply these if the node doesn't already have an age (e.g., from TimeTree)
        changes = {1: 3700, 131567: 3480, 2759: 1898, 33154: 1010}
        for node in changes:
            if node in self.nodes and len(self.nodes[node].nodeages) == 0:
                self.nodes[node].nodeages.update([changes[node]])

        # First, just set every node to the most common value.
        # Skip locked nodes (priority calibrations) - don't overwrite their ages
        for node in self.nodes:
            if len(self.nodes[node].nodeages) > 0:
                # Only update nodeage if it's not locked
                if not self.nodes[node].lock_age:
                    self.nodes[node].nodeage = self.nodes[node].nodeages.most_common(1)[0][0]

        # There are now a lot of nodes that have ages, but we need to make sure that the children are younger than the parents.
        # Alternate between descending (root→tip) and ascending (tip→root) passes
        iteration = 0
        done = False
        max_iterations = 10  # Limit iterations to prevent infinite loops
        while not done and iteration < max_iterations:
            if iteration % 2 == 0:
                print("We are percolating from root to tip")
                self.percolate_descending_root_to_tip()
            else:
                print("We are percolating from tip to root")
                self.percolate_ascending_tip_to_root()
            descending = self.check_nodeages_descending()
            if descending:
                done = True
            print("We are on iteration {}. The descending state is {}".format(iteration, descending))
            iteration += 1
        
        if not descending:
            print(f"Warning: Node ages not fully ordered after {max_iterations} iterations, continuing anyway")

    def percolate_ascending_tip_to_root(self) -> None:
        """
        This does a pass through the tree from leaves to root and ensures parents are older than children.
        Uses reverse BFS (starting from leaves).
        """
        # Start with all leaf nodes
        queue = deque([node for node in self.nodes if len(self.nodes[node].children) == 0])
        visited = set()
        
        while len(queue) > 0:
            child = queue.popleft()
            if child in visited:
                continue
            visited.add(child)
            
            parent = self.nodes[child].parent
            if parent is None:
                continue
                
            # Skip locked nodes - don't modify their ages
            if self.nodes[parent].lock_age:
                # Still add parent to queue to continue traversal
                if parent not in visited:
                    queue.append(parent)
                continue
            
            # Get the maximum child age to ensure parent is older
            child_ages = []
            for sibling in self.nodes[parent].children:
                if len(self.nodes[sibling].nodeages) > 0 and self.nodes[sibling].nodeage is not None:
                    child_ages.append(self.nodes[sibling].nodeage)
            
            if len(child_ages) > 0 and len(self.nodes[parent].nodeages) > 0:
                max_child_age = max(child_ages)
                # If current parent age is younger than any child, pick a new age
                if self.nodes[parent].nodeage is not None and self.nodes[parent].nodeage < max_child_age:
                    # Try to find a parent age that is older than all children
                    for age in sorted(self.nodes[parent].nodeages.keys(), reverse=True):
                        if age >= max_child_age:
                            self.nodes[parent].nodeage = age
                            break
            
            # Add parent to queue for further processing
            if parent not in visited:
                queue.append(parent)
    
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
            # Skip locked nodes (priority calibrations) - don't modify their ages
            if self.nodes[parent].lock_age:
                continue
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

    def correct_missing_nodes(self, priority_node_ages=None) -> None:
        """
        Goes through the tree and interpolates the node ages based on the ages of nodes above and below it.
        This only assigns ages to the nodes that are missing info.
        
        Parameters:
        -----------
        priority_node_ages : dict, optional
            Dictionary mapping node_id -> (age_mya, description) for priority calibration nodes.
            These ages are applied before general interpolation to anchor the tree.
            Example: {33213: (707.0, "Bilateria (Homo sapiens - Drosophila melanogaster)")}
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
        # BUT: Only apply these if the node doesn't already have an age (e.g., from TimeTree)
        changes = {1: 3700, 131567: 3480, 2759: 1898, 33154: 1010}
        for node in changes:
            if node in self.nodes and len(self.nodes[node].nodeages) == 0:
                self.nodes[node].nodeages.update([changes[node]])
                print(f"  Applied fallback age {changes[node]} MYA to node {node}")
        
        # Ensure the root has an age (important for custom topology trees with non-standard root IDs)
        if len(self.nodes[self.root].nodeages) == 0:
            self.nodes[self.root].nodeages.update([3700])  # Default to root of life age
            print(f"  Set root node {self.root} age to 3700 MYA (custom topology root)")

        # Apply priority node ages (from TimeTree) before general interpolation
        # This ensures key calibration points (Bilateria, Cnidaria, etc.) are fixed first
        if priority_node_ages is not None:
            print("\n  Applying priority node ages from TimeTree:")
            for node_id, (age, description) in priority_node_ages.items():
                if node_id in self.nodes:
                    # Lock the node by setting nodeage directly and marking as locked
                    self.nodes[node_id].nodeage = age
                    self.nodes[node_id].nodeageinterpolated = False
                    self.nodes[node_id].nodeages = Counter([age])
                    self.nodes[node_id].lock_age = True  # Lock this node from future modifications
                    node_name = self.nodes[node_id].name if self.nodes[node_id].name else f"Node{node_id}"
                    print(f"    {node_name} [{node_id}]: {age:.2f} MYA ({description}) [LOCKED]")
                else:
                    print(f"    WARNING: Priority node {node_id} not found in tree, skipping ({description})")
            print()

        # There are now a lot of nodes that have ages, but we need to make sure that the children are older than the parents.
        # Do a BFS from the root.
        queue = [self.root]
        while len(queue) > 0:
            parent = queue.pop(0)
            # Skip locked nodes (priority calibrations) - don't modify their ages
            if self.nodes[parent].lock_age:
                for child in self.nodes[parent].children:
                    queue.append(child)
                continue
            acceptable_parent_ages = [x for x in self.nodes[parent].nodeages]
            all_child_ages = []
            #print("The parent is: {}. The parent ages are {}.".format(parent, self.nodes[parent].nodeages))
            for child in self.nodes[parent].children:
                # Skip age comparison for locked children - they have fixed calibration ages
                if self.nodes[child].lock_age:
                    queue.append(child)
                    continue
                if len(acceptable_parent_ages) > 0:
                    # This only works if there are some parent ages to compare to
                    # if the child has an age, we need to make sure that it is older than the parent.
                    child_ages = [x for x in self.nodes[child].nodeages]
                    all_child_ages.append(child_ages)
                    # remove the acceptable parent ages that are smaller than all of the child ages
                    # Allow for floating point tolerance (0.01 MYA = 10,000 years)
                    if len(child_ages) > 0:
                        # we need this condition because the logic doesn't work if the child ages list is empty
                        remove_parent_age = [p for p in acceptable_parent_ages if all([(p + 0.01) < c for c in child_ages])]
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
            # Handle floating point precision - if difference is tiny, treat as zero
            if abs(distance) < 1e-10:  # 1e-10 MYA = essentially zero
                distance = 0
            # check that the distance is positive (after rounding)
            elif distance < 0:
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

    def fix_zero_length_branches(self, tolerance=0.01, max_iterations=5) -> None:
        """
        Fixes zero-length branches (nodes with same age as parent within tolerance).
        
        Strategy:
        1. Pass 1: For childless nodes equal to parent, push PARENT up (toward grandparent)
        2. Pass 2: For internal nodes equal to parent, push NODE down (toward oldest child)
        
        Iterates up to max_iterations times to resolve all zero-length branches.
        
        Parameters:
        -----------
        tolerance : float
            Age difference threshold in MYA (default 0.01 = 10,000 years)
        max_iterations : int
            Maximum number of correction passes (default 5)
        """
        print(f"\nFixing zero-length branches (tolerance={tolerance} MYA, max_iterations={max_iterations})...")
        
        for iteration in range(max_iterations):
            changes_made = 0
            
            # PASS 1: Fix childless nodes - push parent UP (older)
            childless_fixed = 0
            for node_id in self.nodes:
                node = self.nodes[node_id]
                parent_id = node.parent
                
                # Skip root, locked nodes, and nodes with children
                if parent_id is None or node.lock_age or len(node.children) > 0:
                    continue
                
                parent = self.nodes[parent_id]
                if parent.lock_age:  # Don't modify locked parents
                    continue
                    
                # Check if node age equals parent age (within tolerance)
                if node.nodeage is not None and parent.nodeage is not None:
                    diff = abs(parent.nodeage - node.nodeage)
                    
                    if diff <= tolerance:
                        # Push parent UP toward grandparent
                        grandparent_id = parent.parent
                        
                        if grandparent_id is not None:
                            grandparent = self.nodes[grandparent_id]
                            if grandparent.nodeage is not None:
                                # Set parent to midpoint between grandparent and current parent
                                new_parent_age = (grandparent.nodeage + parent.nodeage) / 2.0
                                # Ensure we actually increase the age
                                if new_parent_age > parent.nodeage:
                                    parent.nodeage = new_parent_age
                                    parent.nodeages = Counter([new_parent_age])
                                    childless_fixed += 1
                                    changes_made += 1
            
            # PASS 2: Fix internal nodes with children - push node DOWN (younger)
            internal_fixed = 0
            for node_id in self.nodes:
                node = self.nodes[node_id]
                parent_id = node.parent
                
                # Skip root, locked nodes, and childless nodes
                if parent_id is None or node.lock_age or len(node.children) == 0:
                    continue
                
                parent = self.nodes[parent_id]
                
                # Check if node age equals parent age (within tolerance)
                if node.nodeage is not None and parent.nodeage is not None:
                    diff = abs(parent.nodeage - node.nodeage)
                    
                    if diff <= tolerance:
                        # Push node DOWN toward oldest child
                        # Find the oldest child (child with maximum age)
                        oldest_child_age = None
                        for child_id in node.children:
                            child = self.nodes[child_id]
                            if child.nodeage is not None:
                                if oldest_child_age is None or child.nodeage > oldest_child_age:
                                    oldest_child_age = child.nodeage
                        
                        if oldest_child_age is not None:
                            # Set node to midpoint between parent and oldest child
                            new_node_age = (parent.nodeage + oldest_child_age) / 2.0
                            # Ensure we actually decrease the age and maintain parent > node > child
                            if new_node_age < parent.nodeage and new_node_age > oldest_child_age:
                                node.nodeage = new_node_age
                                node.nodeages = Counter([new_node_age])
                                internal_fixed += 1
                                changes_made += 1
            
            print(f"  Iteration {iteration + 1}: Fixed {childless_fixed} childless nodes, {internal_fixed} internal nodes (total: {changes_made} changes)")
            
            # If no changes were made, we're done
            if changes_made == 0:
                print(f"  Converged after {iteration + 1} iteration(s)")
                break
        
        # Final report
        remaining_zero_length = 0
        for node_id in self.nodes:
            node = self.nodes[node_id]
            parent_id = node.parent
            if parent_id is not None:
                parent = self.nodes[parent_id]
                if node.nodeage is not None and parent.nodeage is not None:
                    diff = abs(parent.nodeage - node.nodeage)
                    if diff <= tolerance:
                        remaining_zero_length += 1
        
        print(f"  After fixing: {remaining_zero_length} zero-length branches remain")

    def analyze_zero_length_branches(self, tolerance=0.01, label="") -> dict:
        """
        Analyzes and reports on zero-length branches in the tree.
        
        Parameters:
        -----------
        tolerance : float
            Age difference threshold in MYA (default 0.01 = 10,000 years)
        label : str
            Optional label for the report (e.g., "BEFORE FIXING", "AFTER FIXING")
        
        Returns:
        --------
        dict : Statistics including counts and examples
        """
        if label:
            print(f"\n{'='*80}")
            print(f"{label}: Checking for nodes with essentially same age as parent...")
            print(f"{'='*80}")
        else:
            print(f"\nAnalyzing zero-length branches (tolerance={tolerance} MYA)...")
        
        same_age_count = 0
        same_age_examples = []
        child_branch_lengths = []
        zero_length_nodes_with_children = 0
        
        # First pass: count zero-length branches and collect examples
        for node_id in self.nodes:
            parent_id = self.nodes[node_id].parent
            if parent_id is not None:
                parent = self.nodes[parent_id]
                node = self.nodes[node_id]
                if parent.nodeage is not None and node.nodeage is not None:
                    diff = abs(parent.nodeage - node.nodeage)
                    if diff <= tolerance:
                        same_age_count += 1
                        if len(same_age_examples) < 5:
                            same_age_examples.append({
                                'node_name': node.name if node.name else f"Node{node_id}",
                                'node_id': node_id,
                                'node_age': node.nodeage,
                                'parent_name': parent.name if parent.name else f"Node{parent_id}",
                                'parent_id': parent_id,
                                'parent_age': parent.nodeage,
                                'diff': diff
                            })
                        
                        # If this zero-length node has children, collect child branch lengths
                        if len(node.children) > 0:
                            zero_length_nodes_with_children += 1
                            for child_id in node.children:
                                child_node = self.nodes[child_id]
                                if child_node.nodeage is not None:
                                    child_branch_length = node.nodeage - child_node.nodeage
                                    if child_branch_length > 0:
                                        child_branch_lengths.append(child_branch_length)
        
        # Report findings
        print(f"  Found {same_age_count} nodes with age within {tolerance} MYA of their parent")
        print(f"  This is {same_age_count}/{len(self.nodes)} = {100*same_age_count/len(self.nodes):.2f}% of all nodes")
        
        if same_age_examples:
            print(f"\n  First {len(same_age_examples)} examples:")
            for ex in same_age_examples:
                print(f"    Node {ex['node_name']} [{ex['node_id']}] age={ex['node_age']:.4f} MYA")
                print(f"      Parent {ex['parent_name']} [{ex['parent_id']}] age={ex['parent_age']:.4f} MYA")
                print(f"      Difference: {ex['diff']:.6f} MYA")
        
        print(f"\n  Zero-length nodes with children: {zero_length_nodes_with_children}")
        print(f"  Total child branches collected: {len(child_branch_lengths)}")
        
        if label:
            print("="*80)
        
        return {
            'total_zero_length': same_age_count,
            'zero_length_with_children': zero_length_nodes_with_children,
            'child_branches': child_branch_lengths,
            'examples': same_age_examples
        }
    
    def generate_tree_report(self, output_file, tolerance=0.01):
        """
        Generate a comprehensive report on the tree structure and statistics.
        Saves to a text file.
        
        Parameters
        ----------
        output_file : str
            Path to output report file
        tolerance : float
            Tolerance in MYA for identifying zero-length branches (default 0.01)
        """
        import numpy as np
        from collections import Counter
        
        with open(output_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("PHYLOGENETIC TREE STRUCTURE REPORT\n")
            f.write("="*80 + "\n\n")
            
            # Basic tree statistics
            f.write("TREE STRUCTURE SUMMARY\n")
            f.write("-"*80 + "\n")
            total_nodes = len(self.nodes)
            leaves = [n for n in self.nodes.values() if len(n.children) == 0]
            internal_nodes = [n for n in self.nodes.values() if len(n.children) > 0]
            
            f.write(f"Total nodes:     {total_nodes:>8}\n")
            f.write(f"Leaf nodes:      {len(leaves):>8} (species)\n")
            f.write(f"Internal nodes:  {len(internal_nodes):>8}\n")
            
            # Root information
            root_id = self.find_root()
            if root_id and root_id in self.nodes:
                root_node = self.nodes[root_id]
                root_age = root_node.nodeage if root_node.nodeage else "MISSING"
                f.write(f"\nRoot node:       {root_id}\n")
                if root_node.name:
                    f.write(f"Root name:       {root_node.name}\n")
                if isinstance(root_age, (int, float)):
                    f.write(f"Root age:        {root_age:.4f} MYA\n")
                else:
                    f.write(f"Root age:        {root_age}\n")
            
            # Node degree distribution (number of children)
            f.write("\n" + "="*80 + "\n")
            f.write("NODE DEGREE DISTRIBUTION (Number of Children)\n")
            f.write("="*80 + "\n")
            
            degree_counts = Counter(len(n.children) for n in self.nodes.values())
            
            # Ensure all degrees from 0 to max are included (even if count is 0)
            # This makes it easy to see the distribution and any gaps
            if degree_counts:
                max_degree = max(degree_counts.keys())
                # Create a complete range from 0 to max_degree
                all_degrees = list(range(0, max_degree + 1))
                # Fill in any missing degrees with count 0
                for degree in all_degrees:
                    if degree not in degree_counts:
                        degree_counts[degree] = 0
            else:
                all_degrees = []
            
            # Header with percentage column
            f.write(f"\n{'Degree':<12} {'Count':<10} {'%':<8} {'Type':<25} Distribution\n")
            f.write("-"*80 + "\n")
            
            max_count = max(degree_counts.values()) if degree_counts else 1
            bar_width = 35  # Slightly reduced to fit percentage column
            
            for degree in all_degrees:
                count = degree_counts[degree]
                pct = (count / total_nodes) * 100 if total_nodes > 0 else 0
                bar_length = int((count / max_count) * bar_width) if max_count > 0 else 0
                bar = '█' * bar_length
                
                if degree == 0:
                    node_type = "Species (leaves)"
                elif degree == 1:
                    node_type = "Unbranched internal"
                elif degree == 2:
                    node_type = "Bifurcating"
                else:
                    node_type = f"Polytomy ({degree}-way)"
                
                f.write(f"{degree:<12} {count:<10} {pct:<7.2f} {node_type:<25} {bar}\n")
            
            # Node age distribution - ALL NODES
            f.write("\n" + "="*80 + "\n")
            f.write("NODE AGE DISTRIBUTION (ALL NODES)\n")
            f.write("="*80 + "\n")
            
            ages = [n.nodeage for n in self.nodes.values() if n.nodeage is not None]
            
            if ages:
                ages_array = np.array(ages)
                f.write(f"\nNodes with ages:  {len(ages)} / {total_nodes}\n")
                f.write(f"Min age:          {np.min(ages_array):.4f} MYA\n")
                f.write(f"Max age:          {np.max(ages_array):.4f} MYA\n")
                f.write(f"Mean age:         {np.mean(ages_array):.4f} MYA\n")
                f.write(f"Median age:       {np.median(ages_array):.4f} MYA\n")
                f.write(f"Std dev:          {np.std(ages_array):.4f} MYA\n")
                
                # Age histogram with 10 MYA bins
                f.write("\nAge Distribution Histogram (10 MYA bins):\n")
                f.write(f"\n{'Age Range (MYA)':<25} {'Count':<10} Distribution\n")
                f.write("-"*80 + "\n")
                
                # Create bins with 10 MYA resolution
                max_age = np.max(ages_array)
                bin_width = 10.0  # 10 million years
                num_bins = int(np.ceil(max_age / bin_width))
                bins = np.arange(0, (num_bins + 1) * bin_width, bin_width)
                
                counts, bin_edges = np.histogram(ages_array, bins=bins)
                max_count_hist = max(counts) if len(counts) > 0 else 1
                
                for i in range(len(counts)):
                    if counts[i] > 0:
                        bin_start = bin_edges[i]
                        bin_end = bin_edges[i + 1]
                        count = counts[i]
                        bar_length = int((count / max_count_hist) * bar_width)
                        bar = '█' * bar_length
                        f.write(f"{bin_start:>7.2f} - {bin_end:<7.2f}    {count:<10} {bar}\n")
            else:
                f.write("\nNo nodes with age information found.\n")
            
            # Internal nodes (non-extant) age distribution
            f.write("\n" + "="*80 + "\n")
            f.write("NODE AGE DISTRIBUTION (INTERNAL NODES ONLY - EXCLUDING EXTANT SPECIES)\n")
            f.write("="*80 + "\n")
            
            # Get ages for nodes WITH children (internal nodes, not extant species)
            internal_ages = [n.nodeage for n in self.nodes.values() 
                            if n.nodeage is not None and len(n.children) > 0]
            
            if internal_ages:
                internal_ages_array = np.array(internal_ages)
                f.write(f"\nInternal nodes:   {len(internal_ages)}\n")
                f.write(f"Min age:          {np.min(internal_ages_array):.4f} MYA\n")
                f.write(f"Max age:          {np.max(internal_ages_array):.4f} MYA\n")
                f.write(f"Mean age:         {np.mean(internal_ages_array):.4f} MYA\n")
                f.write(f"Median age:       {np.median(internal_ages_array):.4f} MYA\n")
                f.write(f"Std dev:          {np.std(internal_ages_array):.4f} MYA\n")
                
                f.write("\nInternal Node Age Distribution (10 MYA bins):\n")
                f.write(f"\n{'Age Range (MYA)':<25} {'Count':<10} Distribution\n")
                f.write("-"*80 + "\n")
                
                # Create bins with 10 MYA resolution
                max_internal_age = np.max(internal_ages_array)
                bin_width = 10.0  # 10 million years
                num_bins = int(np.ceil(max_internal_age / bin_width))
                bins = np.arange(0, (num_bins + 1) * bin_width, bin_width)
                
                counts, bin_edges = np.histogram(internal_ages_array, bins=bins)
                max_count_hist = max(counts) if len(counts) > 0 else 1
                
                for i in range(len(counts)):
                    if counts[i] > 0:
                        bin_start = bin_edges[i]
                        bin_end = bin_edges[i + 1]
                        count = counts[i]
                        bar_length = int((count / max_count_hist) * bar_width)
                        bar = '█' * bar_length
                        f.write(f"{bin_start:>7.2f} - {bin_end:<7.2f}    {count:<10} {bar}\n")
            else:
                f.write("\nNo internal nodes with age information found.\n")
            
            # Branch length distribution
            f.write("\n" + "="*80 + "\n")
            f.write("BRANCH LENGTH DISTRIBUTION\n")
            f.write("="*80 + "\n")
            
            branch_lengths = []
            for node in self.nodes.values():
                if node.parent is not None and node.parent in self.nodes:
                    parent_node = self.nodes[node.parent]
                    if node.nodeage is not None and parent_node.nodeage is not None:
                        branch_len = parent_node.nodeage - node.nodeage
                        if branch_len >= 0:  # Only positive branches
                            branch_lengths.append(branch_len)
            
            if branch_lengths:
                branch_array = np.array(branch_lengths)
                f.write(f"\nTotal branches:   {len(branch_lengths)}\n")
                f.write(f"Min length:       {np.min(branch_array):.4f} MYA\n")
                f.write(f"Max length:       {np.max(branch_array):.4f} MYA\n")
                f.write(f"Mean length:      {np.mean(branch_array):.4f} MYA\n")
                f.write(f"Median length:    {np.median(branch_array):.4f} MYA\n")
                f.write(f"Std dev:          {np.std(branch_array):.4f} MYA\n")
                
                # Branch length histogram (log scale for better visualization)
                f.write("\nBranch Length Distribution (log-scaled bins):\n")
                f.write(f"\n{'Length Range (MYA)':<25} {'Count':<10} Distribution\n")
                f.write("-"*80 + "\n")
                
                # Use log-spaced bins if range is large
                min_len = np.min(branch_array[branch_array > 0]) if np.any(branch_array > 0) else 0.001
                max_len = np.max(branch_array)
                
                if max_len / min_len > 100:  # Large range, use log scale
                    bins = np.logspace(np.log10(min_len), np.log10(max_len), 16)
                else:  # Small range, use linear scale
                    bins = 15
                
                counts, bin_edges = np.histogram(branch_array, bins=bins)
                max_count_hist = max(counts)
                
                for i in range(len(counts)):
                    if counts[i] > 0:
                        bin_start = bin_edges[i]
                        bin_end = bin_edges[i + 1]
                        count = counts[i]
                        bar_length = int((count / max_count_hist) * bar_width)
                        bar = '█' * bar_length
                        f.write(f"{bin_start:>8.4f} - {bin_end:<8.4f}  {count:<10} {bar}\n")
            else:
                f.write("\nNo valid branch lengths found.\n")
            
            # Tree depth statistics (distance from root to leaves)
            f.write("\n" + "="*80 + "\n")
            f.write("TREE DEPTH STATISTICS\n")
            f.write("="*80 + "\n")
            
            root_id = self.find_root()
            if root_id and root_id in self.nodes:
                depths = []  # Time from root to each leaf
                
                for leaf in leaves:
                    if leaf.nodeage is not None:
                        root_age = self.nodes[root_id].nodeage
                        if root_age is not None:
                            depth = root_age - leaf.nodeage
                            depths.append(depth)
                
                if depths:
                    depth_array = np.array(depths)
                    f.write(f"\nLeaves analyzed:  {len(depths)}\n")
                    f.write(f"Min depth:        {np.min(depth_array):.4f} MYA\n")
                    f.write(f"Max depth:        {np.max(depth_array):.4f} MYA\n")
                    f.write(f"Mean depth:       {np.mean(depth_array):.4f} MYA\n")
                    f.write(f"Median depth:     {np.median(depth_array):.4f} MYA\n")
                    f.write(f"Std dev:          {np.std(depth_array):.4f} MYA\n")
                    f.write(f"\n(Depth = time from root to leaf)\n")
                else:
                    f.write("\nCould not calculate depths - missing age information.\n")
            else:
                f.write("\nCould not find root node.\n")
            
            # Zero-length branch analysis
            f.write("\n" + "="*80 + "\n")
            f.write("ZERO-LENGTH BRANCH ANALYSIS\n")
            f.write("="*80 + "\n")
            
            zero_length_count = 0
            zero_length_with_children = 0
            
            for node in self.nodes.values():
                if node.parent is not None and node.parent in self.nodes:
                    parent_node = self.nodes[node.parent]
                    if node.nodeage is not None and parent_node.nodeage is not None:
                        diff = abs(parent_node.nodeage - node.nodeage)
                        if diff <= tolerance:
                            zero_length_count += 1
                            if len(node.children) > 0:
                                zero_length_with_children += 1
            
            f.write(f"\nTolerance:                    {tolerance} MYA\n")
            f.write(f"Zero-length branches:         {zero_length_count}\n")
            f.write(f"  - with children:            {zero_length_with_children}\n")
            f.write(f"  - without children (tips):  {zero_length_count - zero_length_with_children}\n")
            
            if zero_length_count > 0:
                pct = (zero_length_count / len([n for n in self.nodes.values() if n.parent is not None])) * 100
                f.write(f"Percentage of branches:       {pct:.2f}%\n")
            
            # Locked nodes (priority calibrations)
            f.write("\n" + "="*80 + "\n")
            f.write("CALIBRATION INFORMATION\n")
            f.write("="*80 + "\n")
            
            locked_nodes = [n for n in self.nodes.values() if n.lock_age]
            f.write(f"\nLocked nodes (priority calibrations): {len(locked_nodes)}\n")
            
            if locked_nodes:
                f.write("\nLocked node examples (first 10):\n")
                for i, node in enumerate(locked_nodes[:10]):
                    node_name = node.name if node.name else f"Node{node.taxid}"
                    node_age = node.nodeage if node.nodeage else "MISSING"
                    if isinstance(node_age, (int, float)):
                        f.write(f"  {i+1}. {node_name} (taxid {node.taxid}): {node_age:.4f} MYA\n")
                    else:
                        f.write(f"  {i+1}. {node_name} (taxid {node.taxid}): {node_age}\n")
            
            f.write("\n" + "="*80 + "\n")
            f.write("END OF REPORT\n")
            f.write("="*80 + "\n")
        
        print(f"Tree structure report saved to: {output_file}")

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

    def write_newick(self, output_file):
        """
        Write the tree with calibrated node ages to newick format.
        Branch lengths are set to the time difference between parent and child nodes.
        Node names include the clade/species name with taxid in brackets: Name[taxid]
        Uses ete4 PhyloTree.
        """
        # Create root node (ete4 style: create empty node then set name)
        ete_tree = PhyloTree()
        # Format: Name[taxid]
        root_name = self.nodes[self.root].name if self.nodes[self.root].name else f"Node{self.root}"
        ete_tree.name = f"{root_name}[{self.root}]"
        node_map = {self.root: ete_tree}
        
        # Build the tree recursively using BFS to ensure proper structure
        queue = [self.root]
        while queue:
            taxid = queue.pop(0)
            ete_node = node_map[taxid]
            parent_age = self.nodes[taxid].nodeage
            
            for child_taxid in self.nodes[taxid].children:
                child_age = self.nodes[child_taxid].nodeage
                # Branch length is the time difference (parent age - child age)
                branch_length = parent_age - child_age if (parent_age is not None and child_age is not None) else 0
                
                # Add child node with branch length (distance)
                child_ete_node = ete_node.add_child(dist=branch_length)
                # Format node name as: Name[taxid]
                child_name = self.nodes[child_taxid].name if self.nodes[child_taxid].name else f"Node{child_taxid}"
                # Remove single quotes if present in name
                child_name = child_name.replace("'", "")
                child_ete_node.name = f"{child_name}[{child_taxid}]"
                node_map[child_taxid] = child_ete_node
                queue.append(child_taxid)
        
        # Write to file in newick format with branch lengths
        # In ete4, parser=1 for write() includes internal node names in the output
        with open(output_file, 'w') as f:
            f.write(ete_tree.write(parser=1, format_root_node=True))
        print(f"Wrote calibrated tree to: {output_file}")

    def print_edge_information(self, outfile) -> None:
        """
        Prints the edge information of the tree. Use __slots__ to determine what to print.
        """
        # get the fields. Don't include fusions or losses
        fields = [x for x in self.edges[(self.root, list(self.nodes[self.root].children)[0])].__slots__
                  if x not in ["fusions", "losses"]]
        fusion_fields = [x for x in self.edges[(self.root, list(self.nodes[self.root].children)[0])].fusions]
        losses_fields = [x for x in self.edges[(self.root, list(self.nodes[self.root].children)[0])].losses]
        with open(outfile, "w") as outhandle:
            print("\t".join(fields + fusion_fields + losses_fields), file = outhandle)
            for edge in self.edges:
                outstring = "\t".join([str(getattr(self.edges[edge], x))
                                 for x in self.edges[edge].__slots__])
                if len(fusion_fields) > 0:
                    outstring += "\t" + "\t".join([str(self.edges[edge].fusions[x]) for x in fusion_fields])
                if len(losses_fields) > 0:
                    outstring += "\t" + "\t".join([str(self.edges[edge].losses[x]) for x in losses_fields])
                print(outstring, file = outhandle)

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
    NCBI = NCBITaxa()
    
    ## Load the TimeTree.org calibrated newick file (has divergence times)
    print(f"\nLoading TimeTree.org calibrated tree from: {args.time_newick}")
    timetree = PhyloTree(open(args.time_newick).read(), parser=1)
    
    # Extract the TimeTree root age (Metazoa for metazoan-only trees)
    metazoa_age_from_timetree = extract_timetree_root_age_as_metazoa(timetree)
    
    # Load the custom topology tree
    print(f"Loading custom topology tree from: {args.topology_newick}")
    tree = PhyloTree(open(args.topology_newick).read(), parser=1)
    print("  Using custom topology with TimeTree divergence times")
    
    # DIAGNOSTIC: Check for None leaf names
    all_leaves = list(tree.leaves())
    none_leaves = [leaf for leaf in all_leaves if leaf.name is None]
    if none_leaves:
        print(f"\nWARNING: Found {len(none_leaves)} leaves with None names out of {len(all_leaves)} total leaves!")
        print(f"\nINVESTIGATING None-named leaves:")
        print(f"{'='*80}")
        
        for i, leaf in enumerate(none_leaves):
            print(f"\n--- None Leaf #{i+1} ---")
            print(f"  name: {leaf.name}")
            print(f"  dist (branch length): {leaf.dist if hasattr(leaf, 'dist') else 'N/A'}")
            print(f"  support: {leaf.support if hasattr(leaf, 'support') else 'N/A'}")
            print(f"  is_leaf: {leaf.is_leaf}")
            
            # Check parent
            if hasattr(leaf, 'up') and leaf.up:
                parent = leaf.up
                print(f"  parent.name: {parent.name if hasattr(parent, 'name') else 'N/A'}")
                print(f"  parent.dist: {parent.dist if hasattr(parent, 'dist') else 'N/A'}")
                
                # Check siblings
                siblings = [child for child in parent.children if child != leaf]
                print(f"  num_siblings: {len(siblings)}")
                if siblings:
                    print(f"  sibling names: {[s.name for s in siblings[:3]]}")
            else:
                print(f"  parent: None (root?)")
            
            # Check if it has any features/attributes
            if hasattr(leaf, 'features'):
                print(f"  features: {leaf.features}")
            
            # Try to get newick representation
            try:
                newick_str = leaf.write(format=1)
                print(f"  newick (first 100 chars): {newick_str[:100]}")
            except:
                print(f"  newick: <error writing newick>")
        
        print(f"\n{'='*80}")
        print(f"\nPOSSIBLE CAUSES:")
        print(f"   1. Empty nodes from malformed Newick: (,) or ():")
        print(f"   2. Internal nodes mistakenly parsed as leaves")
        print(f"   3. Nodes with only whitespace names that got stripped")
        print(f"   4. Artifacts from tree manipulation in taxids_to_newick.py")
        print(f"\nNEXT STEPS:")
        print(f"   - Check the Newick file for empty leaf labels")
        print(f"   - Look at parent/sibling context to locate position in tree")
        print(f"   - These nodes are being skipped (won't break processing)")
    
    # Annotate the custom tree with times from TimeTree
    timetree_divergences, timetree_name_to_taxid, custom_name_to_taxid = \
        annotate_custom_tree_with_timetree_ages(tree, timetree, NCBI)
    
    print(f"\nUsing custom topology tree with {len(list(tree.leaves()))} species")
    
    # CRITICAL: Build CFtree from custom topology FIRST, before extracting lineages
    # This preserves your custom phylogeny (Myriazoa, Ctenophora placement, etc.)
    print("\nBuilding custom phylogeny tree structure...")
    CFtree = TaxIDtree()
    CFtree.build_from_newick_tree(tree, NCBI)
    print(f"  Custom tree built with {len(CFtree.nodes)} nodes")
    
    # Add lineage information to CFtree based on custom topology
    CFtree.add_lineage_info()
    print("  Added lineage information based on custom topology")
    
    entries = []
    skipped_species = []
    # iterate through the leaves
    all_leaves_list = list(tree.leaves())
    print(f"\nProcessing {len(all_leaves_list)} species to extract taxids and custom lineages...")
    for leaf_idx, thisleaf in enumerate(all_leaves_list):
        if leaf_idx % 500 == 0 and leaf_idx > 0:
            print(f"  Processed {leaf_idx}/{len(all_leaves_list)} species...", end="\r")
        # Skip leaves with None names
        if thisleaf.name is None:
            print(f"  Warning: Skipping leaf with None name")
            skipped_species.append("<None>")
            continue
        # Parse taxid from node name if it's in brackets, e.g., "Homo_sapiens[9606]"
        if '[' in thisleaf.name and thisleaf.name.endswith(']'):
            # Extract taxid from brackets
            species_name = thisleaf.name[:thisleaf.name.rfind('[')]
            taxid = int(thisleaf.name[thisleaf.name.rfind('[')+1:-1])
            taxname = species_name.replace("_", " ")
        else:
            # Fallback to old method if no brackets found
            taxname = thisleaf.name.replace("_", " ")
            try:
                taxid = NCBI.get_name_translator([taxname])[taxname][0]
            except (KeyError, IndexError):
                skipped_species.append(thisleaf.name)
                print(f"  Warning: Could not find taxid for species: {thisleaf.name}, skipping...")
                continue
        
        try:
            # Extract lineage from custom tree instead of NCBI to preserve custom topology
            if taxid in CFtree.nodes:
                custom_lineage = CFtree.get_lineage(taxid)
            else:
                # Fallback: if taxid not in custom tree, use NCBI lineage
                print(f"  Warning: taxid {taxid} ({thisleaf.name}) not in custom tree, using NCBI lineage")
                custom_lineage = NCBI.get_lineage(taxid)
            
            # Store node name in clean format Name[TaxID] for consistency
            # If the name already has brackets (from parser=1), use as-is
            # Otherwise, construct the bracketed format
            if '[' in thisleaf.name and thisleaf.name.endswith(']'):
                clean_nodename = thisleaf.name  # Already in Name[TaxID] format
            else:
                clean_nodename = f"{thisleaf.name}[{taxid}]"  # Construct Name[TaxID] format
            
            entry = {"node": thisleaf,
                     "nodename": clean_nodename,
                     "taxname": taxname,
                     "taxid": taxid,
                     "lineage": custom_lineage}
            entries.append(entry)
        except Exception as e:
            skipped_species.append(thisleaf.name)
            print(f"  Warning: Could not process species: {thisleaf.name}, error: {e}")
            continue
    
    print(f"  Processed {len(all_leaves_list)}/{len(all_leaves_list)} species... Done!")
    if skipped_species:
        print(f"\nSkipped {len(skipped_species)} species due to errors")
        print(f"  Total species processed: {len(entries)}")

    # the ttdf is the timetree df
    ttdf = pd.DataFrame(entries)
    ttdf_taxid_to_sp = {row["taxid"]: row["nodename"] for i, row in ttdf.iterrows()}
    print(ttdf)
    
    # Show example lineage to verify custom topology is preserved
    if len(entries) > 0:
        example_entry = entries[0]
        # Extract display name from Name[TaxID] format for cleaner output
        display_name = example_entry['nodename']
        if '[' in display_name and display_name.endswith(']'):
            display_name = display_name[:display_name.rfind('[')].replace('_', ' ')
        print(f"\nExample lineage for {display_name} (custom topology):")
        print(f"  Taxid path: {example_entry['lineage']}")
        # Show the names for the lineage nodes in Name[TaxID] format
        lineage_names = []
        for tid in example_entry['lineage']:
            if tid in CFtree.nodes:
                node_name = CFtree.nodes[tid].name if CFtree.nodes[tid].name else f"Node{tid}"
                lineage_names.append(f"{node_name}[{tid}]")
            else:
                lineage_names.append(f"Unknown[{tid}]")
        print(f"  Names: {' -> '.join(lineage_names)}")

    # Build TTtree from custom topology for representative matching
    # This uses the same species set as the config, ensuring find_closest_relative() can find matches
    print("\nBuilding TTtree from custom topology for representative matching...")
    TTtree = TaxIDtree()
    for i, row in ttdf.iterrows():
        lineage = row["lineage"]
        for j in range(len(lineage)-1):
            TTtree.add_edge(lineage[j], lineage[j+1])
    print(f"  TTtree built with {len(TTtree.nodes)} nodes from custom topology")

    # Note: We use timetree_divergences dict (taxid-based) instead of species-name-based lookups
    # This ensures all TimeTree calibration ages are properly added to node Counters

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
    print("\nProcessing config file species with custom topology...")
    subspecies_converted = 0
    for sp in config["species"]:
        original_taxid = config["species"][sp]["taxid"]
        taxid = original_taxid
        
        # Convert subspecies to species-level taxid to match tree
        if is_subspecies_or_below(taxid, NCBI):
            taxid = get_species_level_taxid(taxid, NCBI)
            if taxid != original_taxid:
                subspecies_converted += 1
                print(f"  Converted subspecies {original_taxid} -> species {taxid} for {sp}")
        
        # get the taxname using NCBI
        spname = NCBI.get_taxid_translator([taxid])[taxid]
        # Use custom tree lineage instead of NCBI lineage
        if taxid in CFtree.nodes:
            lineage = CFtree.get_lineage(taxid)
        else:
            print(f"  Warning: Config species {sp} (taxid {taxid}) not in custom tree, using NCBI lineage")
            lineage = NCBI.get_lineage(taxid)
        # get the closest species in timetree
        tt_representative = TTtree.find_closest_relative(NCBI, taxid)
        entry = {"node": None,
                 "nodename": sp,
                 "taxname": spname,
                 "taxid": taxid,
                 "taxid_current": lineage[-1],
                 "taxid_in_timetree": tt_representative,
                 "lineage": lineage}
        entries.append(entry)
        if tt_representative is None:
            print(f"  Warning: No TimeTree representative found for {sp} (taxid {taxid})")
    configdf = pd.DataFrame(entries)
    print(f"  Processed {len(entries)} config species")
    if subspecies_converted > 0:
        print(f"  Converted {subspecies_converted} subspecies entries to species-level")

    # CFtree was already built earlier with custom topology
    # Just set the leaf ages to zero
    CFtree.set_leaf_ages_to_zero()

    outfile = f"{args.prefix}.node_ages_for_config.tsv"
    if not os.path.exists(outfile):
        # if the weights file does not yet exist, we need to make it.
        # This step takes the longest of anything in the paper, so it is better to load if possible.
        # Now go through and add all the ages to the config file tree
        
        # DIAGNOSTIC: Check how many config species have TimeTree representatives
        config_tt_taxids = set(configdf["taxid_in_timetree"])
        # Remove None from the set if present
        config_tt_taxids.discard(None)
        timetree_taxids = set()
        for (t1, t2) in timetree_divergences.keys():
            timetree_taxids.add(t1)
            timetree_taxids.add(t2)
        
        config_in_timetree = config_tt_taxids.intersection(timetree_taxids)
        num_none_representatives = (configdf["taxid_in_timetree"] == None).sum()
        print(f"\n    TimeTree divergences dict has {len(timetree_divergences)//2} pairwise ages covering {len(timetree_taxids)} unique taxids")
        print(f"    Config has {len(config_tt_taxids)} unique TimeTree representative taxids (excluding {num_none_representatives} with no representative)")
        print(f"    {len(config_in_timetree)} config representatives are actually in TimeTree")
        print(f"    {len(config_tt_taxids - config_in_timetree)} config representatives are NOT in TimeTree (will have no ages)")
        
        # Show a few examples of missing taxids
        if len(config_tt_taxids - config_in_timetree) > 0:
            missing_examples = list(config_tt_taxids - config_in_timetree)[:10]
            print(f"    Example missing taxids: {missing_examples}")
            # Show what those taxids map to
            for taxid in missing_examples[:5]:
                try:
                    name = NCBI.get_taxid_translator([taxid])[taxid]
                    print(f"      {taxid}: {name}")
                except:
                    print(f"      {taxid}: UNKNOWN")
        
        successful_pairs = 0
        missing_pairs = 0
        for i in range(len(configdf["taxid"])-1):
            # make a print statement that stays on one line
            print("    - Iterating through species, on sample {}/{}".format(i, len(configdf["taxid"])), end="\r")
            for j in range(i+1, len(configdf["taxid"])):
                taxid1_config = configdf["taxid_current"][i]
                taxid2_config = configdf["taxid_current"][j]
                
                # Skip if either taxid is not in CFtree (e.g., config species not in topology)
                if taxid1_config not in CFtree.nodes or taxid2_config not in CFtree.nodes:
                    continue
                
                taxid1_tt     = configdf["taxid_in_timetree"][i]
                taxid2_tt     = configdf["taxid_in_timetree"][j]
                #print("Comparing {} and {}".format(taxid1_tt, taxid2_tt))
                #print("  - config taxids are {} and {}".format(taxid1_config, taxid2_config))
                
                # Skip if either representative is None
                if taxid1_tt is None or taxid2_tt is None:
                    missing_pairs += 1
                    continue
                
                if taxid1_tt != taxid2_tt:
                    # Use taxid-based lookup from timetree_divergences (both directions)
                    div_time = None
                    if (taxid1_tt, taxid2_tt) in timetree_divergences:
                        div_time = timetree_divergences[(taxid1_tt, taxid2_tt)]
                    elif (taxid2_tt, taxid1_tt) in timetree_divergences:
                        div_time = timetree_divergences[(taxid2_tt, taxid1_tt)]
                    
                    if div_time is not None:
                        common_ancestor = CFtree.find_LCA(taxid1_config, taxid2_config)
                        # add the node age to the node in CFtree
                        CFtree.nodes[common_ancestor].nodeages.update([div_time])
                        successful_pairs += 1
                    else:
                        # Pair not in TimeTree - will be interpolated later
                        missing_pairs += 1
        
        print(f"\n    - Collected {successful_pairs} age calibrations from TimeTree")
        print(f"    - {missing_pairs} species pairs not in TimeTree (will be interpolated)")
        
        # Diagnostic: Show how many nodes got ages
        nodes_with_ages = sum(1 for node in CFtree.nodes if len(CFtree.nodes[node].nodeages) > 0)
        print(f"    - {nodes_with_ages}/{len(CFtree.nodes)} nodes have age estimates after pairwise comparison")
        
        # Show a few examples of nodes with ages
        print(f"\n    Example nodes with ages:")
        count = 0
        for node_id in CFtree.nodes:
            if len(CFtree.nodes[node_id].nodeages) > 0 and count < 5:
                node = CFtree.nodes[node_id]
                print(f"      Node {node_id} ({node.name}): {len(node.nodeages)} age estimates, most common = {node.nodeages.most_common(1)[0] if node.nodeages else 'N/A'}")
                count += 1

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

    # Apply the TimeTree root age to the custom phylogeny root AFTER loading from file
    # This prevents the root from getting the 3700 MYA fallback for "root of all life"
    # Must be done AFTER file loading (if file exists) to avoid being overwritten
    if metazoa_age_from_timetree is not None:
        root_node_id = CFtree.find_root()
        CFtree.nodes[root_node_id].nodeages.update([metazoa_age_from_timetree])
        CFtree.nodes[root_node_id].nodeage = metazoa_age_from_timetree
        print(f"  Applied TimeTree root age ({metazoa_age_from_timetree:.2f} MYA) to CFtree root node {root_node_id}")
        print(f"    Root nodeages Counter: {CFtree.nodes[root_node_id].nodeages}")
        print(f"    Root nodeage value: {CFtree.nodes[root_node_id].nodeage}")

    # Extract priority node ages from TimeTree for major clades
    # These will be applied before general interpolation to anchor the tree
    print("\nExtracting priority node ages from TimeTree...")
    priority_node_ages = {}
    
    # 1. BILATERIA (critical anchor point) - Use human-fly divergence
    human_taxid = 9606
    fly_taxid = 7227
    
    # Check if human and fly are in the CFtree
    if human_taxid in CFtree.nodes and fly_taxid in CFtree.nodes:
        # Find the Bilateria node (LCA of human and fly)
        bilateria_node = CFtree.find_LCA(human_taxid, fly_taxid)
        
        # The TimeTree ages were already applied to CFtree nodes during the pairwise comparison loop
        # Extract the age from the node's nodeages Counter
        if bilateria_node in CFtree.nodes and len(CFtree.nodes[bilateria_node].nodeages) > 0:
            # Use the most common age that was already assigned from TimeTree
            bilateria_age = CFtree.nodes[bilateria_node].nodeages.most_common(1)[0][0]
            priority_node_ages[bilateria_node] = (bilateria_age, f"Bilateria (human taxid {human_taxid} - fly taxid {fly_taxid})")
            print(f"  Bilateria node {bilateria_node}: {bilateria_age:.2f} MYA (from existing TimeTree calibration)")
        else:
            print(f"  WARNING: Bilateria node {bilateria_node} has no TimeTree ages, will be interpolated")
    else:
        print(f"  WARNING: Human (taxid {human_taxid}) or Fly (taxid {fly_taxid}) not found in CFtree, Bilateria will be interpolated")
    
    # 2. NON-BILATERIAN PHYLA (secondary calibration points)
    # Cnidaria (6073), Porifera (6040), Ctenophora (10197), Placozoa (10226)
    phyla_to_calibrate = {
        6073: "Cnidaria",
        6040: "Porifera",
        10197: "Ctenophora",
        10226: "Placozoa"
    }
    
    for phylum_taxid, phylum_name in phyla_to_calibrate.items():
        # Find species in this phylum
        phylum_species = []
        for node_id in CFtree.nodes:
            if len(CFtree.nodes[node_id].children) == 0:  # Is a leaf
                lineage = CFtree.get_lineage(node_id)
                if phylum_taxid in lineage:
                    phylum_species.append(node_id)
        
        if len(phylum_species) >= 2:
            # Take first two available species and find their LCA
            sp1_taxid = phylum_species[0]
            sp2_taxid = phylum_species[1]
            phylum_node = CFtree.find_LCA(sp1_taxid, sp2_taxid)
            
            # Extract age from the node's nodeages Counter (already populated from TimeTree)
            if phylum_node in CFtree.nodes and len(CFtree.nodes[phylum_node].nodeages) > 0:
                phylum_age = CFtree.nodes[phylum_node].nodeages.most_common(1)[0][0]
                priority_node_ages[phylum_node] = (phylum_age, f"{phylum_name} (taxid {phylum_taxid})")
                print(f"  {phylum_name} node {phylum_node}: {phylum_age:.2f} MYA (from existing TimeTree calibration)")
            else:
                print(f"  {phylum_name}: Node {phylum_node} has no TimeTree ages, will be interpolated")
        else:
            print(f"  {phylum_name}: Less than 2 species found in tree, will be interpolated")
    
    # 3. MYRIAZOA (custom clade) - Set to 10 MYA younger than root
    # Myriazoa has custom taxid -67 in the tree
    myriazoa_taxid = -67
    if myriazoa_taxid in CFtree.nodes and metazoa_age_from_timetree is not None:
        myriazoa_age = metazoa_age_from_timetree - 10.0  # 10 MYA younger than root
        priority_node_ages[myriazoa_taxid] = (myriazoa_age, f"Myriazoa (10 MYA after Ctenophora split)")
        print(f"  Myriazoa node {myriazoa_taxid}: {myriazoa_age:.2f} MYA (root - 10 MYA)")
    else:
        if myriazoa_taxid not in CFtree.nodes:
            print(f"  WARNING: Myriazoa node (taxid {myriazoa_taxid}) not found in tree")
        else:
            print(f"  WARNING: Cannot set Myriazoa age - no root age available")
    
    print(f"\nTotal priority nodes for calibration: {len(priority_node_ages)}")

    # now we correct the missing nodes
    print("\nNow we are correcting the missing nodes")
    CFtree.correct_missing_nodes(priority_node_ages=priority_node_ages)
    CFtree.ensure_all_leaves_have_age_zero()
    
    # Analyze BEFORE fixing
    stats_before = CFtree.analyze_zero_length_branches(tolerance=0.01, label="BEFORE FIXING")
    
    # Fix zero-length branches
    CFtree.fix_zero_length_branches(tolerance=0.01, max_iterations=5)
    
    # Analyze AFTER fixing
    stats_after = CFtree.analyze_zero_length_branches(tolerance=0.01, label="AFTER FIXING")
    
    # Print improvement summary
    print(f"\nIMPROVEMENT SUMMARY:")
    print(f"  Total zero-length branches: {stats_before['total_zero_length']} → {stats_after['total_zero_length']} "
          f"({stats_before['total_zero_length'] - stats_after['total_zero_length']} fixed)")
    print(f"  Zero-length nodes with children: {stats_before['zero_length_with_children']} → {stats_after['zero_length_with_children']} "
          f"({stats_before['zero_length_with_children'] - stats_after['zero_length_with_children']} fixed)")
    print()
    
    # calculate the dist crown
    CFtree.calc_dist_crown()
    # calculate the lineage info
    CFtree.add_lineage_info()
    # add chromosome information
    if args.chromosome_sizes:
        CFtree.add_chromosome_info_file(args.chromosome_sizes)
    else:
        print("No chromosome sizes file provided, skipping chromosome info")

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
    
    # Check for nodes with same age as their parent (within tolerance)
    print("\nChecking for nodes with essentially same age as parent...")
    tolerance = 0.01  # 0.01 MYA = 10,000 years
    same_age_count = 0
    same_age_examples = []
    
    for node_id in CFtree.nodes:
        parent_id = CFtree.nodes[node_id].parent
        if parent_id is not None:  # Skip root
            node_age = CFtree.nodes[node_id].nodeage
            parent_age = CFtree.nodes[parent_id].nodeage
            
            if node_age is not None and parent_age is not None:
                age_diff = abs(parent_age - node_age)
                if age_diff < tolerance:
                    same_age_count += 1
                    # Collect first 10 examples
                    if len(same_age_examples) < 10:
                        node_name = CFtree.nodes[node_id].name if CFtree.nodes[node_id].name else f"Node{node_id}"
                        parent_name = CFtree.nodes[parent_id].name if CFtree.nodes[parent_id].name else f"Node{parent_id}"
                        same_age_examples.append({
                            'node_id': node_id,
                            'node_name': node_name,
                            'node_age': node_age,
                            'parent_id': parent_id,
                            'parent_name': parent_name,
                            'parent_age': parent_age,
                            'diff': age_diff
                        })
    
    # Note: The detailed zero-length branch analysis is now done earlier
    # using the analyze_zero_length_branches() method before and after fixing
    
    # Generate comprehensive tree structure report
    print("\n" + "="*80)
    print("GENERATING TREE STRUCTURE REPORT")
    print("="*80)
    CFtree.generate_tree_report(f"{args.prefix}.tree_report.txt", tolerance=0.01)

    # now we print the edge information
    CFtree.print_edge_information(f"{args.prefix}.edge_information.tsv")
    CFtree.print_node_information(f"{args.prefix}.node_information.tsv")
    
    # Generate pairwise divergence times for all species
    print("\nGenerating pairwise divergence times file...")
    report_divergence_time_all_vs_all(CFtree, args.prefix)
    
    # Write the calibrated tree to newick format
    CFtree.write_newick(f"{args.prefix}.calibrated_tree.nwk")
    
    # Print final verification of calibrated ages for key nodes
    print("\n" + "="*80)
    print("FINAL CALIBRATED NODE AGES (verification):")
    print("="*80)
    
    # Show root (Metazoa) age
    root_node_id = CFtree.find_root()
    root_name = CFtree.nodes[root_node_id].name if CFtree.nodes[root_node_id].name else f"Node{root_node_id}"
    root_age = CFtree.nodes[root_node_id].nodeage if CFtree.nodes[root_node_id].nodeage else "MISSING"
    print(f"  Metazoa root [{root_name}] (node {root_node_id}): {root_age:.2f} MYA" if isinstance(root_age, (int, float)) else f"  Metazoa root [{root_name}] (node {root_node_id}): {root_age}")
    
    # Show Bilateria age if it was calibrated
    human_taxid = 9606
    fly_taxid = 7227
    if human_taxid in ttdf_taxid_to_sp and fly_taxid in ttdf_taxid_to_sp:
        bilateria_node = CFtree.find_LCA(human_taxid, fly_taxid)
        if bilateria_node in CFtree.nodes:
            bil_name = CFtree.nodes[bilateria_node].name if CFtree.nodes[bilateria_node].name else f"Node{bilateria_node}"
            bil_age = CFtree.nodes[bilateria_node].nodeage if CFtree.nodes[bilateria_node].nodeage else "MISSING"
            print(f"  Bilateria [{bil_name}] (node {bilateria_node}): {bil_age:.2f} MYA" if isinstance(bil_age, (int, float)) else f"  Bilateria [{bil_name}] (node {bilateria_node}): {bil_age}")
    
    # Show non-bilaterian phyla if they were calibrated
    phyla_to_check = {
        6073: "Cnidaria",
        6040: "Porifera",
        10197: "Ctenophora",
        10226: "Placozoa"
    }
    
    for phylum_taxid, phylum_name in phyla_to_check.items():
        # Find species in this phylum
        phylum_species = []
        for node_id in CFtree.nodes:
            if len(CFtree.nodes[node_id].children) == 0:  # Is a leaf
                lineage = CFtree.get_lineage(node_id)
                if phylum_taxid in lineage:
                    phylum_species.append(node_id)
        
        if len(phylum_species) >= 2:
            sp1_taxid = phylum_species[0]
            sp2_taxid = phylum_species[1]
            if sp1_taxid in CFtree.nodes and sp2_taxid in CFtree.nodes:
                phylum_node = CFtree.find_LCA(sp1_taxid, sp2_taxid)
                if phylum_node in CFtree.nodes:
                    ph_name = CFtree.nodes[phylum_node].name if CFtree.nodes[phylum_node].name else f"Node{phylum_node}"
                    ph_age = CFtree.nodes[phylum_node].nodeage if CFtree.nodes[phylum_node].nodeage else "MISSING"
                    print(f"  {phylum_name} [{ph_name}] (node {phylum_node}): {ph_age:.2f} MYA" if isinstance(ph_age, (int, float)) else f"  {phylum_name} [{ph_name}] (node {phylum_node}): {ph_age}")
    
    print("="*80)
    print("\nsuccess")

if __name__ == "__main__":
    main()