#!/usr/bin/env python

"""
This script takes in a list of NCBI taxids and writes a Newick file.

PURPOSE:
========
Generates a phylogenetic tree in Newick format from NCBI taxonomy that can be
visualized in software like FigTree or used for downstream phylogenetic analyses.

USAGE:
======
Basic usage with taxid file (NCBI topology):
    python taxids_to_newick.py -t taxid_list.txt -o ncbi_tree.nwk

Basic usage with config file (NCBI topology):
    python taxids_to_newick.py -c config.yaml -o ncbi_tree.nwk

With custom Ctenophora placement:
    python taxids_to_newick.py -t taxid_list.txt -o ncbi_tree.nwk --custom_phylogeny
    python taxids_to_newick.py -c config.yaml -o ncbi_tree.nwk --custom_phylogeny

INPUT:
======
Either (but not both):
- Text file with one NCBI TaxID per line (integers only), OR
- YAML config file with species entries containing 'taxid' fields
- Requires initialized ete4 NCBI taxonomy database

OUTPUT:
=======
- Newick format tree file (.nwk)
- Species names replace TaxIDs at tips (spaces → underscores)

CUSTOM PHYLOGENY:
=================
When --custom_phylogeny flag is used, Ctenophora (10197) is placed as sister
to all other animals, rather than nested within Eumetazoa as in NCBI taxonomy.

This reflects the phylogenomic hypothesis from Schultz et al. (2023) Nature:
https://doi.org/10.1038/s41586-023-05936-6

The modified topology is:
  Metazoa (33208)
  ├─ Ctenophora (10197)
  └─ Myriazoa (-67) [CUSTOM NODE]
     ├─ Porifera (6040)
     └─ Eumetazoa (6072)
        └─ [Cnidaria, Bilateria, etc.]

Note: Myriazoa is represented by fake taxid -67 (negative to avoid conflicts).

REQUIREMENTS:
=============
- ete4 (updated from ete3)
- NCBI taxonomy database initialized

To initialize NCBI taxonomy (first time only):
    python -c "from ete4 import NCBITaxa; ncbi = NCBITaxa(); ncbi.update_taxonomy_database()"

AUTHORS: Darrin T. Schultz
DATE: 2023-2025
"""

import argparse
import sys
import yaml
from ete4 import NCBITaxa, PhyloTree

def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a Newick tree from a list of NCBI TaxIDs or from a config file",
        epilog="Example: python taxids_to_newick.py -t taxids.txt -o tree.nwk --custom_phylogeny\n"
               "         python taxids_to_newick.py -c config.yaml -o tree.nwk")
    
    # Create mutually exclusive group for input source
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-t', '--taxid_file', type=str,
                             help='File containing one NCBI TaxID per line')
    input_group.add_argument('-c', '--config_file', type=str,
                             help='YAML config file with species entries containing taxid fields')
    
    parser.add_argument('-o', '--output_file', type=str, default='ncbi_tree.nwk', 
                        help='Output Newick file (default: ncbi_tree.nwk)')
    parser.add_argument('--custom_phylogeny', action='store_true',
                        help='Use custom phylogeny with Ctenophora as sister to all other animals (Myriazoa=-67)')
    return parser.parse_args()


def read_taxids_from_config(config_file):
    """
    Read taxids from a YAML config file.
    
    Expected format:
    ----------------
    species:
      SpeciesName-taxid-assembly:
        taxid: 12345
        genus: Genus
        species: species
        ...
    
    Parameters:
    -----------
    config_file : str
        Path to YAML config file
        
    Returns:
    --------
    set : Set of taxids (integers) extracted from config file
    """
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Config file not found: {config_file}")
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing YAML config file: {e}")
    
    if 'species' not in config:
        raise ValueError("Config file must contain a 'species' section")
    
    taxids = set()
    for species_name, species_data in config['species'].items():
        if 'taxid' not in species_data:
            print(f"Warning: No 'taxid' field found for species: {species_name}, skipping...")
            continue
        
        try:
            taxid = int(species_data['taxid'])
            taxids.add(taxid)
        except (ValueError, TypeError):
            print(f"Warning: Invalid taxid for species {species_name}: {species_data.get('taxid')}, skipping...")
            continue
    
    if len(taxids) == 0:
        raise ValueError("No valid taxids found in config file")
    
    print(f"Extracted {len(taxids)} unique taxids from config file")
    return taxids


def build_custom_topology_tree(taxids, ncbi):
    """
    Builds a tree with custom topology by creating subtrees and stitching them together.
    
    This implements the custom topology:
      Metazoa (33208)
      ├─   Ctenophora (10197)
      └─   Myriazoa (-67) [Porifera + Eumetazoa]
       ├─  Porifera (6040)
       └─  Parahoxozoa (-68) [[Cnidaria and Placozoa], Bilateria]
        ├─ Cnidaria (6073) and Placozoa (10226)
        └─ Bilateria (33213)

    In other words, the above in Newick format:
    (Ctenophora,(Porifera,((Cnidaria,Placozoa),Bilateria)Parahoxozoa)Myriazoa)Metazoa;

    Parameters:
    -----------
    taxids : list
        List of taxids to include in tree
    ncbi : NCBITaxa
        NCBI taxonomy database object
        
    Returns:
    --------
    ete4.Tree : Tree with custom topology
    """
    METAZOA_TAXID    = 33208
    CTENOPHORA_TAXID = 10197
    PORIFERA_TAXID   = 6040
    CNIDARIA_TAXID   = 6073
    PLACOZOA_TAXID   = 10226
    BILATERIA_TAXID  = 33213
    MYRIAZOA_TAXID   = -67
    PARAHOXOZOA_TAXID = -68
    
    print("Building custom topology tree by assembling clades...")
    
    # Categorize all taxids into their major clades
    ctenophora_taxa = []
    porifera_taxa = []
    cnidaria_taxa = []
    placozoa_taxa = []
    bilateria_taxa = []
    other_taxa = []  # Non-animal taxa
    
    for taxid in taxids:
        lineage = ncbi.get_lineage(taxid)
        
        if CTENOPHORA_TAXID in lineage:
            ctenophora_taxa.append(taxid)
        elif PORIFERA_TAXID in lineage:
            porifera_taxa.append(taxid)
        elif CNIDARIA_TAXID in lineage:
            cnidaria_taxa.append(taxid)
        elif PLACOZOA_TAXID in lineage:
            placozoa_taxa.append(taxid)
        elif BILATERIA_TAXID in lineage:
            bilateria_taxa.append(taxid)
        elif METAZOA_TAXID in lineage:
            # Some other metazoan group
            other_taxa.append(taxid)
        else:
            # Non-metazoan
            other_taxa.append(taxid)
    
    print(f"  Ctenophora: {len(ctenophora_taxa)} taxa")
    print(f"  Porifera: {len(porifera_taxa)} taxa")
    print(f"  Cnidaria: {len(cnidaria_taxa)} taxa")
    print(f"  Placozoa: {len(placozoa_taxa)} taxa")
    print(f"  Bilateria: {len(bilateria_taxa)} taxa")
    print(f"  Other: {len(other_taxa)} taxa")
    
    # Build subtrees for each major clade using NCBI topology
    subtrees = {}
    
    if ctenophora_taxa:
        if len(ctenophora_taxa) == 1:
            node = PhyloTree()
            node.name = str(ctenophora_taxa[0])
            subtrees['ctenophora'] = node
        else:
            subtrees['ctenophora'] = ncbi.get_topology(ctenophora_taxa)
    
    if porifera_taxa:
        if len(porifera_taxa) == 1:
            node = PhyloTree()
            node.name = str(porifera_taxa[0])
            subtrees['porifera'] = node
        else:
            subtrees['porifera'] = ncbi.get_topology(porifera_taxa)
    
    if cnidaria_taxa:
        if len(cnidaria_taxa) == 1:
            node = PhyloTree()
            node.name = str(cnidaria_taxa[0])
            subtrees['cnidaria'] = node
        else:
            subtrees['cnidaria'] = ncbi.get_topology(cnidaria_taxa)
    
    if placozoa_taxa:
        if len(placozoa_taxa) == 1:
            node = PhyloTree()
            node.name = str(placozoa_taxa[0])
            subtrees['placozoa'] = node
        else:
            subtrees['placozoa'] = ncbi.get_topology(placozoa_taxa)
    
    if bilateria_taxa:
        if len(bilateria_taxa) == 1:
            node = PhyloTree()
            node.name = str(bilateria_taxa[0])
            subtrees['bilateria'] = node
        else:
            subtrees['bilateria'] = ncbi.get_topology(bilateria_taxa)
    
    # Handle other taxa (non-metazoan or other metazoan groups)
    if other_taxa:
        if len(other_taxa) == 1:
            node = PhyloTree()
            node.name = str(other_taxa[0])
            subtrees['other'] = node
        else:
            subtrees['other'] = ncbi.get_topology(other_taxa)
    
    # Now stitch them together according to custom topology
    # Structure: (Ctenophora,(Porifera,((Cnidaria,Placozoa),Bilateria)Parahoxozoa)Myriazoa)Metazoa
    
    # If we have Cnidaria and/or Placozoa, group them
    cnid_plac_node = None
    if 'cnidaria' in subtrees and 'placozoa' in subtrees:
        cnid_plac_node = PhyloTree()
        cnid_plac_node.name = "CnidariaPlacozoa"
        cnid_plac_node.add_child(subtrees['cnidaria'])
        cnid_plac_node.add_child(subtrees['placozoa'])
    elif 'cnidaria' in subtrees:
        cnid_plac_node = subtrees['cnidaria']
    elif 'placozoa' in subtrees:
        cnid_plac_node = subtrees['placozoa']
    
    # Create Parahoxozoa node (Cnidaria+Placozoa sister to Bilateria)
    parahoxozoa_node = None
    if cnid_plac_node and 'bilateria' in subtrees:
        parahoxozoa_node = PhyloTree()
        parahoxozoa_node.name = str(PARAHOXOZOA_TAXID)
        parahoxozoa_node.add_child(cnid_plac_node)
        parahoxozoa_node.add_child(subtrees['bilateria'])
    elif cnid_plac_node:
        parahoxozoa_node = cnid_plac_node
    elif 'bilateria' in subtrees:
        parahoxozoa_node = subtrees['bilateria']
    
    # Create Myriazoa node (Porifera sister to Parahoxozoa)
    myriazoa_node = None
    if 'porifera' in subtrees and parahoxozoa_node:
        myriazoa_node = PhyloTree()
        myriazoa_node.name = str(MYRIAZOA_TAXID)
        myriazoa_node.add_child(subtrees['porifera'])
        myriazoa_node.add_child(parahoxozoa_node)
    elif 'porifera' in subtrees:
        myriazoa_node = subtrees['porifera']
    elif parahoxozoa_node:
        myriazoa_node = parahoxozoa_node
    
    # Create Metazoa root (Ctenophora sister to Myriazoa)
    if 'ctenophora' in subtrees and myriazoa_node:
        metazoa_node = PhyloTree()
        metazoa_node.name = str(METAZOA_TAXID)
        metazoa_node.add_child(subtrees['ctenophora'])
        metazoa_node.add_child(myriazoa_node)
    elif 'ctenophora' in subtrees:
        metazoa_node = subtrees['ctenophora']
    elif myriazoa_node:
        metazoa_node = myriazoa_node
    else:
        raise ValueError("No metazoan taxa found!")
    
    # If there are non-metazoan taxa, add them at the root
    if 'other' in subtrees:
        root = PhyloTree()
        root.name = "root"
        root.add_child(metazoa_node)
        root.add_child(subtrees['other'])
        return root
    
    return metazoa_node

def main():
    args = parse_args()

    # Read taxids from either taxid_file or config_file
    if args.taxid_file:
        print(f"Reading taxids from file: {args.taxid_file}")
        taxids = set()
        with open(args.taxid_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.isdigit():
                    taxids.add(int(line))  # ensure integers
        print(f"Found {len(taxids)} taxids in file")
    elif args.config_file:
        print(f"Reading taxids from config file: {args.config_file}")
        taxids = read_taxids_from_config(args.config_file)
    else:
        # This should never happen due to mutually_exclusive_group(required=True)
        raise ValueError("Must provide either --taxid_file or --config_file")

    taxids = [int(x) for x in sorted(taxids)]
    print(f"First few taxids: {taxids[:5]}")

    if len(taxids) < 2:
        raise ValueError("You need at least two unique taxids to construct a tree.")

    # Initialize NCBI Taxonomy object
    ncbi = NCBITaxa()

    # Optionally update the local database (slow!)
    # ncbi.update_taxonomy_database()

    # Build tree topology
    if args.custom_phylogeny:
        print("Building tree with custom Ctenophora phylogeny...")
        print("  - Ctenophora (10197) as sister to all other animals")
        print("  - Myriazoa (-67) clade for Porifera + Eumetazoa")
        print("  - Parahoxozoa (-68) clade for Cnidaria+Placozoa sister to Bilateria")
        
        # Build tree by creating subtrees and stitching them together
        tree = build_custom_topology_tree(taxids, ncbi)
    else:
        print("Building tree with NCBI taxonomy...")
        tree = ncbi.get_topology(taxids)

    # Optional: replace tip labels to make FigTree-friendly
    name_dict = ncbi.get_taxid_translator(taxids)
    for leaf in tree:
        if leaf.is_leaf:
            tid = int(leaf.name)
            # Just the species name, with spaces replaced by underscores
            leaf.name = name_dict.get(tid, f"taxid_{tid}").replace(" ", "_")

    # Write to Newick
    # ete4 changed the write() API - no longer uses format/outfile keywords
    tree.write(outfile=args.output_file)
    print(f"Newick tree written to: {args.output_file}")
    
    if args.custom_phylogeny:
        print("\nNote: Custom phylogeny flag enabled.")
        print("      For full custom topology support, consider using the lineage-based")
        print("      approach in plot_ALG_fusions_v3.py before tree construction.")

if __name__ == '__main__':
    main()
