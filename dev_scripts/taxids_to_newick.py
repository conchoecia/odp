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
    parser.add_argument('--timetree_list', type=str,
                        help='Output file for TimeTree.org compatible species list (one "Genus species" per line)')
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


def export_timetree_list(taxids, ncbi, output_file):
    """
    Export a list of species names in TimeTree.org compatible format.
    
    TimeTree.org expects one species name per line in "Genus species" format.
    
    Parameters:
    -----------
    taxids : list
        List of taxids to export
    ncbi : NCBITaxa
        NCBI taxonomy database object
    output_file : str
        Path to output file
    """
    print(f"\nExporting TimeTree-compatible species list to: {output_file}")
    
    # Get scientific names for all taxids
    name_dict = ncbi.get_taxid_translator(taxids)
    
    # Filter to only keep proper binomial names (genus + species)
    # TimeTree expects species-level names in "Genus species" format
    species_names = []
    skipped = 0
    
    for taxid in taxids:
        name = name_dict.get(taxid, "")
        # Check if it looks like a binomial (has exactly 2 words)
        parts = name.split()
        if len(parts) >= 2:
            # Take first two parts for binomial
            binomial = f"{parts[0]} {parts[1]}"
            species_names.append(binomial)
        else:
            skipped += 1
            print(f"  Warning: Skipping non-binomial name: {name} (taxid: {taxid})")
    
    # Write to file
    with open(output_file, 'w') as f:
        for name in sorted(species_names):
            f.write(f"{name}\n")
    
    print(f"  Exported {len(species_names)} species names")
    if skipped > 0:
        print(f"  Skipped {skipped} non-species-level taxa")
    print(f"  File ready for upload to TimeTree.org")


def build_subtree_with_labels(taxids, ncbi, root_taxid, root_name):
    """
    Recursively builds a subtree with proper internal node labels.
    
    Parameters:
    -----------
    taxids : list
        List of taxids to include
    ncbi : NCBITaxa
        NCBI taxonomy object
    root_taxid : int
        The taxid of the root of this subtree
    root_name : str
        The name to use for the root node
        
    Returns:
    --------
    PhyloTree : Subtree with labeled internal nodes
    """
    if len(taxids) == 0:
        return None
    
    if len(taxids) == 1:
        # Leaf node
        node = PhyloTree()
        node.name = str(taxids[0])
        return node
    
    # Get the full NCBI tree topology for these taxa
    tree = ncbi.get_topology(taxids)
    
    # Label the root node with the specified root_taxid and root_name
    tree.name = f"{root_name}[{root_taxid}]"
    
    # Now traverse and label all internal nodes based on their descendants
    # Skip the root since we already labeled it
    internal_count = 1  # count the root
    for node in tree.traverse():
        if not node.is_leaf and node != tree:  # Skip root, we already labeled it
            # Get all leaf taxids under this node using traverse
            leaf_taxids = [int(leaf.name) for leaf in node.traverse() if leaf.is_leaf]
            
            if len(leaf_taxids) >= 2:
                # Find common ancestor by getting full lineages (lists, not sets)
                # We need to preserve order to find the most specific common ancestor
                lineages = [ncbi.get_lineage(tid) for tid in leaf_taxids]
                
                # Find common ancestors by intersecting lineages
                common_set = set(lineages[0])
                for lineage in lineages[1:]:
                    common_set = common_set & set(lineage)
                
                if common_set:
                    # Now find the most specific (deepest) common ancestor
                    # by looking for the last common taxid in the first lineage
                    # (lineages go from root to leaf, so last common = most specific)
                    most_specific = None
                    for taxid in reversed(lineages[0]):
                        if taxid in common_set:
                            # Make sure it's not more ancestral than our root
                            if root_taxid > 0:
                                root_lineage_set = set(ncbi.get_lineage(root_taxid))
                                # Only use this if it's within our clade (root_taxid is in its lineage or it equals root_taxid)
                                if taxid == root_taxid or root_taxid in ncbi.get_lineage(taxid):
                                    most_specific = taxid
                                    break
                            else:
                                most_specific = taxid
                                break
                    
                    if most_specific:
                        try:
                            taxon_name = ncbi.get_taxid_translator([most_specific]).get(most_specific, f"taxid_{most_specific}")
                            node.name = f"{taxon_name.replace(' ', '_')}[{most_specific}]"
                            internal_count += 1
                        except Exception as e:
                            node.name = f"node[{most_specific}]"
                            internal_count += 1
            elif len(leaf_taxids) == 1:
                # Single leaf under this node - shouldn't happen but handle it
                node.name = f"internal[{leaf_taxids[0]}]"
                internal_count += 1
    
    print(f"    Labeled {internal_count} internal nodes in {root_name} subtree")
    return tree


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
    ete4.Tree : Tree with custom topology and labeled nodes
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
    
    # Build subtrees for each major clade with labeled internal nodes
    subtrees = {}
    
    if ctenophora_taxa:
        subtrees['ctenophora'] = build_subtree_with_labels(ctenophora_taxa, ncbi, CTENOPHORA_TAXID, "Ctenophora")
    
    if porifera_taxa:
        subtrees['porifera'] = build_subtree_with_labels(porifera_taxa, ncbi, PORIFERA_TAXID, "Porifera")
    
    if cnidaria_taxa:
        subtrees['cnidaria'] = build_subtree_with_labels(cnidaria_taxa, ncbi, CNIDARIA_TAXID, "Cnidaria")
    
    if placozoa_taxa:
        subtrees['placozoa'] = build_subtree_with_labels(placozoa_taxa, ncbi, PLACOZOA_TAXID, "Placozoa")
    
    if bilateria_taxa:
        subtrees['bilateria'] = build_subtree_with_labels(bilateria_taxa, ncbi, BILATERIA_TAXID, "Bilateria")
    
    if other_taxa:
        subtrees['other'] = build_subtree_with_labels(other_taxa, ncbi, -1, "Other")
    
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
        parahoxozoa_node.name = f"Parahoxozoa[{PARAHOXOZOA_TAXID}]"
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
        myriazoa_node.name = f"Myriazoa[{MYRIAZOA_TAXID}]"
        myriazoa_node.add_child(subtrees['porifera'])
        myriazoa_node.add_child(parahoxozoa_node)
    elif 'porifera' in subtrees:
        myriazoa_node = subtrees['porifera']
    elif parahoxozoa_node:
        myriazoa_node = parahoxozoa_node
    
    # Create Metazoa root (Ctenophora sister to Myriazoa)
    if 'ctenophora' in subtrees and myriazoa_node:
        metazoa_node = PhyloTree()
        metazoa_node.name = f"Metazoa[{METAZOA_TAXID}]"
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
    
    # Set up logging to both console and file
    import sys
    log_file = args.output_file.replace('.nwk', '.log').replace('.newick', '.log')
    if not log_file.endswith('.log'):
        log_file = args.output_file + '.log'
    
    class Tee:
        """Write to both stdout and a file"""
        def __init__(self, *files):
            self.files = files
        def write(self, obj):
            for f in self.files:
                f.write(obj)
                f.flush()
        def flush(self):
            for f in self.files:
                f.flush()
    
    log_handle = open(log_file, 'w')
    original_stdout = sys.stdout
    sys.stdout = Tee(sys.stdout, log_handle)
    
    print(f"Logging output to: {log_file}")
    print(f"Output tree file: {args.output_file}")
    print()

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

    # Label leaf nodes with species names (internal nodes already labeled during tree construction)
    print("\nLabeling leaf nodes with species names...")
    name_dict = ncbi.get_taxid_translator(taxids)
    
    leaf_count = 0
    for node in tree.traverse():
        if node.is_leaf:
            tid = int(node.name)
            node.name = name_dict.get(tid, f"taxid_{tid}").replace(" ", "_")
            leaf_count += 1
    
    print(f"  Labeled {leaf_count} leaf nodes")
    print(f"  Total nodes in tree: {len(list(tree.traverse()))}")
    print(f"  - {len([n for n in tree.traverse() if n.is_leaf])} leaf nodes (species)")
    print(f"  - {len([n for n in tree.traverse() if not n.is_leaf])} internal nodes")
    
    # Count internal nodes by number of children
    print("\nInternal node children counts:")
    children_counts = {}
    for node in tree.traverse():
        if not node.is_leaf:
            num_children = len(node.children)
            children_counts[num_children] = children_counts.get(num_children, 0) + 1
    
    for num_children in sorted(children_counts.keys()):
        count = children_counts[num_children]
        print(f"  Internal nodes with {num_children} children: {count}")
    
    # Collapse internal nodes with only one child
    if 1 in children_counts:
        print(f"\nCollapsing {children_counts[1]} internal nodes with single children...")
        nodes_to_check = [tree]
        collapsed_count = 0
        
        while nodes_to_check:
            node = nodes_to_check.pop(0)
            
            # Check each child
            for child in list(node.children):  # Use list() to avoid modifying during iteration
                if not child.is_leaf and len(child.children) == 1:
                    # This child has only one child - collapse it
                    grandchild = child.children[0]
                    node.remove_child(child)
                    node.add_child(grandchild)
                    collapsed_count += 1
                    # Check this node again in case there are multiple levels
                    nodes_to_check.append(node)
                else:
                    # Add child to check queue if it's not a leaf
                    if not child.is_leaf:
                        nodes_to_check.append(child)
        
        print(f"  Collapsed {collapsed_count} single-child internal nodes")
        
        # Recount after collapsing
        print("\nInternal node children counts after collapsing:")
        children_counts_after = {}
        for node in tree.traverse():
            if not node.is_leaf:
                num_children = len(node.children)
                children_counts_after[num_children] = children_counts_after.get(num_children, 0) + 1
        
        for num_children in sorted(children_counts_after.keys()):
            count = children_counts_after[num_children]
            print(f"  Internal nodes with {num_children} children: {count}")
    
    # Debug: check some internal node names
    print("\nDebug: Checking internal node names...")
    internal_with_names = 0
    internal_without_names = 0
    for node in tree.traverse():
        if not node.is_leaf:
            if node.name and node.name.strip():
                internal_with_names += 1
                if internal_with_names <= 5:  # Show first 5
                    print(f"  Sample internal node: '{node.name}'")
            else:
                internal_without_names += 1
    print(f"  Internal nodes WITH names: {internal_with_names}")
    print(f"  Internal nodes WITHOUT names: {internal_without_names}")

    # Write to Newick with internal node names
    # In ete4, we need to manually create the newick string to include internal node names
    # The write() method in ete4 doesn't have the same format parameter as ete3
    
    def needs_quotes(name):
        """Check if a name needs quotes in Newick format."""
        special_chars = set('():;,[] ')
        return any(c in special_chars for c in name)
    
    def to_newick_with_internal_names(node):
        """Recursively build newick string including internal node names."""
        if node.is_leaf:
            # Quote leaf names if they contain special characters
            if needs_quotes(node.name):
                return f"'{node.name}'"
            return node.name
        else:
            # Get newick strings for all children
            children_newick = ','.join([to_newick_with_internal_names(child) for child in node.children])
            # Include the node name if it exists
            if node.name:
                # Quote internal node names if they contain special characters
                if needs_quotes(node.name):
                    return f"({children_newick})'{node.name}'"
                return f"({children_newick}){node.name}"
            else:
                return f"({children_newick})"
    
    newick_string = to_newick_with_internal_names(tree) + ";"
    
    # Validate parentheses are balanced
    open_count = newick_string.count('(')
    close_count = newick_string.count(')')
    if open_count != close_count:
        print(f"WARNING: Unbalanced parentheses in newick string!")
        print(f"  Open: {open_count}, Close: {close_count}")
    else:
        print(f"  Newick validation: {open_count} balanced parentheses")
    
    with open(args.output_file, 'w') as f:
        f.write(newick_string)
    
    print(f"\nNewick tree written to: {args.output_file}")
    print(f"  (Includes internal node labels)")
    
    # Export TimeTree-compatible species list if requested
    if args.timetree_list:
        export_timetree_list(taxids, ncbi, args.timetree_list)
    
    if args.custom_phylogeny:
        print("\nNote: Custom phylogeny flag enabled.")
        print("      For full custom topology support, consider using the lineage-based")
        print("      approach in plot_ALG_fusions_v3.py before tree construction.")
    
    # Close log file and restore stdout
    sys.stdout = original_stdout
    log_handle.close()
    print(f"Log saved to: {log_file}")

if __name__ == '__main__':
    main()
