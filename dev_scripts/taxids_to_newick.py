#!/usr/bin/env python

"""
This script takes in a list of NCBI taxids and writes a Newick file.
The point of this program is to be able to open the Newick file in software like FigTree.
"""

import argparse
from ete3 import NCBITaxa

def parse_args():
    parser = argparse.ArgumentParser(description="Create a Newick tree from a list of NCBI TaxIDs")
    parser.add_argument('-t', '--taxid_file', type=str, required=True, help='File containing one NCBI TaxID per line')
    parser.add_argument('-o', '--output_file', type=str, default='ncbi_tree.nwk', help='Output Newick file')
    return parser.parse_args()

def main():
    args = parse_args()

    # Read and clean taxids
    taxids = set()
    with open(args.taxid_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.isdigit():
                taxids.add(int(line))  # ensure integers

    taxids = [int(x) for x in sorted(taxids)]
    print(taxids[:5])

    if len(taxids) < 2:
        raise ValueError("You need at least two unique taxids to construct a tree.")

    # Initialize NCBI Taxonomy object
    ncbi = NCBITaxa()

    # Optionally update the local database (slow!)
    # ncbi.update_taxonomy_database()

    # Build tree topology
    #tree = ncbi.get_topology(taxids, intermediate_nodes=True)
    tree = ncbi.get_topology(taxids)


    # Optional: replace tip labels to make FigTree-friendly
    name_dict = ncbi.get_taxid_translator(taxids)
    for leaf in tree.iter_leaves():
        tid = int(leaf.name)
        # Just the species name, with spaces replaced by underscores
        leaf.name = name_dict.get(tid, f"taxid_{tid}").replace(" ", "_")

    # Write to Newick
    tree.write(format=1, outfile=args.output_file)
    print(f"Newick tree written to: {args.output_file}")

if __name__ == '__main__':
    main()
