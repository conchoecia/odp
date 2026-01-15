#!/usr/bin/env python3
"""
Create a collapsed phylogenetic tree visualization showing major animal clades
with branch statistics (fusions and dispersals) labeled on branches.

This generates a figure similar to simplified phylogenies
with major clades collapsed and event counts shown on branches.

Usage:
    python plot_collapsed_tree.py -n node_stats.tsv -e edge_stats.tsv -o output.pdf
"""

import argparse
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Polygon

# Import tree handling
from Newick_to_common_ancestors import TaxIDtree
from odp_plotting_functions import format_matplotlib


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Create collapsed phylogenetic tree with branch statistics')
    parser.add_argument('-n', '--node_stats', type=str, required=True, help='Node stats TSV file')
    parser.add_argument('-e', '--edge_stats', type=str, required=True, help='Edge stats TSV file')
    parser.add_argument('-o', '--output', type=str, default='collapsed_tree.pdf', help='Output PDF filename')
    parser.add_argument('--clades', type=str, help='Optional: comma-separated list of taxids to display')
    return parser.parse_args()


# Define major clades to display
MAJOR_CLADES = {
    10197: {'name': 'Ctenophora', 'color': '#E8B4D4'},
    6040: {'name': 'Porifera', 'color': '#D4E8B4'},
    6073: {'name': 'Cnidaria', 'color': '#B4D4E8'},
    6340: {'name': 'Annelida', 'color': '#FFD4B4'},
    6447: {'name': 'Mollusca', 'color': '#D4B4FF'},
    6960: {'name': 'Hexapoda', 'color': '#FFB4D4'},
    6657: {'name': 'Crustacea', 'color': '#B4FFD4'},
    6843: {'name': 'Myriapoda', 'color': '#FFE8B4'},
    50557: {'name': 'Insecta', 'color': '#D4FFB4'},
    7742: {'name': 'Vertebrata', 'color': '#B4E8FF'},
    7586: {'name': 'Echinodermata', 'color': '#FFB4E8'},
    33213: {'name': 'Bilateria', 'color': '#E8FFB4'},
    33317: {'name': 'Protostomia', 'color': '#B4FFE8'},
    33511: {'name': 'Deuterostomia', 'color': '#FFD4E8'},
}


def get_clade_mrca(tree, clade_taxid):
    """
    Find the MRCA node for a given clade taxid.
    
    Args:
        tree: TaxIDtree instance
        clade_taxid: Taxid of the clade to find
        
    Returns:
        TaxNode representing the MRCA, or None if not found
    """
    # Check if this taxid is directly in the tree
    if clade_taxid in tree.nodes:
        return tree.nodes[clade_taxid]
    
    # Otherwise, find the first node whose lineage contains this taxid
    for node_id, node in tree.nodes.items():
        if hasattr(node, 'lineage') and clade_taxid in node.lineage:
            return node
    
    return None


def aggregate_clade_statistics(tree, clade_taxid):
    """
    Aggregate all fusions and dispersals within a clade.
    
    Args:
        tree: TaxIDtree instance
        clade_taxid: Taxid of the clade to aggregate
        
    Returns:
        Dictionary with 'fusions', 'dispersals', 'num_branches', 'num_tips'
    """
    total_fusions = 0
    total_dispersals = 0
    num_branches = 0
    tips_in_clade = set()
    
    # Find all edges within this clade
    for edge_key, edge in tree.edges.items():
        # Check if this edge is within the clade using child node's lineage
        child_node = tree.nodes.get(edge.child_taxid)
        if child_node and hasattr(child_node, 'lineage') and clade_taxid in child_node.lineage:
            num_branches += 1
            
            # Sum statistics
            if hasattr(edge, 'num_fusions_this_branch'):
                val = edge.num_fusions_this_branch
                if isinstance(val, (int, float)) and not np.isnan(val):
                    total_fusions += val
            if hasattr(edge, 'num_losses_this_branch'):
                val = edge.num_losses_this_branch
                if isinstance(val, (int, float)) and not np.isnan(val):
                    total_dispersals += val
            
            # Track tips
            if edge.child_taxid in tree.nodes:
                child_node = tree.nodes[edge.child_taxid]
                if hasattr(child_node, 'isleaf') and child_node.isleaf:
                    tips_in_clade.add(edge.child_taxid)
    
    return {
        'fusions': int(total_fusions),
        'dispersals': int(total_dispersals),
        'num_branches': num_branches,
        'num_tips': len(tips_in_clade)
    }


def get_stem_branch_statistics(tree, clade_taxid):
    """
    Get statistics for the branch leading TO a clade (stem branch).
    
    Args:
        tree: TaxIDtree instance
        clade_taxid: Taxid of the clade
        
    Returns:
        Dictionary with 'fusions', 'dispersals' for the stem branch, or None
    """
    # Find the MRCA node
    mrca_node = get_clade_mrca(tree, clade_taxid)
    if mrca_node is None:
        return None
    
    # Find the edge leading to this node
    for edge_key, edge in tree.edges.items():
        if edge.child_taxid == mrca_node.taxid:
            fusions = edge.num_fusions_this_branch if hasattr(edge, 'num_fusions_this_branch') else 0
            dispersals = edge.num_losses_this_branch if hasattr(edge, 'num_losses_this_branch') else 0
            
            # Safely convert to int, handling NaN and non-numeric values
            if isinstance(fusions, (int, float)) and not np.isnan(fusions):
                fusions = int(fusions)
            else:
                fusions = 0
            
            if isinstance(dispersals, (int, float)) and not np.isnan(dispersals):
                dispersals = int(dispersals)
            else:
                dispersals = 0
            
            return {
                'fusions': fusions,
                'dispersals': dispersals,
                'parent_taxid': edge.parent_taxid,
                'parent_age': edge.parent_age if hasattr(edge, 'parent_age') else None,
                'child_age': edge.child_age if hasattr(edge, 'child_age') else None
            }
    
    return None


def plot_collapsed_tree_custom(tree, clade_definitions, ax, show_within_clade_stats=True):
    """
    Custom plotting function for collapsed phylogenetic tree.
    
    Args:
        tree: TaxIDtree instance
        clade_definitions: Dictionary mapping taxid to display info
        ax: Matplotlib axis
        show_within_clade_stats: If True, show aggregated statistics within collapsed clades
    """
    format_matplotlib()
    
    # First, plot the base tree structure to get positioning
    ax = tree.plot_tree(ax, sort="ascending", draw_horizontal_bars=True, text_older_than=1e12)
    
    # Now overlay branch statistics as text labels
    y_positions = {}  # Track y positions of nodes for label placement
    
    # Get node positions from tree
    for node_id, node in tree.nodes.items():
        if hasattr(node, 'y'):
            y_positions[node_id] = node.y
    
    # Add labels for each clade
    for clade_taxid, clade_info in clade_definitions.items():
        # Get stem branch statistics
        stem_stats = get_stem_branch_statistics(tree, clade_taxid)
        
        # Get within-clade statistics
        if show_within_clade_stats:
            clade_stats = aggregate_clade_statistics(tree, clade_taxid)
        else:
            clade_stats = None
        
        # Find node position
        mrca_node = get_clade_mrca(tree, clade_taxid)
        if mrca_node is None:
            continue
        
        # Add stem branch label if available
        if stem_stats and stem_stats['parent_age'] and stem_stats['child_age']:
            # Position label at midpoint of branch
            x_pos = (stem_stats['parent_age'] + stem_stats['child_age']) / 2
            y_pos = y_positions.get(mrca_node.taxid, 0)
            
            label = f"{stem_stats['fusions']} / {stem_stats['dispersals']}"
            ax.text(x_pos, y_pos, label, fontsize=8, ha='center', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.7))
        
        # Add within-clade label if available
        if clade_stats and show_within_clade_stats:
            # Position near the clade tip
            x_pos = 0 if hasattr(mrca_node, 'nodeage') else -50
            y_pos = y_positions.get(mrca_node.taxid, 0)
            
            clade_label = f"{clade_info['name']}\n({clade_stats['fusions']} / {clade_stats['dispersals']})"
            ax.text(x_pos + 5, y_pos, clade_label, fontsize=9, ha='left', va='center',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor=clade_info.get('color', 'lightgray'), 
                            edgecolor='black', alpha=0.8))
    
    ax.set_xlabel('Age (Ma)', fontsize=12)
    ax.set_ylabel('Taxa', fontsize=12)
    ax.set_title('Collapsed Phylogenetic Tree: Fusions / Dispersals', fontsize=14, fontweight='bold')
    
    # Add legend
    legend_text = 'Branch labels: Fusions / Dispersals (ALG losses)'
    ax.text(0.02, 0.98, legend_text, transform=ax.transAxes, fontsize=9,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    return ax


def main():
    """Main execution function."""
    args = parse_args()
    
    # Load data
    print(f"Loading edge stats from: {args.edge_stats}")
    edgedf = pd.read_csv(args.edge_stats, sep='\t')
    
    print(f"Loading node stats from: {args.node_stats}")
    nodedf = pd.read_csv(args.node_stats, sep='\t')
    
    # Drop per-clade summary columns (from plot_branch_stats_vs_time.py Phase 1)
    # These aren't needed for tree visualization and TaxNode doesn't recognize them
    clade_cols = [col for col in nodedf.columns if 'in_this_clade' in col]
    if clade_cols:
        nodedf = nodedf.drop(columns=clade_cols)
        print(f"Dropped {len(clade_cols)} per-clade summary columns from nodedf")
    
    # Drop columns from edgedf that TaxEdge doesn't recognize
    # Keep only the expected TaxEdge attributes
    expected_edge_cols = ['parent_taxid', 'child_taxid', 'parent_age', 'child_age', 
                          'branch_length', 'dist_crown_plus_this_edge',
                          'parent_lineage', 'child_lineage',
                          'num_fusions_this_branch', 'num_losses_this_branch',
                          'num_fusions_per_my_this_branch', 'num_losses_per_my_this_branch',
                          'num_dispersals_per_my_this_branch',
                          'fusions', 'losses']
    extra_edge_cols = [col for col in edgedf.columns if col not in expected_edge_cols]
    if extra_edge_cols:
        edgedf = edgedf.drop(columns=extra_edge_cols)
        print(f"Dropped {len(extra_edge_cols)} extra columns from edgedf")
    
    # Build tree (don't add helper columns like child_lineage_list yet)
    print("Building phylogenetic tree...")
    tree = TaxIDtree()
    tree.ingest_node_edge(nodedf, edgedf)
    
    print(f"Tree built: {len(tree.nodes)} nodes, {len(tree.edges)} edges")
    
    # Determine which clades to show
    if args.clades:
        clade_taxids = [int(x.strip()) for x in args.clades.split(',')]
        clade_defs = {tid: MAJOR_CLADES.get(tid, {'name': f'Clade_{tid}', 'color': 'lightgray'}) 
                     for tid in clade_taxids}
    else:
        clade_defs = MAJOR_CLADES
    
    # Aggregate statistics
    print(f"\nAggregating statistics for {len(clade_defs)} major clades...")
    for clade_taxid, clade_info in clade_defs.items():
        stats = aggregate_clade_statistics(tree, clade_taxid)
        stem = get_stem_branch_statistics(tree, clade_taxid)
        
        print(f"  {clade_info['name']} (taxid={clade_taxid}):")
        print(f"    Within clade: {stats['fusions']} fusions, {stats['dispersals']} dispersals "
              f"({stats['num_tips']} tips, {stats['num_branches']} branches)")
        if stem:
            print(f"    Stem branch: {stem['fusions']} fusions, {stem['dispersals']} dispersals")
    
    # Create figure
    print(f"\nGenerating collapsed tree plot...")
    fig, ax = plt.subplots(1, 1, figsize=(16, 12))
    
    ax = plot_collapsed_tree_custom(tree, clade_defs, ax, show_within_clade_stats=True)
    
    # Save
    print(f"Saving to: {args.output}")
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    
    print("Done!")


if __name__ == '__main__':
    main()
