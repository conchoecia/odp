#!/usr/bin/env python

from __future__ import annotations

"""
This script makes plots related to plotting the changes on phylogenetic branches over time.

Takes in a file called "Jun240604.edge_information.tsv" that has the following columns:
    - source: the parent node
    - target: the child node
    - source_age: the age of the parent node
    - target_age: the age of the child node
    - branch_length: the branch length between the parent and child node.
    - source_ages: the age counter of the parent node
    - target_ages: the age counter of the child node

Also takes in the df_stats file that has fusions or dispersals.
"""

import argparse
import sys

# Force unbuffered output for real-time monitoring in SLURM logs
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

import os
import pickle
import hashlib
import copy
from multiprocessing import Pool, cpu_count
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
from scipy.stats import pearsonr, spearmanr, kendalltau
import warnings
import logging

from odp_plotting_functions import format_matplotlib

# Configure matplotlib for Illustrator-editable PDFs
plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts for editable text in Illustrator
plt.rcParams['ps.fonttype'] = 42   # Also for EPS output
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['savefig.edgecolor'] = 'none'

def weighted_quantile(values, weights, quantile):
    """
    Compute weighted quantile accounting for phylogenetic weights.
    
    Args:
        values: Array of values (e.g., rates)
        weights: Array of weights (e.g., phylogenetic weights)
        quantile: Desired quantile (0-1)
    
    Returns:
        The value at the specified weighted quantile
    """
    if len(values) == 0:
        return np.nan
    
    values = np.array(values)
    weights = np.array(weights)
    
    # Sort by values
    indices = np.argsort(values)
    sorted_values = values[indices]
    sorted_weights = weights[indices]
    
    # Cumulative weights
    cum_weights = np.cumsum(sorted_weights)
    cum_weights = cum_weights / cum_weights[-1]  # Normalize to [0, 1]
    
    # Interpolate to find value at desired quantile
    return np.interp(quantile, cum_weights, sorted_values)

def renormalize_weights_for_clade(global_phylo_weighter, clade_tip_taxids, time_slice):
    """
    Load global weights for a time slice and renormalize for a specific clade.
    
    Args:
        global_phylo_weighter: PhyloWeighting object from full tree analysis
        clade_tip_taxids: Set of taxids that are tips in this clade
        time_slice: Time slice to load
    
    Returns:
        Dict {taxid: renormalized_weight} for tips in this clade
    """
    # Load global weights from cache
    global_weights = global_phylo_weighter._load_cached_weights(time_slice=time_slice)
    
    if global_weights is None:
        return None
    
    # Filter to only clade tips
    clade_weights = {taxid: weight for taxid, weight in global_weights.items() 
                     if taxid in clade_tip_taxids}
    
    if not clade_weights:
        return {}
    
    # Renormalize so weights average to 1.0 within this clade
    mean_weight = sum(clade_weights.values()) / len(clade_weights)
    renormalized_weights = {taxid: weight / mean_weight 
                           for taxid, weight in clade_weights.items()}
    
    return renormalized_weights

def calculate_clade_event_totals(taxid, edgedf):
    """
    Fast calculation of UNWEIGHTED event totals for a clade.
    Does NOT generate detailed temporal breakdown or plots.
    Weighted calculations happen later in Phase 2 per-clade analyses.
    
    Args:
        taxid: The taxid to filter branches by
        edgedf: DataFrame with branch information
    
    Returns:
        Tuple of (unweighted_fusions, unweighted_dispersals)
    """
    # Filter edgedf to this clade
    edgedf_clade = edgedf[edgedf['child_lineage_list'].apply(lambda x: taxid in x)].copy()
    
    if len(edgedf_clade) == 0:
        return (0.0, 0.0)
    
    # Calculate unweighted totals
    unweighted_fusions = 0.0
    unweighted_dispersals = 0.0
    
    for _, row in edgedf_clade.iterrows():
        num_fusions = row['num_fusions_per_my_this_branch']
        num_dispersals = row['num_dispersals_per_my_this_branch']
        branch_length = row['parent_age'] - row['child_age']
        
        if np.isnan(num_fusions) or np.isinf(num_fusions):
            num_fusions = 0
        if np.isnan(num_dispersals) or np.isinf(num_dispersals):
            num_dispersals = 0
        
        unweighted_fusions += branch_length * num_fusions
        unweighted_dispersals += branch_length * num_dispersals
    
    return (unweighted_fusions, unweighted_dispersals)

class PhyloWeighting:
    """
    Simplified phylogenetic weighting for time trees.
    Computes weights based on phylogenetic distinctiveness to account for oversampled clades.
    """
    def __init__(self, edgedf, nodedf, cache_dir=None, verbose=True):
        """Build DAG from edge structure using time-calibrated ages."""
        self.dag = {}  # node -> {child: branch_length}
        self.tips = set()
        self.cache_dir = cache_dir
        self.edgedf = edgedf  # Store for time-slice queries
        self.nodedf = nodedf  # Store for getting node ages
        self.verbose = verbose  # Control output messages
        
        # Build graph from edges
        for _, row in edgedf.iterrows():
            parent = row['parent_taxid']
            child = row['child_taxid']
            branch_length = row['branch_length']
            
            if parent not in self.dag:
                self.dag[parent] = {}
            self.dag[parent][child] = branch_length
        
        # Identify tips (nodes with no children)
        all_nodes = set(self.dag.keys())
        all_children = set()
        for children in self.dag.values():
            all_children.update(children.keys())
        self.tips = all_children - all_nodes
        
        # Add empty dict for tips
        for tip in self.tips:
            if tip not in self.dag:
                self.dag[tip] = {}
        
        # Find root
        self.root = list(all_nodes - all_children)[0]
        if self.verbose:
            print(f"Phylogenetic weighting: root={self.root}, {len(self.tips)} tips, {len(self.dag)} nodes")
        
        # Compute tree hash for cache invalidation
        self.tree_hash = self._compute_tree_hash(edgedf)
        if self.verbose:
            print(f"Tree hash: {self.tree_hash}")
    
    def _compute_tree_hash(self, edgedf):
        """Compute hash of tree structure to detect changes."""
        hash_data = edgedf[['parent_taxid', 'child_taxid', 'branch_length']].copy()
        hash_data = hash_data.sort_values(['parent_taxid', 'child_taxid']).reset_index(drop=True)
        hash_string = hash_data.to_csv(index=False)
        return hashlib.md5(hash_string.encode()).hexdigest()[:16]
    
    def _get_cache_paths(self, time_slice=None):
        """Get cache file paths."""
        if self.cache_dir is None:
            return None, None
        cache_subdir = os.path.join(self.cache_dir, 'phylo_cache')
        os.makedirs(cache_subdir, exist_ok=True)
        
        if time_slice is not None:
            # Time-specific cache
            time_str = f"t{int(time_slice)}"
            distances_cache = os.path.join(cache_subdir, f'tip_distances_{self.tree_hash}_{time_str}.pkl')
            weights_cache = os.path.join(cache_subdir, f'edge_weights_{self.tree_hash}_{time_str}.pkl')
        else:
            # Present-day cache
            distances_cache = os.path.join(cache_subdir, f'tip_distances_{self.tree_hash}.pkl')
            weights_cache = os.path.join(cache_subdir, f'edge_weights_{self.tree_hash}.pkl')
        return distances_cache, weights_cache
    
    def _load_cached_distances(self, time_slice=None):
        """Load cached tip distances if available."""
        distances_cache, _ = self._get_cache_paths(time_slice)
        if distances_cache and os.path.exists(distances_cache):
            try:
                with open(distances_cache, 'rb') as f:
                    return pickle.load(f)
            except Exception as e:
                print(f"Warning: Failed to load cached distances: {e}")
                return None
        return None
    
    def _save_cached_distances(self, tip_distances, time_slice=None):
        """Save tip distances to cache."""
        distances_cache, _ = self._get_cache_paths(time_slice)
        if distances_cache:
            try:
                with open(distances_cache, 'wb') as f:
                    pickle.dump(tip_distances, f)
            except Exception as e:
                print(f"Warning: Failed to save cached distances: {e}")
    
    def _load_cached_weights(self, time_slice=None):
        """Load cached weights if available."""
        _, weights_cache = self._get_cache_paths(time_slice)
        if weights_cache and os.path.exists(weights_cache):
            try:
                with open(weights_cache, 'rb') as f:
                    return pickle.load(f)
            except Exception as e:
                print(f"Warning: Failed to load cached weights: {e}")
                return None
        return None
    
    def _save_cached_weights(self, weights, time_slice=None):
        """Save weights to cache."""
        _, weights_cache = self._get_cache_paths(time_slice)
        if weights_cache:
            try:
                with open(weights_cache, 'wb') as f:
                    pickle.dump(weights, f)
            except Exception as e:
                print(f"Warning: Failed to save cached weights: {e}")
    
    def _path_to_root(self, node, parent_map):
        """Return path from node to root."""
        path = []
        current = node
        while current in parent_map:
            parent = parent_map[current]
            path.append((current, parent))
            current = parent
        return path
    
    def compute_tip_distances(self, pseudo_tips=None, time_slice=None):
        """Compute pairwise distances between all tips (or pseudo-tips at a time slice)."""
        # Try to load from cache first
        cached_distances = self._load_cached_distances(time_slice)
        if cached_distances is not None:
            return cached_distances
        
        # Use provided pseudo_tips or default to actual tips
        tips_to_use = pseudo_tips if pseudo_tips is not None else self.tips
        
        # Only show detailed progress for non-time-slice computations
        if time_slice is None:
            print(f"Computing tip distances from scratch (this may take a few minutes)...")
        
        # Build parent map
        parent_map = {}
        for parent, children in self.dag.items():
            for child, dist in children.items():
                parent_map[child] = (parent, dist)
        
        # Precompute all paths to root (major optimization)
        tip_list = list(tips_to_use)
        total_tips = len(tip_list)
        paths_to_root = {}
        
        # Suppress per-worker progress when running in parallel (time_slice mode)
        # The main process will show overall progress instead
        
        for tip in tip_list:
            path = {}
            current = tip
            total_dist = 0
            while current in parent_map:
                parent, dist = parent_map[current]
                path[parent] = total_dist
                total_dist += dist
                current = parent
            path[current] = total_dist  # Add root
            paths_to_root[tip] = path
        
        # Compute pairwise distances using precomputed paths
        distances = {}
        
        for i, tip1 in enumerate(tip_list):
            # Only show progress for non-parallel mode (time_slice is None)
            if time_slice is None:
                print(f"\rComputing tip distances: {i+1}/{total_tips} ({100*(i+1)/total_tips:.1f}%)", end='', flush=True)
            
            path1 = paths_to_root[tip1]
            
            for tip2 in tip_list[i:]:
                if tip1 == tip2:
                    distances[(tip1, tip2)] = 0
                    continue
                
                path2 = paths_to_root[tip2]
                
                # Find MRCA - the closest common ancestor
                mrca = None
                for node in path2:
                    if node in path1:
                        mrca = node
                        break
                
                if mrca is not None:
                    distance = path1[mrca] + path2[mrca]
                    distances[(tip1, tip2)] = distance
                    distances[(tip2, tip1)] = distance
                else:
                    distances[(tip1, tip2)] = 0
                    distances[(tip2, tip1)] = 0
        
        print()  # New line after progress complete
        
        # Save to cache
        self._save_cached_distances(distances, time_slice)
        
        return distances
    
    def compute_branch_weights(self, edgedf, time_slice=None, pruned_dag=None):
        """
        Compute phylogenetic weight for each branch.
        Weight based on sum of phylogenetic distances to all tips (or pseudo-tips at a time slice).
        
        Args:
            edgedf: DataFrame with branch information
            time_slice: If provided, compute weights using only branches crossing this time point
            pruned_dag: If provided, use this DAG for descendant lookups (for time-slicing)
        """
        # Determine which branches to use as "tips" for weighting
        if time_slice is not None:
            pseudo_tips = self._get_branches_at_time(time_slice)
            if not pseudo_tips:
                # Silently return all 1.0 weights for empty time slices (no output to avoid spam)
                return [1.0] * len(edgedf)
            # Removed verbose print - will show progress at higher level
        else:
            pseudo_tips = None
            print("Computing weights using present-day tips...")
        
        # Try to load from cache first
        cached_data = self._load_cached_weights(time_slice)
        if cached_data is not None:
            # Handle both dict and array formats from cache
            if isinstance(cached_data, dict):
                # For time-slicing, return dict directly
                if time_slice is not None:
                    # Print statistics for cached time-slice weights (only if verbose)
                    if self.verbose:
                        weights_values = list(cached_data.values())
                        if weights_values:
                            print(f"Phylogenetic weights: min={min(weights_values):.3f}, mean=1.000, max={max(weights_values):.3f}")
                    return cached_data
                # For non-time-slicing, convert dict back to array matching edgedf order
                cached_weights = [cached_data.get(child_id, 1.0) for child_id in edgedf['child_taxid']]
                if len(cached_weights) == len(edgedf):
                    return cached_weights
            elif len(cached_data) == len(edgedf):
                return cached_data
            # Silently recompute if cache size mismatch
        
        # Only print this for non-time-slice or first time slice (to avoid spam)
        if time_slice is None or time_slice == 0:
            print("Computing phylogenetic distances between tips...")
        tip_distances = self.compute_tip_distances(pseudo_tips, time_slice)
        
        # For time-slice weighting: compute weights directly from tip distances
        # No need for descendant traversal since tips have no descendants
        if time_slice is not None:
            tips_to_use = pseudo_tips
            tip_list = list(tips_to_use)
            num_tips = len(tip_list)
            
            # Verbose output suppressed - main process shows overall progress
            
            # Vectorized weight computation - build distance matrix
            tip_to_idx = {tip: idx for idx, tip in enumerate(tip_list)}
            dist_matrix = np.zeros((num_tips, num_tips))
            
            for (tip1, tip2), dist in tip_distances.items():
                if tip1 in tip_to_idx and tip2 in tip_to_idx:
                    i, j = tip_to_idx[tip1], tip_to_idx[tip2]
                    dist_matrix[i, j] = dist
                    dist_matrix[j, i] = dist
            
            # Sum distances for each tip (vectorized)
            distance_sums = dist_matrix.sum(axis=1)
            
            # Create weights dict
            weights_dict = {}
            for idx, tip in enumerate(tip_list):
                avg_dist = distance_sums[idx] / max(num_tips - 1, 1)
                weights_dict[tip] = avg_dist
            
            # Normalize so mean = 1
            weights_array = np.array(list(weights_dict.values()))
            mean_weight = np.mean(weights_array)
            weights_dict = {tip: weight / mean_weight for tip, weight in weights_dict.items()}
            
            if self.verbose:
                print(f"Phylogenetic weights: min={min(weights_dict.values()):.3f}, mean=1.000, max={max(weights_dict.values()):.3f}")
            
            # Save as dictionary for diagnostics
            self._save_cached_weights(weights_dict, time_slice)
            
            # Return dictionary (not array) for time-sliced weights
            return weights_dict
        
        # Original code for non-time-sliced weighting (full tree)
        weights = []
        tips_to_use = self.tips
        tip_list = list(tips_to_use)
        total_edges = len(edgedf)
        
        for counter, (idx, row) in enumerate(edgedf.iterrows(), 1):
            print(f"\rComputing branch weights: {counter}/{total_edges} ({100*counter/total_edges:.1f}%)", end='', flush=True)
            child = row['child_taxid']
            
            # Find all descendant tips
            descendants = self._get_descendants(child)
            
            if not descendants:
                # Leaf branch - use distances to all other tips
                if child in tip_list:
                    total_dist = sum(tip_distances.get((child, other), 0) 
                                    for other in tip_list if other != child)
                    avg_dist = total_dist / max(len(tip_list) - 1, 1)
                else:
                    avg_dist = 1.0
            else:
                # Internal branch - average distinctiveness of descendant tips
                dist_sums = []
                for desc in descendants:
                    if desc in tip_list:
                        total_dist = sum(tip_distances.get((desc, other), 0)
                                       for other in tip_list if other != desc)
                        dist_sums.append(total_dist / max(len(tip_list) - 1, 1))
                avg_dist = np.mean(dist_sums) if dist_sums else 1.0
            
            weights.append(avg_dist)
        
        # Normalize so mean = 1 (preserves total event counts)
        weights = np.array(weights)
        weights = weights / np.mean(weights)
        
        print()  # New line after progress
        print(f"Phylogenetic weights: min={weights.min():.3f}, mean={weights.mean():.3f}, max={weights.max():.3f}")
        
        # Save to cache
        self._save_cached_weights(weights, time_slice)
        
        return weights
    
    def _get_descendants(self, node, dag=None, tips=None):
        """Get all descendant tips of a node."""
        if dag is None:
            dag = self.dag
            tips = self.tips
        elif tips is None:
            # For pruned DAG, compute tips once if not provided
            tips = {n for n in dag if not dag[n]}
        
        if node in tips:
            return {node}
        
        descendants = set()
        if node in dag:
            for child in dag[node]:
                descendants.update(self._get_descendants(child, dag, tips))
        return descendants
    
    def _get_branches_at_time(self, time_mya):
        """Find all branches that cross a given time point."""
        # Ages are stored in edgedf as parent_age and child_age
        branches_at_time = []
        for _, row in self.edgedf.iterrows():
            parent_age = row['parent_age']
            child_age = row['child_age']
            
            # Branch crosses this time if parent_age >= time >= child_age
            if parent_age >= time_mya >= child_age:
                branches_at_time.append(row['child_taxid'])
        
        return branches_at_time
    
    def create_pruned_tree_at_time(self, time_mya, working_dag=None, working_edgedf=None):
        """
        Create pruned tree at time T by cutting branches and deleting descendants.
        
        CRITICAL: ALL tips must be lined up exactly at time T (same distance from root).
        
        Args:
            time_mya: The time point (in Mya) at which to prune the tree
            working_dag: Optional DAG to prune (for progressive pruning). If None, uses self.dag.
            working_edgedf: Optional edgedf to prune (for progressive pruning). If None, uses self.edgedf.
        
        Result: ONE connected component, all tips exactly at time T.
        """
        
        # Use provided structures or defaults
        if working_dag is None:
            working_dag = self.dag
        if working_edgedf is None:
            working_edgedf = self.edgedf
        
        # Step 1: Identify branches crossing T and branches entirely older than T
        tips_at_this_time = set()
        branches_to_include = []
        
        for idx, row in working_edgedf.iterrows():
            parent_age = row['parent_age']
            child_age = row['child_age']
            
            # Case 1: Branch crosses T (will become a tip)
            if parent_age >= time_mya >= child_age:
                tips_at_this_time.add(row['child_taxid'])
                branches_to_include.append(idx)
            # Case 2: Branch entirely older than T (keep unchanged)
            elif child_age >= time_mya:
                branches_to_include.append(idx)
            # Case 3: Branch entirely younger than T (exclude)
        
        edgedf_subset = working_edgedf.iloc[branches_to_include].copy()
        
        # Step 2: Build pruned DAG with adjusted branch lengths
        pruned_dag = {}
        
        for _, row in working_edgedf.iterrows():
            parent = row['parent_taxid']
            child = row['child_taxid']
            parent_age = row['parent_age']
            child_age = row['child_age']
            original_length = row['branch_length']
            
            # Skip if parent doesn't exist at T
            if parent_age < time_mya:
                continue
            
            # Initialize parent
            if parent not in pruned_dag:
                pruned_dag[parent] = {}
            
            # Case 1: Branch crosses T → CUT IT at T
            if child in tips_at_this_time:
                # Shorten branch to end exactly at T
                adjusted_length = parent_age - time_mya
                pruned_dag[parent][child] = adjusted_length
                # Child becomes tip with no children
                if child not in pruned_dag:
                    pruned_dag[child] = {}
            
            # Case 2: Branch entirely older than T → KEEP unchanged
            elif child_age >= time_mya:
                pruned_dag[parent][child] = original_length
            
            # Case 3: Branch entirely younger than T → DELETE
            # (implicit - we don't add it)
        
        return pruned_dag, tips_at_this_time, edgedf_subset

# Worker function for parallel time-slice weight computation (must be at module level for pickling)
def _compute_time_slice_weights_worker(args_tuple):
    """Worker function to compute phylogenetic weights for a single time slice."""
    time_slice, edgedf_copy, phylo_weighter_args = args_tuple
    
    # Recreate phylo_weighter in worker process
    from plot_branch_stats_vs_time import PhyloWeighting
    nodedf_copy, cache_dir = phylo_weighter_args
    worker_weighter = PhyloWeighting(edgedf_copy, nodedf_copy, cache_dir=cache_dir, verbose=False)
    
    # Build pruned tree at this time
    pruned_dag, tips_at_this_time, edgedf_at_time = worker_weighter.create_pruned_tree_at_time(
        time_slice, working_dag=worker_weighter.dag, working_edgedf=edgedf_copy
    )
    
    # Compute weights only for branches that cross this time
    if len(tips_at_this_time) > 0:
        tip_branches = edgedf_at_time[edgedf_at_time['child_taxid'].isin(tips_at_this_time)]
        weights_dict = worker_weighter.compute_branch_weights(
            tip_branches, time_slice=time_slice, pruned_dag=pruned_dag
        )
        return (time_slice, weights_dict)
    else:
        return (time_slice, {})

def parse_args():
    """
    We need two files:
      - the edge_information file.
      - the node information file
      - the chromosome information file
        - From this we want the median number of chromosomes for each taxid
      - the statsdf file.
      - the intensity of extinction file.
      - an optional list of taxids to omit, anything in the clade that has one of these taxids as a parent (or as the actual taxid) will be omitted.
      - a flag R to just rerun the incomplete parts of the analysis. Doesn't bother to output a node file.

    We will check that both files exist.
    """
    # Check if no arguments provided, show help
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    
    parser = argparse.ArgumentParser(
        description='Plot branch stats vs time',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''This script analyzes chromosome fusion/fission events along phylogenetic tree branches over evolutionary time.
It takes edge information (parent-child nodes with ages/branch lengths), node information, and a df_stats file
with fusion/dispersal events. It maps chromosome structural changes to specific branches in the phylogenetic tree,
analyzes the rate of these events over time (in million years), and optionally correlates events with
extinction/origination intensity data.''')
    parser.add_argument("-e", "--edge_information", type=str, required=True, 
                        help='The edge information file (TSV format). Must have columns: source, target, source_age, target_age, branch_length, source_ages, target_ages')
    parser.add_argument("-n", "--node_information", type=str, required=True, 
                        help='The node information file (TSV format). Contains information about each node in the phylogenetic tree including taxid, parent, dist_crown, dist_crown_plus_root')
    parser.add_argument("-s", "--statsdf", type=str, required=True, 
                        help='The df_stats file (TSV format). Contains fusion and dispersal events with taxid information to map events to tree branches')
    parser.add_argument("-S", "--suppress_plotting", action='store_true', 
                        help='Suppress plotting - only generate TSV output files without creating PDF visualizations')
    parser.add_argument("-i", "--intensity_of_extinction", type=str, 
                        help="The intensity of extinction file (TSV format, optional). Must have columns: 'Time (Ma)', 'Diversity All Genera', 'Diversity Short-Lived', 'Diversity Long-Lived', 'Diversity Well-Resolved', 'Extinction Intensity (%%)', 'Origination Intensity(%%)'")
    parser.add_argument("-o", "--omit_taxids", type=str, 
                        help='A comma-separated list of taxids to omit from analysis. Any clade that has one of these taxids as a parent or as the actual taxid will be excluded. Example: 9606,10090,7227')
    parser.add_argument("--analyze_single_clade", type=int,
                        help='Analyze a single clade (taxid) with optional subclade exclusions. Use with --exclude_subclades. Output will be saved to custom_clade_analyses/ with descriptive naming.')
    parser.add_argument("--exclude_subclades", type=str,
                        help='Comma-separated list of taxids to exclude when using --analyze_single_clade. Example: 32443,7777 to exclude Teleostei and Chondrichthyes from Vertebrata analysis.')
    parser.add_argument("-O", "--outdir", type=str, default="branch_stats_output",
                        help='Output directory for all generated files (default: branch_stats_output)')
    parser.add_argument("-P", "--phylogenetic_weighting", action='store_true',
                        help='Apply phylogenetic weighting to account for oversampled clades. Branches in species-rich recent radiations are downweighted relative to phylogenetically isolated lineages.')
    parser.add_argument("-t", "--threads", type=int, default=None,
                        help='Number of parallel threads for phylogenetic weighting (default: auto-detect, max 16)')
    args = parser.parse_args()

    if not os.path.exists(args.edge_information):
        raise ValueError('The edge information file does not exist')

    if not os.path.exists(args.node_information):
        raise ValueError('The node information file does not exist')

    if not os.path.exists(args.statsdf):
        raise ValueError('The statsdf file does not exist')

    # we can optionally plot the intensity of extinction
    if args.intensity_of_extinction:
        if not os.path.exists(args.intensity_of_extinction):
            raise ValueError('The intensity of extinction file does not exist')

    # parse the omit_taxids into a list
    if args.omit_taxids:
        if type(args.omit_taxids) == str:
            args.omit_taxids = [int(x) for x in args.omit_taxids.split(',')]
    
    # parse the exclude_subclades into a list
    if args.exclude_subclades:
        if type(args.exclude_subclades) == str:
            args.exclude_subclades = [int(x) for x in args.exclude_subclades.split(',')]
    else:
        args.exclude_subclades = []
    
    # validate special mode arguments
    if args.analyze_single_clade and not args.exclude_subclades:
        print("WARNING: --analyze_single_clade specified without --exclude_subclades. This will analyze the full clade without exclusions.")
    if args.exclude_subclades and not args.analyze_single_clade:
        raise ValueError("--exclude_subclades requires --analyze_single_clade to be specified")
    
    return args

def plot_phylogenetic_weight_diagnostics(edgedf, outdir, taxid):
    """
    Create diagnostic plots showing phylogenetic weight distribution across the full tree.
    Demonstrates overall phylogenetic balance.
    
    Parameters:
    - edgedf: DataFrame with phylo_weight, branch_length, parent_age, child_age columns
    - outdir: Output directory for plot
    - taxid: Identifier for output filename (typically 'ALL' for global diagnostic)
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Phylogenetic Weight Diagnostics ({taxid})', fontsize=14, fontweight='bold')
    
    weights = edgedf['phylo_weight'].values
    branch_lengths = edgedf['branch_length'].values
    parent_ages = edgedf['parent_age'].values
    
    # Top-left: Histogram of weights
    ax = axes[0, 0]
    ax.hist(weights, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(weights.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean = {weights.mean():.3f}')
    ax.axvline(np.median(weights), color='orange', linestyle='--', linewidth=2, label=f'Median = {np.median(weights):.3f}')
    ax.set_xlabel('Phylogenetic Weight', fontsize=11)
    ax.set_ylabel('Number of Branches', fontsize=11)
    ax.set_title('Weight Distribution', fontsize=12, fontweight='bold')
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(alpha=0.3)
    
    # Add statistics text box
    stats_text = f'N = {len(weights):,}\nMin = {weights.min():.3f}\nMax = {weights.max():.3f}\nStd = {weights.std():.3f}\nCV = {weights.std()/weights.mean():.3f}'
    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Top-right: Cumulative distribution
    ax = axes[0, 1]
    sorted_weights = np.sort(weights)
    cumulative = np.arange(1, len(sorted_weights) + 1) / len(sorted_weights)
    ax.plot(sorted_weights, cumulative, color='steelblue', linewidth=2)
    ax.axvline(1.0, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Weight = 1.0')
    ax.set_xlabel('Phylogenetic Weight', fontsize=11)
    ax.set_ylabel('Cumulative Probability', fontsize=11)
    ax.set_title('Cumulative Distribution', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Bottom-left: Weight vs Branch Length (colored by age)
    ax = axes[1, 0]
    scatter = ax.scatter(branch_lengths, weights, c=parent_ages, cmap='viridis', 
                        alpha=0.5, s=10, edgecolors='none')
    ax.axhline(1.0, color='red', linestyle='--', linewidth=1, alpha=0.7)
    ax.set_xlabel('Branch Length (My)', fontsize=11)
    ax.set_ylabel('Phylogenetic Weight', fontsize=11)
    ax.set_title('Weight vs Branch Length', fontsize=12, fontweight='bold')
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Parent Age (Mya)', fontsize=10)
    ax.grid(alpha=0.3)
    
    # Bottom-right: Weight vs Parent Age
    ax = axes[1, 1]
    ax.scatter(parent_ages, weights, alpha=0.3, s=10, color='steelblue', edgecolors='none')
    ax.axhline(1.0, color='red', linestyle='--', linewidth=1, alpha=0.7, label='Mean weight')
    ax.set_xlabel('Parent Age (Mya)', fontsize=11)
    ax.set_ylabel('Phylogenetic Weight', fontsize=11)
    ax.set_title('Weight vs Time', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    outfile = os.path.join(outdir, f'{taxid}_phylogenetic_weights_diagnostic.pdf')
    plt.savefig(outfile, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved phylogenetic weight diagnostic plot: {outfile}")

def plot_phylogenetic_weight_temporal_heatmap(time_slice_weights, outdir, oldest_age, output_prefix='phylogenetic_weights'):
    """
    Create a heatmap showing how phylogenetic weight distribution evolves over time.
    Each vertical slice shows the weight distribution at that time point.
    
    Parameters:
    - time_slice_weights: Dict mapping time_slice (Mya) -> dict of {taxid: weight}
    - outdir: Output directory for plot
    - oldest_age: Maximum age in the tree
    - output_prefix: Prefix for output filename (e.g., 'ALL' or clade name)
    """
    import matplotlib.colors as mcolors
    
    print("Generating temporal weight distribution heatmap...")
    
    # Define weight bins
    weight_bins = np.arange(0.5, 2.51, 0.05)  # From 0.5 to 2.5 in 0.05 increments
    weight_bin_centers = (weight_bins[:-1] + weight_bins[1:]) / 2
    
    # Create 2D histogram: rows = weight bins, columns = time slices
    time_slices = sorted(time_slice_weights.keys())
    heatmap_data = np.zeros((len(weight_bin_centers), len(time_slices)))
    
    for j, time_slice in enumerate(time_slices):
        weights_dict = time_slice_weights[time_slice]
        if weights_dict:
            weights = list(weights_dict.values())
            counts, _ = np.histogram(weights, bins=weight_bins)
            heatmap_data[:, j] = counts
    
    # Flip data horizontally so present (0 Mya) is on right and oldest is on left
    heatmap_data = np.flip(heatmap_data, axis=1)
    
    # Create log2-space histogram for bottom panel
    # Use evenly-spaced bins in log2 space
    log2_weight_bins = np.arange(np.log2(0.5), np.log2(2.5) + 0.05, 0.05)
    log2_heatmap_data = np.zeros((len(log2_weight_bins) - 1, len(time_slices)))
    
    for j, time_slice in enumerate(time_slices):
        weights_dict = time_slice_weights[time_slice]
        if weights_dict:
            weights = list(weights_dict.values())
            # Transform weights to log2 space before binning
            log2_weights = np.log2(weights)
            counts, _ = np.histogram(log2_weights, bins=log2_weight_bins)
            log2_heatmap_data[:, j] = counts
    
    # Flip log2 data horizontally to match top panel
    log2_heatmap_data = np.flip(log2_heatmap_data, axis=1)
    
    # Create figure with two panels (linear and log2) with white background
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12), sharex=True, facecolor='white')
    
    # Top panel: Log scale (better for diffuse counts)
    im1 = ax1.imshow(heatmap_data, aspect='auto', origin='lower', 
                     cmap='YlOrRd', 
                     norm=mcolors.LogNorm(vmin=1, vmax=heatmap_data.max() + 1),
                     extent=[oldest_age, 0, weight_bins[0], weight_bins[-1]])
    
    cbar1 = plt.colorbar(im1, ax=ax1, label='Number of Tips')
    ax1.axhline(1.0, color='blue', linestyle='--', linewidth=2, alpha=0.7, label='Weight = 1.0')
    ax1.set_ylabel('Phylogenetic Weight', fontsize=12, fontweight='bold')
    ax1.set_title('Phylogenetic Weight Distribution Over Time (Log Scale)', 
                  fontsize=12, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(False)
    ax1.set_facecolor('white')
    
    # Bottom panel: Log2 scale (data rebinned in log2 space)
    im2 = ax2.imshow(log2_heatmap_data, aspect='auto', origin='lower', 
                     cmap='YlOrRd',  # Same colormap as top panel
                     norm=mcolors.LogNorm(vmin=1, vmax=log2_heatmap_data.max() + 1),
                     extent=[oldest_age, 0, log2_weight_bins[0], log2_weight_bins[-1]])
    
    cbar2 = plt.colorbar(im2, ax=ax2, label='Number of Tips')
    ax2.axhline(0.0, color='blue', linestyle='--', linewidth=2, alpha=0.7, label='log₂(Weight) = 0 (Weight = 1.0)')
    ax2.set_xlabel('Time (Mya)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('log₂(Phylogenetic Weight)', fontsize=12, fontweight='bold')
    ax2.set_title('Phylogenetic Weight Distribution Over Time (Log₂ Scale)', 
                  fontsize=12, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(False)
    ax2.set_facecolor('white')
    
    # Add text annotation to top panel
    ax1.text(0.02, 0.98, 
            f'Time slices: {len(time_slices)}\nWeight range: [{weight_bins[0]:.1f}, {weight_bins[-1]:.1f}]',
            transform=ax1.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    outfile = os.path.join(outdir, f'{output_prefix}_temporal_heatmap.pdf')
    plt.savefig(outfile, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved temporal weight heatmap: {outfile}")

def plot_phylogenetic_weights_temporal(edgedf_subset, outdir, taxid):
    """
    Create temporal analysis showing how phylogenetic weights vary across time periods.
    Shows whether phylogenetic correction is consistent across deep vs shallow branches.
    
    Parameters:
    - edgedf_subset: Filtered DataFrame for this clade/subset
    - outdir: Output directory for plot
    - taxid: Identifier for output filename
    """
    if len(edgedf_subset) == 0:
        print(f"Skipping temporal weight plot for {taxid}: no branches")
        return
    
    # Define time bins (100 Mya intervals)
    max_age = edgedf_subset['parent_age'].max()
    if max_age <= 0:
        print(f"Skipping temporal weight plot for {taxid}: all ages <= 0")
        return
    
    bin_size = 100
    n_bins = int(np.ceil(max_age / bin_size))
    bins = [(i * bin_size, (i + 1) * bin_size) for i in range(n_bins)]
    
    # Assign branches to time bins
    bin_labels = []
    bin_weights = []
    bin_counts = []
    
    for bin_start, bin_end in bins:
        mask = (edgedf_subset['parent_age'] >= bin_start) & (edgedf_subset['parent_age'] < bin_end)
        weights_in_bin = edgedf_subset.loc[mask, 'phylo_weight'].values
        
        if len(weights_in_bin) > 0:
            bin_labels.append(f"{bin_start}-{bin_end}")
            bin_weights.append(weights_in_bin)
            bin_counts.append(len(weights_in_bin))
    
    if len(bin_weights) == 0:
        print(f"Skipping temporal weight plot for {taxid}: no bins with data")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(10, len(bin_labels) * 1.5), 6))
    
    # Violin plots
    positions = np.arange(len(bin_labels))
    parts = ax.violinplot(bin_weights, positions=positions, widths=0.7,
                          showmeans=True, showmedians=True)
    
    # Color the violin plots
    for pc in parts['bodies']:
        pc.set_facecolor('steelblue')
        pc.set_alpha(0.6)
    
    # Overlay scatter points
    for i, weights in enumerate(bin_weights):
        x = np.random.normal(i, 0.04, size=len(weights))
        ax.scatter(x, weights, alpha=0.3, s=20, color='darkblue')
    
    # Add horizontal line at mean=1.0
    ax.axhline(1.0, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Mean = 1.0')
    
    # Labels and formatting
    ax.set_xticks(positions)
    ax.set_xticklabels(bin_labels, rotation=45, ha='right')
    ax.set_xlabel('Age Bin (Mya)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Phylogenetic Weight', fontsize=12, fontweight='bold')
    ax.set_title(f'Temporal Variation in Phylogenetic Weights ({taxid})', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3, axis='y')
    
    # Add statistics panel
    cv_values = [np.std(w) / np.mean(w) if len(w) > 0 else 0 for w in bin_weights]
    stats_text = "\n".join([f"{bin_labels[i]}: n={bin_counts[i]:,}, CV={cv_values[i]:.3f}" 
                             for i in range(len(bin_labels))])
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
            family='monospace')
    
    plt.tight_layout()
    outfile = os.path.join(outdir, f'{taxid}_phylogenetic_weights_temporal.pdf')
    plt.savefig(outfile, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved temporal weight analysis plot: {outfile}")

def plot_weighted_vs_unweighted_comparison(result_df, result_df_unweighted, outdir, output_prefix='ALL'):
    """
    Compare weighted vs unweighted rates to verify phylogenetic weighting is working.
    
    Creates a 3-panel figure showing:
    - Top: Fusion rates (weighted vs unweighted overlaid)
    - Middle: Dispersal rates (weighted vs unweighted overlaid)
    - Bottom: Difference plot (weighted - unweighted)
    
    Args:
        output_prefix: Prefix for output filename (default 'ALL')
    """
    format_matplotlib()
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 9))
    
    # Merge dataframes on age to ensure alignment (they may have different time slices)
    merged_df = pd.merge(result_df[['age', 'fusion_rate_at_this_age_mean', 'dispersal_rate_at_this_age_mean']], 
                         result_df_unweighted[['age', 'fusion_rate_at_this_age_mean', 'dispersal_rate_at_this_age_mean']], 
                         on='age', suffixes=('_weighted', '_unweighted'))
    
    # Prepare data - use mean rates (not ratios) to show actual weighting effect
    ages = merged_df['age'].values
    fusion_weighted = merged_df['fusion_rate_at_this_age_mean_weighted'].values
    fusion_unweighted = merged_df['fusion_rate_at_this_age_mean_unweighted'].values
    dispersal_weighted = merged_df['dispersal_rate_at_this_age_mean_weighted'].values
    dispersal_unweighted = merged_df['dispersal_rate_at_this_age_mean_unweighted'].values
    
    # Panel 1: Fusion rates comparison
    ax = axes[0]
    ax.plot(ages, fusion_weighted, 'b-', linewidth=2, label='Weighted', alpha=0.8)
    ax.plot(ages, fusion_unweighted, 'b--', linewidth=2, label='Unweighted', alpha=0.6)
    ax.set_ylabel('Fusion Rate\n(events/My)', fontsize=11)
    ax.set_title('Phylogenetic Weighting Verification: Fusion Rates', fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(alpha=0.3)
    ax.set_xlim([ages.min(), ages.max()])
    
    # Panel 2: Dispersal rates comparison
    ax = axes[1]
    ax.plot(ages, dispersal_weighted, 'r-', linewidth=2, label='Weighted', alpha=0.8)
    ax.plot(ages, dispersal_unweighted, 'r--', linewidth=2, label='Unweighted', alpha=0.6)
    ax.set_ylabel('Dispersal Rate\n(events/My)', fontsize=11)
    ax.set_title('Dispersal Rates', fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(alpha=0.3)
    ax.set_xlim([ages.min(), ages.max()])
    
    # Panel 3: Difference plot (weighted - unweighted)
    ax = axes[2]
    fusion_diff = fusion_weighted - fusion_unweighted
    dispersal_diff = dispersal_weighted - dispersal_unweighted
    ax.plot(ages, fusion_diff, 'b-', linewidth=2, label='Fusion difference', alpha=0.7)
    ax.plot(ages, dispersal_diff, 'r-', linewidth=2, label='Dispersal difference', alpha=0.7)
    ax.axhline(0, color='black', linestyle=':', linewidth=1)
    ax.set_xlabel('Age (Ma)', fontsize=11)
    ax.set_ylabel('Rate Difference\n(Weighted - Unweighted)', fontsize=11)
    ax.set_title('Difference (Weighted - Unweighted)', fontsize=12, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(alpha=0.3)
    ax.set_xlim([ages.min(), ages.max()])
    
    # Add summary statistics
    fusion_max_diff = np.abs(fusion_diff).max()
    dispersal_max_diff = np.abs(dispersal_diff).max()
    fusion_mean_diff = np.mean(np.abs(fusion_diff))
    dispersal_mean_diff = np.mean(np.abs(dispersal_diff))
    
    stats_text = (f'Max |diff| - Fusion: {fusion_max_diff:.4f}, Dispersal: {dispersal_max_diff:.4f}\n'
                  f'Mean |diff| - Fusion: {fusion_mean_diff:.4f}, Dispersal: {dispersal_mean_diff:.4f}')
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    outfile = os.path.join(outdir, f'{output_prefix}_phylogenetic_weighting_verification.pdf')
    plt.savefig(outfile, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved weighted vs unweighted comparison plot: {outfile}")

def plot_clade_weight_distribution(edgedf, outdir, target_clades):
    """
    Show weight distributions by major taxonomic clade to verify oversampled clades are downweighted.
    
    Args:
        edgedf: DataFrame with phylogenetic weights and lineage information
        outdir: Output directory
        target_clades: List of clade names to analyze
    """
    format_matplotlib()
    
    # Define clade name to taxid mapping (common NCBI taxids)
    clade_taxids = {
        'Lepidoptera': 7088,
        'Cnidaria': 6073,
        'Hexapoda': 6960,
        'Spiralia': 1206794,
        'Vertebrata': 7742,
        'Bilateria': 33213,
        'Panarthropoda': 88770,
        'Arthropoda': 6656,
        'Brachiopoda': 7565,
        'Bryozoa': 7509,
        'Chaetognatha': 7067,
        'Chordata': 7711,
        'Ctenophora': 10197,
        'Cycliophora': 50939,
        'Echinodermata': 7586,
        'Entoprocta': 51939,
        'Gastrotrich': 51028,
        'Gnathostomulida': 51023,
        'Hemichordata': 7696,
        'Kinorhyncha': 51032,
        'Loricifera': 51027,
        'Mollusca': 6447,
        'Nematoda': 6231,
        'Nematomorpha': 51031,
        'Nemertea': 51293,  # Note: user wrote "Memertea" but correct spelling is "Nemertea"
        'Onychophora': 51233,
        'Orthonectida': 51024,
        'Phoronida': 51502,
        'Placozoa': 10226,
        'Platyhelminthes': 6157,
        'Porifera': 6040,
        'Priapulida': 51025,
        'Rotifera': 10190,
        'Tardigrada': 42242,
        'Xenacoelomorpha': 1206795,
        'Protostomia': 33317,
        'Deuterostomia': 33511,
    }
    
    # Identify which clades are present in the data
    clade_weights = {}
    clade_branch_counts = {}
    
    for clade_name, clade_taxid in clade_taxids.items():
        if clade_name not in target_clades:
            continue
        
        # Find branches where this clade is in the lineage
        mask = edgedf['child_lineage_list'].apply(lambda x: clade_taxid in x if isinstance(x, (list, set)) else False)
        clade_branches = edgedf[mask]
        
        if len(clade_branches) > 0:
            clade_weights[clade_name] = clade_branches['phylo_weight'].values
            clade_branch_counts[clade_name] = len(clade_branches)
    
    if len(clade_weights) == 0:
        print("No target clades found in the dataset")
        return
    
    # Sort clades by median weight
    sorted_clades = sorted(clade_weights.keys(), 
                          key=lambda x: np.median(clade_weights[x]))
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(6, len(sorted_clades) * 0.4)))
    
    # Create violin plot
    positions = range(len(sorted_clades))
    parts = ax.violinplot([clade_weights[c] for c in sorted_clades],
                          positions=positions,
                          vert=False,
                          widths=0.7,
                          showmeans=True,
                          showmedians=True)
    
    # Color by whether median weight is above or below 1.0
    for i, clade in enumerate(sorted_clades):
        median_weight = np.median(clade_weights[clade])
        color = 'red' if median_weight < 1.0 else 'steelblue'
        parts['bodies'][i].set_facecolor(color)
        parts['bodies'][i].set_alpha(0.6)
    
    # Set labels
    ax.set_yticks(positions)
    ax.set_yticklabels([f'{c} (n={clade_branch_counts[c]})' for c in sorted_clades], fontsize=9)
    ax.set_xlabel('Phylogenetic Weight', fontsize=11)
    ax.set_title('Phylogenetic Weights by Major Clade', fontsize=12, fontweight='bold')
    ax.axvline(1.0, color='black', linestyle='--', linewidth=1, label='Mean weight = 1.0')
    ax.grid(alpha=0.3, axis='x')
    ax.legend(fontsize=9)
    
    # Add text box explaining colors
    color_text = 'Red: Oversampled clades (downweighted)\nBlue: Undersampled clades (upweighted)'
    ax.text(0.98, 0.02, color_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    outfile = os.path.join(outdir, 'ALL_phylogenetic_weights_by_clade.pdf')
    plt.savefig(outfile, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved clade weight distribution plot: {outfile}")
    print(f"Found {len(sorted_clades)} clades: {', '.join(sorted_clades)}")

def plot_clade_weight_distribution_from_cache(phylo_weighter, edgedf, outdir, target_clades):
    """
    Show weight distributions by major taxonomic clade using cached time-slice weights at t=0.
    This version works with time-slicing mode.
    
    Args:
        phylo_weighter: PhyloWeighting object with cached weights
        edgedf: DataFrame with lineage information
        outdir: Output directory
        target_clades: List of clade names to analyze
    """
    format_matplotlib()
    
    # Load weights from cache at t=0 (present day, where all tips exist)
    import pickle
    cache_subdir = os.path.join(phylo_weighter.cache_dir, 'phylo_cache')
    cache_file = os.path.join(cache_subdir, f'edge_weights_{phylo_weighter.tree_hash}_t0.pkl')
    
    if not os.path.exists(cache_file):
        print(f"Cannot generate clade weight plot: cache file not found at {cache_file}")
        return
    
    with open(cache_file, 'rb') as f:
        weights_dict = pickle.load(f)
    
    # Define clade name to taxid mapping (common NCBI taxids)
    clade_taxids = {
        'Lepidoptera': 7088,
        'Cnidaria': 6073,
        'Hexapoda': 6960,
        'Spiralia': 1206794,
        'Vertebrata': 7742,
        'Bilateria': 33213,
        'Panarthropoda': 88770,
        'Arthropoda': 6656,
        'Brachiopoda': 7565,
        'Bryozoa': 7509,
        'Chaetognatha': 7067,
        'Chordata': 7711,
        'Ctenophora': 10197,
        'Cycliophora': 50939,
        'Echinodermata': 7586,
        'Entoprocta': 51939,
        'Gastrotrich': 51028,
        'Gnathostomulida': 51023,
        'Hemichordata': 7696,
        'Kinorhyncha': 51032,
        'Loricifera': 51027,
        'Mollusca': 6447,
        'Nematoda': 6231,
        'Nematomorpha': 51031,
        'Nemertea': 51293,
        'Onychophora': 51233,
        'Orthonectida': 51024,
        'Phoronida': 51502,
        'Placozoa': 10226,
        'Platyhelminthes': 6157,
        'Porifera': 6040,
        'Priapulida': 51025,
        'Rotifera': 10190,
        'Tardigrada': 42242,
        'Xenacoelomorpha': 1206795,
        'Protostomia': 33317,
        'Deuterostomia': 33511,
    }
    
    # Identify which clades are present and collect their weights
    clade_weights = {}
    clade_branch_counts = {}
    
    # Get tip branches (child_age == 0)
    tip_branches = edgedf[edgedf['child_age'] == 0].copy()
    
    for clade_name, clade_taxid in clade_taxids.items():
        if clade_name not in target_clades:
            continue
        
        # Find tip branches where this clade is in the lineage
        mask = tip_branches['child_lineage_list'].apply(lambda x: clade_taxid in x if isinstance(x, (list, set)) else False)
        clade_tips = tip_branches[mask]
        
        if len(clade_tips) > 0:
            # Get weights for these tips from weights_dict
            tip_weights = []
            for taxid in clade_tips['child_taxid']:
                if taxid in weights_dict:
                    tip_weights.append(weights_dict[taxid])
            
            if len(tip_weights) > 0:
                clade_weights[clade_name] = np.array(tip_weights)
                clade_branch_counts[clade_name] = len(tip_weights)
    
    if len(clade_weights) == 0:
        print("No target clades found in the dataset")
        return
    
    # Create violin plot
    fig, ax = plt.subplots(figsize=(10, max(8, len(clade_weights) * 0.5)))
    
    # Sort clades by mean weight
    sorted_clades = sorted(clade_weights.keys(), key=lambda c: np.mean(clade_weights[c]))
    positions = np.arange(len(sorted_clades))
    
    # Prepare data for violin plot
    plot_data = [clade_weights[c] for c in sorted_clades]
    
    # Create violin plots
    parts = ax.violinplot(plot_data, positions=positions, vert=False, widths=0.7,
                          showmeans=True, showmedians=True)
    
    # Color violins by whether they're over or under-weighted
    for i, pc in enumerate(parts['bodies']):
        mean_weight = np.mean(plot_data[i])
        if mean_weight > 1.0:
            pc.set_facecolor('lightcoral')  # Oversampled (downweighted)
        else:
            pc.set_facecolor('skyblue')  # Undersampled (upweighted)
        pc.set_alpha(0.6)
    
    # Set labels
    ax.set_yticks(positions)
    ax.set_yticklabels([f'{c} (n={clade_branch_counts[c]})' for c in sorted_clades], fontsize=9)
    ax.set_xlabel('Phylogenetic Weight', fontsize=11)
    ax.set_title('Phylogenetic Weights by Major Clade', fontsize=12, fontweight='bold')
    ax.axvline(1.0, color='black', linestyle='--', linewidth=1, label='Mean weight = 1.0')
    ax.grid(alpha=0.3, axis='x')
    ax.legend(fontsize=9)
    
    # Add text box explaining colors
    color_text = 'Red: Oversampled clades (downweighted)\nBlue: Undersampled clades (upweighted)'
    ax.text(0.98, 0.02, color_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    outfile = os.path.join(outdir, 'ALL_phylogenetic_weights_by_clade.pdf')
    plt.savefig(outfile, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved clade weight distribution plot: {outfile}")
    print(f"Found {len(sorted_clades)} clades: {', '.join(sorted_clades)}")

def plot_event_count_conservation(result_df, result_df_unweighted, edgedf, outdir, output_prefix='ALL'):
    """
    Verify that phylogenetic weighting preserves total event counts (normalization check).
    
    Shows three bars each for fusions and dispersals:
    1. Raw counts from branches
    2. Unweighted total across time bins
    3. Weighted total across time bins
    
    Args:
        output_prefix: Prefix for output filename (default 'ALL')
    """
    format_matplotlib()
    
    # Calculate raw counts from edgedf
    raw_fusions = edgedf['num_fusions_this_branch'].sum()
    raw_dispersals = edgedf['num_dispersals_this_branch'].sum()
    
    # Calculate totals from analysis (sum across all time bins)
    weighted_fusions = result_df['total_fusions_at_this_age'].sum()
    weighted_dispersals = result_df['total_dispersals_at_this_age'].sum()
    unweighted_fusions = result_df_unweighted['total_fusions_at_this_age'].sum()
    unweighted_dispersals = result_df_unweighted['total_dispersals_at_this_age'].sum()
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # Panel 1: Fusions
    ax = axes[0]
    bars = ax.bar(['Raw\nCounts', 'Unweighted\nAnalysis', 'Weighted\nAnalysis'],
                  [raw_fusions, unweighted_fusions, weighted_fusions],
                  color=['gray', 'skyblue', 'steelblue'],
                  alpha=0.8,
                  edgecolor='black')
    ax.set_ylabel('Total Fusion Events', fontsize=12)
    ax.set_title('Fusion Event Count Conservation', fontsize=13, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Add percentage differences
    pct_diff_unweighted = 100 * (unweighted_fusions - raw_fusions) / raw_fusions
    pct_diff_weighted = 100 * (weighted_fusions - raw_fusions) / raw_fusions
    stats_text = f'Unweighted: {pct_diff_unweighted:+.2f}%\nWeighted: {pct_diff_weighted:+.2f}%'
    ax.text(0.98, 0.50, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Panel 2: Dispersals
    ax = axes[1]
    bars = ax.bar(['Raw\nCounts', 'Unweighted\nAnalysis', 'Weighted\nAnalysis'],
                  [raw_dispersals, unweighted_dispersals, weighted_dispersals],
                  color=['gray', 'lightcoral', 'red'],
                  alpha=0.8,
                  edgecolor='black')
    ax.set_ylabel('Total Dispersal Events', fontsize=12)
    ax.set_title('Dispersal Event Count Conservation', fontsize=13, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Add percentage differences
    pct_diff_unweighted = 100 * (unweighted_dispersals - raw_dispersals) / raw_dispersals
    pct_diff_weighted = 100 * (weighted_dispersals - raw_dispersals) / raw_dispersals
    stats_text = f'Unweighted: {pct_diff_unweighted:+.2f}%\nWeighted: {pct_diff_weighted:+.2f}%'
    ax.text(0.98, 0.50, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    outfile = os.path.join(outdir, f'{output_prefix}_event_count_conservation.pdf')
    plt.savefig(outfile, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    print(f"Saved event count conservation plot: {outfile}")
    
    # Print summary
    print(f"\nEvent Count Conservation Summary:")
    print(f"  Fusions - Raw: {raw_fusions:.1f}, Unweighted: {unweighted_fusions:.1f}, Weighted: {weighted_fusions:.1f}")
    print(f"  Dispersals  - Raw: {raw_dispersals:.1f}, Unweighted: {unweighted_dispersals:.1f}, Weighted: {weighted_dispersals:.1f}")

def plot_fusions_per_branch_vs_time(outprefix, resultsdf_weighted, resultsdf_unweighted, intensity_of_extinction_filepath = None):
    """
    Makes a plot comparing phylogenetically weighted vs unweighted fusion/dispersal rates over time.
    
    Columns 0-1: Phylogenetically weighted rates (linear and log2)
    Columns 2-3: Unweighted rates (linear and log2)
    """
    # call the function to properly format the text
    format_matplotlib()

    # MAGIC NUMBERS
    fs = 8
    # make the axis limits go from -1200 to 0
    xmin = -900
    xmax = 0

    # Now make three plots, each top of one another.
    # The 0th row shows both fusions and fissions ratios together in the same axis.
    # The 1st row shows just the fusion ratio on the axis.
    # The 2nd row shows just the fission ratio on the axis.
    # Use two different panels, plot the fusions in the top panel, and the dispersals in the bottom panel. Color the fusions blue and the dispersals red.
    # We also make a second axis which is the log2 of the ratio of fusions to dispersals.

    # The third column uses the mean of the fusion rate at each age.
    # The fourth column uses the log2 of the fusion rate at each age.

    num_rows = 3
    if intensity_of_extinction_filepath is not None:
        num_rows = 5
        # read in the intensity of extinction file
        intensity_of_extinction_df = pd.read_csv(intensity_of_extinction_filepath, sep='\t')

    red = "#D22C16"
    blue = "#3054A3"

    fig, ax = plt.subplots(num_rows, 4, figsize=(12, 12))
    
    # Columns 0-1: Weighted, Columns 2-3: Unweighted
    for coli in [0,1,2,3]:
        # Select which dataframe to use
        if coli < 2:
            resultsdf = resultsdf_weighted
        else:
            resultsdf = resultsdf_unweighted
            
        # The age axis is the same in all of the plots
        for rowi in range(num_rows):
            ax[rowi][coli].set_xlabel("Million Years ago (Mya)")
        
        # Always plot fusions_ratio and dispersals_ratio (the biologically meaningful weighted metrics)
        yfusion = "fusions_ratio"
        yloss   = "dispersals_ratio"
        
        # Determine if we're doing the log2 of the values
        if coli % 2 == 0:
            # we use the native values
            ax[0][coli].plot(resultsdf["age"], resultsdf[yfusion], color=blue)
            ax[0][coli].plot(resultsdf["age"], resultsdf[yloss]  ,   color=red)
            ax[1][coli].plot(resultsdf["age"], resultsdf[yfusion], color=blue)
            ax[2][coli].plot(resultsdf["age"], resultsdf[yloss]  , color=red)
            
            # Add shaded confidence intervals for linear plots
            if 'fusion_rate_at_this_age_10th' in resultsdf.columns and 'fusion_rate_at_this_age_90th' in resultsdf.columns:
                ax[1][coli].fill_between(resultsdf["age"], 
                                         resultsdf['fusion_rate_at_this_age_10th'].astype('float64'),
                                         resultsdf['fusion_rate_at_this_age_90th'].astype('float64'),
                                         color=blue, alpha=0.2)
            if 'dispersal_rate_at_this_age_10th' in resultsdf.columns and 'dispersal_rate_at_this_age_90th' in resultsdf.columns:
                ax[2][coli].fill_between(resultsdf["age"],
                                         resultsdf['dispersal_rate_at_this_age_10th'].astype('float64'),
                                         resultsdf['dispersal_rate_at_this_age_90th'].astype('float64'),
                                         color=red, alpha=0.2)
            label_append = ""
        else:
            # we use the log2 of the values
            ax[0][coli].plot(resultsdf["age"], np.log2(resultsdf[yfusion].astype('float64')), color=blue )
            ax[0][coli].plot(resultsdf["age"], np.log2(resultsdf[yloss].astype('float64')  ), color=red  )
            ax[1][coli].plot(resultsdf["age"], np.log2(resultsdf[yfusion].astype('float64')), color=blue )
            ax[2][coli].plot(resultsdf["age"], np.log2(resultsdf[yloss].astype('float64')  ), color=red  )
            
            # Add shaded confidence intervals for log2 plots (filter out invalid values)
            if 'fusion_rate_at_this_age_10th' in resultsdf.columns and 'fusion_rate_at_this_age_90th' in resultsdf.columns:
                # Only plot where both quantiles are positive
                valid_mask = (resultsdf['fusion_rate_at_this_age_10th'] > 0) & (resultsdf['fusion_rate_at_this_age_90th'] > 0)
                if valid_mask.any():
                    ax[1][coli].fill_between(resultsdf["age"][valid_mask], 
                                             np.log2(resultsdf['fusion_rate_at_this_age_10th'][valid_mask].astype('float64')),
                                             np.log2(resultsdf['fusion_rate_at_this_age_90th'][valid_mask].astype('float64')),
                                             color=blue, alpha=0.2)
            if 'dispersal_rate_at_this_age_10th' in resultsdf.columns and 'dispersal_rate_at_this_age_90th' in resultsdf.columns:
                # Only plot where both quantiles are positive
                valid_mask = (resultsdf['dispersal_rate_at_this_age_10th'] > 0) & (resultsdf['dispersal_rate_at_this_age_90th'] > 0)
                if valid_mask.any():
                    ax[2][coli].fill_between(resultsdf["age"][valid_mask],
                                             np.log2(resultsdf['dispersal_rate_at_this_age_10th'][valid_mask].astype('float64')),
                                             np.log2(resultsdf['dispersal_rate_at_this_age_90th'][valid_mask].astype('float64')),
                                             color=red, alpha=0.2)
            label_append = " (log2)"

        # Set titles based on weighting
        if coli < 2:
            ylabel0 = "Changes/My"
            ylabel1 = "Fusions/My"
            ylabel2 = "Dispersals/My"
            ylabel3 = "Extinction\nIntensity (%)"
            ylabel4 = "Origination\nIntensity (%)"
            title0  = "Phylo-weighted: Fusions or Dispersals/My"
            title1  = "Phylo-weighted: Fusions/My"
            title2  = "Phylo-weighted: Dispersals/My"
            title3  = "Rohde & Muller (2005)\nExtinction Intensity"
            title4  = "Rohde & Muller (2005)\nOrigination Intensity"
        else:
            ylabel0 = "Changes/My"
            ylabel1 = "Fusions/My"
            ylabel2 = "Dispersals/My"
            ylabel3 = "Extinction\nIntensity (%)"
            ylabel4 = "Origination\nIntensity (%)"
            title0  = "Unweighted: Fusions or Dispersals/My"
            title1  = "Unweighted: Fusions/My"
            title2  = "Unweighted: Dispersals/My"
            title3  = "Rohde & Muller (2005)\nExtinction Intensity"
            title4  = "Rohde & Muller (2005)\nOrigination Intensity"

        ax[0][coli].set_title( title0                , fontsize = fs)
        ax[0][coli].set_ylabel(ylabel0 + label_append, fontsize = fs)
        # second, just fusions
        ax[1][coli].set_title( title1                , fontsize = fs)
        ax[1][coli].set_ylabel(ylabel1 + label_append, fontsize = fs)
        # third, just dispersals
        ax[2][coli].set_title( title2                , fontsize = fs)
        ax[2][coli].set_ylabel(ylabel2 + label_append, fontsize = fs)

        if intensity_of_extinction_filepath is not None:
            # This is intensity of extinction
            for i, row in intensity_of_extinction_df.iterrows():
                left_x = -1 * row['Time (Ma)']
                left_y = 0
                height = row['Extinction Intensity (%)']
                width = 1
                rectangle = patches.Rectangle((left_x, left_y), width, height, fill="#555555", edgecolor=None)
                ax[3][coli].add_patch(rectangle)
            # This is intensity of origination
            for i, row in intensity_of_extinction_df.iterrows():
                left_x = -1 * row['Time (Ma)']
                left_y = 0
                height = row['Origination Intensity (%)']
                width = 1
                rectangle = patches.Rectangle((left_x, left_y), width, height, fill="#555555", edgecolor=None)
                ax[4][coli].add_patch(rectangle)
            ax[3][coli].set_title( title3 , fontsize = fs)
            ax[3][coli].set_ylabel(ylabel3, fontsize = fs)
            ax[4][coli].set_title( title4 , fontsize = fs)
            ax[4][coli].set_ylabel(ylabel4, fontsize = fs)

            # WE NOW SCALE THE Y-AXES
            # get the values where 'Time (Ma)' is between xmin and xmax
            subdf = intensity_of_extinction_df[(intensity_of_extinction_df['Time (Ma)'] >= xmax) & (intensity_of_extinction_df['Time (Ma)'] <= -1 * xmin)]
            # the ylim is 1.1 times the maximum value
            ymax = 1.1 * subdf['Extinction Intensity (%)'].max()
            ax[3][coli].set_ylim(0, ymax)
            # Now do the same for origination intensity
            ymax = 1.1 * subdf['Origination Intensity (%)'].max()
            ax[4][coli].set_ylim(0, ymax)

        for axi in range(num_rows):
            ax[axi][coli].set_xlim(xmin, xmax)

        # change the fontsize of the axes and the titles
        for axi in range(num_rows):
            ax[axi][coli].tick_params(axis='both', which='major', labelsize=fs)
            ax[axi][coli].set_title( ax[axi][coli].get_title(),    fontsize=fs)
            ax[axi][coli].set_xlabel(ax[axi][coli].get_xlabel(),   fontsize=fs)
            ax[axi][coli].set_ylabel(ax[axi][coli].get_ylabel(),   fontsize=fs)

    # increase the horizontal and vertical space between the panels
    plt.subplots_adjust(hspace=0.6, wspace=0.6)

    # save the plot as a pdf
    pdfout = outprefix.rstrip(".pdf") + ".pdf"
    plt.savefig(pdfout, facecolor='white', edgecolor='none')
    # close the figure to free up memory
    plt.close(fig)

def _subdf_no_missing(df, col1name, col2name) -> pd.DataFrame:
    """
    Given two column names and a df, extract those two columns, then remove rows with missing values, or values that are -inf or inf.
    A view of the original df is fine since we do not modify the subdf.
    Returns a view of the original df.
    """
    statsdf = df[[col1name, col2name]]
    return statsdf[~statsdf.isin([np.nan, np.inf, -np.inf]).any(axis=1)]

def _full_sk_stats(df, col1name, col2name, spearman_or_kendall) -> tuple:
    """
    Given a df and two column names,
    Returns a 3-element tuple of correlation, the p-value, and n (sample size)
    """
    statsdf = _subdf_no_missing(df, col1name, col2name)
    if spearman_or_kendall == "spearman":
        corr, p = spearmanr(statsdf[col1name], statsdf[col2name])
    elif spearman_or_kendall == "kendall":
        corr, p = kendalltau(statsdf[col1name], statsdf[col2name])
    elif spearman_or_kendall == "pearson":
        statsdf = statsdf[[col1name, col2name]]
        statsdf[col1name] = np.log2(statsdf[col1name])
        statsdf[col2name] = np.log2(statsdf[col2name])
        statsdf = _subdf_no_missing(statsdf, col1name, col2name)
        if len(statsdf) >= 2:
            corr, p = pearsonr(statsdf[col1name], statsdf[col2name])
        else:
            corr, p = np.nan, np.nan
    else:
        raise ValueError(f"Unknown correlation type {spearman_or_kendall}")
    return (corr, p, len(statsdf))

def plot_intensity_of_extinction(outprefix, count_df, intensity_of_extinction_filepath, suppress_plotting = False, fontsize = 8):
    """
    This makes another pdf of the intensity of extinction plus a df of the counts.
    This dataframe should be extracted from the supplementary information of this paper: https://www.nature.com/articles/nature03339

    Rohde, Robert A., and Richard A. Muller. "Cycles in fossil diversity." Nature 434.7030 (2005): 208-210.
    """
    # read in the intensity of extinction file
    intensity_of_extinction_df = pd.read_csv(intensity_of_extinction_filepath, sep='\t')

    # here we will make 6 plots.
    # Top row is fusions and dispersals plotted together
    # Second row is fusions
    # Third row is dispersals
    # For the first column, the y-axis (ratio of changes) is in linear scale
    # For the second column, the y-axis (ratio of changes) is in log2 scale
    # The x-axis will be the intensity of extinction

    # first, we need to get the intensity of extinction for each age.
    mya_to_extinction_intensity  = {row['Time (Ma)']: row['Extinction Intensity (%)'] for i, row in intensity_of_extinction_df.iterrows()}
    mya_to_origination_intensity = {row['Time (Ma)']: row['Origination Intensity (%)'] for i, row in intensity_of_extinction_df.iterrows()}

    # map the ages, but if the ages aren't present, give a value of -1
    count_df["age_positive"] = count_df["age"].apply(lambda x: -1 * x)
    count_df["intensity_of_extinction"]  = count_df["age_positive"].map(mya_to_extinction_intensity).fillna(-1)
    count_df["intensity_of_origination"] = count_df["age_positive"].map(mya_to_origination_intensity).fillna(-1)

    # only plot the values where intensity of extinction is not -1
    plotdf = count_df[count_df["intensity_of_extinction"] != -1].copy()
    plotdf["extinction_intensity_normalized"]  = plotdf["intensity_of_extinction"] / 100
    plotdf["origination_intensity_normalized"] = plotdf["intensity_of_origination"] / 100
    # set the correct columns to type float64
    for thiscol in ["fusions_ratio", "dispersals_ratio", "fusion_rate_at_this_age_mean", "dispersal_rate_at_this_age_mean"]:
        plotdf[thiscol] = plotdf[thiscol].astype('float64')
    # We want to measure the spearman correlation between the intensity of extinction and the fusions and dispersals. Use log for the fusions and dispersals.
    #  - It is not appropriate to use Pearson correlation, because the data are not linear (the percent, ranging between 0-100 or 0-1, cannot be linearly related to the number of fusions or dispersals).
    # get the values that have no nans or inf values
    fusion_spearman_corr,   fusion_spearman_p,   fusion_spearman_n   = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusions_ratio", "spearman")
    fusion_kendalltau_corr, fusion_kendalltau_p, fusion_kendalltau_n = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusions_ratio", "kendall" )
    fusion_pearson_corr,    fusion_pearson_p,    fusion_pearson_n    = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusions_ratio", "pearson" )
    loss_spearman_corr,     loss_spearman_p,     loss_spearman_n     = _full_sk_stats(plotdf, "extinction_intensity_normalized", "dispersals_ratio",  "spearman" )
    loss_kendalltau_corr,   loss_kendalltau_p,   loss_kendalltau_n   = _full_sk_stats(plotdf, "extinction_intensity_normalized", "dispersals_ratio",  "kendall"  )
    loss_pearson_corr,      loss_pearson_p,      loss_pearson_n      = _full_sk_stats(plotdf, "extinction_intensity_normalized", "dispersals_ratio",  "pearson"  )
    statsresults = {}
    statsresults["extinction_fusion_numpbranch_spearman_r"]   = fusion_spearman_corr
    statsresults["extinction_fusion_numpbranch_spearman_p"]   = fusion_spearman_p
    statsresults["extinction_fusion_numpbranch_spearman_n"]   = fusion_spearman_n
    statsresults["extinction_fusion_numpbranch_kendalltau_r"]  = fusion_kendalltau_corr
    statsresults["extinction_fusion_numpbranch_kendalltau_p"]  = fusion_kendalltau_p
    statsresults["extinction_fusion_numpbranch_kendalltau_n"]  = fusion_kendalltau_n
    statsresults["extinction_fusion_numpbranch_pearson_r"]      = fusion_pearson_corr
    statsresults["extinction_fusion_numpbranch_pearson_p"]      = fusion_pearson_p
    statsresults["extinction_fusion_numpbranch_pearson_n"]      = fusion_pearson_n
    statsresults["extinction_losses_numpbranch_spearman_r"]   = loss_spearman_corr
    statsresults["extinction_losses_numpbranch_spearman_p"]   = loss_spearman_p
    statsresults["extinction_losses_numpbranch_spearman_n"]   = loss_spearman_n
    statsresults["extinction_losses_numpbranch_kendalltau_r"]  = loss_kendalltau_corr
    statsresults["extinction_losses_numpbranch_kendalltau_p"]  = loss_kendalltau_p
    statsresults["extinction_losses_numpbranch_kendalltau_n"]  = loss_kendalltau_n
    statsresults["extinction_losses_numpbranch_pearson_r"]      = loss_pearson_corr
    statsresults["extinction_losses_numpbranch_pearson_p"]      = loss_pearson_p
    statsresults["extinction_losses_numpbranch_pearson_n"]      = loss_pearson_n


    ofusion_spearman_corr,   ofusion_spearman_p,   ofusion_spearman_n   = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusions_ratio", "spearman")
    ofusion_kendalltau_corr, ofusion_kendalltau_p, ofusion_kendalltau_n = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusions_ratio", "kendall")
    ofusion_pearson_corr,    ofusion_pearson_p,    ofusion_pearson_n    = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusions_ratio", "pearson")
    oloss_spearman_corr,     oloss_spearman_p,     oloss_spearman_n     = _full_sk_stats(plotdf, "origination_intensity_normalized", "dispersals_ratio" , "spearman")
    oloss_kendalltau_corr,   oloss_kendalltau_p,   oloss_kendalltau_n   = _full_sk_stats(plotdf, "origination_intensity_normalized", "dispersals_ratio" , "kendall")
    oloss_pearson_corr,      oloss_pearson_p,      oloss_pearson_n      = _full_sk_stats(plotdf, "origination_intensity_normalized", "dispersals_ratio" , "pearson")

    statsresults["origination_fusion_numpbranch_spearman_r"]   = ofusion_spearman_corr
    statsresults["origination_fusion_numpbranch_spearman_p"]   = ofusion_spearman_p
    statsresults["origination_fusion_numpbranch_spearman_n"]   = ofusion_spearman_n
    statsresults["origination_fusion_numpbranch_kendalltau_r"]  = ofusion_kendalltau_corr
    statsresults["origination_fusion_numpbranch_kendalltau_p"]  = ofusion_kendalltau_p
    statsresults["origination_fusion_numpbranch_kendalltau_n"]  = ofusion_kendalltau_n
    statsresults["origination_fusion_numpbranch_pearson_r"]      = ofusion_pearson_corr
    statsresults["origination_fusion_numpbranch_pearson_p"]      = ofusion_pearson_p
    statsresults["origination_fusion_numpbranch_pearson_n"]      = ofusion_pearson_n
    statsresults["origination_losses_numpbranch_spearman_r"]   = oloss_spearman_corr
    statsresults["origination_losses_numpbranch_spearman_p"]   = oloss_spearman_p
    statsresults["origination_losses_numpbranch_spearman_n"]   = oloss_spearman_n
    statsresults["origination_losses_numpbranch_kendalltau_r"]  = oloss_kendalltau_corr
    statsresults["origination_losses_numpbranch_kendalltau_p"]  = oloss_kendalltau_p
    statsresults["origination_losses_numpbranch_kendalltau_n"]  = oloss_kendalltau_n
    statsresults["origination_losses_numpbranch_pearson_r"]      = oloss_pearson_corr
    statsresults["origination_losses_numpbranch_pearson_p"]      = oloss_pearson_p
    statsresults["origination_losses_numpbranch_pearson_n"]      = oloss_pearson_n

    # Now these are on the rates at this age mean
    rate_fusion_spearman_corr,   rate_fusion_spearman_p,   rate_fusion_spearman_n   = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusion_rate_at_this_age_mean", "spearman")
    rate_fusion_kendalltau_corr, rate_fusion_kendalltau_p, rate_fusion_kendalltau_n  = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusion_rate_at_this_age_mean", "kendall" )
    rate_fusion_pearson_corr,    rate_fusion_pearson_p,    rate_fusion_pearson_n      = _full_sk_stats(plotdf, "extinction_intensity_normalized", "fusion_rate_at_this_age_mean", "pearson" )
    rate_loss_spearman_corr,     rate_loss_spearman_p,     rate_loss_spearman_n     = _full_sk_stats(plotdf, "extinction_intensity_normalized", "dispersal_rate_at_this_age_mean",   "spearman" )
    rate_loss_kendalltau_corr,   rate_loss_kendalltau_p,   rate_loss_kendalltau_n    = _full_sk_stats(plotdf, "extinction_intensity_normalized", "dispersal_rate_at_this_age_mean",   "kendall"  )
    rate_loss_pearson_corr,      rate_loss_pearson_p,      rate_loss_pearson_n        = _full_sk_stats(plotdf, "extinction_intensity_normalized", "dispersal_rate_at_this_age_mean",   "pearson"  )
    statsresults["extinction_fusion_ratepmya_spearman_r"]   = rate_fusion_spearman_corr
    statsresults["extinction_fusion_ratepmya_spearman_p"]   = rate_fusion_spearman_p
    statsresults["extinction_fusion_ratepmya_spearman_n"]   = rate_fusion_spearman_n
    statsresults["extinction_fusion_ratepmya_kendalltau_r"]  = rate_fusion_kendalltau_corr
    statsresults["extinction_fusion_ratepmya_kendalltau_p"]  = rate_fusion_kendalltau_p
    statsresults["extinction_fusion_ratepmya_kendalltau_n"]  = rate_fusion_kendalltau_n
    statsresults["extinction_fusion_ratepmya_pearson_r"]      = rate_fusion_pearson_corr
    statsresults["extinction_fusion_ratepmya_pearson_p"]      = rate_fusion_pearson_p
    statsresults["extinction_fusion_ratepmya_pearson_n"]      = rate_fusion_pearson_n
    statsresults["extinction_losses_ratepmya_spearman_r"]   = rate_loss_spearman_corr
    statsresults["extinction_losses_ratepmya_spearman_p"]   = rate_loss_spearman_p
    statsresults["extinction_losses_ratepmya_spearman_n"]   = rate_loss_spearman_n
    statsresults["extinction_losses_ratepmya_kendalltau_r"]  = rate_loss_kendalltau_corr
    statsresults["extinction_losses_ratepmya_kendalltau_p"]  = rate_loss_kendalltau_p
    statsresults["extinction_losses_ratepmya_kendalltau_n"]  = rate_loss_kendalltau_n
    statsresults["extinction_losses_ratepmya_pearson_r"]      = rate_loss_pearson_corr
    statsresults["extinction_losses_ratepmya_pearson_p"]      = rate_loss_pearson_p
    statsresults["extinction_losses_ratepmya_pearson_n"]      = rate_loss_pearson_n


    rate_ofusion_spearman_corr,   rate_ofusion_spearman_p,   rate_ofusion_spearman_n   = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusion_rate_at_this_age_mean", "spearman")
    rate_ofusion_kendalltau_corr, rate_ofusion_kendalltau_p, rate_ofusion_kendalltau_n  = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusion_rate_at_this_age_mean", "kendall")
    rate_ofusion_pearson_corr,    rate_ofusion_pearson_p,    rate_ofusion_pearson_n      = _full_sk_stats(plotdf, "origination_intensity_normalized", "fusion_rate_at_this_age_mean", "pearson")
    rate_oloss_spearman_corr,     rate_oloss_spearman_p,     rate_oloss_spearman_n     = _full_sk_stats(plotdf, "origination_intensity_normalized", "dispersal_rate_at_this_age_mean",   "spearman")
    rate_oloss_kendalltau_corr,   rate_oloss_kendalltau_p,   rate_oloss_kendalltau_n    = _full_sk_stats(plotdf, "origination_intensity_normalized", "dispersal_rate_at_this_age_mean",   "kendall")
    rate_oloss_pearson_corr,      rate_oloss_pearson_p,      rate_oloss_pearson_n        = _full_sk_stats(plotdf, "origination_intensity_normalized", "dispersal_rate_at_this_age_mean",   "pearson")

    statsresults["origination_fusion_ratepmya_spearman_r"]   = rate_ofusion_spearman_corr
    statsresults["origination_fusion_ratepmya_spearman_p"]   = rate_ofusion_spearman_p
    statsresults["origination_fusion_ratepmya_spearman_n"]   = rate_ofusion_spearman_n
    statsresults["origination_fusion_ratepmya_kendalltau_r"]  = rate_ofusion_kendalltau_corr
    statsresults["origination_fusion_ratepmya_kendalltau_p"]  = rate_ofusion_kendalltau_p
    statsresults["origination_fusion_ratepmya_kendalltau_n"]  = rate_ofusion_kendalltau_n
    statsresults["origination_fusion_ratepmya_pearson_r"]      = rate_ofusion_pearson_corr
    statsresults["origination_fusion_ratepmya_pearson_p"]      = rate_ofusion_pearson_p
    statsresults["origination_fusion_ratepmya_pearson_n"]      = rate_ofusion_pearson_n
    statsresults["origination_losses_ratepmya_spearman_r"]   = rate_oloss_spearman_corr
    statsresults["origination_losses_ratepmya_spearman_p"]   = rate_oloss_spearman_p
    statsresults["origination_losses_ratepmya_spearman_n"]   = rate_oloss_spearman_n
    statsresults["origination_losses_ratepmya_kendalltau_r"]  = rate_oloss_kendalltau_corr
    statsresults["origination_losses_ratepmya_kendalltau_p"]  = rate_oloss_kendalltau_p
    statsresults["origination_losses_ratepmya_kendalltau_n"]  = rate_oloss_kendalltau_n
    statsresults["origination_losses_ratepmya_pearson_r"]      = rate_oloss_pearson_corr
    statsresults["origination_losses_ratepmya_pearson_p"]      = rate_oloss_pearson_p
    statsresults["origination_losses_ratepmya_pearson_n"]      = rate_oloss_pearson_n

    def gen_stat_string(spearman_r, spearman_p, spearman_n, kendall_r, kendall_p, kendall_n, pearson_r, pearson_p, pearson_n) -> str:
        """
        generates a string with the statistics
        """
        outstring  = f"Spearman r={spearman_r:.2f}\n"
        outstring += f"Spearman p={spearman_p:.4e}\n"
        outstring += f"Spearman n={spearman_n}\n\n"
        outstring += f"Kendall r={kendall_r:.2f}\n"
        outstring += f"Kendall p={kendall_p:.4e}\n"
        outstring += f"Kendall n={kendall_n}\n\n"
        outstring += f"Pearson (log/log) r={pearson_r:.2f}\n"
        outstring += f"Pearson (log/log) p={pearson_p:.4e}\n"
        outstring += f"Pearson (log/log) n={pearson_n}"
        return outstring

    if not suppress_plotting:
        # EXTINCTION INTENSITY
        # now make the plots
        red = "#D22C16"
        blue = "#3054A3"
        num_rows = 6
        fig, ax = plt.subplots(num_rows, 4, figsize=(20, 20))
        fontsize = 8
        for coli in [0,1,2,3]:
            ylabel1 = "Changes/million years"
            ylabel2 = "Fusions/million years"
            ylabel3 = "Dispersals/million years"
            ylabel4 = "Changes/million years"
            ylabel5 = "Fusions/million years"
            ylabel6 = "Dispersals/million years"

            yfusiontop = plotdf["fusions_ratio"]
            yfusionbot = plotdf["fusion_rate_at_this_age_mean"]
            ylosstop   = plotdf["dispersals_ratio"]
            ylossbot   = plotdf["dispersal_rate_at_this_age_mean"]
            if coli < 2:
                xlabel = "Extinction Intensity (%)"
                x      = plotdf["intensity_of_extinction"]
                if coli % 2 == 0:
                    title1 = "Fusions or Dispersals/My vs Extinction Intensity"
                    title2 = "Fusions/My vs Extinction Intensity"
                    title3 = "Dispersals/My vs Extinction Intensity"
                    title4 = "Fusions or Dispersals/My vs Extinction Intensity"
                    title5 = "Fusions/My vs Extinction Intensity"
                    title6 = "Dispersals/My vs Extinction Intensity"
                else:
                    title1 = "Fusions or Dispersals/My vs Extinction Intensity (log2)"
                    title2 = "Fusions/My vs Extinction Intensity (log2)"
                    title3 = "Dispersals/My vs Extinction Intensity (log2)"
                    title4 = "Fusions or Dispersals/My vs Extinction Intensity (log2)"
                    title5 = "Fusions/My vs Extinction Intensity (log2)"
                    title6 = "Dispersals/My vs Extinction Intensity (log2)"
                stattext1 = gen_stat_string(fusion_spearman_corr,      fusion_spearman_p,      fusion_spearman_n,      fusion_kendalltau_corr,      fusion_kendalltau_p,      fusion_kendalltau_n     ,      fusion_pearson_corr,      fusion_pearson_p,      fusion_pearson_n)
                stattext2 = gen_stat_string(loss_spearman_corr,        loss_spearman_p,        loss_spearman_n,        loss_kendalltau_corr,        loss_kendalltau_p,        loss_kendalltau_n       ,        loss_pearson_corr,        loss_pearson_p,        loss_pearson_n)
                stattext3 = gen_stat_string(rate_fusion_spearman_corr, rate_fusion_spearman_p, rate_fusion_spearman_n, rate_fusion_kendalltau_corr, rate_fusion_kendalltau_p, rate_fusion_kendalltau_n, rate_fusion_pearson_corr, rate_fusion_pearson_p, rate_fusion_pearson_n)
                stattext4 = gen_stat_string(rate_loss_spearman_corr,   rate_loss_spearman_p,   rate_loss_spearman_n,   rate_loss_kendalltau_corr,   rate_loss_kendalltau_p,   rate_loss_kendalltau_n  ,   rate_loss_pearson_corr,   rate_loss_pearson_p,   rate_loss_pearson_n)
            else:
                xlabel = "Origination Intensity (%)"
                x      = plotdf["intensity_of_origination"]
                if coli % 2 == 0:
                    title1 = "Fusions or Dispersals/My vs Origination Intensity"
                    title2 = "Fusions/My vs Origination Intensity"
                    title3 = "Dispersals/My vs Origination Intensity"
                    title4 = "Fusions or Dispersals/My vs Origination Intensity"
                    title5 = "Fusions/My vs Origination Intensity"
                    title6 = "Dispersals/My vs Origination Intensity"
                else:
                    title1 = "Fusions or Dispersals/My vs Origination Intensity (log2)"
                    title2 = "Fusions/My vs Origination Intensity (log2)"
                    title3 = "Dispersals/My vs Origination Intensity (log2)"
                    title4 = "Fusions or Dispersals/My vs Origination Intensity (log2)"
                    title5 = "Fusions/My vs Origination Intensity (log2)"
                    title6 = "Dispersals/My vs Origination Intensity (log2)"
                stattext1 = gen_stat_string(ofusion_spearman_corr,      ofusion_spearman_p,      ofusion_spearman_n,      ofusion_kendalltau_corr,      ofusion_kendalltau_p,      ofusion_kendalltau_n     ,      ofusion_pearson_corr,      ofusion_pearson_p,      ofusion_pearson_n )
                stattext2 = gen_stat_string(oloss_spearman_corr,        oloss_spearman_p,        oloss_spearman_n,        oloss_kendalltau_corr,        oloss_kendalltau_p,        oloss_kendalltau_n       ,        oloss_pearson_corr,        oloss_pearson_p,        oloss_pearson_n )
                stattext3 = gen_stat_string(rate_ofusion_spearman_corr, rate_ofusion_spearman_p, rate_ofusion_spearman_n, rate_ofusion_kendalltau_corr, rate_ofusion_kendalltau_p, rate_ofusion_kendalltau_n, rate_ofusion_pearson_corr, rate_ofusion_pearson_p, rate_ofusion_pearson_n )
                stattext4 = gen_stat_string(rate_oloss_spearman_corr,   rate_oloss_spearman_p,   rate_oloss_spearman_n,   rate_oloss_kendalltau_corr,   rate_oloss_kendalltau_p,   rate_oloss_kendalltau_n  ,   rate_oloss_pearson_corr,   rate_oloss_pearson_p,   rate_oloss_pearson_n )
            for rowi in range(num_rows):
                ax[rowi][coli].set_xlabel(xlabel, fontsize=fontsize)
            if coli % 2 == 1:
                yfusiontop = np.log2(yfusiontop)
                yfusionbot = np.log2(yfusionbot)
                ylosstop   = np.log2(ylosstop)
                ylossbot   = np.log2(ylossbot)
            
            # Get age for coloring (negative values for proper timeline orientation)
            age_colors = -plotdf["age"]  # Negative so older = larger number
            
            # first, combined
            ax[0][coli].set_title(title1, fontsize=fontsize)
            ax[0][coli].scatter(x, yfusiontop, c=age_colors, cmap='viridis', vmin=age_colors.min(), vmax=age_colors.max())
            ax[0][coli].scatter(x, ylosstop, c=age_colors, cmap='viridis', vmin=age_colors.min(), vmax=age_colors.max())
            ax[0][coli].set_ylabel(ylabel1, fontsize=fontsize)
            # second, just fusions
            ax[1][coli].scatter(x, yfusiontop, c=age_colors, cmap='viridis', vmin=age_colors.min(), vmax=age_colors.max())
            ax[1][coli].set_title(title2,   fontsize=fontsize)
            ax[1][coli].set_ylabel(ylabel2, fontsize=fontsize)
            ax[1][coli].text(0.5, 0.5, stattext1, fontsize=fontsize, transform=ax[1][coli].transAxes, va = "center")
            # third, just dispersals
            ax[2][coli].scatter(x, ylosstop, c=age_colors, cmap='viridis', vmin=age_colors.min(), vmax=age_colors.max())
            ax[2][coli].set_title(title3,   fontsize=fontsize)
            ax[2][coli].set_ylabel(ylabel3, fontsize=fontsize)
            ax[2][coli].text(0.5, 0.5, stattext2, fontsize=fontsize, transform=ax[2][coli].transAxes, va = "center")

            # fourth, combined, but for the rate rather than for the ratio
            ax[3][coli].scatter(x, yfusionbot, c=age_colors, cmap='viridis', vmin=age_colors.min(), vmax=age_colors.max())
            ax[3][coli].scatter(x, ylossbot, c=age_colors, cmap='viridis', vmin=age_colors.min(), vmax=age_colors.max())
            ax[3][coli].set_title(title4,   fontsize=fontsize)
            ax[3][coli].set_ylabel(ylabel4, fontsize=fontsize)
            # second, just fusions
            ax[4][coli].scatter(x, yfusionbot, c=age_colors, cmap='viridis', vmin=age_colors.min(), vmax=age_colors.max())
            ax[4][coli].set_title(title5,   fontsize=fontsize)
            ax[4][coli].set_ylabel(ylabel5, fontsize=fontsize)
            ax[4][coli].text(0.5, 0.5, stattext3, fontsize=fontsize, transform=ax[4][coli].transAxes, va = "center")
            # third, just dispersals
            ax[5][coli].scatter(x, ylossbot, c=age_colors, cmap='viridis', vmin=age_colors.min(), vmax=age_colors.max())
            ax[5][coli].set_title(title6,   fontsize=fontsize)
            ax[5][coli].set_ylabel(ylabel6, fontsize=fontsize)
            ax[5][coli].text(0.5, 0.5, stattext4, fontsize=fontsize, transform=ax[5][coli].transAxes, va = "center")

        # change the fontsize of the axes and the titles
        for axi in range(num_rows):
            for axj in range(4):
                ax[axi][axj].tick_params(axis='both', which='major', labelsize=fontsize)
                ax[axi][axj].set_title(ax[axi][axj].get_title(), fontsize=fontsize)

        # force change the dot size
        dotsize = 1
        for axi in range(num_rows):
            for axj in [0,1,2,3]:
                for coll in ax[axi][axj].collections:
                    coll.set_sizes([dotsize])

        # change the opacity to 0.25
        opacity = 0.25
        for axi in range(num_rows):
            for axj in [0,1,2,3]:
                for coll in ax[axi][axj].collections:
                    coll.set_alpha(opacity)

        # Add colorbar for age mapping (use the last scatter plot collection)
        # Get age colors for colorbar
        age_colors_for_cbar = -plotdf["age"]
        # Add colorbar to the right side of the figure
        sm = plt.cm.ScalarMappable(cmap='viridis', 
                                   norm=plt.Normalize(vmin=age_colors_for_cbar.min(), 
                                                     vmax=age_colors_for_cbar.max()))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax.ravel().tolist(), location='right', pad=0.05, fraction=0.05)
        cbar.set_label('Age (Mya)', fontsize=10)
        # Invert the colorbar labels to show actual ages (not negative)
        cbar_ticks = cbar.get_ticks()
        cbar.set_ticklabels([f'{int(-t)}' for t in cbar_ticks])

        # increase the horizontal and vertical space between the panels
        plt.subplots_adjust(hspace=0.5, wspace=0.5)

        # save to a pdf with tight layout to prevent overlap
        outpdf =  outprefix.rstrip(".pdf").rstrip(".tsv") + ".pdf"
        plt.savefig(outpdf, facecolor='white', edgecolor='none', bbox_inches='tight')
        # close the figure to free up memory
        plt.close(fig)
    return statsresults

def get_edge_stats_single_taxid(taxid, edgedf, nodedf=None, phylo_weighter=None, use_time_slicing=True, num_threads=None, outdir=None, clade_tips_for_renorm=None, output_prefix=None) -> tuple:
    """
    This function takes the edgedf and builds up statistics about the edges.
    The things we will calculate:
      - The rate of fusions on this specific branch
      - The rate of fusions on this and child branches

    Args:
        taxid: The taxid to filter branches by
        edgedf: DataFrame with branch information
        nodedf: DataFrame with node information (needed for parallel processing)
        phylo_weighter: PhyloWeighting object for computing time-specific weights
        use_time_slicing: If True, compute weights per time slice; if False, use edgedf['phylo_weight']
        num_threads: Number of parallel threads for weight computation
        outdir: Output directory for diagnostic plots (if None, uses phylo_weighter cache dir)
        clade_tips_for_renorm: Set of tip taxids to filter and renormalize weights for clade analysis
        output_prefix: Prefix for output filenames (defaults to 'taxid{taxid}')

    Returns:
      - A tuple of (resultsdf_weighted, resultsdf_unweighted, (total_fusions, total_losses))
        where resultsdf_weighted uses phylogenetic weights and resultsdf_unweighted treats all branches equally
    """
    # Set output prefix
    if output_prefix is None:
        output_prefix = f'taxid{taxid}'
    
    # Use pre-processed lineage list (much faster than eval in loop)
    edgedf = edgedf[edgedf['child_lineage_list'].apply(lambda x: taxid in x)].copy()
    edgedf.reset_index(inplace=True, drop=True)

    # build a structured vector to hold the counts at each age.
    # get the unique source ages
    # make a datastructure with the counts of the source ages
    columns = {"total_branches_at_this_age":      0,
               "total_fusions_at_this_age":       0,
               "total_dispersals_at_this_age":        0,
               "fusion_rate_at_this_age_mean":   -1,
               "fusion_rate_at_this_age_median": -1,
               "fusion_rate_at_this_age_list":   [],
               "fusion_rate_at_this_age_weight_list": [],
               "fusion_rate_at_this_age_10th":   -1,
               "fusion_rate_at_this_age_90th":   -1,
               "dispersal_rate_at_this_age_mean":     -1,
               "dispersal_rate_at_this_age_median":   -1,
               "dispersal_rate_at_this_age_list":     [],
               "dispersal_rate_at_this_age_weight_list": [],
               "dispersal_rate_at_this_age_10th":     -1,
               "dispersal_rate_at_this_age_90th":     -1,
               }

    # get the oldest age in source_age in the edgedf
    oldest_age = int(edgedf['parent_age'].max())
    # populate the data structure with one row per integer
    count_struct_weighted = {i: copy.deepcopy(columns) for i in range(oldest_age + 1)}
    count_struct_unweighted = {i: copy.deepcopy(columns) for i in range(oldest_age + 1)}

    # Pre-compute time-specific weights using parallel processing or renormalization
    time_slice_weights = {}
    if use_time_slicing and phylo_weighter is not None:
        # Check if this is a clade analysis (renormalize from global) or root analysis (compute from scratch)
        if clade_tips_for_renorm is not None:
            # Clade analysis: load and renormalize global weights
            print(f"\n{'='*80}")
            print(f"LOADING AND RENORMALIZING GLOBAL WEIGHTS FOR CLADE")
            print(f"{'='*80}")
            print(f"Clade has {len(clade_tips_for_renorm)} tips")
            print(f"Loading global weights for {oldest_age + 1} time slices and renormalizing...")
            print(f"{'='*80}\n")
            
            import time
            start_time = time.time()
            
            # Set up diagnostics directory
            if outdir is not None:
                diagnostics_dir = outdir
            else:
                diagnostics_dir = os.path.join(phylo_weighter.cache_dir, 'global_diagnostics')
            os.makedirs(diagnostics_dir, exist_ok=True)
            
            # Load and renormalize weights for each time slice
            all_weights_for_diagnostic = []
            for time_slice in range(oldest_age + 1):
                renormed_weights = renormalize_weights_for_clade(
                    phylo_weighter, clade_tips_for_renorm, time_slice
                )
                
                if renormed_weights is not None:
                    time_slice_weights[time_slice] = renormed_weights
                    all_weights_for_diagnostic.extend(renormed_weights.values())
                else:
                    # Global weights not in cache - shouldn't happen if root was run first
                    print(f"WARNING: Global weights for time slice {time_slice} not found in cache!")
                    time_slice_weights[time_slice] = {}
                
                # Progress indicator every 100 slices
                if (time_slice + 1) % 100 == 0:
                    print(f"  Processed {time_slice + 1}/{oldest_age + 1} time slices")
            
            total_time = time.time() - start_time
            
            # Print summary statistics
            if all_weights_for_diagnostic:
                all_weights_arr = np.array(all_weights_for_diagnostic)
                print(f"\n{'='*80}")
                print(f"Weight renormalization complete in {total_time:.2f} seconds")
                print(f"Renormalized weight distribution across all time slices:")
                print(f"  Min:    {all_weights_arr.min():.3f}")
                print(f"  Mean:   {all_weights_arr.mean():.3f} (expect ~1.0 after renormalization)")
                print(f"  Median: {np.median(all_weights_arr):.3f}")
                print(f"  Max:    {all_weights_arr.max():.3f}")
                print(f"{'='*80}\n")
                
                # Generate diagnostic plot for time 0
                if 0 in time_slice_weights and time_slice_weights[0]:
                    print("Generating phylogenetic weight diagnostic plot for t=0...")
                    tips_at_t0 = set(time_slice_weights[0].keys())
                    tips_edgedf = edgedf[edgedf['child_taxid'].isin(tips_at_t0)].copy()
                    tips_edgedf['phylo_weight'] = tips_edgedf['child_taxid'].map(time_slice_weights[0])
                    plot_phylogenetic_weight_diagnostics(tips_edgedf, diagnostics_dir, output_prefix)
                    print("Diagnostic plot saved.\n")
                
                # Generate temporal heatmap
                plot_phylogenetic_weight_temporal_heatmap(time_slice_weights, diagnostics_dir, oldest_age, output_prefix)
        else:
            # Root analysis: compute weights from scratch (existing code)
            print(f"\n{'='*80}")
            print(f"PARALLELIZED TIME-SLICE PHYLOGENETIC WEIGHTING")
            print(f"{'='*80}")
            print(f"Computing weights for {oldest_age + 1} time slices")
            num_workers = num_threads if num_threads else min(cpu_count(), 16)
            print(f"Using {num_workers} parallel workers")
            print(f"Results will be cached for future runs.")
            print(f"{'='*80}\n")
        
            import time
            start_time = time.time()
            
            # Set up diagnostics directory
            if outdir is not None:
                diagnostics_dir = outdir
            else:
                diagnostics_dir = os.path.join(phylo_weighter.cache_dir, 'global_diagnostics')
            os.makedirs(diagnostics_dir, exist_ok=True)
            
            # Check which time slices are already cached
            time_slices = list(range(oldest_age + 1))
            print(f"Checking cache for {len(time_slices)} time slices...")
            cached_slices = []
            uncached_slices = []
            
            for i, time_slice in enumerate(time_slices):
                if (i + 1) % 100 == 0 or i == len(time_slices) - 1:
                    print(f"\rChecking cache: {i+1}/{len(time_slices)} ({100*(i+1)/len(time_slices):.1f}%)", end='', flush=True)
                
                _, weights_cache = phylo_weighter._get_cache_paths(time_slice)
                if weights_cache and os.path.exists(weights_cache):
                    # Load cached data
                    try:
                        with open(weights_cache, 'rb') as f:
                            cached_data = pickle.load(f)
                        time_slice_weights[time_slice] = cached_data
                        cached_slices.append(time_slice)
                    except:
                        uncached_slices.append(time_slice)
                else:
                    uncached_slices.append(time_slice)
            
            print(f"\n\nCache check complete:")
            print(f"  - Found {len(cached_slices)} cached time slices")
            print(f"  - Need to compute {len(uncached_slices)} time slices\n")
            
            # Collect all weights for statistics (whether cached or newly computed)
            all_weights_for_diagnostic = []
            for time_slice in time_slices:
                if time_slice in time_slice_weights and time_slice_weights[time_slice]:
                    all_weights_for_diagnostic.extend(time_slice_weights[time_slice].values())
            
            if len(uncached_slices) == 0:
                print("All time slices found in cache!\n")
                
                # Generate diagnostic plot for time 0 (when loaded from cache)
                if 0 in time_slice_weights and time_slice_weights[0]:
                    print("Generating phylogenetic weight diagnostic plot for t=0...")
                    tips_at_t0 = set(time_slice_weights[0].keys())
                    tips_edgedf = edgedf[edgedf['child_taxid'].isin(tips_at_t0)].copy()
                    tips_edgedf['phylo_weight'] = tips_edgedf['child_taxid'].map(time_slice_weights[0])
                    plot_phylogenetic_weight_diagnostics(tips_edgedf, diagnostics_dir, output_prefix)
                    print("Diagnostic plot saved.\n")
            else:
                # Prepare arguments for parallel processing (only uncached slices)
                phylo_args = (nodedf, phylo_weighter.cache_dir)
                args_list = [(t, edgedf, phylo_args) for t in uncached_slices]
                total_jobs = len(args_list)
                
                # Process in parallel with progress tracking
                print(f"Launching {total_jobs} parallel jobs across {num_workers} workers...\n")
                
                with Pool(num_workers) as pool:
                    # Use imap_unordered to get results as they complete
                    result_iterator = pool.imap_unordered(_compute_time_slice_weights_worker, args_list, chunksize=1)
                    
                    # Collect results with progress tracking
                    completed = 0
                    
                    for time_slice, weights_dict in result_iterator:
                        completed += 1
                        time_slice_weights[time_slice] = weights_dict
                        if weights_dict:
                            all_weights_for_diagnostic.extend(weights_dict.values())
                        
                        # Calculate progress statistics
                        elapsed = time.time() - start_time
                        percent = 100.0 * completed / total_jobs
                        running = min(num_workers, total_jobs - completed)
                        
                        if completed > 0:
                            rate = completed / elapsed  # jobs per second
                            remaining_jobs = total_jobs - completed
                            eta_seconds = remaining_jobs / rate if rate > 0 else 0
                            eta_minutes = eta_seconds / 60
                            
                            print(f"Completed: {completed}/{total_jobs} ({percent:.1f}%) | "
                                  f"Running: {running} | "
                                  f"Rate: {rate*60:.1f} jobs/min | "
                                  f"ETA: {eta_minutes:.1f} min")
                        
                        # Generate diagnostic plot for time 0
                        if time_slice == 0 and weights_dict:
                            print("\nGenerating phylogenetic weight diagnostic plot for t=0...")
                            
                            tips_at_t0 = set(weights_dict.keys())
                            tips_edgedf = edgedf[edgedf['child_taxid'].isin(tips_at_t0)].copy()
                            tips_edgedf['phylo_weight'] = tips_edgedf['child_taxid'].map(weights_dict)
                            plot_phylogenetic_weight_diagnostics(tips_edgedf, diagnostics_dir, output_prefix)
                            print("Diagnostic plot saved.")
                
                print()  # New line after progress complete
            
            total_time = time.time() - start_time
            # Print summary statistics
            if all_weights_for_diagnostic:
                all_weights_arr = np.array(all_weights_for_diagnostic)
                print(f"\n{'='*80}")
                if len(uncached_slices) > 0:
                    print(f"Time-slice weight computation complete in {total_time/60:.1f} minutes")
                else:
                    print(f"All weights loaded from cache in {total_time/60:.1f} minutes")
                print(f"Overall weight distribution across all time slices:")
                print(f"  Min:    {all_weights_arr.min():.3f}")
                print(f"  Mean:   {all_weights_arr.mean():.3f} (expect ~1.0 with no temporal bias)")
                print(f"  Median: {np.median(all_weights_arr):.3f}")
                print(f"  Max:    {all_weights_arr.max():.3f}")
                print(f"{'='*80}\n")
                
                # Generate temporal heatmap
                plot_phylogenetic_weight_temporal_heatmap(time_slice_weights, diagnostics_dir, oldest_age, output_prefix)
            else:
                if len(uncached_slices) > 0:
                    print(f"Time-slice weight computation complete in {total_time/60:.1f} minutes\n")
                else:
                    print(f"All weights loaded from cache in {total_time/60:.1f} minutes\n")

    total_fusions_in_this_branch = 0
    total_dispersals_in_this_branch = 0
    # Now go through and record the total number of branches at each age. This is done just using the edgedf.
    # We compute both weighted (using phylo_weight) and unweighted (all weights = 1) versions
    for i, row in edgedf.iterrows():
        # We just need to calculate this once
        num_fusions_my_this_branch = row["num_fusions_per_my_this_branch"]
        if np.isnan(num_fusions_my_this_branch) or np.isinf(num_fusions_my_this_branch):
            num_fusions_my_this_branch = 0
        num_dispersals_my_this_branch = row["num_dispersals_per_my_this_branch"]
        if np.isnan(num_dispersals_my_this_branch) or np.isinf(num_dispersals_my_this_branch):
            num_dispersals_my_this_branch = 0
        # now update the count_struct for each date
        # I need to iterate from the child age to the parent age, includsive.
        start = row["child_age"]
        end   = row["parent_age"]
        start_bin = int(row['child_age'])
        end_bin = int(row['parent_age'])
        child_taxid = row['child_taxid']
        
        # Use global phylo_weight if not using time-slicing
        if not use_time_slicing or phylo_weighter is None:
            phylo_weight = row['phylo_weight']
        
        # iterate through each time bin that this edge spans
        for bin_index in range(start_bin, end_bin + 1):
            # Get time-specific weight if using time-slicing
            if use_time_slicing and phylo_weighter is not None:
                phylo_weight = time_slice_weights[bin_index].get(child_taxid, 1.0)
            
            # Calculate the overlapping segment length within the bin
            overlap_start = max(start, bin_index)
            overlap_end   = min(end,   bin_index + 1)
            overlap_length = overlap_end - overlap_start

            # accumulate changes and lengths - WEIGHTED version
            if overlap_length > 0:
                count_struct_weighted[bin_index]['total_branches_at_this_age'] += overlap_length * phylo_weight
                count_struct_weighted[bin_index]['total_fusions_at_this_age']  += overlap_length * num_fusions_my_this_branch * phylo_weight
                count_struct_weighted[bin_index]['total_dispersals_at_this_age']   += overlap_length * num_dispersals_my_this_branch * phylo_weight
            elif overlap_length > 1:
                raise IOError("The overlap length is greater than 1. This should not happen. We are doing bins of 1 million years.")

            # accumulate changes and lengths - UNWEIGHTED version (weight = 1.0)
            if overlap_length > 0:
                count_struct_unweighted[bin_index]['total_branches_at_this_age'] += overlap_length
                count_struct_unweighted[bin_index]['total_fusions_at_this_age']  += overlap_length * num_fusions_my_this_branch
                count_struct_unweighted[bin_index]['total_dispersals_at_this_age']   += overlap_length * num_dispersals_my_this_branch

            # Append to lists for quantile calculations
            count_struct_weighted[bin_index]['fusion_rate_at_this_age_list'].append(num_fusions_my_this_branch)
            count_struct_weighted[bin_index]['fusion_rate_at_this_age_weight_list'].append(phylo_weight)
            count_struct_weighted[bin_index]['dispersal_rate_at_this_age_list'].append(num_dispersals_my_this_branch)
            count_struct_weighted[bin_index]['dispersal_rate_at_this_age_weight_list'].append(phylo_weight)
            
            count_struct_unweighted[bin_index]['fusion_rate_at_this_age_list'].append(num_fusions_my_this_branch)
            count_struct_unweighted[bin_index]['fusion_rate_at_this_age_weight_list'].append(1.0)
            count_struct_unweighted[bin_index]['dispersal_rate_at_this_age_list'].append(num_dispersals_my_this_branch)
            count_struct_unweighted[bin_index]['dispersal_rate_at_this_age_weight_list'].append(1.0)

    # now calculate the quantiles for each position - WEIGHTED version
    for i in range(oldest_age + 1):
        # Fusion rate statistics with phylogenetic weighting
        fusion_values = count_struct_weighted[i]['fusion_rate_at_this_age_list']
        fusion_weights = count_struct_weighted[i]['fusion_rate_at_this_age_weight_list']
        
        if len(fusion_values) > 0:
            count_struct_weighted[i]['fusion_rate_at_this_age_mean'] = np.average(fusion_values, weights=fusion_weights)
            count_struct_weighted[i]['fusion_rate_at_this_age_median'] = weighted_quantile(fusion_values, fusion_weights, 0.5)
        else:
            count_struct_weighted[i]['fusion_rate_at_this_age_mean'] = 0.0
            count_struct_weighted[i]['fusion_rate_at_this_age_median'] = 0.0
            
        count_struct_weighted[i]['fusion_rate_at_this_age_10th']   = weighted_quantile(
            fusion_values, fusion_weights, 0.1
        )
        count_struct_weighted[i]['fusion_rate_at_this_age_90th']   = weighted_quantile(
            fusion_values, fusion_weights, 0.9
        )
        
        # Dispersal rate statistics with phylogenetic weighting
        dispersal_values = count_struct_weighted[i]['dispersal_rate_at_this_age_list']
        dispersal_weights = count_struct_weighted[i]['dispersal_rate_at_this_age_weight_list']
        
        if len(dispersal_values) > 0:
            count_struct_weighted[i]['dispersal_rate_at_this_age_mean'] = np.average(dispersal_values, weights=dispersal_weights)
            count_struct_weighted[i]['dispersal_rate_at_this_age_median'] = weighted_quantile(dispersal_values, dispersal_weights, 0.5)
        else:
            count_struct_weighted[i]['dispersal_rate_at_this_age_mean'] = 0.0
            count_struct_weighted[i]['dispersal_rate_at_this_age_median'] = 0.0
            
        count_struct_weighted[i]['dispersal_rate_at_this_age_10th']     = weighted_quantile(
            dispersal_values, dispersal_weights, 0.1
        )
        count_struct_weighted[i]['dispersal_rate_at_this_age_90th']     = weighted_quantile(
            dispersal_values, dispersal_weights, 0.9
        )

    # now calculate the quantiles for each position - UNWEIGHTED version
    for i in range(oldest_age + 1):
        # Fusion rate statistics (unweighted)
        fusion_values = count_struct_unweighted[i]['fusion_rate_at_this_age_list']
        fusion_weights = count_struct_unweighted[i]['fusion_rate_at_this_age_weight_list']
        
        if len(fusion_values) > 0:
            count_struct_unweighted[i]['fusion_rate_at_this_age_mean'] = np.mean(fusion_values)
            count_struct_unweighted[i]['fusion_rate_at_this_age_median'] = np.median(fusion_values)
        else:
            count_struct_unweighted[i]['fusion_rate_at_this_age_mean'] = 0.0
            count_struct_unweighted[i]['fusion_rate_at_this_age_median'] = 0.0
            
        count_struct_unweighted[i]['fusion_rate_at_this_age_10th']   = weighted_quantile(
            fusion_values, fusion_weights, 0.1
        )
        count_struct_unweighted[i]['fusion_rate_at_this_age_90th']   = weighted_quantile(
            fusion_values, fusion_weights, 0.9
        )
        
        # Dispersal rate statistics (unweighted)
        dispersal_values = count_struct_unweighted[i]['dispersal_rate_at_this_age_list']
        dispersal_weights = count_struct_unweighted[i]['dispersal_rate_at_this_age_weight_list']
        
        if len(dispersal_values) > 0:
            count_struct_unweighted[i]['dispersal_rate_at_this_age_mean'] = np.mean(dispersal_values)
            count_struct_unweighted[i]['dispersal_rate_at_this_age_median'] = np.median(dispersal_values)
        else:
            count_struct_unweighted[i]['dispersal_rate_at_this_age_mean'] = 0.0
            count_struct_unweighted[i]['dispersal_rate_at_this_age_median'] = 0.0
            
        count_struct_unweighted[i]['dispersal_rate_at_this_age_10th']     = weighted_quantile(
            dispersal_values, dispersal_weights, 0.1
        )
        count_struct_unweighted[i]['dispersal_rate_at_this_age_90th']     = weighted_quantile(
            dispersal_values, dispersal_weights, 0.9
        )

    # all of the data should be in the struct now. Turn it into pandas dfs
    resultsdf_weighted = pd.DataFrame(count_struct_weighted).T
    resultsdf_unweighted = pd.DataFrame(count_struct_unweighted).T
    
    # remove the rows where total_branches_at_this_age is 0
    resultsdf_weighted = resultsdf_weighted[resultsdf_weighted['total_branches_at_this_age'] > 0]
    resultsdf_unweighted = resultsdf_unweighted[resultsdf_unweighted['total_branches_at_this_age'] > 0]

    for event in ["fusions", "dispersals"]:
        resultsdf_weighted[f"{event}_ratio"] = resultsdf_weighted[f"total_{event}_at_this_age"] / resultsdf_weighted["total_branches_at_this_age"]
        resultsdf_unweighted[f"{event}_ratio"] = resultsdf_unweighted[f"total_{event}_at_this_age"] / resultsdf_unweighted["total_branches_at_this_age"]
    
    resultsdf_weighted["age"] = -1 * resultsdf_weighted.index
    resultsdf_unweighted["age"] = -1 * resultsdf_unweighted.index
    
    return resultsdf_weighted, resultsdf_unweighted, (total_fusions_in_this_branch, total_dispersals_in_this_branch)

def add_events_to_edge_df(edgedf, eventdf) -> pd.DataFrame:
    """
    All that this does is tracks the types and quantity of changes on each edge.
    This information will later be used in the other functions for plotting.
    """
    # we need a taxid_to_parent dictionary to get the mapping to id the branch length
    # Actually, we can just get the index of this row, because we now know that the target column is unique,
    #   and we will want to work with multiple values from this row.
    taxid_to_parent = {row['child_taxid']: row['parent_taxid'] for i, row in edgedf.iterrows()}

    # Initialize a dictionary to hold the new columns
    new_columns = {"num_fusions_this_branch": 0,
                   "num_dispersals_this_branch": 0,
                   "num_fusions_per_my_this_branch": 0,
                   "num_dispersals_per_my_this_branch": 0}
    # get all of the unique categories without the parentheses from the eventdf
    ALGs         = sorted(set([x for x in eventdf['change'].unique() if x[0] != '(']))
    for thisALG in ALGs:
        new_columns["num_" + thisALG + "_this_branch"] = 0

    # Get all of the possible combinations
    Combinations = sorted(set([x for x in eventdf['change'].unique() if x[0] == '(']))
    for thisCombo in Combinations:
        combostr = "+".join(eval(thisCombo))
        new_columns["num_" + combostr + "_this_branch"] = 0
    
    # Drop columns that already exist to avoid duplicates
    # Also explicitly drop old "losses" columns that are now "dispersals"
    cols_to_drop = [col for col in new_columns.keys() if col in edgedf.columns]
    legacy_cols = ['num_losses_this_branch', 'num_losses_per_my_this_branch']
    cols_to_drop.extend([col for col in legacy_cols if col in edgedf.columns])
    if cols_to_drop:
        edgedf = edgedf.drop(columns=cols_to_drop)
        print(f"Dropped {len(cols_to_drop)} existing columns to repopulate them fresh: {cols_to_drop}")
    
    new_columns_df = pd.DataFrame(new_columns, index=edgedf.index)
    # Concatenate the new columns to the original DataFrame
    edgedf = pd.concat([edgedf, new_columns_df], axis=1)

    missing_events = 0
    for i, row in eventdf.iterrows():
        # This is a single event. We need to know the branch on which it occurred.
        event_branch_id = row["target_taxid"]
        if event_branch_id not in taxid_to_parent:
            missing_events += 1
            continue
        else:
            parent_taxid = taxid_to_parent[event_branch_id]
            edge = (parent_taxid, event_branch_id)

            # get the indices of the edgedf that this matches,
            edge_row = edgedf[(edgedf['parent_taxid'] == parent_taxid) & (edgedf['child_taxid'] == event_branch_id)]
            if len(edge_row) > 1:
                raise ValueError('For some reason we found this edge more than once. There is probably an error with the output software. {} {}'.format(
                    parent_taxid, event_branch_id))
            elif len(edge_row) == 0:
                missing_events += 1
                continue
            else:
                # determine whether the event was a fusion or a dispersal
                colnames = []
                changetype = None
                if row["change"][0] == "(":
                    changetype = "fusions"
                    colnames.append("num_fusions_this_branch")
                    combostr = "+".join(eval(row["change"]))
                    colnames.append("num_" + combostr + "_this_branch")
                else:
                    changetype = "dispersals"
                    colnames.append("num_dispersals_this_branch")
                    colnames.append("num_" + row["change"] + "_this_branch")
                for thiscol in colnames:
                    edgedf.loc[edge_row.index, thiscol] += 1
    print("There were {} missing events".format(missing_events))

    edgedf["num_fusions_per_my_this_branch"] = edgedf["num_fusions_this_branch"] / edgedf["branch_length"]
    edgedf["num_dispersals_per_my_this_branch"] = edgedf["num_dispersals_this_branch"] / edgedf["branch_length"]
    return edgedf

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    print(f"Output directory: {args.outdir}\n")
    
    # Set up logging to capture warnings in a file
    logging.basicConfig(
        filename=os.path.join(args.outdir, 'warnings.log'),
        level=logging.WARNING,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.captureWarnings(True)
    
    # Suppress console output of specific warnings
    warnings.filterwarnings('ignore', category=RuntimeWarning, message='divide by zero encountered in log2')
    warnings.filterwarnings('ignore', message='An input array is constant')

    # read in the dataframes
    edgedf  = pd.read_csv(args.edge_information, sep='\t', index_col=False)
    nodedf  = pd.read_csv(args.node_information, sep='\t', index_col=False)
    # enforce that node parent is int. First, convert NaNs to type None
    nodedf['parent'] = nodedf['parent'].apply(lambda x: -1 if pd.isna(x) else int(x)).astype(int)
    eventdf = pd.read_csv(args.statsdf, sep='\t', index_col=False)
    
    # OPTIMIZATION: Pre-convert string lineages to lists once (huge speedup)
    print("Pre-processing lineage columns...")
    edgedf['child_lineage_list'] = edgedf['child_lineage'].apply(eval)
    print("Done pre-processing.\n")

    print("This is nodedf")
    print(nodedf)
    print("This is edgedf")
    print(edgedf)
    print("\n=== DEBUG: edgedf structure ===")
    print(f"Index: {edgedf.index}")
    print(f"Index names: {edgedf.index.names}")
    print(f"Columns: {list(edgedf.columns)}")
    print(f"Shape: {edgedf.shape}")
    print("\nFirst few rows with column alignment:")
    print(edgedf.head(10).to_string())
    print("=== END DEBUG ===\n")

    edgedf = add_events_to_edge_df(edgedf, eventdf)
    print("Edgedf has been marked up")
    # save this to a new file, as it has the number of changes per branch
    edgedf.to_csv(os.path.join(args.outdir, "modified_edge_list.tsv"), sep='\t', index=False)

    # print the columns that have _this_branch

    # verify that there every value in the target column is unique
    if len(edgedf['child_taxid'].unique()) != len(edgedf):
        raise ValueError('The target column (child_taxid) is not unique')
    else:
        print("Every value in the target (child_taxid) column is unique. This is a legal DAG.")

    # Initialize phylogenetic weighter for time-slice weighting
    print("Initializing phylogenetic weighting system...")
    phylo_weighter = PhyloWeighting(edgedf, nodedf, cache_dir=args.outdir)
    print("Phylogenetic weighting system ready.\n")
    
    # Create subdirectories for output
    per_clade_dir = os.path.join(args.outdir, 'per_clade_analyses')
    os.makedirs(per_clade_dir, exist_ok=True)
    global_diagnostics_dir = os.path.join(args.outdir, 'global_diagnostics')
    os.makedirs(global_diagnostics_dir, exist_ok=True)

    #for node in nodedf['taxid']:
    nodedf["fusions_in_this_clade"] = -1
    nodedf["losses_in_this_clade"] = -1
    # change the types of the previous two to int
    nodedf["fusions_in_this_clade"] = nodedf["fusions_in_this_clade"].astype(int)
    nodedf["losses_in_this_clade"] = nodedf["losses_in_this_clade"].astype(int)
    nodedf["fusions_in_this_clade_div_dist_crown"] = -1.00000000001
    nodedf["losses_in_this_clade_div_dist_crown"]  = -1.00000000001
    nodedf["fusions_in_this_clade_div_dist_crown_plus_root"] = -1.00000000001
    nodedf["losses_in_this_clade_div_dist_crown_plus_root"]  = -1.00000000001

    # ========================================================================
    # Helper function for building descriptive names for custom analyses
    # ========================================================================
    def build_custom_clade_names(target_taxid, exclude_taxids, nodedf):
        """
        Build descriptive directory and output prefix names for custom clade analysis.
        
        Returns:
            tuple: (directory_name, output_prefix)
                   e.g., ("Vertebrata_7742_minus_Teleostei_32443", "Vertebrata_minus_Teleostei")
        """
        # Get target clade name
        target_row = nodedf[nodedf['taxid'] == target_taxid]
        if len(target_row) > 0 and pd.notna(target_row.iloc[0]['name']) and target_row.iloc[0]['name'] != "":
            target_name = target_row.iloc[0]['name']
        else:
            target_name = str(target_taxid)
        
        # Build directory name with taxids
        dir_name = f"{target_name}_{target_taxid}"
        
        # Build output prefix with names only
        prefix = target_name
        
        if exclude_taxids:
            dir_name += "_minus"
            prefix += "_minus"
            
            for excl_taxid in exclude_taxids:
                excl_row = nodedf[nodedf['taxid'] == excl_taxid]
                if len(excl_row) > 0 and pd.notna(excl_row.iloc[0]['name']) and excl_row.iloc[0]['name'] != "":
                    excl_name = excl_row.iloc[0]['name']
                else:
                    excl_name = str(excl_taxid)
                
                dir_name += f"_{excl_name}_{excl_taxid}"
                prefix += f"_{excl_name}"
        
        return dir_name, prefix

    # Mapping dict will have strings as keys that will become column names.
    #  The values will be another dict in which the keys are the ages(negative integers)
    #  and the values are the counts. These dicts will then be iterated through to, so we
    #  can add multiple columns to the dataframe in case there are multiple things we want
    #  to measure.
    mapping_dict = {}
    # intensity of extinction
    if args.intensity_of_extinction is not None:
        IOEdf = pd.read_csv(args.intensity_of_extinction, sep='\t')
        # First do Rohde and Muller Intensities
        TIME        = "Time (Ma)"
        IOEdf[TIME] = IOEdf[TIME].astype(int) * -1
        EXTINCTION  = "Extinction Intensity (%)"
        ORIGINATION = "Origination Intensity (%)"
        # zip together time as the key and the intensity as the value
        mapping_dict["rohde_extinction"]  = dict(zip(IOEdf[TIME], IOEdf[EXTINCTION]))
        mapping_dict["rohde_origination"] = dict(zip(IOEdf[TIME], IOEdf[ORIGINATION]))

    # Process each node in the tree
    ## just do Metazoa now
    #for node in [33208]:
    
    # Define major clades for weight distribution analysis
    target_clades = ['Lepidoptera', 'Cnidaria', 'Hexapoda', 'Spiralia', 'Vertebrata', 
                     'Bilateria', 'Panarthropoda', 'Arthropoda', 'Brachiopoda', 'Bryozoa', 
                     'Chaetognatha', 'Chordata', 'Ctenophora', 'Cycliophora', 'Echinodermata', 
                     'Entoprocta', 'Gastrotrich', 'Gnathostomulida', 'Hemichordata', 'Kinorhyncha', 
                     'Loricifera', 'Mollusca', 'Nematoda', 'Nematomorpha', 'Nemertea', 'Onychophora', 
                     'Orthonectida', 'Phoronida', 'Placozoa', 'Platyhelminthes', 'Porifera', 
                     'Priapulida', 'Rotifera', 'Tardigrada', 'Xenacoelomorpha', 'Protostomia', 
                     'Deuterostomia']
    
    # ========================================================================
    # SPECIAL MODE: Single clade analysis with subclade exclusions
    # ========================================================================
    if args.analyze_single_clade:
        print("\n" + "="*80)
        print("SPECIAL MODE: Single clade analysis with subclade exclusions")
        print("="*80 + "\n")
        
        target_taxid = args.analyze_single_clade
        exclude_taxids = args.exclude_subclades
        
        # Build descriptive names
        dir_name, output_prefix = build_custom_clade_names(target_taxid, exclude_taxids, nodedf)
        custom_dir = os.path.join(args.outdir, "custom_clade_analyses", dir_name)
        os.makedirs(custom_dir, exist_ok=True)
        
        print(f"Target clade: {target_taxid}")
        print(f"Excluding subclades: {exclude_taxids}")
        print(f"Output directory: {custom_dir}")
        print(f"Output prefix: {output_prefix}\n")
        
        # Validate exclusion relationships
        print("Validating exclusion relationships...")
        for excl_taxid in exclude_taxids:
            total_excl_edges = len(edgedf[edgedf['child_lineage_list'].apply(lambda x: excl_taxid in x)])
            excl_within_target = len(edgedf[edgedf['child_lineage_list'].apply(
                lambda x: excl_taxid in x and target_taxid in x)])
            
            if total_excl_edges > 0:
                percentage = 100.0 * excl_within_target / total_excl_edges
                print(f"  Taxid {excl_taxid}: {percentage:.1f}% ({excl_within_target}/{total_excl_edges}) of its branches are within target taxid {target_taxid}")
            else:
                print(f"  WARNING: Taxid {excl_taxid} has no branches in the tree")
        print()
        
        # Filter edges: include target, exclude subclades
        print("Filtering edges and nodes...")
        edgedf_filtered = edgedf[edgedf['child_lineage_list'].apply(
            lambda x: target_taxid in x and not any(excl in x for excl in exclude_taxids))].copy()
        
        print(f"  Original edges in target clade: {len(edgedf[edgedf['child_lineage_list'].apply(lambda x: target_taxid in x)])}")
        print(f"  Filtered edges (after exclusions): {len(edgedf_filtered)}")
        
        if len(edgedf_filtered) == 0:
            print("\nERROR: No edges remaining after filtering. The exclusions removed the entire clade.")
            sys.exit(1)
        
        # Filter nodes to match filtered edges
        nodes_in_filtered_edges = set(edgedf_filtered['parent_taxid'].tolist() + edgedf_filtered['child_taxid'].tolist())
        nodedf_filtered = nodedf[nodedf['taxid'].isin(nodes_in_filtered_edges)].copy()
        print(f"  Filtered nodes: {len(nodedf_filtered)}\n")
        
        # Compute clade tips
        clade_tips = set()
        for _, row in edgedf_filtered.iterrows():
            child = row['child_taxid']
            if child not in edgedf_filtered['parent_taxid'].values:
                clade_tips.add(child)
        
        num_tips = len(clade_tips)
        print(f"Clade has {num_tips} species after filtering")
        
        if num_tips < 50:
            print(f"  WARNING: Clade has fewer than 50 species (minimum recommended: 50)")
            print(f"  Proceeding anyway since this is an explicitly requested analysis\n")
        
        # Create a dedicated PhyloWeighting object for this filtered tree
        # This ensures weights are computed based on the filtered topology
        # and cached separately in the custom analysis directory
        print("Creating PhyloWeighting object for filtered tree...")
        custom_phylo_weighter = PhyloWeighting(edgedf_filtered, nodedf_filtered, 
                                              cache_dir=custom_dir, verbose=True)
        print()
        
        # Always compute from scratch for filtered trees (don't use global cache)
        # The filtered tree has a different topology, so weights must be recomputed
        print("Computing phylogenetic weights from scratch for filtered tree")
        print(f"(Cache will be stored in {custom_dir}/phylo_cache/ for future runs)\n")
        clade_tips_for_renorm = None
        
        # Run analysis
        print("Running branch statistics analysis...")
        resultsdf_weighted, resultsdf_unweighted, (fusions, dispersals) = get_edge_stats_single_taxid(
            target_taxid, edgedf_filtered, nodedf=nodedf_filtered, 
            phylo_weighter=custom_phylo_weighter, use_time_slicing=True,
            num_threads=args.threads, outdir=custom_dir, 
            clade_tips_for_renorm=clade_tips_for_renorm,
            output_prefix=output_prefix
        )
        
        # Generate verification plots
        print(f"\nGenerating verification plots for {output_prefix}...")
        plot_weighted_vs_unweighted_comparison(resultsdf_weighted, resultsdf_unweighted, 
                                               custom_dir, output_prefix=output_prefix)
        plot_event_count_conservation(resultsdf_weighted, resultsdf_unweighted, 
                                      edgedf_filtered, custom_dir, output_prefix=output_prefix)
        
        # Add extinction intensity mapping if provided (to both dataframes)
        if args.intensity_of_extinction is not None:
            for thiscol in mapping_dict:
                resultsdf_weighted[thiscol] = resultsdf_weighted["age"].map(mapping_dict[thiscol])
                resultsdf_unweighted[thiscol] = resultsdf_unweighted["age"].map(mapping_dict[thiscol])
        
        # Generate temporal weight analysis plot (if phylo_weight exists)
        if 'phylo_weight' in edgedf_filtered.columns:
            print(f"  Generating temporal weight analysis for {output_prefix}...")
            plot_phylogenetic_weights_temporal(edgedf_filtered, custom_dir, output_prefix)
        else:
            print(f"  Skipping temporal weight analysis for {output_prefix} (not compatible with time-slicing)")
        
        # Save both weighted and unweighted results
        print(f"\nSaving results to TSV files...")
        outprefix_weighted = os.path.join(custom_dir, f"{output_prefix}_changes_vs_age_weighted")
        resultsdf_weighted.to_csv(outprefix_weighted + ".tsv", sep='\t', index=False)
        outprefix_unweighted = os.path.join(custom_dir, f"{output_prefix}_changes_vs_age_unweighted")
        resultsdf_unweighted.to_csv(outprefix_unweighted + ".tsv", sep='\t', index=False)
        print(f"  Saved: {outprefix_weighted}.tsv")
        print(f"  Saved: {outprefix_unweighted}.tsv")
        
        # Plot changes vs time
        outprefix_time = os.path.join(custom_dir, f"{output_prefix}_changes_vs_time")
        if not args.suppress_plotting:
            print(f"\n  Plotting changes vs time...")
            if args.intensity_of_extinction is not None:
                plot_fusions_per_branch_vs_time(outprefix_time, resultsdf_weighted, resultsdf_unweighted, 
                                               intensity_of_extinction_filepath=args.intensity_of_extinction)
            else:
                plot_fusions_per_branch_vs_time(outprefix_time, resultsdf_weighted, resultsdf_unweighted)
            print(f"  Saved: {outprefix_time}.pdf")
        
        # Plot changes vs intensity of extinction
        if args.intensity_of_extinction is not None:
            print(f"\n  Plotting changes vs extinction intensity...")
            outprefix_intensity = os.path.join(custom_dir, f"{output_prefix}_changes_vs_intensity")
            plot_intensity_of_extinction(outprefix_intensity, resultsdf_weighted, 
                                        args.intensity_of_extinction, 
                                        suppress_plotting=args.suppress_plotting)
            if not args.suppress_plotting:
                print(f"  Saved: {outprefix_intensity}.pdf")
        
        print("\n" + "="*80)
        print("SPECIAL MODE: Analysis complete")
        print(f"Results saved to: {custom_dir}")
        print("="*80 + "\n")
        
        # Exit early - don't run normal per-clade analysis
        return
    
    # ========================================================================
    # PHASE 1: Fast calculation of event totals for ALL clades
    # ========================================================================
    print("\n" + "="*80)
    print("PHASE 1: Calculating unweighted event totals for all clades")
    print("="*80 + "\n")
    
    # Initialize columns (unweighted only)
    nodedf['fusions_in_this_clade_unweighted'] = 0.0
    nodedf['dispersals_in_this_clade_unweighted'] = 0.0
    
    for i, row in nodedf.iterrows():
        node = row["taxid"]
        clade_name = row["name"]
        if pd.isna(clade_name) or clade_name == "":
            clade_name = f"Node{node}"
        
        # Print progress every 100 clades, plus first and last
        if i == 0 or (i+1) % 100 == 0 or i == len(nodedf) - 1:
            print(f"  Processing clade {i+1}/{len(nodedf)}: {clade_name} (taxid={node})")
        
        # Get clade tips
        edgedf_clade = edgedf[edgedf['child_lineage_list'].apply(lambda x: node in x)].copy()
        clade_tips = set()
        for _, edge_row in edgedf_clade.iterrows():
            child = edge_row['child_taxid']
            if child not in edgedf_clade['parent_taxid'].values:
                clade_tips.add(child)
        
        # Calculate unweighted totals
        unweighted_fusions, unweighted_dispersals = calculate_clade_event_totals(node, edgedf)
        
        # Store in nodedf
        nodedf.loc[i, 'fusions_in_this_clade_unweighted'] = unweighted_fusions
        nodedf.loc[i, 'dispersals_in_this_clade_unweighted'] = unweighted_dispersals
    
    print("\n" + "="*80)
    print("Phase 1 complete! Calculating rates...")
    print("="*80 + "\n")
    
    # Calculate rates using unweighted values
    nodedf["fusions_in_this_clade_div_dist_crown"] = nodedf["fusions_in_this_clade_unweighted"] / nodedf["dist_crown"]
    nodedf["dispersals_in_this_clade_div_dist_crown"] = nodedf["dispersals_in_this_clade_unweighted"] / nodedf["dist_crown"]
    nodedf["fusions_in_this_clade_div_dist_crown_plus_root"] = nodedf["fusions_in_this_clade_unweighted"] / nodedf["dist_crown_plus_root"]
    nodedf["dispersals_in_this_clade_div_dist_crown_plus_root"] = nodedf["dispersals_in_this_clade_unweighted"] / nodedf["dist_crown_plus_root"]
    
    # Save the summary file
    summary_file = os.path.join(args.outdir, "modified_node_list.tsv")
    nodedf.to_csv(summary_file, sep='\t', index=False)
    print(f"Saved summary file: {summary_file}\n")
    
    # ========================================================================
    # PHASE 2: Detailed analysis and plotting for selected clades
    # ========================================================================
    print("\n" + "="*80)
    print("PHASE 2: Generating detailed plots for clades with ≥30 species")
    print("="*80 + "\n")
    
    verification_plots_generated = False
    
    for i, row in nodedf.iterrows():
        node = row["taxid"]
        # Get clade name from the node information file (already handles custom topology nodes)
        clade_name = row["name"]
        if pd.isna(clade_name) or clade_name == "":
            clade_name = f"Node{node}"  # Fallback name if missing
        
        print("  - We are processing node {} / {}".format(i, len(nodedf)))
        print(f"    Clade: {clade_name} (taxid={node})")
        
        # Create clade-specific PhyloWeighting object for recentered weights
        # Filter edgedf to this clade first
        edgedf_clade = edgedf[edgedf['child_lineage_list'].apply(lambda x: node in x)].copy()
        print(f"    Clade has {len(edgedf_clade)} branches")
        
        # Count number of tips (species) in this clade
        # A tip is a child that's not a parent in this clade
        clade_tips = set()
        for _, row in edgedf_clade.iterrows():
            child = row['child_taxid']
            if child not in edgedf_clade['parent_taxid'].values:
                clade_tips.add(child)
        
        num_tips = len(clade_tips)
        
        # Skip entire analysis for clades with <30 species (except root)
        required_num_tips = 50
        if i > 0 and num_tips < required_num_tips:
            print(f"    Clade has only {num_tips} species - skipping analysis (minimum: {required_num_tips})")
            continue
        
        print(f"    Clade has {num_tips} species - generating full analysis")
        
        # For the first node (root), compute weights from scratch for global analysis
        # For all nodes (including root), also run per-clade analysis with renormalized weights
        if i == 0:
            # Root node: compute from scratch and save to global cache
            print(f"    Computing phylogenetic weights from scratch (root analysis)")
            resultsdf_weighted, resultsdf_unweighted, (fusions, dispersals) = get_edge_stats_single_taxid(
                node, edgedf, nodedf=nodedf, phylo_weighter=phylo_weighter, use_time_slicing=True, 
                num_threads=args.threads, outdir=global_diagnostics_dir, clade_tips_for_renorm=None,
                output_prefix='ALL'
            )
            # Also run per-clade analysis for root with filtered tips
            print(f"    Running per-clade analysis for {clade_name} (renormalized for clade tips)")
            resultsdf_weighted, resultsdf_unweighted, (fusions, dispersals) = get_edge_stats_single_taxid(
                node, edgedf, nodedf=nodedf, phylo_weighter=phylo_weighter, use_time_slicing=True, 
                num_threads=args.threads, outdir=per_clade_dir, clade_tips_for_renorm=clade_tips,
                output_prefix=f'{clade_name}_{node}'
            )
        else:
            # Per-clade: renormalize from global weights
            print(f"    Reusing global phylogenetic weights, renormalized for this clade")
            resultsdf_weighted, resultsdf_unweighted, (fusions, dispersals) = get_edge_stats_single_taxid(
                node, edgedf, nodedf=nodedf, phylo_weighter=phylo_weighter, use_time_slicing=True, 
                num_threads=args.threads, outdir=per_clade_dir, clade_tips_for_renorm=clade_tips,
                output_prefix=f'{clade_name}_{node}'
            )
        
        # Generate global diagnostic plots only for the first (root) node
        if not verification_plots_generated:
            print("\n  Generating global diagnostic plots...")
            
            # Save global verification plots to global_diagnostics_dir
            print("    - Weighted vs unweighted comparison plot")
            plot_weighted_vs_unweighted_comparison(resultsdf_weighted, resultsdf_unweighted, global_diagnostics_dir)
            
            print("    - Event count conservation plot")
            plot_event_count_conservation(resultsdf_weighted, resultsdf_unweighted, edgedf, global_diagnostics_dir)
            
            # Clade weight distribution plot - use time-slicing compatible version
            print("    - Clade-specific weight distribution plot")
            plot_clade_weight_distribution_from_cache(phylo_weighter, edgedf, global_diagnostics_dir, target_clades)
            
            verification_plots_generated = True
            print("  Global diagnostic plots complete.\n")
        
        # Generate verification plots for each clade (including root)
        print(f"\n  Generating verification plots for {clade_name}...")
        
        # Weighted vs unweighted comparison
        plot_weighted_vs_unweighted_comparison(resultsdf_weighted, resultsdf_unweighted, per_clade_dir, output_prefix=f'{clade_name}_{node}')
        
        # Event count conservation
        plot_event_count_conservation(resultsdf_weighted, resultsdf_unweighted, edgedf_clade, per_clade_dir, output_prefix=f'{clade_name}_{node}')
        
        # Filter edgedf to this clade for temporal weight analysis
        edgedf_subset = edgedf[edgedf['child_lineage_list'].apply(lambda x: node in x)].copy()
        
        # Add extinction intensity mapping if provided (to both dataframes)
        if args.intensity_of_extinction is not None:
            for thiscol in mapping_dict:
                resultsdf_weighted[thiscol] = resultsdf_weighted["age"].map(mapping_dict[thiscol])
                resultsdf_unweighted[thiscol] = resultsdf_unweighted["age"].map(mapping_dict[thiscol])
        
        # Generate temporal weight analysis plot for this subset (if phylo_weight exists)
        if 'phylo_weight' in edgedf_subset.columns:
            print(f"  Generating temporal weight analysis for {clade_name}_{node}...")
            plot_phylogenetic_weights_temporal(edgedf_subset, per_clade_dir, f"{clade_name}_{node}")
        else:
            print(f"  Skipping temporal weight analysis for {clade_name}_{node} (not compatible with time-slicing)")
        
        # Save both weighted and unweighted results
        outprefix0_weighted = os.path.join(per_clade_dir, f"{clade_name}_{node}_changes_vs_age_weighted")
        resultsdf_weighted.to_csv(outprefix0_weighted + ".tsv", sep='\t', index=False)
        outprefix0_unweighted = os.path.join(per_clade_dir, f"{clade_name}_{node}_changes_vs_age_unweighted")
        resultsdf_unweighted.to_csv(outprefix0_unweighted + ".tsv", sep='\t', index=False)

        outprefix1 = os.path.join(per_clade_dir, f"{clade_name}_{node}_changes_vs_time")
        if not args.suppress_plotting:
            # plot the changes per time - now passing both weighted and unweighted dataframes
            if args.intensity_of_extinction is not None:
                plot_fusions_per_branch_vs_time(outprefix1, resultsdf_weighted, resultsdf_unweighted, 
                                               intensity_of_extinction_filepath = args.intensity_of_extinction)
            else:
                plot_fusions_per_branch_vs_time(outprefix1, resultsdf_weighted, resultsdf_unweighted)

        outprefix2 = os.path.join(per_clade_dir, f"{clade_name}_{node}_changes_vs_intensity")
        # now plot the intensity of extinction with the change types depending on whether they are present.
        # Use weighted results for intensity correlation analysis
        if args.intensity_of_extinction is not None and i > 0:
            statsdf = plot_intensity_of_extinction(outprefix2, resultsdf_weighted, args.intensity_of_extinction, suppress_plotting = args.suppress_plotting)
            for key in statsdf:
                nodedf.loc[nodedf['taxid'] == node, key] = statsdf[key]
    
    print("\n" + "="*80)
    print("PHASE 2 COMPLETE!")
    print("="*80)
    print(f"\nAll analyses finished. Summary file: {summary_file}")
    print(f"Detailed plots and TSVs saved to: {per_clade_dir}\n")

if __name__== '__main__':
    main()