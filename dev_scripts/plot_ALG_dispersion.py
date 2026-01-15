#!/usr/bin/env python

"""
    Filename:   plot_ALG_dispersion.py
   File type:   python script (.py)
      Author:   darrin t schultz (github: @conchoecia)
Date created:   January 1st, 2026

Description:
  - This script creates dispersion plots showing ALG conservation categorized by
    median ALG conservation levels (0-20%, 40-60%, 80-100%).
  - For each genome, calculates the median conservation across all ALGs to determine
    which dispersion bin it belongs to.
  - Each datapoint in the plot is one ALG's conservation percentage in one genome.
  - ALGs are sorted by size (smallest to largest) on the x-axis.
  - Uses box-and-whisker plots colored by ALG.

Usage:
  python plot_ALG_dispersion.py -d /path/to/species_vs_ALG_rbh_files/ \\
                                 -a /path/to/ALG_database.rbh \\
                                 -n BCnS \\
                                 -o output_directory/
"""

# plotting options
import matplotlib.pyplot as plt

import argparse
import pandas as pd
import os
import sys
from glob import glob
from ete4 import NCBITaxa

# Add source directory to path for rbh_tools import
script_path = os.path.dirname(os.path.abspath(__file__))
source_path = os.path.join(script_path, "../source")
sys.path.insert(1, source_path)
import rbh_tools

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot ALG conservation dispersion categorized by median conservation levels.",
        epilog="This script analyzes species vs ALG RBH files to create dispersion plots.")
    
    parser.add_argument("-d", "--directory", required=True,
        help="Directory containing species vs ALG RBH files (e.g., species_vs_BCnS.rbh)")
    parser.add_argument("-a", "--alg_rbh", required=True,
        help="Path to ALG database RBH file (e.g., BCnSSimakov2022.rbh)")
    parser.add_argument("-n", "--algname", default="BCnS",
        help="Name of ALG database (default: BCnS)")
    parser.add_argument("-o", "--outdir", default="./alg_dispersion_plots",
        help="Output directory for plots (default: ./alg_dispersion_plots)")
    parser.add_argument("-m", "--minsig", type=float, default=0.05,
        help="Fisher's Exact Test p-value threshold for significance (default: 0.05)")
    parser.add_argument("-s", "--species",
        help="Specific species to plot (default: plot all species found in directory)")
    parser.add_argument("--metadata",
        help="TSV/CSV file with species metadata (e.g., per_species_ALG_presence_fusions.tsv)")
    parser.add_argument("--species_col", default="species",
        help="Column name for species in metadata file (default: 'species')")
    parser.add_argument("--lineage_col", default="taxidstring",
        help="Column name for lineage in metadata file (default: 'taxidstring')")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.directory):
        print(f"ERROR: Directory {args.directory} does not exist", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.alg_rbh):
        print(f"ERROR: ALG RBH file {args.alg_rbh} does not exist", file=sys.stderr)
        sys.exit(1)
    
    return args

def find_species_rbh_files(directory, algname):
    """Find all species vs ALG RBH files in the directory."""
    # Try the ODP naming pattern: BCnSSimakov2022_<species>_xy_reciprocal_best_hits.plotted.rbh
    pattern1 = os.path.join(directory, f"{algname}*_*_xy_reciprocal_best_hits.plotted.rbh")
    files = glob(pattern1)
    
    if len(files) == 0:
        # Try alternative pattern: <species>_vs_<algname>.rbh
        pattern2 = os.path.join(directory, f"*_vs_{algname}.rbh")
        files = glob(pattern2)
    
    species_to_file = {}
    for filepath in files:
        basename = os.path.basename(filepath)
        
        # Parse ODP format: BCnSSimakov2022_<species>_xy_reciprocal_best_hits.plotted.rbh
        if basename.startswith(algname) or basename.startswith("BCnS"):
            # Remove prefix and suffix
            species = basename.replace(f"{algname}Simakov2022_", "")
            species = species.replace("BCnSSimakov2022_", "")
            species = species.replace("_xy_reciprocal_best_hits.plotted.rbh", "")
            species = species.replace("_yx_reciprocal_best_hits.plotted.rbh", "")
        else:
            # Parse simpler format: <species>_vs_BCnS.rbh
            species = basename.replace(f"_vs_{algname}.rbh", "")
        
        species_to_file[species] = filepath
    
    return species_to_file

def taxidstring_to_lineage(taxidstring, ncbi):
    """
    Convert semicolon-separated taxid string to human-readable lineage.
    Handles custom taxids: -67 (Myriazoa), -68 (Parahoxozoa)
    
    Args:
        taxidstring: String like "1;131567;2759;33154;33208;..."
        ncbi: NCBITaxa instance
    
    Returns:
        String like "root;cellular organisms;Eukaryota;Opisthokonta;Metazoa;..."
    """
    custom_taxid_names = {
        -67: "Myriazoa",
        -68: "Parahoxozoa"
    }
    
    # Parse taxids from string
    taxids = [int(x) for x in taxidstring.split(";")]
    
    # Separate valid NCBI taxids from custom ones
    valid_taxids = [t for t in taxids if t > 0]
    
    # Get names for valid taxids
    taxid_to_name = {}
    if valid_taxids:
        taxid_to_name = ncbi.get_taxid_translator(valid_taxids)
    
    # Build lineage string
    lineage_parts = []
    for taxid in taxids:
        if taxid in custom_taxid_names:
            lineage_parts.append(custom_taxid_names[taxid])
        elif taxid in taxid_to_name:
            lineage_parts.append(taxid_to_name[taxid])
        else:
            lineage_parts.append(str(taxid))  # fallback to numeric
    
    return ";".join(lineage_parts)

def parse_metadata(metadata_file, species_col="species", lineage_col="taxidstring"):
    """
    Parse metadata file to extract species and lineage information.
    
    Args:
        metadata_file: Path to TSV or CSV file
        species_col: Column name containing species names
        lineage_col: Column name containing lineage (e.g., 'taxidstring')
    
    Returns:
        Dictionary mapping species name -> lineage string
    """
    # Detect file format by extension
    if metadata_file.endswith('.csv'):
        df = pd.read_csv(metadata_file)
    else:
        df = pd.read_csv(metadata_file, sep='\t')
    
    # Check required columns exist
    if species_col not in df.columns:
        print(f"ERROR: Column '{species_col}' not found in metadata file", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        sys.exit(1)
    if lineage_col not in df.columns:
        print(f"ERROR: Column '{lineage_col}' not found in metadata file", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        sys.exit(1)
    
    # Convert taxidstrings to human-readable lineages
    print("Converting taxidstrings to human-readable lineages...", file=sys.stderr)
    ncbi = NCBITaxa()
    species_to_lineage = {}
    
    for idx, row in df.iterrows():
        species = row[species_col]
        taxidstring = row[lineage_col]
        
        if pd.isna(taxidstring):
            species_to_lineage[species] = "Unknown"
        else:
            species_to_lineage[species] = taxidstring_to_lineage(str(taxidstring), ncbi)
    
    print(f"Loaded lineage data for {len(species_to_lineage)} species", file=sys.stderr)
    return species_to_lineage

def calculate_alg_conservation_per_species(rbh_file, algname, sorted_algs, minsig=0.05):
    """
    Calculate conservation percentage for each ALG in a species.
    
    Returns:
        dict: {alg_name: conservation_fraction}
    """
    try:
        df = pd.read_csv(rbh_file, sep="\t")
    except Exception as e:
        print(f"ERROR reading {rbh_file}: {e}", file=sys.stderr)
        return {}
    
    if len(df) == 0:
        return {}
    
    # Filter by significance if whole_FET column exists
    if 'whole_FET' in df.columns:
        df = df[df['whole_FET'] <= minsig]
    
    if len(df) == 0:
        return {}
    
    # Get species name from columns
    species_cols = [col for col in df.columns if col.endswith('_gene') and algname not in col]
    if len(species_cols) == 0:
        return {}
    
    species = species_cols[0].replace('_gene', '')
    
    # Calculate conservation per ALG
    alg_conservation = {}
    
    for alg in sorted_algs:
        # Filter rows for this ALG
        alg_df = df[df['gene_group'] == alg]
        
        if len(alg_df) > 0:
            # Conservation = number of genes found / total genes in this ALG
            # This is the fraction of ALG genes that are conserved in this species
            alg_conservation[alg] = len(alg_df)
        else:
            alg_conservation[alg] = 0
    
    return alg_conservation

def plot_dispersion_by_alg(species_to_rbh, alg_df, algname, minsig, outdir, species_to_lineage=None):
    """
    Create ALG dispersion plot for all species.
    
    For each species:
    1. Calculate conservation for each ALG
    2. Use MEDIAN ALG conservation to determine bin (0-20%, 40-60%, 80-100%)
    3. Each ALG's conservation becomes a datapoint in that bin
    
    Args:
        species_to_lineage: Optional dict mapping species -> lineage string
    """
    
    # Sort ALGs by size (smallest to largest)
    sorted_algs = alg_df['ALGname'].tolist()
    alg_sizes = dict(zip(alg_df['ALGname'], alg_df['Size']))
    alg_colors = dict(zip(alg_df['ALGname'], alg_df['Color']))
    
    # Total genes per ALG (from the ALG database)
    total_genes_per_alg = alg_sizes.copy()
    
    # Define dispersion bins based on MEDIAN ALG conservation
    # High dispersion = low conservation, low dispersion = high conservation
    dispersion_bins = {
        '0%-20% dispersion': (0.8, 1.0),    # Low dispersion = high conservation
        '40%-60% dispersion': (0.4, 0.6),   # Medium dispersion = medium conservation
        '80%-100% dispersion': (0.0, 0.2)   # High dispersion = low conservation
    }
    
    # Data structure: {bin_name: {alg: [conservation fractions]}}
    bin_to_alg_data = {bin_name: {alg: [] for alg in sorted_algs} for bin_name in dispersion_bins.keys()}
    
    # Process each species
    num_species = len(species_to_rbh)
    num_binned = 0
    
    # Track all species with their median conservation for ranking
    species_conservation_data = []
    
    print(f"Processing {num_species} species...", file=sys.stderr)
    
    for species, rbh_file in species_to_rbh.items():
        # Get raw gene counts per ALG for this species
        alg_gene_counts = calculate_alg_conservation_per_species(
            rbh_file, algname, sorted_algs, minsig)
        
        if not alg_gene_counts or sum(alg_gene_counts.values()) == 0:
            continue
        
        # Convert to fractions (conserved genes / total genes in ALG)
        alg_conservation_fractions = {}
        for alg in sorted_algs:
            count = alg_gene_counts.get(alg, 0)
            total = total_genes_per_alg.get(alg, 1)
            alg_conservation_fractions[alg] = count / total if total > 0 else 0.0
        
        # Calculate MEDIAN conservation across all ALGs
        conservation_values = [v for v in alg_conservation_fractions.values() if v > 0]
        if len(conservation_values) == 0:
            continue
        
        median_conservation = pd.Series(conservation_values).median()
        
        # Store species data for ranking (including per-ALG conservation)
        species_data = {
            'species': species,
            'median_conservation': median_conservation,
            'num_algs_detected': len(conservation_values),
            'mean_conservation': pd.Series(conservation_values).mean(),
            'std_conservation': pd.Series(conservation_values).std()
        }
        # Add lineage if available
        if species_to_lineage and species in species_to_lineage:
            species_data['lineage'] = species_to_lineage[species]
        
        # Add per-ALG conservation values
        for alg in sorted_algs:
            species_data[f'ALG_{alg}'] = alg_conservation_fractions.get(alg, 0.0)
        
        species_conservation_data.append(species_data)
        
        # Determine which bin this species belongs to
        target_bin = None
        for bin_name, (lower, upper) in dispersion_bins.items():
            if lower <= median_conservation <= upper:
                target_bin = bin_name
                break
        
        if target_bin is None:
            continue
        
        num_binned += 1
        
        # Add each ALG's conservation as a datapoint
        for alg in sorted_algs:
            frac = alg_conservation_fractions.get(alg, 0.0)
            bin_to_alg_data[target_bin][alg].append(frac)
    
    print(f"Binned {num_binned}/{num_species} species", file=sys.stderr)
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(22.5, 6))
    fig.suptitle(f"ALG conservation dispersion (FET p<{minsig})", fontsize=14)
    
    # Plot each bin
    for idx, bin_name in enumerate(['0%-20% dispersion', '40%-60% dispersion', '80%-100% dispersion']):
        ax = axes[idx]
        ax.set_title(bin_name)
        ax.set_xlabel(f"{algname} ALG")
        if idx == 0:
            ax.set_ylabel("Fraction of ALG genes conserved")
        
        # Prepare data for box plot
        plot_data = []
        plot_positions = []
        plot_labels = []
        plot_colors = []
        
        for pos, alg in enumerate(sorted_algs):
            data = bin_to_alg_data[bin_name][alg]
            if len(data) > 0:
                plot_data.append(data)
                plot_positions.append(pos)
                plot_labels.append(alg)
                plot_colors.append(alg_colors.get(alg, '#1f77b4'))
        
        if len(plot_data) > 0:
            # Create box-and-whisker plot
            bp = ax.boxplot(plot_data, positions=plot_positions,
                           widths=0.6, patch_artist=True,
                           showfliers=True, notch=False)
            
            # Color boxes by ALG
            for patch, color in zip(bp['boxes'], plot_colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            ax.set_xticks(plot_positions)
            ax.set_xticklabels(plot_labels, rotation=90, fontsize=8, ha='center')
        else:
            ax.text(0.5, 0.5, 'No data in this bin', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
        
        ax.set_ylim(-0.05, 1.05)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"ALG_dispersion_{algname}.pdf")
    plt.savefig(outfile, format='pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved plot to {outfile}", file=sys.stderr)
    
    # Create ranking dataframe and save extremes
    if len(species_conservation_data) > 0:
        ranking_df = pd.DataFrame(species_conservation_data)
        ranking_df = ranking_df.sort_values('median_conservation', ascending=False)
        
        # Save top 10 most conserved and 10 most dispersed
        outfile_ranking = os.path.join(outdir, f"ALG_conservation_ranking_{algname}.tsv")
        ranking_df.to_csv(outfile_ranking, sep="\t", index=False)
        print(f"Saved full ranking to {outfile_ranking}", file=sys.stderr)
        
        # Print summary
        print("\n=== TOP 10 MOST CONSERVED GENOMES ===", file=sys.stderr)
        print("(High median ALG conservation = low dispersion)", file=sys.stderr)
        for idx, row in ranking_df.head(10).iterrows():
            lineage_str = f" [{row['lineage']}]" if 'lineage' in row and pd.notna(row['lineage']) else ""
            print(f"  {row['species']}: {row['median_conservation']:.3f} (n={int(row['num_algs_detected'])} ALGs){lineage_str}", file=sys.stderr)
        
        print("\n=== TOP 10 MOST DISPERSED GENOMES ===", file=sys.stderr)
        print("(Low median ALG conservation = high dispersion)", file=sys.stderr)
        for idx, row in ranking_df.tail(10).iterrows():
            lineage_str = f" [{row['lineage']}]" if 'lineage' in row and pd.notna(row['lineage']) else ""
            print(f"  {row['species']}: {row['median_conservation']:.3f} (n={int(row['num_algs_detected'])} ALGs){lineage_str}", file=sys.stderr)
        print("", file=sys.stderr)

def main():
    args = parse_args()
    
    # Parse ALG database
    print(f"Loading ALG database: {args.alg_rbh}", file=sys.stderr)
    try:
        alg_df = rbh_tools.parse_ALG_rbh_to_colordf(args.alg_rbh)
        print(f"  Found {len(alg_df)} ALGs", file=sys.stderr)
    except Exception as e:
        print(f"ERROR parsing ALG database: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Find species RBH files
    print(f"Searching for species vs {args.algname} RBH files in {args.directory}", file=sys.stderr)
    species_to_rbh = find_species_rbh_files(args.directory, args.algname)
    
    if len(species_to_rbh) == 0:
        print(f"ERROR: No RBH files found matching pattern *_vs_{args.algname}.rbh", file=sys.stderr)
        sys.exit(1)
    
    print(f"  Found {len(species_to_rbh)} species", file=sys.stderr)
    
    # Filter to specific species if requested
    if args.species:
        if args.species in species_to_rbh:
            species_to_rbh = {args.species: species_to_rbh[args.species]}
            print(f"  Filtering to species: {args.species}", file=sys.stderr)
        else:
            print(f"ERROR: Species {args.species} not found in directory", file=sys.stderr)
            print(f"  Available species: {', '.join(sorted(species_to_rbh.keys()))}", file=sys.stderr)
            sys.exit(1)
    
    # Load metadata if provided
    species_to_lineage = None
    if args.metadata:
        if not os.path.exists(args.metadata):
            print(f"ERROR: Metadata file {args.metadata} does not exist", file=sys.stderr)
            sys.exit(1)
        print(f"\nLoading metadata from {args.metadata}...", file=sys.stderr)
        species_to_lineage = parse_metadata(args.metadata, args.species_col, args.lineage_col)
    
    # Create dispersion plot
    plot_dispersion_by_alg(species_to_rbh, alg_df, args.algname, args.minsig, args.outdir, species_to_lineage)
    
    print("Done!", file=sys.stderr)

if __name__ == "__main__":
    main()
