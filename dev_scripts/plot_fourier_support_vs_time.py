#!/usr/bin/env python
"""
Plot periodicity support vs time window from fourier analysis results.

This script reads *_rsims.tsv files from fourier analyses run with different
time windows and plots how the support for different periods changes as the
time window is expanded.

Usage:
    python plot_fourier_support_vs_time.py -d /path/to/fourier_analysis/ -o output.pdf
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Plot periodicity support vs time window')
    parser.add_argument('-d', '--directory', type=str, required=True,
                        help='Path to directory containing *_rsims.tsv files')
    parser.add_argument('-o', '--output', type=str, default='support_vs_time_window.pdf',
                        help='Output PDF filename (default: support_vs_time_window.pdf)')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Find all TSV files
    files = glob.glob(os.path.join(args.directory, '*_rsims.tsv'))
    print(f"Found {len(files)} total files")
    
    dispersal_files = [f for f in files if 'dispersal' in f]
    fusion_files = [f for f in files if 'fusion' in f]
    print(f"  - {len(fusion_files)} fusion files")
    print(f"  - {len(dispersal_files)} dispersal files")
    
    if len(files) == 0:
        print("ERROR: No *_rsims.tsv files found!")
        exit(1)
    
    # Read all TSV files
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        dfs.append(df)
    
    all_data = pd.concat(dfs, ignore_index=True)
    print(f"\nTotal rows: {len(all_data)}")
    
    # Extract change_type from the clade column (filename) instead of the change_type column
    # which may be NaN for dispersals due to a bug
    all_data['change_type_fixed'] = all_data['clade'].apply(
        lambda x: 'dispersals' if 'dispersal' in str(x).lower() else 'fusions'
    )
    
    # Extract max_time from clade column
    all_data['max_time'] = all_data['clade'].str.extract(r'(\d+)_(?:padded|unpadded)', expand=False)
    
    # Drop old files without time window numbers
    before = len(all_data)
    all_data = all_data.dropna(subset=['max_time'])
    print(f"Dropped {before - len(all_data)} rows without time numbers")
    all_data['max_time'] = all_data['max_time'].astype(int)
    
    # Check change_type
    print(f"\nchange_type_fixed value counts:")
    print(all_data['change_type_fixed'].value_counts())
    
    # Filter for observed values only
    observed = all_data[all_data['obs_exp'] == 'observed'].copy()
    print(f"\nObserved rows: {len(observed)}")
    print(f"Observed change_type_fixed counts:")
    print(observed['change_type_fixed'].value_counts())
    
    # Bin periods into 10 Myr bins
    observed['period_bin'] = (observed['period'] // 10) * 10
    observed['period_bin_label'] = observed['period_bin'].astype(int).astype(str) + '-' + (observed['period_bin'] + 10).astype(int).astype(str) + ' Myr'
    
    # Group by change_type, period_bin, and max_time, then average support
    binned = observed.groupby(['change_type_fixed', 'period_bin', 'period_bin_label', 'max_time'])['support_mean'].mean().reset_index()
    
    # Separate by change type
    fusions = binned[binned['change_type_fixed'] == 'fusions']
    dispersals = binned[binned['change_type_fixed'] == 'dispersals']
    
    print(f"\nFusions rows: {len(fusions)}")
    print(f"Dispersals rows: {len(dispersals)}")
    
    if len(dispersals) > 0:
        print(f"\nDispersal support range: {dispersals['support_mean'].min():.3f} - {dispersals['support_mean'].max():.3f}")
        print(f"Dispersal periods: {sorted(dispersals['period_bin_label'].unique())}")
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot fusions
    for period_bin_label in sorted(fusions['period_bin_label'].unique()):
        period_data = fusions[fusions['period_bin_label'] == period_bin_label].sort_values('max_time')
        ax1.plot(period_data['max_time'], period_data['support_mean'], 
                marker='o', linewidth=2, markersize=6, label=period_bin_label)
    
    ax1.set_xlabel('Maximum Time Window (MYA)', fontsize=12)
    ax1.set_ylabel('Support (mean)', fontsize=12)
    ax1.set_title('Fusions: Periodicity Support vs Time Window (10 Myr bins)', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Plot dispersals
    if len(dispersals) > 0:
        for period_bin_label in sorted(dispersals['period_bin_label'].unique()):
            period_data = dispersals[dispersals['period_bin_label'] == period_bin_label].sort_values('max_time')
            ax2.plot(period_data['max_time'], period_data['support_mean'], 
                    marker='o', linewidth=2, markersize=6, label=period_bin_label)
        ax2.legend()
    else:
        ax2.text(0.5, 0.5, 'No dispersal data found', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=14)
    
    ax2.set_xlabel('Maximum Time Window (MYA)', fontsize=12)
    ax2.set_ylabel('Support (mean)', fontsize=12)
    ax2.set_title('Dispersals: Periodicity Support vs Time Window (10 Myr bins)', fontsize=14, fontweight='bold')
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {args.output}")

if __name__ == '__main__':
    main()
