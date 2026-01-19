#!/usr/bin/env python

"""
Join two TSV files on Assembly Accession.

Joins:
- subsample_allsamples.neighbors_250.mind_1.0.missing_large.supplemented.df (sample column, 3rd '-'-delimited field)
- merged.tsv (Assembly Accession column)

Output has subsample file columns first, then merged.tsv columns.
"""

import argparse
import pandas as pd
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--subsample', required=True,
                        help='Path to subsample .df file')
    parser.add_argument('--merged', required=True,
                        help='Path to merged.tsv file')
    parser.add_argument('--output', required=True,
                        help='Path to output joined TSV file')
    parser.add_argument('--join-type', default='left',
                        choices=['left', 'right', 'inner', 'outer'],
                        help='Type of join (default: left)')
    return parser.parse_args()


def main():
    args = parse_args()
    
    print(f"Reading subsample file: {args.subsample}")
    subsample_df = pd.read_csv(args.subsample, sep='\t')
    print(f"  {len(subsample_df)} rows, {len(subsample_df.columns)} columns")
    
    print(f"Reading merged file: {args.merged}")
    merged_df = pd.read_csv(args.merged, sep='\t')
    print(f"  {len(merged_df)} rows, {len(merged_df.columns)} columns")
    
    # Extract 3rd field from sample column (0-indexed, so field 2)
    print("\nExtracting Assembly Accession from sample column...")
    subsample_df['assembly_accession_join'] = subsample_df['sample'].str.split('-').str[2]
    
    # Add underscore after 3-letter prefix (GCF or GCA)
    # Convert GCF964019385.1 -> GCF_964019385.1 to match format in merged.tsv
    subsample_df['assembly_accession_join'] = subsample_df['assembly_accession_join'].str.replace(
        r'^(GC[AF])([0-9])', r'\1_\2', regex=True
    )
    
    # Check for missing values
    missing_count = subsample_df['assembly_accession_join'].isna().sum()
    if missing_count > 0:
        print(f"  WARNING: {missing_count} rows have missing assembly accession")
    
    print(f"  Example sample: {subsample_df['sample'].iloc[0]}")
    print(f"  Extracted accession: {subsample_df['assembly_accession_join'].iloc[0]}")
    
    # Perform join
    print(f"\nPerforming {args.join_type} join...")
    joined_df = subsample_df.merge(
        merged_df,
        left_on='assembly_accession_join',
        right_on='Assembly Accession',
        how=args.join_type,
        suffixes=('', '_merged')
    )
    
    # Drop the temporary join column
    joined_df = joined_df.drop(columns=['assembly_accession_join'])
    
    print(f"  Result: {len(joined_df)} rows, {len(joined_df.columns)} columns")
    
    # Report join statistics
    if args.join_type == 'left':
        unmatched = joined_df['Assembly Accession'].isna().sum()
        if unmatched > 0:
            print(f"  WARNING: {unmatched} rows from subsample file had no match in merged file")
    
    # Save output
    print(f"\nSaving output to: {args.output}")
    joined_df.to_csv(args.output, sep='\t', index=False)
    print("Done!")


if __name__ == '__main__':
    main()
