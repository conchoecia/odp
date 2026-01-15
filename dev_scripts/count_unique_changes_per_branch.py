#!/usr/bin/env python3
"""
Quick script to count unique changes per branch from changestrings.

Usage:
    python count_unique_changes_per_branch.py per_species_ALG_presence_fusions.tsv > unique_changes_summary.tsv
"""

import sys
import pandas as pd
from collections import defaultdict

def parse_changestring(changestring):
    """Parse a changestring and return list of (source, target, fusions, losses, splits)."""
    # Handle negative taxids (replace -- with -~)
    changestring_safe = changestring.replace('--', '-~')
    parts = changestring_safe.split('-')
    
    branches = []
    
    for i in range(0, len(parts)-1, 2):
        # Parse taxids
        source = parts[i].replace('~', '-')
        try:
            source_taxid = int(source)
        except:
            continue
            
        if i+2 >= len(parts):
            break
            
        target = parts[i+2].replace('~', '-')
        try:
            target_taxid = int(target)
        except:
            continue
        
        # Parse changes
        if i+1 >= len(parts):
            break
            
        changes_str = parts[i+1].lstrip('(').rstrip(')')
        
        try:
            components = changes_str.split(']|[')
            fusions = eval('[' + components[0].lstrip('[') + ']')
            losses = eval('[' + components[1] + ']')
            splits = eval('[' + components[2].rstrip(']') + ']')
            
            branches.append((source_taxid, target_taxid, fusions, losses, splits))
        except:
            continue
    
    return branches


def main():
    if len(sys.argv) < 2:
        print("Usage: python count_unique_changes_per_branch.py <per_species_ALG_presence_fusions.tsv>", file=sys.stderr)
        sys.exit(1)
    
    filename = sys.argv[1]
    df = pd.read_csv(filename, sep='\t')
    
    # Dictionary: (source_taxid, target_taxid) -> {'fusions': set(), 'losses': set(), 'splits': set()}
    branch_changes = defaultdict(lambda: {'fusions': set(), 'losses': set(), 'splits': set()})
    
    print(f"Parsing {len(df)} changestrings...", file=sys.stderr)
    for idx, row in df.iterrows():
        if idx % 100 == 0:
            print(f"  Processed {idx}/{len(df)} species...", file=sys.stderr)
        
        changestring = row['changestrings']
        if pd.isna(changestring):
            continue
        
        branches = parse_changestring(changestring)
        
        for source_taxid, target_taxid, fusions, losses, splits in branches:
            key = (source_taxid, target_taxid)
            
            for fusion in fusions:
                branch_changes[key]['fusions'].add(tuple(sorted(fusion)))
            
            for loss in losses:
                branch_changes[key]['losses'].add(loss)
            
            for split in splits:
                branch_changes[key]['splits'].add(split)
    
    # Output
    print("source_taxid\ttarget_taxid\tnum_unique_fusions\tnum_unique_losses\tnum_unique_splits\tunique_fusions\tunique_losses\tunique_splits")
    
    for (source, target) in sorted(branch_changes.keys()):
        info = branch_changes[(source, target)]
        
        num_fusions = len(info['fusions'])
        num_losses = len(info['losses'])
        num_splits = len(info['splits'])
        
        fusions_str = '; '.join(['+'.join(f) for f in sorted(info['fusions'])])
        losses_str = '; '.join(sorted(info['losses']))
        splits_str = '; '.join(sorted(info['splits']))
        
        print(f"{source}\t{target}\t{num_fusions}\t{num_losses}\t{num_splits}\t{fusions_str}\t{losses_str}\t{splits_str}")


if __name__ == '__main__':
    main()