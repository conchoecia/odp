#!/usr/bin/env python

"""
The purpose of this script is to time the rbh_to_distance_gbgz function
 and to see if the disk locks up.

When run with snakemake, all of the jobs stop running (seemingly due to IO issues).
"""

import argparse

from PhyloTreeUMAP import rbh_to_distance_gbgz

def parse_args():
    """
    Args we need:
      - rbh_file
      - ALGname
      - samplename
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--rbh_file", type=str, required=True)
    parser.add_argument("--ALGname", type=str, required=True)
    parser.add_argument("--samplename", type=str, required=True)
    args = parser.parse_args()
    ## if any of the required args are missing, print an error
    #if not args.rbh_file:
    #    parser.error("Please provide a rbh_file")
    #if not args.ALGname:
    #    parser.error("Please provide an ALGname")
    #if not args.samplename:
    #    parser.error("Please provide a samplename")
    return args

def main():
    args = parse_args()
    outfile = f"{args.samplename}.gb.gz"
    print(f"processing {args.rbh_file}")
    rbh_to_distance_gbgz(args.rbh_file,
                         outfile,
                         args.ALGname)

if __name__ == '__main__':
    main()