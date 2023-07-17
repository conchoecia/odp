#!/usr/bin/env python

"""
Reads in a genome and performs an emergent cluster analysis
"""

import argparse
import pandas as pd
# import biopython module for reading fasta files
from Bio import SeqIO

# set up cmd line arguments to read in these files:
#  - a .rbh file with the groupby and color columns
#  - chrom file of sp1
#  - chrom file of sp2
#  - assembly fasta file of sp1
#  - assembly fasta file of sp2

# also add arguments for these parameters
#  - sp1 name
#  - sp2 name
#  - sp1 chrom of interest
#  - sp2 chrom of interest
#  - output prefix
#  - comma-separated groups to include in PCA

def parse_args():
    """
    parse the command line arguments
    """
    parser = argparse.ArgumentParser(description="Run PCA on a set of orthologs")
    parser.add_argument("-r",  "--rbh", required=True, help="Input .rbh file")
    parser.add_argument("-s1", "--species1", required=True, help="Species 1 name")
    parser.add_argument("-s2", "--species2", required=True, help="Species 2 name")
    parser.add_argument(       "--chrfile1", required=True, help="Species 1 chrom file.")
    parser.add_argument(       "--chrfile2", required=True, help="Species 2 chrom file.")  
    parser.add_argument("-f1", "--fasta1", required=True, help="Species 1 assembly fasta file")
    parser.add_argument("-f2", "--fasta2", required=True, help="Species 1 assembly fasta file")
    parser.add_argument("-c1", "--chrom1", required=True, help="Species 1 chromosome of interest")
    parser.add_argument("-c2", "--chrom2", required=True, help="Species 2 chromosome of interest")
    parser.add_argument("-g",  "--groups", required=True, help="Comma-separated groups to include in PCA")
    parser.add_argument("-o", "--output", required=True, help="Output prefix")
    args = parser.parse_args()
    return args

def get_chromosome_lengths(fasta_file):
    """
    Use biopython to read in the fasta file
     and return a dictionary of chromosome lengths
    """
    chrom_lengths = {}
    with open(fasta_file) as handle:
        records = SeqIO.parse(handle, "fasta")
        for record in records:
            chrom_lengths[record.id] = len(record.seq)
    return chrom_lengths

def get_chrom_lens(sp, fastafile):
    """
    reads in a fasta file and returns a dict
    of the chrom length. Also handles the case
    where the chrom lengths have already been seen.

    In this case, we will read in the chrom lengths
    from an existing file.
    """
    # get the chrom lengths of the two species unless they exist in the pwd
    rbh1file = "{}_chromsizes.tsv".format(sp)
    sp1_chrom_lengths = {}
    if not os.path.exists(rbh1file):
        sp1_chrom_lengths = get_chromosome_lengths(fastafile)
        with open(rbh1file, "w") as outfile:
            for entry in sp1_chrom_lengths:
                outfile.write("{}\t{}\n".format(entry, sp1_chrom_lengths[entry]))
    else:
        with open(rbh1file) as infile:
            for line in infile:
                line = line.strip().split("\t")
                sp1_chrom_lengths[line[0]] = int(line[1])
    return sp1_chrom_lengths

def main():
    args = parse_args()
    # get the chrom lengths of the two species unless they exist in the pwd
    rbh1file = "{}_chromsizes.tsv.format(args.species1)"
    sp1_chrom_lengths = get_chrom_lens(args.species1, args.fasta1)
    sp2_chrom_lengths = get_chrom_lens(args.species2, args.fasta2)

    # get the rbh files as pandas dataframes
    rbh1 = pd.read_csv(args.rbh, sep="\t")
    rbh1 = rbh1[rbh1["{}_scaf".format(args.species1)] == args.chrom1]
    # only keep the rows that are the chromosome of interest for this species
    rbh2 = pd.read_csv(args.rbh, sep="\t")
    rbh2 = rbh2[rbh2["{}_scaf".format(args.species2)] == args.chrom2]
    print(rbh1)

    # read in the chrom file as a dataframe
    # read in the rbh file as a dataframe
    pass

if __name__ == "__main__":
    main()
