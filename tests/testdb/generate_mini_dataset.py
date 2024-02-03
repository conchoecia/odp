#!/usr/bin/env python

"""
Takes in a dataset and generates a version that is only 10% as large, in terms of proteins, as the original.
The resulting genome assembly fasta file will be as small as possible.
"""

# This block imports fasta-parser as fasta
import os
import sys
this_path = os.path.dirname(os.path.abspath(__file__))
dependencies_path = os.path.join(this_path, "../../dependencies/fasta-parser")
sys.path.insert(1, dependencies_path)
import fasta

import argparse
import gzip
import pandas as pd
import random

def parse_args():
    """
    In this, we need to get the file paths for:
      1. The protein fasta file, -p --protein
      2. The genome assembly fasta file -g --genome
      3. The chrom file -c --chrom
    We also need:
      - The output prefix --prefix
    """
    parser = argparse.ArgumentParser(description="Takes in a dataset and generates a version that is only 10% as large, in terms of proteins, as the original. The resulting genome assembly fasta file will be as small as possible.")
    parser.add_argument("-p", "--protein", help="The protein fasta file", required=True)
    parser.add_argument("-g", "--genome", help="The genome assembly fasta file", required=True)
    parser.add_argument("-c", "--chrom", help="The chrom file", required=True)
    parser.add_argument("--prefix", help="The output prefix", required=True)

    args = parser.parse_args()
    # check that protein file exists
    if not os.path.exists(args.protein):
        parser.error(f"The protein file {args.protein} does not exist.")
    # check that genome file exists
    if not os.path.exists(args.genome):
        parser.error(f"The genome file {args.genome} does not exist.")
    # check that the chrom file exists
    if not os.path.exists(args.chrom):
        parser.error(f"The chrom file {args.chrom} does not exist.")
    return args


def write_new_protein_file(protein_file, output_protein_path, protein_sample):
    """
    Takes in a sample of proteins and writes them to a new, gzipped file
    """
    # write out the sampled proteins to a new, gzipped file. If the file exists already, overwrite it.
    # open the output file as a gzip file
    outhandle = gzip.open(output_protein_path, "wt")
    for record in fasta.parse(protein_file):
        if record.id in protein_sample:
            outhandle.write(record.format(wrap=80))
    outhandle.close()

def write_new_fasta_file(chromdf, output_genome_path):
    """
    Write a pseudo-genome assembly file
    """
    # get all of the scaffolds in the chrom file, and a dictionary of the max length for that scaf
    chrom_to_size = chromdf.groupby("scaf")["stop"].max().to_dict()
    # open the genome as a gzip file
    outhandle = gzip.open(output_genome_path, "wt")
    for this_scaf, this_size in chrom_to_size.items():
        # make a string that is slightly longer than the size of the max scaffold size
        outstring = "GATTACA" * ((this_size // 7) + 2)
        record = fasta.Record(id=this_scaf,
                seq=outstring)
        outhandle.write(record.format(wrap=80))
    outhandle.close()

def generate_and_write_new_chrom_file(chrom_file, output_chrom_path, protein_sample) -> pd.DataFrame:
    """
    Generates a new chrom file and returns a pandas dataframe
    """
    # read in the chrom file with pandas. use gzip to open it if needed
    # columns are: ["protein", "scaf", "strand", "start", stop"]
    df = pd.read_csv(chrom_file, sep="\t", header=None)
    df.columns = ["protein", "scaf", "strand", "start", "stop"]
    # only get the rows that are in the protein sample
    df = df[df["protein"].isin(protein_sample)].reset_index(drop=True)
    # make a column that is length, and divide it by 10, then round it to an integer
    df["length"] = (df["stop"] - df["start"]) // 100
    # if the value of length is zero, make it 1000
    df.loc[df["length"] == 0, "length"] = 10
    # Groupby the chromosome names and do a running sum of the lengths.
    # The running sum should be a new start and a new stop.
    df["stop"] = df.groupby("scaf")["length"].cumsum()
    df["start"] = df["stop"] - df["length"]
    # drop the length column
    df = df.drop(columns=["length"])
    # write the new chrom file
    df.to_csv(output_chrom_path, sep="\t", header=False, index=False)
    return df

def main():
    args = parse_args()

    # read in all the protein sequence IDs
    protein_ids = set()
    for record in fasta.parse(args.protein):
        if record.id in protein_ids:
            raise IOerror("Duplicate protein ID: {record.id} --- We saw this more than once.")
        protein_ids.add(record.id)

    # sample 10% of the proteins
    protein_sample = set(random.sample(protein_ids, int(0.1 * len(protein_ids))))

    # write out the sampled proteins to a new file.
    output_protein_path = f"{args.prefix}.pep.gz"
    write_new_protein_file(args.protein, output_protein_path, protein_sample)

    # generate a new chrom file. Results in a df.
    output_chrom_path = f"{args.prefix}.chrom.gz"
    chromdf = generate_and_write_new_chrom_file(args.chrom, output_chrom_path, protein_sample)

    # write out a fake genome assembly fasta file
    output_genome_path = f"{args.prefix}.fa.gz"
    write_new_fasta_file(chromdf, output_genome_path)

if __name__ == "__main__":
    main()