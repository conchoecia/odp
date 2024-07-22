#!/usr/bin/env python

import argparse
from collections import Counter
import ete3


# get the path of this script, so we know where to look for the plotdfs file
# This block imports fasta-parser as fasta
import os
import sys

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, lil_matrix, save_npz, load_npz
import time

from PhyloTreeUMAP import (algcomboix_file_to_dict)

def parse_args():
    """
    Here, we need:
      - a path to a coo file
      - a path to a sample df
      - a path to the coo combination file
      - a list of the taxids for which we want to save the unique pairs tsv file
    """
    parser = argparse.ArgumentParser(description="Define features" )
    parser.add_argument("--coo_path",             type=str, help="Path to the coo file", required=True)
    parser.add_argument("--sample_df_path",       type=str, help="Path to the sample df", required=True)
    parser.add_argument("--coo_combination_path", type=str, help="Path to the coo combination file", required=True)
    parser.add_argument("--taxid_list",           type=str, help="Comma-separated list of taxids for which we want to save the unique pairs tsv file", required=False, default = "")

    args = parser.parse_args()
    # check that both the coo file and the df file exist. Same with coo combination file
    if not os.path.exists(args.coo_path):
        raise FileNotFoundError(f"{args.coo_path} does not exist")
    if not os.path.exists(args.sample_df_path):
        raise FileNotFoundError(f"{args.sample_df_path} does not exist")
    if not os.path.exists(args.coo_combination_path):
        raise FileNotFoundError(f"{args.coo_combination_path} does not exist")

    # The taxid list will be provided to us as a list of strings.
    # Here, we clean up the comma-separated list of taxids and convert it to a list of integers
    # If the taxid_list is empty, we will return an empty list.
    outlist = []
    if args.taxid_list != "":
        for taxid in args.taxid_list.split(","):
            if not taxid.isdigit():
                raise ValueError(f"taxid {taxid} is not an integer. Exiting.")
            outlist.append(int(taxid))
        args.taxid_list = outlist
    else:
        args.taxid_list = []
    return args

def load_coo(cdf, coofile, ALGcomboix, missing_value_as):
    """
    This loads a coo file and converts the missing values to the missing_value_as.
    Assumes that none of the values in the matrix will be -1, as this matrix should
      only contain positive integers.
    """
    print("loading the coo file")
    lil = load_npz(coofile).tolil()
    # check that the largest row index of the lil matrix is less than the largest index of cdf - 1
    if lil.shape[0] > max(cdf.index) + 1:
        raise ValueError(f"The largest row index of the lil matrix, {lil.shape[0]}, is greater than the largest index of cdf, {max(cdf.index)}. Exiting.")
    # check that the largest value of the ALGcomboix is less than the number of columns of the lil matrix - 1
    if max(ALGcomboix.values()) > lil.shape[1] - 1:
        raise ValueError(f"The largest value of the ALGcomboix, {max(ALGcomboix.values())}, is greater than the number of columns of the lil matrix, {lil.shape[1]}. Exiting.")

    # If the matrix is large, we have to convert the real zeros to -1 before we change to csf
    # we have to flip the values of the lil matrix
    print("setting zeros to -1")
    lil.data[lil.data == 0] = -1
    # We have to convert this to a dense matrix now. There is no way to modify the large values in a sparse matrix.
    print("Converting to a dense matrix. RAM will increase now.")
    # Goodbye, RAM.
    matrix = lil.toarray().astype(float)
    del lil
    # if the missing_values is "large", then we have to convert the 0 to the missing_value_as
    # Here we switch the representation, namely we don't have to access the data with .data now that this
    #  is a dense matrix.
    print(f"setting zeros to {missing_value_as}")
    matrix[matrix == 0] = missing_value_as
    # now we convert the -1s to 0
    print("converting -1s to 0")
    matrix[matrix == -1] = 0
    return matrix

def compute_statistics(col, inindex, outindex):
    """
    This function is applied pairwise to the columns of the matrix.
    It computes the statistics for the in and out samples.
    This could probably be done more efficiently in the future by using a DFS, but who knows really without optimizing.
        This is fine for now. :)
    """
    len_inindex = len(inindex)
    len_outindex = len(outindex)
    invals    = col.loc[inindex]
    outvals   = col.loc[outindex]
    notna_in  = invals.notna().sum()
    notna_out = outvals.notna().sum()

    mean_in   = invals.mean()
    sd_in     = invals.std()

    mean_out  = outvals.mean()
    sd_out    = outvals.std()

    return {
        "pair":          col.name,
        "notna_in":      notna_in,
        "notna_out":     notna_out,
        "mean_in":       mean_in,
        "sd_in":         sd_in,
        "mean_out":      mean_out,
        "sd_out":        sd_out,
        "occupancy_in":  notna_in  / len_inindex,
        "occupancy_out": notna_out / len_outindex}

def process_coo_file(sampledffile, ALGcomboixfile, coofile,
                     dfoutfilepath, missing_value_as = np.nan,
                     taxid_list = []):
    """
    Handles loading in the coo file and transforms it to a matrix that we can work with.

    Required args:
      - sampledffile: str, path to the sample df
        - ALGcomboixfile: str, path to the coo combination file
        - coofile: str, path to the coo file
        - dfoutfilepath: str, path to the output df file that we will write to
    Optional args:
      - missing_value_as: int, the value that we will use to represent missing values
    """
    ###make sure missing_value_as is an integer
    ## ensure that missing_value_as is an integer or np.nan
    #if not isinstance(missing_value_as, (int, np.nan)):
    #    raise ValueError(f"missing_value_as must be an integer or np.nan. Got {missing_value_as} instead. Exiting.")

    # check that all of the relevant files are actually present
    for filepath in [sampledffile, ALGcomboixfile, coofile]:
        if not os.path.exists(filepath):
            raise ValueError(f"The filepath {filepath} does not exist. Exiting.")

    # check that the file ending for the df outfile is .df
    if not dfoutfilepath.endswith(".df"):
        raise ValueError(f"The dfoutfilepath {dfoutfilepath} does not end with '.df'. Exiting.")

    # check that the type of taxid_list is a list
    if not isinstance(taxid_list, list):
        raise ValueError(f"taxid_list must be a list. Got {taxid_list} instead. Exiting.")

    # read in the sample dataframe. We will need this later
    cdf = pd.read_csv(sampledffile, sep = "\t", index_col = 0)
    # print the columns of cdf
    # Read in the ALGcomboixfile
    ALGcomboix = algcomboix_file_to_dict(ALGcomboixfile)
    ALGcomboix_inverse = {v: k for k, v in ALGcomboix.items()}
    matrix = load_coo(cdf, coofile, ALGcomboix, missing_value_as)
    print("done loading the matrix")
    print(f"dimensions of the matrix are: {matrix.shape}")
    print(matrix)
    # convert the matrix to a df, where species are labeled with the cdf samples, and the combinations are labeled with the ALGcomboix
    # we have to transpose the matrix to get the species as the rows
    print(cdf)
    df = pd.DataFrame(matrix, index=cdf.index,
                      columns=[i for i in range(matrix.shape[1])])

    # now find the things that define the species
    # convert the column "taxid_list" to a list of ints. use eval
    cdf["taxid_list"] = cdf["taxid_list"].apply(eval)
    # first collect all of the ncbi taxids from the cdf dict
    all_taxids = set()
    for i, row in cdf.iterrows():
        all_taxids.update(row["taxid_list"])

    # we should delete the columns from the df that are all nan for each taxid
    print("We are deleting the columns that are full of nans for each taxid. shape: ", df.shape)
    df = df.dropna(axis = 1, how = "all")
    print("New shape: ", df.shape)

    # If the user specified some taxids, we will only iterate over those taxids
    # Otherwise, go through all of the taxids in the cdf
    if len(taxid_list) > 0:
        iterate_taxids = taxid_list
        # make sure that all of the taxids that we want to iterate through are in all_taxids.
        for taxid in iterate_taxids:
            if taxid not in all_taxids:
                raise ValueError(f"taxid {taxid} is not in the taxids in the cdf. Exiting.")
    else:
        iterate_taxids = all_taxids

    entries = []
    # precompute the "in" and "out" indices
    # Precompute the inindex and outindex for each taxid
    # innan_dict is a dictionary that stores the columns that are nan for each taxid
    inindex_dict = {}
    outindex_dict = {}
    innan_dict   = {}

    start = time.time()
    print("Starting the precomputation of inindex and outindex")
    for taxid in iterate_taxids:
        inindex  = cdf[cdf["taxid_list"].apply(lambda x: taxid in x)].index
        outindex = cdf.index.difference(inindex)
        inindex_dict[ taxid ] = inindex
        outindex_dict[taxid ] = outindex
        # give me all of the colnames that are only nan for this taxid
        innan_dict[   taxid]  = df.loc[inindex,  df.columns].isna().all()
    print("  - Time to precompute inindex and outindex: ", time.time() - start)

    # load NCBI now that we know we will use it. Past the "taxid not in dataset" check
    NCBI = ete3.NCBITaxa()
    # For each taxid, get the pairs that are unique to this taxid
    counter = 1
    for taxid in iterate_taxids:
        nodename = NCBI.get_taxid_translator([taxid])[taxid]
        print(f"Starting taxid {taxid}, {nodename} ({counter} of {len(iterate_taxids)})")
        replace_dict = {" ": "", ",": "", ";": "", "(": "", ")": "", ".": "", "-": "", "_": ""}
        for k, v in replace_dict.items():
            nodename = nodename.replace(k, v)
        outprefix = f"{nodename}_{taxid}"
        outfile = outprefix + "_unique_pair_df.tsv.gz"
        if os.path.exists(outfile):
            print(f"{outfile} already exists. Skipping.")
            counter += 1
            continue

        # get the samples that have this taxid
        inindex  = inindex_dict[taxid]
        # if there are at least two samples, we should continue to analyze this
        if len(inindex) < 2:
            continue
        # get everything that isn't in inindex
        outindex = outindex_dict[taxid]
        # if the size of this is zero, we don't analyze it
        if len(outindex) == 0:
            continue
        # Apply function to each column. Only do it on the columns that are not all nan
        ignore_columns = innan_dict[taxid][innan_dict[taxid]].index
        # set everything after the first 50 columns to True
        #ignore_columns = df.columns[50:]
        print("  - We are filtering out ", len(ignore_columns), " columns.")
        stats = df.loc[:, df.columns.difference(ignore_columns)].apply(
            lambda col: compute_statistics(col, inindex, outindex)
        )
        print("This is stats")
        # turn this list of dictionaries into a df
        unique_pair_df = pd.DataFrame(stats.tolist())
        #unique_pair_df["mean_in_out_ratio"] = unique_pair_df["mean_in"] / unique_pair_df["mean_out"]
        # The following distribution is close to normal in log space, if not a little skewed.
        #  Because it has this property, for each pair we can measure where it falls in the distribution.
        #  If sd_in_out_ratio is nan, that is because the value does not appear in the out samples, meaning this is a new feature for this taxid.
        #  If sd_in_out_ratio is 0, this means that the variance was 0 for the in samples. This probably means that this pair was only measured twice?
        #  So from this number we can rank the pairs based on that ratio, then filter even more for well-represented pairs.
        #  From the other information we can get the things that do not occur in other samples, and that have a small SD or a small mean.
        #unique_pair_df["sd_in_out_ratio"]   = unique_pair_df["sd_in"]   / unique_pair_df["sd_out"]
        print("This is unique_pair_df")
        print(unique_pair_df)
        # save this so we can play with it
        unique_pair_df.to_csv(outfile, sep = "\t", index = False, compression = "gzip")
        counter += 1

    # convert entries to a df
    entries_df = pd.DataFrame(entries)
    # save this to a gzipped tsv
    entries_df.to_csv("unique_pairs.tsv.gz", sep = "\t", index = False, compression = "gzip")

def main():
    args = parse_args()
    print(args)
    #process_coo_file(sampledffile, ALGcomboixfile, coofile,
    #                 dfoutfilepath, missing_value_as = 9999999999)
    if len(args.taxid_list) > 0:
        process_coo_file(args.sample_df_path, args.coo_combination_path, args.coo_path, "test.df", taxid_list = args.taxid_list)
    else:
        process_coo_file(args.sample_df_path, args.coo_combination_path, args.coo_path, "test.df")

if __name__ == "__main__":
    main()