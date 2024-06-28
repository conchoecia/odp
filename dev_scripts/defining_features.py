#!/usr/bin/env python

from collections import Counter
import argparse

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
    """
    parser = argparse.ArgumentParser(description="Define features" )
    parser.add_argument("--coo_path",             type=str, help="Path to the coo file", required=True)
    parser.add_argument("--sample_df_path",       type=str, help="Path to the sample df", required=True)
    parser.add_argument("--coo_combination_path", type=str, help="Path to the coo combination file", required=True)

    args = parser.parse_args()
    # check that both the coo file and the df file exist. Same with coo combination file
    if not os.path.exists(args.coo_path):
        raise FileNotFoundError(f"{args.coo_path} does not exist")
    if not os.path.exists(args.sample_df_path):
        raise FileNotFoundError(f"{args.sample_df_path} does not exist")
    if not os.path.exists(args.coo_combination_path):
        raise FileNotFoundError(f"{args.coo_combination_path} does not exist")
    return parser.parse_args()

# Function to get column names of the lowest 1% values for each row
def get_lowest_1_percent_columns(row, percentage=1):
    threshold = np.nanpercentile(row, percentage)
    return row.index[row < threshold].tolist()

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

def process_coo_file(sampledffile, ALGcomboixfile, coofile,
                     dfoutfilepath, missing_value_as = np.nan):
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
    df = pd.DataFrame(matrix, index = cdf.index,
                      columns = [i for i in range(matrix.shape[1])])
    #columns = [ALGcomboix_inverse[i] for i in range(matrix.shape[1])])
    # I think this is messing up the indices
    ## remove the columns in which all values are nan
    #print(f"dropping empties. dimensions: {df.shape}")
    #df = df.dropna(axis=1, how='all')
    #print(f"new dimensions: {df.shape}")

    # get the lowest 1% columns
    # time this
    start = time.time()
    print("starting the lowest 1% columns in apply")
    lowest_1_percent_columns = df.apply(get_lowest_1_percent_columns, axis=1)
    # the above is a single column, convert to a df with the "sample" column of the cdf
    lowest_1_percent_columns = pd.DataFrame(lowest_1_percent_columns, columns = ["pairs"])
    lowest_1_percent_columns["sample"] = cdf["sample"]
    # add this as a column to the cdf
    cdf["lowest_1_percent_columns"] = lowest_1_percent_columns["pairs"]
    print(cdf)
    # change the single column name to "pairs"
    print(f"Time to get lowest 1% columns: {time.time() - start}")
    # convert the indices of the lowest 1% columns to the sample names
    pair_to_sample = {}
    for i, row in lowest_1_percent_columns.iterrows():
        sample = row["sample"]
        for j in row["pairs"]:
            if j not in pair_to_sample:
                pair_to_sample[j] = []
            pair_to_sample[j].append(sample)
    # save pair_to_sample to a file, gzip it
    with open("pair_to_sample.txt", "w") as f:
        for k, v in pair_to_sample.items():
            f.write(f"{k}\t{v}\n")

    # save a pdf with the distribution of the values for the first 10 species. Make one pdf, and stack the plots on top of one another
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(10, 1, figsize = (10, 30))
    for i, row in df.iterrows():
        if i > 9:
            break
        # plot the z-scores
        axs[i].hist(row, bins = 100)
        axs[i].set_title(f"Species {i}")
    # add more vertical space between each plot. do not use tight layout
    plt.subplots_adjust(hspace = 0.5)
    plt.savefig("values.pdf")
    # close
    plt.close()

    # now find the things that define the species
    # convert the column "taxid_list" to a list of ints. use eval
    cdf["taxid_list"] = cdf["taxid_list"].apply(eval)
    # first collect all of the ncbi taxids from the cdf dict
    all_taxids = set()
    for i, row in cdf.iterrows():
        all_taxids.update(row["taxid_list"])

    entries = []
    # For each taxid, get the pairs that are unique to this taxid
    for taxid in all_taxids:
        # get the samples that have this taxid
        inindex  = cdf[cdf["taxid_list"].apply(lambda x: taxid in x)].index
        # if there are at least two samples, we should continue to analyze this
        if len(inindex) < 2:
            continue
        # get everything that isn't in inindex
        outindex = cdf[~cdf.index.isin(inindex)].index
        # if the size of this is zero, we don't analyze it
        if len(outindex) == 0:
            continue
        # get all the pairs in the "in" samples. We should make a set with the cdf["lowest_1_percent_columns"] column
        inpairs = set()
        for index in inindex:
            inpairs.update(cdf.loc[index, "lowest_1_percent_columns"])
        outpairs = set()
        for index in outindex:
            outpairs.update(cdf.loc[index, "lowest_1_percent_columns"])
        # get the things that are unique to this taxid
        unique_pairs = inpairs - outpairs
        # for each unique pair, get the median value of the pairs for the inindex rows
        unique_pair_entries = [{"pair": pair, "mean_value": df.loc[inindex, pair].mean(),
                                "median_value": df.loc[inindex, pair].median(),
                                # for occupancy, count the things that aren't nan
                                "occupancy": df.loc[inindex, pair].notna().sum() / len(inindex)} \
                               for pair in unique_pairs]
        # make this into a df so we can sort
        unique_pair_df = pd.DataFrame(unique_pair_entries)
        # log of low occupancy is a smaller number. 1/mean value means larger numbers are smaller.
        # This means we need to sort descending, to get larger values first
        occweight = 1
        unique_pair_df['composite_score'] = np.log1p((unique_pair_df['occupancy']+0.000000000000000001) * occweight) + np.log1p(1/(unique_pair_df['mean_value'] + 1))
        # sort by composite score
        unique_pair_df = unique_pair_df.sort_values("composite_score", ascending = False)

        # sort these by the median value
        entries.append({"taxid": taxid,
                        "num_samples_in_taxid":      len(inindex),
                        "num_samples_outside_taxid": len(outindex),
                        "num_unique_pairs": len(unique_pairs),
                        "unique_pairs": [(int(row["pair"]), int(row["mean_value"]), row["occupancy"]) \
                                           for i, row in unique_pair_df.iterrows()]
                        })
    # convert entries to a df
    entries_df = pd.DataFrame(entries)
    # save this to a gzipped tsv
    entries_df.to_csv("unique_pairs.tsv.gz", sep = "\t", index = False, compression = "gzip")

def main():
    args = parse_args()
    #process_coo_file(sampledffile, ALGcomboixfile, coofile,
    #                 dfoutfilepath, missing_value_as = 9999999999)
    process_coo_file(args.sample_df_path, args.coo_combination_path, args.coo_path, "test.df")

if __name__ == "__main__":
    main()