#!/usr/bin/env python

"""
This takes in a list of RBH files from multiple species and constructs a rbh database file.
This is useful in the case where we have limited .rbh files that only have a few columns.

The reason for writing this script was to accommodate alternative ways of introducing homologous loci into various odp pipelines.
In one case, the user has CNEs all identified with the common coordinates of one genome, but has no predefined orthologs.

rbh                 gene_group  color   locus_id            scaffold_id  position
grp10:86431_86432   grp1        #F5733  grp10:86431_86432   scaf13       39380201
grp23:73015_73016   grp2        #3A02C  grp23:73015_73016   scaf24       56209189
grp12:62466_62467   grp3        #178B4  grp12:62466_62467   scaf23       21323678
grp16:31990_31990   grp4        #E31AC  grp16:31990_31990   scaf19       41377561
grp17:44757_44757   grp5        #6A39A  grp17:44757_44757   scaf17       76107270
grp13:75073_75074   grp6        #F7F00  grp13:75073_75074   scaf12       136956227

The program will:
1. Load in all of the RBH files
2. Determine the sample names in them
3. Construct a new database rbh file with the appropriate columns
"""

import argparse
import pandas as pd
import os
import sys

def parse_args():
    """
    Args to parse:
    - rbh_direc:  Directory containing all the rbh files
    - ALGsetName: Name of the ALGset, just something unique like "BCnSSimakov2022", "Schultz2022", etc.
                   This name will be used as a prefix for the output filenames as well.
    """
    parser = argparse.ArgumentParser(description='Construct a RBH database file from multiple RBH files.')
    parser.add_argument('-d', '--rbh_direc',
                        type=str,
                        help='Directory containing all the rbh files')
    parser.add_argument('-A', '--ALGsetName',
                        type=str,
                        help='Name of the ALGset')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    print("Looking in the following directory for .rbh files: {}".format(args.rbh_direc), file = sys.stderr)

    # 1. Load in the absolute paths of all of the RBH files
    rbh_files = [os.path.join(args.rbh_direc, f) for f in os.listdir(args.rbh_direc) if f.endswith('.rbh')]
    print("We have found the following rbh files:", file = sys.stderr)
    for f in rbh_files:
        print(f, file = sys.stderr)
    # load all of these into pandas dataframes
    rbh_dfs = [pd.read_csv(f, sep='\t') for f in rbh_files]

    # 2. Determine the sample names in them
    rbh_df_dict = {}
    for i in range(len(rbh_dfs)):
        columns = rbh_dfs[i].columns
        # get the unique entries of the columns that end in "_gene", "_scaf", "_pos"
        cols_of_interest = [x for x in columns if x.endswith('_gene') or x.endswith('_scaf') or x.endswith('_pos')]
        # make sure there is only 1 _ chat in the column name. If not, there is something wrong
        if any([col.count('_') != 1 for col in cols_of_interest]):
            print(f"Error: Column names in {rbh_files[i]} are not formatted correctly. They can only have one '_' char. Exiting.")
            print(cols_of_interest)
            sys.exit()
        # For this to work, there can only be one '_' char in the string
        uniques = set([col.split('_')[0] for col in cols_of_interest])
        # There should only be one unique entry in the set - that is the intention of this script
        # See the docstring for more info.
        if len(uniques) > 1:
            print(f"Error: More than one unique sample name found in {rbh_files[i]}. Exiting.")
        # We now use this information to to construct a new database file
        unique = uniques.pop()
        rbh_df_dict[unique] = rbh_dfs[i]
    del rbh_dfs

    # 3. Construct a new database rbh file with the appropriate columns
    # Start by trying to join two of the dataframes
    # get two of the keys for the dataframe
    sample_keys = list(rbh_df_dict.keys())
    # merge all of the columns on "rbh" "gene_group", "color". Allow missing values
    merged_df = pd.merge(rbh_df_dict[sample_keys[0]], rbh_df_dict[sample_keys[1]], on=["rbh", "gene_group", "color"], how='outer')
    # now merge the rest of the dataframes
    for i in range(2, len(sample_keys)):
        merged_df = pd.merge(merged_df, rbh_df_dict[sample_keys[i]], on=["rbh", "gene_group", "color"], how='outer')
    # the only columns should be "rbh", "gene_group", "color",  and a set of columns ending in "_scaf", "_pos", "_gen"
    for thiscol in merged_df.columns:
        if (thiscol not in ["rbh", "gene_group", "color"]) and (thiscol.split("_")[1] not in ["scaf", "pos", "gene"]):
            print(f"Error: Column {thiscol} should not be in the dataframe. Should only be 'rbh', 'gene_group', 'color', and columns ending in '_scaf', '_pos', '_gene'. Exiting.")
            sys.exit()

    # This merged df is the dabase rbh file that we need
    # Check that all the columns in rbh are unique
    if len(merged_df["rbh"].unique()) != len(merged_df):
        print("Error: The 'rbh' column is not unique. Exiting.")
        sys.exit()

    # Check that the gene_group has the same color for each row
    # Do this by getting unique entries of gene_group and color
    unique_color_gene_group = merged_df.groupby("gene_group")["color"].unique()
    if any([len(x) > 1 for x in unique_color_gene_group]):
        print("Error: The 'gene_group' column has more than one color for a single gene_group. Exiting.")
        print("Here are the gene_groups with more than one color:")
        print(unique_color_gene_group[unique_color_gene_group.apply(len) > 1])
        sys.exit()

    # Now we are sure that this file is legal and can be used by other programs
    # Write this to a file. Do not allow columns ending in _pos to write floats
    output_file = f"{args.ALGsetName}_rbh_database.rbh"
    # This is fine, but usually we would want to do this with all of the columns.
    # Right now we don't because pandas is having trouble with large dataframes.
    merged_df[["rbh", "gene_group", "color"]].to_csv(output_file, sep='\t', index=False, float_format='%.0f')

    # Now we need to generate a pseudo genome file
    merged_chrom = merged_df.copy()
    merged_chrom = merged_chrom[["rbh", "gene_group"]]
    merged_chrom["strand"] = "."
    merged_chrom["start"]  = -1
    merged_chrom["stop"]   = -1

    # randomize the index
    merged_chrom = merged_chrom.sample(frac=1).reset_index(drop=True)
    gb = merged_chrom.groupby("gene_group")
    # for each group
    # use the randomized order to add indices to the start and stop.
    # Increment by 1, stop is 1 more than start
    for name, group in gb:
        start = 1
        for i, row in group.iterrows():
            merged_chrom.at[i, "start"] = start
            merged_chrom.at[i, "stop"] = start + 1
            start += 1
    # convert this back to a df, sorted by gene_group then start+stop
    merged_chrom = merged_chrom.sort_values(by=["gene_group", "start", "stop"]).reset_index(drop=True)
    # Replace all the missing values with empty strings
    merged_chrom = merged_chrom.fillna('')
    # save this as a .chrom file, named the same as the previous rbh file, no header
    chrom_file = f"{args.ALGsetName}_rbh_database.chrom"
    merged_chrom.to_csv(chrom_file, sep='\t', index=False, header=False)

    # make a directory using the {args.ALGsetName}_rbh
    # save the rbh and chrom files in this directory
    os.makedirs(f"{args.ALGsetName}_rbh", exist_ok=True)
    # Now that we have these coordinates for all of the gene groups, we can use this to add
    #  the ALGsetName_scaf, ALGsetName_pos, ALGsetName_gene columns to the original rbh files
    #  This will allow us to use the original rbh files in the pipeline
    for sample in sample_keys:
        # tempdf to format the column names correctly
        tempdf = merged_chrom.copy()
        tempdf = tempdf[["rbh", "gene_group", "start"]]
        tempdf = tempdf.rename(columns={"gene_group": f"{args.ALGsetName}_scaf", "start": f"{args.ALGsetName}_pos"})
        # merge this with the original df, rbh in original and ALGsetName_gene in tempdf
        original_df = rbh_df_dict[sample]
        original_df = pd.merge(original_df, tempdf, on="rbh", how='left')
        original_df[f"{args.ALGsetName}_gene"] = original_df["rbh"]
        # make sure there are no NaN values. There should be none
        if any(original_df[f"{args.ALGsetName}_scaf"].isna()):
            print(f"Error: There are NaN values in the {args.ALGsetName}_scaf column. Exiting.")
            sys.exit()
        # save this to the new directory using the name {args.ALGsetName}_{sample}.rbh
        original_df.to_csv(f"{args.ALGsetName}_rbh/{args.ALGsetName}_{sample}_customGenerated.rbh", sep='\t', index=False)

if __name__ == '__main__':
    main()